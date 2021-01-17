import logging
from typing import TYPE_CHECKING, List, Tuple, Union

import numpy as np
from pint.quantity import Quantity
from scipy import integrate
from scipy.interpolate import interp1d

from wind_stats.constants import AIR_DENSITY
from wind_stats.gwa_reader import get_gwc_data, get_weibull_parameters
from wind_stats.stats import kde_distribution, weibull
from wind_stats.units import units
from wind_stats.utils import vertical_wind_profile, wind_power

if TYPE_CHECKING:  # pragma: no cover
    import xarray as xr

logger = logging.getLogger(__name__)


class PowerCurve:
    """Power curve.

    Parameters
    ----------
    wind_speed: `pint.Quantity`
        wind speed data
    wind_speed: `pint.Quantity`
        power data
    """

    @units.check(None, "[speed]", " [power]")
    def __init__(self, wind_speed: Quantity, power: Quantity) -> None:
        self.wind_speed = wind_speed
        self.power = power

    def __call__(self, x: Union[float, np.ndarray, Quantity]) -> Quantity:
        """Linear interpolation on the power curve.

        Parameters
        ----------
        x: array_like
            wind speed(s) at which to compute power.

        Returns
        -------
        power

        Notes
        -----
        Values outside the defined curve range will return 0W.
        """
        if isinstance(x, Quantity):
            wind_speed = self.wind_speed.m_as(x)
            x = x.m
        else:
            wind_speed = self.wind_speed.m

        interpolate = interp1d(
            wind_speed, self.power.m, bounds_error=False, fill_value=0
        )
        power = interpolate(x) * self.power.units
        return power


class WindDistribution:
    """Wind distribution."""

    def __init__(self, distribution) -> None:
        self.distribution = distribution

    def __repr__(self) -> str:
        mean = self.mean_wind_speed
        dist_type = self.distribution.dist.name
        return f"<{self.__class__.__name__}>(type: {dist_type}, mean: {mean})"

    @classmethod
    def from_data(
        cls,
        wind_speed_data,
        roughness_length: Union[float, List[float]],
        measurement_height: float,
        height: float,
    ) -> "WindDistribution":
        """Create KDE WindDistribution based on measurement data."""
        # scale wind date with vertical log profile
        if measurement_height == height:
            data = wind_speed_data
        else:
            data = vertical_wind_profile(
                height, roughness_length, measurement_height, wind_speed_data
            )

        distribution = kde_distribution(data, name="KDE")
        # freeze distribution to have the same signature as other dists
        frozen_dist = distribution.freeze()
        return cls(frozen_dist)

    @classmethod
    def from_gwc(
        cls,
        gwc_dataset: "xr.Dataset",
        roughness_length: Union[float, List[float]],
        height: float,
    ) -> "WindDistribution":
        """Create Weibull WindDistribution based on GWC file dataset."""
        A, k, _ = get_weibull_parameters(gwc_dataset, roughness_length, height)
        distribution = weibull(A, k)
        return cls(distribution)

    @classmethod
    def weibull(cls, A: float, k: float) -> "WindDistribution":
        """Create Weibull WindDistribution."""
        distribution = weibull(A, k)
        return cls(distribution)

    def pdf(self, x: float) -> float:
        """Probability density function."""
        return self.distribution.pdf(x)

    @property
    def mean_wind_speed(self) -> Quantity:
        return self.distribution.mean() * units("m/s")

    def moment(self, n: int) -> float:
        """Get n-raw moment of the distribution."""
        return self.distribution.moment(n)


class Site:
    """Site location.

    Attributes
    ----------
    latitude
    longitude
    distribution: WindDistribution
    """

    def __init__(
        self,
        latitude: float,
        longitude: float,
        distribution: WindDistribution,
        air_density: Union[Quantity, float, None] = None,
    ) -> None:
        self.latitude = latitude
        self.longitude = longitude
        self.distribution = distribution

        if air_density is not None:
            if isinstance(air_density, Quantity):
                self.air_density = air_density.to("kg/m**3")
            else:
                self.air_density = units.Quantity(air_density, "kg/m**3")
        self.air_density = air_density or AIR_DENSITY

    def __repr__(self) -> str:
        return (
            f"<{self.__class__.__name__}>\n"
            f"GPS Coordinates: latitude:{self.latitude}, longitude: {self.longitude}"
        )

    @classmethod
    def create_gwa_data(
        cls,
        latitude: float,
        longitude: float,
        roughness_length: Union[float, List[float]],
        height: float,
    ) -> "Site":
        """Create Site with Global Wind Atlas Data.

        Retrieve GWA data & initiate Site with Wind distribution.
        """
        gwc_data = get_gwc_data(latitude, longitude)
        wind_distribution = WindDistribution.from_gwc(
            gwc_data, roughness_length, height
        )
        site = cls(latitude, longitude, wind_distribution)
        return site

    @property
    def mean_wind(self) -> Quantity:
        return self.distribution.mean_wind_speed

    @property
    def mean_power_density(self) -> Quantity:
        return self.get_mean_power_density()

    def get_mean_power_density(self) -> Quantity:
        r"""Get mean power density in W/m².

        Notes
        -----
        Mean power density :

        .. math::
            \bar{P}_{density} = \frac{1}{2} \cdot \rho \cdot \int_{0}^{\infty}[v^3 pdf(v)] dv

        """
        return (0.5 * self.air_density.m * self.distribution.moment(3)) * units(
            "W/m**2"
        )


class WindTurbine:
    def __init__(
        self,
        name: str,
        power_curve: Tuple[Quantity, Quantity],
        diameter: float,
        height: float,
    ) -> None:
        self.name = name
        self.power_curve = PowerCurve(*power_curve)
        self.diameter = diameter * units.m
        self.rated_power = max(self.power_curve.power)
        self.hub_height = height * units.m

    def __repr__(self) -> str:
        return (
            f"<{self.__class__.__name__}>("
            f"{self.name}, {self.rated_power.to_compact()}, "
            f"height:{self.hub_height}, diameter:{self.diameter})"
        )

    @property
    def rotor_area(self) -> Quantity:
        """Swept rotor area."""
        return np.pi * self.diameter ** 2 / 4

    def get_power_coefficients(self) -> Tuple[np.ndarray, np.ndarray]:
        r"""Get Cp coefficients for wind speeds defined in the power curve.

        Notes
        -----
        The power coefficient :math:`C_p` is the ratio between the power
        extracted & the theoretical wind power available going through the
        swept area.

        .. math::
            C_p = \frac{P}{\frac{1}{2}\rho A V^{2}}

        """
        wind_speeds = self.power_curve.wind_speed
        available_wind_power = wind_power(self.rotor_area, self.power_curve.wind_speed)
        cp = self.power_curve.power / available_wind_power

        return wind_speeds.m_as("m/s"), cp.m_as(units.dimensionless)

    def get_mean_power(self, site: Site) -> Quantity:
        r"""Mean power output.

        .. math::
           \bar{P} = \int_{0}^{\infty}[P(v) \cdot pdf(v)] dv

        Parameters
        ----------
        site: Site

        Returns
        -------
        `pint.Quantity`

        """
        distribution = site.distribution
        wind_speeds = self.power_curve.wind_speed

        def f(wind_speed: float) -> float:
            return distribution.pdf(wind_speed) * self.power_curve(wind_speed).m

        mean_power = integrate.quad(
            f,
            min(wind_speeds).m,
            max(wind_speeds).m,
            points=wind_speeds.m,
            limit=max(50, len(wind_speeds)),
        )[0]
        return mean_power * self.power_curve.power.u

    @units.check(None, None, "[time]")
    def get_energy_production(self, site: Site, time: Quantity) -> Quantity:
        r"""Calculate energy output over a period of time.

        Parameters
        ----------
        site: Site
            Site where the wind turbine is located.
        time: `pint.Quantity`
            period of time the result energy is produced.

        Notes
        -----
        The energy produced is the integration of the power production on
        the wind speed probability density function over time:

        .. math::
           E = \bar{P} \cdot t

        .. math::
           E = \int_{0}^{\infty}[P(v) \cdot pdf(v) \cdot t] dv

        See Also
        --------
        get_annual_energy_production:
            Get energy output over a year.
        """
        mean_power = self.get_mean_power(site)

        energy = time * mean_power

        return energy.to("Wh")

    def get_annual_energy_production(self, site: Site) -> Quantity:
        """Get annual energy production.

        Integrate wind frequency with power curve & annual hours.

        Parameters
        ----------
        site: Site
            Site where the wind turbine is located.

        Returns
        -------
        energy: `pint.Quantity`
            Annual energy to be expected in Wh.

        See Also
        --------
        get_energy_production:
            Get energy production over any period of time.
        """
        return self.get_energy_production(site, 1 * units.year)
