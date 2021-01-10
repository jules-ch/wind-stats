import logging
from typing import TYPE_CHECKING, List, Tuple, Union

import numpy as np
from pint.quantity import Quantity
from scipy import stats
from scipy import integrate
from scipy.interpolate import interp1d

from wind_stats.gwa_reader import get_gwc_data, get_weibull_parameters
from wind_stats.stats import kde_distribution
from wind_stats.utils import vertical_wind_profile

from .constants import AIR_DENSITY
from .units import units

if TYPE_CHECKING:  # pragma: no cover
    import xarray as xr

logger = logging.getLogger(__name__)


@units.check("[area]", "[speed]")
def wind_power(area: Quantity, wind_speed: Quantity) -> Quantity:
    """Calculate available wind power

    Parameters
    ----------
    area : `pint.Quantity`
        swept surface
    wind_speed: `pint.Quantity`
        wind speed

    Returns
    -------
    `pint.Quantity`
        The available wind power
    """
    return (0.5 * AIR_DENSITY * area * wind_speed ** 3).to("W")


def weibull(A: float, k: float) -> stats.weibull_min:
    r"""Weibull distribution with A, k parameters.

    .. math::
        f(x) = c \frac{x}{A}^{c-1} \exp(-\frac{x}{A}^c)
    """
    return stats.weibull_min(k, scale=A)


class PowerCurve:
    """Power curve

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

    def __call__(self, x) -> Quantity:
        """Linear interpolation

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
    """Wind distribution"""

    def __init__(self, distribution) -> None:
        self.distribution = distribution

    def __repr__(self) -> str:
        mean = self.mean_wind_speed
        dist_type = self.distribution.dist.name
        return f"{self.__class__.__name__}(type: {dist_type}, mean: {mean})"

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
        """Create Weibull WindDistribution based on GWC file dataset"""
        A, k, f = get_weibull_parameters(gwc_dataset, roughness_length, height)
        distribution = weibull(A, k)
        return cls(distribution)

    @classmethod
    def weibull(cls, A: float, k: float) -> "WindDistribution":
        """Create Weibull WindDistribution"""
        distribution = weibull(A, k)
        return cls(distribution)

    def pdf(self, x: float) -> float:
        """Probability density function"""
        return self.distribution.pdf(x)

    @property
    def mean_wind_speed(self) -> Quantity:
        return self.distribution.mean() * units("m/s")

    def get_power_density(self, wind_speed: float) -> Quantity:
        return 0.5 * AIR_DENSITY * self.pdf(wind_speed)

    def moment(self, n: int) -> float:
        return self.distribution.moment(n)


class Site:
    def __init__(
        self, latitude: float, longitude: float, distribution: WindDistribution
    ) -> None:
        self.latitude = latitude
        self.longitude = longitude
        self.distribution = distribution

    def __repr__(self):
        return f"<{self.__class__.__name__}\n GPS Coordinates: latitude:{self.latitude}, longitude: {self.longitude}>"

    @classmethod
    def create_gwa_data(
        cls,
        latitude: float,
        longitude: float,
        roughness_length: Union[float, List[float]],
        height: float,
    ) -> "Site":
        """Create Site with Global Wind Atlas Data

        Retrieve GWA data & initiate Site with Wind distribution.
        """

        gwc_data = get_gwc_data(latitude, longitude)
        wind_distribution = WindDistribution.from_gwc(
            gwc_data, roughness_length, height
        )
        return cls(latitude, longitude, wind_distribution)

    @property
    def mean_wind(self) -> Quantity:
        return self.distribution.mean_wind_speed

    @property
    def mean_power_density(self) -> Quantity:
        return self.get_mean_power_density()

    def get_mean_power_density(self) -> Quantity:
        r"""Get mean power density in W/m²

        .. math::
            \bar{P} = \frac{1}{2}\rho A^3 * \Gamma(1+\frac{1}{3})

        """

        return (0.5 * AIR_DENSITY.m * self.distribution.moment(3)) * units("W/m**2")


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
            f"WindTurbine("
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

    def annual_energy_distribution(self, site, wind_speeds=np.linspace(0, 25)):
        distribution = site.distribution

        return (
            distribution.pdf(wind_speeds)
            * self.power_curve(wind_speeds)
            * (1 * units.year).to("hours")
        )

    def get_mean_power(self, site: Site) -> Quantity:
        r"""Mean power output

        .. math::
           \bar{P} = \int_{0}^{\infty}[pdf(v) \cdot P(v)] dv
        """
        distribution = site.distribution
        wind_speeds = self.power_curve.wind_speed

        def f(wind_speed):
            return distribution.pdf(wind_speed) * self.power_curve(wind_speed).m

        mean_power, _ = integrate.quad(
            f,
            min(wind_speeds).m,
            max(wind_speeds).m,
            points=wind_speeds.m,
            limit=max(50, len(wind_speeds)),
        )
        return mean_power * self.power_curve.power.u

    def get_annual_energy_production(self, site: Site) -> Quantity:
        r"""Get annual energy production.

        Integrate wind frequency with power curve & annual hours.

        Parameters
        ----------
        site: Site
            Site where the wind turbine is located.

        Returns
        -------
        energy: `pint.Quantity`
            Annual energy to be expected in Wh.

        Notes
        -----
        The energy producted is the integration of the power production on
        the wind  speed probability density function over time:

        .. math::
           E = \int_{0}^{\infty}[pdf(v) \cdot P(v) \cdot t] dv

        """

        annual_hours = (1 * units.year).to("hours")

        mean_power = self.get_mean_power(site)
        annual_energy_production = annual_hours * mean_power

        return annual_energy_production
