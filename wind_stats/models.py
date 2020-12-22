import logging
from typing import List, Tuple

import numpy as np
from pint.quantity import Quantity
from scipy.integrate import quad
from scipy.interpolate.interpolate import interp1d
from scipy.special import gamma
from scipy.stats import weibull_min

from .constants import AIR_DENSITY
from .units import units

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


def weibull(A: float, k: float) -> weibull_min:
    r"""Weibull distribution with A, k parameters.

    .. math::
        f(x) = c \frac{x}{A}^{c-1} \exp(-\frac{x}{A}^c)
    """
    return weibull_min(k, scale=A)


class PowerCurve:
    """Power curve

    Power in Watts
    Wind speeds in meter/second
    """

    def __init__(self, wind_speed: List[float], power: List[float]) -> None:
        self.wind_speed = wind_speed * units("m/s")
        self.power = power * units.W


class WindDistribution:
    """Wind distribution as a Weibull distribution."""

    def __init__(self, A: float, k: float) -> None:
        self.A = A
        self.k = k
        self.weibull = weibull(A, k)

    def pdf(self, x: float) -> float:
        return self.weibull.pdf(x)

    @classmethod
    def from_mean_wind(cls, mean_wind: float, k: float) -> "WindDistribution":
        A = mean_wind / gamma(1 + 1 / k)
        return cls(A, k)

    @property
    def mean_wind_speed(self) -> Quantity:
        return self.weibull.mean() * units("m/s")

    def get_power_density(self, wind_speed: float):
        return 0.5 * AIR_DENSITY * self.pdf(wind_speed)

    def get_max_power_density_wind_speed(self):
        A = self.A
        k = self.k
        return A * ((k + 2) / k) ** (1 / k) * units("m/s")


class Site:

    latitude: float
    longitude: float
    distribution: WindDistribution

    def __init__(self, latitude, longitude, distribution) -> None:
        self.latitude = latitude
        self.longitude = longitude
        self.distribution = distribution

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

        A = self.distribution.A
        k = self.distribution.k
        return (0.5 * AIR_DENSITY.m * A ** 3 * gamma(1 + 3 / k)) * units("W/m**2")


class WindTurbine:
    def __init__(
        self,
        name: str,
        power_curve: Tuple[List[float], List[float]],
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

        return wind_speeds.m, cp.m

    def annual_energy_distribution(self, site, wind_speeds=np.linspace(0, 25)):
        distribution = site.distribution
        power_function = interp1d(
            self.power_curve.wind_speed.m, self.power_curve.power.m, bounds_error=False
        )

        return (
            distribution.pdf(wind_speeds)
            * power_function(wind_speeds)
            * self.power_curve.power.u
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
            power_function = interp1d(
                self.power_curve.wind_speed.m, self.power_curve.power.m
            )
            return distribution.pdf(wind_speed) * power_function(wind_speed)

        mean_power, _ = quad(
            f,
            min(wind_speeds).m,
            max(wind_speeds).m,
            points=wind_speeds.m,
            limit=max(50, len(wind_speeds)),
        )
        return mean_power * units.W

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
