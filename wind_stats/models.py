import logging
from dataclasses import dataclass
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


@units.wraps("W", ("m**2", "m/s"))
def wind_power(area, wind_speed):
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
    return 0.5 * AIR_DENSITY.m * area * wind_speed ** 3


def weibull(A, k):
    """Weibull distribution with A, k parameters.

    .. math::
        f(x) = c \frac{x}{A}^{c-1} \exp(-\frac{x}{A}^c)
    """
    return weibull_min(k, scale=A)


class PowerCurve:
    """Power curve

    Power in Watts
    Wind speeds in meter/second
    """

    def __init__(self, wind_speed, power) -> None:
        self.wind_speed = wind_speed * units("m/s")
        self.power = power * units.W


class WindDistribution:
    """Wind distribution as a Weibull distribution."""

    def __init__(self, A, k) -> None:
        self.A = A
        self.k = k
        self.weibull = weibull(A, k)

    @classmethod
    def from_mean_wind(cls, mean_wind, k):
        A = mean_wind / gamma(1 + 1 / k)
        return cls(A, k)

    def power_density(self, wind_speed):
        return 0.5 * AIR_DENSITY * self.weibull.pdf(wind_speed)


@dataclass
class Site:

    latitude: float
    longitude: float
    distribution: WindDistribution

    @property
    def mean_wind(self):
        return self.distribution.weibull.mean() * units("m/s")

    @property
    def mean_power_density(self):
        return self.get_mean_power_density()

    @units.wraps("W/m**2", None)
    def get_mean_power_density(self) -> Quantity:
        """Get mean power density in W/m²

        .. math::
            \bar{P} = \frac{1}{2}\rho A^3 * \Gamma(1+\frac{1}{3})

        """

        A = self.distribution.A
        k = self.distribution.k
        return 0.5 * AIR_DENSITY.m * A ** 3 * gamma(1 + 3 / k)


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

    @property
    def rotor_area(self) -> Quantity:
        """Swept rotor area."""
        return np.pi * self.diameter ** 2 / 4

    def get_power_coefficients(self) -> Tuple[np.ndarray]:
        """Get Cp coefficients for wind speeds defined in the power curve."""

        wind_speeds = self.power_curve.wind_speed
        available_wind_power = wind_power(self.rotor_area, self.power_curve.wind_speed)
        cp = self.power_curve.power / available_wind_power

        return wind_speeds.m, cp.m

    def get_mean_power(self, site: Site) -> Quantity:
        """Mean power output

        ..math
           \int_{0}^{\infty}(pdf(v) \cdot P(v)) dv
        """
        distribution = site.distribution
        wind_speeds = self.power_curve.wind_speed

        def f(wind_speed):
            power_function = interp1d(
                self.power_curve.wind_speed.m, self.power_curve.power.m
            )
            return distribution.weibull.pdf(wind_speed) * power_function(wind_speed)

        mean_power, _ = quad(
            f, min(wind_speeds).m, max(wind_speeds).m, points=wind_speeds.m
        )
        return mean_power * units.W

    def get_annual_energy_production(self, site: Site) -> Quantity:
        """
        Get annual energy production.

        Integrate wind frequency with power curve & annual hours.

        ..math
           \int_{0}^{\infty}(pdf(v) \cdot P(v) \cdot t) dv

        Returns
        -------
            Annual energy to be expected in Wh.
        """

        annual_hours = (1 * units.year).to("hours")

        mean_power = self.get_mean_power(site)
        annual_energy_production = annual_hours * mean_power

        return annual_energy_production
