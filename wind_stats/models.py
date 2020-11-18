from dataclasses import dataclass, field
from typing import List

import numpy as np
from numpy.core.fromnumeric import mean
from scipy.integrate import quad
from scipy.interpolate.interpolate import interp1d
from scipy.special import gamma
from scipy.stats import rv_continuous, weibull_min

from wind_stats.constants import AIR_DENSITY

from .constants import UREG


@UREG.wraps("W", ("m**2", "m/s"))
def wind_power(area, wind_speed):
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
        self.wind_speed = wind_speed * UREG("m/s")
        self.power = power * UREG.W


class WindDistribution:
    def __init__(self, A, k) -> None:
        self.A = A
        self.k = k
        self.weibull = weibull(A, k)


@dataclass
class Site:

    latitude: float
    longitude: float
    distribution: WindDistribution
    turbines: List["WindTurbine"] = field(default_factory=list)

    @property
    def mean_wind(self):
        return self.distribution.weibull.mean() * UREG("m/s")

    @property
    def mean_power_density(self):
        return self.get_mean_power_density()

    @UREG.wraps("W/m**2", None)
    def get_mean_power_density(self):
        """Get mean power density in W/m²

        .. math::
            \bar{P} = \frac{1}{2}\rho A^3 * \Gamma(1+\frac{1}{3})

        """

        A = self.distribution.A
        k = self.distribution.k
        return 0.5 * AIR_DENSITY.m * A ** 3 * gamma(1 + 3 / k)


class WindTurbine:
    def __init__(self, power_curve, diameter) -> None:
        self.power_curve = power_curve
        self.diameter = diameter * UREG.m
        self.site = None
        self.rated_power = max(self.power_curve.power)

    def set_site(self, site: Site):
        self.site = site

    @property
    def rotor_area(self):
        return np.pi * self.diameter ** 2 / 4

    def get_power_coefficients(self):
        return self.power_curve.power / wind_power(
            self.rotor_area, self.power_curve.wind_speed
        )

    def get_mean_power(self):
        """Mean power output

        ..math
           \int_{0}^{\infty}(pdf(v) \cdot P(v)) dv
        """

        distribution = self.site.distribution
        wind_speeds = self.power_curve.wind_speed

        def f(wind_speed):
            power_function = interp1d(
                self.power_curve.wind_speed.m, self.power_curve.power.m
            )
            return distribution.weibull.pdf(wind_speed) * power_function(wind_speed)

        mean_power, _ = quad(
            f, min(wind_speeds).m, max(wind_speeds).m, points=wind_speeds.m
        )
        return mean_power * UREG.W

    def get_annual_energy_production(self):
        """
        Get annual energy production.

        Integrate wind frequency with power curve & annual hours.

        ..math
           \int_{0}^{\infty}(pdf(v) \cdot P(v) \cdot t) dv

        Returns:
            Annual energy to be expected in Wh.
        """

        if not self.site:
            raise ValueError("No site has been set")
        annual_hours = (1 * UREG.year).to("hours")

        mean_power = self.get_mean_power()
        annual_energy_production = annual_hours * mean_power

        return annual_energy_production
