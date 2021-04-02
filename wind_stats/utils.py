"""Utilities functions."""


import numpy as np
from pint import Quantity

from wind_stats.constants import ISA_AIR_DENSITY, Rd, Rv
from wind_stats.units import units


def vertical_wind_profile(
    height,
    roughness_length,
    reference_height: float,
    reference_wind_speed: float,
    zero_displacement_plane: float = 0.0,
) -> float:
    """Scale wind speed with Log Wind profile.

    Parameters
    ----------
    height:
        height where wind speed is evaluated.
    roughness_length:
        roughness length parameter in meters.
    reference_height: float
        reference height where reference wind speed is measured.
    reference_wind_speed: float
        reference wind speed measured at `reference_height`.
    zero_displacement_plane: float
        the height at which the mean velocity is zero due to large
        obstacles such as buildings/canopy.

    Returns
    -------
    wind_speed
    """
    wind_speed = (
        reference_wind_speed
        * (np.log(height - zero_displacement_plane) / roughness_length)
        / (np.log(reference_height - zero_displacement_plane) / roughness_length)
    )

    return wind_speed


@units.check("[area]", "[speed]", "[density]")
def wind_power(
    area: Quantity, wind_speed: Quantity, air_density: Quantity = ISA_AIR_DENSITY
) -> Quantity:
    """Calculate available wind power.

    Parameters
    ----------
    area : `pint.Quantity`
        swept surface
    wind_speed: `pint.Quantity`
        wind speed
    air_density: `pint.Quantity`
        air_density (default to ISA air density)

    Returns
    -------
    `pint.Quantity`
        The available wind power
    """
    return (0.5 * air_density * area * wind_speed ** 3).to("W")


@units.check("[temperature]", "[pressure]", None)
def calculate_air_density(
    temperature: Quantity, pressure: Quantity, relative_humidity: float
):

    """Air density function based on revised formula for the density of moist air."""
    sat_pressure_0c = 6.112 * units.millibar

    saturation_vapor_pressure = sat_pressure_0c * np.exp(
        17.67
        * (temperature - 273.15 * units.kelvin)
        / (temperature - 29.65 * units.kelvin)
    )

    vapor_pressure = saturation_vapor_pressure * relative_humidity
    dry_air_pressure = pressure - vapor_pressure

    air_density = (dry_air_pressure / (Rd * temperature)) + (
        vapor_pressure / (Rv * temperature)
    )

    return air_density.to("kg/m**3")
