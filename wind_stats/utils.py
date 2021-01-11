import numpy as np
from pint import Quantity

from wind_stats.units import units
from wind_stats.constants import AIR_DENSITY


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


@units.check("[area]", "[speed]")
def wind_power(area: Quantity, wind_speed: Quantity) -> Quantity:
    """Calculate available wind power.

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
