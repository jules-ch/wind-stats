import numpy as np


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
