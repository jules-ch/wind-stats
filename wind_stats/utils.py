from typing import Tuple, Union

import numpy as np

from wind_stats.constants import KARMAN_CONSTANT


def wind_speed_reduction(
    x_L: float,
    y_L: float,
    x_R: float,
    y_R: float,
    h: float,
    z: float,
    z_0: float,
    C_h: float,
    wind_angle: Union[float, np.ndarray] = 0.0,
):
    """Shelter model

    WEMOD algorithm outlined in Taylor, P. A. and J. R. Salmon, 1993.

    Vectorized for wind angles.

    Parameters
    ----------

    x_L: coordinate on x-axis of the extreme left point of the obstacle.
    y_L: coordinate on y-axis of the extreme left point of the obstacle.
    x_R: coordinate on x-axis of the extreme right point of the obstacle.
    x_R: coordinate on y-axis of the extreme right point of the obstacle.
    h: height of the obstacle.
    z: height where the speed reduction is measured.
    z_0: roughness length.
    C_h: the normalized wake moment coefficient.
    wind_angle: wind direction (0 means wind from North)

    Returns
    -------

    velocity reduction ratio.


    Notes
    -----

    Typical values for the normalized wake coefficient :math:`\tilde{C_h}` :

          Nature of the obstacle        :math:`\tilde{C_h}`
    ================================    ===================
    2D fences or dense line of trees    0.8(1-φ)
    Long low square buildings           0.25 - 0.4
    3D cubical buildings                0.2 - 0.35
    Vertical cylinder                   0.2
    """

    # is_not_scalar = isinstance(wind_angle, (list, tuple, np.ndarray))

    wind_angle = np.atleast_1d(wind_angle)

    alpha_L = np.arctan2(x_L, y_L)
    alpha_R = np.arctan2(x_R, y_R)

    alpha_R += 2 * np.pi if alpha_R == 0 else 0

    # At most 0.1° intervals
    intervals = int(abs(alpha_R - alpha_L) // np.radians(0.1))
    delta_beta = abs(alpha_R - alpha_L) / intervals

    beta_i = np.linspace(alpha_L + delta_beta / 2, alpha_R - delta_beta / 2, intervals)
    psi = np.radians(wind_angle)

    psi, beta_i = np.meshgrid(psi, beta_i)

    R_i = ((x_R * y_L) - (x_L * y_R)) / (
        (x_R - x_L) * np.cos(beta_i) - (y_R - y_L) * np.sin(beta_i)
    )
    Γ = 12.1875

    if (x_R - x_L) == 0:
        gamma = np.pi / 2 if (y_L > y_R) > 0 else -np.pi / 2
    elif (y_L - y_R) == 0:
        gamma = np.pi if (x_L > x_R) > 0 else 0
    else:
        gamma = np.arctan((y_L - y_R) / (x_R - x_L))

    width_i = (
        2
        * R_i
        * np.sin(delta_beta / 2)
        * (1 / np.cos(beta_i - gamma))
        * np.abs(np.cos(psi - gamma))
    )

    x_i, y_i = R_i * np.cos(psi - beta_i), R_i * np.sin(psi - beta_i)

    n = 1 / 7

    c_a = np.power(np.log((h + z_0) / z_0) / (2 * (KARMAN_CONSTANT ** 2)), 1 / (n + 2))

    a_f = 0.5
    a_g = 0.67 * c_a ** 1.5

    with np.errstate(invalid="ignore"):
        eta_i = (z / h) * np.power(x_i / h, -1 / (n + 2))
        lambda_i = (y_i / h) * np.power(x_i / h, -0.5)

    # Spreading coefficient as a gaussian
    F = 1 / (np.sqrt(2 * np.pi) * a_f) * np.exp(-(lambda_i ** 2) / (2 * a_f ** 2))
    G = c_a * eta_i * np.exp(-a_g * eta_i ** 1.5)

    with np.errstate(invalid="ignore"):
        wind_reductions = (
            Γ
            * C_h
            * (width_i / h)
            * np.power(x_i / h, -1.5)
            * G
            * F
            * (np.log((h + z_0) / z_0) / np.log((z + z_0) / z_0))
        )

    return np.nansum(wind_reductions, axis=0)
