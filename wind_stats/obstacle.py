"""Module implementing Shelter model for basic obstacles."""

import abc
from typing import Tuple, Union

import numpy as np

from .constants import KARMAN_CONSTANT
from .geometry import azimuth_to_cartesian_angle, rotate, translate


class Obstacle:
    """Obstacle base class

    Attributes
    ----------

    coords: array_like
        coordinates of the obstacle as x, y
    height:
        height of the obstacle
    """

    height: float
    coords: Tuple[np.ndarray, np.ndarray]

    @property
    @abc.abstractmethod
    def normalized_wake_coefficient(self) -> float:
        raise NotImplementedError()

    def get_coordinates(self, x_offset=0.0, y_offset=0.0):
        coords = self.coords
        x, y = self.coords

        centroid_x = np.mean(x)
        centroid_y = np.mean(y)

        rotation = -np.arctan2(centroid_y - y_offset, centroid_x - x_offset)

        new_coords = rotate(
            translate(coords, -x_offset, -y_offset), rotation, (0, 0), use_radians=True
        )

        angles = np.arctan2(new_coords[1], new_coords[0])

        x_L, y_L = coords[0][np.argmax(angles)], coords[1][np.argmax(angles)]
        try:
            x_R, y_R = (
                coords[0][np.argmax(angles) + 2],
                coords[1][np.argmax(angles) + 2],
            )
        except IndexError:
            x_R, y_R = (
                coords[0][np.argmax(angles) - 2],
                coords[1][np.argmax(angles) - 2],
            )

        return x_L, y_L, x_R, y_R, angles

    def get_wind_reduction(self, z, z_0, wind_angle, x_offset=0.0, y_offset=0.0):
        c_h = self.normalized_wake_coefficient
        coords = self.coords
        x, y = self.coords

        centroid_x = np.mean(x)
        centroid_y = np.mean(y)

        rotation = -np.arctan2(centroid_y - y_offset, centroid_x - x_offset)

        new_coords = rotate(
            translate(coords, -x_offset, -y_offset), rotation, (0, 0), use_radians=True
        )

        angles = np.arctan2(new_coords[1], new_coords[0])

        x_L, y_L = coords[0][np.argmax(angles)], coords[1][np.argmax(angles)]
        x_R, y_R = coords[0][np.argmin(angles)], coords[1][np.argmin(angles)]

        return wemod_wind_speed_reduction(
            x_L,
            y_L,
            x_R,
            y_R,
            self.height,
            z,
            z_0,
            c_h,
            wind_angle=wind_angle,
            x_offset=x_offset,
            y_offset=y_offset,
        )


class Building(Obstacle):
    def __init__(self, center, length, width, height, orientation) -> None:
        super().__init__()
        x_center, y_center = center
        self.length = length
        self.width = width
        self.height = height
        self.orientation = orientation

        x = np.array([-length / 2, length / 2, length / 2, -length / 2])
        y = np.array([width / 2, width / 2, -width / 2, -width / 2])

        angle = azimuth_to_cartesian_angle(orientation)
        self.coords = rotate(
            translate((x, y), x_center, y_center), angle, (x_center, y_center)
        )

    @property
    def normalized_wake_coefficient(self):
        if self.width / self.height <= 2:
            return 0.4
        else:
            return 0.35


class Tree(Obstacle):
    def __init__(self, center, radius, height) -> None:
        super().__init__()
        theta = np.linspace(0, 2 * np.pi, 100)
        x, y = radius * np.cos(theta), radius * np.sin(theta)
        x_center, y_center = center
        self.height = height

        self.coords = translate((x, y), x_center, y_center)

    @property
    def normalized_wake_coefficient(self):
        return 0.4 * (1 - 0.5)


def wemod_wind_speed_reduction(
    x_L: float,
    y_L: float,
    x_R: float,
    y_R: float,
    h: float,
    z: float,
    z_0: float,
    C_h: float,
    wind_angle: Union[float, np.ndarray] = 0.0,
    x_offset: float = 0.0,
    y_offset: float = 0.0,
) -> float:
    r"""Shelter model

    WEMOD algorithm outlined in Taylor, P. A. and J. R. Salmon, 1993.

    Vectorized for wind angles.

    Parameters
    ----------

    x_L: float
        coordinate on x-axis of the extreme left point of the obstacle.
    y_L: float
        coordinate on y-axis of the extreme left point of the obstacle.
    x_R: float
        coordinate on x-axis of the extreme right point of the obstacle.
    x_R: float
        coordinate on y-axis of the extreme right point of the obstacle.
    h: float
        height of the obstacle.
    z: height where the speed reduction is measured.
    z_0: roughness length.
    C_h: the normalized wake moment coefficient.
    wind_angle: azimuth wind direction (0 means wind from North)
    x_offset: float
        offset on x-axis to the location where reduction is computed.
    y_offset: float
        offset on x-axis to the location where reduction is computed

    Returns
    -------

    velocity reduction ratio : float


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

    # Offset coordinates
    x_L = x_L - x_offset
    x_R = x_R - x_offset
    y_L = y_L - y_offset
    y_R = y_R - y_offset

    wind_angle = np.atleast_1d(wind_angle)

    alpha_L = np.arctan2(x_L, y_L)
    alpha_R = np.arctan2(x_R, y_R)

    if alpha_R < alpha_L:
        alpha_R += 2 * np.pi

    alpha_R += 2 * np.pi if alpha_R == 0 else 0

    # At most 0.1° intervals
    intervals = int(abs(alpha_R - alpha_L) // np.radians(0.1))
    delta_beta = abs(alpha_R - alpha_L) / intervals

    beta_i = np.linspace(
        alpha_L + delta_beta / 2, alpha_R - delta_beta / 2, intervals
    ) % (2 * np.pi)
    psi = np.radians(wind_angle)

    psi, beta_i = np.meshgrid(psi, beta_i)

    radii = ((x_R * y_L) - (x_L * y_R)) / (
        (x_R - x_L) * np.cos(beta_i) - (y_R - y_L) * np.sin(beta_i)
    )
    Γ = 12.1875

    if (x_R - x_L) == 0:
        gamma = np.pi / 2 if (y_L > y_R) > 0 else 3 * np.pi / 2
    elif (y_L - y_R) == 0:
        gamma = np.pi if (x_L > x_R) > 0 else 0
    else:
        gamma = np.arctan((y_L - y_R) / (x_R - x_L)) % (2 * np.pi)

    width_i = (
        2
        * radii
        * np.sin(delta_beta / 2)
        * (1 / np.abs(np.cos(beta_i - gamma)))
        * np.abs(np.cos(psi - gamma))
    )

    x_i, y_i = radii * np.cos(psi - beta_i), radii * np.sin(psi - beta_i)

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

