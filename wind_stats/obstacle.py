"""Module implementing Shelter model for basic obstacles."""

import abc
from typing import List, Tuple

import numpy as np

from .geometry import azimuth_to_cartesian_angle, rotate, translate
from .utils import wemod_wind_speed_reduction


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
    coords: Tuple[List[float], List[float]]

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

        x = [-length / 2, length / 2, length / 2, -length / 2]
        y = [width / 2, width / 2, -width / 2, -width / 2]

        angle = azimuth_to_cartesian_angle(orientation)
        self.coords = rotate(
            translate(np.array([x, y]), x_center, y_center), angle, (x_center, y_center)
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

        self.coords = translate(np.array([x, y]), x_center, y_center)

    @property
    def normalized_wake_coefficient(self):
        return 0.4 * (1 - 0.5)
