"""Geometry operations module.

Basic affine 2D transformations."""
from math import cos, pi, sin
from typing import Tuple

import numpy as np


def affine_2d_transformation(
    coordinates: Tuple[np.ndarray, np.ndarray], matrix: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:

    """Numpy based affine 2D transformation.

    [x']   | a  b xoff | [x]
    [y'] = | d  e yoff | [y]
    [1 ]   | 0  0   1  | [1]

    Parameters
    ----------
    coordinates: tuple
        x, y coordinates
    matrix:
        transformation matrix [3, 3]

    Returns
    -------
    x, y = coordinates
    """
    if matrix.shape != (3, 3):
        raise ValueError("2D transformation matrix must be of shape [3, 3]")

    x, y = coordinates
    vector = np.array([x, y, np.ones(len(x))])
    x_new, y_new, _ = np.matmul(matrix, vector)
    return x_new, y_new


def rotate(
    coordinates,
    angle: float = 0,
    origin: Tuple[float, float] = (0, 0),
    use_radians=False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Returns rotated coordinates.
    The angle of rotation can be specified in either degrees (default) or
    radians by setting ``use_radians=True``.
    Positive angles are counter-clockwise and negative are clockwise rotations.
    The point of origin can be a keyword 'center' for the bounding box
    center (default), 'centroid' for the geometry's centroid, a Point object
    or a coordinate tuple (x0, y0).

    The affine transformation matrix for 2D rotation is:
      | cos(θ) -sin(θ) xoff |
      | sin(θ)  cos(θ) yoff |
      |   0       0      1  |

    where the offsets are calculated from the origin Point(x0, y0):
        xoff = x0 - x0 * cos(r) + y0 * sin(r)
        yoff = y0 - x0 * sin(r) - y0 * cos(r)
    """
    if not use_radians:
        angle = angle * pi / 180.0
    cosp = cos(angle)
    sinp = sin(angle)

    x0, y0 = origin

    xoff = x0 - x0 * cosp + y0 * sinp
    yoff = y0 - x0 * sinp - y0 * cosp

    matrix = np.array(
        [
            [cosp, -sinp, xoff],
            [sinp, cosp, yoff],
            [0, 0, 1],
        ]
    )
    return affine_2d_transformation(coordinates, matrix)


def translate(
    coordinates: Tuple[np.ndarray, np.ndarray], xoff=0.0, yoff=0.0
) -> Tuple[np.ndarray, np.ndarray]:
    r"""Returns a translated geometry shifted by offsets along each dimension.
    The general 2D affine transformation matrix for translation is:
        | 1  0  xoff |
        | 0  1  yoff |
        | 0  0   1   |
    """
    matrix = np.array(
        [
            [1, 0, xoff],
            [0, 1, yoff],
            [0, 0, 1],
        ]
    )

    return affine_2d_transformation(coordinates, matrix)


def azimuth_to_cartesian_angle(azimuth: float, radians: bool = False) -> float:
    """Convert cartographical azimuth angle to cartesian angle.

    Parameters
    ----------
    azimuth: float
        cartographical azimuth angle 0 is North, 180 is South when using
        degrees.

    Returns
    -------

    angle: float
        cartesian angle between [0...360°] or [0...2π] based on `radians`
        parameter.

    Examples
    --------
    >>> azimuth_to_regular_angle(90)
    0
    >>> azimuth_to_regular_angle(0)
    90
    >>> azimuth_to_regular_angle(270)
    180
    >>> azimuth_to_regular_angle(270, radians=True)
    3.141592653589793
    """

    cartesian_angle = (90 - azimuth) % 360

    if radians:
        cartesian_angle = np.deg2rad(cartesian_angle)

    return cartesian_angle
