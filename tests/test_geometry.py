import numpy as np
from pytest import approx

from wind_stats.geometry import azimuth_to_cartesian_angle, rotate, translate


def test_translate():
    x = [0]
    y = [0]
    assert translate((x, y)) == (x, y)
    assert translate((x, y), xoff=5, yoff=-5) == ([5], [-5])


def test_rotate():
    x = [5]
    y = [0]
    assert rotate((x, y)) == (x, y)
    assert rotate((x, y), angle=90)[0] == approx([0])
    assert rotate((x, y), angle=90)[1] == approx([5])


def test_azimuth_to_cartesian_angle():
    assert azimuth_to_cartesian_angle(0.0) == 90.0
    assert azimuth_to_cartesian_angle(45.0) == 45.0
    assert azimuth_to_cartesian_angle(90.0) == 0.0
    assert azimuth_to_cartesian_angle(0.0, radians=True) == np.pi / 2
    assert azimuth_to_cartesian_angle(180.0) == 270.0
    assert azimuth_to_cartesian_angle(270.0) == 180.0
    assert azimuth_to_cartesian_angle(360.0) == 90.0
