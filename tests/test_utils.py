from wind_stats.constants import ISA_PRESSURE, ISA_TEMPERATURE
from wind_stats.utils import calculate_air_density
from pytest import approx


def test_calculate_air_density():
    air_density = calculate_air_density(ISA_TEMPERATURE, ISA_PRESSURE, 0)
    assert air_density.m == approx(1.225, abs=1e-4)
