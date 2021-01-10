from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal
from pytest import approx
from scipy import stats

from wind_stats import GWAReader
from wind_stats.models import PowerCurve, Site, WindDistribution, WindTurbine, weibull
from wind_stats.units import units

TEST_DATA = Path(__file__).parent


class TestWindDistribution:
    def test_mean_wind_speed(self):
        distribution = weibull(6, 2)
        wind_distribution = WindDistribution(distribution)

        assert wind_distribution.mean_wind_speed == distribution.mean() * units("m/s")

    def test_from_gwc(self):
        with (TEST_DATA / "gwa3_gwc_test_file.lib").open() as f:
            dataset = GWAReader.load(f)
        wind_distribution = WindDistribution.from_gwc(dataset, 0.5, 100.0)
        assert wind_distribution

    def test_from_data(self):
        test_dist = stats.norm(loc=6)
        data = test_dist.rvs(10000)
        wind_distribution = WindDistribution.from_data(data, 0.5, 100.0, 100.0)
        assert approx(test_dist.mean(), wind_distribution.distribution.mean)
        # with scaled wind speed
        wind_distribution = WindDistribution.from_data(data, 0.5, 50.0, 100.0)
        assert wind_distribution

    def test_weibull(self):
        wind_distribution = WindDistribution.weibull(A=6, k=2)
        assert wind_distribution.distribution.dist.name == "weibull_min"
        assert wind_distribution.distribution.kwds == {"scale": 6}
        assert wind_distribution.distribution.args == (2,)

    def test_repr(self):
        distribution = weibull(6, 2)
        wind_distribution = WindDistribution(distribution)
        assert (
            repr(wind_distribution)
            == "WindDistribution(type: weibull_min, mean: 5.3174 m/s)"
        )

    def test_moment(self):
        # first moment is mean
        distribution = weibull(6, 2)
        wind_distribution = WindDistribution(distribution)
        assert wind_distribution.moment(1) == wind_distribution.distribution.mean()


class TestPowerCurve:
    def test_call(self):
        windspeed = units.Quantity(np.linspace(0, 25), "m/s")
        power = units.Quantity(2 * np.linspace(0, 25), "kW")
        powercurve = PowerCurve(windspeed, power)

        assert_array_equal(
            powercurve([1.0, 2.0, 3.0]), units.Quantity([2.0, 4.0, 6.0], "kW")
        )
        assert_array_almost_equal(
            powercurve(units.Quantity([3.6, 7.2, 10.8], "km/h")),
            units.Quantity([2.0, 4.0, 6.0], "kW"),
        )


class TestWindTurbine:

    # https://en.wind-turbine-models.com/turbines/1467-siemens-swt-3.3-130-ln
    # Power curve data
    @pytest.fixture
    def turbine(self):

        wind_speed = [
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
        ] * units("m/s")
        power = [
            43.0,
            184.0,
            421.0,
            778.0,
            1270.0,
            1905.0,
            2593.0,
            3096.0,
            3268.0,
            3297.0,
            3300.0,
            3300.0,
            3300.0,
            3300.0,
            3300.0,
            3300.0,
            3300.0,
            3300.0,
            3300.0,
            3300.0,
            3300.0,
            3300.0,
            3300.0,
        ] * units.kW

        wind_turbine = WindTurbine(
            "Siemens SWT-3.3-130 LN", (wind_speed, power), 130, 135
        )
        return wind_turbine

    def test_rotor_area(self, turbine):
        assert turbine.rotor_area.m == approx(13273.23)
        assert turbine.rotor_area.u == units("m**2")

    def test_get_mean_power(self, turbine):
        distribution = weibull(6, 2)
        wind_distribution = WindDistribution(distribution)

        site = Site(0, 0, wind_distribution)
        assert turbine.get_mean_power(site).m == approx(852.943)

    def test_get_power_coefficients(self, turbine):
        ws, cp = turbine.get_power_coefficients()

        assert_array_almost_equal(
            cp,
            np.array(
                [
                    0.196,
                    0.354,
                    0.414,
                    0.443,
                    0.455,
                    0.458,
                    0.438,
                    0.381,
                    0.302,
                    0.235,
                    0.185,
                    0.148,
                    0.120,
                    0.099,
                    0.083,
                    0.070,
                    0.059,
                    0.051,
                    0.044,
                    0.038,
                    0.033,
                    0.029,
                    0.026,
                ]
            ),
            decimal=3,
        )

    def test_get_annual_energy_production(self, turbine):
        distribution = weibull(6, 2)
        wind_distribution = WindDistribution(distribution)

        site = Site(0, 0, wind_distribution)

        assert turbine.get_annual_energy_production(site).m_as("MWh") == approx(
            7476.9041
        )
