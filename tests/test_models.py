import pytest
from pytest import approx

from wind_stats.models import Site, WindDistribution, WindTurbine, weibull
from wind_stats.units import units


class TestWindDistribution:
    def test_mean_wind_speed(self):
        distribution = weibull(6, 2)
        wind_distribution = WindDistribution(distribution)

        assert wind_distribution.mean_wind_speed == distribution.mean() * units("m/s")


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

    def test_get_power_coefficients(self, turbine):
        pass

    def test_get_mean_power(self, turbine):
        distribution = weibull(6, 2)
        wind_distribution = WindDistribution(distribution)

        site = Site(0, 0, wind_distribution)
        assert turbine.get_mean_power(site).m == approx(852.943)

    def test_get_annual_energy_production(self, turbine):
        distribution = weibull(6, 2)
        wind_distribution = WindDistribution(distribution)

        site = Site(0, 0, wind_distribution)

        assert turbine.get_annual_energy_production(site).m_as("MWh") == approx(
            7476.9041
        )
