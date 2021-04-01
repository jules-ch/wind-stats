"""Constants module.

This module provides basic access to constants.
"""

from wind_stats.units import units

KARMAN_CONSTANT = 0.4

R = 8.314462618 * units("J / mol / K")
Md = dry_air_molecular_weight = 28.96546e-3 * units("kg / mol")
Rd = dry_air_gas_constant = R / Md
Mw = water_molecular_weight = 18.015268 * units("g / mol")
Rv = water_gas_constant = R / Mw

# ISA (International Atmosphere)

ISA_TEMPERATURE = 15 * units("°C")
ISA_PRESSURE = 101.325 * units("kPa")
ISA_AIR_DENSITY = 1.225 * units("kg/m**3")