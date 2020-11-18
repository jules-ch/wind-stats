from pint import UnitRegistry

UREG = UnitRegistry()
UREG.default_format = ".3f~P"

AIR_DENSITY = 1.225 * UREG("kg/m**3")