"""Wind stats units support module.

This module add units support for wind-stats with `pint`.

Attributes
----------
units: `pint.UnitRegistry`

"""
from pint import UnitRegistry

units = UnitRegistry()
units.default_format = ".5g~P"

units.define("@alias hour = h")
