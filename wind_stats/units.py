"""Wind stats units support module."""

from pint import UnitRegistry

units = UnitRegistry()
units.default_format = ".5g~P"

units.define("@alias hour = h")
