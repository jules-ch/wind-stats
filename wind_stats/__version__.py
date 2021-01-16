try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:  # Can remove when we require Python > 3.7
    from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version(__package__)
except PackageNotFoundError:
    __version__ = "unknown"
