try:
    from importlib.metadata import PackageNotFoundError, version  # type: ignore
except ImportError:
    from importlib_metadata import PackageNotFoundError, version  # type: ignore

try:
    __version__ = version(__package__)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
