from importlib.metadata import version as _pkg_version, PackageNotFoundError

from .ehframe import Location, System
from .config import config, config2d
from .eh2d import System2D, Section2D, Fill, HollowBlock, Slab

try:
    __version__ = _pkg_version("enerhabitat")
except PackageNotFoundError:        # not installed (e.g. running from source)
    __version__ = "0.3.0"
