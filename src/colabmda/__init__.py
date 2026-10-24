"""ColabMDA public package metadata."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("colabmda")
except PackageNotFoundError:
    __version__ = "0+uninstalled"
