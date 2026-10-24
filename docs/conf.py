"""Sphinx configuration for the README and two reference pages."""

from importlib.metadata import metadata as package_metadata
from importlib.metadata import version as package_version

project = "ColabMDA"
author = package_metadata("colabmda").get("Author", "Shamrat Kumar Paul")
copyright = f"2026, {author}"
release = package_version("colabmda")

extensions = ["myst_parser"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
myst_enable_extensions = ["colon_fence"]
myst_heading_anchors = 3
