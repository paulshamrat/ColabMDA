"""Fetch the fkbp example; follow the main ColabMDA README."""

from pathlib import Path

from colabmda.core.proteins import fetch

fetch(Path(__file__).resolve().parents[3], "fkbp")
