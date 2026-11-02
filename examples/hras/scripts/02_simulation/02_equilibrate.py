"""Equilibrate one reviewed hras variant on a GPU runtime."""

import argparse
from pathlib import Path

from colabmda.core.proteins import simulate

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--variant", default="WT")
args = parser.parse_args()
simulate(Path(__file__).resolve().parents[3], "hras", "equilibrate", args.variant)
