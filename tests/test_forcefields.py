"""Force field selection and padding resolution unit tests."""

import pytest
from openmm import unit

from colabmda.openmm.modular.em import _resolve_padding, _select_forcefield


def test_forcefield_selection():
    # Valid pairings
    ff_opc = _select_forcefield("ff19SB", "opc")
    assert ff_opc["protein_xml"] == "amber19-all.xml"
    assert ff_opc["water_xml"] == "amber19/opc.xml"

    ff_tip3p = _select_forcefield("ff14SB", "tip3p")
    assert ff_tip3p["protein_xml"] == "amber14-all.xml"
    assert ff_tip3p["water_xml"] == "amber14/tip3p.xml"

    # Invalid pairing raises ValueError
    with pytest.raises(ValueError):
        _select_forcefield("ff19SB", "nonexistent_water")


def test_resolve_padding_explicit():
    positions = [(0.0, 0.0, 0.0) * unit.nanometer, (1.0, 1.0, 1.0) * unit.nanometer]
    # Explicit valid padding
    pad = _resolve_padding(positions, "1.2")
    assert pad == 1.2

    # Under 1.0 nm cutoff raises ValueError
    with pytest.raises(ValueError):
        _resolve_padding(positions, "0.8")


def test_resolve_padding_auto():
    # Compact solute (< 8 nm span)
    positions = [(0.0, 0.0, 0.0) * unit.nanometer, (2.0, 2.0, 2.0) * unit.nanometer]
    pad = _resolve_padding(positions, "auto")
    assert pad == 1.0
