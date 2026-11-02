from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class IonDefinition:
    """Standardized chemical definition of an ion in molecular modeling."""

    symbol: str  # wwPDB 2-character chemical code (e.g. "MG", "MN")
    element: str  # Periodic table element name (e.g. "Mg", "Mn")
    charge: int  # Formal charge (e.g. +2, +1, -1)
    amber_unit: str  # Amber tleap unit name (e.g. "Mg+", "Mn2+", "Na+")
    openmm_name: str  # OpenMM residue/atom identifier
    radius_nm: float  # Approximate ionic radius in nanometers


# Central Ion Registry for catalytic metals and bulk solvent salts
ION_REGISTRY: dict[str, IonDefinition] = {
    "MG": IonDefinition("MG", "Mg", 2, "Mg+", "MG", 0.072),
    "MN": IonDefinition("MN", "Mn", 2, "Mn2+", "MN", 0.083),
    "CA": IonDefinition("CA", "Ca", 2, "Ca2+", "CA", 0.100),
    "ZN": IonDefinition("ZN", "Zn", 2, "Zn2+", "ZN", 0.074),
    "CO": IonDefinition("CO", "Co", 2, "Co2+", "CO", 0.075),
    "NA": IonDefinition("NA", "Na", 1, "Na+", "NA", 0.102),
    "K": IonDefinition("K", "K", 1, "K+", "K", 0.138),
    "CL": IonDefinition("CL", "Cl", -1, "Cl-", "CL", 0.181),
}


def get_ion(symbol: str) -> IonDefinition:
    """Retrieve an ion definition by its symbol (case-insensitive)."""
    norm = symbol.strip().upper()
    if norm not in ION_REGISTRY:
        supported = ", ".join(ION_REGISTRY.keys())
        raise ValueError(f"Unknown ion symbol '{symbol}'. Supported ions: {supported}")
    return ION_REGISTRY[norm]


def parse_ion_spec(spec: str | dict[str, int] | None) -> dict[str, int]:
    """Parse ion specification strings into a normalized {SYMBOL: count} dictionary.

    Examples:
        >>> parse_ion_spec("MG:2")
        {'MG': 2}
        >>> parse_ion_spec("MG:1,MN:1")
        {'MG': 1, 'MN': 1}
        >>> parse_ion_spec("0")
        {}
        >>> parse_ion_spec(None)
        {}
    """
    if spec is None:
        return {}
    if isinstance(spec, dict):
        return {k.strip().upper(): int(v) for k, v in spec.items() if int(v) > 0}

    s = spec.strip()
    if not s or s == "0" or s.lower() == "none":
        return {}

    result: dict[str, int] = {}
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            ion, count_str = part.split(":", 1)
            ion_sym = ion.strip().upper()
            get_ion(ion_sym)  # Validate against registry
            result[ion_sym] = int(count_str.strip())
        else:
            # Single bare ion implies count = 1
            ion_sym = part.upper()
            get_ion(ion_sym)
            result[ion_sym] = result.get(ion_sym, 0) + 1

    return {k: v for k, v in result.items() if v > 0}
