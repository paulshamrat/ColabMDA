from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from colabmda.core.ions import parse_ion_spec

# Known kinase substrate peptide consensus presets
SUBSTRATE_PRESETS: dict[str, dict[str, str]] = {
    "eb2": {
        "name": "eb2",
        "seq": "KRPSSAPSG",  # Native human EB2 Ser222 consensus motif (RPxSA)
        "description": "End-Binding Protein 2 (Ser222 phospho-acceptor)",
    },
    "eb1": {
        "name": "eb1",
        "seq": "KRPTSAPSG",  # Thr/Ser alternative
        "description": "End-Binding Protein 1 consensus motif",
    },
    "amph1": {
        "name": "amph1",
        "seq": "KRPSSGPPS",  # Amphiphysin-1 synaptic vesicle motif
        "description": "Amphiphysin-1 consensus motif",
    },
}


@dataclass
class SystemSpec:
    """Democratized, multi-component specification of any kinase assembly."""

    protein: str = "cdkl5"
    variant: str = "WT"
    ligand: str | Path | None = None  # "ATP", "ADP", "GDP", "ANP", or Path to .sdf
    ions: dict[str, int] = field(default_factory=dict)  # {"MG": 2}
    phos_residues: list[str] = field(default_factory=list)  # ["Y171"]
    peptide_name: str | None = None  # "eb2", "eb1", "amph1"
    peptide_seq: str | None = None  # e.g. "KRPSSAPSG"

    def __post_init__(self):
        # Normalize variant
        if self.variant:
            self.variant = self.variant.strip()

        # Parse ion dictionary if passed as string
        if isinstance(self.ions, str):
            self.ions = parse_ion_spec(self.ions)

        # Normalize ligand string
        if isinstance(self.ligand, str):
            self.ligand = self.ligand.strip().upper()
            if self.ligand in {"NONE", "0"}:
                self.ligand = None

        # Auto-populate substrate sequence from preset if missing
        if self.peptide_name and not self.peptide_seq:
            p_norm = self.peptide_name.lower().strip()
            if p_norm in SUBSTRATE_PRESETS:
                self.peptide_seq = SUBSTRATE_PRESETS[p_norm]["seq"]

    @property
    def is_apo(self) -> bool:
        """True if kinase is unliganded (no nucleotide and no active-site ions)."""
        return self.ligand is None and not self.ions and not self.peptide_name

    @property
    def is_phosphorylated(self) -> bool:
        """True if any residues are phosphorylated (e.g. pY171)."""
        return bool(self.phos_residues)

    @property
    def state_tag(self) -> str:
        """Generate a standardized, self-describing directory tag."""
        if self.is_apo:
            return "apo"

        parts: list[str] = []

        # 1. Ligand part
        if self.ligand:
            if isinstance(self.ligand, Path):
                parts.append(self.ligand.stem.lower())
            else:
                parts.append(str(self.ligand).lower())

        # 2. Ion part
        if self.ions:
            ion_tokens = []
            for ion_sym, count in sorted(self.ions.items()):
                ion_tokens.append(f"{count}{ion_sym.lower()}")
            parts.append("_".join(ion_tokens))

        # 3. Peptide part
        if self.peptide_name:
            parts.append(self.peptide_name.lower())
        elif self.peptide_seq:
            parts.append("pep")

        return "_".join(parts) if parts else "apo"

    @classmethod
    def from_tag(cls, tag: str, variant: str = "WT", protein: str = "cdkl5") -> SystemSpec:
        """Construct a SystemSpec from a standard biological tag.

        Supported tags:
            - "apo": Unliganded monomer
            - "atpmg" / "atp_1mg": Ground holo state (ATP + 1 Mg)
            - "atp_2mg": Pre-catalytic complex (ATP + 2 Mg)
            - "adp_1mg": Post-catalytic product (ADP + 1 Mg)
            - "gdp_1mg": Crystal resting state (GDP + 1 Mg)
            - "atp_2mg_eb2": Active ternary complex (pY171 + ATP + 2 Mg + EB2)
            - "atp_2mn_eb2": Active ternary complex with Manganese (pY171 + ATP + 2 Mn + EB2)
        """
        clean_tag = tag.strip().lower()

        if clean_tag == "apo":
            return cls(protein=protein, variant=variant, ligand=None, ions={})

        if clean_tag in {"atpmg", "atp_mg", "atp_1mg"}:
            return cls(protein=protein, variant=variant, ligand="ATP", ions={"MG": 1})

        if clean_tag == "atp_2mg":
            return cls(protein=protein, variant=variant, ligand="ATP", ions={"MG": 2})

        if clean_tag == "adp_1mg":
            return cls(protein=protein, variant=variant, ligand="ADP", ions={"MG": 1})

        if clean_tag == "gdp_1mg":
            return cls(protein=protein, variant=variant, ligand="GDP", ions={"MG": 1})

        if clean_tag in {"atp_2mg_eb2", "atp_2mg_py171_eb2", "atp_2mg_py171_peb2"}:
            return cls(
                protein=protein,
                variant=variant,
                ligand="ATP",
                ions={"MG": 2},
                phos_residues=["Y171"],
                peptide_name="eb2",
                peptide_seq=SUBSTRATE_PRESETS["eb2"]["seq"],
            )

        if clean_tag == "atp_2mn_eb2":
            return cls(
                protein=protein,
                variant=variant,
                ligand="ATP",
                ions={"MN": 2},
                phos_residues=["Y171"],
                peptide_name="eb2",
                peptide_seq=SUBSTRATE_PRESETS["eb2"]["seq"],
            )

        # Dynamic parsing for composed tags: e.g. "atp_2mg_eb1", "adp_1mg_amph1", "adp_2mg_py171"
        tokens = clean_tag.split("_")
        ligand_val = None
        ions_val: dict[str, int] = {}
        peptide_val = None
        phos_val: list[str] = []

        for token in tokens:
            if token in {"atp", "adp", "gdp", "gtp", "anp"}:
                ligand_val = token.upper()
            elif re.match(r"^\d+[a-z]{1,2}$", token):
                # Ion token like "2mg" or "1mn"
                m = re.match(r"^(\d+)([a-z]{1,2})$", token)
                if m:
                    count, sym = int(m.group(1)), m.group(2).upper()
                    ions_val[sym] = count
            elif token in {"py171", "y171"}:
                if "Y171" not in phos_val:
                    phos_val.append("Y171")
            elif token in SUBSTRATE_PRESETS:
                peptide_val = token
                if "Y171" not in phos_val:
                    phos_val.append("Y171")  # Substrate complex implies active phosphorylated kinase
            elif token.startswith("p") and token[1:] in SUBSTRATE_PRESETS:
                peptide_val = token[1:]
                if "Y171" not in phos_val:
                    phos_val.append("Y171")

        return cls(
            protein=protein,
            variant=variant,
            ligand=ligand_val,
            ions=ions_val,
            phos_residues=phos_val,
            peptide_name=peptide_val,
        )
