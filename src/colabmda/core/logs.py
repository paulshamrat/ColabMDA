"""Verified stage records in execution logs, with read support for older runs."""

import hashlib
import os
import tempfile
from pathlib import Path

_HEADER = "\n# ColabMDA output verification\n"


def log_directory(folder):
    """Mirror an owning data folder beneath its project's logs directory."""
    folder = Path(folder).resolve()
    for ancestor in [folder, *folder.parents]:
        if ancestor.name == "data":
            return ancestor.parent / "logs" / folder.relative_to(ancestor)
    return folder / "logs"


def _sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class StageLog:
    """Store status and metadata in one stage log; never create a DONE file.

    Path-like methods keep the validation at existing workflow boundaries.
    Older flags remain readable until their owning stage is explicitly migrated.
    """

    def __init__(self, folder, legacy_name):
        self.folder = Path(folder).resolve()
        self.legacy = self.folder / legacy_name
        self.stage = legacy_name.removesuffix("_DONE").lower()
        self.path = log_directory(self.folder) / (self.stage + ".log")

    def __fspath__(self):
        return str(self.path if self.path.exists() or not self.legacy.exists() else self.legacy)

    def __str__(self):
        return os.fspath(self)

    def read_text(self, *args, **kwargs):
        return Path(os.fspath(self)).read_text(*args, **kwargs)

    def exists(self):
        if not self.path.exists():
            return self.legacy.is_file()
        try:
            record = self.path.read_text().rsplit(_HEADER, 1)[1]
            if record.splitlines()[-1] != "status: complete":
                return False
            artifacts = [
                line.split("\t", 2) for line in record.splitlines() if line.startswith("artifact\t")
            ]
            if not artifacts:
                return False
            for _, digest, name in artifacts:
                relative = Path(name)
                if relative.is_absolute() or ".." in relative.parts:
                    return False
                path = self.folder / relative
                if not path.resolve().is_relative_to(self.folder):
                    return False
                if (
                    path.is_symlink()
                    or not path.is_file()
                    or not path.stat().st_size
                    or _sha(path) != digest
                ):
                    return False
            return True
        except (OSError, IndexError, ValueError):
            return False

    is_file = exists

    def _outputs(self):
        folder, stage = self.folder, self.stage

        def first(candidates):
            return next((p for p in candidates if p.is_file()), candidates[0])

        if stage == "structure":
            if (folder / "WT.pdb").is_file() and (folder / "input.json").is_file():
                return [folder / "WT.pdb", folder / "input.json"]
            return [
                first([folder / "raw/template.pdb", folder / "template.pdb"]),
                first(
                    [
                        folder / "raw/uniprot.fasta",
                        folder / "uniprot.fasta",
                        *sorted(folder.glob("*.fasta")),
                    ]
                ),
            ]
        if stage == "model":
            direct = folder / (folder.name + ".pdb")
            return [direct] if direct.is_file() else sorted(folder.glob("models/*/*.pdb"))
        if stage == "assembly":
            return [first([folder / (folder.name + ".pdb"), folder / "complex.pdb"])]
        if stage == "prep":
            return [folder / name for name in ["system.prmtop", "system.inpcrd", "solvated.pdb"]]
        if stage == "eq":
            return [
                folder / "system.xml",
                first([folder / "eq.pdb", folder / "equilibrated.pdb"]),
                first([folder / "npt.state.xml", folder / "eq_final.chk"]),
            ]
        if stage == "merge":
            from colabmda.openmm.merge import find_trajectory_pair

            pair, ca = find_trajectory_pair(folder), find_trajectory_pair(folder, ca=True)
            return (
                [*pair, *ca, first([folder / "md.log", folder / "prod_full.log"])]
                if pair and ca
                else []
            )
        if stage == "analysis":
            return [
                folder / name
                for name in [
                    "rmsd.csv",
                    "rg.csv",
                    "rmsf.csv",
                    "figures/rmsd_vs_time.png",
                    "figures/rg_vs_time.png",
                    "figures/rmsf_per_residue.png",
                ]
            ]
        if stage == "aggregate":
            return [
                folder / name
                for name in [
                    "rmsd_avg.csv",
                    "rg_avg.csv",
                    "rmsf_avg.csv",
                    "figures/rmsd_aggregate.png",
                    "figures/rg_aggregate.png",
                    "figures/rmsf_aggregate.png",
                ]
            ]
        if stage == "visualize":
            return [
                folder / "visualize.pml",
                first([folder / "view.pdb", *sorted(folder.glob("*.pdb"))]),
            ]
        if stage in {"run", "md", "smoke"}:
            return sorted(folder.rglob("prod.chk")) + sorted(folder.rglob("*.dcd"))
        return [p for p in folder.glob("*.png") if p.stem.lower() == stage]

    def _write(self, content):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", dir=self.path.parent, delete=False
        ) as stream:
            temporary = Path(stream.name)
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            temporary.replace(self.path)
        finally:
            temporary.unlink(missing_ok=True)

    def write_text(self, content, *args, **kwargs):
        outputs = self._outputs()
        if not outputs or any(
            p.is_symlink() or not p.is_file() or not p.stat().st_size for p in outputs
        ):
            raise ValueError(f"Cannot complete {self.stage}: required outputs are missing or empty")
        details = str(content).split(_HEADER, 1)[0].rstrip()
        entries = "".join(f"artifact\t{_sha(p)}\t{p.relative_to(self.folder)}\n" for p in outputs)
        self._write(details + _HEADER + entries + "status: complete\n")
        if not self.exists():
            raise ValueError(
                f"Outputs changed while recording {self.stage}; previous metadata preserved"
            )
        self.legacy.unlink(missing_ok=True)
        return len(content)

    def unlink(self, missing_ok=False):
        if not self.path.exists() and not self.legacy.exists():
            if not missing_ok:
                raise FileNotFoundError(self.path)
            return
        details = self.read_text().split(_HEADER, 1)[0].rstrip()
        self._write(details + _HEADER + "status: incomplete\n")
        self.legacy.unlink(missing_ok=True)


def migrate_project_logs(root):
    """Move completed legacy flags and visual provenance into mirrored stage logs."""
    root = Path(root).resolve()
    data = root / "data"
    if any(data.rglob("*.dcd.tmp")):
        raise ValueError(
            "Unfinished production chunks remain; finish the run before migrating logs"
        )
    records = []
    for marker in sorted(data.rglob("*_DONE")):
        if marker.is_symlink():
            raise ValueError("Refusing a symlinked completion record")
        record = StageLog(marker.parent, marker.name)
        if record.path.exists() and _HEADER in record.path.read_text() and not record.exists():
            raise ValueError(f"Current stage log is incomplete or outputs changed: {record.path}")
        outputs = record._outputs()
        if not outputs or any(not p.is_file() or not p.stat().st_size for p in outputs):
            raise ValueError(f"Incomplete outputs for {marker}")
        records.append((marker, record))
    for source, record in records:
        details = source.read_text()
        if record.path.exists():
            details = record.read_text().split(_HEADER, 1)[0] + "\n" + details
        record.write_text(details)
    moved = [str(record.path.relative_to(root)) for _, record in records]
    for source in sorted(data.rglob("*.log")):
        if source.name not in {"visualize.log", "structure_snapshots.log", "model_snapshots.log"}:
            continue
        record = StageLog(source.parent, source.stem.upper() + "_DONE")
        details = source.read_text()
        if record.path.exists():
            details = record.read_text().split(_HEADER, 1)[0] + "\n" + details
        record.write_text(details)
        source.unlink()
        moved.append(str(record.path.relative_to(root)))
    return moved
