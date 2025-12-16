from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class DockingConfig:
    # input files (PDBQT paths, relative to working dir, like original C++)
    receptor: Path
    ligand: Path
    target: Optional[Path]

    # search box
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float

    # outputs
    out: Path
    log: Optional[Path]

    # misc options
    cpu: int = 1
    seed: Optional[int] = None
    exhaustiveness: int = 8
    num_modes: int = 9
    energy_range: float = 3.0


def _parse_kv_line(line: str) -> Optional[tuple[str, str]]:
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    if "#" in line:
        line = line.split("#", 1)[0].strip()
    if not line:
        return None
    if "=" not in line:
        return None
    key, value = line.split("=", 1)
    return key.strip(), value.strip()


def load_config(path: Path) -> DockingConfig:
    """
    Parse a Vina-style config file, mirroring the semantics in main.cpp.
    Example (from example/1bzc/1bzc_config.txt):

        receptor = ../coreset_pdbqt/1bzc/1bzc_protein.pdbqt
        ligand   = ../coreset_pdbqt/1bzc/1bzc_ligand.pdbqt
        target   = ../coreset_pdbqt/1bzc/1bzc_ligand.pdbqt
        center_x = -18.643
        ...
    """
    kv: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            parsed = _parse_kv_line(raw)
            if parsed is None:
                continue
            key, value = parsed
            kv[key] = value

    def req(name: str) -> str:
        if name not in kv:
            raise ValueError(f"Missing required config key: {name}")
        return kv[name]

    base_dir = path.parent

    receptor = (base_dir / req("receptor")).resolve()
    ligand = (base_dir / req("ligand")).resolve()
    target = kv.get("target")
    target_path = (base_dir / target).resolve() if target else None

    center_x = float(req("center_x"))
    center_y = float(req("center_y"))
    center_z = float(req("center_z"))
    size_x = float(req("size_x"))
    size_y = float(req("size_y"))
    size_z = float(req("size_z"))

    out_str = kv.get("out")
    if out_str:
        out_path = (base_dir / out_str).resolve()
    else:
        # mimic default_output(ligand_name) in main.cpp: strip ".pdbqt" and add "_out.pdbqt"
        lig_name = ligand.name
        if lig_name.endswith(".pdbqt"):
            lig_name = lig_name[:-6]
        out_path = (ligand.parent / f"{lig_name}_out.pdbqt").resolve()

    log_str = kv.get("log")
    log_path = (base_dir / log_str).resolve() if log_str else None

    cpu = int(kv.get("cpu", "1"))
    exhaustiveness = int(kv.get("exhaustiveness", "8"))
    num_modes = int(kv.get("num_modes", "9"))
    energy_range = float(kv.get("energy_range", "3.0"))

    seed_val: Optional[int]
    if "seed" in kv:
        try:
            seed_val = int(kv["seed"])
        except ValueError:
            seed_val = None
    else:
        seed_val = None

    return DockingConfig(
        receptor=receptor,
        ligand=ligand,
        target=target_path,
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        out=out_path,
        log=log_path,
        cpu=cpu if cpu > 0 else 1,
        seed=seed_val,
        exhaustiveness=exhaustiveness if exhaustiveness > 0 else 1,
        num_modes=num_modes if num_modes > 0 else 1,
        energy_range=energy_range,
    )



