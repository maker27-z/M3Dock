from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

from .atom import Atom
from .model import Model, Ligand
from .types import Vec3


def _parse_float_field(line: str, start: int, end: int) -> float:
    """Parse a fixed-width float field (1-based indices as in C++ code)."""
    # convert to Python 0-based slicing
    s = max(0, start - 1)
    e = min(len(line), end)
    field = line[s:e].strip()
    return float(field) if field else 0.0


def _parse_string_field(line: str, start: int, end: int) -> str:
    s = max(0, start - 1)
    e = min(len(line), end)
    return line[s:e].strip()


def _parse_pdbqt_atom_line(line: str) -> Atom:
    """
    Python counterpart of parse_pdbqt_atom_string in parse_pdbqt.cpp,
    but simplified: we parse coords, charge and AD type.
    """
    x = _parse_float_field(line, 31, 38)
    y = _parse_float_field(line, 39, 46)
    z = _parse_float_field(line, 47, 54)
    charge = 0.0
    # charge is in columns 69-76 if present
    try:
        charge = _parse_float_field(line, 69, 76)
    except ValueError:
        charge = 0.0

    # AD type is typically in columns 78-79 (like " C", "OA", etc.)
    ad_type = _parse_string_field(line, 78, 79)
    name = _parse_string_field(line, 13, 16)
    resname = _parse_string_field(line, 18, 20)
    resid = _parse_string_field(line, 23, 26)
    chain = _parse_string_field(line, 22, 22)

    return Atom(
        coords=Vec3(x=x, y=y, z=z),
        ad_type=ad_type,
        charge=charge,
        name=name,
        residue_name=resname,
        residue_id=resid,
        chain_id=chain,
    )


def parse_ligand_pdbqt(path: Path) -> Model:
    """
    Simplified ligand PDBQT parser:
    - constructs a Model with one ligand
    - all atoms from ROOT/BRANCH blocks are flattened into a single linear ligand
    - torsions and rotatable bonds are NOT yet reconstructed (will be added later)
    """
    atoms: List[Atom] = []
    with path.open("r", encoding="utf-8") as f:
        in_model = False
        for line in f:
            if line.startswith("MODEL"):
                in_model = True
                continue
            if line.startswith("ENDMDL"):
                break
            if not in_model and (line.startswith("ATOM  ") or line.startswith("HETATM")):
                # Vina input usually has a single MODEL, but we are tolerant here
                in_model = True
            if not in_model:
                continue
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                atoms.append(_parse_pdbqt_atom_line(line))

    if not atoms:
        raise ValueError(f"No atoms parsed from ligand PDBQT: {path}")

    model = Model(atoms=atoms, grid_atoms=[], ligands=[])
    # For now treat whole atom list as one ligand, and approximate degrees_of_freedom
    dof = 6  # translation (3) + rotation (3); torsions will be added later
    lig = Ligand(begin=0, end=len(atoms), degrees_of_freedom=dof)
    model.ligands.append(lig)
    model.finalize()
    return model


def parse_receptor_pdbqt(path: Path) -> Model:
    """
    Simplified receptor PDBQT parser:
    - all ATOM/HETATM records are treated as rigid grid_atoms
    - no flex residues (flexible side chains) are handled yet
    """
    grid_atoms: List[Atom] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                grid_atoms.append(_parse_pdbqt_atom_line(line))

    if not grid_atoms:
        raise ValueError(f"No atoms parsed from receptor PDBQT: {path}")

    model = Model(atoms=[], grid_atoms=grid_atoms, ligands=[])
    model.finalize()
    return model


def parse_complex(receptor_path: Path, ligand_path: Path) -> Model:
    """
    Helper that mirrors C++ parse_bundle(rigid_name, flex_name_opt, ligand_names)
    for the simplest case: rigid receptor + single ligand.
    """
    rec = parse_receptor_pdbqt(receptor_path)
    lig = parse_ligand_pdbqt(ligand_path)

    # merge into one Model: receptor grid_atoms + ligand atoms
    model = Model()
    model.grid_atoms = rec.grid_atoms
    model.atoms = lig.atoms
    model.ligands = lig.ligands
    model.finalize()
    return model



