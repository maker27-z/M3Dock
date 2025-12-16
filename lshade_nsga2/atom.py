from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

from .types import Vec3


@dataclass
class AtomIndex:
    """Python equivalent of C++ atom_index."""

    i: int
    in_grid: bool = False


@dataclass
class Bond:
    """Simplified Python equivalent of C++ bond."""

    connected_atom_index: AtomIndex
    length: float
    rotatable: bool = False


@dataclass
class Atom:
    """
    Simplified atom combining the roles of C++ atom_base + atom.

    We keep:
    - coords       : current 3D coordinates
    - ad_type      : AutoDock atom type string (e.g. 'C', 'OA', 'HD')
    - charge       : partial atomic charge
    - name, res    : PDB/PDBQT identifiers where available
    - bonds        : list of Bond (will be populated later when we port bonding logic)
    """

    coords: Vec3
    ad_type: str
    charge: float = 0.0
    name: str = ""
    residue_name: str = ""
    residue_id: str = ""
    chain_id: str = ""
    bonds: List[Bond] = field(default_factory=list)

    def position_tuple(self) -> Tuple[float, float, float]:
        return self.coords.as_tuple()



