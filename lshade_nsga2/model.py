from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .atom import Atom, AtomIndex, Bond
from .types import Vec3, GridDims


@dataclass
class InteractingPair:
    """Python counterpart of C++ interacting_pair (simplified)."""

    a: int
    b: int


@dataclass
class Ligand:
    """
    Highly simplified ligand structure.
    In C++ it inherits from flexible_body and atom_range and contains:
      - degrees_of_freedom
      - interacting_pairs
      - context for writing PDBQT
    Here we keep only what we need immediately and will extend later.
    """

    begin: int
    end: int
    degrees_of_freedom: int
    pairs: List[InteractingPair] = field(default_factory=list)

    def atom_indices(self) -> range:
        return range(self.begin, self.end)


@dataclass
class Model:
    """
    Python analogue of the core of C++ model, but starting with a reduced
    feature set sufficient to:
      - hold receptor (grid_atoms) and ligand (atoms)
      - expose number of movable atoms
      - compute simple bounding box for movable atoms

    As we port more code, this class will gain:
      - bonds & mobility-based distance types
      - internal vs other interacting_pairs, grids, etc.
    """

    # movable atoms: ligand + flex (for now we treat all as movable)
    atoms: List[Atom] = field(default_factory=list)
    # rigid receptor atoms used as grid_atoms in C++
    grid_atoms: List[Atom] = field(default_factory=list)
    # ligands description (currently assuming one ligand)
    ligands: List[Ligand] = field(default_factory=list)

    # cached number of movable atoms
    _num_movable_atoms: int = 0

    def finalize(self) -> None:
        """
        After constructing Model (adding grid_atoms, atoms, ligands),
        call this to compute cached counts etc.
        """
        # For now we treat all non-grid atoms as movable
        self._num_movable_atoms = len(self.atoms)
        # Ensure ligand ranges are valid
        for lig in self.ligands:
            if lig.end > len(self.atoms):
                raise ValueError("Ligand atom range exceeds atom list size")

    # ---- basic API used later by optimizers ----

    @property
    def num_movable_atoms(self) -> int:
        return self._num_movable_atoms

    def ligand_degrees_of_freedom(self, ligand_index: int) -> int:
        return self.ligands[ligand_index].degrees_of_freedom

    def movable_coords(self) -> List[Vec3]:
        return [a.coords for a in self.atoms[: self._num_movable_atoms]]

    def movable_atoms_box(self, add_padding: float = 0.0, granularity: float = 0.375) -> GridDims:
        """
        Python version of model::movable_atoms_box with a slightly simplified
        interface: returns GridDims.
        """
        if self._num_movable_atoms == 0:
            raise ValueError("No movable atoms in model")

        coords = self.movable_coords()
        xs = [c.x for c in coords]
        ys = [c.y for c in coords]
        zs = [c.z for c in coords]

        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        min_z, max_z = min(zs), max(zs)

        min_x -= add_padding / 2
        min_y -= add_padding / 2
        min_z -= add_padding / 2
        max_x += add_padding / 2
        max_y += add_padding / 2
        max_z += add_padding / 2

        def _make_dim(lo: float, hi: float) -> Tuple[float, float, int]:
            span = hi - lo
            n = max(1, int((span / granularity) + 0.9999))
            real_span = granularity * n
            center = (lo + hi) * 0.5
            begin = center - real_span / 2.0
            end = begin + real_span
            return begin, end, n

        bx, ex, nx = _make_dim(min_x, max_x)
        by, ey, ny = _make_dim(min_y, max_y)
        bz, ez, nz = _make_dim(min_z, max_z)

        from .types import GridDim, GridDims

        return GridDims(
            x=GridDim(begin=bx, end=ex, n=nx),
            y=GridDim(begin=by, end=ey, n=ny),
            z=GridDim(begin=bz, end=ez, n=nz),
        )

    # ---- simple pair construction (placeholder until full bond/topology logic is ported) ----

    def build_intra_ligand_pairs(self) -> List[InteractingPair]:
        """
        Build all unique pairs within each ligand (i < j), skipping hydrogens
        by heuristic (AD type 'H' or 'HD').
        """
        pairs: List[InteractingPair] = []
        for lig in self.ligands:
            for i in range(lig.begin, lig.end):
                ai = self.atoms[i]
                if ai.ad_type in {"H", "HD"}:
                    continue
                for j in range(i + 1, lig.end):
                    aj = self.atoms[j]
                    if aj.ad_type in {"H", "HD"}:
                        continue
                    pairs.append(InteractingPair(a=i, b=j))
        return pairs

    def iter_ligand_atoms(self):
        for lig in self.ligands:
            for idx in range(lig.begin, lig.end):
                yield idx, self.atoms[idx]




