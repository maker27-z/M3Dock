from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class Vec3:
    x: float
    y: float
    z: float

    def as_tuple(self) -> Tuple[float, float, float]:
        return self.x, self.y, self.z


@dataclass
class GridDim:
    begin: float
    end: float
    n: int


@dataclass
class GridDims:
    """Python equivalent of C++ grid_dims (3 dimensions)."""

    x: GridDim
    y: GridDim
    z: GridDim

    def as_list(self) -> List[GridDim]:
        return [self.x, self.y, self.z]



