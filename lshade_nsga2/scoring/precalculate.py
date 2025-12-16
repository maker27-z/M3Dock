from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

from ..atom import Atom
from ..model import InteractingPair
from .weighted_terms import WeightedTerms


def _distance_squared(a: Atom, b: Atom) -> float:
    dx = a.coords.x - b.coords.x
    dy = a.coords.y - b.coords.y
    dz = a.coords.z - b.coords.z
    return dx * dx + dy * dy + dz * dz


@dataclass
class Precalculate:
    """
    Very lightweight placeholder for C++ precalculate:
    - stores weighted_terms
    - exposes eval_fast(type_pair_index, r2)
    For now, type_pair_index is ignored; later we will map AD/XS types.
    """

    wt: WeightedTerms
    cutoff_sqr: float = field(init=False)

    def __post_init__(self) -> None:
        c = self.wt.cutoff()
        self.cutoff_sqr = c * c

    def eval_pair(self, a: Atom, b: Atom) -> float:
        r2 = _distance_squared(a, b)
        if r2 >= self.cutoff_sqr:
            return 0.0
        from math import sqrt

        r = sqrt(r2)
        return self.wt.eval_pair(a, b, r)


def eval_interacting_pairs(p: Precalculate, pairs: List[InteractingPair], atoms: List[Atom]) -> float:
    e = 0.0
    for ip in pairs:
        e += p.eval_pair(atoms[ip.a], atoms[ip.b])
    return e


