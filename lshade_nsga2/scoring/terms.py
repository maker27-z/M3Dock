from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List

from ..atom import Atom
from ..model import Model

# ---- basic math utilities ----


def gaussian(x: float, width: float) -> float:
    from math import exp

    return exp(-(x / width) ** 2)


def slope_step(x_bad: float, x_good: float, x: float) -> float:
    if x_bad < x_good:
        if x <= x_bad:
            return 0.0
        if x >= x_good:
            return 1.0
    else:
        if x >= x_bad:
            return 0.0
        if x <= x_good:
            return 1.0
    return (x - x_bad) / (x_good - x_bad)


# ---- term definitions (simplified) ----


@dataclass
class Term:
    name: str


@dataclass
class Usable(Term):
    cutoff: float
    eval_fn: Callable[[Atom, Atom, float], float]

    def eval(self, a: Atom, b: Atom, r: float) -> float:
        if r >= self.cutoff:
            return 0.0
        return self.eval_fn(a, b, r)


@dataclass
class ConfIndependent(Term):
    size: int
    eval_fn: Callable[[Model, float, List[float]], float]

    def eval(self, m: Model, x: float, weights: List[float]) -> float:
        return self.eval_fn(m, x, weights)


@dataclass
class Terms:
    usable_terms: List[Usable] = field(default_factory=list)
    conf_independent_terms: List[ConfIndependent] = field(default_factory=list)

    def max_cutoff(self) -> float:
        return max((t.cutoff for t in self.usable_terms), default=0.0)


# ---- predefined simple vina-like terms (approximation of everything.cpp defaults) ----


def make_default_vina_terms() -> Terms:
    """
    Build a minimal set of Vina-like usable terms based on everything.cpp defaults:
      - gauss(0, 0.5)
      - gauss(3, 2.0)
      - repulsion(0.0)
      - hydrophobic(0.5, 1.5)
      - non_dir_h_bond(-0.7, 0)

    NOTE: This is a simplified implementation; hydrophobic / H-bond detection
    uses heuristic based on AD types.
    """

    def hydrophobic_type(ad: str) -> bool:
        return ad in {"C", "A", "N", "NA", "SA", "P", "S", "Cl", "Br", "I", "F"}

    def hbond_possible(ad1: str, ad2: str) -> bool:
        donors = {"HD"}
        acceptors = {"OA", "NA"}
        return (ad1 in donors and ad2 in acceptors) or (ad2 in donors and ad1 in acceptors)

    def gauss(offset: float, width: float, cutoff: float) -> Usable:
        def fn(a: Atom, b: Atom, r: float) -> float:
            return gaussian(r - offset, width)

        return Usable(name=f"gauss(o={offset},w={width})", cutoff=cutoff, eval_fn=fn)

    def repulsion(offset: float, cutoff: float) -> Usable:
        def fn(a: Atom, b: Atom, r: float) -> float:
            d = r - offset
            return 0.0 if d > 0 else d * d

        return Usable(name=f"repulsion(o={offset})", cutoff=cutoff, eval_fn=fn)

    def hydrophobic_term(good: float, bad: float, cutoff: float) -> Usable:
        def fn(a: Atom, b: Atom, r: float) -> float:
            if hydrophobic_type(a.ad_type) and hydrophobic_type(b.ad_type):
                return slope_step(bad, good, r)
            return 0.0

        return Usable(name="hydrophobic", cutoff=cutoff, eval_fn=fn)

    def non_dir_h_bond(good: float, bad: float, cutoff: float) -> Usable:
        def fn(a: Atom, b: Atom, r: float) -> float:
            if hbond_possible(a.ad_type, b.ad_type):
                return slope_step(bad, good, r)
            return 0.0

        return Usable(name="non_dir_h_bond", cutoff=cutoff, eval_fn=fn)

    cutoff = 8.0
    terms = Terms()
    terms.usable_terms.extend(
        [
            gauss(0.0, 0.5, cutoff),
            gauss(3.0, 2.0, cutoff),
            repulsion(0.0, cutoff),
            hydrophobic_term(0.5, 1.5, cutoff),
            non_dir_h_bond(-0.7, 0.0, cutoff),
        ]
    )
    return terms



