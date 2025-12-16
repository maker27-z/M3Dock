from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

from ..atom import Atom
from ..model import Model
from .terms import Terms


@dataclass
class WeightedTerms:
    """
    Simplified Python equivalent of C++ weighted_terms:
    - holds a Terms instance
    - holds weights (external + conf-independent)
    """

    terms: Terms
    weights: List[float] = field(default_factory=list)

    def cutoff(self) -> float:
        return self.terms.max_cutoff()

    def eval_pair(self, a: Atom, b: Atom, r: float) -> float:
        acc = 0.0
        usable = self.terms.usable_terms
        usable_weights = self.weights[: len(usable)]
        for w, term in zip(usable_weights, usable):
            acc += w * term.eval(a, b, r)
        return acc

    def conf_independent(self, m: Model, energy: float) -> float:
        offset = len(self.terms.usable_terms)
        weights_ci = self.weights[offset:]
        x = energy
        for term, w in zip(self.terms.conf_independent_terms, weights_ci):
            x = term.eval(m, x, [w])
        return x


