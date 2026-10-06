"""Which processes couple two radiation bands directly, i.e. make the band-band block of the Jacobian non-diagonal."""

from dataclasses import dataclass

import sympy as sp

from ..symbols import n_, x_
from .band import BandSet


@dataclass(frozen=True)
class BandStructure:
    """couplings: (process, band species whose density the row depends on, band species of the row)"""

    bands: tuple
    couplings: tuple

    @property
    def diagonal(self):
        """True if no band's row depends on another band: the bands then couple only through the matter variables,
        and the band-band Jacobian block is diagonal"""
        return not self.couplings

    @property
    def pairs(self):
        return sorted({(src, dst) for _, src, dst in self.couplings})

    def metadata(self):
        return {"band_jacobian_diagonal": self.diagonal, "band_couplings": [list(p) for p in self.pairs]}


def _closure(symbols, deps):
    """symbols with everything they depend on through deps (symbol -> free symbols of its definition)"""
    out, todo = set(symbols), list(symbols)
    while todo:
        for d in deps.get(todo.pop(), ()):
            if d not in out:
                out.add(d)
                todo.append(d)
    return out


def band_structure(processes, bands):
    """BandStructure of a Model (its processes after its rules, with its derived quantities and intermediates
    expanded) or of a sequence of processes.

    bands: BandSet or sequence of band species names
    """
    species = tuple(bands.species if isinstance(bands, BandSet) else bands)
    deps = {}
    if hasattr(processes, "derived"):
        deps.update({sp.Symbol(k) if isinstance(k, str) else k: sp.sympify(v).free_symbols
                     for k, v in processes.derived.items()})
        deps.update({w: sp.sympify(W).free_symbols for w, W in processes.intermediates})
        atoms = processes.subprocesses
    else:
        atoms = [a for p in processes for a in p.subprocesses]
    symbols = {s: {n_(s), x_(s)} for s in species}
    couplings = []
    for p in atoms:
        net = p._network
        for dst in species:
            if dst not in net:
                continue
            free = _closure(sp.sympify(dict.__getitem__(net, dst).rhs).free_symbols, deps)
            for src in species:
                if src != dst and free & symbols[src]:
                    couplings.append((p.name, src, dst))
    return BandStructure(species, tuple(couplings))
