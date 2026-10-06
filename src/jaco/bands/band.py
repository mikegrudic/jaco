"""Radiation bands and the ordered band set that numbers them."""

from dataclasses import dataclass, replace

from ..declarations import Species

UNITS = ("photons", "energy")
SHAPES = ("ppl", "tracked")


@dataclass(frozen=True)
class Band:
    """A photon-energy interval [E_lo, E_hi) [eV] transported as one radiation species.

    Parameters
    ----------
    name: str
        Identifier; the species is photon_<name>.
    E_lo, E_hi: float
        Edges [eV].
    unit: str
        What the species counts: "photons" (photons per H nucleus) or "energy" (eV per H nucleus).
    shape: str
        The spectrum assumed inside the band:

        - "ppl": u_E proportional to E^slope, fixed (the piecewise-power-law "fixed slope" spectrum); slope -1 is
          E u_E = const;
        - "tracked": a dilute blackbody at a radiation temperature T_rad the host carries per cell; its projections are
          tables in T_rad.
    slope: float
        Index of u_E for a "ppl" band.
    """

    name: str
    E_lo: float
    E_hi: float
    unit: str = "energy"
    shape: str = "ppl"
    slope: float = -1.0

    def __post_init__(self):
        if not 0 < self.E_lo < self.E_hi:
            raise ValueError(f"band {self.name}: need 0 < E_lo < E_hi, got [{self.E_lo}, {self.E_hi}]")
        if self.unit not in UNITS:
            raise ValueError(f"band {self.name}: unit must be one of {UNITS}, not {self.unit!r}")
        if self.shape not in SHAPES:
            raise ValueError(f"band {self.name}: shape must be one of {SHAPES}, not {self.shape!r}")

    @property
    def species(self):
        return f"photon_{self.name}"

    def contains(self, E):
        return self.E_lo <= E < self.E_hi

    def with_slope(self, slope):
        return replace(self, slope=float(slope))

    def declaration(self):
        """The radiation Species of the band"""
        what = "photons" if self.unit == "photons" else "energy [eV]"
        return Species(self.species, "radiation", f"{what} per H nucleus in {self.E_lo:g}-{self.E_hi:g} eV", floor=0.0)


class BandSet:
    """Bands in index order; their energy intervals may not overlap, and may leave gaps (energy emitted into a gap
    escapes every band)."""

    def __init__(self, bands):
        self.bands = tuple(bands)
        names = [b.name for b in self.bands]
        if len(set(names)) != len(names):
            raise ValueError(f"duplicate band names in {names}")
        ordered = sorted(self.bands, key=lambda b: b.E_lo)
        for lo, hi in zip(ordered, ordered[1:]):
            if hi.E_lo < lo.E_hi:
                raise ValueError(f"bands {lo.name} and {hi.name} overlap")

    def __iter__(self):
        return iter(self.bands)

    def __len__(self):
        return len(self.bands)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.bands[self.index(key)]
        return self.bands[key]

    def __eq__(self, other):
        return isinstance(other, BandSet) and self.bands == other.bands

    def __hash__(self):
        return hash(self.bands)

    def __repr__(self):
        return "BandSet(" + ", ".join(f"{b.name}[{b.E_lo:g},{b.E_hi:g})" for b in self.bands) + ")"

    @property
    def names(self):
        return tuple(b.name for b in self.bands)

    @property
    def species(self):
        return tuple(b.species for b in self.bands)

    def index(self, name):
        for i, b in enumerate(self.bands):
            if b.name == name:
                return i
        raise KeyError(f"no band {name!r} in {self.names}")

    def band_at(self, E):
        """The band containing photon energy E, or None"""
        return next((b for b in self.bands if b.contains(E)), None)

    def gaps(self, E_lo, E_hi):
        """[(a, b)] sub-intervals of [E_lo, E_hi] that no band covers"""
        out, x = [], E_lo
        for b in sorted(self.bands, key=lambda b: b.E_lo):
            if b.E_hi <= x or b.E_lo >= E_hi:
                continue
            if b.E_lo > x:
                out.append((x, b.E_lo))
            x = max(x, b.E_hi)
        if x < E_hi:
            out.append((x, E_hi))
        return out

    def replace(self, name, *new):
        """The set with band name replaced by the bands new, in its place"""
        i = self.index(name)
        return BandSet(self.bands[:i] + tuple(new) + self.bands[i + 1:])

    def split(self, name, edges, names=None):
        """The set with band name cut at the interior energies edges; the parts keep its unit, shape and slope and
        are named names (default <name>_0, <name>_1, ...)"""
        b = self[name]
        cuts = [b.E_lo, *sorted(edges), b.E_hi]
        if any(not lo < hi for lo, hi in zip(cuts, cuts[1:])):
            raise ValueError(f"split edges {edges} are not inside band {name} [{b.E_lo}, {b.E_hi}]")
        names = names or [f"{name}_{i}" for i in range(len(cuts) - 1)]
        if len(names) != len(cuts) - 1:
            raise ValueError(f"{len(cuts) - 1} parts need as many names")
        return self.replace(name, *(replace(b, name=n, E_lo=lo, E_hi=hi) for n, lo, hi in zip(names, cuts, cuts[1:])))

    def declarations(self):
        """The radiation Species of the bands"""
        return [b.declaration() for b in self.bands]
