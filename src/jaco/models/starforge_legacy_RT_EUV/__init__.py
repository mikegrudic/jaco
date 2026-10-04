"""GIZMO's legacy cooling module with the ionizing band alone (M1 RADTRANSFER with RT_CHEM_PHOTOION and no dust bands,
e.g. HII_region_simple and the Iliev test), as one jaco model whose unknowns include the band's photons.

starforge_legacy plus, from jaco.models.starforge.radiation: photoionization of H, each event taking one photon of the
band at c_tilde and heating the gas by eps_HI (no donation: without RT_OPTICAL_NIR GIZMO's kick loses the absorbed
energy), and Compton heating off the band. GIZMO keeps transport only; the host hands the model the band's
post-transport photons as the step's initial value and writes the solved ones back. Without RT_INFRARED GIZMO's dust
temperature is its non-RT estimate (the parameter Td), and the emission corrections take the CMB.
"""

from jaco.declarations import Species, Output
from ..starforge.starforge import GIZMO_FAMILY  # noqa: F401  (the host-code family)
from ..starforge import radiation as rt
from ..starforge.radiation import EUV
from ..starforge_legacy import make_model as make_legacy, GIZMO_CLUMPING
from ..starforge_legacy_RT import GIZMO_DENSITY, GIZMO_ABUNDANCES, BAND_SPECIES, RT_PARAMETERS


def make_model():
    """The STARFORGE_LEGACY_RT_EUV model"""
    base = make_legacy()
    return base.evolve(
        list(base.processes.values()) + [rt.photoionization(), rt.compton_off_bands([EUV])],
        solve_vars=base.solve_vars + (EUV,),
        time_dependent=("T", "H+", "H_2", EUV),
        rules=[GIZMO_CLUMPING, GIZMO_DENSITY, GIZMO_ABUNDANCES],
        parameters=[*base.parameters, *RT_PARAMETERS],
        species=[*base.species, BAND_SPECIES[EUV]],
        outputs=[Output("photoionization_rate", heat_of={(p, EUV): -1 / rt.rsol for p in rt.EUV_SINKS},
                        units="cm^-3 s^-1", doc="photoionizations by the ionizing band per unit volume (each takes one "
                                                "photon)")],
    )
