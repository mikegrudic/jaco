"""Outputs for an energy ledger of the H2 chemistry: what the collisional dissociations took from the gas and how many
there were, the same for the formation channels that give the gas heat, the energy the formed H2 radiated (capped), the
net H2 production of every process, the whole chemical heat of the H2/H- network, and the heat of every other process
but PdV work (the radiation the optically-thick limit caps, the formation radiation included, and the gas-dust
exchange, which passes the dust's capped radiation to the gas), with the limit's factor. With PdV work these heats are
every process, so their sum is the energy row's rho du/dt. The solve is unchanged; the host integrates these over each
step (rate * Delta_t under backward Euler)."""

from jaco.declarations import Output
from .h2_chemistry import h2_chemistry_processes
from .h2_chemistry.chemical_heat import RADIATION_OF
from .optically_thick import LIMIT_ID, F_THICK

DISSOCIATION = tuple(f"Collisional dissociation of H_2 by {c}" for c in ("H+", "e-", "H_2", "H", "He"))
FORMATION_CHANNELS = ("Formation of H_2 on dust grains", "3-body formation of H_2", "Associative detachment of H with H-")


def ledger_outputs(model):
    """The ledger's Outputs for a starforge model"""
    chemistry = [i for i in model.processes if i in {p.name for p in h2_chemistry_processes()}
                 or i.startswith("3-body formation of H_2")]
    formation = [i for i in model.processes if i in FORMATION_CHANNELS or i.startswith("3-body formation of H_2")]
    h2_rows = [i for i, p in model.processes.items() if "H_2" in p.network]
    radiation = [RADIATION_OF + i for i in formation if RADIATION_OF + i in model]
    outputs = [
        Output("h2_diss_heat", heat_of=DISSOCIATION, units="erg cm^-3 s^-1",
               doc="heat the collisional dissociations of H2 give the gas (negative)"),
        Output("h2_diss_events", heat_of={(i, "H_2"): -1 for i in DISSOCIATION}, units="cm^-3 s^-1",
               doc="collisional dissociations of H2"),
        Output("h2_form_heat", heat_of=formation, units="erg cm^-3 s^-1",
               doc="heat of H2 formation on grains, by three-body collisions and by associative detachment"),
        Output("h2_form_events", heat_of={(i, "H_2"): 1 for i in formation}, units="cm^-3 s^-1",
               doc="H2 formed on grains, by three-body collisions and by associative detachment"),
        Output("h2_form_radiated", heat_of=radiation, units="erg cm^-3 s^-1",
               doc="heat of the radiation of the H2 formed, capped like the other cooling (negative)"),
        Output("h2_net", heat_of={(i, "H_2"): 1 for i in h2_rows}, units="cm^-3 s^-1",
               doc="net production of H2 by every process, d n_H2 / dt"),
        Output("h2_chem_heat", heat_of=chemistry, units="erg cm^-3 s^-1",
               doc="heat of every H2/H- reaction"),
    ]
    if any(r.name == LIMIT_ID for r in model.rules):
        radiative = [i for i, p in model.processes.items() if i not in set(chemistry) | {"PdV work"} and "heat" in p.network]
        outputs += [
            Output("q_capped", heat_of=radiative, units="erg cm^-3 s^-1",
                   doc="heat of every process but the H2/H- reactions and PdV work: the capped radiation and the "
                       "gas-dust exchange"),
            Output("f_thick", F_THICK, doc="the optically-thick limit's factor on their heat"),
        ]
    return outputs


def with_ledger(model):
    return model.evolve(outputs=(*model.outputs, *ledger_outputs(model)))
