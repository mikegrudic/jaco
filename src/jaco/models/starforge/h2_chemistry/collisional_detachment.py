"""Implementation of collisional detachment for the H- ion"""

from ..symbols import sp, T
from jaco.processes import ChemicalReaction


def Hminus_collisional_detachment(collider):
    """Collisional dissociation of H- by colliding species"""

    lnTeV = sp.ln(T / 8.617e-5)
    match collider:
        case "H":
            k_lo = 2.5634e-9 * T**1.78186
            k_hi = sp.exp(
                -2.0372609e1
                + 1.13944933e0 * lnTeV
                - 1.4210135e-1 * lnTeV * lnTeV
                + 8.4644554e-3 * lnTeV * lnTeV * lnTeV
                - 1.4327641e-3 * lnTeV * lnTeV * lnTeV * lnTeV
                + 2.0122503e-4 * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV
                + 8.6639632e-5 * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV
                - 2.5850097e-5 * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV
                + 2.4555012e-6 * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV
                - 8.0683825e-8 * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV
            )
            k = sp.Piecewise((k_lo, T < 1160.45), (k_hi, T >= 1160.45))
        case "e-":
            k = sp.exp(
                -1.801849334e1
                + 2.36085220e0 * lnTeV
                - 2.827443e-1 * lnTeV**2
                + 1.62331664e-2 * lnTeV**3
                - 3.36501203e-2 * lnTeV**4
                + 1.17832978e-2 * lnTeV**5
                - 1.65619470e-3 * lnTeV**6
                + 1.06827520e-4 * lnTeV**7
                - 2.63128581e-6 * lnTeV**8
            )
        case _:
            raise NotImplementedError(f"Collisional dissociation rate of H- with {collider} not implemented.")

    return ChemicalReaction(
        f"H- + {collider} -> H + {collider} + e-",
        k,
        name=f"Collisional detachment of H- by {collider}",
        bibliography=["1987ephh.book.....J"],
    )


model_process = Hminus_collisional_detachment("H")
