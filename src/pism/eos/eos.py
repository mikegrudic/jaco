"""Should eventually contain routines to generate the symbolic expression for internal energy and pressure, given
the abundances available in a network.
"""

import sympy as sp
from astropy.constants import k_B, m_p
from astropy import units as u
from pism.symbols import n_, T

boltzmann_cgs = k_B.to(u.erg / u.K).value
protonmass_cgs = m_p.to(u.g).value
# write down internal energy density in terms of number densities
egy_density = boltzmann_cgs * T * (1.5 * (n_("e-") + n_("H") + n_("H+") + n_("He") + n_("He+") + n_("He++")))
rho = protonmass_cgs * (n_("H") + n_("H+") + 4 * (n_("He") + n_("He+") + n_("He++")))
u = sp.simplify(egy_density / rho)

# def u_from_specieslist(specieslist):
