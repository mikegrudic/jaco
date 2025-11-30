"""Atomic data"""

import periodictable

atomic_weights = {str(el): el.mass for el in periodictable.elements}
atomic_numbers = {str(el): el.number for el in periodictable.elements}
atoms = list(atomic_weights)
bibliography = "https://iupac.qmul.ac.uk/AtWt/"
