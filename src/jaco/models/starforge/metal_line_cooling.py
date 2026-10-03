"""Metal line cooling for the STARFORGE model.

Loads tabulated cooling rate coefficients from spcool_tables.hdf5 and builds
cooling processes of the form:

    rate per volume = k(n_Htot, T) * n_e * n_X,tot * C_2

where n_X,tot is the total number density of element X summed over every
network species carrying it (the tables already sum over ionization stages).
The table values are normalized such that multiplying by n_e * n_X,tot yields
the cooling rate density; the 3D (z, n_H, T) table is sliced to a 2D
(n_Htot, T) grid at fixed redshift. As in GIZMO, the tables are tapered below
their 100 K edge instead of being held at the edge value, and are applied at
all T: negative entries (UV-background photoheating) heat the gas, and the
CMB-bath factor multiplies the summed rate only where it is net cooling.

The HDF5 file lives alongside this module but is not tracked in git (8.7 MB); build it once with
``python -m jaco.models.starforge.convert_spcool_tables``.
"""

from importlib.resources import files
import numpy as np
import h5py

import sympy as sp
from jaco.processes import ThermalTerm
from jaco.symbols import T, table_interp_2d, n_, x_
from .symbols import n_Htot, log_T, cmb_bath_factor


def _hdf5_path():
    """Locate spcool_tables.hdf5 in the package directory."""
    path = files("jaco.models.starforge").joinpath("spcool_tables.hdf5")
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} not found. Generate it with: python -m jaco.models.starforge.convert_spcool_tables"
        )
    return path

# (HDF5 dataset name, species symbol used in the chemical network)
METAL_SPECIES = [
    ("Carbon_cooling", "C"),
    ("Nitrogen_cooling", "N"),
    ("Oxygen_cooling", "O"),
    ("Neon_cooling", "Ne"),
    ("Magnesium_cooling", "Mg"),
    ("Silicon_cooling", "Si"),
    ("Sulfur_cooling", "S"),
    ("Calcium_cooling", "Ca"),
    ("Iron_cooling", "Fe"),
]


# network species carrying each element, with multiplicity; elements not listed are carried by their atom alone
ELEMENT_CARRIERS = {"C": {"C": 1, "C+": 1, "CO": 1}, "O": {"O": 1, "CO": 1}}

# GIZMO: LambdaMetal *= exp(-min((2 - log10 T)^2 / 0.1, 40)) below the tables' 100 K edge
low_T_taper = sp.Piecewise((sp.exp(-sp.Min((2 - log_T) ** 2 / 0.1, 40)), T < 100), (1, True))


def element_number_density(element):
    """Total number density of an element summed over the network species that carry it.

    Written in abundances so that atom conservation reduces the sum to the total-abundance parameter.
    """
    return n_Htot * sum(count * x_(s) for s, count in ELEMENT_CARRIERS.get(element, {element: 1}).items())


def _load_redshift_slice(dataset_name, z=0.0, hdf5_path=None):
    """Read a 2D (n_Htot, T) slice of the cooling table at the closest redshift to z.

    Returns
    -------
    table_2d : np.ndarray, shape (n_nH, n_T)
        Rate coefficient values (erg cm^3 / s).
    nH_axis, T_axis : np.ndarray
        Linear-space axis values (the table is uniform in log space).
    """
    if hdf5_path is None:
        hdf5_path = _hdf5_path()
    with h5py.File(str(hdf5_path), "r") as f:
        log_T = f["log_T"][:]
        log_nH = f["log_nH"][:]
        log_zp1 = f["log_one_plus_z"][:]
        iz = int(np.argmin(np.abs(log_zp1 - np.log10(1.0 + z))))
        table_2d = f[dataset_name][iz].astype(np.float64)
    return table_2d, 10.0**log_nH, 10.0**log_T


def metal_line_cooling_rate(species, dataset, z=0.0):
    """Tabulated volumetric cooling rate of one element (erg cm^-3 s^-1; negative = net heating), without clumping
    or CMB-bath factors: `k(n_Htot, T) * n_e * n_X,tot`, tapered below 100 K.

    Parameters
    ----------
    species : str
        Element symbol used in the chemical network (e.g. "C", "O", "Fe").
    dataset : str
        HDF5 dataset name in spcool_tables.hdf5 (e.g. "Carbon_cooling").
    z : float
        Redshift slice to use.
    """
    table_2d, nH_axis, T_axis = _load_redshift_slice(dataset, z=z)
    table_name = f"{dataset}_z{int(round(z))}"
    k = table_interp_2d(
        table_name, table_2d, [nH_axis, T_axis], n_Htot, T,
        log_axes=[True, True],
    )
    return k * low_T_taper * n_("e-") * element_number_density(species)


def cmb_corrected(Lambda):
    """GIZMO applies the CMB-bath factor to net metal-line cooling only, not to net photoheating."""
    return sp.Max(Lambda, 0) * cmb_bath_factor + sp.Min(Lambda, 0)


def metal_line_cooling_process(species, dataset, z=0.0):
    """Cooling process for a single tabulated element (see metal_line_cooling_rate)."""
    return ThermalTerm(
        -cmb_corrected(metal_line_cooling_rate(species, dataset, z)),
        name=f"{species} line cooling",
        bibliography=["2009MNRAS.393...99W"],  # Wiersma+ tables used in GIZMO
        clumping=sp.Symbol("C_2"),
    )


def metal_line_cooling(z=0.0):
    """Combined metal-line cooling of all tabulated elements; the CMB-bath factor acts on the sum, as in GIZMO.

    Scaled by the switch parameter f_metal: GIZMO applies these tables only when its UV background is loaded (J_UV != 0).
    """
    Lambda = sum(metal_line_cooling_rate(sp_, ds, z=z) for ds, sp_ in METAL_SPECIES)
    return ThermalTerm(
        -sp.Symbol("f_metal") * cmb_corrected(Lambda),
        name="Metal line cooling",
        bibliography=["2009MNRAS.393...99W"],  # Wiersma+ tables used in GIZMO
        clumping=sp.Symbol("C_2"),
    )
