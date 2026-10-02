"""Convert GIZMO's per-species metal-line cooling tables into spcool_tables.hdf5.

GIZMO's COOL_METAL_LINES_BY_SPECIES tables (spcool_0 ... spcool_48, from Wiersma, Schaye & Smith 2009,
distributed at https://users.flatironinstitute.org/~mgrudic/gizmo_tests/spcool_tables.tgz) are raw
native-endian float32 arrays of shape (10, 41, 176): index 0 is the tabulated n_e/n_H, indices 1-9 are
Lambda_X/n_H^2 (erg cm^3/s) for C, N, O, Ne, Mg, Si, S, Ca, Fe at the table's solar abundances, on
log10(n_H) = -8..0 and log10(T) = 2..9. File i holds log10(1+z) = i/48.

The output stores, per element, the 2-body coefficient

    k_X(n_H, T) = Lambda_X / ((n_e/n_H)_table * x_X,table)

so that the cooling rate density is k_X * n_e * n_X for arbitrary electron and element abundances.
x_X,table are the Wiersma+09 (Table 1) solar number abundances that the tables assume.

Usage:
    python -m jaco.models.starforge.convert_spcool_tables [/path/to/spcool_tables] [-o out.hdf5]

With no directory argument the tables are downloaded from SPCOOL_URL.
"""

import argparse
import os
import tarfile
import tempfile
import urllib.request

import h5py
import numpy as np

SPCOOL_URL = "https://users.flatironinstitute.org/~mgrudic/gizmo_tests/spcool_tables.tgz"
N_NH, N_T, N_Z = 41, 176, 49
LOG_NH = np.linspace(-8.0, 0.0, N_NH)
LOG_T = np.linspace(2.0, 9.0, N_T)
LOG_ONE_PLUS_Z = np.arange(N_Z) / 48.0

# (dataset name, Wiersma+09 Table 1 solar abundance n_X/n_H), in GIZMO table order (species index 1-9)
SPECIES = [
    ("Carbon_cooling", 2.46e-4),
    ("Nitrogen_cooling", 8.51e-5),
    ("Oxygen_cooling", 4.90e-4),
    ("Neon_cooling", 1.00e-4),
    ("Magnesium_cooling", 3.47e-5),
    ("Silicon_cooling", 3.47e-5),
    ("Sulfur_cooling", 1.86e-5),
    ("Calcium_cooling", 2.29e-6),
    ("Iron_cooling", 2.82e-5),
]


def read_gizmo_table(path):
    raw = np.fromfile(path, dtype=np.float32)
    expected = (1 + len(SPECIES)) * N_NH * N_T
    if raw.size != expected:
        raise ValueError(f"{path}: {raw.size} floats, expected {expected}")
    return raw.reshape(1 + len(SPECIES), N_NH, N_T).astype(np.float64)


def convert(spcool_dir, output):
    tables = np.array([read_gizmo_table(os.path.join(spcool_dir, f"spcool_{i}")) for i in range(N_Z)])
    ne_over_nH = tables[:, 0]
    if np.any(ne_over_nH <= 0):
        raise ValueError("non-positive tabulated n_e/n_H")
    with h5py.File(output, "w") as f:
        f.attrs["source"] = "GIZMO COOL_METAL_LINES_BY_SPECIES tables (Wiersma, Schaye & Smith 2009)"
        f.attrs["normalization"] = "Lambda_X / ((n_e/n_H)_table * x_X,table) [erg cm^3 / s]; rate density = k * n_e * n_X"
        f.create_dataset("log_nH", data=LOG_NH)
        f.create_dataset("log_T", data=LOG_T)
        f.create_dataset("log_one_plus_z", data=LOG_ONE_PLUS_Z)
        for k, (name, x_solar) in enumerate(SPECIES, start=1):
            coeff = tables[:, k] / (ne_over_nH * x_solar)
            dset = f.create_dataset(name, data=coeff.astype(np.float32), compression="gzip",
                                    compression_opts=9, shuffle=True)
            dset.attrs["x_solar"] = x_solar


def main():
    default_out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "spcool_tables.hdf5")
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("spcool_dir", nargs="?", help="directory containing spcool_0 ... spcool_48")
    parser.add_argument("-o", "--output", default=default_out, help=f"output file (default: {default_out})")
    args = parser.parse_args()
    if args.spcool_dir:
        convert(args.spcool_dir, args.output)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            tgz = os.path.join(tmp, "spcool_tables.tgz")
            print(f"downloading {SPCOOL_URL}")
            urllib.request.urlretrieve(SPCOOL_URL, tgz)
            with tarfile.open(tgz) as tar:
                tar.extractall(tmp, filter="data")
            convert(os.path.join(tmp, "spcool_tables"), args.output)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
