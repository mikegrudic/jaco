"""Free-fall collapse of one parcel in three models: nelson (Nelson & Langer 1999), starforge, and popsicle (Menon et
al. 2025's KROME network), with jaco.onezone.

    python examples/freefall_three_models.py -o freefall_three_models.png --save tracks.npz

The same collapse for all three (d rho/dt = rho / t_ff from n_H = 0.1 cm^-3 and 300 K), the same external FUV field
(chi = 1 Draine = 1.69 Habing), solar metallicity, shielding columns of half a Jeans length; each model keeps its own
abundances and host contract:

- starforge runs with the GIZMO host closures (the default host).
- popsicle runs with its own closures and start state (host="plain").
- nelson's host inputs (A_V, tau_uv, the CO and H2 columns, the velocity gradient) have no closures in jaco, so they
  are written here from the same Jeans-length column. Its G_0 is the incident field, so it also runs with
  host="plain". The model assumes molecular hydrogen and has no optically thick physics: outside n_H ~ 1e2-1e5 cm^-3
  its track is an extrapolation.

One process per model: the systems take minutes to build (starforge and popsicle) and seconds to integrate.
"""
import argparse
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np

N0, T0 = 0.1, 300.0
LABELS = {"nelson": "Nelson & Langer 1999", "starforge": "starforge", "popsicle": "popsicle (Menon+25)"}
COLORS = {"nelson": "#1b6e8a", "starforge": "#c0392b", "popsicle": "#444444"}


def run(name, n_stop):
    """the free-fall track of one model, as a dict of arrays"""
    warnings.simplefilter("ignore")
    from jaco.onezone import OneZone
    t0 = time.time()
    if name == "starforge":
        from jaco.models.starforge import make_model
        zone, x0 = OneZone(make_model()), None
    elif name == "popsicle":
        from jaco.models.popsicle import make_model
        from jaco.models.popsicle.closures import PopsicleClosures, parameters, INITIAL_ABUNDANCES
        zone = OneZone(make_model(), params=parameters(1.0), closures=PopsicleClosures.DEFAULT, host="plain")
        x0 = dict(INITIAL_ABUNDANCES)
    elif name == "nelson":
        from jaco.models.nelson import make_model, PAPER_PARAMETERS

        def N_H(st):  # half a Jeans length of H nuclei, the shielding column of the other two
            return st["n_Htot"] * st["lambda_J"] / 2

        closures = {
            "A_V": lambda st: N_H(st) / 1.87e21,
            "tau_uv": lambda st: 1.25e-21 * N_H(st),  # NL97's UV optical depth per column
            "N_H2": lambda st: 0.5 * N_H(st),  # the model's hydrogen is molecular
            "N_CO": lambda st: max(st["x"]["CO"] * N_H(st), 1.0),
            "dv_dr": lambda st: abs(st["dlnrho_dt"]) / 3,  # homologous collapse
        }
        # these closures depend on the state, which the analytic Jacobian does not chain: difference it instead
        zone = OneZone(make_model(), params={**PAPER_PARAMETERS, "G_0": 1.69}, closures=closures, host="plain",
                       analytic_jacobian=False)
        # NL99's start: carbon as C+, oxygen atomic, the metals ionized
        x0 = {"C+": PAPER_PARAMETERS["x_C,tot"], "Mg+": PAPER_PARAMETERS["x_Mg,tot"], "He+": 1e-9, "H3+": 1e-9,
              "CH": 0.0, "OH": 0.0, "CO": 0.0, "HCO+": 0.0}
    else:
        raise ValueError(f"unknown model {name}")
    t1 = time.time()
    tr = zone.free_fall(n0=N0, T0=T0, n_stop=n_stop, x0=x0, n_points=600)
    print(f"{name}: built in {t1 - t0:.0f} s, integrated in {time.time() - t1:.0f} s; minimum {tr['T'].min():.3g} K "
          f"at n_H = {tr['n_Htot'][np.argmin(tr['T'])]:.3g} cm^-3", flush=True)
    return {k: v for k, v in tr.items() if isinstance(v, np.ndarray)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("models", nargs="*", default=list(LABELS), help="any of " + ", ".join(LABELS))
    ap.add_argument("--n-stop", type=float, default=1e15, help="final n_H [cm^-3]")
    ap.add_argument("-o", "--output", default="freefall_three_models.png", help="figure path")
    ap.add_argument("--save", default=None, help="write the tracks to this .npz")
    args = ap.parse_args()

    with ProcessPoolExecutor(len(args.models)) as pool:
        tracks = dict(zip(args.models, pool.map(run, args.models, [args.n_stop] * len(args.models))))
    if args.save:
        np.savez(args.save, **{f"{m}/{k}": v for m, tr in tracks.items() for k, v in tr.items()})

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    amp = r"\&" if plt.rcParams["text.usetex"] else "&"
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    for m, tr in tracks.items():
        ax.loglog(tr["n_Htot"], tr["T"], color=COLORS[m], lw=1.6, label=LABELS[m].replace("&", amp))
    ax.set_xlabel(r"$n_{\rm H}$ [cm$^{-3}$]")
    ax.set_ylabel(r"$T$ [K]")
    ax.set_xlim(N0, args.n_stop)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
