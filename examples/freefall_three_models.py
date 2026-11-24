"""Free-fall collapse of one parcel in three models, at several metallicities: nelson (Nelson & Langer 1999),
starforge, and popsicle (Menon et al. 2026's KROME network), with jaco.onezone.

    python examples/freefall_three_models.py -o freefall_three_models.png --save tracks.npz
    python examples/freefall_three_models.py --Z 1 starforge popsicle        # one metallicity, two of the models

The same collapse for all three (d rho/dt = rho / t_ff from n_H = 0.1 cm^-3 and 300 K), the same external FUV field
(chi = 1 Draine = 1.69 Habing), shielding columns of half a Jeans length, metallicity Z in solar units (default 1,
1e-2, 1e-4 and 0); each model keeps its own abundances and host contract:

- starforge runs with the GIZMO host closures (the default host); its element totals and dust scale with Z.
- popsicle runs with its own closures, start state and per-Z parameters (host="plain").
- nelson's host inputs (A_V, tau_uv, the CO and H2 columns, the velocity gradient) have no closures in jaco, so they
  are written here from the same Jeans-length column. Its G_0 is the incident field, so it also runs with
  host="plain". The model assumes molecular hydrogen and has no optically thick physics: outside n_H ~ 1e2-1e5 cm^-3
  its track is an extrapolation. Its C, O and metal abundances and its dust extinction are scaled with Z here, but its
  thermal fits (the molecular-line cooling, the gas-dust exchange, the photoelectric heating) carry no metallicity:
  below solar its track is not a low-metallicity model, and it is drawn dashed.

One process per model and metallicity: the systems take minutes to build (starforge and popsicle) and seconds to
minutes to integrate.
"""
import argparse
import os
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np

N0, T0 = 0.1, 300.0
LABELS = {"nelson": "Nelson & Langer 1999", "starforge": "starforge", "popsicle": "popsicle (Menon+26)"}
COLORS = {"nelson": "#1b6e8a", "starforge": "#c0392b", "popsicle": "#444444"}


def run(name, Z, n_stop):
    """the free-fall track of one model at metallicity Z (solar units), as a dict of arrays; None if it fails"""
    warnings.simplefilter("ignore")
    from jaco.onezone import OneZone
    t0 = time.time()
    if name == "starforge":
        from jaco.models.starforge import make_model
        zone, x0 = OneZone(make_model()), None
        zone.params = zone.metallicity(Z)
    elif name == "popsicle":
        from jaco.models.popsicle import make_model
        from jaco.models.popsicle.closures import PopsicleClosures, parameters, INITIAL_ABUNDANCES
        zone = OneZone(make_model(), params=parameters(Z), closures=PopsicleClosures.DEFAULT, host="plain")
        x0 = dict(INITIAL_ABUNDANCES)
    elif name == "nelson":
        from jaco.models.nelson import make_model, PAPER_PARAMETERS

        def N_H(st):  # half a Jeans length of H nuclei, the shielding column of the other two
            return st["n_Htot"] * st["lambda_J"] / 2

        closures = {
            "A_V": lambda st: Z * N_H(st) / 1.87e21,
            "tau_uv": lambda st: Z * 1.25e-21 * N_H(st),  # NL97's UV optical depth per column
            "N_H2": lambda st: 0.5 * N_H(st),  # the model's hydrogen is molecular
            "N_CO": lambda st: max(st["x"]["CO"] * N_H(st), 1.0),
            "dv_dr": lambda st: abs(st["dlnrho_dt"]) / 3,  # homologous collapse
        }
        totals = {k: Z * PAPER_PARAMETERS[k] for k in ("x_C,tot", "x_O,tot", "x_Mg,tot")}
        # these closures depend on the state, which the analytic Jacobian does not chain: difference it instead
        zone = OneZone(make_model(), params={**PAPER_PARAMETERS, **totals, "G_0": 1.69}, closures=closures,
                       host="plain", analytic_jacobian=False)
        # NL99's start: carbon as C+, oxygen atomic, the metals ionized
        x0 = {"C+": totals["x_C,tot"], "Mg+": totals["x_Mg,tot"], "He+": 1e-9, "H3+": 1e-9,
              "CH": 0.0, "OH": 0.0, "CO": 0.0, "HCO+": 0.0}
    else:
        raise ValueError(f"unknown model {name}")
    t1 = time.time()
    try:
        tr = zone.free_fall(n0=N0, T0=T0, n_stop=n_stop, x0=x0, n_points=600)
    except Exception as e:  # one failed track should not cost the others
        print(f"{name} Z = {Z:g}: FAILED after {time.time() - t0:.0f} s: {type(e).__name__}: {e}", flush=True)
        return None
    print(f"{name} Z = {Z:g}: built in {t1 - t0:.0f} s, integrated in {time.time() - t1:.0f} s; minimum "
          f"{tr['T'].min():.3g} K at n_H = {tr['n_Htot'][np.argmin(tr['T'])]:.3g} cm^-3", flush=True)
    return {k: v for k, v in tr.items() if isinstance(v, np.ndarray)}


def z_label(Z):
    if Z in (0, 1):
        return f"$Z = {Z:g}$"
    e = np.log10(Z)
    return rf"$Z = 10^{{{e:g}}}$" if e == int(e) else f"$Z = {Z:g}$"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("models", nargs="*", default=list(LABELS), help="any of " + ", ".join(LABELS))
    ap.add_argument("--Z", default="1,1e-2,1e-4,0", help="metallicities in solar units, comma-separated")
    ap.add_argument("--n-stop", type=float, default=1e15, help="final n_H [cm^-3]")
    ap.add_argument("--jobs", type=int, default=None, help="processes (default: one per track, up to the CPU count)")
    ap.add_argument("-o", "--output", default="freefall_three_models.png", help="figure path")
    ap.add_argument("--save", default=None, help="write the tracks to this .npz")
    ap.add_argument("--load", default=None, help="plot the tracks of a saved .npz instead of running them")
    args = ap.parse_args()
    Zs = [float(z) for z in args.Z.split(",")]
    jobs = [(m, Z) for Z in Zs for m in args.models]

    if args.load:
        saved = np.load(args.load)
        tracks = {(m, Z): {k.split("/", 2)[2]: saved[k] for k in saved.files if k.startswith(f"{m}/Z={Z:g}/")}
                  for m, Z in jobs}
        tracks = {job: tr for job, tr in tracks.items() if tr}
    else:
        with ProcessPoolExecutor(args.jobs or min(len(jobs), os.cpu_count() or 1)) as pool:
            results = pool.map(run, *zip(*[(m, Z, args.n_stop) for m, Z in jobs]))
            tracks = {job: tr for job, tr in zip(jobs, results) if tr is not None}
    if args.save:
        np.savez(args.save, **{f"{m}/Z={Z:g}/{k}": v for (m, Z), tr in tracks.items() for k, v in tr.items()})

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    amp = r"\&" if plt.rcParams["text.usetex"] else "&"
    ncols = 2 if len(Zs) > 1 else 1
    nrows = -(-len(Zs) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols + 0.6, 3.1 * nrows + 0.5), sharex=True, sharey=True,
                             squeeze=False)
    note = None  # the first dashed line, for its legend entry
    for ax, Z in zip(axes.flat, Zs):
        for m in args.models:
            if (m, Z) not in tracks:
                continue
            tr = tracks[(m, Z)]
            unscaled = m == "nelson" and Z != 1  # its thermal fits carry no metallicity
            line, = ax.loglog(tr["n_Htot"], tr["T"], color=COLORS[m], lw=1.6, ls="--" if unscaled else "-",
                              label=LABELS[m].replace("&", amp))
            if unscaled and note is None:
                note = (ax, line)
        ax.text(0.5, 0.93, z_label(Z), transform=ax.transAxes, ha="center", va="top")
        ax.set_xlim(N0, args.n_stop)
        ax.tick_params(which="both", direction="in", top=True, right=True)
    for ax in axes.flat[len(Zs):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel(r"$n_{\rm H}$ [cm$^{-3}$]")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$T$ [K]")
    models_legend = axes.flat[0].legend(frameon=False, fontsize=8, loc="center")
    if note is not None:
        ax, line = note
        if ax is axes.flat[0]:
            ax.add_artist(models_legend)  # a second legend on the same axes replaces the first otherwise
        ax.legend([line], [f"Nelson {amp} Langer: thermal fits not scaled with $Z$"], frameon=False, fontsize=8,
                  loc="lower left")
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
