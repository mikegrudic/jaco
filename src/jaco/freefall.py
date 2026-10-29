"""One-zone free-fall collapse of a model: the density-temperature track of a parcel whose density grows at the
free-fall rate while its thermochemistry evolves (Omukai 2000; Glover & Jappsen 2007; Menon et al. 2025, Fig. 2).

    track = free_fall(make_model(), n0=0.1, T0=300, n_stop=1e15)
    plt.loglog(track["n_Htot"], track["T"])

or, from the command line, python -m jaco.freefall starforge starforge_legacy -o track.png.

The parcel: d rho/dt = rho / t_ff with t_ff = sqrt(3 pi / (32 G rho)), compressional heating pdv_work = P / t_ff per
unit volume, velocity gradient 1/t_ff (trace-free part zero: homologous compression does not clump), and the
shielding columns of half a Jeans length, N_H = n_Htot lambda_J / 2 with lambda_J = sqrt(pi c_s^2 / (G rho)).
Everything else the model reads (radiation field, cosmic rays, metal abundances) is the model's default or the
caller's `params`; `closures` replaces any of the built-in prescriptions by a function of the state.

The integration variable is s = ln(n_Htot / n0): the thermochemistry's own rates set the stiffness, which the
integrator (scipy's LSODA by default) handles; the exponential density growth does not add to it.
"""

import argparse
import sys

import numpy as np
from scipy.integrate import solve_ivp

from .ode import OdeSystem

__all__ = ["free_fall", "FreeFallClosures"]

G_CGS = 6.674e-8
PROTONMASS = 1.6726e-24


def t_freefall(rho):
    return np.sqrt(3 * np.pi / (32 * G_CGS * rho))


class FreeFallClosures:
    """The built-in prescriptions for the parameters a host would supply, as functions of the state dict the
    integrator passes (n_Htot, T, rho, P, x, t_ff, lambda_J)."""

    @staticmethod
    def pdv_work(st):
        return st["P"] / st["t_ff"]

    @staticmethod
    def rho(st):
        return st["rho"]

    @staticmethod
    def grad_v(st):
        return 1.0 / st["t_ff"]

    @staticmethod
    def grad_v_tf(st):
        return 0.0

    @staticmethod
    def N_H(st):
        return st["n_Htot"] * st["lambda_J"] / 2

    N_H_cool = N_H

    @staticmethod
    def cell_size(st):
        return st["lambda_J"]

    DEFAULT = {"pdv_work": pdv_work.__func__, "rho": rho.__func__, "∇v": grad_v.__func__, "∇v_tf": grad_v_tf.__func__,
               "N_H": N_H.__func__, "N_H_cool": N_H.__func__, "Δx": cell_size.__func__}


def free_fall(model, n0=0.1, T0=300.0, x0=None, n_stop=1e15, params=None, closures=None, n_points=400,
              method="LSODA", rtol=1e-5, atol=1e-12, system=None, verbose=False):
    """Integrate a model's thermochemistry along a free-fall density history.

    Parameters
    ----------
    model: jaco.model.Model
    n0, T0: float
        Initial H nucleus density [cm^-3] and temperature [K].
    x0: dict, optional
        Initial abundances per H nucleus of the solve species; the rest start at their floors (atomic, neutral gas).
    n_stop: float
        Final density [cm^-3].
    params: dict, optional
        Parameter values by expression name (e.g. {"G_0": 0, "ISRF": 0}) for the whole track.
    closures: dict, optional
        Parameter name -> function of the state dict, replacing or adding to FreeFallClosures.DEFAULT (only the
        parameters the model reads are used).
    n_points: int
        Output points, uniform in ln n.
    method, rtol, atol: solve_ivp options. atol applies to T in K and to the abundances per H.
    system: OdeSystem, optional
        A prepared system for this model (its construction is the slow part: a minute for starforge).
    verbose: bool
        Report progress every decade of density.

    Returns
    -------
    dict of arrays: t [s], n_Htot, T, rho, P, u, the abundances by species name, the balance variables (e.g. Td),
    heat [erg cm^-3 s^-1] and pdv_work, plus "system" (the OdeSystem) and "sol" (the solve_ivp result).
    """
    sys_ = system or OdeSystem(model, params)
    sys_.reset()
    P_fixed = dict(params or {})
    closures = {**FreeFallClosures.DEFAULT, **(closures or {})}
    active = {k: f for k, f in closures.items() if k in sys_.parameters}
    y0 = sys_.state(T0, x0 or {})
    s_end = np.log(n_stop / n0)
    state = {"n_Htot": n0}
    last_decade = [None]

    def params_at(n, T, x):
        p = dict(P_fixed)
        p["n_Htot"] = n
        # a first pass for rho and P at this state (the closures need them); then the closures
        base = sys_.rates(T, x, {**p, **{k: 0.0 for k in active if k in ("pdv_work",)}}) if active else None
        if base is not None:
            st = {"n_Htot": n, "T": T, "x": x, "rho": base["rho"], "P": base["P"], "u": base["u"]}
            st["t_ff"] = t_freefall(base["rho"])
            st["lambda_J"] = np.sqrt(np.pi * base["P"] / base["rho"] / (G_CGS * base["rho"]))
            for k, f in active.items():
                p[k] = f(st)
        return p

    def rhs(s, y):
        n = n0 * np.exp(s)
        T, x = y[0], {sp_: v for sp_, v in zip(sys_.species, y[1:])}
        p = params_at(n, T, x)
        r = sys_.rates(T, x, p)
        tff = t_freefall(r["rho"])
        if verbose:
            d = int(np.floor(np.log10(n)))
            if last_decade[0] is None or d > last_decade[0]:  # the integrator's trial steps cross decades both ways
                last_decade[0] = d
                print(f"  n = 1e{d} cm^-3: T = {T:.4g} K", file=sys.stderr)
        return tff * np.array([r["dTdt"]] + [r["dxdt"][sp_] for sp_ in sys_.species])

    s_eval = np.linspace(0, s_end, n_points)
    sol = solve_ivp(rhs, (0, s_end), y0, method=method, t_eval=s_eval, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"free_fall: {sol.message}")
    n = n0 * np.exp(sol.t)
    out = {"n_Htot": n, "T": sol.y[0]}
    for i, sp_ in enumerate(sys_.species):
        out[sp_] = sol.y[1 + i]
    # the derived quantities and the time along the track
    cols = {k: [] for k in ["rho", "P", "u", "heat", "pdv_work"] + sys_.balances}
    sys_.reset()
    for k in range(len(n)):
        x = {sp_: out[sp_][k] for sp_ in sys_.species}
        p = params_at(n[k], out["T"][k], x)
        r = sys_.rates(out["T"][k], x, p)
        for c in cols:
            cols[c].append(r[c] if c in r else p.get(c, np.nan))
    for c, v in cols.items():
        out[c] = np.array(v)
    tff = t_freefall(out["rho"])
    out["t"] = np.concatenate([[0], np.cumsum(0.5 * (tff[1:] + tff[:-1]) * np.diff(sol.t))])
    out["system"], out["sol"] = sys_, sol
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="density-temperature track of a free-fall collapse for jaco models")
    ap.add_argument("models", nargs="+", help="model names (jaco.models.<name>)")
    ap.add_argument("--n0", type=float, default=0.1)
    ap.add_argument("--T0", type=float, default=300.0)
    ap.add_argument("--n-stop", type=float, default=1e15)
    ap.add_argument("--param", action="append", default=[], metavar="NAME=VALUE", help="e.g. G_0=0, ISRF=0, Z_d=0.01")
    ap.add_argument("--rtol", type=float, default=1e-5)
    ap.add_argument("-o", "--output", default="free_fall_track.png")
    ap.add_argument("--save", default=None, help="write the tracks to this .npz")
    args = ap.parse_args(argv)
    params = {}
    for kv in args.param:
        k, v = kv.split("=", 1)
        params[k] = float(v)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from importlib import import_module
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    tracks = {}
    for name in args.models:
        model = import_module(f"jaco.models.{name}").make_model()
        print(f"{name}: building the system", file=sys.stderr)
        tr = free_fall(model, n0=args.n0, T0=args.T0, n_stop=args.n_stop, params=params, rtol=args.rtol, verbose=True)
        tracks[name] = {k: v for k, v in tr.items() if isinstance(v, np.ndarray)}
        ax.loglog(tr["n_Htot"], tr["T"], lw=1.6, label=name)
        print(f"{name}: T = {tr['T'][-1]:.4g} K at n = {tr['n_Htot'][-1]:.3g}; {tr['sol'].nfev} evaluations", file=sys.stderr)
    ax.set_xlabel(r"$n_{\rm H}$ [cm$^{-3}$]")
    ax.set_ylabel("T [K]")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(args.output)
    if args.save:
        np.savez(args.save, **{f"{m}/{k}": v for m, d in tracks.items() for k, v in d.items()})


if __name__ == "__main__":
    main()
