"""One-zone models: a parcel of gas with the model's thermochemistry and a prescribed density history.

    zone = OneZone(make_model(), params={"G_0": 0, "ISRF": 0})
    ff = zone.free_fall(n0=0.1, T0=300, n_stop=1e15)             # collapse at the free-fall rate
    eq = zone.equilibrium(np.logspace(-2, 6, 33), T0=1e4)          # thermochemical equilibrium at each density
    ev = zone.evolve(n=100, T0=1e4, t_end=1e14)                    # time evolution at fixed density
    tr = zone.track(dlnrho_dt=lambda st: 1e-13, n0=1, T0=1e4, t_end=1e14)   # any compression law

or, from the command line, python -m jaco.onezone free_fall starforge starforge_legacy -o track.png and
python -m jaco.onezone equilibrium starforge --n-grid 1e-2,1e6,33.

Every mode integrates the same right-hand side (jaco.ode.OdeSystem): d ln rho/dt is the mode's, the compressional
heating is pdv_work = P d ln rho/dt per unit volume, the velocity gradient is |d ln rho/dt| with no trace-free part
(homologous compression does not clump), and the shielding columns are half a Jeans length, N_H = n_Htot lambda_J / 2
with lambda_J = sqrt(pi c_s^2 / (G rho)). Everything else the model reads (radiation field, cosmic rays, metal
abundances, the dust's ISRF absorption, the optically-thick opacity) is the model's default or the caller's `params`;
`closures` replaces any of the built-in prescriptions by a function of the state. The free-fall mode (Omukai 2000;
Glover & Jappsen 2007; Menon et al. 2025, Fig. 2) integrates in s = ln(n/n0), so the exponential density growth adds
no stiffness; the others integrate in time.
"""

import argparse
import sys

import numpy as np
from scipy.integrate import solve_ivp

from .ode import OdeSystem

__all__ = ["OneZone", "Closures", "t_freefall", "jeans_length"]

G_CGS = 6.674e-8


def t_freefall(rho):
    return np.sqrt(3 * np.pi / (32 * G_CGS * rho))


def jeans_length(P, rho):
    return np.sqrt(np.pi * P / rho / (G_CGS * rho))


class Closures:
    """The built-in prescriptions for the parameters a host would supply, as functions of the state dict
    (n_Htot, T, x, rho, P, u, dlnrho_dt, t_ff, lambda_J)."""

    @staticmethod
    def pdv_work(st):
        return st["P"] * st["dlnrho_dt"]

    @staticmethod
    def rho(st):
        return st["rho"]

    @staticmethod
    def grad_v(st):
        return abs(st["dlnrho_dt"])

    @staticmethod
    def grad_v_tf(st):
        return 0.0

    @staticmethod
    def N_H(st):
        return st["n_Htot"] * st["lambda_J"] / 2

    @staticmethod
    def cell_size(st):
        return st["lambda_J"]

    DEFAULT = {"pdv_work": pdv_work.__func__, "rho": rho.__func__, "∇v": grad_v.__func__, "∇v_tf": grad_v_tf.__func__,
               "N_H": N_H.__func__, "N_H_cool": N_H.__func__, "Δx": cell_size.__func__}


class OneZone:
    """A model's thermochemistry in a parcel with a prescribed density history.

    Parameters
    ----------
    model: jaco.model.Model
    params: dict, optional
        Parameter values by expression name (e.g. {"G_0": 0, "ISRF": 0}) for every evaluation.
    closures: dict, optional
        Parameter name -> function of the state dict, replacing or adding to Closures.DEFAULT (only the parameters the
        model reads are used).
    system: OdeSystem, optional
        A prepared system for this model; building it is the slow part (two minutes for starforge).
    method, rtol, atol: solve_ivp options; atol applies to T in K and to the abundances per H.

    Every mode returns a dict of arrays: t [s], n_Htot, T, rho, P, u, the abundances by species name, the balance
    variables (e.g. Td), heat [erg cm^-3 s^-1], pdv_work and dlnrho_dt, plus "sol" (the solve_ivp result).
    """

    def __init__(self, model, params=None, closures=None, system=None, method="LSODA", rtol=1e-5, atol=1e-12):
        self.model = model
        self.system = system or OdeSystem(model, params)
        self.params = dict(params or {})
        closures = {**Closures.DEFAULT, **(closures or {})}
        self.closures = {k: f for k, f in closures.items() if k in self.system.parameters}
        self.method, self.rtol, self.atol = method, rtol, atol

    @property
    def species(self):
        return self.system.species

    # ------------------------------------------------------------------ the right-hand side

    def _params_at(self, n, T, x, dlnrho_dt):
        p = dict(self.params)
        p["n_Htot"] = n
        if not self.closures:
            return p
        base = self.system.rates(T, x, {**p, "pdv_work": 0.0})  # rho and P at this state, for the closures
        st = {"n_Htot": n, "T": T, "x": x, "rho": base["rho"], "P": base["P"], "u": base["u"], "dlnrho_dt": dlnrho_dt,
              "t_ff": t_freefall(base["rho"]), "lambda_J": jeans_length(base["P"], base["rho"])}
        for k, f in self.closures.items():
            p[k] = f(st)
        return p

    def rates(self, n, T, x, dlnrho_dt=0.0):
        """the system's rates at (n, T, x) under a compression rate d ln rho/dt, with the closures applied"""
        p = self._params_at(n, T, x, dlnrho_dt)
        r = self.system.rates(T, x, p)
        r["params"] = p
        return r

    def _unpack(self, y):
        return y[0], {s: v for s, v in zip(self.species, y[1:])}

    # ------------------------------------------------------------------ modes

    def track(self, dlnrho_dt, n0, T0, t_end, x0=None, n_points=400, verbose=False):
        """Time evolution under a compression law d ln rho/dt = dlnrho_dt(state), the state dict as for the closures
        plus "t". The density follows from it; use free_fall() for the free-fall law."""
        sys_ = self.system
        sys_.reset()
        y0 = np.concatenate([[np.log(n0)], sys_.state(T0, x0 or {})])  # ln n, T, x
        progress = [None]

        def rhs(t, y):
            n = np.exp(y[0])
            T, x = self._unpack(y[1:])
            base = self._params_at(n, T, x, 0.0)
            r0 = sys_.rates(T, x, {**base, "pdv_work": 0.0})
            st = {"t": t, "n_Htot": n, "T": T, "x": x, "rho": r0["rho"], "P": r0["P"], "u": r0["u"],
                  "t_ff": t_freefall(r0["rho"]), "lambda_J": jeans_length(r0["P"], r0["rho"])}
            d = float(dlnrho_dt(st))
            r = sys_.rates(T, x, self._params_at(n, T, x, d))
            self._progress(verbose, progress, n, T, t)
            return np.array([d, r["dTdt"]] + [r["dxdt"][s] for s in self.species])

        t_eval = np.linspace(0, t_end, n_points)
        sol = solve_ivp(rhs, (0, t_end), y0, method=self.method, t_eval=t_eval, rtol=self.rtol, atol=self.atol)
        if not sol.success:
            raise RuntimeError(f"OneZone.track: {sol.message}")
        return self._collect(sol.t, np.exp(sol.y[0]), sol.y[1:], sol, lambda st: dlnrho_dt({**st, "t": st["t"]}))

    def free_fall(self, n0=0.1, T0=300.0, n_stop=1e15, x0=None, n_points=400, verbose=False):
        """Collapse at the free-fall rate, d rho/dt = rho / t_ff, from n0 to n_stop; integrated in ln n."""
        sys_ = self.system
        sys_.reset()
        y0 = sys_.state(T0, x0 or {})
        s_end = np.log(n_stop / n0)
        progress = [None]

        def rhs(s, y):
            n = n0 * np.exp(s)
            T, x = self._unpack(y)
            r0 = sys_.rates(T, x, {**self._params_at(n, T, x, 0.0), "pdv_work": 0.0})
            tff = t_freefall(r0["rho"])
            r = sys_.rates(T, x, self._params_at(n, T, x, 1.0 / tff))
            self._progress(verbose, progress, n, T)
            return tff * np.array([r["dTdt"]] + [r["dxdt"][sp] for sp in self.species])

        s_eval = np.linspace(0, s_end, n_points)
        sol = solve_ivp(rhs, (0, s_end), y0, method=self.method, t_eval=s_eval, rtol=self.rtol, atol=self.atol)
        if not sol.success:
            raise RuntimeError(f"OneZone.free_fall: {sol.message}")
        n = n0 * np.exp(sol.t)
        out = self._collect(None, n, sol.y, sol, lambda st: 1.0 / st["t_ff"])
        tff = t_freefall(out["rho"])
        out["t"] = np.concatenate([[0], np.cumsum(0.5 * (tff[1:] + tff[:-1]) * np.diff(sol.t))])
        return out

    def evolve(self, n, T0, t_end, x0=None, n_points=400, verbose=False):
        """Time evolution at fixed density n from (T0, x0)."""
        return self.track(lambda st: 0.0, n, T0, t_end, x0=x0, n_points=n_points, verbose=verbose)

    def equilibrium(self, n_grid, T0=1e4, x0=None, tol=1e-6, t_max_ff=1e4, verbose=False):
        """Thermochemical equilibrium at each density of n_grid: the state relaxed in time at fixed density until no
        rate would change it by more than tol over a free-fall time (or t_max_ff free-fall times have passed), each
        density warm-started from the previous one's answer, so a bistable model follows the branch the previous
        density was on (as a slow compression would). The returned dict has one entry per density plus "relaxed"
        (whether the tolerance was met) and "t" (the relaxation time used)."""
        sys_ = self.system
        n_grid = np.asarray(n_grid, dtype=float)
        T, x = T0, dict(x0 or {})
        rows = []
        for n in n_grid:
            sys_.reset()
            y0 = sys_.state(T, x)
            r0 = self.rates(n, T, x)
            tff = t_freefall(r0["rho"])
            scale = np.array([None] + [max(sys_.floors[s], 1e-6) for s in self.species], dtype=object)

            def rhs(t, y):
                Ti, xi = self._unpack(y)
                r = sys_.rates(Ti, xi, self._params_at(n, Ti, xi, 0.0))
                return np.array([r["dTdt"]] + [r["dxdt"][s] for s in self.species])

            def settled(t, y):
                d = rhs(t, y)
                rel = [abs(d[0]) * tff / max(y[0], 1.0)] + [abs(d[k]) * tff / max(abs(y[k]), scale[k])
                                                           for k in range(1, len(y))]
                return np.log10(max(rel) + 1e-300) - np.log10(tol)
            settled.terminal, settled.direction = True, -1

            sol = solve_ivp(rhs, (0, t_max_ff * tff), y0, method=self.method, rtol=self.rtol, atol=self.atol,
                            events=settled)
            if not sol.success:
                raise RuntimeError(f"OneZone.equilibrium at n = {n:g}: {sol.message}")
            T, x = self._unpack(sol.y[:, -1])
            x = {s: float(v) for s, v in x.items()}
            rows.append((float(T), x, float(sol.t[-1]), bool(sol.t_events[0].size)))
            if verbose:
                print(f"  n = {n:.3g} cm^-3: T = {T:.4g} K after {sol.t[-1] / tff:.3g} t_ff"
                      f"{'' if rows[-1][3] else ' (not relaxed)'}", file=sys.stderr)
        Ts = np.array([r[0] for r in rows])
        X = np.array([[r[1][s] for r in rows] for s in self.species]).reshape(len(self.species), -1)
        out = self._collect(np.array([r[2] for r in rows]), n_grid, np.vstack([Ts, X]) if len(self.species) else Ts[None, :],
                            None, lambda st: 0.0)
        out["relaxed"] = np.array([r[3] for r in rows])
        return out

    # ------------------------------------------------------------------ helpers

    def _progress(self, verbose, progress, n, T, t=None):
        if not verbose:
            return
        d = int(np.floor(np.log10(n)))
        if progress[0] is None or d > progress[0]:  # trial steps cross decades both ways
            progress[0] = d
            where = f"n = 1e{d} cm^-3" if t is None else f"t = {t:.3g} s, n = {n:.3g} cm^-3"
            print(f"  {where}: T = {T:.4g} K", file=sys.stderr)

    def _collect(self, t, n, Y, sol, dlnrho_dt_of):
        """the track dict: the integrated state plus the derived quantities re-evaluated at every output point"""
        sys_ = self.system
        out = {"n_Htot": np.asarray(n, dtype=float), "T": np.asarray(Y[0], dtype=float)}
        for i, s in enumerate(self.species):
            out[s] = np.asarray(Y[1 + i], dtype=float)
        cols = {k: [] for k in ["rho", "P", "u", "heat", "pdv_work", "dlnrho_dt"] + sys_.balances}
        sys_.reset()
        for k in range(len(out["n_Htot"])):
            x = {s: out[s][k] for s in self.species}
            r0 = sys_.rates(out["T"][k], x, {**self._params_at(out["n_Htot"][k], out["T"][k], x, 0.0), "pdv_work": 0.0})
            st = {"t": t[k] if t is not None else np.nan, "n_Htot": out["n_Htot"][k], "T": out["T"][k], "x": x,
                  "rho": r0["rho"], "P": r0["P"], "u": r0["u"], "t_ff": t_freefall(r0["rho"]),
                  "lambda_J": jeans_length(r0["P"], r0["rho"])}
            d = float(dlnrho_dt_of(st))
            p = self._params_at(out["n_Htot"][k], out["T"][k], x, d)
            r = sys_.rates(out["T"][k], x, p)
            for c in cols:
                cols[c].append(r[c] if c in r else (d if c == "dlnrho_dt" else p.get(c, np.nan)))
        for c, v in cols.items():
            out[c] = np.array(v, dtype=float)
        if t is not None:
            out["t"] = np.asarray(t, dtype=float)
        out["sol"] = sol
        return out


def _parse_grid(text):
    lo, hi, num = text.split(",")
    return np.logspace(np.log10(float(lo)), np.log10(float(hi)), int(num))


def main(argv=None):
    ap = argparse.ArgumentParser(description="one-zone density-temperature tracks of jaco models")
    ap.add_argument("mode", choices=["free_fall", "equilibrium"])
    ap.add_argument("models", nargs="+", help="model names (jaco.models.<name>)")
    ap.add_argument("--n0", type=float, default=0.1, help="free_fall: initial density")
    ap.add_argument("--n-stop", type=float, default=1e15, help="free_fall: final density")
    ap.add_argument("--n-grid", default="1e-2,1e6,33", metavar="LO,HI,NUM", help="equilibrium: log-spaced densities")
    ap.add_argument("--T0", type=float, default=None, help="initial temperature (300 K for free_fall, 1e4 K for equilibrium)")
    ap.add_argument("--param", action="append", default=[], metavar="NAME=VALUE", help="e.g. G_0=0, ISRF=0, Z_d=0.01")
    ap.add_argument("--rtol", type=float, default=1e-5)
    ap.add_argument("-o", "--output", default=None, help="figure path (default <mode>_track.png)")
    ap.add_argument("--save", default=None, help="write the tracks to this .npz")
    args = ap.parse_args(argv)
    params = {}
    for kv in args.param:
        k, v = kv.split("=", 1)
        params[k] = float(v)
    T0 = args.T0 if args.T0 is not None else (300.0 if args.mode == "free_fall" else 1e4)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from importlib import import_module
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    tracks = {}
    for name in args.models:
        model = import_module(f"jaco.models.{name}").make_model()
        print(f"{name}: building the system", file=sys.stderr)
        zone = OneZone(model, params=params, rtol=args.rtol)
        if args.mode == "free_fall":
            tr = zone.free_fall(n0=args.n0, T0=T0, n_stop=args.n_stop, verbose=True)
        else:
            tr = zone.equilibrium(_parse_grid(args.n_grid), T0=T0, verbose=True)
        tracks[name] = {k: v for k, v in tr.items() if isinstance(v, np.ndarray)}
        ax.loglog(tr["n_Htot"], tr["T"], lw=1.6, label=name)
        print(f"{name}: T = {tr['T'][-1]:.4g} K at n = {tr['n_Htot'][-1]:.3g}", file=sys.stderr)
    ax.set_xlabel(r"$n_{\rm H}$ [cm$^{-3}$]")
    ax.set_ylabel("T [K]")
    ax.set_title({"free_fall": "free-fall collapse", "equilibrium": "thermochemical equilibrium"}[args.mode], fontsize=10)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(top=False, right=False, which="both")
    fig.tight_layout()
    out = args.output or f"{args.mode}_track.png"
    fig.savefig(out, dpi=150)
    print(out)
    if args.save:
        np.savez(args.save, **{f"{m}/{k}": v for m, d in tracks.items() for k, v in d.items()})


if __name__ == "__main__":
    main()
