"""One-zone models: a parcel of gas with the model's thermochemistry and a prescribed density history.

    zone = OneZone(make_model(), params={"G_0": 0, "ISRF": 0})
    ff = zone.free_fall(n0=0.1, T0=300, n_stop=1e15)             # collapse at the free-fall rate
    eq = zone.equilibrium(np.logspace(-2, 6, 33), T0=1e4)          # thermochemical equilibrium at each density
    ev = zone.evolve(n=100, T0=1e4, t_end=1e14)                    # time evolution at fixed density
    tr = zone.track(dlnrho_dt=lambda st: 1e-13, n0=1, T0=1e4, t_end=1e14)   # any compression law

    steps = OneZone(make_model(), host="menon", step=True)          # implicit steps, as GIZMO's one-zone driver
    st = steps.free_fall_steps(n0=0.1, T0=300, n_stop=1e15)         # takes them (a model written as a step needs it)

or, from the command line, python -m jaco.onezone free_fall starforge starforge_legacy -o track.png and
python -m jaco.onezone equilibrium starforge --n-grid 1e-2,1e6,33.

host="menon" closes the parcel of the one-zone collapse test of Menon et al. (2026, Sec. 3.2, Fig. 2; MenonClosures):
a static parcel a Jeans length across with the columns of their KROME driver, no turbulence;
examples/one_zone_collapse.ipynb runs it against their curves.

Every mode integrates the same right-hand side (jaco.ode.OdeSystem): d ln rho/dt is the mode's, the compressional
heating is pdv_work = P d ln rho/dt per unit volume, the velocity gradient is |d ln rho/dt| with no trace-free part
(homologous compression does not clump), and the shielding columns are half a Jeans length, N_H = n_Htot lambda_J / 2
with lambda_J = sqrt(pi c_s^2 / (G rho)). Everything else the model reads (radiation field, cosmic rays, metal
abundances, the dust's ISRF absorption, the optically-thick opacity) is the model's default or the caller's `params`;
`closures` replaces any of the built-in prescriptions by a function of the state. The free-fall mode (Omukai 2000;
Glover & Jappsen 2007; Menon et al. 2025, Fig. 2) integrates in s = ln(n/n0), so the exponential density growth adds
no stiffness; the others integrate in time.

free_fall_steps instead takes one backward-Euler step of the model per Delta_t (OdeSystem(step=True).solve_step), in
the loop of GIZMO's ONE_ZONE_COLLAPSE driver: first order in Delta_t, and the only mode for a model whose rows use the
step itself (starforge_legacy_RT).
"""

import argparse
import sys

import numpy as np
from scipy.integrate import solve_ivp

from .ode import OdeSystem

__all__ = ["OneZone", "Closures", "GizmoHostClosures", "MenonClosures", "t_freefall", "jeans_length"]

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

    GRAD_V_FLOOR = 1e-3  # of the free-fall rate: a parcel at rest keeps a residual gradient (GIZMO's cells never have
    # an exactly zero one), so the LVG escape-probability columns and their derivatives stay finite

    @classmethod
    def grad_v(cls, st):
        return max(abs(st["dlnrho_dt"]), cls.GRAD_V_FLOOR / st["t_ff"])

    @staticmethod
    def grad_v_tf(st):
        return 0.0

    @staticmethod
    def N_H(st):
        return st["n_Htot"] * st["lambda_J"] / 2

    @staticmethod
    def cell_size(st):
        return st["lambda_J"]



class GizmoHostClosures(Closures):
    """The prescriptions GIZMO's host layer (cooling/jaco.cc) would supply for the same parcel, for GIZMO-family
    models: the column is the slab's mass column rho lambda_J/2 in nucleons; the FUV field is the ISRF attenuated by
    dust, 1.7 ISRF exp(-500 Z_d Sigma) (get_FUV_G0); the velocity gradient is GIZMO's Frobenius norm of homologous
    compression, |d ln rho/dt| / sqrt(3); the dust absorbs the CMB and the attenuated ISRF (dust_absorption_rate:
    0.262 eV/cc CMB, 0.31 ISRF IR, 0.66 ISRF optical attenuated by 180 Z_d Sigma and re-emitted into the IR); the
    optically-thick opacity is the Semenov Planck mean at T times the surviving dust, floored at 0.1 cm^2/g, plus
    GIZMO's high-temperature gas opacities above 1500 K (optically_thick_cooling_opacity). Z_d and ISRF are read from
    the parameters in force."""

    PROTONMASS = 1.6726e-24
    EV = 1.602176634e-12
    C_LIGHT = 2.99792458e10
    _dust = None

    @classmethod
    def _dust_functions(cls):
        """kappa_P(T, T) and the dust survival fraction, compiled like the model's expressions"""
        if cls._dust is None:
            import sympy as sp
            from .interpolation import tables_in
            from .models.starforge.dust_opacity import semenov_planck_mean, dust_survival
            from .ode_codegen import emit_source, compile_source
            T = sp.Symbol("T_")
            exprs = [sp.sympify(semenov_planck_mean(T, T)), sp.sympify(dust_survival(T))]
            f = compile_source(emit_source(exprs, ["T_"], [], [], tables_in(exprs), func_name="dust"), "numba", "dust")
            no_params = np.empty(0)
            cls._dust = (lambda t: f(np.array([float(t)]), no_params)[0], lambda t: f(np.array([float(t)]), no_params)[1])
        return cls._dust

    @staticmethod
    def sigma(st):
        """the slab's mass column [g cm^-2]: half a Jeans length of gas at the parcel's density"""
        return st["rho"] * st["lambda_J"] / 2

    @classmethod
    def N_H(cls, st):
        return cls.sigma(st) / cls.PROTONMASS

    @classmethod
    def grad_v(cls, st):
        return max(abs(st["dlnrho_dt"]) / np.sqrt(3.0), cls.GRAD_V_FLOOR / st["t_ff"])

    @classmethod
    def G_0(cls, st):
        p = st["params"]
        return p.get("ISRF", 1.0) * 1.7 * np.exp(-max(p.get("Z_d", 1.0), 1e-4) * cls.sigma(st) * 500.0)

    @classmethod
    def A_dust(cls, st):
        kappa_P, f_surv = cls._dust_functions()
        p = st["params"]
        Zd, isrf = p.get("Z_d", 1.0), p.get("ISRF", 1.0)
        fac = cls.C_LIGHT * st["rho"] * cls.EV  # absorption per unit volume per (eV cm^-3 of radiation) per (cm^2/g)
        e_cmb, e_ir, e_opt = 0.262, 0.31 * isrf, 0.66 * isrf + 7.8e-3
        kappa_cmb = float(kappa_P(2.73)) * Zd * float(f_surv(2.73))
        kappa_ir = float(kappa_P(20.0)) * Zd * float(f_surv(20.0))
        kappa_opt = 180.0 * Zd
        tau_opt = kappa_opt * cls.sigma(st)
        return fac * (e_cmb * kappa_cmb + kappa_opt * e_opt * np.exp(max(-tau_opt, -100.0))
                      + kappa_ir * (-0.5 * np.expm1(max(-tau_opt, -100.0)) * e_opt + e_ir))

    @classmethod
    def kappa_thick(cls, st):
        kappa_P, f_surv = cls._dust_functions()
        p, T, rho = st["params"], st["T"], st["rho"]
        Zd = p.get("Z_d", 1.0)
        kappa = max(0.1, float(kappa_P(T)) * Zd * float(f_surv(T)))
        if T > 1500.0:
            Z_mass, X = 0.0142 * Zd, 0.76
            x_e = st["x"].get("H+", 0.0) + st["x"].get("He+", 0.0) + 2 * st["x"].get("He++", 0.0)
            k_electron = 0.2 * (1 + X)
            k_molecular = 0.1 * Z_mass
            k_hminus = 1.1e-25 * np.sqrt(Z_mass * rho) * T**7.7
            k_kramers = 4.0e25 * (1 + X) * (Z_mass + 0.001) * rho / T**3.5
            k_radiative = k_molecular + 1.0 / (1.0 / max(k_hminus, 1e-300) + 1.0 / (k_electron + k_kramers))
            k_conductive = 2.6e-7 * x_e * T * T / (rho * rho)
            kappa += 1.0 / (1.0 / k_radiative + 1.0 / max(k_conductive, 1e-300))
        return kappa


# grad_v is a classmethod: bound here, after the class exists (its __func__ would take the state as cls)
Closures.DEFAULT = {"pdv_work": Closures.pdv_work, "rho": Closures.rho, "∇v": Closures.grad_v,
                    "∇v_tf": Closures.grad_v_tf, "N_H": Closures.N_H, "N_H_cool": Closures.N_H, "Δx": Closures.cell_size}


GizmoHostClosures.DEFAULT = {
    "pdv_work": Closures.pdv_work, "rho": Closures.rho, "∇v": GizmoHostClosures.grad_v, "∇v_tf": Closures.grad_v_tf,
    "N_H": GizmoHostClosures.N_H, "N_H_cool": GizmoHostClosures.N_H, "Δx": Closures.cell_size,
    "G_0": GizmoHostClosures.G_0, "A_dust": GizmoHostClosures.A_dust, "kappa_thick": GizmoHostClosures.kappa_thick}


class MenonClosures(GizmoHostClosures):
    """The one-zone collapse test of Menon et al. (2026, Sec. 3.2 and Fig. 2), after Omukai (2000): a static parcel
    a Jeans length across, as GIZMO's ONE_ZONE_COLLAPSE driver also sets it up, with the columns of the KROME driver
    that made their figures (tests/popsicle_semenov_photo_cr_full_collapseTest/test.f90 in Sharda's KROME fork):

    - free fall and Jeans length at the gas density, as for every host (Menon's Figs. 2-3 follow from the total gas
      density in t_ff: their cooling then balances the compressional heating where their gas is isothermal; GIZMO's
      driver takes n_H m_H instead);
    - the shielding and line columns over half a Jeans length (KROME's num2col, n lambda_J / 2): N_H, and the FUV
      attenuation and the dust's absorption of the background on it, GizmoHostClosures' rho lambda_J / 2;
    - the dust's and the gas continuum's optical depth over the full Jeans length (KROME's tau = rho kappa lambda_J):
      N_H_cool = rho lambda_J / m_p;
    - no flow: no turbulence or clumping (grad_v_tf = 0, and the cell size Delta_x = 0, since Delta_x enters these
      models only multiplied by a velocity gradient, as the sub-grid turbulence of the clumping factor and the H2
      line width). The one velocity gradient left is the CO line's: Menon's CO cooling (Omukai et al. 2010 fits in
      N(CO) / v_th,CO, N(CO) from num2col) is a static column with a thermal line width, which an LVG rate reproduces
      with the Sobolev length set to that column's lambda_J / 2, |grad v| = 2 v_th,CO / lambda_J,
      v_th,CO = sqrt(2 k T / m_CO);
    - the FUV field, the dust's absorption of the background and the thick opacity are GizmoHostClosures' (with
      ISRF = 0, Menon's chi_FUV = 0, only the CMB heats the dust)."""

    BOLTZMANN = 1.380649e-16
    M_CO = 28.0 * GizmoHostClosures.PROTONMASS

    @classmethod
    def N_H_cool(cls, st):
        """the dust's column [cm^-2]: a full Jeans length of gas, in nucleons"""
        return st["rho"] * st["lambda_J"] / cls.PROTONMASS

    @classmethod
    def grad_v(cls, st):
        return 2 * np.sqrt(2 * cls.BOLTZMANN * st["T"] / cls.M_CO) / st["lambda_J"]

    @staticmethod
    def cell_size(st):
        return 0.0


MenonClosures.DEFAULT = {  # rebound, so that the classmethods see MenonClosures
    "pdv_work": Closures.pdv_work, "rho": Closures.rho, "∇v": MenonClosures.grad_v, "∇v_tf": Closures.grad_v_tf,
    "N_H": MenonClosures.N_H, "N_H_cool": MenonClosures.N_H_cool, "Δx": MenonClosures.cell_size,
    "G_0": MenonClosures.G_0, "A_dust": MenonClosures.A_dust, "kappa_thick": MenonClosures.kappa_thick}



class OneZone:
    """A model's thermochemistry in a parcel with a prescribed density history.

    Parameters
    ----------
    model: jaco.model.Model
    params: dict, optional
        Parameter values by expression name (e.g. {"G_0": 0, "ISRF": 0}) for every evaluation.
    closures: dict, optional
        Parameter name -> function of the state dict, replacing or adding to the host's set (only the parameters the
        model reads are used).
    backend: "numba", "numpy" or "jax"
        How the model's expressions are evaluated (jaco.ode_codegen); numba is ~100x faster than the others per call.
    analytic_jacobian: bool
        Hand the integrator the system's Jacobian (OdeSystem.jacobian) instead of letting it difference the
        right-hand side.
    host: "gizmo", "menon" or "plain"
        Which prescriptions close the parameters a host would supply: GizmoHostClosures (what cooling/jaco.cc
        computes for the cell: attenuated FUV field, dust absorption of the CMB and ISRF, the thick opacity, the
        nucleon column), MenonClosures (the same, in the static parcel of Menon et al. 2026, with their columns) or
        the bare Closures (constant fields, no dust absorption, the model's default opacity).
    system: OdeSystem, optional
        A prepared system for this model; building it is the slow part (two minutes for starforge).
    method, rtol, atol: solve_ivp options; atol applies to T in K and to the abundances per H.
    step: bool
        Build the system in step mode (OdeSystem(step=True)) for free_fall_steps(), the implicit stepping that a model
        written as one backward-Euler step needs (starforge_legacy_RT); the integrating modes need step=False.

    Every mode returns a dict of arrays: t [s], n_Htot, T, rho, P, u, the abundances by species name, the balance
    variables (e.g. Td), heat [erg cm^-3 s^-1], pdv_work and dlnrho_dt, plus "sol" (the solve_ivp result).
    """

    def __init__(self, model, params=None, closures=None, system=None, method="LSODA", rtol=1e-5, atol=1e-12,
                 host="gizmo", backend="numba", analytic_jacobian=True, step=False):
        self.model = model
        self.system = system or OdeSystem(model, params, backend=backend, step=step)
        self.analytic_jacobian = analytic_jacobian
        self.params = dict(params or {})
        base = {"gizmo": GizmoHostClosures, "menon": MenonClosures, "plain": Closures}[host].DEFAULT
        closures = {**base, **(closures or {})}
        self.closures = {k: f for k, f in closures.items() if k in self.system.parameters}
        self.method, self.rtol, self.atol = method, rtol, atol

    @property
    def species(self):
        return self.system.species

    ELEMENT_TOTALS = ("x_C,tot", "x_O,tot", "x_N", "x_Ne", "x_Mg", "x_Si", "x_S", "x_Ca", "x_Fe", "Z", "Z_d")

    def metallicity(self, Z):
        """parameter overrides for a metallicity Z in solar units: every heavy-element total the model reads and the
        dust abundance Z_d scaled from the system's defaults (dust-to-gas proportional to Z)"""
        d = self.system.defaults
        return {k: d[k] * Z for k in self.ELEMENT_TOTALS if k in d}

    # ------------------------------------------------------------------ the right-hand side

    def _params_at(self, n, T, x, dlnrho_dt):
        """the parameters at a state: the fixed ones, n_Htot and the closures (from rho and P at the state as given)"""
        p = dict(self.params)
        p["n_Htot"] = n
        if not self.closures:
            return p
        th = self.system.thermo(T, x, p)
        st = {"n_Htot": n, "T": T, "x": th["x"], "rho": th["rho"], "P": th["P"], "u": th["u"], "dlnrho_dt": dlnrho_dt,
              "t_ff": t_freefall(th["rho"]), "lambda_J": jeans_length(th["P"], th["rho"]),
              "params": {**self.system.defaults, **p}}
        for k, f in self.closures.items():
            p[k] = f(st)
        return p

    def rates(self, n, T, x, dlnrho_dt=0.0):
        """the system's rates at (n, T, x) under a compression rate d ln rho/dt, with the closures applied"""
        p = self._params_at(n, T, x, dlnrho_dt)
        r = self.system.rates(T, x, p)
        r["params"] = p
        return r

    # The integrator's state is w = (ln u, time-dependent abundances): a trial step cannot make u negative, which has
    # no temperature. OdeSystem works with y = (u, abundances).
    @staticmethod
    def _y(w):
        return np.concatenate([[np.exp(w[0])], w[1:]])

    def _state_at(self, n, w, law, t=None):
        """the full state behind an integrator state w = (ln u, time-dependent abundances) at density n: T and the
        steady-state species solved from u, the compression rate from the law, the parameters from the closures.
        The closures depend on T (through the Jeans length), so the solve is repeated once with them updated.
        Returns T, x, params, d ln rho/dt."""
        sys_ = self.system
        u, xi = float(np.exp(w[0])), {s: float(v) for s, v in zip(sys_.integrated, w[1:])}
        T = sys_._last_T if sys_._last_T is not None else 1e3
        x = xi
        d = 0.0
        for _ in range(2):
            p = self._params_at(n, T, x, d)
            T, x = sys_.temperature(u, xi, p)
            th = sys_.thermo(T, x, p)
            st = {"t": t, "n_Htot": n, "T": T, "x": x, "rho": th["rho"], "P": th["P"], "u": th["u"],
                  "t_ff": t_freefall(th["rho"]), "lambda_J": jeans_length(th["P"], th["rho"]),
                  "params": {**sys_.defaults, **p}}
            d = float(law(st))
        p = self._params_at(n, T, x, d)
        return T, x, p, d

    def _jacobian_at(self, n, w, law, t=None):
        """the Jacobian of the integrator's right-hand side (in w = (ln u, x)) at a state, the closures' dependence on
        the state chained in from their partial derivatives in T and the abundances (central differences of the closure
        formulas at fixed density; cheap, no solves)"""
        sys_ = self.system
        y = self._y(w)
        T, x, p, d = self._state_at(n, w, law, t)
        closed = [q for q in sys_._deriv_params if q in self.closures]
        dparams = {}
        if closed:
            dparams = {q: {} for q in closed}
            for c in ["T"] + [f"x_{s}" for s in sys_.species]:
                if c == "T":
                    h = 1e-6 * T
                    pp, pm = self._params_at(n, T + h, x, d), self._params_at(n, T - h, x, d)
                else:
                    v = x[c[2:]]
                    h = 1e-6 * max(v, 1e-12)
                    pp = self._params_at(n, T, {**x, c[2:]: v + h}, d)
                    pm = self._params_at(n, T, {**x, c[2:]: max(v - h, 0.0)}, d)
                    h = (v + h) - max(v - h, 0.0)
                    h = h / 2
                for q in closed:
                    dparams[q][c] = (pp[q] - pm[q]) / (2 * h)
        J = sys_.jacobian(y, p, dparams)
        if np.all(np.isfinite(J)):
            u = y[0]  # d(du/dt / u)/d(ln u) = J_uu - (du/dt)/u; the rest rescales rows and columns
            J[0, :] /= u
            J[:, 0] *= u
            J[0, 0] -= sys_.rates(T, x, p, solved=True)["dudt"] / u
        else:  # a singular parameter derivative at this state: difference the whole thing
            J = np.empty((len(w), len(w)))
            for k in range(len(w)):
                h = 1e-6 * max(abs(w[k]), 1e-30)
                e = np.zeros(len(w))
                e[k] = h
                J[:, k] = (self._rhs_at(n, w + e, law, t) - self._rhs_at(n, w - e, law, t)) / (2 * h)
            sys_.temperature(y[0], {s: v for s, v in zip(sys_.integrated, y[1:])}, p)
        return J, d

    def _rhs_at(self, n, w, law, t=None):
        """d ln u/dt, dx/dt of the integrator's state at density n"""
        sys_ = self.system
        T, x, p, d = self._state_at(n, w, law, t)
        r = sys_.rates(T, x, p, solved=True)
        return np.array([r["dudt"] / r["u"]] + [r["dxdt"][s] for s in sys_.integrated])

    def _initial(self, n, T0, x0):
        """the integrator state w = (ln u, x) at (n, T0, x0)"""
        sys_ = self.system
        sys_.reset()
        p = self._params_at(n, T0, x0 or {}, 0.0)
        y = sys_.state(T0, x0 or {}, p)
        return np.concatenate([[np.log(y[0])], y[1:]])

    # ------------------------------------------------------------------ modes

    def track(self, dlnrho_dt, n0, T0, t_end, x0=None, n_points=400, verbose=False):
        """Time evolution under a compression law d ln rho/dt = dlnrho_dt(state), the state dict as for the closures
        plus "t". The density follows from it; use free_fall() for the free-fall law."""
        sys_ = self.system
        y0 = np.concatenate([[np.log(n0)], self._initial(n0, T0, x0)])  # ln n, ln u, x
        progress = [None]

        def rhs(t, y):
            n = np.exp(y[0])
            T, x, p, d = self._state_at(n, y[1:], dlnrho_dt, t)
            r = sys_.rates(T, x, p, solved=True)
            self._progress(verbose, progress, n, T, t)
            return np.array([d, r["dudt"] / r["u"]] + [r["dxdt"][s] for s in sys_.integrated])

        def jac(t, y):
            """analytic in (u, x); the ln n column and row by finite differences (the law and the density enter the
            rates as parameters)"""
            m = len(y)
            J = np.zeros((m, m))
            n = np.exp(y[0])
            T, x, p, d = self._state_at(n, y[1:], dlnrho_dt, t)
            J[1:, 1:] = self._jacobian_at(n, y[1:], dlnrho_dt, t)[0]
            h = 1e-6
            f0, f1 = rhs(t, y), rhs(t, y + np.array([h] + [0.0] * (m - 1)))
            J[:, 0] = (f1 - f0) / h
            return J

        t_eval = np.linspace(0, t_end, n_points)
        sol = solve_ivp(rhs, (0, t_end), y0, method=self.method, t_eval=t_eval, rtol=self.rtol, atol=self.atol,
                        jac=jac if self.analytic_jacobian else None)
        if not sol.success:
            raise RuntimeError(f"OneZone.track: {sol.message}")
        return self._collect(sol.t, np.exp(sol.y[0]), sol.y[1:], sol, dlnrho_dt)

    def free_fall(self, n0=0.1, T0=300.0, n_stop=1e15, x0=None, n_points=400, verbose=False, max_dlnn=0.5):
        """Collapse at the free-fall rate, d rho/dt = rho / t_ff, from n0 to n_stop; integrated in ln n.

        max_dlnn caps the integrator's step in ln n. Where nothing but compression acts (no electrons, metals or
        radiation, e.g. equilibrium ions at Z = 0) the right-hand side is so smooth that the integrator's first trial
        step otherwise reaches n_stop, where the adiabat lies far outside the temperatures the state solve accepts."""
        sys_ = self.system
        y0 = self._initial(n0, T0, x0)
        s_end = np.log(n_stop / n0)
        progress = [None]
        law = lambda st: 1.0 / st["t_ff"]

        def rhs(s, y):
            n = n0 * np.exp(s)
            T, x, p, d = self._state_at(n, y, law)
            r = sys_.rates(T, x, p, solved=True)
            self._progress(verbose, progress, n, T)
            return np.array([r["dudt"] / r["u"]] + [r["dxdt"][sp] for sp in sys_.integrated]) / d

        def jac(s, y):
            J, d = self._jacobian_at(n0 * np.exp(s), y, law)
            return J / d

        s_eval = np.linspace(0, s_end, n_points)
        sol = solve_ivp(rhs, (0, s_end), y0, method=self.method, t_eval=s_eval, rtol=self.rtol, atol=self.atol,
                        jac=jac if self.analytic_jacobian else None, max_step=max_dlnn)
        if not sol.success:
            raise RuntimeError(f"OneZone.free_fall: {sol.message}")
        n = n0 * np.exp(sol.t)
        out = self._collect(None, n, sol.y, sol, law)
        tff = t_freefall(out["rho"])
        out["t"] = np.concatenate([[0], np.cumsum(0.5 * (tff[1:] + tff[:-1]) * np.diff(sol.t))])
        return out

    # ------------------------------------------------------------------ implicit stepping (step mode)

    SUBSTEP_MAX = 2000  # adaptive substeps of a failed step, and their smallest fraction of it (GIZMO's JACO solver's
    SUBSTEP_MIN_FRAC = 1e-6  # tier 3: JACO_SUBCYCLE_MAXSTEPS, JACO_SUBCYCLE_MIN_FRAC)

    def _bands(self):
        kinds = self.model.network.species_kinds
        return [s for s in self.system.species if kinds.get(s) == "radiation"]

    def _band_absorption(self):
        """band -> function of the parameter dict giving the band's absorption coefficient kappa_abs rho [cm^-1], as
        the model's absorption processes take it (jaco.models.starforge.radiation): the dust's at f_abs = 1/2 for the
        photoelectric, NUV and optical bands, the dust's and the gas's absorption opacities for the IR band, sigma_HI on
        the neutral H for the ionizing band"""
        if getattr(self, "_absorption", None) is not None:
            return self._absorption
        import sympy as sp
        from .interpolation import tables_in
        from .ode_codegen import emit_source, compile_source
        from .symbols import sanitize_name, sanitize_symbols
        from .models.starforge import radiation as rt
        from .models.starforge.dust_opacity import T_dust, T_rad, rho, ir_dust_opacity, ir_gas_opacity, band_dust_opacity
        exprs = {b: 0.5 * band_dust_opacity(k, f, T_d=T_dust) * rho for b, (k, f) in rt.DUST_BAND_OPACITY.items()}
        exprs[rt.IR] = (ir_dust_opacity(T_dust, T_rad) + ir_gas_opacity(T_rad, T_dust)) * rho
        exprs = {b: sp.sympify(e) for b, e in exprs.items() if b in self._bands()}
        names = sorted({str(s) for e in exprs.values() for s in e.free_symbols})
        f = compile_source(emit_source([sanitize_symbols(e) for e in exprs.values()], [],
                                       [sanitize_name(q) for q in names], [], tables_in(list(exprs.values())),
                                       func_name="band_absorption"), "numba", "band_absorption")
        order = list(exprs)

        defaults = self.system.defaults

        def absorption(p):
            """p: parameters and the state (T, Td, x_<species>); the eliminated abundances the opacities read are
            filled from the solved ones: x_H = 1 - x_H+ - 2 x_H2, x_He = y - x_He+ - x_He++ and the electrons as the
            ions give them (x_H+ + x_He+ + 2 x_He++; the metals' electrons left out: they matter to the IR opacity
            only in ionized gas, where the ions dominate)"""
            p = {**defaults, **p}
            g = lambda k: float(p.get(k, 0.0))  # noqa: E731
            p.setdefault("x_H", max(1.0 - g("x_H+") - 2 * g("x_H_2"), 0.0))
            p.setdefault("x_He", max(g("y") - g("x_He+") - g("x_He++"), 0.0))
            p.setdefault("x_e-", g("x_H+") + g("x_He+") + 2 * g("x_He++"))
            out = f(np.empty(0), np.array([float(p[q]) for q in names]))
            a = dict(zip(order, (float(v) for v in out)))
            if rt.EUV in self._bands():
                a[rt.EUV] = float(p["sigma_HI"]) * float(p["n_Htot"]) * max(1.0 - p["x_H+"], 0.0)
            return a
        self._absorption = absorption
        return absorption

    def _solve_with_substeps(self, T, x, balance, p):
        """one implicit step over p["Δt"]; if Newton fails, adaptive substeps as GIZMO's JACO solver takes them (tier 3:
        half the step first, the substep doubled after a success and halved after a failure; between substeps the
        start-of-step energy and time-dependent abundances move to the substep's answer, everything else stays).
        Returns (T, x, balance, outputs, iterations, substeps)."""
        sys_ = self.system
        res = sys_.solve_step(T, x, balance, p)
        if res is not None:
            return (*res, 0)
        dt = p["Δt"]
        done, sub, n_sub, its = 0.0, 0.5 * dt, 0, 0
        cur = (T, dict(x), dict(balance), {})
        ps = dict(p)
        while done < dt * (1 - 1e-12):
            if n_sub >= self.SUBSTEP_MAX or sub < self.SUBSTEP_MIN_FRAC * dt:
                raise RuntimeError(f"OneZone.free_fall_steps: the step failed at n_Htot = {p['n_Htot']:.4g}, "
                                   f"T = {T:.4g} K, Delta_t = {dt:.4g} s, even in substeps")
            ps["Δt"] = min(sub, dt - done)
            if n_sub > 0:
                ps["u_initial"] = sys_.thermo(cur[0], cur[1], ps)["u"]
                for s in sys_.integrated:
                    ps[f"x_{s}_initial"] = cur[1][s]
            r = sys_.solve_step(cur[0], cur[1], cur[2], ps)
            if r is None:
                sub *= 0.5
                continue
            cur = r[:4]
            its += r[4]
            done += ps["Δt"]
            n_sub += 1
            sub *= 2
        return (*cur, its, n_sub)

    def free_fall_steps(self, n0=0.1, T0=300.0, n_stop=1e15, x0=None, Td0=None, T_rad0=None, f_ff=0.05, f_cool=0.5,
                        f_rt=0.4, every=1, verbose=False, max_steps=None):
        """Collapse at the free-fall rate in implicit steps, the loop of GIZMO's one-zone driver (ONE_ZONE_COLLAPSE,
        cooling/one_zone.cc, one_zone_run_track), so that the two take the same steps. Needs OneZone(step=True).

        Each step: Delta_t = min(f_ff t_ff, f_cool t_cool, f_rt lambda_J / c_tilde) (the last with radiation bands,
        c_tilde = rsol c), t_cool = u / |net heating without the PdV work| at the current state over the previous step;
        the density by exact free fall over Delta_t, rho_new = rho / (1 - Delta_t / (2 t_ff))^2; the closures at the new
        density, the current temperature and abundances (so lambda_J at the current P/rho, and the compression work
        pdv_work = rho_new (P/rho) d ln rho / Delta_t, as GIZMO's DtInternalEnergy); with radiation bands, each band's
        excess over its starting background (an energy, or photon, density held fixed) decays by exp(-Delta_t / t_esc),
        t_esc = lambda_J (1 + tau) / c_tilde, tau = kappa_abs rho lambda_J, and the IR band's T_rad is re-weighted, as
        one_zone_relax_radiation does; then one backward-Euler step of the model (OdeSystem.solve_step) from the current
        state, in adaptive substeps if it fails (_solve_with_substeps). With bands, the host inputs follow GIZMO's host
        layer (cooling/jaco.cc, jaco_pack_bands and jaco_write_bands): T_rad, the start-of-step dust temperature
        Td_initial, the adiabatic index gamma_eos = 1 + P/(rho u) of the current state, and T_rad_new carried to the
        next step. rsol, Z_metals and the band constants are the caller's params.

        x0 holds the starting abundances (bands included, per H nucleus); Td0 and T_rad0 the starting dust and IR band
        temperatures (by default the dust's balance at the start, and T_CMB). every: keep every so many steps (the
        last always). max_steps: stop after so many steps, wherever the density is.

        Returns a dict of arrays, one row per kept step: t [s], n_Htot, T, u, rho, P, every abundance and balance
        variable, T_rad, and t_ff, t_cool, lambda_J and t_esc (inf without bands) at the row, with dt, the Newton
        iterations and the substeps of the step that led to it; and "steps", the number of steps taken."""
        sys_ = self.system
        if not sys_.step_mode:
            raise RuntimeError("free_fall_steps needs OneZone(..., step=True)")
        bands = self._bands()
        x0 = dict(x0 or {})
        rsol_c = None
        # the starting state: steady-state species and balance variables solved at T0
        p = self._params_at(n0, T0, x0, 0.0)
        sys_.reset()
        r = sys_.rates(T0, x0, p)
        T, x = T0, dict(r["x"])
        balance = {b: r[b] for b in sys_.balances}
        if Td0 is not None and "Td" in balance:
            balance["Td"] = Td0
        T_rad = T_rad0 if T_rad0 is not None else p.get("T_CMB", 2.73)
        background = {b: x[b] * n0 for b in bands}  # per unit volume, held fixed
        T_bg = T_rad
        n, t, dt, step = n0, 0.0, None, 0
        cols = {k: [] for k in ["t", "n_Htot", "T", "u", "rho", "P", "dt", "t_ff", "t_cool", "lambda_J", "t_esc",
                                "iterations", "substeps", "T_rad"] + list(sys_.species) + list(sys_.balances)}
        its = n_sub = 0
        progress = [None]
        while True:
            th = sys_.thermo(T, x, p)
            rho, P, u = th["rho"], th["P"], th["u"]
            t_ff, lam = t_freefall(rho), jeans_length(P, rho)
            # the cooling time: the net heating without the compression, over the previous step (GIZMO's
            # GetNetCoolingRate with OneZoneTimestep = the last step; f_ff t_ff before the first)
            pc = dict(p)
            pc.update({"Δt": dt if dt else f_ff * t_ff, "u_initial": u, "pdv_work": 0.0})
            pc.update({f"x_{s}_initial": x[s] for s in sys_.integrated})
            if bands:
                pc.update({"T_rad": T_rad, "Td_initial": balance.get("Td", T)})
            heats = sys_.heat_breakdown(T, x, pc, solved=True)
            net = sum(v for k, v in heats.items() if k != "PdV work")
            t_cool = rho * u / abs(net) if net != 0 else np.inf
            t_esc = np.inf
            if bands:
                rsol_c = float(pc["rsol"]) * 2.99792458e10
                a = self._band_absorption()({**pc, "Td": balance.get("Td", T), "T": T,
                                             **{f"x_{s}": v for s, v in x.items()}})
                t_esc = min(lam * (1 + a[b] * lam) / rsol_c for b in bands)
            finished = n >= n_stop or (max_steps is not None and step >= max_steps)
            if step % every == 0 or finished:
                for k, v in (("t", t), ("n_Htot", n), ("T", T), ("u", u), ("rho", rho), ("P", P),
                             ("dt", dt or 0.0), ("t_ff", t_ff), ("t_cool", t_cool), ("lambda_J", lam),
                             ("t_esc", t_esc), ("iterations", its), ("substeps", n_sub), ("T_rad", T_rad)):
                    cols[k].append(v)
                for s in sys_.species:
                    cols[s].append(x[s])
                for b in sys_.balances:
                    cols[b].append(balance[b])
            self._progress(verbose, progress, n, T, t)
            if finished:
                break

            dt = min(f_ff * t_ff, f_cool * t_cool)
            if bands:
                dt = min(dt, f_rt * lam / rsol_c)
            contraction = 1.0 - 0.5 * dt / t_ff
            rho_new = rho / contraction**2
            dlnrho = np.log(rho_new / rho)
            n_new = n * rho_new / rho
            p = self._params_at(n_new, T, x, dlnrho / dt)
            p.update({"Δt": dt, "u_initial": u})
            p.update({f"x_{s}_initial": x[s] for s in sys_.integrated})
            if bands:
                # the transport's stand-in: what the zone holds above the background escapes over t_esc
                th_new = sys_.thermo(T, x, p)
                lam_new = jeans_length(th_new["P"], th_new["rho"])
                a = self._band_absorption()({**p, "T_rad": T_rad, "Td": balance.get("Td", T), "T": T,
                                             **{f"x_{s}": v for s, v in x.items()}})
                for b in bands:
                    x_bg = background[b] / n_new
                    excess = (x[b] - x_bg) * np.exp(-dt / (lam_new * (1 + a[b] * lam_new) / rsol_c))
                    if b == "photon_IR":
                        T_rad = (x_bg + excess) / (x_bg / T_bg + excess / max(T_rad, 1e-300)) if excess > 0 else T_bg
                    x[b] = x_bg + excess
                    p[f"x_{b}_initial"] = x[b]
                p.update({"T_rad": T_rad, "Td_initial": balance.get("Td", T), "gamma_eos": 1.0 + P / (rho * u)})
            T, x, balance, outputs, its, n_sub = self._solve_with_substeps(T, x, balance, p)
            if bands and "T_rad_new" in outputs:
                T_rad = max(float(outputs["T_rad_new"]), p.get("T_CMB", 2.73))
            n, t = n_new, t + dt
            step += 1
        out = {k: np.array(v, dtype=float) for k, v in cols.items()}
        out["steps"] = step
        return out

    def evolve(self, n, T0, t_end, x0=None, n_points=400, verbose=False):
        """Time evolution at fixed density n from (T0, x0)."""
        return self.track(lambda st: 0.0, n, T0, t_end, x0=x0, n_points=n_points, verbose=verbose)

    def equilibrium(self, n_grid, T0=1e4, x0=None, tol=1e-4, t_max_ff=100.0, x_floor=1e-6, verbose=False):
        """Thermochemical equilibrium at each density of n_grid: the state relaxed in time at fixed density until no
        rate would change u by more than tol of itself, or a time-dependent abundance by more than tol of itself or
        x_floor per H (whichever is larger), over a free-fall time; or until t_max_ff free-fall times have passed.
        Each density is warm-started from the previous one's answer, so a bistable model follows the branch the
        previous density was on (as a slow compression would). The returned dict has one entry per density plus
        "relaxed" (whether the tolerance was met) and "t" (the relaxation time used). Trace species jitter at the
        integrator's absolute tolerance; x_floor keeps that from counting as change."""
        sys_ = self.system
        n_grid = np.asarray(n_grid, dtype=float)
        T, x = T0, dict(x0 or {})
        rows = []
        still = lambda st: 0.0
        for n in n_grid:
            y0 = self._initial(n, T, x)
            tff = t_freefall(sys_.thermo(T, x, {**self.params, "n_Htot": n})["rho"])
            scale = [None] + [max(sys_.floors[s], x_floor) for s in sys_.integrated]

            def rhs(t, y):
                Ti, xi, p, _ = self._state_at(n, y, still)
                r = sys_.rates(Ti, xi, p, solved=True)
                return np.array([r["dudt"] / r["u"]] + [r["dxdt"][s] for s in sys_.integrated])

            def jac(t, y):
                return self._jacobian_at(n, y, still)[0]

            def settled(t, y):
                d = rhs(t, y)
                rel = [abs(d[0]) * tff] + [abs(d[k]) * tff / max(abs(y[k]), scale[k]) for k in range(1, len(y))]
                return np.log10(max(rel) + 1e-300) - np.log10(tol)
            settled.terminal, settled.direction = True, -1

            sol = solve_ivp(rhs, (0, t_max_ff * tff), y0, method=self.method, rtol=self.rtol, atol=self.atol,
                            events=settled, jac=jac if self.analytic_jacobian else None)
            if not sol.success:
                raise RuntimeError(f"OneZone.equilibrium at n = {n:g}: {sol.message}")
            T, x, _, _ = self._state_at(n, sol.y[:, -1], still)
            rows.append((sol.y[:, -1].copy(), float(sol.t[-1]), bool(sol.t_events[0].size)))
            if verbose:
                print(f"  n = {n:.3g} cm^-3: T = {T:.4g} K after {sol.t[-1] / tff:.3g} t_ff"
                      f"{'' if rows[-1][2] else ' (not relaxed)'}", file=sys.stderr)
        out = self._collect(np.array([r[1] for r in rows]), n_grid, np.array([r[0] for r in rows]).T, None, still)
        out["relaxed"] = np.array([r[2] for r in rows])
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

    def _collect(self, t, n, Y, sol, law):
        """the track dict: the integrated state (u, time-dependent abundances) resolved into T and every abundance,
        plus the derived quantities re-evaluated at every output point"""
        sys_ = self.system
        n = np.asarray(n, dtype=float)
        cols = {k: [] for k in ["T", "u", "rho", "P", "heat", "pdv_work", "dlnrho_dt"] + list(self.species) + sys_.balances}
        sys_.reset()
        for k in range(len(n)):
            T, x, p, d = self._state_at(n[k], Y[:, k], law, t[k] if t is not None else None)
            r = sys_.rates(T, x, p, solved=True)
            for c in cols:
                cols[c].append(r[c] if c in r else r["x"][c] if c in r["x"] else T if c == "T" else
                               d if c == "dlnrho_dt" else p.get(c, np.nan))
        out = {"n_Htot": n}
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
    ap.add_argument("--Z", default=None, metavar="Z1,Z2,...", help="metallicities in solar units: one track per value "
                    "(heavy-element totals and Z_d scaled from the model's defaults)")
    ap.add_argument("--jobs", type=int, default=1, help="tracks run concurrently (each builds its own system)")
    ap.add_argument("--host", default="gizmo", choices=["gizmo", "menon", "plain"],
                    help="closures: GIZMO's host layer, Menon et al. (2026)'s Jeans-column parcel, or the bare set")
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
    Zs = [float(z) for z in args.Z.split(",")] if args.Z else [None]
    jobs = [(name, Z) for name in args.models for Z in Zs]
    common = dict(mode=args.mode, params=params, rtol=args.rtol, T0=T0, n0=args.n0, n_stop=args.n_stop, n_grid=args.n_grid,
                  host=args.host)
    if args.jobs > 1 and len(jobs) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            results = list(ex.map(_run_job, jobs, [common] * len(jobs)))
    else:
        results = [_run_job(j, common) for j in jobs]
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    tracks = {}
    cmap = plt.get_cmap("viridis")
    for (name, Z), tr in zip(jobs, results):
        if tr is None:
            continue
        key = name if Z is None else f"{name}/Z={Z:g}"
        tracks[key] = tr
        if Z is None:
            ax.loglog(tr["n_Htot"], tr["T"], lw=1.6, label=name)
        else:
            c = cmap(0.1 + 0.8 * (np.log10(max(Z, 1e-7)) + 7) / 7)
            ax.loglog(tr["n_Htot"], tr["T"], lw=1.4, color=c, label=f"Z = {Z:g} Z$_\\odot$" if Z > 0 else "Z = 0")
    ax.set_xlabel(r"$n_{\rm H}$ [cm$^{-3}$]")
    ax.set_ylabel("T [K]")
    title = {"free_fall": "free-fall collapse", "equilibrium": "thermochemical equilibrium"}[args.mode]
    ax.set_title(title + (f": {args.models[0]}" if args.Z and len(args.models) == 1 else ""), fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    fig.tight_layout()
    out = args.output or f"{args.mode}_track.png"
    fig.savefig(out, dpi=150)
    print(out)
    if args.save:
        np.savez(args.save, **{f"{m}/{k}": v for m, d in tracks.items() for k, v in d.items()})


def _run_job(job, common):
    """one track (a worker builds its own system); returns the arrays or None on failure"""
    from importlib import import_module
    name, Z = job
    try:
        model = import_module(f"jaco.models.{name}").make_model()
        zone = OneZone(model, params=common["params"], rtol=common["rtol"], host=common.get("host", "gizmo"))
        params = dict(common["params"])
        if Z is not None:
            params.update(zone.metallicity(Z))
            zone.params = params
        label = name if Z is None else f"{name} Z={Z:g}"
        print(f"{label}: system built", file=sys.stderr, flush=True)
        if common["mode"] == "free_fall":
            tr = zone.free_fall(n0=common["n0"], T0=common["T0"], n_stop=common["n_stop"])
        else:
            tr = zone.equilibrium(_parse_grid(common["n_grid"]), T0=common["T0"])
        print(f"{label}: T = {tr['T'][-1]:.4g} K at n = {tr['n_Htot'][-1]:.3g}", file=sys.stderr, flush=True)
        return {k: v for k, v in tr.items() if isinstance(v, np.ndarray)}
    except Exception as e:  # one failed track must not lose the others
        print(f"{job}: FAILED: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
        return None


if __name__ == "__main__":
    main()
