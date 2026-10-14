"""
Linear theory for a sinusoidal perturbation in a uniformly COOLING background.

Unlike Das, Choudhury & Sharma (2021; heating = cooling on average), the
background here is not in equilibrium: it cools on the same timescale on which
the perturbation evolves.  We therefore integrate the linearised equations
along the time-dependent background instead of looking for exp(-i w t) modes.

Lagrangian formulation (mass coordinate m, fixed mass wavenumber k_m, physical
wavenumber k = k_m rho0(t)).  The perturbation is

    delta rho / rho0 = s(t) cos(k_m m),  delta p / p0 = q(t) cos(k_m m),
    delta v          = w(t) sin(k_m m),

and with  D = dln rho0/dt,  r = (gamma-1) n_e n_H Lambda / p0  (inverse
isochoric cooling time),  h = (gamma-1) H / p0  (heating, uniform per volume),

    ds/dt = s D - k w
    dw/dt = k (p0/rho0) q
    dq/dt = q (r - h) + gamma s D - gamma k w - r [2 s + Lambda_T (q - s)]

Two backgrounds:

  'isochoric'  : infinite uniform medium (the "infinite cloud").  D = 0,
                 dp0/dt = -(gamma-1)(L - H).  k is constant.
  'isobaric'   : cloud much smaller than c_s t_cool, held at the constant
                 pressure of the hot medium.  D = (r - h)/gamma, p0 = const,
                 so the cloud is compressed and k grows with rho0.

Useful limits (pure cooling, h = 0)
  * short wavelength (isobaric mode):  ds/dt = (2 - Lambda_T) r s / gamma.
    The INSTANTANEOUS rate is identical to Das et al.'s isobaric rate, but the
    background evolves, so the total amplification is finite:
        isobaric cloud :  s/s0 = [t_cool,P(T0)/t_cool,P(T)]
                              = (T0/T)^2 Lambda(T)/Lambda(T0)
        infinite medium:  s/s0 = [(T0/T)^2 Lambda(T)/Lambda(T0)]^(1/gamma)
    i.e. the perturbation is just a *time shift* of the cooling history.
  * long wavelength (isochoric mode): s = const (no density growth at all);
    only dq/dt = r(1 - Lambda_T) q - r (2 - Lambda_T) s (pressure/temperature
    perturbation, compare Das: -Lambda_T r).
"""
import numpy as np
from scipy.integrate import solve_ivp

import cooling as C
import cooling_models as CM


def integrate(n0, T0, lam_kpc, Z=None, table=None, background="isochoric",
              amp=1.0, heat=0.0, Tfloor=1.0e4, t_end=None, Tstop=None,
              npts=400, rtol=1e-8, t_out=None, model=None, pert=None,
              stop_factor=1.05):
    """
    Integrate the linearised equations.

    n0      : initial cloud density in code units (rho / m_p, i.e. cm^-3)
    T0      : initial cloud temperature [K]
    lam_kpc : initial wavelength [kpc]
    model   : a cooling model from cooling_models.py (net loss L(rho, T)).
              Default: CIEFloor(Z) (or a CoolingTable passed as `table`),
              i.e. the original collisional cooling with a floor at Tfloor.
    heat    : extra uniform volumetric heating [erg cm^-3 s^-1]
    amp     : perturbation amplitude A
    pert    : type of seed (default "density"):
                "density"     delta rho/rho = A at uniform pressure (s=A)
                "temperature" delta T/T = A at fixed density (q=A)
                "sound"       standing sound wave, s=A, q=gamma A
                "velocity"    delta v = A c_s sin(k m) (w=A c_s)
                "metal"       delta Z/Z = A at uniform rho, p; Z is carried
                              with the gas, so this is a constant forcing
                              -(gamma-1)/p (dL/dlnZ) A in the energy equation
                              (needs a model with dL_dlnZ, e.g. "pie_z")
    Tstop   : stop when the background reaches this temperature.  Default:
              stop_factor x the floor / equilibrium temperature.  Tstop=0 integrates
              to t_end, through the floor (the post-cooling phase).
    t_out   : optional output times [Myr] (default: npts uniform samples)
    Returns a dict of arrays (t in Myr).
    """
    if model is None:
        model = CM.CIEFloor(Z, Tfloor=Tfloor)
        if table is not None:
            model.tab = table
    g = C.GAMMA
    rho0 = n0 * C.mp
    p0 = C.p_of(rho0, T0)
    km = 2 * np.pi / (lam_kpc * C.kpc) / rho0          # mass wavenumber
    Tf_model = getattr(model, "Tfloor", Tfloor)

    def Lnet(rho, T):
        L, Lr, LT = model.derivs(rho, T)
        if heat:
            f, dlnf = CM.cutoff(T, Tf_model)
            L, LT = L - heat * f, LT - heat * f * dlnf
        return L, Lr, LT

    tc0 = p0 / (g - 1) / Lnet(rho0, T0)[0]
    if t_end is None:
        t_end = 50.0 * tc0
    stop_fac = stop_factor if Tstop is None else None
    cs0 = C.cs_of(T0)

    # y = [ln rho0, ln p0, s, q, w/cs0]
    def rhs(t, y):
        lr, lp, s, q, wn = y
        rho, p = np.exp(lr), np.exp(lp)
        w = wn * cs0
        T = C.T_of(rho, p)
        L, Lr, LT = Lnet(rho, T)
        r = (g - 1) * L / p
        if background == "isochoric":
            D, dlp = 0.0, -r
        elif background == "isobaric":
            D, dlp = r / g, 0.0
        else:
            raise ValueError(background)
        k = km * rho
        ds = s * D - k * w
        dw = k * p / rho * q
        dq = q * r + g * s * D - g * k * w - (g - 1) / p * (Lr * s + LT * (q - s))
        if zeta:
            dq -= (g - 1) / p * model.dL_dlnZ(rho, T) * zeta
        return [D, dlp, ds, dq, dw / cs0]

    def hit_floor(t, y):
        rho, p = np.exp(y[0]), np.exp(y[1])
        Ts = stop_fac * model.T_eq(rho) if stop_fac else Tstop
        return C.T_of(rho, p) - Ts
    hit_floor.terminal = True

    pert = pert or "density"
    s0, q0, w0, zeta = {"density": (amp, 0.0, 0.0, 0.0),
                        "temperature": (0.0, amp, 0.0, 0.0),
                        "sound": (amp, g * amp, 0.0, 0.0),
                        "velocity": (0.0, 0.0, amp, 0.0),
                        "metal": (0.0, 0.0, 0.0, amp)}[pert]
    y0 = [np.log(rho0), np.log(p0), s0, q0, w0]
    sol = solve_ivp(rhs, (0, t_end), y0, method="LSODA", rtol=rtol,
                    atol=[1e-10, 1e-10] + [1e-12 * max(abs(amp), 1e-30)] * 3,
                    events=hit_floor, dense_output=True)
    if t_out is not None:
        t = np.asarray(t_out, float) * C.Myr
    else:
        t = np.linspace(0, sol.t[-1], npts)
    if sol.t_events[0].size:                        # append the stopping point
        te = sol.t_events[0][0]
        t = np.append(t[t < te], te) if t_out is None else t[t <= te]
    Y = sol.sol(t)
    rho, p = np.exp(Y[0]), np.exp(Y[1])
    T = C.T_of(rho, p)
    k = km * rho
    Ls = np.array([Lnet(a, b) for a, b in zip(rho, T)])
    tc = p / (g - 1) / np.where(np.abs(Ls[:, 0]) > 0, np.abs(Ls[:, 0]), 1e-300)
    out = dict(t=t / C.Myr, T=T, n=rho / C.mp, p=p, s=Y[2], q=Y[3],
               dT=Y[3] - Y[2], dK=Y[3] - g * Y[2],      # delta T/T, delta ln(p/rho^gamma)
               w=Y[4] * cs0, k=k, lam_kpc=2 * np.pi / k / C.kpc,
               kcs_tcool=k * C.cs_of(T) * tc, tcool0_Myr=tc0 / C.Myr,
               cs_tcool0_kpc=cs0 * tc0 / C.kpc,
               LambdaT=Ls[:, 2] / np.where(Ls[:, 0] != 0, Ls[:, 0], 1e-300))
    return out


T_END_FACTOR = 1.5   # cooling counts as finished once T < 1.5 T_eq (or floor)


def t_cool_end(n0, T0, model, background="isochoric", heat=0.0):
    """Time [Myr] for unperturbed gas starting at (n0, T0) to cool to
    T_END_FACTOR x T_eq.  Used to normalise all time axes ("t_floor").
    (Reaching 1.05 T_eq can take 2-3x longer, because cooling just above the
    floor / equilibrium is slow, so it is a poor marker of the end of cooling.)"""
    o = integrate(n0, T0, 1.0, model=model, background=background, amp=0.0,
                  heat=heat, stop_factor=T_END_FACTOR)
    return o["t"][-1]


def timeshift_amplification(T, T0, tab, background):
    """Analytic short-wavelength (isobaric-mode) amplification s(T)/s0."""
    x = (T0 / T) ** 2 * tab(T) / tab(T0)
    return x if background == "isobaric" else x ** (1.0 / C.GAMMA)


def das_rates(T, tab):
    """Das et al. (thermal-balance) asymptotic growth rates in units of 1/t_cool
    (t_cool = p/((gamma-1) n_e n_H Lambda)):  isobaric (2-L_T)/gamma, isochoric -L_T."""
    LT = tab.dlnL_dlnT(T)
    return (2 - LT) / C.GAMMA, -LT


def cooling_rates(T, tab):
    """Same, for a purely cooling background: isobaric (2-L_T)/gamma for
    delta rho/rho; isochoric (1-L_T) for delta p/p (delta rho/rho frozen)."""
    LT = tab.dlnL_dlnT(T)
    return (2 - LT) / C.GAMMA, 1 - LT
