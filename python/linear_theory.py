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


def _cutoff(T, Tfloor):
    """PLUTO's smooth cooling cutoff, f = 1 - sech((T/Tf)^12), and dln f/dln T."""
    y = (T / Tfloor) ** 12
    if y > 700:
        return 1.0, 0.0
    sech = 1.0 / np.cosh(y)
    f = 1.0 - sech
    dlnf = 12.0 * y * sech * np.tanh(y) / max(f, 1e-300)
    return f, dlnf


def integrate(n0, T0, lam_kpc, Z=None, table=None, background="isochoric",
              amp=1.0, heat=0.0, Tfloor=1.0e4, t_end=None, Tstop=None,
              npts=400, rtol=1e-8, t_out=None):
    """
    Integrate the linearised equations.

    n0      : initial cloud density in code units (rho / m_p, i.e. cm^-3)
    T0      : initial cloud temperature [K]
    lam_kpc : initial wavelength [kpc]
    Z       : metallicity (solar units) -- or pass a CoolingTable as `table`
    heat    : uniform volumetric heating rate [erg cm^-3 s^-1]
    amp     : initial delta rho / rho (isobaric seed: q = w = 0)
    t_out   : optional output times [Myr] (default: npts uniform samples)
    Returns a dict of arrays (t in Myr).
    """
    tab = table if table is not None else C.CoolingTable.for_Z(Z)
    g = C.GAMMA
    rho0 = n0 * C.mp
    p0 = C.p_of(rho0, T0)
    km = 2 * np.pi / (lam_kpc * C.kpc) / rho0          # mass wavenumber
    tc0 = p0 / (g - 1) / (C.nenH_of_rho(rho0) * tab(T0))
    if t_end is None:
        t_end = 50.0 * tc0
    # Tstop: stop when the background reaches this temperature (default just
    # above the floor).  Tstop=0 integrates to t_end, i.e. through the floor,
    # which captures the ballistic (post-cooling) phase.
    Tstop = Tstop if Tstop is not None else 1.05 * Tfloor
    cs0 = C.cs_of(T0)

    # y = [ln rho0, ln p0, s, q, w/cs0]
    def rhs(t, y):
        lr, lp, s, q, wn = y
        rho, p = np.exp(lr), np.exp(lp)
        w = wn * cs0
        T = C.T_of(rho, p)
        f, dlnf = _cutoff(T, Tfloor)
        L = C.nenH_of_rho(rho) * tab(T) * f
        LT = tab.dlnL_dlnT(T) + dlnf
        r = (g - 1) * L / p
        hh = (g - 1) * heat / p * f
        if background == "isochoric":
            D = 0.0
            dlp = -(r - hh)
        elif background == "isobaric":
            D = (r - hh) / g
            dlp = 0.0
        else:
            raise ValueError(background)
        k = km * rho
        ds = s * D - k * w
        dw = k * p / rho * q
        dq = q * (r - hh) + g * s * D - g * k * w - r * (2 * s + LT * (q - s))
        return [D, dlp, ds, dq, dw / cs0]

    def hit_floor(t, y):
        return C.T_of(np.exp(y[0]), np.exp(y[1])) - Tstop
    hit_floor.terminal = True

    y0 = [np.log(rho0), np.log(p0), amp, 0.0, 0.0]
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
    tc = p / (g - 1) / (C.nenH_of_rho(rho) * np.maximum(tab(T), 1e-40))
    out = dict(t=t / C.Myr, T=T, n=rho / C.mp, p=p, s=Y[2], q=Y[3],
               w=Y[4] * cs0, k=k, lam_kpc=2 * np.pi / k / C.kpc,
               kcs_tcool=k * C.cs_of(T) * tc, tcool0_Myr=tc0 / C.Myr,
               cs_tcool0_kpc=cs0 * tc0 / C.kpc, LambdaT=tab.dlnL_dlnT(T))
    return out


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
