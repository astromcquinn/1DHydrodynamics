"""
Cooling function, composition constants and helper rates.

All conventions here mirror pluto_setup/radiat.cpp so that the linear theory
and the simulation analysis use exactly the same microphysics:

    cooling rate per volume  = n_e n_H Lambda(T)     [erg cm^-3 s^-1]
    n_H = rho X / (m_u A_H),  n_e = n_H (1 + 2 f_He + A_Z f_Z / 2)   (fully ionised)
    T   = mu m_u p / (rho k_B)

Cooling function
----------------
Lambda(T, Z) = Lambda_HHe(T) + (Z/Zsun) Lambda_metal,sun(T)

* Lambda_HHe: primordial H+He collisional-ionisation-equilibrium cooling computed
  from the Katz, Weinberg & Hernquist (1996) rates (collisional excitation,
  ionisation, recombination, dielectronic recombination, bremsstrahlung).
* Lambda_metal,sun: a smooth (PCHIP in log-log) curve approximating the metal
  contribution to the Sutherland & Dopita (1993) solar CIE curve.

IMPORTANT: the metal part is an approximation typed in from the published
curves (the S&D93 tables could not be downloaded in the environment where this
was written).  It reproduces the characteristic features (C/O peak at
~1-2e5 K, Fe/Ne shoulder near 1e6 K, decline to the bremsstrahlung floor) to
~0.1-0.2 dex.  To use the real tables, write them with `write_pluto_table`
from your own (T, Lambda) arrays, or call `load_table(path)`.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator

# ---------------------------------------------------------------------------
# constants (cgs), identical to gPLUTO's gpluto.hpp
# ---------------------------------------------------------------------------
kB = 1.3806505e-16
amu = 1.66053886e-24
mp = 1.67262171e-24
pc = 3.0856775807e18
kpc = 1.0e3 * pc
Myr = 3.15576e13
A_H, A_He, A_Z = 1.008, 4.004, 30.0
X_H, Y_He = 0.7110, 0.2741           # gPLUTO defaults for COOLING = TABULATED
Z_m = 1.0 - X_H - Y_He
f_He = Y_He / A_He * A_H / X_H
f_Z = Z_m / A_Z * A_H / X_H
MU = (A_H + f_He * A_He + f_Z * A_Z) / (2.0 + 3.0 * f_He + f_Z * (1.0 + 0.5 * A_Z))
NE_OVER_NH = 1.0 + 2.0 * f_He + 0.5 * A_Z * f_Z
GAMMA = 5.0 / 3.0

# code units used in pluto_setup/init.cpp
UNIT_RHO = mp
UNIT_L = kpc
UNIT_V = 1.0e7
UNIT_T = UNIT_L / UNIT_V            # 9.78 Myr


def nH_of_rho(rho_cgs):
    return rho_cgs * X_H / (amu * A_H)


def nenH_of_rho(rho_cgs):
    nH = nH_of_rho(rho_cgs)
    return NE_OVER_NH * nH * nH


def T_of(rho_cgs, p_cgs):
    return MU * amu * p_cgs / (rho_cgs * kB)


def p_of(rho_cgs, T):
    return rho_cgs * kB * T / (MU * amu)


def cs_of(T):
    return np.sqrt(GAMMA * kB * T / (MU * amu))


# ---------------------------------------------------------------------------
# primordial CIE cooling (Katz, Weinberg & Hernquist 1996)
# ---------------------------------------------------------------------------
def lambda_HHe(T):
    T = np.asarray(T, dtype=float)
    T3, T5, T6 = T / 1e3, T / 1e5, T / 1e6
    s = 1.0 / (1.0 + np.sqrt(T5))
    # rates
    aHp = 8.4e-11 / np.sqrt(T) * T3**-0.2 / (1 + T6**0.7)
    aHep = 1.5e-10 * T**-0.6353
    ad = 1.9e-3 * T**-1.5 * np.exp(-470000.0 / T) * (1 + 0.3 * np.exp(-94000.0 / T))
    aHepp = 3.36e-10 / np.sqrt(T) * T3**-0.2 / (1 + T6**0.7)
    GH0 = 5.85e-11 * np.sqrt(T) * np.exp(-157809.1 / T) * s
    GHe0 = 2.38e-11 * np.sqrt(T) * np.exp(-285335.4 / T) * s
    GHep = 5.68e-12 * np.sqrt(T) * np.exp(-631515.0 / T) * s
    # ionisation equilibrium (per n_H)
    xH0 = aHp / (aHp + GH0)
    xHp = 1 - xH0
    r1 = GHe0 / (aHep + ad)          # He+/He0
    r2 = GHep / aHepp                # He++/He+
    y0 = 1.0 / (1 + r1 + r1 * r2)
    y1, y2 = r1 * y0, r1 * r2 * y0
    nHe0, nHep, nHepp = f_He * y0, f_He * y1, f_He * y2
    ne = xHp + nHep + 2 * nHepp
    # cooling, per n_H^2
    c = 0.0
    c += 7.50e-19 * s * np.exp(-118348.0 / T) * ne * xH0
    c += 5.54e-17 * T**-0.397 * s * np.exp(-473638.0 / T) * ne * nHep
    c += 1.27e-21 * np.sqrt(T) * s * np.exp(-157809.1 / T) * ne * xH0
    c += 9.38e-22 * np.sqrt(T) * s * np.exp(-285335.4 / T) * ne * nHe0
    c += 4.95e-22 * np.sqrt(T) * s * np.exp(-631515.0 / T) * ne * nHep
    c += 8.70e-27 * np.sqrt(T) * T3**-0.2 / (1 + T6**0.7) * ne * xHp
    c += 1.55e-26 * T**0.3647 * ne * nHep
    c += 3.48e-26 * np.sqrt(T) * T3**-0.2 / (1 + T6**0.7) * ne * nHepp
    c += 1.24e-13 * T**-1.5 * np.exp(-470000.0 / T) * (1 + 0.3 * np.exp(-94000.0 / T)) * ne * nHep
    gff = 1.1 + 0.34 * np.exp(-(5.5 - np.log10(T)) ** 2 / 3.0)
    c += 1.42e-27 * gff * np.sqrt(T) * (xHp + nHep + 4 * nHepp) * ne
    # normalise to n_e(fully ionised) n_H, the convention used in radiat.cpp
    return c / NE_OVER_NH


# metal contribution for solar abundances, log10 Lambda / (erg cm^3 s^-1),
# normalised to n_e n_H (approximation to S&D93 solar CIE minus primordial)
_LOGT_MET = np.array([4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0,
                      6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.5])
_LOGL_MET = np.array([-24.0, -22.7, -22.2, -21.9, -21.55, -21.25, -21.2, -21.22,
                      -21.42, -21.68, -21.8, -21.85, -22.05, -22.3, -22.55,
                      -22.75, -22.95, -23.15, -23.35, -23.5, -23.6, -23.9])
_met_interp = PchipInterpolator(_LOGT_MET, _LOGL_MET, extrapolate=True)


def lambda_metal_sun(T):
    return 10.0 ** _met_interp(np.log10(T))


def lambda_total(T, Z):
    """Lambda(T, Z) in erg cm^3 s^-1 (n_e n_H normalisation); Z in solar units."""
    return lambda_HHe(T) + Z * lambda_metal_sun(T)


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------
def make_table(Z, logTmin=3.8, logTmax=8.5, dlogT=0.005):
    T = 10.0 ** np.arange(logTmin, logTmax + 0.5 * dlogT, dlogT)
    return T, lambda_total(T, Z)


def write_pluto_table(path, T, L):
    with open(path, "w") as f:
        for t, l in zip(T, L):
            f.write(f"{t:.8e}  {l:.8e}\n")


def load_table(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1]


class CoolingTable:
    """Lambda(T) with the same linear-in-T interpolation that PLUTO uses."""

    def __init__(self, T, L):
        self.T = np.asarray(T, float)
        self.L = np.asarray(L, float)
        self.lnT = np.log(self.T)
        self.lnL = np.log(self.L)

    @classmethod
    def for_Z(cls, Z):
        return cls(*make_table(Z))

    def __call__(self, T):
        return np.interp(T, self.T, self.L, left=0.0, right=0.0)

    def dlnL_dlnT(self, T):
        """Logarithmic slope Lambda_T = dlnLambda/dlnT (centred differences)."""
        s = np.gradient(self.lnL, self.lnT)
        return np.interp(np.log(T), self.lnT, s)

    # ---- characteristic times (cgs) ----
    def tcool_isochoric(self, n_or_rho, T, rho=True):
        """e/(n_e n_H Lambda), e = p/(gamma-1)."""
        rho_c = n_or_rho
        p = p_of(rho_c, T)
        return p / (GAMMA - 1) / (nenH_of_rho(rho_c) * self(T))
