"""
Cooling models shared by the linear theory and the gPLUTO driver.

Every model returns the NET radiative loss per unit volume L(rho, T)
[erg cm^-3 s^-1] (cooling minus heating) and its logarithmic derivatives, with
exactly the interpolation used in pluto_setup/radiat.cpp.

  CIEFloor(Z)            -- original setup: approximate S&D93-like CIE curve
                            (cooling.py), n_e n_H normalisation, cut off below
                            T_floor = 1e4 K.  COOL_MODE = 0 in PLUTO.
  Cloudy2D("noUVB", Z)   -- Cloudy CIE (Grackle tables), n_H^2 normalisation,
                            same 1e4 K floor.  COOL_MODE = 1.
  Cloudy2D("hm2012", Z)  -- Cloudy photoionisation equilibrium in the
                            Haardt & Madau (2012) z = 0 UV background:
                            photoheating included, no floor; the gas settles at
                            the density-dependent T_eq(n).  COOL_MODE = 1.

Note the photoionised model breaks the n -> a n, (x, t) -> (x, t)/a scaling of
pure collisional cooling: results now depend on the density itself.
"""
import os

import numpy as np

import cooling as C

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
LOG_ZERO = -99.0


def cutoff(T, Tfloor):
    """PLUTO's smooth cooling cutoff f = 1 - sech((T/Tf)^12) and dln f/dln T."""
    y = (T / Tfloor) ** 12
    if y > 700:
        return 1.0, 0.0
    sech = 1.0 / np.cosh(y)
    f = 1.0 - sech
    return f, 12.0 * y * sech * np.tanh(y) / max(f, 1e-300)


class CIEFloor:
    """The original model: Lambda(T) n_e n_H, linear-in-T table, 1e4 K floor."""
    cool_mode = 0

    def __init__(self, Z, Tfloor=1.0e4):
        self.Z, self.Tfloor = Z, Tfloor
        self.tab = C.CoolingTable.for_Z(Z)
        self.name = f"CIE (approx.), Z={Z}"
        self.label = "cie"

    def derivs(self, rho, T):
        f, dlnf = cutoff(T, self.Tfloor)
        L = C.nenH_of_rho(rho) * self.tab(T) * f
        return L, 2.0 * L, L * (self.tab.dlnL_dlnT(T) + dlnf)

    def net(self, rho, T):
        return C.nenH_of_rho(rho) * self.tab(T) * cutoff(T, self.Tfloor)[0]

    def T_eq(self, rho):
        return self.Tfloor

    def write_pluto(self, d):
        C.write_pluto_table(os.path.join(d, "cooltable.dat"), self.tab.T, self.tab.L)
        return dict(cool_mode=0, T_floor=self.Tfloor)


class Cloudy2D:
    """Grackle/Cloudy tables in (log n_H, log T); rates are Lambda/n_H^2."""
    cool_mode = 1

    def __init__(self, kind, Z, Tfloor=None):
        d = np.load(os.path.join(DATA, f"cloudy_{kind}_z0.npz"))
        self.kind, self.Z = kind, Z
        self.logn, self.logT = d["logn"], d["logT"]
        cool = d["prim_cool"] + Z * d["metal_cool"]
        heat = d["prim_heat"] + Z * d["metal_heat"]
        with np.errstate(divide="ignore"):
            self.lc = np.where(cool > 0, np.log10(np.maximum(cool, 1e-300)), LOG_ZERO)
            self.lh = np.where(heat > 0, np.log10(np.maximum(heat, 1e-300)), LOG_ZERO)
        # default: floor at 1e4 K without a UV background, none with it
        if Tfloor is None:
            Tfloor = 1.0e4 if kind == "noUVB" else 50.0
        self.Tfloor = Tfloor
        self.name = {"noUVB": "CIE (Cloudy)", "hm2012": "UV bkg (HM12, z=0)"}[kind] + f", Z={Z}"
        self.label = {"noUVB": "ciecl", "hm2012": "pie"}[kind]

    # -- bilinear interpolation in (log n, log T), identical to radiat.cpp --
    def _cell(self, x, y):
        i = int(np.clip(np.searchsorted(self.logn, x, side="right") - 1, 0, len(self.logn) - 2))
        j = int(np.clip(np.searchsorted(self.logT, y, side="right") - 1, 0, len(self.logT) - 2))
        dx = self.logn[i + 1] - self.logn[i]
        dy = self.logT[j + 1] - self.logT[j]
        tx, ty = (x - self.logn[i]) / dx, (y - self.logT[j]) / dy
        inx, iny = 0.0 <= tx <= 1.0, 0.0 <= ty <= 1.0      # clamped -> no slope
        return i, j, min(max(tx, 0.0), 1.0), min(max(ty, 0.0), 1.0), dx, dy, inx, iny

    def _interp_d(self, tab, x, y):
        """value and partial derivatives (per dex) of the bilinear interpolant."""
        i, j, tx, ty, dx, dy, inx, iny = self._cell(x, y)
        f00, f10 = tab[i, j], tab[i + 1, j]
        f01, f11 = tab[i, j + 1], tab[i + 1, j + 1]
        v = (1 - tx) * (1 - ty) * f00 + tx * (1 - ty) * f10 + (1 - tx) * ty * f01 + tx * ty * f11
        fx = ((1 - ty) * (f10 - f00) + ty * (f11 - f01)) / dx if inx else 0.0
        fy = ((1 - tx) * (f01 - f00) + tx * (f11 - f10)) / dy if iny else 0.0
        return v, fx, fy

    def rates(self, nH, T):
        x, y = np.log10(nH), np.log10(T)
        return 10.0 ** self._interp_d(self.lc, x, y)[0], 10.0 ** self._interp_d(self.lh, x, y)[0]

    def _L(self, lognH, logT):
        c = 10.0 ** self._interp_d(self.lc, lognH, logT)[0]
        h = 10.0 ** self._interp_d(self.lh, lognH, logT)[0]
        f = cutoff(10.0 ** logT, self.Tfloor)[0]
        return (10.0 ** lognH) ** 2 * (c - h) * f

    def net(self, rho, T):
        return self._L(np.log10(C.nH_of_rho(rho)), np.log10(T))

    def derivs(self, rho, T):
        """L, dL/dln(rho)|_T, dL/dln(T)|_rho -- exact slopes of the interpolant
        PLUTO uses (piecewise constant within a table cell)."""
        nH = C.nH_of_rho(rho)
        x, y = np.log10(nH), np.log10(T)
        fc, fcx, fcy = self._interp_d(self.lc, x, y)
        fh, fhx, fhy = self._interp_d(self.lh, x, y)
        c, h = 10.0 ** fc, 10.0 ** fh
        f, dlnf = cutoff(T, self.Tfloor)
        n2 = nH * nH * f
        L = n2 * (c - h)
        Lr = n2 * (2.0 * (c - h) + c * fcx - h * fhx)
        LT = n2 * (c * fcy - h * fhy) + L * dlnf
        return L, Lr, LT

    def T_eq(self, rho):
        """Highest temperature below 10^5.5 K at which net cooling = 0."""
        if self.kind == "noUVB":
            return self.Tfloor
        x = round(float(np.log10(C.nH_of_rho(rho))), 3)
        cache = self.__dict__.setdefault("_teq_cache", {})
        if x not in cache:
            cache[x] = self._T_eq(x)
        return cache[x]

    def _T_eq(self, x):
        ys = np.linspace(2.5, 5.5, 601)
        Ls = np.array([self._L(x, y) for y in ys])
        idx = np.where((Ls[:-1] <= 0) & (Ls[1:] > 0))[0]
        if len(idx) == 0:
            return 10.0 ** ys[0]
        k = idx[-1]
        y0 = ys[k] - Ls[k] * (ys[k + 1] - ys[k]) / (Ls[k + 1] - Ls[k])
        return 10.0 ** y0

    def write_pluto(self, d):
        # PLUTO still reads cooltable.dat at start-up; give it the 1D CIE table
        T1, L1 = C.make_table(self.Z)
        C.write_pluto_table(os.path.join(d, "cooltable.dat"), T1, L1)
        with open(os.path.join(d, "cooltable2d.dat"), "w") as f:
            f.write(f"{len(self.logn)} {len(self.logT)}\n")
            f.write(" ".join(f"{v:.8f}" for v in self.logn) + "\n")
            f.write(" ".join(f"{v:.8f}" for v in self.logT) + "\n")
            for tab in (self.lc, self.lh):
                for row in tab:
                    f.write(" ".join(f"{v:.8f}" for v in row) + "\n")
        return dict(cool_mode=1, T_floor=self.Tfloor)


def make_model(cooling, Z):
    if cooling == "cie":
        return CIEFloor(Z)
    if cooling == "ciecl":
        return Cloudy2D("noUVB", Z)
    if cooling == "pie":
        return Cloudy2D("hm2012", Z)
    raise ValueError(cooling)
