"""
Figures comparing the original collisional cooling (CIE + 1e4 K floor) with a
photoionising UV background (Cloudy, Haardt & Madau 2012, z = 0).

    python make_figures_uvb.py
"""
import glob
import os

import numpy as np

import cooling as C
import cooling_models as CM
import linear_theory as LT
import plutorun as P
from metrics import t_cool as tc
from make_figures import (FIG, RES, S1, S2, S3, S4, S5, INK2, GRID, plt, save,
                          cloud_metrics)

NCOL = {1e-4: S1, 1e-3: S2, 1e-2: S3}
CIECOL = "#52514e"


def runs(sweep, pattern="*.npz"):
    return [P.load(f) for f in sorted(glob.glob(os.path.join(RES, "runs", sweep, pattern)))]


def fig_cooling_models():
    T = np.logspace(3.5, 7.5, 600)
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9))
    Z = 0.3
    old = CM.CIEFloor(Z)
    cl = CM.Cloudy2D("noUVB", Z, Tfloor=1.0)
    pie = CM.Cloudy2D("hm2012", Z)
    for n, c in NCOL.items():
        rho = n * C.mp
        nH = C.nH_of_rho(rho)
        net = np.array([pie.net(rho, t) for t in T]) / nH**2
        ax[0].loglog(T, np.where(net > 0, net, np.nan), color=c, label=f"UV bkg, n = {n:g}")
        ax[0].loglog(T, np.where(net < 0, -net, np.nan), color=c, ls=":", lw=1.2)
        Teq = pie.T_eq(rho)
        ax[0].axvline(Teq, color=c, lw=0.8, ls="--")
        tc = np.array([C.p_of(rho, t) / (C.GAMMA - 1) / abs(pie.net(rho, t)) for t in T]) / C.Myr
        ax[1].loglog(T, tc * n / 1e-3, color=c, label=f"UV bkg, n = {n:g}")
    ce = np.array([cl.net(1e-3 * C.mp, t) for t in T]) / C.nH_of_rho(1e-3 * C.mp) ** 2
    ax[0].loglog(T, ce, color=CIECOL, lw=1.5, label="CIE (Cloudy)")
    oa = np.array([C.nenH_of_rho(1e-3 * C.mp) * old.tab(t) for t in T]) / C.nH_of_rho(1e-3 * C.mp) ** 2
    ax[0].loglog(T, oa, color=CIECOL, lw=1.0, ls="--", label="CIE (approx., original runs)")
    for lab, mod in [("CIE (Cloudy)", cl)]:
        tc = np.array([C.p_of(1e-3 * C.mp, t) / (C.GAMMA - 1) / mod.net(1e-3 * C.mp, t) for t in T]) / C.Myr
        ax[1].loglog(T, tc, color=CIECOL, lw=1.5, label=lab)
    ax[0].set_ylim(1e-24, 1e-21)
    ax[0].set_xlabel("T [K]")
    ax[0].set_ylabel(r"net $\Lambda / n_H^2$ [erg cm$^3$ s$^{-1}$]")
    ax[0].set_title("Net cooling, Z = 0.3 (dotted: net heating)")
    ax[0].legend(fontsize=7.5, loc="upper right")
    ax[1].set_xlabel("T [K]")
    ax[1].set_ylabel(r"$t_{cool}\times(n/10^{-3})$ [Myr]")
    ax[1].set_title("Cooling time, scaled to n = 10$^{-3}$ cm$^{-3}$")
    ax[1].set_ylim(1, 1e4)
    save(fig, "fig12_cooling_models.png")


def fig_validation_uvb():
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9))
    for c, x in zip([S1, S2, S3], [0.1, 1.0, 10.0]):
        f = os.path.join(RES, "validation_pie", f"valpie_x{x:g}.npz")
        if not os.path.exists(f):
            continue
        d = P.load(f)
        m = CM.make_model("pie", 0.3)
        o = LT.integrate(1e-3, 1e6, d["params"]["lam"], model=m, background="isochoric",
                         amp=1e-4, Tstop=0, t_end=d["t"][-1] * C.Myr * 1.001, t_out=d["t"])
        n = len(o["t"])
        ax[0].plot(o["t"], o["s"] / 1e-4, color=c, lw=1.5, label=fr"$\lambda/c_s t_{{cool}}$ = {x:g}")
        ax[0].plot(d["t"][:n:5], d["amp_mode"][:n:5] / 1e-4, "o", ms=4, mfc="none", color=c)
        ax[1].semilogy(o["t"], o["T"], color=c, lw=1.5)
        ax[1].semilogy(d["t"][:n:5], d["Tmean"][:n:5], "o", ms=4, mfc="none", color=c)
    ax[0].set_yscale("symlog", linthresh=1)
    ax[0].set_xlabel("t [Myr]")
    ax[0].set_ylabel(r"$\delta\rho/\rho / A_0$ (seeded mode)")
    ax[0].set_title("UV bkg, n = 10$^{-3}$: gPLUTO (circles) vs linear theory")
    ax[0].legend(fontsize=8)
    ax[1].set_xlabel("t [Myr]")
    ax[1].set_ylabel("mean T [K]")
    ax[1].set_title("Gas settles at T$_{eq}$ ≈ 1.3×10$^4$ K instead of a floor")
    save(fig, "fig13_validation_uvb.png")


def fig_linear_uvb():
    f = os.path.join(RES, "linear_maps_uvb.npz")
    if not os.path.exists(f):
        return
    m = dict(np.load(f))
    mo = dict(np.load(os.path.join(RES, "linear_maps.npz")))
    for T0, tag in [(1e6, ""), (5e5, "_T5e5")]:
        fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0), sharey=True)
        for a, bg in [(ax[0], "isochoric"), (ax[1], "isobaric")]:
            xmax = 1.0 if bg == "isobaric" else 1e9
            sel = (mo["bg"] == bg) & (mo["Z"] == 0.3) & (mo["T0"] == T0) & (mo["x"] <= xmax)
            o = np.argsort(mo["x"][sel])
            a.loglog(mo["x"][sel][o], mo["G_max_cool"][sel][o], color=CIECOL, lw=1.0, ls="--",
                     label="CIE approx. + floor (original)")
            for cool, n, c, lab in [("ciecl", 1e-2, CIECOL, "CIE (Cloudy) + floor")] + \
                    [("pie", n, NCOL[n], f"UV bkg, n = {n:g}") for n in (1e-4, 1e-3, 1e-2)]:
                sel = ((m["cooling"] == cool) & (m["n0"] == n) & (m["bg"] == bg) & (m["Z"] == 0.3)
                       & (m["T0"] == T0) & (m["x"] <= xmax))
                o = np.argsort(m["x"][sel])
                a.loglog(m["x"][sel][o], m["G_max_cool"][sel][o], color=c, lw=1.6, label=lab)
                a.loglog(m["x"][sel][o], m["G_3tf"][sel][o], color=c, lw=0.9, ls=":")
            a.set_xlabel(r"$\lambda_0/(c_s t_{cool})_0$")
            a.set_title(("Infinite medium" if bg == "isochoric" else r"Small cloud ($L \ll c_s t_{cool}$)")
                        + f", Z = 0.3, T$_0$ = {T0:.0e} K".replace("e+0", "×10$^") + "$")
        ax[0].set_ylabel(r"amplification $|\delta\rho/\rho|/A_0$")
        ax[0].legend(fontsize=7.5, loc="upper right")
        ax[1].text(0.04, 1.2, "solid: max while cooling (to 1.5 T$_{eq}$)\ndotted: max by 3 t$_{cool}$",
                   fontsize=8, color=INK2)
        save(fig, f"fig14_linear_uvb{tag}.png")


def fig_infinite_uvb():
    R = runs("infinite_uvb")
    if not R:
        return
    f = os.path.join(RES, "linear_maps_uvb.npz")
    m = dict(np.load(f)) if os.path.exists(f) else None
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9), sharey=True)
    for a, Z in zip(ax, [0.1, 0.3, 1.0]):
        for cool, n, c, lab in [("ciecl", 1e-2, CIECOL, "CIE (Cloudy)")] + \
                [("pie", n, NCOL[n], f"UV bkg, n = {n:g}") for n in (1e-4, 1e-3, 1e-2)]:
            pts = []
            for d in R:
                p = d["params"]
                if p["cooling"] != cool or p["Z"] != Z or p["rho_bg"] != n or p["amp"] != 1e-3:
                    continue
                tf = tc(d)
                cool_ph = d["t"] <= tf
                pts.append((p["lam"] / d["cs_tcool0"], np.abs(d["amp_mode"][cool_ph]).max() / 1e-3))
            pts = np.array(sorted(pts))
            if len(pts):
                a.loglog(pts[:, 0], pts[:, 1], "o", color=c, ms=5, label=lab)
            if m is not None:
                sel = ((m["cooling"] == cool) & (m["n0"] == n) & (m["bg"] == "isochoric")
                       & (m["Z"] == Z) & (m["T0"] == 1e6))
                o = np.argsort(m["x"][sel])
                a.loglog(m["x"][sel][o], m["G_max_cool"][sel][o], color=c, lw=1, ls="--")
        a.set_xlabel(r"$\lambda_0/(c_s t_{cool})_0$")
        a.set_title(f"Infinite medium, Z = {Z}, T$_0$ = 10$^6$ K")
    ax[0].set_ylabel(r"max $|\delta\rho/\rho|/A_0$ while cooling")
    ax[0].legend(fontsize=8)
    ax[1].text(0.035, 1.15, "points: gPLUTO, A$_0$ = 10$^{-3}$;  dashed: linear theory", fontsize=8, color=INK2)
    save(fig, "fig15_infinite_uvb.png")


def fig_cloud_uvb():
    Rp = [d for d in runs("cloud_uvb") if d["params"]["chi"] == 2.0]
    Rc = {d["params"]["name"].replace("cl_", "clpie_"): d for d in runs("cloud") if d["params"]["chi"] == 2.0}
    if not Rp:
        return
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.1))
    for d in Rp:
        mp = cloud_metrics(d)
        c = NCOL[mp["n"]]
        ok = mp["cells_per_lam_final"] >= 10
        ax[0].loglog(mp["x"], mp["G_cool"], "o", color=c, ms=5.5, mfc=c if ok else "none")
        dc = Rc.get(d["params"]["name"])
        if dc is not None:
            mc = cloud_metrics(dc)
            ax[1].loglog(mc["G_cool"], mp["G_cool"], "o", color=c, ms=5.5, mfc=c if ok else "none")
            ax[2].loglog(mc["G_end"], mp["G_end"], "o", color=c, ms=5.5, mfc=c if ok else "none")
    for a in ax[1:]:
        lim = [0.8, 3e3]
        a.plot(lim, lim, color=INK2, lw=0.8, ls="--")
        a.set_xlabel("CIE + floor (original runs)")
        a.set_ylabel("UV background")
    ax[0].set_xlabel(r"$\lambda_0/(c_s t_{cool})_0$ (cloud gas, UV bkg)")
    ax[0].set_ylabel(r"max $|\delta\rho/\rho|/A_0$ while cooling")
    ax[0].set_title("Cloud runs with a UV background")
    ax[1].set_title("Same clouds: growth while cooling")
    ax[2].set_title("Same clouds: growth by 2.5 t$_{cool}$")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=c, marker="o", ls="none", label=f"n$_{{bg}}$ = {n:g}") for n, c in NCOL.items()]
    h += [Line2D([], [], color=INK2, marker="o", mfc="none", ls="none", label="< 10 cells/λ after collapse")]
    ax[0].legend(handles=h, fontsize=7.5)
    save(fig, "fig16_cloud_uvb.png")


def fig_amplitude_uvb():
    Rp, Rc = runs("amplitude_uvb"), runs("amplitude")
    if not Rp:
        return
    fig, ax = plt.subplots(2, 3, figsize=(13, 6.8), sharey=True)
    amps = [1e-3, 1e-2, 0.1, 0.3]
    cols = [S1, S2, S3, S4]
    for row, (R, lab) in enumerate([(Rc, "CIE + floor"), (Rp, "UV background")]):
        for a, n in zip(ax[row], [1e-4, 1e-3, 1e-2]):
            for c, A in zip(cols, amps):
                for d in R:
                    p = d["params"]
                    if p["rho_bg"] == n and p["amp"] == A:
                        tf = tc(d)
                        k = d["t"] <= 3 * tf
                        a.semilogy(d["t"][k] / tf, np.maximum(d["amp_max"][k], 1e-4), color=c, lw=1.4,
                                   label=f"A$_0$ = {A:g}")
                        Lx = p["L_cloud"] / d["cs_tcool0"]
            a.set_title(f"{lab}: n$_{{bg}}$ = {n:g} (L/c$_s$t$_{{cool}}$ = {Lx:.2g})")
            a.set_ylim(1e-4, 5e2)
            a.axvline(1, color=GRID, lw=1)
            if row == 1:
                a.set_xlabel(r"t / t$_{cool}$")
        ax[row, 0].set_ylabel(r"interior max $\rho/\bar\rho - 1$")
    ax[0, 0].legend(fontsize=7.5)
    save(fig, "fig17_amplitude_uvb.png")


if __name__ == "__main__":
    fig_cooling_models()
    fig_validation_uvb()
    fig_linear_uvb()
    fig_infinite_uvb()
    fig_cloud_uvb()
    fig_amplitude_uvb()
