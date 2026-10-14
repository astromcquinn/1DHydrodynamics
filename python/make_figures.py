"""
Make the figures in ../results/figures from the sweep outputs.

    python make_figures.py
"""
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cooling as C
import linear_theory as LT
import plutorun as P
from metrics import t_cool as tc

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

# categorical slots (fixed order) and text/grid colours
S1, S2, S3, S4, S5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
CAT = [S1, S2, S3, S4, S5, "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
ZCOL = {0.05: S5, 0.1: S1, 0.3: S2, 1.0: S3}

plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "savefig.facecolor": "#fcfcfb", "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 2.0, "font.size": 10, "legend.frameon": False,
    "axes.titlesize": 10.5, "axes.titleweight": "bold",
})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name), dpi=150)
    plt.close(fig)
    print("wrote", name)


def runs(sweep):
    return [P.load(f) for f in sorted(glob.glob(os.path.join(RES, "runs", sweep, "*.npz")))]


# ---------------------------------------------------------------------------
def fig_cooling():
    T = np.logspace(4.0, 8.0, 800)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for Z in [0.05, 0.1, 0.3, 1.0]:
        tab = C.CoolingTable.for_Z(Z)
        ax[0].loglog(T, tab(T), color=ZCOL[Z], label=f"Z = {Z} Z$_\\odot$")
        ax[1].semilogx(T, tab.dlnL_dlnT(T), color=ZCOL[Z])
    ax[0].set_ylim(10**-23.6, 1e-21)
    ax[0].set_xlabel("T [K]")
    ax[0].set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]  ($n_e n_H$ norm.)")
    ax[0].set_title("Cooling function used (approximate S&D93 CIE)")
    ax[0].legend()
    for y, lab in [(2, r"isobaric unstable below ($\Lambda_T<2$), both cases"),
                   (1, r"cooling gas: $\delta p$ grows below ($\Lambda_T<1$)"),
                   (0, r"Das et al. isochoric unstable below ($\Lambda_T<0$)")]:
        ax[1].axhline(y, color=INK2, lw=0.8, ls="--")
        ax[1].text(9.5e7, y + 0.1, lab, color=INK2, fontsize=7.5, ha="right",
                   bbox=dict(fc="#fcfcfb", ec="none", pad=0.5))
    ax[1].set_ylim(-3, 4)
    ax[1].set_xlim(1e4, 1e8)
    ax[1].set_xlabel("T [K]")
    ax[1].set_ylabel(r"$\Lambda_T = d\ln\Lambda/d\ln T$")
    ax[1].set_title("Logarithmic slope sets the growth rates")
    save(fig, "fig1_cooling_function.png")


# ---------------------------------------------------------------------------
def fig_validation():
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for c, lam in zip([S1, S2, S3], [0.5, 5.0, 50.0]):
        f = os.path.join(RES, "validation", f"val_lam{lam:g}.npz")
        if not os.path.exists(f):
            continue
        d = P.load(f)
        o = LT.integrate(0.01, 1e6, lam, Z=0.3, background="isochoric",
                         amp=1e-4, t_out=d["t"])
        n = len(o["t"])
        x = lam / o["cs_tcool0_kpc"]
        ax[0].plot(o["t"], np.abs(o["s"]) / 1e-4, color=c, lw=1.5,
                   label=fr"$\lambda$ = {lam:g} kpc ($\lambda/c_s t_{{cool}}$ = {x:.2g})")
        ax[0].plot(d["t"][:n:4], np.abs(d["amp_mode"][:n:4]) / 1e-4, "o", ms=4,
                   mfc="none", color=c)
        ax[1].plot(o["T"], np.abs(o["s"]) / 1e-4, color=c, lw=1.5)
        ax[1].plot(d["Tmean"][:n:4], np.abs(d["amp_mode"][:n:4]) / 1e-4, "o",
                   ms=4, mfc="none", color=c)
    tab = C.CoolingTable.for_Z(0.3)
    TT = np.logspace(np.log10(3e4), 6, 200)
    ax[1].plot(TT, LT.timeshift_amplification(TT, 1e6, tab, "isochoric"),
               color=INK2, lw=1, ls="--", label="short-$\\lambda$ limit (time shift)")
    for a in ax:
        a.set_yscale("log")
        a.set_ylabel(r"$|\delta\rho/\rho| / A_0$")
    ax[0].set_xlabel("t [Myr]")
    ax[1].set_xscale("log")
    ax[1].invert_xaxis()
    ax[1].set_xlabel("mean T [K] (cooling →)")
    ax[0].legend(fontsize=8)
    ax[1].legend(fontsize=8)
    ax[0].set_title("gPLUTO (circles) vs linear theory (lines)")
    ax[1].set_title("Same, against the background temperature")
    save(fig, "fig2_validation.png")


# ---------------------------------------------------------------------------
def fig_linear_maps():
    f = os.path.join(RES, "linear_maps.npz")
    if not os.path.exists(f):
        return
    m = dict(np.load(f))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.0), sharey=True)
    for a, bg, title in [(ax[0], "isochoric", "Infinite cooling medium"),
                         (ax[1], "isobaric", "Small cloud ($L \\ll c_s t_{cool}$)")]:
        for Z in [0.05, 0.1, 0.3, 1.0]:
            sel = (m["bg"] == bg) & (m["Z"] == Z) & (m["T0"] == 1e6)
            if bg == "isobaric":          # lambda < L < c_s t_cool only
                sel &= m["x"] <= 1.0
            o = np.argsort(m["x"][sel])
            a.loglog(m["x"][sel][o], m["G_max_cool"][sel][o], color=ZCOL[Z],
                     label=f"Z = {Z}")
            a.loglog(m["x"][sel][o], m["G_3tf"][sel][o], color=ZCOL[Z], lw=1, ls=":")
        sel = (m["bg"] == bg) & (m["Z"] == 0.3) & (m["T0"] == 1e6)
        if bg == "isobaric":
            sel &= m["x"] <= 1.0
        o = np.argsort(m["x"][sel])
        a.loglog(m["x"][sel][o], np.where(m["x"][sel][o] < 1, m["das_isobaric"][sel][o],
                                          m["das_isochoric"][sel][o]),
                 color=INK2, lw=1, ls="--", label="heated medium (Das et al.), same elapsed time")
        a.axvline(1, color=GRID, lw=1)
        a.set_xlabel(r"$\lambda_0 / (c_s t_{cool})_0$")
        a.set_title(title + r", $T_0=10^6$ K")
    ax[0].set_ylabel(r"amplification $|\delta\rho/\rho|/A_0$")
    ax[0].legend(fontsize=8, loc="upper right")
    ax[1].text(0.04, 3e3, "solid: max while cooling (to 1.5 T$_{floor}$)\ndotted: max by $3\\,t_{cool}$ (incl. ballistic phase)",
               fontsize=8, color=INK2, va="top")
    for a in ax:
        a.set_xlim(0.02, 400)
    save(fig, "fig3_linear_amplification.png")

    # T0 dependence for Z = 0.3
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.0), sharey=True)
    for a, bg in [(ax[0], "isochoric"), (ax[1], "isobaric")]:
        for c, T0 in zip([S1, S2, S3], [1e6, 5e5, 2.5e5]):
            sel = (m["bg"] == bg) & (m["Z"] == 0.3) & (m["T0"] == T0)
            if bg == "isobaric":
                sel &= m["x"] <= 1.0
            o = np.argsort(m["x"][sel])
            a.loglog(m["x"][sel][o], m["G_max_cool"][sel][o], color=c,
                     label=f"T$_0$ = {T0:.2g} K")
        a.set_xlabel(r"$\lambda_0 / (c_s t_{cool})_0$")
        a.set_title(f"{bg} background, Z = 0.3")
    ax[0].set_ylabel(r"max amplification during cooling")
    ax[0].legend(fontsize=8)
    save(fig, "fig3b_linear_amplification_T0.png")

    # minimum amplitude to go nonlinear, in physical units for several n
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.0), sharey=True)
    for a, bg in [(ax[0], "isochoric"), (ax[1], "isobaric")]:
        sel = (m["bg"] == bg) & (m["Z"] == 0.3) & (m["T0"] == 1e6)
        if bg == "isobaric":
            sel &= m["x"] <= 1.0
        o = np.argsort(m["x"][sel])
        x, G = m["x"][sel][o], m["G_3tf"][sel][o]
        for c, n in zip([S1, S2, S3], [1e-2, 1e-3, 1e-4]):
            # x is measured with the initial (cloud) gas; for n=0.01 cs*tcool:
            cst = sweep_cs_tcool(n, 1e6, 0.3)
            a.loglog(x * cst, np.minimum(1.0 / G, 1.0), color=c, label=f"n = {n:g} cm$^{{-3}}$")
        a.set_xlabel(r"initial wavelength $\lambda_0$ [kpc]")
        a.set_title(f"{bg} background, Z = 0.3, T$_0$ = 10$^6$ K")
    ax[0].set_ylabel(r"min. $\delta\rho/\rho$ that reaches order unity")
    ax[0].legend(fontsize=8)
    save(fig, "fig4_min_amplitude_nonlinear.png")


def sweep_cs_tcool(n, T, Z):
    tab = C.CoolingTable.for_Z(Z)
    rho = n * C.mp
    tc = C.p_of(rho, T) / (C.GAMMA - 1) / (C.nenH_of_rho(rho) * tab(T))
    return C.cs_of(T) * tc / C.kpc


# ---------------------------------------------------------------------------
def fig_infinite_sims():
    R = runs("infinite")
    if not R:
        return
    f = os.path.join(RES, "linear_maps.npz")
    m = dict(np.load(f)) if os.path.exists(f) else None
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9), sharey=True)
    for a, Z in zip(ax, [0.1, 0.3, 1.0]):
        for c, A in zip([S1, S2, S3], [1e-3, 1e-2, 1e-1]):
            pts = []
            for d in R:
                p = d["params"]
                if p["Z"] != Z or p["amp"] != A:
                    continue
                x = p["lam"] / d["cs_tcool0"]
                tf = tc(d)
                cool = d["t"] <= tf
                pts.append((x, np.abs(d["amp_mode"][cool]).max() / A,
                            np.abs(d["amp_mode"]).max() / A,
                            (d["amp_max"].max() + 1) / (1 + A)))
            pts = np.array(sorted(pts))
            if len(pts):
                a.loglog(pts[:, 0], pts[:, 1], "o-", color=c, ms=5, lw=1.2,
                         label=f"sim, A$_0$ = {A:g}")
        if m is not None:
            sel = (m["bg"] == "isochoric") & (m["Z"] == Z) & (m["T0"] == 1e6)
            o = np.argsort(m["x"][sel])
            a.loglog(m["x"][sel][o], m["G_max_cool"][sel][o], color=INK2, lw=1,
                     ls="--", label="linear theory")
        a.set_title(f"Infinite medium, Z = {Z}, n = 0.01, T$_0$ = 10$^6$ K")
        a.set_xlabel(r"$\lambda_0/(c_s t_{cool})_0$")
    ax[0].set_ylabel(r"max $|\delta\rho/\rho|/A_0$ while cooling (mode)")
    ax[0].legend(fontsize=8)
    save(fig, "fig5_infinite_sims.png")


# ---------------------------------------------------------------------------
def cloud_metrics(d):
    p = d["params"]
    tf = tc(d)
    A = p["amp"]
    t = d["t"]
    cool = t <= tf
    dx = p["box"] / p["nx"]
    compression = d["nmean"].max() / d["nmean"][0]
    return dict(
        Lx=p["L_cloud"] / d["cs_tcool0"], x=p["lam"] / d["cs_tcool0"],
        G_cool=np.abs(d["amp_mode"][cool]).max() / A,
        G_end=np.abs(d["amp_mode"][t <= 2.5 * tf]).max() / A,
        student=d["nmax_over_nmean"].max(),
        inner_max=d["amp_max"].max(), Z=p["Z"], n=p["rho_bg"], L=p["L_cloud"],
        nw=p["L_cloud"] / p["lam"], chi=p["chi"],
        cells_per_lam_final=p["lam"] / dx / compression, compression=compression)


def fig_cloud():
    R = runs("cloud")
    if not R:
        return
    M = [cloud_metrics(d) for d in R]
    f = os.path.join(RES, "linear_maps.npz")
    m = dict(np.load(f)) if os.path.exists(f) else None
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.1))
    bins = [(0, 1.5, S1, "L/c$_s$t$_{cool}$ < 1.5"), (1.5, 10, S2, "1.5–10"),
            (10, 100, S3, "10–100"), (100, 1e9, S4, "> 100")]

    def col(Lx):
        for lo, hi, c, _ in bins:
            if lo <= Lx < hi:
                return c

    for mm in M:
        ok = mm["cells_per_lam_final"] >= 10
        c = col(mm["Lx"])
        ax[0].loglog(mm["x"], mm["G_cool"], "o", color=c, ms=5.5,
                     mfc=c if ok else "none")
        ax[1].loglog(mm["x"], mm["G_end"], "o", color=c, ms=5.5,
                     mfc=c if ok else "none")
        ax[2].loglog(max(mm["inner_max"], 1e-2), mm["student"], "o", color=c, ms=5.5,
                     mfc=c if ok else "none")
    if m is not None:
        for bg, ls, lab in [("isochoric", "--", "linear theory: infinite medium"),
                            ("isobaric", ":", "linear theory: small isobaric cloud")]:
            sel = (m["bg"] == bg) & (m["Z"] == 0.3) & (m["T0"] == 5e5)
            if bg == "isobaric":
                sel &= m["x"] <= 1.0
            o = np.argsort(m["x"][sel])
            ax[0].loglog(m["x"][sel][o], m["G_max_cool"][sel][o], color=INK2, lw=1.2, ls=ls, label=lab)
            ax[1].loglog(m["x"][sel][o], m["G_3tf"][sel][o], color=INK2, lw=1.2, ls=ls)
    for a in ax[:2]:
        a.set_xlabel(r"$\lambda_0 / (c_s t_{cool})_0$  (cloud gas)")
        a.set_ylim(0.8, 3e3)
    ax[0].set_ylabel(r"max $|\delta\rho/\rho|/A_0$ of seeded mode")
    ax[0].set_title("Cloud runs: growth while cooling (to 1.5 T$_{floor}$)")
    ax[1].set_title("Including post-cooling phase (to 2.5 $t_{cool}$)")
    ax[2].set_xlabel(r"interior max $\rho/\bar\rho - 1$ (edges excluded)")
    ax[2].set_ylabel(r"whole-cloud max$(\rho)/\bar\rho - 1$ (incl. edges)")
    lim = [1e-2, 1e2]
    ax[2].plot(lim, lim, color=INK2, lw=0.8, ls="--")
    ax[2].set_title("A whole-cloud maximum is set by the edges")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=c, marker="o", ls="none", label=lab) for _, _, c, lab in bins]
    h += [Line2D([], [], color=INK2, marker="o", mfc="none", ls="none",
                 label="< 10 cells/λ after collapse")]
    leg1 = ax[0].legend(handles=h, fontsize=7.5, loc="upper right")
    ax[0].add_artist(leg1)
    ax[0].legend(fontsize=7.5, loc="lower left")
    save(fig, "fig6_cloud_sweep.png")


def fig_isochoric():
    """Growth in the isochoric regime (x > 1.5, L/c_s t_cool > 1.5)."""
    R = runs("cloud")
    if not R:
        return
    pts = []
    for d in R:
        m = cloud_metrics(d)
        if m["Lx"] < 1.5 or m["x"] < 1.5:
            continue
        tf = tc(d)
        cool = d["t"] <= 0.95 * tf
        pts.append((m["x"], m["G_cool"] - 1, m["G_end"] - 1,
                    np.abs(d["amp_mode_p"][cool]).max() / m_amp(d)))
    pts = np.array(sorted(pts))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
    ax[0].loglog(pts[:, 0], pts[:, 3], "o", color=S2, ms=5, label=r"$|\delta T/T|$ = $|\delta p/p|$, max before end of cooling")
    ax[0].loglog(pts[:, 0], np.maximum(pts[:, 1], 1e-3) + 1, "o", color=S1, ms=5, mfc="none",
                 label=r"$|\delta\rho/\rho|$, max while cooling")
    ax[0].axhline(1, color=INK2, lw=0.8, ls="--")
    ax[0].set_xlabel(r"$\lambda_0/(c_s t_{cool})_0$")
    ax[0].set_ylabel(r"amplitude / initial $\delta\rho/\rho$")
    ax[0].set_title("While cooling: T and p grow, density does not")
    ax[0].legend(fontsize=8)
    xx = np.logspace(np.log10(1.5), np.log10(500), 50)
    ax[1].loglog(pts[:, 0], pts[:, 2], "o", color=S1, ms=5, label="cloud runs, by 2.5 t$_{cool}$")
    ax[1].loglog(xx, 25 / xx**2, color=INK2, lw=1, ls="--", label=r"$25\,(c_s t_{cool}/\lambda)^2$")
    ax[1].set_xlabel(r"$\lambda_0/(c_s t_{cool})_0$")
    ax[1].set_ylabel(r"$G - 1$ (density growth)")
    ax[1].set_title("Afterwards: density grows from the pressure imprint")
    ax[1].legend(fontsize=8)
    save(fig, "fig11_isochoric_regime.png")


def m_amp(d):
    return d["params"]["amp"]


def fig_profiles():
    picks = [("cloud", "cl_Z0p3_n0p01_L100_nw4p5_chi2", "Large cloud: L/(c$_s$t$_{cool}$) ≈ 190 (isochoric)"),
             ("cloud", "cl_Z0p3_n0p0001_L30_nw4p5_chi2", "Small cloud: L/(c$_s$t$_{cool}$) ≈ 0.6 (≈isobaric)")]
    fig, ax = plt.subplots(2, 3, figsize=(13, 6.6))
    for row, (sw, name, title) in enumerate(picks):
        f = os.path.join(RES, "runs", sw, name + ".npz")
        if not os.path.exists(f):
            continue
        d = P.load(f)
        prof, tp, x = d["profiles"], d["profile_t"], d["x"]
        L = d["params"]["L_cloud"]
        cmap = plt.get_cmap("viridis")
        for i in range(len(tp)):
            c = cmap(i / max(len(tp) - 1, 1))
            rho, prs, T = prof[i][0], prof[i][1], prof[i][2]
            ax[row, 0].semilogy(x / L, rho, color=c, lw=1)
            ax[row, 1].semilogy(x / L, prs / prof[0][1].max(), color=c, lw=1)
            ax[row, 2].semilogy(x / L, T, color=c, lw=1,
                                label=f"{tp[i]:.3g} Myr" if i % 3 == 0 or i == len(tp) - 1 else None)
        for a in ax[row]:
            a.set_xlim(-0.9, 0.9)
            a.set_xlabel("x / L")
        ax[row, 0].set_ylabel(r"n [cm$^{-3}$]")
        ax[row, 1].set_ylabel(r"p / p$_{hot}$")
        ax[row, 2].set_ylabel("T [K]")
        ax[row, 0].set_title(title)
        ax[row, 2].legend(fontsize=7, loc="lower right")
    save(fig, "fig7_profiles.png")


def fig_amplitude():
    R = runs("amplitude")
    if not R:
        return
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9), sharey=True)
    amps = [1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3]
    for a, n in zip(ax, [1e-4, 1e-3, 1e-2]):
        for c, A in zip(CAT, amps):
            for d in R:
                p = d["params"]
                if p["rho_bg"] == n and p["amp"] == A:
                    tf = tc(d)
                    k = d["t"] <= 3 * tf
                    a.semilogy(d["t"][k] / tf, d["amp_rms"][k] / (A / np.sqrt(2)), color=c,
                               lw=1.5, label=f"A$_0$ = {A:g}")
                    Lx = p["L_cloud"] / d["cs_tcool0"]
        a.set_title(f"n$_{{bg}}$ = {n:g}: L/(c$_s$t$_{{cool}}$) = {Lx:.2g}")
        a.set_xlabel(r"t / t$_{cool}$")
        a.axvline(1, color=GRID, lw=1)
    ax[0].set_ylabel(r"rms $\delta\rho/\rho$ (cloud interior) / initial rms")
    ax[0].set_ylim(1e-2, 1e4)
    ax[0].legend(fontsize=8)
    ax[0].text(1.25, 2e3, "after ~t$_{cool}$ the collapsed\ncloud spans < 10 cells:\nsee resolution study", fontsize=7.5, color=INK2)
    save(fig, "fig8_amplitude_scan.png")


def fig_heating():
    R = runs("heating")
    if not R:
        return
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for a, n in zip(ax, [1e-3, 1e-2]):
        for c, mode in zip([S1, S2, S3], ["mask", "balance", "both"]):
            for d in R:
                if d["params"]["rho_bg"] == n and d["params"]["name"].endswith(mode):
                    tf = tc(d)
                    a.semilogy(d["t"] / tf, np.abs(d["amp_mode"]) / 0.02, color=c, label={
                        "mask": "cooling off outside cloud",
                        "balance": "heating = background cooling",
                        "both": "both"}[mode])
        a.set_xlabel(r"t / t$_{cool}$")
        a.set_title(f"n$_{{bg}}$ = {n:g}, L = 50 kpc, Z = 0.3, χ = 2")
    ax[0].set_ylabel(r"$|\delta\rho/\rho|/A_0$ (seeded mode)")
    ax[0].legend(fontsize=8)
    save(fig, "fig9_background_thermal_model.png")


def fig_resolution():
    R = runs("resolution")
    if not R:
        return
    A_runs = {d["params"]["amp"]: d for d in runs("amplitude") if d["params"]["rho_bg"] == 1e-4}
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9), sharey=True)
    for a, A in zip(ax, [1e-3, 1e-2, 1e-1]):
        sets = [("2048 cells / 150 kpc", A_runs.get(A), S1)]
        for d in R:
            if d["params"]["amp"] == A:
                sets.append((f"{d['params']['nx']} cells / 100 kpc", d,
                             S2 if d["params"]["nx"] == 4096 else S3))
        for lab, d, c in sets:
            if d is None:
                continue
            tf = tc(d)
            k = d["t"] <= 3 * tf
            a.semilogy(d["t"][k] / tf, np.maximum(d["amp_max"][k], 1e-4), color=c, lw=1.5, label=lab)
        a.set_title(f"Small cloud (L/c$_s$t$_{{cool}}$ ≈ 0.9), A$_0$ = {A:g}")
        a.set_xlabel(r"t / t$_{cool}$")
        a.set_ylim(1e-4, 1e2)
    ax[0].set_ylabel(r"interior max $\rho/\bar\rho - 1$")
    ax[0].legend(fontsize=8)
    save(fig, "fig10_resolution.png")


if __name__ == "__main__":
    fig_cooling()
    fig_validation()
    fig_linear_maps()
    fig_infinite_sims()
    fig_cloud()
    fig_profiles()
    fig_amplitude()
    # fig_heating()  # heated-background runs were numerically pathological (see REPORT)
    fig_resolution()
    fig_isochoric()
