"""
Figures for the chi = 1 version of the write-up (cooling switched on in a
region of the 10^6 K gas, no initial contrast; the infinite medium is the
L -> infinity limit).

    python make_figures_paper.py
"""
import glob
import os

import numpy as np
from matplotlib.lines import Line2D

import cooling as C
import cooling_models as CM
import linear_theory as LT
import plutorun as P
from metrics import t_cool as tc
from sweep import _fmt
from make_figures import (FIG, RES, S1, S2, S3, S4, S5, CAT, INK2, GRID, plt, save,
                          cloud_metrics)

CHI1_SETS = ["chi1", "chi1_large"]


def runs(*sweeps):
    out = []
    for sw in sweeps:
        out += [P.load(f) for f in sorted(glob.glob(os.path.join(RES, "runs", sw, "*.npz")))]
    return out


def lt_long(cool, n, x, t_max):
    """Linear theory for an infinite medium out to t_max (in units of t_cool)."""
    m = CM.make_model(cool, 0.3)
    tf = LT.t_cool_end(n, 1e6, m)
    ref = LT.integrate(n, 1e6, 1.0, model=m, background="isochoric", amp=0.0,
                       stop_factor=LT.T_END_FACTOR)
    t_out = np.linspace(0, min(t_max, 3) * tf, 900)
    if t_max > 3:
        t_out = np.unique(np.concatenate([t_out, np.logspace(np.log10(3 * tf),
                                                              np.log10(t_max * tf), 20000)]))
    o = LT.integrate(n, 1e6, x * ref["cs_tcool0_kpc"], model=m, background="isochoric",
                     amp=1.0, Tstop=0, t_end=t_max * tf * C.Myr, t_out=t_out)
    o["tf"] = tf
    return o


def envelope(t, s, period, t0=1.0):
    """|s| up to t0 (end of cooling); afterwards the maximum of |s| over the
    preceding oscillation period (but not before t0), i.e. the amplitude of the
    sound wave the perturbation turns into."""
    s = np.abs(s)
    out = s.copy()
    lo = np.searchsorted(t, np.maximum(t0, t - period))
    for i in np.nonzero(t > t0)[0]:
        out[i] = s[lo[i]:i + 1].max()
    return out


# ---------------------------------------------------------------------------
def fig_validation_long():
    """Fig. 2: infinite medium followed to 200 t_cool, sims vs linear theory."""
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.0))
    for a, (cool, n, title) in zip(ax[:2], [("cie", 1e-2, "CIE + 10$^4$ K floor, n = 10$^{-2}$"),
                                             ("pie", 1e-3, "UV background, n = 10$^{-3}$")]):
        for c, x in zip([S1, S2, S3], [0.1, 1.0, 10.0]):
            f = os.path.join(RES, "runs", "longval", f"lv_{cool}_x{_fmt(x)}.npz")
            o = lt_long(cool, n, x, 1000.0)
            tf = o["tf"]
            t = o["t"] / tf
            # envelope over one cold-gas sound period
            cs_cold = C.cs_of(o["T"][-1]) / C.kpc * C.Myr        # kpc / Myr
            period = x * o["cs_tcool0_kpc"] / cs_cold / tf
            env = envelope(t, o["s"], period)
            a.loglog(t, env, color=c, lw=1.8,
                     label=fr"$\lambda/c_st_{{cool}}$ = {x:g}")
            if os.path.exists(f):
                d = P.load(f)
                td = d["t"] / tc(d)
                env_d = envelope(td, d["amp_mode"] / d["params"]["amp"], period)
                pick = np.unique(np.searchsorted(td, np.logspace(-1, np.log10(td[-1]), 40)))
                pick = pick[pick < len(td)]
                a.loglog(td[pick], env_d[pick], "o", ms=4, mfc="none", color=c)
        a.axvline(1, color=GRID, lw=1)
        t_H = 13.8e3 / LT.t_cool_end(n, 1e6, CM.make_model(cool, 0.3))
        a.axvline(t_H, color=INK2, lw=0.8, ls="--")
        a.text(t_H * 0.85, 0.015, "Hubble time", rotation=90, fontsize=7.5, color=INK2,
               ha="right", va="bottom")
        a.set_xlim(0.1, 1000)
        a.set_ylim(0.01, 300)
        a.set_xlabel(r"t / t$_{cool}$")
        a.set_title(title)
    ax[0].set_ylabel(r"$|\delta\rho/\rho| / A_0$")
    ax[0].legend(fontsize=8, loc="upper left")
    h = [Line2D([], [], color=INK2, lw=1.8, label="linear theory"),
         Line2D([], [], color=INK2, marker="o", mfc="none", ls="none", label="gPLUTO (to 200 t$_{cool}$)")]
    ax[1].legend(handles=h, fontsize=8, loc="upper left")

    # right: against background temperature while cooling, with the time-shift limit
    a = ax[2]
    tab = C.CoolingTable.for_Z(0.3)
    for c, x in zip([S1, S2, S3], [0.1, 1.0, 10.0]):
        o = lt_long("cie", 1e-2, x, 1.0)
        a.loglog(o["T"], np.abs(o["s"]), color=c, lw=1.8)
        f = os.path.join(RES, "runs", "longval", f"lv_cie_x{_fmt(x)}.npz")
        if os.path.exists(f):
            d = P.load(f)
            k = d["t"] <= tc(d)
            a.loglog(d["Tmean"][k], np.abs(d["amp_mode"][k]) / d["params"]["amp"], "o",
                     ms=4, mfc="none", color=c)
    TT = np.logspace(np.log10(1.5e4), 6, 200)
    a.loglog(TT, LT.timeshift_amplification(TT, 1e6, tab, "isochoric"), color=INK2, lw=1,
             ls="--", label="short-λ limit (time shift)")
    a.invert_xaxis()
    a.set_xlabel("background T [K]  (cooling →)")
    a.set_title("While cooling, against temperature (CIE)")
    a.legend(fontsize=8, loc="upper left")
    save(fig, "figP02_validation_long.png")


# ---------------------------------------------------------------------------
LBINS = [(0, 0.5, S1, "L/c$_s$t$_{cool}$ < 0.5"), (0.5, 1.5, S2, "0.5–1.5"),
         (1.5, 5, S3, "1.5–5"), (5, 50, S4, "5–50"), (50, 1e9, S5, "> 50")]


def lcol(Lx):
    for lo, hi, c, _ in LBINS:
        if lo <= Lx < hi:
            return c


def _lt_lines(a, cool, which):
    """Linear theory: infinite medium (dashed) and small isobaric cloud (dotted)."""
    if cool == "cie":
        m = dict(np.load(os.path.join(RES, "linear_maps.npz")))
        base = (m["Z"] == 0.3) & (m["T0"] == 1e6)
    else:
        m = dict(np.load(os.path.join(RES, "linear_maps_uvb.npz")))
        base = (m["cooling"] == "pie") & (m["n0"] == 1e-3) & (m["Z"] == 0.3) & (m["T0"] == 1e6)
    key = "G_max_cool" if which == "cool" else "G_2tf"
    for bg, ls, lab in [("isochoric", "--", "linear theory: infinite medium"),
                        ("isobaric", ":", "linear theory: small cloud (isobaric)")]:
        sel = base & (m["bg"] == bg)
        if bg == "isobaric":
            sel &= m["x"] <= 1.0
        o = np.argsort(m["x"][sel])
        a.loglog(m["x"][sel][o], m[key][sel][o], color=INK2, lw=1.1, ls=ls, label=lab)


def fig_chi1_clouds():
    R = runs(*CHI1_SETS)
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.0), sharex=True, sharey=True)
    for col, (cool, title) in enumerate([("cie", "CIE + 10$^4$ K floor"), ("pie", "UV background")]):
        for row, which in enumerate(["cool", "end"]):
            a = ax[row, col]
            for d in R:
                if d["params"]["cooling"] != cool:
                    continue
                mm = cloud_metrics(d)
                ok = mm["cells_per_lam_final"] >= 10
                c = lcol(mm["Lx"])
                a.loglog(mm["x"], mm["G_cool"] if which == "cool" else mm["G_end"], "o",
                         color=c, ms=5.5, mfc=c if ok else "none")
            _lt_lines(a, cool, which)
            if which == "cool":
                # small seeds: the same small clouds with A = 1e-3 (amplitude scan)
                for d in runs("chi1_amplitude"):
                    p = d["params"]
                    if p["cooling"] != cool or p["amp"] != 1e-3 or p["rho_bg"] > 1e-3:
                        continue
                    mm = cloud_metrics(d)
                    a.loglog(mm["x"], mm["G_cool"], "*", color=lcol(mm["Lx"]), ms=12,
                             mec="k", mew=0.6)
                a.axhline(1.3 / 0.02, color=INK2, lw=0.8, ls="-.")
                a.text(20, 1.3 / 0.02 * 1.15, r"$\delta\rho/\rho\approx1.3$ for A$_0$ = 0.02",
                       fontsize=7.5, color=INK2)
            a.set_title(f"{title}: " + ("max while cooling" if which == "cool"
                                         else "max by 2.5 t$_{cool}$ (2 t$_{cool}$ for lines)"))
            if row == 1:
                a.set_xlabel(r"$\lambda_0/(c_s t_{cool})_0$")
        ax[0, col].legend(fontsize=7.5, loc="lower left")
    for a in ax[:, 0]:
        a.set_ylabel(r"$|\delta\rho/\rho|/A_0$ of the seeded mode")
    h = [Line2D([], [], color=c, marker="o", ls="none", label=lab) for _, _, c, lab in LBINS]
    h += [Line2D([], [], color=INK2, marker="o", mfc="none", ls="none",
                 label="< 10 cells/λ after collapse"),
          Line2D([], [], color=INK2, marker="*", ms=11, mec="k", ls="none",
                 label="A$_0$ = 10$^{-3}$ (others: 0.02)")]
    ax[0, 1].legend(handles=h, fontsize=7.5, loc="upper right")
    ax[0, 0].set_ylim(0.5, 5e3)
    save(fig, "figP_chi1_clouds.png")


def fig_chi1_isochoric():
    """Long-wavelength (isochoric) regime: x >= 1.5 and L >= 1.5 c_s t_cool."""
    R = runs(*CHI1_SETS)
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0))
    pref = []
    for d in R:
        mm = cloud_metrics(d)
        if mm["Lx"] < 1.5 or mm["x"] < 1.5:
            continue
        tf = tc(d)
        k = d["t"] <= 0.95 * tf
        A = d["params"]["amp"]
        c = S1 if d["params"]["cooling"] == "cie" else S2
        ax[0].loglog(mm["x"], np.abs(d["amp_mode_T"][k]).max() / A, "o", color=c, ms=5)
        ax[0].loglog(mm["x"], mm["G_cool"], "s", color=c, ms=5, mfc="none")
        ax[1].loglog(mm["x"], mm["G_end"] - 1, "o", color=c, ms=5)
        pref.append((d["params"]["cooling"], mm["x"], (mm["G_end"] - 1) * mm["x"] ** 2))
    xx = np.logspace(np.log10(1.5), np.log10(600), 50)
    p = np.array([q[2] for q in pref])
    P0 = np.median(p)
    ax[1].loglog(xx, P0 / xx ** 2, color=INK2, lw=1, ls="--",
                 label=fr"${P0:.0f}\,(c_s t_{{cool}}/\lambda)^2$")
    ax[0].axhline(1, color=INK2, lw=0.8, ls="--")
    ax[0].set_ylabel(r"amplitude / A$_0$")
    ax[0].set_title("While cooling: δT/T grows (filled), δρ/ρ does not (open)")
    ax[1].set_ylabel(r"$G-1$ by 2.5 t$_{cool}$")
    ax[1].set_title("Afterwards: density grows from the pressure imprint")
    for a in ax:
        a.set_xlabel(r"$\lambda_0/(c_s t_{cool})_0$")
    h = [Line2D([], [], color=S1, marker="o", ls="none", label="CIE + floor"),
         Line2D([], [], color=S2, marker="o", ls="none", label="UV background")]
    ax[0].legend(handles=h, fontsize=8)
    ax[1].legend(fontsize=8)
    save(fig, "figP_chi1_isochoric.png")
    q = np.percentile(p, [10, 50, 90])
    print("isochoric prefactor 10/50/90%%: %.1f %.1f %.1f  (N=%d)" % (q[0], q[1], q[2], len(p)))
    for cool in ["cie", "pie"]:
        pp = [r[2] for r in pref if r[0] == cool]
        print("  ", cool, "median %.1f N=%d" % (np.median(pp), len(pp)))


def fig_chi1_profiles():
    picks = [("chi1_large", "chi1L_cie_n0p01_L1000_nw4p5", "CIE, n = 10$^{-2}$, L = 1000 kpc"),
             ("chi1", "chi1_cie_n0p001_L30_nw4p5", "CIE, n = 10$^{-3}$, L = 30 kpc"),
             ("chi1", "chi1_pie_n0p001_L100_nw4p5", "UV, n = 10$^{-3}$, L = 100 kpc")]
    fig, ax = plt.subplots(len(picks), 3, figsize=(13, 3.3 * len(picks)))
    for row, (sw, name, title) in enumerate(picks):
        f = os.path.join(RES, "runs", sw, name + ".npz")
        if not os.path.exists(f):
            continue
        d = P.load(f)
        prof, tp, x = d["profiles"], d["profile_t"], d["x"]
        L = d["params"]["L_cloud"]
        t_c = tc(d)
        cmap = plt.get_cmap("viridis")
        for i in range(len(tp)):
            c = cmap(i / max(len(tp) - 1, 1))
            rho, prs, T = prof[i][0], prof[i][1], prof[i][2]
            ax[row, 0].semilogy(x / L, rho, color=c, lw=1)
            ax[row, 1].semilogy(x / L, prs / prof[0][1].max(), color=c, lw=1)
            ax[row, 2].semilogy(x / L, T, color=c, lw=1,
                                label=f"{tp[i] / t_c:.2g} t$_{{cool}}$" if i % 3 == 0 or i == len(tp) - 1 else None)
        for a in ax[row]:
            a.set_xlim(-0.9, 0.9)
        ax[row, 0].set_ylabel(r"n [cm$^{-3}$]")
        ax[row, 1].set_ylabel(r"p / p$_{0}$")
        ax[row, 2].set_ylabel("T [K]")
        ax[row, 0].set_title(title + f"  (L/c$_s$t$_{{cool}}$ = {L / d['cs_tcool0']:.2g})")
        ax[row, 2].legend(fontsize=7, loc="lower right")
    for a in ax[-1]:
        a.set_xlabel("x / L")
    save(fig, "figP_chi1_profiles.png")


def fig_chi1_amplitude():
    R = runs("chi1_amplitude")
    fig, ax = plt.subplots(2, 3, figsize=(13, 6.8), sharey=True, sharex=True)
    amps = [1e-3, 1e-2, 3e-2, 0.1, 0.3]
    for row, (cool, lab) in enumerate([("cie", "CIE + floor"), ("pie", "UV")]):
        for a, n in zip(ax[row], [1e-4, 1e-3, 1e-2]):
            Lx = None
            for c, A in zip(CAT, amps):
                for d in R:
                    p = d["params"]
                    if p["rho_bg"] == n and p["amp"] == A and p["cooling"] == cool:
                        tf = tc(d)
                        a.semilogy(d["t"] / tf, np.maximum(d["amp_max"], 1e-4), color=c, lw=1.4,
                                   label=f"A$_0$ = {A:g}")
                        Lx = p["L_cloud"] / d["cs_tcool0"]
            if Lx is not None:
                a.set_title(f"{lab}: n = {n:g} (L/c$_s$t$_{{cool}}$ = {Lx:.2g})")
            a.set_ylim(1e-4, 5e2)
            a.axvline(1, color=GRID, lw=1)
            if row == 1:
                a.set_xlabel(r"t / t$_{cool}$")
        ax[row, 0].set_ylabel(r"interior max $\rho/\bar\rho - 1$")
    ax[0, 0].legend(fontsize=7.5)
    save(fig, "figP_chi1_amplitude.png")


def fig_chi1_resolution():
    R = runs("chi1_resolution")
    base = {d["params"]["amp"]: d for d in runs("chi1_amplitude")
            if d["params"]["rho_bg"] == 1e-4 and d["params"]["cooling"] == "cie"}
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9), sharey=True)
    for a, A in zip(ax, [1e-3, 1e-2, 1e-1]):
        sets = [(base.get(A), S1)] + [(d, S2 if d["params"]["nx"] == 4096 else S3)
                                      for d in R if d["params"]["amp"] == A]
        for d, c in sets:
            if d is None:
                continue
            tf = tc(d)
            k = d["t"] <= 3 * tf
            a.semilogy(d["t"][k] / tf, np.maximum(d["amp_max"][k], 1e-4), color=c, lw=1.5,
                       label=f"{d['params']['nx']} cells ({d['params']['nx'] // 3} across the cloud)")
            Lx = d["params"]["L_cloud"] / d["cs_tcool0"]
        a.set_title(f"Small cloud (L/c$_s$t$_{{cool}}$ ≈ {Lx:.2g}), A$_0$ = {A:g}")
        a.set_xlabel(r"t / t$_{cool}$")
        a.set_ylim(1e-4, 3e2)
    ax[0].set_ylabel(r"interior max $\rho/\bar\rho - 1$")
    ax[0].legend(fontsize=8)
    save(fig, "figP_chi1_resolution.png")


def fig_chi1_uv_vs_cie():
    R = runs(*CHI1_SETS)
    by = {d["params"]["name"]: d for d in R}
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    ratios = []
    for name, dp in by.items():
        if "_pie_" not in name:
            continue
        dc = by.get(name.replace("_pie_", "_cie_"))
        if dc is None:
            continue
        mp, mc = cloud_metrics(dp), cloud_metrics(dc)
        ok = min(mp["cells_per_lam_final"], mc["cells_per_lam_final"]) >= 10
        c = lcol(mc["Lx"])
        ax[0].loglog(mc["G_cool"], mp["G_cool"], "o", color=c, ms=5.5, mfc=c if ok else "none")
        ax[1].loglog(mc["G_end"], mp["G_end"], "o", color=c, ms=5.5, mfc=c if ok else "none")
        if ok:
            ratios.append((mp["G_cool"] / mc["G_cool"], mp["G_end"] / mc["G_end"]))
    for a, t in zip(ax, ["growth while cooling", "growth by 2.5 t$_{cool}$"]):
        lim = [0.8, 3e3]
        a.plot(lim, lim, color=INK2, lw=0.8, ls="--")
        a.set_xlabel("CIE + floor")
        a.set_ylabel("UV background")
        a.set_title("Same χ = 1 clouds: " + t)
    h = [Line2D([], [], color=c, marker="o", ls="none", label=lab) for _, _, c, lab in LBINS]
    h += [Line2D([], [], color=INK2, marker="o", mfc="none", ls="none", label="< 10 cells/λ after collapse")]
    ax[0].legend(handles=h, fontsize=7.5)
    save(fig, "figP_chi1_uv_vs_cie.png")
    r = np.array(ratios)
    for j, lab in enumerate(["cool", "2.5tc"]):
        print("UV/CIE", lab, "N=%d median %.2f 10-90%% %.2f-%.2f" % (
            len(r), np.median(r[:, j]), *np.percentile(r[:, j], [10, 90])))


def fig_chi1_seeds():
    import make_figures_pert as MP
    # n = 1e-3: the L = 200 kpc clouds (the L = 50 kpc ones collapse)
    MP.fig_clouds(["chi1_pert", "chi1_pert_big"], "figP_chi1_seeds.png",
                  keep=lambda d: d["params"]["rho_bg"] == 1e-2 or d["params"]["L_cloud"] == 200.0)


# ---------------------------------------------------------------------------
NL_TCOL = {3e5: S1, 1e6: S2, 3e6: S3}
NL_MARK = {0.3: "v", 1.0: "o", 3.0: "s", 10.0: "D"}


def nl_metrics(d):
    """Peak and late-time contrast between cold cloud gas and the hot ambient gas."""
    p = d["params"]
    t = d["t"] / tc(d)
    chi = np.asarray(d["nl_rho_cold"]) / np.asarray(d["nl_rho_hot"])
    chimax = np.asarray(d["nl_rho_max"]) / np.asarray(d["nl_rho_hot"])
    ok = np.isfinite(chi)
    late = ok & (t >= 0.8 * t[-1])
    i = np.nanargmax(np.where(ok, chi, np.nan)) if ok.any() else 0
    defic = 1 - np.asarray(d["nl_p_min"]) / np.asarray(d["nl_p_hot"])
    return dict(T=p["T_bg"], A=p["amp"], cool=p["cooling"],
                Lx=round(p["L_cloud"] / d["cs_tcool0"], 1) if p["L_cloud"] > 0 else np.inf,
                x=p["lam"] / d["cs_tcool0"],
                chi_peak=np.nanmax(chi) if ok.any() else np.nan, t_peak=t[i],
                chimax_peak=np.nanmax(chimax[ok]) if ok.any() else np.nan,
                chi_late=np.nanmedian(chi[late]) if late.any() else np.nan,
                cells_peak=int(np.asarray(d["nl_cells_densest"])[i]),
                ncold_late=int(np.median(np.asarray(d["nl_ncold"])[late])) if late.any() else 0,
                cold_frac=float(np.asarray(d["nl_cold_mass"])[-1] / np.asarray(d["nl_mass"])[0]),
                deficit=np.nanmax(defic[t <= 1.5]))


def fig_nonlinear():
    R = [d for d in runs("nonlinear") if d["params"]["L_cloud"] > 0]
    if not R:
        return
    M = [nl_metrics(d) for d in R]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for col, cool, title in [(0, "cie", "CIE + 10$^4$ K floor"), (1, "pie", "UV background")]:
        a = ax[col]
        for mm in M:
            if mm["cool"] != cool:
                continue
            c = NL_TCOL[mm["T"]]
            jit = {0.3: 0.85, 1.0: 0.95, 3.0: 1.05, 10.0: 1.15}[mm["Lx"]]
            a.loglog(mm["A"] * jit, mm["chi_peak"], NL_MARK[mm["Lx"]], color=c, ms=7)
            a.loglog(mm["A"] * jit, mm["chi_late"], NL_MARK[mm["Lx"]], color=c, ms=7, mfc="none")
        for T, c in NL_TCOL.items():
            a.axhline(T / 1e4, color=c, lw=0.8, ls=":")
        a.axhline(300, color=INK2, lw=1.2, ls="--")
        a.text(0.026, 330, "Gronke & Oh χ$_f$ ≈ 300", fontsize=8, color=INK2)
        a.set_xlabel("seed amplitude A")
        a.set_title(title)
        a.set_ylim(1, 3e4)
    ax[0].set_ylabel(r"cold / hot density, $\langle\rho_{\rm cold}\rangle/\rho_{\rm hot}$")
    h = [Line2D([], [], color=c, marker="o", ls="none", label=f"T$_{{hot}}$ = {T:.0e} K")
         for T, c in NL_TCOL.items()]
    h += [Line2D([], [], color=INK2, marker=m, ls="none", label=f"L/c$_s$t$_{{cool}}$ = {Lx:g}")
          for Lx, m in NL_MARK.items()]
    h += [Line2D([], [], color=INK2, marker="o", ls="none", label="peak"),
          Line2D([], [], color=INK2, marker="o", mfc="none", ls="none", label="late (last 20%)"),
          Line2D([], [], color=INK2, ls=":", label="T$_{hot}$/10$^4$ K")]
    ax[1].legend(handles=h, fontsize=7, loc="upper left", ncol=2)
    # pressure deficit vs late contrast
    a = ax[2]
    for mm in M:
        c = NL_TCOL[mm["T"]]
        a.semilogy(mm["deficit"], mm["chi_peak"], NL_MARK[mm["Lx"]], color=c, ms=7,
                   mfc=c if mm["cool"] == "cie" else "none")
    a.axhline(300, color=INK2, lw=1.2, ls="--")
    a.set_xlabel(r"max pressure deficit $1-p_{\min}/p_{\rm hot}$ (t ≤ 1.5 t$_{cool}$)")
    a.set_ylabel("peak cold / hot density")
    a.set_title("Contrast vs pressure deficit (filled: CIE, open: UV)")
    save(fig, "figP_nonlinear.png")
    print("%-5s %-8s %5s %5s %8s %8s %8s %6s %6s %6s %6s" % ("cool", "T", "Lx", "A", "chi_pk", "t_pk",
          "chi_lt", "cells", "ncold", "defic", "fcold"))
    for mm in sorted(M, key=lambda q: (q["cool"], q["T"], q["Lx"], q["A"])):
        print("%-5s %-8.0e %5g %5g %8.1f %8.2f %8.1f %6d %6d %6.2f %6.2f" % (
            mm["cool"], mm["T"], mm["Lx"], mm["A"], mm["chi_peak"], mm["t_peak"], mm["chi_late"],
            mm["cells_peak"], mm["ncold_late"], mm["deficit"], mm["cold_frac"]))


if __name__ == "__main__":
    import sys
    want = sys.argv[1:]
    for name, fn in list(globals().items()):
        if name.startswith("fig_") and callable(fn) and (not want or name in want):
            fn()
