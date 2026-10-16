"""
Figures for chi = 1 clouds: no initial density/temperature contrast, the
'cloud' is background gas in which cooling is switched on.

    python make_figures_chi1.py
"""
import glob
import os

import numpy as np

import plutorun as P
from metrics import t_end as tc  # legacy script: old time unit (time to cool)
from make_figures import FIG, RES, S1, S2, S3, S4, INK2, GRID, plt, save, cloud_metrics


def runs():
    return [P.load(f) for f in sorted(glob.glob(os.path.join(RES, "runs", "chi1", "*.npz")))]


def fig_chi1():
    R = runs()
    if not R:
        return
    mo = dict(np.load(os.path.join(RES, "linear_maps.npz")))
    mu = dict(np.load(os.path.join(RES, "linear_maps_uvb.npz")))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.1), sharey=True)
    bins = [(0, 0.5, S1, "L/c$_s$t$_{cool}$ < 0.5"), (0.5, 1.5, S2, "0.5–1.5"),
            (1.5, 5, S3, "1.5–5"), (5, 1e9, S4, "> 5")]

    def col(Lx):
        for lo, hi, c, _ in bins:
            if lo <= Lx < hi:
                return c

    for a, cool, title in [(ax[0], "cie", "Original CIE + floor"), (ax[1], "pie", "UV background")]:
        for d in R:
            p = d["params"]
            if p["cooling"] != cool:
                continue
            m = cloud_metrics(d)
            ok = m["cells_per_lam_final"] >= 10
            c = col(m["Lx"])
            a.loglog(m["x"], m["G_cool"], "o", color=c, ms=6, mfc=c if ok else "none")
            a.loglog(m["x"], m["G_end"], "^", color=c, ms=6, mfc=c if ok else "none")
        if cool == "cie":
            sel = (mo["bg"] == "isochoric") & (mo["Z"] == 0.3) & (mo["T0"] == 1e6)
            o = np.argsort(mo["x"][sel])
            a.loglog(mo["x"][sel][o], mo["G_max_cool"][sel][o], color=INK2, lw=1.2, ls="--",
                     label="infinite medium, while cooling")
            a.loglog(mo["x"][sel][o], mo["G_2tf"][sel][o], color=INK2, lw=0.9, ls=":",
                     label="infinite medium, by 2 t$_{cool}$")
        else:
            sel = ((mu["cooling"] == "pie") & (mu["n0"] == 1e-3) & (mu["bg"] == "isochoric")
                   & (mu["Z"] == 0.3) & (mu["T0"] == 1e6))
            o = np.argsort(mu["x"][sel])
            a.loglog(mu["x"][sel][o], mu["G_max_cool"][sel][o], color=INK2, lw=1.2, ls="--",
                     label="infinite medium (n = 10$^{-3}$), while cooling")
            a.loglog(mu["x"][sel][o], mu["G_2tf"][sel][o], color=INK2, lw=0.9, ls=":",
                     label="infinite medium, by 2 t$_{cool}$")
        a.set_xlabel(r"$\lambda_0/(c_s t_{cool})_0$")
        a.set_title(f"χ = 1 clouds, Z = 0.3, {title}")
        a.legend(fontsize=7.5, loc="lower left")
    ax[0].set_ylabel(r"$|\delta\rho/\rho|/A_0$ of seeded mode")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=c, marker="s", ls="none", label=lab) for _, _, c, lab in bins]
    h += [Line2D([], [], color=INK2, marker="o", ls="none", label="while cooling"),
          Line2D([], [], color=INK2, marker="^", ls="none", label="by 2.5 t$_{cool}$"),
          Line2D([], [], color=INK2, marker="o", mfc="none", ls="none", label="< 10 cells/λ after collapse")]
    ax[1].legend(handles=h, fontsize=7.5, loc="upper right")
    save(fig, "fig22_chi1_clouds.png")


def fig_chi1_profiles():
    picks = [("chi1_pie_n0p001_L100_nw4p5", "χ = 1, UV bkg, n$_{bg}$ = 10$^{-3}$, L = 100 kpc"),
             ("chi1_pie_n0p01_L30_nw4p5", "χ = 1, UV bkg, n$_{bg}$ = 10$^{-2}$, L = 30 kpc")]
    fig, ax = plt.subplots(2, 3, figsize=(13, 6.6))
    for row, (name, title) in enumerate(picks):
        f = os.path.join(RES, "runs", "chi1", name + ".npz")
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
            a.set_xlabel("x / L")
        ax[row, 0].set_ylabel(r"n [cm$^{-3}$]")
        ax[row, 1].set_ylabel(r"p / p$_{hot}$")
        ax[row, 2].set_ylabel("T [K]")
        ax[row, 0].set_title(title + f"  (L/c$_s$t$_{{cool}}$ = {L / d['cs_tcool0']:.2g})")
        ax[row, 2].legend(fontsize=7, loc="lower right")
    save(fig, "fig23_chi1_profiles.png")


if __name__ == "__main__":
    fig_chi1()
    fig_chi1_profiles()
