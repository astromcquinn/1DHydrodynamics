"""
Figures comparing different kinds of seed perturbation (density, temperature,
sound wave, velocity, metallicity).

    python make_figures_pert.py
"""
import glob
import os

import numpy as np

import cooling as C
import cooling_models as CM
import linear_theory as LT
import plutorun as P
from metrics import t_cool as tc
from make_figures import FIG, RES, S1, S2, S3, S4, S5, INK2, GRID, plt, save

PERTS = ["density", "temperature", "sound", "velocity", "metal"]
PCOL = dict(zip(PERTS, [S1, S2, S3, S4, S5]))
PLAB = {"density": "density (isobaric)", "temperature": "temperature (isochoric)",
        "sound": "sound wave", "velocity": "velocity", "metal": "metallicity"}
CLAB = {"ciecl_z": "CIE + floor", "pie_z": "UV background"}


def runs(sweep):
    return [P.load(f) for f in sorted(glob.glob(os.path.join(RES, "runs", sweep, "*.npz")))]


def fig_histories():
    R = {(d["params"]["cooling"], d["params"]["pert"], round(d["params"]["lam"] / d["cs_tcool0"], 2), d["params"]["amp"]): d
         for d in runs("pert_infinite")}
    fig, ax = plt.subplots(3, 2, figsize=(11, 9.5), sharex="row")
    for row, x in enumerate([0.1, 1.0, 10.0]):
        for pt in PERTS:
            d = R.get(("pie_z", pt, x, 1e-3))
            if d is None:
                continue
            m = CM.make_model("pie_z", 0.3)
            tf = tc(d)
            tt = np.linspace(0, d["t"][-1], 600)
            o = LT.integrate(1e-3, 1e6, d["params"]["lam"], model=m, amp=1.0, pert=pt,
                             Tstop=0, t_end=tt[-1] * C.Myr * 1.001, t_out=tt)
            c = PCOL[pt]
            ax[row, 0].plot(o["t"] / tf, o["s"], color=c, lw=1.5, label=PLAB[pt])
            ax[row, 0].plot(d["t"][::6] / tf, d["amp_mode"][::6] / 1e-3, "o", color=c, ms=3.5, mfc="none")
            ax[row, 1].plot(o["t"] / tf, o["dK"], color=c, lw=1.5)
            ax[row, 1].plot(d["t"][::6] / tf, d["amp_mode_K"][::6] / 1e-3, "o", color=c, ms=3.5, mfc="none")
        for a in ax[row]:
            a.set_yscale("symlog", linthresh=1)
            a.axvline(1, color=GRID, lw=1)
            a.axhline(0, color=INK2, lw=0.6)
        ax[row, 0].set_ylabel(r"$\delta\rho/\rho$ / A")
        ax[row, 1].set_ylabel(r"$\delta\ln K$ / A  (entropy)")
        ax[row, 0].set_title(fr"$\lambda = {x:g}\,c_s t_{{cool}}$: density response")
        ax[row, 1].set_title(fr"$\lambda = {x:g}\,c_s t_{{cool}}$: entropy perturbation")
    for a in ax[2]:
        a.set_xlabel(r"t / t$_{cool}$  (gas reaches T$_{eq}$ at 1)")
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=5, fontsize=8.5, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle("Infinite medium, UV background, n = 10$^{-3}$, Z = 0.3, T$_0$ = 10$^6$ K: "
                 "lines = linear theory, circles = gPLUTO (A = 10$^{-3}$)", fontsize=10, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(FIG, "fig18_pert_histories.png"), dpi=150)
    plt.close(fig)
    print("wrote fig18_pert_histories.png")


def fig_linear_survey():
    f = os.path.join(RES, "linear_maps_pert.npz")
    if not os.path.exists(f):
        return
    m = dict(np.load(f))
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.6), sharey=True, sharex=True)
    for row, cool in enumerate(["ciecl_z", "pie_z"]):
        for col, (key, ttl) in enumerate([("s_max_cool", "max while cooling"),
                                          ("s_max", "max by 3 t$_{cool}$")]):
            a = ax[row, col]
            for pt in PERTS:
                sel = (m["cooling"] == cool) & (m["pert"] == pt) & (m["bg"] == "isochoric")
                o = np.argsort(m["x"][sel])
                a.loglog(m["x"][sel][o], np.maximum(m[key][sel][o], 1e-3), color=PCOL[pt], lw=1.6,
                         label=PLAB[pt])
            a.axhline(1, color=INK2, lw=0.6, ls="--")
            a.set_title(f"{CLAB[cool]}: |δρ/ρ|/A, {ttl}")
            a.set_ylim(1e-2, 3e2)
        ax[row, 0].set_ylabel(r"$|\delta\rho/\rho|$ / A")
    for a in ax[1]:
        a.set_xlabel(r"$\lambda_0/(c_s t_{cool})_0$")
    ax[0, 0].legend(fontsize=8)
    fig.suptitle("Linear theory, infinite medium, n = 10$^{-3}$, Z = 0.3, T$_0$ = 10$^6$ K", fontsize=10)
    save(fig, "fig19_pert_linear_survey.png")


def fig_clouds(sweep="pert_cloud", out="fig20_pert_clouds.png", keep=None):
    R = [d for sw in ([sweep] if isinstance(sweep, str) else sweep) for d in runs(sw)]
    if keep is not None:
        R = [d for d in R if keep(d)]
    if not R:
        return
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), sharey=True)
    for a, nbg in zip(ax, [1e-3, 1e-2]):
        for i, pt in enumerate(PERTS):
            for j, cool in enumerate(["ciecl_z", "pie_z"]):
                for d in R:
                    p = d["params"]
                    if p["pert"] != pt or p["cooling"] != cool or p["rho_bg"] != nbg:
                        continue
                    tf = tc(d)
                    cool_ph = d["t"] <= tf
                    gc = np.abs(d["amp_mode"][cool_ph]).max() / p["amp"]
                    ge = np.abs(d["amp_mode"][d["t"] <= 2.5 * tf]).max() / p["amp"]
                    xpos = i + (j - 0.5) * 0.28
                    a.semilogy([xpos, xpos], [max(gc, 1e-2), max(ge, 1e-2)], color=PCOL[pt], lw=1.2)
                    a.semilogy(xpos, max(gc, 1e-2), "o", color=PCOL[pt], ms=6,
                               mfc="none" if cool == "ciecl_z" else PCOL[pt])
                    a.semilogy(xpos, max(ge, 1e-2), "^", color=PCOL[pt], ms=6,
                               mfc="none" if cool == "ciecl_z" else PCOL[pt])
                    Lx = p["L_cloud"] / d["cs_tcool0"]
                    x = p["lam"] / d["cs_tcool0"]
                    Lkpc, chi = p["L_cloud"], p["chi"]
        a.set_xticks(range(len(PERTS)))
        a.set_xticklabels([PLAB[p].split(" (")[0] for p in PERTS], fontsize=8.5)
        a.axhline(1, color=INK2, lw=0.6, ls="--")
        a.set_title(f"Cloud, n$_{{bg}}$ = {nbg:g}, L = {Lkpc:g} kpc, χ = {chi:g}: "
                    f"L ≈ {Lx:.2g}, λ ≈ {x:.2g} c$_s$t$_{{cool}}$")
    ax[0].set_ylabel(r"$|\delta\rho/\rho|$ / A of the seeded mode")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=INK2, marker="o", mfc="none", ls="none", label="max while cooling, CIE + floor"),
         Line2D([], [], color=INK2, marker="o", ls="none", label="max while cooling, UV background"),
         Line2D([], [], color=INK2, marker="^", ls="none", label="max by 2.5 t$_{cool}$")]
    ax[1].legend(handles=h, fontsize=7.5, loc="upper right")
    save(fig, out)


def fig_nonlinear():
    R = [d for d in runs("pert_infinite") if d["params"]["amp"] == 0.1]
    if not R:
        return
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), sharey=True)
    for a, cool in zip(ax, ["ciecl_z", "pie_z"]):
        for pt in PERTS:
            pts = sorted((d["params"]["lam"] / d["cs_tcool0"], d["amp_max"].max(), d["amp_max"][-1])
                         for d in R if d["params"]["pert"] == pt and d["params"]["cooling"] == cool)
            if not pts:
                continue
            pts = np.array(pts)
            a.loglog(pts[:, 0], np.maximum(pts[:, 1], 1e-3), "o-", color=PCOL[pt], ms=5, lw=1.2, label=PLAB[pt])
            a.loglog(pts[:, 0], np.maximum(pts[:, 2], 1e-3), "x:", color=PCOL[pt], ms=6, lw=0.9)
        a.set_xlabel(r"$\lambda_0/(c_s t_{cool})_0$")
        a.set_title(f"{CLAB[cool]}, A = 0.1: max ρ/ρ̄ − 1 (o) and at 3 t$_{{cool}}$ (x)")
    ax[0].set_ylabel(r"max $\rho/\bar\rho - 1$")
    ax[0].legend(fontsize=8)
    save(fig, "fig21_pert_nonlinear.png")


if __name__ == "__main__":
    fig_histories()
    fig_linear_survey()
    fig_clouds()
    fig_nonlinear()
