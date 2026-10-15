"""
Regenerate the numbers quoted in the chi = 1 write-up (main.tex in
astromcquinn/cooling-perturbations-paper).

    python paper_numbers.py            # everything (linear-theory parts take ~1-2 min)
    python paper_numbers.py clouds     # one section; see SECTIONS at the bottom

Each section prints the numbers with the place in the draft they appear.
Older numbers for chi > 1 clouds are in report_numbers.py (superseded).
"""
import glob
import os
import sys

import numpy as np

import cooling as C
import cooling_models as CM
import linear_theory as LT
import plutorun as P
from metrics import t_cool as tc, t_cooled
from make_figures import cloud_metrics

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "..", "results", "runs")


def runs(*sweeps):
    return [P.load(f) for sw in sweeps for f in sorted(glob.glob(os.path.join(RUNS, sw, "*.npz")))]


def header(t):
    print("\n=== " + t + " ===")


# ---------------------------------------------------------------------------
def long_time():
    """Sec. 4.1 (Fig. 2): envelope of the linear-theory mode out to 1000 t_cool."""
    header("Long-time linear theory (Sec. 'Validation and the long-time evolution')")
    for cool, n in [("cie", 1e-2), ("pie", 1e-3)]:
        m = CM.make_model(cool, 0.3)
        tf = LT.t_cool_end(n, 1e6, m)
        ref = LT.integrate(n, 1e6, 1.0, model=m, background="isochoric", amp=0.0,
                           stop_factor=LT.T_END_FACTOR)
        print(f"{cool} n={n:g}: t_cool = {tf:.1f} Myr; Hubble time = {13.8e3 / tf:.0f} t_cool")
        for x in [0.1, 1, 10]:
            t_out = np.concatenate([np.linspace(0, 3 * tf, 600)[:-1],
                                    np.logspace(np.log10(3 * tf), np.log10(1000 * tf), 3000)])
            o = LT.integrate(n, 1e6, x * ref["cs_tcool0_kpc"], model=m, background="isochoric",
                             amp=1.0, Tstop=0, t_end=1000 * tf * C.Myr, t_out=t_out)
            s, t = np.abs(o["s"]), o["t"] / tf
            out = []
            for a, b in [(0, 1), (1, 3), (3, 10), (10, 30), (30, 100), (100, 300), (300, 1000)]:
                k = (t > a) & (t <= b)
                out.append(f"({a},{b}]:{s[k].max():.3g}")
            print(f"   x={x:<4g} max|s| in t/t_cool windows: " + " ".join(out))


def heated_table():
    """Sec. 2.4, Table 'sigma': t_cool and (2-Lambda_T)/gamma along the cooling."""
    header("Table of instantaneous rates (Sec. 'Comparison with a heated medium')")
    tab = C.CoolingTable.for_Z(0.3)
    rho = 1e-2 * C.mp
    for T in [1e6, 5e5, 3e5, 1e5, 3e4, 2e4]:
        tcool = C.p_of(rho, T) / (C.GAMMA - 1) / (C.nenH_of_rho(rho) * tab(T)) / C.Myr
        LTs = tab.dlnL_dlnT(T)
        print(f"T={T:.0e}: t_cool={tcool:.3g} Myr  Lambda_T={LTs:+.2f}  (2-LT)/g={(2 - LTs) / C.GAMMA:.2f}")


def isobaric_validity():
    """Sec. 2.2: growth of L/(c_s t_cool) along an isobaric track."""
    header("Isobaric-limit validity, L/(c_s t_cool) ~ Lambda/T^1.5 (Sec. 'The two backgrounds')")
    tab = C.CoolingTable.for_Z(0.3)
    for T in [5e5, 3e5, 1e5, 5e4, 2e4]:
        print(f"T={T:.0e}: ratio relative to 1e6 K = {tab(T) / T**1.5 / (tab(1e6) / 1e6**1.5):.3g}")


def thresholds():
    """Sec. 'Growth versus lambda' / Fig. 4: minimum seeds from linear_maps.npz."""
    header("Minimum seed amplitudes (Fig. min_amplitude)")
    m = dict(np.load(os.path.join(HERE, "..", "results", "linear_maps.npz")))
    for bg in ["isochoric", "isobaric"]:
        sel = (m["bg"] == bg) & (m["Z"] == 0.3) & (m["T0"] == 1e6)
        o = np.argsort(m["x"][sel])
        for x, g, g3, d3 in zip(m["x"][sel][o], m["G_max_cool"][sel][o], m["G_3tf"][sel][o],
                                m["dT_max_3tf"][sel][o]):
            print(f"{bg:9s} x={x:7.3g}  G_cool={g:7.3g}  G_3tc={g3:7.3g}  "
                  f"A_min(rho)={1 / g3:.2g}  A_min(T)={1 / d3:.2g}")


def clouds():
    """Sec. 'Clouds': Table 'clouds', small-cloud ranges, density scaling."""
    header("chi = 1 clouds (Table 'clouds' and text)")
    R = runs("chi1", "chi1_large")
    for cool in ["cie", "pie"]:
        M = [cloud_metrics(d) for d in R if d["params"]["cooling"] == cool]
        for lo, hi in [(0.3, 1), (1, 3), (3, 10), (10, 1e4)]:
            s = [mm for mm in M if lo <= mm["x"] < hi and mm["Lx"] >= 1.5 and mm["cells_per_lam_final"] >= 10]
            if s:
                print(f"{cool} x {lo}-{hi}: N={len(s)}  median G_cool={np.median([q['G_cool'] for q in s]):.3g}"
                      f"  median G_2.5={np.median([q['G_end'] for q in s]):.3g}")
        sm = [mm for mm in M if mm["Lx"] < 1.5]
        print(f"{cool} small clouds (L<1.5 c_s t_cool): G_cool {min(q['G_cool'] for q in sm):.0f}-"
              f"{max(q['G_cool'] for q in sm):.0f}, compression {min(q['compression'] for q in sm):.0f}-"
              f"{max(q['compression'] for q in sm):.0f}")
    print("density scaling (n, L) -> G_cool, G_2.5:")
    for d in runs("chi1_scaling"):
        mm = cloud_metrics(d)
        print(f"   {d['params']['name']}: {mm['G_cool']:.5g} {mm['G_end']:.5g}")


def isochoric():
    """Sec. 'Long-wavelength regime': prefactor of G-1 = P x^-2, and the Burkert & Lin check."""
    header("Long-wavelength regime")
    R = runs("chi1", "chi1_large")
    pref = {"cie": [], "pie": []}
    for d in R:
        mm = cloud_metrics(d)
        if mm["Lx"] >= 1.5 and mm["x"] >= 1.5:
            pref[d["params"]["cooling"]].append((mm["G_end"] - 1) * mm["x"] ** 2)
    allp = pref["cie"] + pref["pie"]
    print(f"prefactor: N={len(allp)} median {np.median(allp):.1f}, 10-90% "
          f"{np.percentile(allp, 10):.1f}-{np.percentile(allp, 90):.1f}; "
          f"CIE {np.median(pref['cie']):.1f}, UV {np.median(pref['pie']):.1f}")
    m = CM.make_model("cie", 0.3)
    n = 1e-2
    tf = LT.t_cool_end(n, 1e6, m)
    ref = LT.integrate(n, 1e6, 1.0, model=m, background="isochoric", amp=0.0,
                       stop_factor=LT.T_END_FACTOR)
    tc0 = ref["tcool0_Myr"]
    print(f"Burkert & Lin beta=0 coasting rate 4 pi^2/gamma = {4 * np.pi**2 / C.GAMMA:.1f} x^-2 per tau_c(0)")
    for x in [3, 10, 30, 100]:
        o = LT.integrate(n, 1e6, x * ref["cs_tcool0_kpc"], model=m, background="isochoric", amp=1.0,
                         Tstop=0, t_end=2.5 * tf * C.Myr, t_out=np.linspace(0, 2.5 * tf, 2001))
        i = np.argmin(abs(o["t"] - tf))
        rate = np.gradient(o["s"], o["t"])[i] * tc0 * x**2
        print(f"   x={x}: (G(tc)-1)x^2={(o['s'][i] - 1) * x**2:.1f}  ds/dt*tau_c0*x^2 at tc={rate:.1f}"
              f"  (G(2.5tc)-1)x^2={(o['s'][-1] - 1) * x**2:.1f}")


def uv_vs_cie():
    header("UV / CIE ratios for resolved cloud pairs")
    by = {d["params"]["name"]: d for d in runs("chi1", "chi1_large")}
    r = []
    for name, dp in by.items():
        if "_pie_" not in name or name.replace("_pie_", "_cie_") not in by:
            continue
        mp, mc = cloud_metrics(dp), cloud_metrics(by[name.replace("_pie_", "_cie_")])
        if min(mp["cells_per_lam_final"], mc["cells_per_lam_final"]) >= 10:
            r.append((mp["G_cool"] / mc["G_cool"], mp["G_end"] / mc["G_end"]))
    r = np.array(r)
    for j, lab in enumerate(["while cooling", "by 2.5 t_cool"]):
        print(f"{lab}: N={len(r)} median {np.median(r[:, j]):.2f} "
              f"(10-90% {np.percentile(r[:, j], 10):.2f}-{np.percentile(r[:, j], 90):.2f})")


def amplitude_resolution():
    header("Seed-amplitude scan and resolution study")
    for d in runs("chi1_amplitude", "chi1_resolution"):
        p = d["params"]
        mm = cloud_metrics(d)
        t = d["t"] / tc(d)
        k = t <= 3
        print(f"{p['name']:30s} nx={p['nx']:5d} G_cool={mm['G_cool']:8.4g} A*G={mm['G_cool'] * p['amp']:.2f}"
              f"  peak max={d['amp_max'][k].max():.3g} at {t[k][np.argmax(d['amp_max'][k])]:.2f} tc"
              f"  compression={mm['compression']:.0f}  end-of-cooling/t_cool={t_cooled(d) / tc(d):.2f}")


def small_cloud_history():
    header("Small-cloud compression history (amplitude scan, A = 1e-3)")
    for name in ["chi1A_cie_n0p0001_A0p001", "chi1A_cie_n0p001_A0p001", "chi1A_pie_n0p0001_A0p001"]:
        d = P.load(os.path.join(RUNS, "chi1_amplitude", name + ".npz"))
        t = d["t"] / tc(d)
        print(name)
        for T in [4e5, 2e5, 1e5, 5e4, 2e4]:
            i = np.argmax(d["Tmean"] < T)
            print(f"   T<{T:.0e} at t={t[i]:.2f} tc: n/n0={d['nmean'][i] / d['nmean'][0]:.3g} "
                  f"(isobaric {1e6 / d['Tmean'][i]:.3g}), G={abs(d['amp_mode'][i]) / d['params']['amp']:.3g}")


def seeds():
    header("Seed types in clouds")
    for d in runs("chi1_pert_big") + [d for d in runs("chi1_pert") if d["params"]["rho_bg"] == 1e-2]:
        p = d["params"]
        t = d["t"] / tc(d)
        A = p["amp"]
        print(f"{p['name']:40s} L/cstc={p['L_cloud'] / d['cs_tcool0']:.2g} x={p['lam'] / d['cs_tcool0']:.2g}"
              f"  cooling {np.abs(d['amp_mode'][t <= 1]).max() / A:.3g}  by 2.5tc "
              f"{np.abs(d['amp_mode'][t <= 2.5]).max() / A:.3g}")


def conduction():
    header("Field length vs c_s t_cool (Discussion)")
    tab = C.CoolingTable.for_Z(0.3)
    n = 1e-2
    rho = n * C.mp
    cst = C.cs_of(1e6) * C.p_of(rho, 1e6) / (C.GAMMA - 1) / (C.nenH_of_rho(rho) * tab(1e6))
    for T in [1e6, 5e5, 2e5, 1e5]:
        kappa = 6e-7 * T**2.5          # Spitzer, erg s^-1 cm^-1 K^-1
        lamF = np.sqrt(kappa * T / (C.nenH_of_rho(rho) * tab(T)))
        lamk = 2 * np.pi / np.sqrt(2 - tab.dlnL_dlnT(T)) * lamF
        print(f"T={T:.0e}: lambda_F={lamF / C.kpc:.3g} kpc, lambda_kappa={lamk / C.kpc:.3g} kpc, "
              f"x_kappa={lamk / cst:.3g}")


def runaway():
    """Discussion, Table 'runaway': where t_cool(T) drops below t_ff, and the
    short-wavelength amplification left from there to 1.5e4 K."""
    header("Runaway cooling (Table 'runaway')")
    tab = C.CoolingTable.for_Z(0.3)
    T = np.logspace(np.log10(1.5e4), 6, 4000)
    F = lambda Tx, Tf: (Tx / Tf) ** 2 * tab(Tf) / tab(Tx)
    for R in [3, 10, 30, 100]:
        out = []
        for tc_T, beta in [(T / tab(T), 1 / C.GAMMA), (T**2 / tab(T), 1.0)]:
            tc_T = tc_T / tc_T[-1]
            Tx = T[tc_T <= 1.0 / R].max()
            out.append(f"T_x={Tx:.2g} K amp={F(Tx, 1.5e4) ** beta:.0f}")
        print(f"t_cool0/t_ff={R}: isochoric {out[0]};  isobaric {out[1]}")
    for cool in ["cie", "pie"]:
        m = CM.make_model(cool, 0.3)
        print(cool, "time to cool [Myr]:", [round(LT.t_cool_end(n, 1e6, m)) for n in (1e-4, 1e-3, 1e-2)])
    for r in [50, 100, 200]:
        print(f"t_ff (v_c=200 km/s) at {r} kpc: {np.sqrt(2) * r * C.kpc / 200e5 / C.Myr:.0f} Myr")


SECTIONS = dict(runaway=runaway, long_time=long_time, heated_table=heated_table, isobaric_validity=isobaric_validity,
                thresholds=thresholds, clouds=clouds, isochoric=isochoric, uv_vs_cie=uv_vs_cie,
                amplitude_resolution=amplitude_resolution, small_cloud_history=small_cloud_history,
                seeds=seeds, conduction=conduction)

if __name__ == "__main__":
    for name in (sys.argv[1:] or SECTIONS):
        SECTIONS[name]()
