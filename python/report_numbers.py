"""Print the numbers quoted in results/REPORT.md, from the stored runs."""
import glob
import os

import numpy as np

import plutorun as P
from metrics import t_end as tc  # legacy script: old time unit (time to cool)
from make_figures import cloud_metrics

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")


def runs(sw):
    return [P.load(f) for f in sorted(glob.glob(os.path.join(RES, "runs", sw, "*.npz")))]


def G(d, upto=1.0):
    A = d["params"]["amp"]
    return np.abs(d["amp_mode"][d["t"] <= upto * tc(d)]).max() / A


def xval(d):
    return d["params"]["lam"] / d["cs_tcool0"]


def section(t):
    print("\n==", t)


section("infinite CIE (approx) sims, A=1e-3: G while cooling / by 2.5 t_cool")
for Z in (0.1, 0.3, 1.0):
    rr = sorted((xval(d), G(d), G(d, 2.5)) for d in runs("infinite")
                if d["params"]["Z"] == Z and d["params"]["amp"] == 1e-3)
    print(f"Z={Z}: " + "  ".join(f"x={x:.2g}:{g:.3g}/{g2:.3g}" for x, g, g2 in rr))
rr = sorted((xval(d), G(d), d["amp_max"].max()) for d in runs("infinite")
            if d["params"]["Z"] == 0.3 and d["params"]["amp"] == 0.1)
print("A=0.1 Z=0.3: " + "  ".join(f"x={x:.2g}:G={g:.3g},contrast={c:.3g}" for x, g, c in rr))

section("cloud CIE sims")
M = [cloud_metrics(d) for d in runs("cloud")]
big = [m for m in M if m["Lx"] >= 1.5]
print("N(L/cst>=1.5) =", len(big))
iso = [(m["G_end"] - 1) * m["x"] ** 2 for m in M if m["Lx"] >= 1.5 and m["x"] >= 1.5]
print("(G_2.5-1) x^2: median %.0f, 10-90%% %.0f-%.0f, N=%d" % (np.median(iso), *np.percentile(iso, [10, 90]), len(iso)))
for m in sorted(M, key=lambda m: m["x"]):
    if m["chi"] == 1.3:
        print("chi=1.3 Z=%.1f n=%g: x=%.2g Gcool=%.3g G2.5=%.3g" % (m["Z"], m["n"], m["x"], m["G_cool"], m["G_end"]))
ex = [(m["x"], m["G_cool"], m["G_end"]) for m in M if m["chi"] == 2.0 and m["Lx"] >= 1.5]
for lo, hi in [(0.15, 0.3), (0.8, 1.6), (1.6, 3), (4, 6), (10, 30)]:
    sel = [e for e in ex if lo <= e[0] < hi]
    if sel:
        print("x in [%g,%g): Gcool median %.3g  G2.5 median %.3g  (N=%d)" % (lo, hi, np.median([e[1] for e in sel]), np.median([e[2] for e in sel]), len(sel)))
# isochoric: dT/T growth before end of cooling
dT = []
for d in runs("cloud"):
    m = cloud_metrics(d)
    if m["Lx"] >= 1.5 and m["x"] >= 1.5:
        dT.append(np.abs(d["amp_mode_p"][d["t"] <= tc(d)]).max() / d["params"]["amp"])
print("isochoric dT/T growth before end of cooling: %.0f-%.0f (A0=0.02 -> %.2f-%.2f)" % (min(dT), max(dT), 0.02 * min(dT), 0.02 * max(dT)))

section("amplitude scan (L=50, chi=2, Z=0.3)")
for sw in ("amplitude", "amplitude_uvb"):
    for d in runs(sw):
        p = d["params"]
        t = d["t"] / tc(d)
        i15 = np.argmin(abs(t - 1.5)); i25 = np.argmin(abs(t - 2.5))
        print("%s n=%g A=%g: L/cst=%.2g G_cool=%.4g G2.5=%.4g peak=%.3g at1.5=%.3g at2.5=%.3g end=%.3g t_end/tc=%.1f" % (
            sw, p["rho_bg"], p["amp"], p["L_cloud"] / d["cs_tcool0"], G(d), G(d, 2.5), d["amp_max"].max(),
            d["amp_max"][i15], d["amp_max"][i25], d["amp_max"][-1], t[-1]))

section("resolution")
for d in runs("resolution"):
    p = d["params"]
    t = d["t"] / tc(d); i15 = np.argmin(abs(t - 1.5)); i25 = np.argmin(abs(t - 2.5))
    print("A=%g nx=%d G_cool=%.4g peak=%.3g at1.5=%.3g at2.5=%.3g" % (p["amp"], p["nx"], G(d), d["amp_max"].max(), d["amp_max"][i15], d["amp_max"][i25]))

section("infinite UV sims Z=0.3, A=1e-3: G while cooling / by 2.5 t_cool")
for cool, n in [("ciecl", 1e-2), ("pie", 1e-4), ("pie", 1e-3), ("pie", 1e-2)]:
    rr = sorted((xval(d), G(d), G(d, 2.5)) for d in runs("infinite_uvb")
                if d["params"]["cooling"] == cool and d["params"]["rho_bg"] == n and d["params"]["Z"] == 0.3 and d["params"]["amp"] == 1e-3)
    print(f"{cool} n={n}: " + "  ".join(f"x={x:.2g}:{g:.3g}/{g2:.3g}" for x, g, g2 in rr))
rr = sorted((xval(d), d["amp_max"].max()) for d in runs("infinite_uvb") if d["params"]["amp"] == 0.1)
print("pie A=0.1 contrast: " + "  ".join(f"x={x:.2g}:{c:.3g}" for x, c in rr))

section("cloud UV vs CIE pairs")
Rc = {d["params"]["name"].replace("cl_", "clpie_"): d for d in runs("cloud")}
rc, re, iso = [], [], []
for d in runs("cloud_uvb"):
    mp = cloud_metrics(d); dc = Rc[d["params"]["name"]]; mc = cloud_metrics(dc)
    if mp["cells_per_lam_final"] < 10 or mc["cells_per_lam_final"] < 10:
        continue
    rc.append(mp["G_cool"] / mc["G_cool"]); re.append(mp["G_end"] / mc["G_end"])
    if mp["Lx"] >= 1.5 and mp["x"] >= 1.5:
        iso.append((mp["G_end"] - 1) * mp["x"] ** 2)
print("pairs %d; Gcool ratio median %.2f (10-90%% %.2f-%.2f); G2.5 ratio median %.2f (%.2f-%.2f); pie iso coeff median %.0f (%.0f-%.0f) N=%d" % (
    len(rc), np.median(rc), *np.percentile(rc, [10, 90]), np.median(re), *np.percentile(re, [10, 90]), np.median(iso), *np.percentile(iso, [10, 90]), len(iso)))

section("seed types, infinite, A=1e-3: G while cooling / by 2.5 t_cool / value at 2.5 t_cool (signed)")
for cool in ("ciecl_z", "pie_z"):
    for d in sorted(runs("pert_infinite"), key=lambda d: (xval(d), d["params"]["pert"])):
        p = d["params"]
        if p["cooling"] != cool or p["amp"] != 1e-3:
            continue
        i = np.argmin(abs(d["t"] / tc(d) - 2.5))
        print("%s x=%.2g %-11s %.3g / %.3g / %+.3g" % (cool, xval(d), p["pert"], G(d), G(d, 2.5), d["amp_mode"][i] / p["amp"]))
for d in sorted(runs("pert_infinite"), key=lambda d: (d["params"]["cooling"], xval(d), d["params"]["pert"])):
    p = d["params"]
    if p["amp"] == 0.1:
        print("A=0.1 %s x=%.2g %-11s peak contrast %.3g end %.3g" % (p["cooling"], xval(d), p["pert"], d["amp_max"].max(), d["amp_max"][-1]))

section("seed types, clouds A=0.02")
for d in runs("pert_cloud"):
    p = d["params"]
    print("%s n=%g %-11s x=%.2g L/cst=%.2g G_cool=%.3g G2.5=%.3g peak interior contrast=%.3g size_end=%.1f" % (
        p["cooling"], p["rho_bg"], p["pert"], xval(d), p["L_cloud"] / d["cs_tcool0"], G(d), G(d, 2.5), d["amp_max"].max(), d["size"][-1]))
