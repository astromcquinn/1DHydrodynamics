"""
Linear-theory survey: amplification of an isobaric seed perturbation in a
cooling background, as a function of x = lambda/(c_s t_cool,0), initial
temperature T0 and metallicity Z, for both backgrounds:

  isochoric : infinite uniform cooling medium ("infinite cloud")
  isobaric  : cloud much smaller than c_s t_cool held at fixed pressure

Because the equations are invariant under n -> a n, (x, t) -> (x, t)/a, the
density only enters through x; we fix n = 0.01 cm^-3.

Output: ../results/linear_maps.npz
"""
import itertools
import os
from multiprocessing import Pool

import numpy as np

import cooling as C
import linear_theory as LT

N0 = 0.01
ZS = [0.05, 0.1, 0.3, 1.0]
T0S = [1.0e6, 5.0e5, 2.5e5]
XS = np.logspace(-1.5, 2.5, 17)
BGS = ["isochoric", "isobaric"]


def one(args):
    Z, T0, x, bg = args
    tab = C.CoolingTable.for_Z(Z)
    ref = LT.integrate(N0, T0, 1.0, Z=Z, table=tab, background=bg, amp=0.0,
                       stop_factor=LT.T_END_FACTOR)
    tf = ref["t"][-1]                    # time to cool to 1.5 T_floor/T_eq
    lam = x * ref["cs_tcool0_kpc"]
    t_out = np.linspace(0, 3 * tf, 1500)
    o = LT.integrate(N0, T0, lam, Z=Z, table=tab, background=bg, amp=1.0, Tstop=0,
                     t_end=3 * tf * C.Myr, t_out=t_out)
    s = np.abs(o["s"])
    cool = o["t"] <= tf
    i_f = np.argmin(np.abs(o["t"] - tf))
    i_2 = np.argmin(np.abs(o["t"] - 2 * tf))
    # Das et al. (heating = cooling) comparison: growth over the same time tf
    # at the asymptotic rate appropriate for this x at the initial state
    LT0 = tab.dlnL_dlnT(T0)
    rate_isobaric = (2 - LT0) / C.GAMMA / ref["tcool0_Myr"]
    rate_isochoric = -LT0 / ref["tcool0_Myr"]
    return dict(Z=Z, T0=T0, x=x, bg=bg, tf=tf, tcool0=ref["tcool0_Myr"],
                G_floor=s[i_f], G_max_cool=s[cool].max(), G_2tf=s[: i_2 + 1].max(),
                G_3tf=s.max(), T_at_max=o["T"][cool][np.argmax(s[cool])],
                das_isobaric=np.exp(rate_isobaric * tf),
                das_isochoric=np.exp(rate_isochoric * tf),
                timeshift=LT.timeshift_amplification(o["T"][cool][np.argmax(s[cool])], T0, tab, bg))


def main_cie():
    jobs = list(itertools.product(ZS, T0S, XS, BGS))
    nproc = int(os.environ.get("NPROC", 1))
    with Pool(nproc) as pool:
        rows = pool.map(one, jobs, chunksize=4)
    keys = rows[0].keys()
    out = {k: np.array([r[k] for r in rows]) for k in keys}
    path = os.path.join(os.path.dirname(__file__), "..", "results", "linear_maps.npz")
    np.savez(path, **out)
    print("wrote", path, len(rows), "rows")


# ---------------------------------------------------------------------------
# With the Cloudy tables: UV background (density matters) and Cloudy CIE.
#     NPROC=2 python linear_maps.py uvb   ->  ../results/linear_maps_uvb.npz
# ---------------------------------------------------------------------------
XS_UVB = np.logspace(-1.5, 2.5, 13)
COMBOS_UVB = [("pie", 1e-4), ("pie", 1e-3), ("pie", 1e-2), ("ciecl", 1e-2)]


def one_uvb(args):
    import cooling_models as CM
    cool, n0, Z, T0, x, bg = args
    model = CM.make_model(cool, Z)
    ref = LT.integrate(n0, T0, 1.0, model=model, background=bg, amp=0.0,
                       stop_factor=LT.T_END_FACTOR)
    tf = ref["t"][-1]
    lam = x * ref["cs_tcool0_kpc"]
    t_out = np.linspace(0, 3 * tf, 1500)
    o = LT.integrate(n0, T0, lam, model=model, background=bg, amp=1.0, Tstop=0,
                     t_end=3 * tf * C.Myr, t_out=t_out)
    s = np.abs(o["s"])
    cool_ph = o["t"] <= tf
    i_2 = np.argmin(np.abs(o["t"] - 2 * tf))
    # Das-like rates at the initial state, generalised to L(rho, T):
    rho0 = n0 * C.mp
    L, Lr, LTd = model.derivs(rho0, T0)
    p0 = C.p_of(rho0, T0)
    rate_isobaric = (C.GAMMA - 1) * (Lr - LTd) / (C.GAMMA * p0) * C.Myr
    rate_isochoric = -(C.GAMMA - 1) * LTd / p0 * C.Myr
    return dict(cooling=cool, n0=n0, Z=Z, T0=T0, x=x, bg=bg, tf=tf,
                tcool0=ref["tcool0_Myr"], cs_tcool0=ref["cs_tcool0_kpc"],
                T_end=ref["T"][-1],
                G_max_cool=s[cool_ph].max(), G_2tf=s[: i_2 + 1].max(), G_3tf=s.max(),
                T_at_max=o["T"][cool_ph][np.argmax(s[cool_ph])],
                das_isobaric=np.exp(rate_isobaric * tf),
                das_isochoric=np.exp(rate_isochoric * tf))


def main_uvb():
    jobs = [(c, n) + rest for (c, n) in COMBOS_UVB
            for rest in itertools.product([0.1, 0.3, 1.0], [1.0e6, 5.0e5], XS_UVB, BGS)]
    nproc = int(os.environ.get("NPROC", 1))
    with Pool(nproc) as pool:
        rows = pool.map(one_uvb, jobs, chunksize=2)
    out = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    path = os.path.join(os.path.dirname(__file__), "..", "results", "linear_maps_uvb.npz")
    np.savez(path, **out)
    print("wrote", path, len(rows), "rows")


# ---------------------------------------------------------------------------
# Different kinds of seed perturbation (density, temperature, sound wave,
# velocity, metallicity) with the split Cloudy tables.
#     NPROC=4 python linear_maps.py pert  ->  ../results/linear_maps_pert.npz
# ---------------------------------------------------------------------------
PERTS = ["density", "temperature", "sound", "velocity", "metal"]


def one_pert(args):
    import cooling_models as CM
    cool, pert, n0, Z, T0, x, bg = args
    model = CM.make_model(cool, Z)
    ref = LT.integrate(n0, T0, 1.0, model=model, background=bg, amp=0.0,
                       stop_factor=LT.T_END_FACTOR)
    tf = ref["t"][-1]
    lam = x * ref["cs_tcool0_kpc"]
    t_out = np.linspace(0, 3 * tf, 1500)
    o = LT.integrate(n0, T0, lam, model=model, background=bg, amp=1.0, Tstop=0,
                     t_end=3 * tf * C.Myr, t_out=t_out, pert=pert)
    s, K = np.abs(o["s"]), np.abs(o["dK"])
    cool_ph = o["t"] <= tf
    late = o["t"] >= 2 * tf
    return dict(cooling=cool, pert=pert, n0=n0, Z=Z, T0=T0, x=x, bg=bg, tf=tf,
                s_max_cool=s[cool_ph].max(), s_max=s.max(),
                s_late_max=s[late].max(), s_late_min=s[late].min(),
                s_end=o["s"][-1], K_max=K.max(), K_end=o["dK"][-1],
                s_sign_cool=np.sign(o["s"][cool_ph][np.argmax(s[cool_ph])]))


def main_pert():
    jobs = list(itertools.product(["ciecl_z", "pie_z"], PERTS, [1e-3], [0.3], [1e6],
                                  XS_UVB, BGS))
    nproc = int(os.environ.get("NPROC", 1))
    with Pool(nproc) as pool:
        rows = pool.map(one_pert, jobs, chunksize=2)
    out = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    path = os.path.join(os.path.dirname(__file__), "..", "results", "linear_maps_pert.npz")
    np.savez(path, **out)
    print("wrote", path, len(rows), "rows")


if __name__ == "__main__":
    import sys
    if "uvb" in sys.argv[1:]:
        main_uvb()
    elif "pert" in sys.argv[1:]:
        main_pert()
    else:
        main_cie()
