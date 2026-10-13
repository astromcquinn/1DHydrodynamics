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
    ref = LT.integrate(N0, T0, 1.0, table=tab, background=bg, amp=0.0)
    tf = ref["t"][-1]                                   # time to reach the floor
    lam = x * ref["cs_tcool0_kpc"]
    t_out = np.linspace(0, 3 * tf, 1500)
    o = LT.integrate(N0, T0, lam, table=tab, background=bg, amp=1.0, Tstop=0,
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


if __name__ == "__main__":
    jobs = list(itertools.product(ZS, T0S, XS, BGS))
    nproc = int(os.environ.get("NPROC", 1))
    with Pool(nproc) as pool:
        rows = pool.map(one, jobs, chunksize=4)
    keys = rows[0].keys()
    out = {k: np.array([r[k] for r in rows]) for k in keys}
    path = os.path.join(os.path.dirname(__file__), "..", "results", "linear_maps.npz")
    np.savez(path, **out)
    print("wrote", path, len(rows), "rows")
