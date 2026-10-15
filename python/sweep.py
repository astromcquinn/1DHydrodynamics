"""
Parameter sweeps of gPLUTO runs (parallel, resumable).

    python sweep.py infinite      # periodic "infinite cloud" runs
    python sweep.py cloud         # tophat cloud in a hot background
    python sweep.py amplitude     # initial-amplitude scan (nonlinearity)
    python sweep.py scaling       # density-scaling check
    python sweep.py heating       # cooling-masked vs heated background
    python sweep.py resolution    # resolution study of a small (isobaric) cloud
    python sweep.py infinite_uvb cloud_uvb amplitude_uvb
                                  # same with a photoionising UV background
    python sweep.py chi1          # chi = 1: no initial contrast, cooling switched on
    python sweep.py chi1_large chi1_amplitude chi1_resolution chi1_pert chi1_pert_big chi1_scaling
                                  # more chi = 1 sets
    python sweep.py longval       # infinite medium followed to 200 t_cool
    python sweep.py nonlinear     # nonlinear outcome vs hot-gas temperature
    python sweep.py pert_infinite pert_cloud
                                  # density, temperature, sound, velocity and
                                  # metallicity seeds

Each run is written to ../results/runs/<sweep>/<name>.npz; existing files are
skipped, so an interrupted sweep can simply be restarted.
Set PLUTO_EXE (compiled code) and PLUTO_WORK (scratch dir) in the environment.
"""
import os
import sys
import itertools
import time
from multiprocessing import Pool

import numpy as np

import cooling as C
import cooling_models as CM
import plutorun as P

NPROC = int(os.environ.get("NPROC", os.cpu_count() or 1))


def cs_tcool(n, T, Z):
    tab = C.CoolingTable.for_Z(Z)
    rho = n * C.mp
    tc = C.p_of(rho, T) / (C.GAMMA - 1) / (C.nenH_of_rho(rho) * tab(T))
    return C.cs_of(T) * tc / C.kpc


def _fmt(v):
    return ("%g" % v).replace("+", "").replace(".", "p").replace("-", "m")


def infinite_runs():
    runs = []
    n, T = 0.01, 1e6
    for Z, x, A in itertools.product([0.1, 0.3, 1.0],
                                     [0.03, 0.1, 0.3, 1, 3, 10, 30],
                                     [1e-3, 1e-2, 1e-1]):
        lam = x * cs_tcool(n, T, Z)
        name = f"inf_Z{_fmt(Z)}_x{_fmt(x)}_A{_fmt(A)}"
        runs.append(dict(name=name, rho_bg=n, T_bg=T, chi=1.0, L_cloud=-1,
                         lam=lam, Z=Z, amp=A, nx=512, tstop_factor=3.0,
                         nout=150))
    return runs


def cloud_runs():
    runs = []
    T_bg = 1e6
    for Z, nbg, L, nw, chi in itertools.product(
            [0.1, 0.3, 1.0], [1e-4, 1e-3, 1e-2], [30.0, 100.0],
            [2.5, 4.5, 8.5], [2.0]):
        name = f"cl_Z{_fmt(Z)}_n{_fmt(nbg)}_L{_fmt(L)}_nw{_fmt(nw)}_chi{_fmt(chi)}"
        runs.append(dict(name=name, rho_bg=nbg, T_bg=T_bg, chi=chi, L_cloud=L,
                         lam=L / nw, Z=Z, amp=0.02, nx=2048, box=3 * L,
                         tstop_factor=2.5, nout=120))
    # contrast scan at the poster's fiducial values
    for Z, chi, nbg in itertools.product([0.1, 0.3], [1.3, 5.0], [1e-3, 1e-2]):
        L, nw = 50.0, 4.5
        name = f"cl_Z{_fmt(Z)}_n{_fmt(nbg)}_L{_fmt(L)}_nw{_fmt(nw)}_chi{_fmt(chi)}"
        runs.append(dict(name=name, rho_bg=nbg, T_bg=T_bg, chi=chi, L_cloud=L,
                         lam=L / nw, Z=Z, amp=0.02, nx=2048, box=3 * L,
                         tstop_factor=2.5, nout=120))
    return runs


def amplitude_runs():
    runs = []
    T_bg, chi, Z, L, nw = 1e6, 2.0, 0.3, 50.0, 4.5
    for nbg, A in itertools.product([1e-4, 1e-3, 1e-2],
                                    [1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3]):
        name = f"amp_n{_fmt(nbg)}_A{_fmt(A)}"
        runs.append(dict(name=name, rho_bg=nbg, T_bg=T_bg, chi=chi, L_cloud=L,
                         lam=L / nw, Z=Z, amp=A, nx=2048, box=3 * L,
                         tstop_factor=2.5, nout=120))
    return runs


def scaling_runs():
    """Same L*n: the runs must be identical in scaled units."""
    runs = []
    for nbg, L in [(1e-2, 50.0), (1e-3, 500.0), (1e-4, 5000.0)]:
        name = f"scale_n{_fmt(nbg)}_L{_fmt(L)}"
        runs.append(dict(name=name, rho_bg=nbg, T_bg=1e6, chi=2.0, L_cloud=L,
                         lam=L / 4.5, Z=0.3, amp=0.02, nx=2048, box=3 * L,
                         tstop_factor=2.5, nout=120))
    return runs


def heating_runs():
    runs = []
    for nbg, mode in itertools.product([1e-3, 1e-2], ["mask", "balance", "both"]):
        kw = dict(cool_mask=1, heat_mode="none")
        if mode == "balance":
            kw = dict(cool_mask=0, heat_mode="balance")
        if mode == "both":
            kw = dict(cool_mask=1, heat_mode="balance")
        name = f"heat_n{_fmt(nbg)}_{mode}"
        runs.append(dict(name=name, rho_bg=nbg, T_bg=1e6, chi=2.0, L_cloud=50.0,
                         lam=50 / 4.5, Z=0.3, amp=0.02, nx=2048, box=150.0,
                         tstop_factor=2.5, nout=120, **kw))
    return runs


def resolution_runs():
    """Small (nearly isobaric) cloud: it is compressed ~T_c/T_floor = 50 times,
    so the post-collapse phase needs many more cells than the cooling phase."""
    runs = []
    for A, nx in itertools.product([1e-3, 1e-2, 1e-1], [4096, 16384]):
        name = f"res_A{_fmt(A)}_nx{nx}"
        runs.append(dict(name=name, rho_bg=1e-4, T_bg=1e6, chi=2.0, L_cloud=50.0,
                         lam=50 / 4.5, Z=0.3, amp=A, nx=nx, box=100.0,
                         tstop_factor=2.5, nout=120))
    return runs


def cs_tcool_model(n, T, Z, cooling):
    m = CM.make_model(cooling, Z)   # noqa: same for the "_z" variants
    rho = n * C.mp
    tc = C.p_of(rho, T) / (C.GAMMA - 1) / m.net(rho, T)
    return C.cs_of(T) * tc / C.kpc


XS_INF = [0.03, 0.1, 0.3, 1, 3, 10, 30]


def infinite_uvb_runs():
    """Infinite medium with the Cloudy tables: UV background ("pie") at three
    densities, and Cloudy CIE ("ciecl", density-independent) for reference."""
    runs = []
    T = 1e6
    combos = [("pie", n) for n in (1e-4, 1e-3, 1e-2)] + [("ciecl", 1e-2)]
    for (cool, n), Z, x in itertools.product(combos, [0.1, 0.3, 1.0], XS_INF):
        lam = x * cs_tcool_model(n, T, Z, cool)
        name = f"inf_{cool}_n{_fmt(n)}_Z{_fmt(Z)}_x{_fmt(x)}_A0p001"
        runs.append(dict(name=name, rho_bg=n, T_bg=T, chi=1.0, L_cloud=-1,
                         lam=lam, Z=Z, amp=1e-3, nx=512, tstop_factor=3.0,
                         nout=150, cooling=cool))
    for x in XS_INF:                       # nonlinear seeds
        lam = x * cs_tcool_model(1e-3, T, 0.3, "pie")
        name = f"inf_pie_n0p001_Z0p3_x{_fmt(x)}_A0p1"
        runs.append(dict(name=name, rho_bg=1e-3, T_bg=T, chi=1.0, L_cloud=-1,
                         lam=lam, Z=0.3, amp=0.1, nx=512, tstop_factor=3.0,
                         nout=150, cooling="pie"))
    return runs


def cloud_uvb_runs():
    runs = []
    for kw in cloud_runs():
        kw = dict(kw, cooling="pie", name=kw["name"].replace("cl_", "clpie_"))
        runs.append(kw)
    return runs


def amplitude_uvb_runs():
    return [dict(kw, cooling="pie", name=kw["name"].replace("amp_", "amppie_"))
            for kw in amplitude_runs()]


PERTS = ["density", "temperature", "sound", "velocity", "metal"]


def pert_infinite_runs():
    """Different kinds of seed perturbation, infinite medium, Cloudy tables
    with (pie_z) and without (ciecl_z) a UV background."""
    runs = []
    n, T, Z = 1e-3, 1e6, 0.3
    for cool, pt, x, A in itertools.product(["ciecl_z", "pie_z"], PERTS,
                                            [0.1, 1.0, 10.0], [1e-3, 0.1]):
        lam = x * cs_tcool_model(n, T, Z, cool)
        name = f"pinf_{cool}_{pt}_x{_fmt(x)}_A{_fmt(A)}"
        runs.append(dict(name=name, rho_bg=n, T_bg=T, chi=1.0, L_cloud=-1,
                         lam=lam, Z=Z, amp=A, nx=512, tstop_factor=3.0,
                         nout=150, cooling=cool, pert=pt))
    return runs


def pert_cloud_runs():
    runs = []
    for cool, pt, nbg in itertools.product(["ciecl_z", "pie_z"], PERTS, [1e-3, 1e-2]):
        L, nw = 50.0, 4.5
        name = f"pcl_{cool}_{pt}_n{_fmt(nbg)}"
        runs.append(dict(name=name, rho_bg=nbg, T_bg=1e6, chi=2.0, L_cloud=L,
                         lam=L / nw, Z=0.3, amp=0.02, nx=2048, box=3 * L,
                         tstop_factor=2.5, nout=120, cooling=cool, pert=pt))
    return runs


def chi1_runs():
    """chi = 1: the 'cloud' is ordinary background gas in which cooling is
    switched on (no initial density or temperature contrast)."""
    runs = []
    for cool, nbg, L, nw in itertools.product(["cie", "pie"], [1e-4, 1e-3, 1e-2],
                                              [30.0, 100.0], [2.5, 4.5, 8.5]):
        name = f"chi1_{cool}_n{_fmt(nbg)}_L{_fmt(L)}_nw{_fmt(nw)}"
        runs.append(dict(name=name, rho_bg=nbg, T_bg=1e6, chi=1.0, L_cloud=L,
                         lam=L / nw, Z=0.3, amp=0.02, nx=2048, box=3 * L,
                         tstop_factor=3.0, nout=150, cooling=cool))
    return runs


def _chi1(name, nbg, L, nw, cool="cie", A=0.02, **kw):
    d = dict(name=name, rho_bg=nbg, T_bg=1e6, chi=1.0, L_cloud=L, lam=L / nw,
             Z=0.3, amp=A, nx=2048, box=3 * L, tstop_factor=3.0, nout=150,
             cooling=cool)
    d.update(kw)
    return d


def chi1_large_runs():
    """chi = 1 clouds much larger than c_s t_cool (long-wavelength regime)."""
    return [_chi1(f"chi1L_{cool}_n0p01_L{_fmt(L)}_nw{_fmt(nw)}", 1e-2, L, nw, cool)
            for cool, L, nw in itertools.product(["cie", "pie"], [300.0, 1000.0, 3000.0],
                                                 [2.5, 4.5, 8.5])]


def chi1_amplitude_runs():
    """Seed-amplitude scan of a chi = 1 cloud (L = 50 kpc), run to 6 t_cool."""
    return [_chi1(f"chi1A_{cool}_n{_fmt(n)}_A{_fmt(A)}", n, 50.0, 4.5, cool, A=A,
                  tstop_factor=6.0, nout=300)
            for cool, n, A in itertools.product(["cie", "pie"], [1e-4, 1e-3, 1e-2],
                                                [1e-3, 1e-2, 3e-2, 0.1, 0.3])]


def chi1_resolution_runs():
    """The small chi = 1 cloud of the amplitude scan at 2x and 8x resolution."""
    return [_chi1(f"chi1R_A{_fmt(A)}_nx{nx}", 1e-4, 50.0, 4.5, "cie", A=A, nx=nx,
                  tstop_factor=3.0)
            for A, nx in itertools.product([1e-3, 1e-2, 0.1], [4096, 16384])]


def chi1_pert_runs():
    """The five seed types in a chi = 1 cloud (L = 50 kpc)."""
    return [_chi1(f"chi1P_{cool}_{pt}_n{_fmt(n)}", n, 50.0, 4.5, cool, pert=pt)
            for cool, pt, n in itertools.product(["ciecl_z", "pie_z"], PERTS, [1e-3, 1e-2])]


def chi1_pert_big_runs():
    """The five seeds at n = 1e-3 in a cloud large enough not to collapse
    (L = 200 kpc, same wavelength as chi1_pert)."""
    return [_chi1(f"chi1P_{cool}_{pt}_n0p001_L200", 1e-3, 200.0, 18.0, cool, pert=pt,
                  nx=6144, box=400.0)
            for cool, pt in itertools.product(["ciecl_z", "pie_z"], PERTS)]


def chi1_scaling_runs():
    return [_chi1(f"chi1S_n{_fmt(n)}_L{_fmt(L)}", n, L, 4.5)
            for n, L in [(1e-2, 50.0), (1e-3, 500.0), (1e-4, 5000.0)]]


def longval_runs():
    """Infinite medium followed for 200 t_cool (well past the floor), A = 1e-4."""
    runs = []
    for cool, n in [("cie", 1e-2), ("pie", 1e-3)]:
        for x in [0.1, 1.0, 10.0]:
            lam = x * cs_tcool_model(n, 1e6, 0.3, cool)
            runs.append(dict(name=f"lv_{cool}_x{_fmt(x)}", rho_bg=n, T_bg=1e6, chi=1.0,
                             L_cloud=-1, lam=lam, Z=0.3, amp=1e-4, nx=512,
                             tstop_factor=200.0, nout=2000, cooling=cool))
    return runs


def nonlinear_runs():
    """Nonlinear outcome (contrast with the hot gas, cf. Gronke & Oh 2020): chi = 1
    clouds with L/(c_s t_cool) = 0.3-10 and L/lambda = 4.5, seeds of 3-30%, for hot
    gas at 3e5, 1e6 and 3e6 K, with and without the UV background (n = 1e-3);
    plus infinite-medium runs with 30% seeds."""
    runs = []
    n = 1e-3
    for cool, T, Lx, A in itertools.product(["cie", "pie"], [3e5, 1e6, 3e6],
                                           [0.3, 1.0, 3.0, 10.0], [0.03, 0.1, 0.3]):
        L = Lx * cs_tcool_model(n, T, 0.3, cool)
        runs.append(dict(name=f"nl_{cool}_T{_fmt(T)}_Lx{_fmt(Lx)}_A{_fmt(A)}", rho_bg=n,
                         T_bg=T, chi=1.0, L_cloud=L, lam=L / 4.5, Z=0.3, amp=A, nx=4096,
                         box=3 * L, tstop_factor=4.0, nout=200, cooling=cool))
    for T, x in itertools.product([3e5, 1e6, 3e6], [0.1, 0.3, 1.0, 3.0]):
        lam = x * cs_tcool_model(n, T, 0.3, "cie")
        runs.append(dict(name=f"nlinf_cie_T{_fmt(T)}_x{_fmt(x)}_A0p3", rho_bg=n, T_bg=T,
                         chi=1.0, L_cloud=-1, lam=lam, Z=0.3, amp=0.3, nx=1024,
                         tstop_factor=4.0, nout=200, cooling="cie"))
    return runs


def nonlinear_res_runs():
    """Resolution check of the nonlinear sweep: CIE, A = 0.1, L/(c_s t_cool) = 0.3
    and 1, at 4x the cells (16384)."""
    return [dict(kw, nx=16384, name=kw["name"] + "_nx16384") for kw in nonlinear_runs()
            if kw["cooling"] == "cie" and kw["amp"] == 0.1 and kw["L_cloud"] > 0
            and any(f"_Lx{s}_" in kw["name"] for s in ("0p3", "1"))]


SWEEPS = dict(nonlinear=nonlinear_runs, nonlinear_res=nonlinear_res_runs, chi1=chi1_runs, chi1_large=chi1_large_runs, chi1_amplitude=chi1_amplitude_runs,
              chi1_resolution=chi1_resolution_runs, chi1_pert=chi1_pert_runs, chi1_pert_big=chi1_pert_big_runs,
              chi1_scaling=chi1_scaling_runs, longval=longval_runs, pert_infinite=pert_infinite_runs, pert_cloud=pert_cloud_runs,
              infinite_uvb=infinite_uvb_runs, cloud_uvb=cloud_uvb_runs,
              amplitude_uvb=amplitude_uvb_runs, resolution=resolution_runs, infinite=infinite_runs, cloud=cloud_runs, amplitude=amplitude_runs,
              scaling=scaling_runs, heating=heating_runs)


def _do(args):
    sweep, kw = args
    outdir = os.path.join(os.path.dirname(__file__), "..", "results", "runs", sweep)
    fn = os.path.join(outdir, kw["name"] + ".npz")
    if os.path.exists(fn):
        return kw["name"], "skip", 0.0
    t0 = time.time()
    try:
        r = P.Run(outdir=outdir, **kw)
        res = r.execute()
        status = "ok" if res["returncode"] == 0 else "rc=%d" % res["returncode"]
    except Exception as e:  # keep the sweep going
        status = "ERR " + repr(e)
    return kw["name"], status, time.time() - t0


if __name__ == "__main__":
    for sweep in sys.argv[1:]:
        jobs = [(sweep, kw) for kw in SWEEPS[sweep]()]
        # longest (lowest-density) runs first for better load balancing
        jobs.sort(key=lambda j: j[1]["rho_bg"] * max(j[1]["L_cloud"], 1.0))
        print(f"[{sweep}] {len(jobs)} runs on {NPROC} processes", flush=True)
        with Pool(NPROC) as pool:
            for name, status, dt in pool.imap_unordered(_do, jobs):
                print(f"  {name:45s} {status:8s} {dt:7.1f}s", flush=True)
