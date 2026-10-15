"""
Driver for the gPLUTO cooling-cloud problem in ../pluto_setup.

    run = Run(name="test", rho_bg=0.01, T_bg=1e6, chi=1, L_cloud=-1, lam=10, Z=0.3)
    res = run.execute()          # runs PLUTO, analyses, deletes raw output

`execute` writes <outdir>/<name>.npz with the time series of diagnostics and a
few profiles, so a sweep of hundreds of runs takes little disk space.

Environment:  PLUTO_EXE  path to the compiled gPLUTO executable built from
              ../pluto_setup  (see README.md)
"""
import os
import shutil
import subprocess
import json
import numpy as np

import cooling as C
import cooling_models as CM
import linear_theory as LT

PLUTO_EXE = os.environ.get("PLUTO_EXE", "pluto")
WORK = os.environ.get("PLUTO_WORK", "/tmp/pluto_runs")

INI = """[Grid]

X1-grid   {x0:.8g}   {nx:d}   {x1:.8g}
X2-grid   0.0   1   1.0
X3-grid   0.0   1   1.0

[Time]

CFL            {cfl}
CFL_max_var    1.1
tstop          {tstop:.8g}
first_dt       1.e-8

[Solver]

Solver         hllc

[Boundary]

X1-beg    {bc}
X1-end    {bc}
X2-beg    periodic
X2-end    periodic
X3-beg    periodic
X3-end    periodic

[Static Grid Output]

uservar    0
dbl       {dt_out:.8g}  -1   single_file
flt       -1.0  -1   single_file
vtk       -1.0  -1   single_file
tab       -1.0  -1
ppm       -1.0  -1
png       -1.0  -1
log        1000
analysis  -1.0  -1

[Parameters]

RHO_BG          {rho_bg:.8g}
T_BG            {T_bg:.8g}
CHI             {chi:.8g}
L_CLOUD         {L_cloud:.8g}
EDGE_WIDTH      {edge:.8g}
AMP             {amp:.8g}
LAMBDA          {lam:.8g}
HEAT_CGS        {heat:.8g}
T_FLOOR         {T_floor:.8g}
MAX_COOL_RATE   {max_cool:.8g}
PHASE           {phase:.8g}
COOL_MASK       {cool_mask:d}
COOL_MODE       {cool_mode:d}
PERT_TYPE       {pert_type:d}
ZMET            {zmet:.8g}
"""


PERT_TYPES = {"density": 0, "temperature": 1, "sound": 2, "velocity": 3, "metal": 4}


def read_dbl(path, nx):
    a = np.fromfile(path, dtype="<f8")
    return a.reshape(a.size // nx, nx)  # rho vx1 vx2 vx3 prs tr1 [tr2 = Z]


class Run:
    def __init__(self, name, rho_bg=0.01, T_bg=1e6, chi=1.0, L_cloud=-1.0,
                 lam=10.0, Z=0.3, amp=1e-3, heat_mode="none", cool_mask=None,
                 T_floor=None, cooling="cie", pert="density", nx=None, box=None, edge=0.0, phase=0.0,
                 tstop=None, nout=100, tstop_factor=1.6, cfl=0.4,
                 max_cool=0.05, outdir="../results/runs", keep_raw=False):
        self.p = dict(name=name, rho_bg=rho_bg, T_bg=T_bg, chi=chi,
                      L_cloud=L_cloud, lam=lam, Z=Z, amp=amp,
                      heat_mode=heat_mode, cooling=cooling, pert=pert, edge=edge,
                      phase=phase, cfl=cfl, max_cool=max_cool)
        self.infinite = L_cloud <= 0
        if cool_mask is None:   # default: cooling off outside the cloud
            cool_mask = 0 if self.infinite else 1
        self.p["cool_mask"] = int(cool_mask)
        # cooling model: "cie" (original, approx. CIE + 1e4 K floor),
        # "ciecl" (Cloudy CIE + floor), "pie" (Cloudy, HM12 UV background)
        self.model = CM.make_model(cooling, Z)
        if T_floor is not None:
            self.model.Tfloor = T_floor
        T_floor = self.model.Tfloor
        self.p["T_floor"] = T_floor
        self.tab = C.CoolingTable.for_Z(Z)
        n_c, T_c = chi * rho_bg, T_bg / chi
        # heating = cooling of the hot background (per unit volume)
        if heat_mode == "balance":
            rho = rho_bg * C.mp
            self.p["heat"] = self.model.net(rho, T_bg)
        else:
            self.p["heat"] = 0.0
        # reference (unperturbed) isochoric history of the cloud gas
        if pert == "metal" and not hasattr(self.model, "dL_dlnZ"):
            raise ValueError('metal perturbations need a split table, e.g. cooling="pie_z"')
        ref = LT.integrate(n_c, T_c, lam, model=self.model, background="isochoric",
                           amp=0.0, heat=self.p["heat"], stop_factor=LT.T_END_FACTOR)
        self.t_floor_iso = ref["t"][-1]                       # Myr
        self.tcool0 = ref["tcool0_Myr"]
        self.cs_tcool0 = ref["cs_tcool0_kpc"]
        if tstop is None:
            tstop = tstop_factor * self.t_floor_iso
        self.p["tstop_Myr"] = tstop
        if self.infinite:
            box = box or lam
            nx = nx or 512
            x0, x1, bc = 0.0, box, "periodic"
        else:
            box = box or max(4.0 * L_cloud, 2.0 * L_cloud + 100.0)
            nx = nx or 4096
            x0, x1, bc = -0.5 * box, 0.5 * box, "outflow"
        self.p.update(nx=nx, box=box)
        self.x = x0 + (np.arange(nx) + 0.5) * (x1 - x0) / nx
        self.ini = INI.format(x0=x0, x1=x1, nx=nx, bc=bc, cfl=cfl,
                              tstop=tstop * C.Myr / C.UNIT_T,
                              dt_out=tstop * C.Myr / C.UNIT_T / nout,
                              rho_bg=rho_bg, T_bg=T_bg, chi=chi,
                              L_cloud=L_cloud, edge=edge, amp=amp, lam=lam,
                              heat=self.p["heat"], T_floor=T_floor,
                              max_cool=max_cool, phase=phase,
                              cool_mask=self.p["cool_mask"],
                              cool_mode=self.model.cool_mode,
                              pert_type=PERT_TYPES[pert], zmet=Z)
        self.outdir = outdir
        self.keep_raw = keep_raw

    # ------------------------------------------------------------------
    def execute(self, verbose=False):
        name = self.p["name"]
        d = os.path.join(WORK, name)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        with open(os.path.join(d, "pluto.ini"), "w") as f:
            f.write(self.ini)
        self.model.write_pluto(d)
        with open(os.path.join(d, "run.log"), "w") as log:
            ret = subprocess.run([os.path.abspath(PLUTO_EXE)], cwd=d, stdout=log,
                                 stderr=subprocess.STDOUT)
        res = self.analyse(d)
        res["returncode"] = ret.returncode
        os.makedirs(self.outdir, exist_ok=True)
        np.savez_compressed(os.path.join(self.outdir, name + ".npz"),
                            params=json.dumps(self.p), **res)
        if not self.keep_raw:
            shutil.rmtree(d, ignore_errors=True)
        return res

    # ------------------------------------------------------------------
    def _nonlinear_diag(self, rho, prs, T, trc, sel, dx, T_cold):
        """Diagnostics of the nonlinear outcome at one output.

        cold gas: cloud gas below T_cold = 2 T_eq (2 T_floor);
        hot gas:  for clouds, the non-cooling ambient medium (tracer < 0.5) above
                  0.5 T_bg; for the infinite medium, cells above half the current
                  maximum temperature.
        nl_rho_cold is the mass-weighted mean density of the cold gas, nl_ncold the
        number of separate cold sheets, nl_cells_densest the width (in cells) of
        the cold sheet that contains the densest cell."""
        nan = float("nan")
        if self.infinite:
            hot = T > 0.5 * T.max()
        else:
            hot = (trc < 0.5) & (T > 0.5 * self.p["T_bg"])
        cold = sel & (T < T_cold)
        out = dict(nl_rho_max=float(rho[sel].max()), nl_mass=float(rho[sel].sum() * dx),
                   nl_cold_mass=float(rho[cold].sum() * dx),
                   nl_rho_hot=float(np.median(rho[hot])) if hot.any() else nan,
                   nl_T_hot=float(np.median(T[hot])) if hot.any() else nan,
                   nl_p_hot=float(np.median(prs[hot])) if hot.any() else nan,
                   nl_p_min=float(prs[sel].min()))
        if cold.any():
            out["nl_rho_cold"] = float(np.sum(rho[cold] ** 2) / np.sum(rho[cold]))
            edges = np.diff(np.concatenate([[0], cold.astype(int), [0]]))
            starts, ends = np.nonzero(edges == 1)[0], np.nonzero(edges == -1)[0]
            out["nl_ncold"] = len(starts)
            imax = np.argmax(np.where(cold, rho, -1.0))
            j = np.nonzero((starts <= imax) & (ends > imax))[0]
            out["nl_cells_densest"] = int(ends[j[0]] - starts[j[0]]) if len(j) else 0
        else:
            out.update(nl_rho_cold=nan, nl_ncold=0, nl_cells_densest=0)
        return out

    def analyse(self, d):
        nx = self.p["nx"]
        times = []
        with open(os.path.join(d, "dbl.out")) as f:
            for line in f:
                times.append(float(line.split()[1]))
        times = np.array(times) * C.UNIT_T / C.Myr
        dx = self.p["box"] / nx
        keys = ["t", "Tmean", "pmean", "nmean", "amp_mode", "amp_mode_p",
                "amp_mode_T", "amp_mode_K", "amp_mode_v", "amp_mode_Z",
                "amp_h2", "amp_max", "amp_rms", "size", "nmax_over_nmean",
                # nonlinear outcome (see _nonlinear_diag)
                "nl_rho_max", "nl_rho_cold", "nl_rho_hot", "nl_T_hot", "nl_p_hot",
                "nl_p_min", "nl_mass", "nl_cold_mass", "nl_ncold", "nl_cells_densest"]
        g = C.GAMMA
        cs0 = np.sqrt(g * C.kB * self.p["T_bg"] / self.p["chi"] / (C.MU * C.amu)) / C.UNIT_V
        ts = {k: [] for k in keys}
        prof_idx = np.unique(np.linspace(0, len(times) - 1, 24).astype(int))
        T_cold = 2.0 * self.model.T_eq(self.p["chi"] * self.p["rho_bg"] * C.mp)
        profiles = []
        lam = self.p["lam"]
        rho_c0 = self.p["chi"] * self.p["rho_bg"]
        for i, t in enumerate(times):
            fn = os.path.join(d, "data.%04d.dbl" % i)
            if not os.path.exists(fn):
                break
            V = read_dbl(fn, nx)
            rho, vx, prs, trc = V[0], V[1], V[4], V[5]
            T = prs / rho * (C.UNIT_V**2 * C.amu / C.kB) * C.MU
            if self.infinite:
                sel = np.ones(nx, bool)
            else:
                sel = trc > 0.5
                if sel.sum() < 3:
                    # cloud collapsed to a few cells and diluted by numerical mixing:
                    # fall back to the cells with the most cloud material
                    sel = trc >= 0.5 * trc.max()
            m_c = rho[sel] * dx
            M = m_c.sum()
            # mass coordinate relative to the cloud's centre of mass
            mcoord = np.cumsum(m_c) - 0.5 * m_c - 0.5 * M
            km = 2 * np.pi / (lam * rho_c0)       # initial mass wavenumber
            nmean = M / (sel.sum() * dx)
            Tmean = np.sum(T[sel] * m_c) / M
            pmean = np.mean(prs[sel])
            if self.infinite:
                inner = np.ones(sel.sum(), bool)
            else:  # skip the outer 20% of the cloud mass on each side (edges)
                inner = np.abs(mcoord) < 0.3 * M
            dr = rho[sel][inner]
            pr = prs[sel][inner]
            mm = mcoord[inner]
            # least-squares projection on cos/sin(k_m m) and the 2nd harmonic
            A = np.vstack([np.ones_like(mm), np.cos(km * mm), np.sin(km * mm),
                           np.cos(2 * km * mm), np.sin(2 * km * mm)]).T
            cr, *_ = np.linalg.lstsq(A, dr, rcond=None)
            cp, *_ = np.linalg.lstsq(A, pr, rcond=None)
            # other fields: relative perturbations of T, entropy K = p/rho^g,
            # velocity (in units of the initial cloud sound speed) and Z
            Tin = T[sel][inner]
            cT, *_ = np.linalg.lstsq(A, Tin, rcond=None)
            lnK = np.log(pr) - g * np.log(dr)
            cK, *_ = np.linalg.lstsq(A, lnK - lnK.mean(), rcond=None)
            cv, *_ = np.linalg.lstsq(A, vx[sel][inner], rcond=None)
            if V.shape[0] > 6:
                cZ, *_ = np.linalg.lstsq(A, V[6][sel][inner], rcond=None)
                ampZ = np.hypot(cZ[1], cZ[2]) / cZ[0] * np.sign(cZ[1]) if cZ[0] > 0 else 0.0
            else:
                ampZ = 0.0
            ts["t"].append(t)
            ts["Tmean"].append(Tmean)
            ts["pmean"].append(pmean)
            ts["nmean"].append(nmean)
            ts["amp_mode"].append(np.hypot(cr[1], cr[2]) / cr[0] * np.sign(cr[1]))
            ts["amp_mode_p"].append(np.hypot(cp[1], cp[2]) / cp[0] * np.sign(cp[1]))
            ts["amp_mode_T"].append(np.hypot(cT[1], cT[2]) / cT[0] * np.sign(cT[1]))
            ts["amp_mode_K"].append(np.hypot(cK[1], cK[2]) * np.sign(cK[1]))
            ts["amp_mode_v"].append(np.hypot(cv[1], cv[2]) / cs0 * np.sign(cv[2]))
            ts["amp_mode_Z"].append(ampZ)
            ts["amp_h2"].append(np.hypot(cr[3], cr[4]) / cr[0])
            ts["amp_max"].append(dr.max() / dr.mean() - 1)
            ts["amp_rms"].append(np.std(dr) / dr.mean())
            ts["size"].append(sel.sum() * dx)
            # max over the whole cloud, edges included (for comparison)
            ts["nmax_over_nmean"].append(rho[sel].max() / nmean - 1)
            for k, v in self._nonlinear_diag(rho, prs, T, trc, sel, dx, T_cold).items():
                ts[k].append(v)
            if i in prof_idx:
                profiles.append(np.vstack([rho, prs, T, trc, vx]).astype(np.float32))
        out = {k: np.array(v) for k, v in ts.items()}
        out["profiles"] = np.array(profiles)
        out["profile_t"] = times[prof_idx[:len(profiles)]]
        out["x"] = self.x
        out["t_floor_iso"] = self.t_floor_iso
        out["tcool0"] = self.tcool0
        out["cs_tcool0"] = self.cs_tcool0
        return out


def load(path):
    d = dict(np.load(path, allow_pickle=False))
    d["params"] = json.loads(str(d["params"]))
    return d
