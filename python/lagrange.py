"""
A 1D Lagrangian (mass-coordinate) hydrodynamics code for the cooling-cloud
problem, as a companion to the gPLUTO setup in ../pluto_setup.

Why: on a fixed grid a cooling cloud that collapses by a factor 100-1000 in
density ends up a few cells wide, and hot gas is mixed into it numerically at
the edge.  In the mass coordinate every cell keeps its identity: the cloud edge
is a fixed mass coordinate (so the cooling switch is exact), a collapsing cloud
keeps all its cells, and a constant external pressure is a boundary condition,
which gives the fully isobaric background of the linear theory directly.

Scheme: staggered-grid von Neumann-Richtmyer.  Velocities live on the nodes,
density / internal energy / pressure in the cells.  Leapfrog in time, a
quadratic + linear artificial viscosity for shocks, and an implicit
(time-centred) p dV energy update that follows the adiabat exactly.  Cooling is
operator split and sub-cycled so that no sub-step changes a cell's internal
energy by more than `max_cool` (5%), as in the PLUTO setup.

Configurations (the same as plutorun.Run):
  * infinite medium  : L_cloud <= 0, periodic box of one wavelength, all gas cools
  * embedded cloud   : cloud of length L_cloud in a box of hot gas that does not
                       cool; constant pressure p0 at the outer boundaries
  * isobaric cloud   : box = L_cloud (cloud gas only) with constant pressure p0
                       at both ends: the cloud contracts at constant pressure,
                       the isobaric background of linear_theory.integrate

    import lagrange as LG
    r = LG.LRun("test", rho_bg=1e-2, T_bg=1e6, L_cloud=-1, lam=0.5, Z=0.3, amp=1e-4)
    res = r.execute()          # -> ../results/runs_lagrange/test.npz, same keys as plutorun
"""
import json
import os
import time

import numpy as np

import cooling as C
import cooling_models as CM
import linear_theory as LT

GAMMA = C.GAMMA
T_MIN = 10.0            # never let a cell cool below this


# ---------------------------------------------------------------------------
# vectorised net cooling L(rho, T) [erg cm^-3 s^-1] for the three models
# ---------------------------------------------------------------------------
def _cutoff_vec(T, Tfloor):
    y = np.minimum((T / Tfloor) ** 12, 700.0)
    return 1.0 - 1.0 / np.cosh(y)


class VecCooling:
    """Vectorised wrapper around a cooling model from cooling_models."""

    def __init__(self, model):
        self.m = model
        self.Tfloor = model.Tfloor
        self.kind = "cie" if isinstance(model, CM.CIEFloor) else "cloudy"
        if self.kind == "cloudy":
            self.logn, self.logT, self.lc, self.lh = model.logn, model.logT, model.lc, model.lh

    def _bilinear(self, tab, x, y):
        i = np.clip(np.searchsorted(self.logn, x, side="right") - 1, 0, len(self.logn) - 2)
        j = np.clip(np.searchsorted(self.logT, y, side="right") - 1, 0, len(self.logT) - 2)
        tx = np.clip((x - self.logn[i]) / (self.logn[i + 1] - self.logn[i]), 0.0, 1.0)
        ty = np.clip((y - self.logT[j]) / (self.logT[j + 1] - self.logT[j]), 0.0, 1.0)
        return ((1 - tx) * (1 - ty) * tab[i, j] + tx * (1 - ty) * tab[i + 1, j]
                + (1 - tx) * ty * tab[i, j + 1] + tx * ty * tab[i + 1, j + 1])

    def net(self, rho, T):
        if self.kind == "cie":
            return C.nenH_of_rho(rho) * self.m.tab(T) * _cutoff_vec(T, self.Tfloor)
        nH = C.nH_of_rho(rho)
        x, y = np.log10(nH), np.log10(T)
        c = 10.0 ** self._bilinear(self.lc, x, y)
        h = 10.0 ** self._bilinear(self.lh, x, y)
        return nH * nH * (c - h) * _cutoff_vec(T, self.Tfloor)


# ---------------------------------------------------------------------------
class LRun:
    def __init__(self, name, rho_bg=0.01, T_bg=1e6, L_cloud=-1.0, lam=10.0, Z=0.3,
                 amp=1e-3, cooling="cie", pert="density", T_floor=None,
                 nx=None, box=None, cloud_cells=None, phase=0.0,
                 tstop=None, tstop_factor=1.6, nout=100, cfl=0.4, max_cool=0.05,
                 cq=2.0, cl=0.5, outdir="../results/runs_lagrange", seed_in="cloud",
                 store_all=False):
        """
        nx          : total number of cells (default 512 infinite / 2048 cloud)
        cloud_cells : for an embedded cloud, how many of the nx cells are cloud
                      gas (default: the cloud's share of the box, as on a uniform
                      grid).  More cells in the cloud = better resolution there.
        box         : box length [kpc]; box == L_cloud gives the isobaric cloud
                      (cloud gas only, constant pressure at both ends).
        """
        self.p = dict(name=name, rho_bg=rho_bg, T_bg=T_bg, chi=1.0, L_cloud=L_cloud,
                      lam=lam, Z=Z, amp=amp, heat_mode="none", cooling=cooling,
                      pert=pert, edge=0.0, phase=phase, cfl=cfl, max_cool=max_cool,
                      cq=cq, cl=cl, heat=0.0, lagrangian=True)
        self.infinite = L_cloud <= 0
        self.model = CM.make_model(cooling, Z)
        if T_floor is not None:
            self.model.Tfloor = T_floor
        self.p["T_floor"] = self.model.Tfloor
        self.p["cool_mask"] = 0 if self.infinite else 1
        self.vc = VecCooling(self.model)
        ref = LT.integrate(rho_bg, T_bg, lam, model=self.model, background="isochoric",
                           amp=0.0, stop_factor=LT.T_END_FACTOR)
        self.t_floor_iso = ref["t"][-1]
        self.tcool0 = ref["tcool0_Myr"]
        self.cs_tcool0 = ref["cs_tcool0_kpc"]
        if tstop is None:
            tstop = tstop_factor * self.t_floor_iso
        self.p["tstop_Myr"] = tstop
        self.nout = nout
        # ---- grid (uniform in x initially; cell masses carry the seed) ----
        if self.infinite:
            box = box or lam
            nx = nx or 512
            x0, x1 = 0.0, box
            self.bc = "periodic"
        else:
            box = box or 3.0 * L_cloud
            nx = nx or 2048
            x0, x1 = -0.5 * box, 0.5 * box
            self.bc = "pressure"
        self.isobaric = (not self.infinite) and abs(box - L_cloud) < 1e-9 * L_cloud
        self.p.update(nx=nx, box=box)
        if self.infinite or self.isobaric:
            xn = np.linspace(x0, x1, nx + 1)
            mask = np.ones(nx, bool)
        else:
            # cloud_cells cells across the cloud, the rest split over the two wings
            if cloud_cells is None:
                cloud_cells = int(round(nx * L_cloud / box))
            nw = (nx - cloud_cells) // 2
            cloud_cells = nx - 2 * nw
            xn = np.concatenate([np.linspace(x0, -0.5 * L_cloud, nw + 1)[:-1],
                                 np.linspace(-0.5 * L_cloud, 0.5 * L_cloud, cloud_cells + 1)[:-1],
                                 np.linspace(0.5 * L_cloud, x1, nw + 1)])
            mask = np.zeros(nx, bool)
            mask[nw:nw + cloud_cells] = True
        self.p["cloud_cells"] = int(mask.sum())
        self.x_nodes0 = xn * C.kpc
        xc = 0.5 * (xn[1:] + xn[:-1])
        self.x = xc                                         # initial cell centres [kpc]
        self.mask = mask
        rho0 = rho_bg * C.mp
        p0 = C.p_of(rho0, T_bg)
        self.p0 = p0
        k = 2 * np.pi / (lam * C.kpc)
        xcen = 0.5 * (x0 + x1)                              # seed phase origin: box centre
        xs = (xc - xcen) * C.kpc
        wave = np.cos(k * xs + phase)
        seeded = mask if seed_in == "cloud" else np.ones(nx, bool)
        rho = np.full(nx, rho0)
        prs = np.full(nx, p0)
        u = np.zeros(nx + 1)
        if pert == "density":
            rho = np.where(seeded, rho0 * (1 + amp * wave), rho0)
        elif pert == "temperature":
            prs = np.where(seeded, p0 * (1 + amp * wave), p0)
        elif pert == "sound":
            rho = np.where(seeded, rho0 * (1 + amp * wave), rho0)
            prs = np.where(seeded, p0 * (1 + GAMMA * amp * wave), p0)
        elif pert == "velocity":
            xsn = self.x_nodes0 - xcen * C.kpc
            sn = (np.abs(xsn) <= 0.5 * L_cloud * C.kpc) if not self.infinite else np.ones(nx + 1, bool)
            u = np.where(sn, amp * C.cs_of(T_bg) * np.sin(k * xsn + phase), 0.0)
        else:
            raise ValueError(pert)
        dxc = np.diff(self.x_nodes0)
        self.dm = rho * dxc                                 # cell masses [g/cm^2]
        self.e = prs / ((GAMMA - 1) * rho)                  # specific internal energy
        self.u = u
        self.xn = self.x_nodes0.copy()
        self.outdir = outdir
        self.store_all = store_all                          # profiles at every output
        self.k_m = k / rho0                                 # mass wavenumber (initial)
        self.rho_c0 = rho0

    # ------------------------------------------------------------------
    def _state(self):
        V = np.diff(self.xn) / self.dm
        rho = 1.0 / V
        prs = (GAMMA - 1) * rho * self.e
        return V, rho, prs

    def _ghost_pressure(self, prs, q):
        """Pressure + viscosity in the two ghost cells."""
        if self.bc == "periodic":
            return prs[-1] + q[-1], prs[0] + q[0]
        return self.p0, self.p0

    def _step(self, dt_max):
        V, rho, prs = self._state()
        cs = np.sqrt(GAMMA * prs * V)
        du = np.diff(self.u)                                # node i+1 minus node i, per cell
        comp = du < 0
        q = np.where(comp, rho * (self.p["cq"] ** 2 * du * du + self.p["cl"] * cs * np.abs(du)), 0.0)
        # time step: sound crossing of each cell, and the viscous constraint
        dx = np.diff(self.xn)
        dt = self.p["cfl"] * np.min(dx / (cs + self.p["cq"] * np.abs(du) + 1e-30))
        dt = min(dt, dt_max)
        # momentum: nodes 1..N-1 interior; boundary nodes use ghost pressures
        P = prs + q
        gl, gr = self._ghost_pressure(prs, q)
        Pext = np.concatenate([[gl], P, [gr]])
        if self.bc == "periodic":
            mnode = 0.5 * (np.concatenate([[self.dm[-1]], self.dm]) + np.concatenate([self.dm, [self.dm[0]]]))
        else:
            mnode = 0.5 * (np.concatenate([[self.dm[0]], self.dm]) + np.concatenate([self.dm, [self.dm[-1]]]))
        acc = -(Pext[1:] - Pext[:-1]) / mnode
        self.u = self.u + dt * acc
        if self.bc == "periodic":
            self.u[-1] = self.u[0]
        self.xn = self.xn + dt * self.u
        # energy: implicit time-centred p dV (exact adiabat for an ideal gas)
        Vn = np.diff(self.xn) / self.dm
        dV = Vn - V
        a = (GAMMA - 1) * dV / (V + Vn)
        self.e = ((1 - a) * self.e - q * dV) / (1 + a)
        # cooling, sub-cycled
        self._cool(Vn, dt)
        return dt

    def _cool(self, V, dt):
        rho = 1.0 / V
        m = self.mask
        e = self.e[m]
        rho_m = rho[m]
        e_min = C.p_of(rho_m, T_MIN) / ((GAMMA - 1) * rho_m)

        def rate(e_):
            T = np.maximum(C.T_of(rho_m, (GAMMA - 1) * rho_m * e_), T_MIN)
            return self.vc.net(rho_m, T) / rho_m                        # erg g^-1 s^-1

        r0 = rate(e)
        n_sub = int(np.clip(np.ceil(dt * np.max(np.abs(r0) / e) / self.p["max_cool"]), 1, 20000))
        h = dt / n_sub
        for _ in range(n_sub):
            e1 = np.maximum(e - 0.5 * h * rate(e), e_min)
            e = np.maximum(e - h * rate(e1), e_min)
        self.e[m] = e

    # ------------------------------------------------------------------
    def execute(self, verbose=False):
        name = self.p["name"]
        tstop = self.p["tstop_Myr"] * C.Myr
        t_out = np.linspace(0.0, tstop, self.nout + 1)
        prof_idx = set(range(self.nout + 1)) if self.store_all else \
            set(np.unique(np.linspace(0, self.nout, 24).astype(int)))
        keys = ["t", "Tmean", "pmean", "nmean", "amp_mode", "amp_mode_p", "amp_mode_T",
                "amp_mode_K", "amp_mode_v", "amp_mode_Z", "amp_h2", "amp_max", "amp_rms",
                "size", "nmax_over_nmean", "nl_rho_max", "nl_rho_cold", "nl_rho_hot",
                "nl_T_hot", "nl_p_hot", "nl_p_min", "nl_mass", "nl_cold_mass", "nl_ncold",
                "nl_cells_densest"]
        ts = {k: [] for k in keys}
        profiles, profile_x, profile_t = [], [], []
        t, nstep, wall = 0.0, 0, time.time()
        for i_out, to in enumerate(t_out):
            while t < to - 1e-9 * tstop:
                dt = self._step(to - t)
                t += dt
                nstep += 1
            d = self._diag(t)
            for k in keys:
                ts[k].append(d[k])
            if i_out in prof_idx:
                V, rho, prs = self._state()
                xc = 0.5 * (self.xn[1:] + self.xn[:-1]) / C.kpc
                T = C.T_of(rho, prs)
                uc = 0.5 * (self.u[1:] + self.u[:-1]) / C.UNIT_V
                profiles.append(np.vstack([rho / C.mp, prs / (C.UNIT_RHO * C.UNIT_V ** 2), T,
                                           self.mask.astype(float), uc]).astype(
                                               np.float64 if self.store_all else np.float32))
                profile_x.append(xc.astype(np.float32))
                profile_t.append(t / C.Myr)
            if verbose:
                print(f"{name}: t = {t / C.Myr:8.3f} Myr ({t / C.Myr / self.tcool0:.3f} t_cool), "
                      f"{nstep} steps, Tmean = {d['Tmean']:.3g}, G = {d['amp_mode'] / self.p['amp']:.4g}",
                      flush=True)
        out = {k: np.array(v) for k, v in ts.items()}
        out["profiles"] = np.array(profiles)
        out["profile_x"] = np.array(profile_x)
        out["profile_t"] = np.array(profile_t)
        out["x"] = self.x
        out["t_floor_iso"] = self.t_floor_iso
        out["tcool0"] = self.tcool0
        out["cs_tcool0"] = self.cs_tcool0
        out["nsteps"] = nstep
        out["wall_s"] = time.time() - wall
        out["returncode"] = 0
        os.makedirs(self.outdir, exist_ok=True)
        np.savez_compressed(os.path.join(self.outdir, name + ".npz"), params=json.dumps(self.p), **out)
        return out

    # ------------------------------------------------------------------
    def _diag(self, t):
        """The diagnostics of plutorun.Run.analyse, in the mass coordinate."""
        V, rho, prs = self._state()
        T = C.T_of(rho, prs)
        uc = 0.5 * (self.u[1:] + self.u[:-1])
        sel = self.mask
        m_c = self.dm[sel]
        M = m_c.sum()
        mcoord = np.cumsum(m_c) - 0.5 * m_c - 0.5 * M
        km = self.k_m
        xw = np.diff(self.xn)
        nmean = M / xw[sel].sum()                           # mean density of the cloud gas
        Tmean = np.sum(T[sel] * m_c) / M
        pmean = np.mean(prs[sel])
        inner = np.ones(sel.sum(), bool) if self.infinite else (np.abs(mcoord) < 0.3 * M)
        dr, pr, Tin, mm = rho[sel][inner], prs[sel][inner], T[sel][inner], mcoord[inner]
        A = np.vstack([np.ones_like(mm), np.cos(km * mm), np.sin(km * mm),
                       np.cos(2 * km * mm), np.sin(2 * km * mm)]).T
        def proj(f):
            c, *_ = np.linalg.lstsq(A, f, rcond=None)
            return c
        cr, cp, cT = proj(dr), proj(pr), proj(Tin)
        lnK = np.log(pr) - GAMMA * np.log(dr)
        cK = proj(lnK - lnK.mean())
        cv = proj(uc[sel][inner])
        cs0 = C.cs_of(self.p["T_bg"])
        out = dict(t=t / C.Myr, Tmean=Tmean, pmean=pmean / (C.UNIT_RHO * C.UNIT_V ** 2), nmean=nmean / C.mp,
                   amp_mode=np.hypot(cr[1], cr[2]) / cr[0] * np.sign(cr[1]),
                   amp_mode_p=np.hypot(cp[1], cp[2]) / cp[0] * np.sign(cp[1]),
                   amp_mode_T=np.hypot(cT[1], cT[2]) / cT[0] * np.sign(cT[1]),
                   amp_mode_K=np.hypot(cK[1], cK[2]) * np.sign(cK[1]),
                   amp_mode_v=np.hypot(cv[1], cv[2]) / cs0 * np.sign(cv[2]),
                   amp_mode_Z=0.0, amp_h2=np.hypot(cr[3], cr[4]) / cr[0],
                   amp_max=dr.max() / dr.mean() - 1, amp_rms=np.std(dr) / dr.mean(),
                   size=xw[sel].sum() / C.kpc, nmax_over_nmean=rho[sel].max() / nmean - 1)
        # nonlinear-outcome diagnostics (as plutorun._nonlinear_diag)
        T_cold = 2.0 * self.model.T_eq(self.rho_c0)
        hot = (T > 0.5 * T.max()) if self.infinite else (~sel) & (T > 0.5 * self.p["T_bg"])
        cold = sel & (T < T_cold)
        nan = float("nan")
        out.update(nl_rho_max=rho[sel].max() / C.mp, nl_mass=M, nl_cold_mass=self.dm[cold].sum(),
                   nl_rho_hot=np.median(rho[hot]) / C.mp if hot.any() else nan,
                   nl_T_hot=np.median(T[hot]) if hot.any() else nan,
                   nl_p_hot=np.median(prs[hot]) / (C.UNIT_RHO * C.UNIT_V ** 2) if hot.any() else nan,
                   nl_p_min=prs[sel].min() / (C.UNIT_RHO * C.UNIT_V ** 2))
        if cold.any():
            out["nl_rho_cold"] = np.sum(self.dm[cold] * rho[cold]) / self.dm[cold].sum() / C.mp
            edges = np.diff(np.concatenate([[0], cold.astype(int), [0]]))
            starts, ends = np.nonzero(edges == 1)[0], np.nonzero(edges == -1)[0]
            out["nl_ncold"] = len(starts)
            imax = np.argmax(np.where(cold, rho, -1.0))
            j = np.nonzero((starts <= imax) & (ends > imax))[0]
            out["nl_cells_densest"] = int(ends[j[0]] - starts[j[0]]) if len(j) else 0
        else:
            out.update(nl_rho_cold=nan, nl_ncold=0, nl_cells_densest=0)
        return out


def load(path):
    d = dict(np.load(path, allow_pickle=False))
    d["params"] = json.loads(str(d["params"]))
    return d
