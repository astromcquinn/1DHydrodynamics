"""
Time normalisation for stored runs.

t_cool(d): the time unit [Myr], the initial isochoric e-folding cooling time
p/[(gamma-1) n_e n_H Lambda] of the cloud gas (stored in each run as tcool0).
It is also the time in the length unit c_s t_cool.

t_end(d): time [Myr] for the unperturbed cloud gas of run d to cool
(isochorically) to 1.5 x its floor / equilibrium temperature -- the "end of
cooling" (0.6 t_cool for CIE from 1e6 K).  Runs made before this definition stored
t_floor_iso = time to reach 1.05 T_eq, which can be 2-3x later because the
last approach to the floor is slow; this recomputes it from the run's
parameters with the linear-theory background integrator (cached on disk).
"""
import json
import os

import cooling_models as CM
import linear_theory as LT
from linear_theory import T_AFTER  # noqa: F401  (re-exported)

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "..", "results", "t_cool_cache.json")
_cache = None


def _load():
    global _cache
    if _cache is None:
        try:
            with open(CACHE) as f:
                _cache = json.load(f)
        except (OSError, ValueError):
            _cache = {}
    return _cache


def t_cool(d):
    """Time unit: initial e-folding cooling time [Myr]."""
    return float(d["tcool0"])


def t_end(d):
    """Time [Myr] for the unperturbed cloud gas to cool to 1.5 T_eq (end of cooling)."""
    p = d["params"]
    cooling = p.get("cooling", "cie")
    n_c, T_c = p["chi"] * p["rho_bg"], p["T_bg"] / p["chi"]
    heat = p.get("heat", 0.0)
    key = f"{cooling}|{p['Z']}|{n_c:.6g}|{T_c:.6g}|{heat:.6g}|{p.get('T_floor', 1e4)}|{LT.T_END_FACTOR}"
    c = _load()
    if key not in c:
        model = CM.make_model(cooling, p["Z"])
        if p.get("T_floor") is not None:
            model.Tfloor = p["T_floor"]
        c[key] = LT.t_cool_end(n_c, T_c, model, heat=heat)
        with open(CACHE, "w") as f:
            json.dump(c, f, indent=0)
    return c[key]


def t_cooled(d):
    """End of the cooling phase of run d [Myr]: t_end(d), or, if later, the first
    time the cloud's own mass-weighted mean temperature falls below 3 x the
    floor/equilibrium temperature.  A compressed cloud does not follow the
    isochoric reference history and can finish cooling later.  (The mean is
    used with 3 rather than 1.5 x T_eq because a nonlinear temperature
    perturbation delays the mean of an otherwise isochoric cloud by a few per
    cent; with 3 x T_eq uncompressed clouds end at t_end(d).)"""
    import numpy as np
    import cooling as C
    p = d["params"]
    model = CM.make_model(p.get("cooling", "cie"), p["Z"])
    if p.get("T_floor") is not None:
        model.Tfloor = p["T_floor"]
    T_end = 2 * LT.T_END_FACTOR * model.T_eq(p["chi"] * p["rho_bg"] * C.mp)
    tf = t_end(d)
    i = np.nonzero(d["Tmean"] <= T_end)[0]
    return max(tf, d["t"][i[0]]) if len(i) else tf
