"""
Time normalisation for stored runs.

t_cool(d): time [Myr] for the unperturbed cloud gas of run d to cool
(isochorically) to 1.5 x its floor / equilibrium temperature -- the "end of
cooling" used on every time axis.  Runs made before this definition stored
t_floor_iso = time to reach 1.05 T_eq, which can be 2-3x later because the
last approach to the floor is slow; this recomputes it from the run's
parameters with the linear-theory background integrator (cached on disk).
"""
import json
import os

import cooling_models as CM
import linear_theory as LT

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
