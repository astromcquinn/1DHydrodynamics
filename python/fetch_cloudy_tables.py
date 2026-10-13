"""
Extract the z = 0 Cloudy cooling/heating tables used by Grackle into small
npz files in python/data/.

Source: Grackle data files (https://github.com/grackle-project/grackle_data_files)
  CloudyData_UVB=HM2012.h5 : photoionisation equilibrium in the Haardt & Madau
                             (2012) UV background, Smith et al. (2017)
  CloudyData_noUVB.h5      : collisional ionisation equilibrium
All rates are Lambda / n_H^2 [erg cm^3 s^-1]; the metal tables are for solar
metallicity and scale linearly with Z.

    python fetch_cloudy_tables.py [download_dir]
"""
import os
import sys
import urllib.request

import h5py
import numpy as np

BASE = "https://raw.githubusercontent.com/grackle-project/grackle_data_files/main/input/"
FILES = {"hm2012": "CloudyData_UVB%3DHM2012.h5", "noUVB": "CloudyData_noUVB.h5"}
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _clean(a):
    a = np.array(a, dtype=float)
    a[~np.isfinite(a)] = 0.0
    return np.maximum(a, 0.0)


def extract(path, kind, z=0.0):
    with h5py.File(path, "r") as f:
        attrs = f["CoolingRates/Primordial/Cooling"].attrs
        logn = np.array(attrs["Parameter1"], float)
        logT = np.log10(np.array(attrs["Temperature"], float))
        out = dict(logn=logn, logT=logT)
        names = [("prim_cool", "CoolingRates/Primordial/Cooling"),
                 ("metal_cool", "CoolingRates/Metals/Cooling")]
        if kind == "hm2012":
            names += [("prim_heat", "CoolingRates/Primordial/Heating"),
                      ("metal_heat", "CoolingRates/Metals/Heating")]
            zgrid = np.array(attrs["Parameter2"], float)
            iz = int(np.argmin(np.abs(zgrid - z)))
            out["z"] = zgrid[iz]
        for key, ds in names:
            a = f[ds][:]
            out[key] = _clean(a[:, iz, :] if a.ndim == 3 else a)
        if kind != "hm2012":
            out["prim_heat"] = np.zeros_like(out["prim_cool"])
            out["metal_heat"] = np.zeros_like(out["metal_cool"])
            out["z"] = np.nan
    return out


if __name__ == "__main__":
    dl = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(OUT, exist_ok=True)
    for kind, fn in FILES.items():
        local = os.path.join(dl, fn.replace("%3D", "="))
        if not os.path.exists(local):
            print("downloading", fn)
            urllib.request.urlretrieve(BASE + fn, local)
        d = extract(local, kind)
        np.savez_compressed(os.path.join(OUT, f"cloudy_{kind}_z0.npz"), **d)
        print("wrote", f"data/cloudy_{kind}_z0.npz")
