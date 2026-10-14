# Growth of perturbations in cooling gas (1D, PLUTO)

Tools for studying how a small density perturbation evolves inside gas that is
**cooling** (no heating inside the cloud), either in a finite cloud embedded in a hot
background or in an infinite uniform medium. Compare Das, Choudhury & Sharma (2021,
MNRAS), where heating balances cooling on average.

Results and figures are in [`results/REPORT.md`](results/REPORT.md).

## Layout

| path | what |
|---|---|
| `pluto_setup/` | gPLUTO problem (`definitions.hpp`, `init.cpp`, `radiat.cpp`, example `pluto.ini`, `build_gpluto.sh`) |
| `pluto_setup/pluto4.4_untested/` | the same problem ported to legacy PLUTO 4.4 (C). **Not compiled or tested** |
| `python/cooling.py` | cooling function Λ(T, Z), composition constants (identical to `radiat.cpp`) |
| `python/cooling_models.py` | the three cooling models (approx. CIE + floor; Cloudy CIE + floor; Cloudy with the HM12 UV background), each giving the net loss L(ρ, T), its derivatives, T_eq and the PLUTO tables |
| `python/fetch_cloudy_tables.py`, `python/data/` | z = 0 slices of Grackle's Cloudy tables (Smith et al. 2017; HM2012 UV background and no-UVB) |
| `python/linear_theory.py` | linearised equations integrated along a *cooling* background (infinite medium or small isobaric cloud) |
| `python/plutorun.py` | writes `pluto.ini` and `cooltable.dat`, runs gPLUTO, measures the mode amplitude, saves a compact `.npz` and deletes the raw output |
| `python/sweep.py` | parallel, resumable parameter sweeps (`infinite`, `cloud`, `amplitude`, `scaling`, `heating`) |
| `python/linear_maps.py` | linear-theory survey over λ/(c_s t_cool), T₀, Z |
| `python/make_figures.py` | all figures in `results/figures/` |
| `python/make_cooltable.py` | write a `cooltable.dat` for a given Z |
| `results/` | sweep outputs (`runs/*/*.npz`), `linear_maps.npz`, figures, report |

## The PLUTO problem

Code units: length 1 kpc, velocity 100 km/s (time unit 9.78 Myr), density m_p
(code density = ρ/m_p in cm⁻³). The parameters in `pluto.ini` are read **by position**
(gPLUTO skips the line after `[Parameters]`, so keep that blank line):

| parameter | meaning |
|---|---|
| `RHO_BG`, `T_BG` | density (ρ/m_p) and temperature of the hot background |
| `CHI` | cloud density contrast; the cloud has T = T_BG/CHI (pressure equilibrium) |
| `L_CLOUD` | cloud length [kpc]. **≤ 0 means an infinite cloud**: the whole periodic box is cloud gas |
| `EDGE_WIDTH` | tanh smoothing of the tophat edges [kpc] (0 = sharp) |
| `AMP`, `LAMBDA`, `PHASE` | δρ/ρ = AMP cos(2π(x−x_c)/LAMBDA + PHASE) inside the cloud, at uniform pressure (isobaric seed) |
| `COOL_MASK` | 1: cooling is multiplied by the passive tracer, so only cloud gas cools ("cooling off outside the cloud"); 2: sharp version, cool only where the tracer exceeds 0.5 |
| `HEAT_CGS` | uniform constant volumetric heating [erg cm⁻³ s⁻¹] (e.g. the background's cooling rate) |
| `T_FLOOR` | cooling (and heating) switched off smoothly below this temperature (use ~50 K with a UV background) |
| `MAX_COOL_RATE` | max fractional pressure change per step from cooling (time-step control) |
| `COOL_MODE` | 0: 1D table `cooltable.dat`, Λ(T) n_e n_H. 1: 2D table `cooltable2d.dat` of cooling **and** heating rates vs (n_H, T), Λ n_H², e.g. photoionisation equilibrium in a UV background |

`radiat.cpp` replaces gPLUTO's tabulated cooling. It adds the heating, the tracer
mask and the 2D (density, temperature) table, and uses n_e for a fully ionised
H+He+metal plasma; the distributed file leaves out the helium electrons. With
`COOL_MODE = 1` the mean molecular weight is still that of fully ionised gas,
which slightly misplaces T near 10⁴ K where hydrogen starts to recombine. Unlike PLUTO 4.4, gPLUTO never initialises
`minCoolingTemp`, so `init.cpp` sets it explicitly.

### Build and run (gPLUTO, CPU)

```bash
git clone https://gitlab.com/PLUTO-code/gPLUTO.git      # tested with commit 870aa80
export PLUTO_DIR=$PWD/gPLUTO
pluto_setup/build_gpluto.sh /path/to/build               # g++, no MPI, ~20 s
cd /path/to/build && python ../1DHydrodynamics/python/make_cooltable.py 0.3 && ./pluto
```

### Sweeps

```bash
pip install numpy scipy matplotlib
cd python
export PLUTO_EXE=/path/to/build/pluto PLUTO_WORK=/tmp/pluto_runs NPROC=4
python sweep.py infinite cloud amplitude scaling    # resumable
python sweep.py infinite_uvb cloud_uvb amplitude_uvb # with a UV background
python sweep.py chi1                                 # chi = 1: no initial contrast
python sweep.py pert_infinite pert_cloud             # other kinds of seed
NPROC=4 python linear_maps.py          # original cooling
NPROC=4 python linear_maps.py uvb      # Cloudy tables, with and without UV background
python make_figures.py && python make_figures_uvb.py
```

A single cloud run (2048 cells) takes 1 s to a few minutes on one core. The cost
grows as 1/(n L) because small or low-density clouds need many sound-crossing
times to cool.

### Using the script directly

```python
import plutorun as P, linear_theory as LT
r = P.Run("mytest", rho_bg=1e-3, T_bg=1e6, chi=2, L_cloud=50, lam=50/4.5,
          Z=0.3, amp=0.02, cooling="pie",        # "cie" (default), "ciecl", "pie"
          outdir="../results/runs/mine")
res = r.execute()       # res["t"], res["Tmean"], res["amp_mode"], ...
lt = LT.integrate(2e-3, 5e5, 50/4.5, Z=0.3, background="isobaric", amp=0.02)
```

## Cooling models

* `cie` (default, used for the original survey): Λ(T, Z) from `cooling.py`, with a
  hard floor at 10⁴ K.
* `ciecl`: Cloudy collisional-equilibrium tables from Grackle, same floor.
* `pie`: Cloudy photoionisation-equilibrium tables in the Haardt & Madau (2012)
  UV background at z = 0, including photoheating, with no floor. Gas settles at
  T_eq(n) ≈ 1.9×10⁴ / 1.3×10⁴ / 1.0×10⁴ K for n = 10⁻⁴ / 10⁻³ / 10⁻² cm⁻³.
  This model breaks the exact density scaling of pure collisional cooling.

### Caveat on the approximate CIE curve

The S&D93 tables could not be downloaded in the environment this was built in.
`cooling.py` builds Λ(T, Z) from the exact primordial H+He CIE cooling (Katz et al.
1996 rates) plus a smooth approximation to the S&D93 solar metal contribution scaled
by Z. Compared with the Cloudy CIE table it is ~0.1–0.3 dex too high between
2×10⁵ and 3×10⁶ K, and the detailed wiggles in dlnΛ/dlnT between
10⁵ and 10⁶ K come from the approximation. Swap in your own table with
`cooling.load_table` / `CoolingTable(T, L)`, or write it with `write_pluto_table`.
