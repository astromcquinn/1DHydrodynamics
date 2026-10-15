# Literature notes for the cooling-perturbation draft

Key points of each cited paper as they bear on our results. Notes marked **(R)**
were read in full by the drafting session. The others were summarized by helper
agents from the uploaded PDFs and spot-checked. Equation numbers refer to the
papers. Quotes are verbatim.

## Pure cooling (no heating): direct precedents
- **Schwarz, McCray & Stein 1972, ApJ 175, 673 (R).** The closest precedent.
  - Setup: interstellar gas cooling at constant density from 8200 K, no heating.
  - Method: linearized equations whose coefficients change in time (§IV), integrated numerically, plus 1D simulations.
  - Short-wavelength law: eq. 26 gives d ln n1/dt = (2 − s_T)/(γ τ_c); eqs. 27 and 30 integrate it to n1 ∝ (T0/T)^((2−s_T)/γ).
  - Heated comparison: eq. 29 gives the heated (FGH, per-mass) rate (1 − s_T)/(γ τ_c).
  - Long wavelengths keep growing "long after the instability shuts off" (inertia).
  - Linear growth is ~1e2, nearly independent of λ for 1–10 pc.
  - Threshold (note added in proof): "an initial perturbation of minimum amplitude (between 5 and 10 percent) is needed" to produce a contrast above 100:1, which they call "a key difference between instability in a medium in equilibrium, and instability in a cooling medium."
- **Burkert & Lin 2000, ApJ 537, 270 (R).**
  - Setup: isochoric cooling, no heating, Λ ∝ T^β, linear theory plus 1D hydro.
  - Isobaric growth (eq. 22): ρ1/ρ0 ∝ (T0/T)^((2−β)/Γ). We verified numerically that this equals our time-shift law.
  - The isobaric → isochoric transition occurs when Q = τ_c k c_iso drops below 1 (eqs. 20, 24, 25).
  - After the transition the perturbation grows linearly in time from its velocity (eq. 27).
  - For K ≤ 1 the early solution is eqs. 13–14: V ≈ (β−2)K²A τ²/2. They expected the reversed pressure gradient to erase the fluctuation after cooling.
  - Critical K for going nonlinear before the transition: eq. 28, K_crit = A^((2β−3)Γ/(4−2β)) (316 for A = 0.01, β = 0).
  - Thermal conduction erases λ < λ_κ = 2π λ_F/(2−β)^(1/2) (eq. 35). Heating gives a two-phase medium (§6).
  - Conversion to our variables: K = 2π/(√γ x).
- **Mathews & Bregman 1978, ApJ 224, 308 (R, appendix).**
  - Short-wavelength growth along compressing cooling-flow parcels (appendix eqs. A6–A8).
  - Growth factors e^Γ ≈ 1e4±1.
  - They argue that frozen-in (entropy) perturbations, not travelling waves, produce condensations.
- **Field 1965, ApJ 142, 531 (§VI read).**
  - A static medium with no heat input: relative growth if N < 2.
  - With heating per unit mass independent of ρ and T: N < 1. Hence "a system which is cooling in the unperturbed state is more unstable than one in thermal equilibrium".
  - The sound-wave damping by radiative relaxation in his dispersion relation scales as k² at long wavelength.
- **Balbus 1986, ApJ 303, L79.**
  - Lagrangian criterion with ℒ → ℒ/T.
  - Eq. 23 gives the short-wavelength growth on any homogeneous time-dependent background. For ℒ = n²Λ it reduces to our per-e-fold law.
- **Koyama & Inutsuka 2000, ApJ 532, 980.** Thermal instability in a cooling post-shock layer, with heating. Its linear analysis freezes the coefficients, valid only for t ≪ γτ_net.
- **Hennebelle & Pérault 1999, A&A 351, 309; 2000, A&A 359, 1124.**
  - Heated, bistable interstellar gas: condensation needs a roughly transonic converging flow on scales ≳ c_s t_rad.
  - Self-similar pure-cooling solutions.
  - A transverse magnetic field limits compression.

## Heated media / haloes
- **Das, Choudhury & Sharma 2021, MNRAS 502, 4935 (R, §2–3).**
  - Uniform medium with heating per unit volume that depends on neither ρ nor T (δq+ = 0). That is why our per-e-fold equivalence is exact for their setup.
  - Nonlinear runs use δρ/ρ = 1e-5; they do not map a threshold amplitude.
- **McCourt et al. 2012, MNRAS 419, 3319; Sharma et al. 2012, MNRAS 420, 3174.**
  - Stratified, heated haloes: condensation where t_TI/t_ff ≲ 10 (or ≲ 1 plane-parallel), with buoyancy saturating growth.
  - Unheated cooling flows need large seeds (Sharma appendix A, eq. A6; McCourt §5.1).
- **Choudhury, Sharma & Quataert 2019, MNRAS 488, 3195.** The condensation threshold in t_cool/t_ff rises with seed amplitude. Blob criterion (eq. 14): (t_cool/t_ff)_blob ≲ 10.
- **Mohapatra et al. 2023, MNRAS 525, 3831.** Criterion min[t_ti/min(t_mix, t_ff)] ≲ 1.8 exp(6σ_s). Isochoric gas at 10^5.5–10^6 K, where t_cool < t_cs.
- **Balbus & Soker 1989, ApJ 341, 611.**
  - In cooling flows, buoyancy and advection limit linear growth to modest power laws (×2–3 for oscillatory modes).
  - Eulerian frozen-background growth rates miss amplitude changes.
- **Waters & Proga 2019, ApJ 875, 158.**
  - Classification of non-isobaric linear regimes for heated photoionized gas, including gas off the equilibrium curve.
  - Damped oscillations of large clouds settling onto the equilibrium curve.
  - Non-isobaric behaviour appears only in the nonlinear regime (contrast with our δT/T result).
- **Pizzolato & Soker 2005, ApJ 632, 821.** "Cold feedback": blobs with δρ/ρ ≳ 1–2 cool isobarically against drag and buoyancy.
- **Sánchez-Salcedo, Vázquez-Semadeni & Gazol 2002, ApJ 577, 768.**
  - Nonlinear thermal instability in heated interstellar gas.
  - Small-scale velocity forcing (rms Mach number ≳ 0.3) suppresses condensation.
  - §5.1 contrasts their heated case with Burkert & Lin's unheated one.
- **Gronke & Oh 2020, MNRAS 494, L27.**
  - Shattering requires a final contrast χ_f ≳ 300 and δP/P ~ 1. It happens even for initially linear contrasts (χ_i = 1.01).
  - It has the same loss of sonic contact as our small clouds.

## Cooling flows in galaxy haloes
- **Stern, Fielding, Faucher-Giguère & Quataert 2020, MNRAS 492, 6042 (R, §2–3).**
  This is Paper II; Paper I is Stern et al. 2019, MNRAS 488, 2549 (not read).
  - Steady CGM cooling flows using the Wiersma et al. tables with the HM12 UV background.
  - In the subsonic part, "radiative cooling is balanced via heating by compression
    as the gas flows inwards", and T ≈ T_vir.
  - The sonic point is at t_cool/t_ff ≈ 0.7. Inside it the cooling runs away: T falls
    by more than 10× over about 2× in radius.
  - 3D Athena run seeded with 3% isobaric perturbations: the perturbations grow
    significantly only inside the sonic radius, where t_cool ≪ r/v_r.
  - Cited in the draft's applicability subsection (`sec:gravity`).

## Not yet read (possibly relevant)
Goldsmith 1970, 1972; Defouw 1970; David & Bregman 1989; Kovalenko & Shchekinov
1999; Meerson & Sasorov 1987; McCray, Stein & Kafatos 1975; Loewenstein 1989.
