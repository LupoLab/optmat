# Ar 800 nm / 6 fs Gaussian intensity scan — TDSE vs PPT

Sibling of `../he_intensity_scan` for **argon** (3p, Tong–Lin SAE, Ip = 15.76 eV). Same method,
measures, scripts and PPT comparison; see that README for the metric/numerics rationale. This file
records only the Ar-specific choices and the (different!) result.

```
uv run make_inputs.py     # ar_scan_00..09.inp (gauss) + manifest.csv
uv run make_isurff.py     # ar_isurff_00..06.inp (cos8 + iSurff)
./tRecX ar_scan_NN.inp ; ./tRecX ar_isurff_NN.inp
julia scan_ppt.jl         # -> results_ar_intensity_scan{,_loglog}.png, results_ar_scan_{zoom,relerror}.png
```

## Intensity scale — matched to He by the barrier parameter (the sensible choice)
The interesting variable for a total-yield-vs-intensity curve is the **barrier-suppression parameter**
u = E₀/E_b (E_b = Ip²/4Z): it sets where the curve saturates. Matching u reproduces the *same
yield-curve shape*. Since E_b ∝ Ip², this is **I ∝ Ip⁴**, so

  I_Ar = (Ip_Ar/Ip_He)⁴ · I_He = (15.76/24.59)⁴ = **0.169 × I_He  →  1.7×10¹³ … 8.4×10¹⁴ W/cm²**,

giving the same **u-range [0.26, 1.85]** (multiphoton tail → over-barrier saturation) as He. The
Keldysh γ range differs (**2.8 → 0.40**, vs He's 1.43 → 0.20): Ar reaches *deeper* into the
multiphoton regime at the low end — matching u and matching γ are not the same because γ ∝ √Ip/E₀
while u ∝ E₀/Ip². (Matching γ instead would push the top to ~3×10¹⁵ where Ar is u≈3.6, saturated
across most of the window — less informative.)

## Numerics (vs He)
- **Deep SAE core** (V → −18/r at small r) sets RK stability → **fixStep = 0.015** (vs He's
  laser-limited step). Initial state = 3p valence = propagation `state=5` (eig-index 4 + ECS offset);
  confirmed ⟨H₀⟩ = −0.5797 au.
- **No big-box problem**: quiver α = E₀/ω² ≤ 48 a₀ < BOX = 80 over the whole range, so the He
  box>quiver fix isn't needed — uniform **BOX = 80, standard ECS theta = 0.3**. Eta = 30, cutEnergy = 60.

## Result
| I (W/cm²) | γ | u | chi-loss | norm-loss | iSurff | PPT plain | PPT +T-L | PPT +Zhang |
|-----------|------|------|------|------|------|------|------|------|
| 4.03e13 | 1.81 | 0.40 | 0.0007% | 0.0005% | 0.0005% | 0.0003% | 0.0002% | 0.0003% |
| 6.22e13 | 1.46 | 0.50 | 0.012%  | 0.011%  | 0.0096% | 0.0079% | 0.0048% | 0.0056% |
| 9.61e13 | 1.17 | 0.62 | 0.170%  | 0.147%  | 0.137%  | 0.142%  | 0.078%  | 0.087%  |
| 1.48e14 | 0.94 | 0.78 | 1.60%   | 1.49%   | 1.38%   | 1.80%   | 0.86%   | 0.93%   |
| 2.29e14 | 0.76 | 0.96 | 10.5%   | 10.2%   | 9.68%   | 15.2%   | 6.40%   | 6.69%   |
| 3.54e14 | 0.61 | 1.20 | 40.1%   | 38.3%   | —       | 65.7%   | 29.6%   | 30.1%   |
| 5.46e14 | 0.49 | 1.49 | 82.4%   | 79.9%   | —       | 99.4%   | 73.2%   | 73.3%   |
| 8.44e14 | 0.40 | 1.85 | 98.8%   | 97.3%   | —       | 100.0%  | 97.3%   | 97.3%   |

**Conclusion — the opposite of He.** For He the Tong–Lin/Zhang BSI-corrected PPT *tracked* the TDSE
within ±20%. For **Ar, neither PPT flavour tracks it — the TDSE is bracketed between them:**
- **plain PPT over-predicts** by up to **+57%** in the tunnelling region (and *under*-predicts at the
  lowest, deeply multiphoton points where tunnelling theory fails);
- **Tong–Lin / Zhang BSI over-correct → under-predict by 25–58%** across most of the range (the Ar BSI
  coefficient α=9 is more aggressive than He's α=7), converging only at saturation.

So for Ar the true (iSurff) yield sits *between* plain and BSI-corrected PPT — the same bracketing seen
in the structured-field Ar study, and a clean contrast to the He scan. The chi-vs-iSurff **excitation
gap** is +41% → +9% (larger than He at the low end: more multiphoton Rydberg excitation at γ≈1.8).

Figures: `results_ar_intensity_scan{,_loglog}.png`, `results_ar_scan_zoom.png`,
`results_ar_scan_relerror.png` (the clearest view of the bracketing).

## m=1 orbital and the full 3p⁶ shell  (`make_m1.py`, `scan_m_combined.jl`)
The scan above is the single **m=0** 3p orbital (aligned with the field). The 3p⁶ shell also has the
two **m=±1** orbitals (perpendicular), which ionise much less. We rerun the whole scan for m=1 and
combine into the shell ionisation.

**m=1 setup:** `Phi,1,0.,6.28318530718,expIm[1]` selects a single m=1 azimuthal function;
`assocLegendre{Phi}` then builds P_l¹ (l≥1). The m=1 manifold is 2p,3p,3d,… (no s): 3p(m=1) is
eig-index 1 → `Initial: state=2` (⟨H₀⟩=−0.5797, same 3p energy as m=0). Everything else identical.

**Shell ionisation** (probability that ≥1 of the 6 electrons leaves), treating the orbitals as
independent with spin degeneracy: **P = 1 − p₀² p₁⁴** (p₀ = m=0 survival, p₁ = m=±1 survival;
2 electrons in m=0, 4 in m=±1). Computed for chi[0,20], norm and iSurff survivals, and compared to
**PPT with `msum=true, occupancy=2`** (the full shell: m=0 + m=±1 channels × 2 spins).

| I (W/cm²) | 1−χ (m=0) | 1−χ (m=1) | shell 1−p₀²p₁⁴ (χ) | shell (iSurff) | PPT msum,occ2 |
|-----------|------|------|------|------|------|
| 9.61e13 | 0.170% | 0.0042% | 0.356% | 0.289% | 0.297% |
| 1.48e14 | 1.60%  | 0.049%  | 3.36%  | 2.91%  | 3.76%  |
| 2.29e14 | 10.5%  | 0.43%   | 21.2%  | 19.8%  | 29.6%  |
| 3.54e14 | 40.1%  | 2.75%   | 67.9%  | —      | 90.0%  |
| 5.46e14 | 82.4%  | 12.7%   | 98.2%  | —      | 100.0% |
| 8.44e14 | 98.8%  | 38.5%   | 99.998% | —     | 100.0% |

**Findings:**
- **m=1 ionises far less than m=0** — suppressed by ~3–4× in intensity (e.g. at 1.5e14, 0.049% vs
  1.60%; at the top, 38% vs 99%). The perpendicular 3p lobes barely feel the field.
- The **shell yield exceeds single m=0** (it counts *any* of the 6 electrons; dominated by the two
  m=0 electrons, ≈ 1−p₀²), so it saturates earlier than the single-orbital curve.
- **PPT(msum,occ=2) tracks the TDSE shell very well at low intensity (multiphoton) and at
  saturation**, but **over-predicts the tunnelling rise** (29.6% vs 21.2% at 2.3e14; 90% vs 68% at
  3.5e14) — the same plain-PPT over-prediction as the single-orbital m=0 comparison.

Figures: `results_ar_shell.png` (linear), `results_ar_shell_loglog.png`.
