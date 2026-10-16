# Two-photon TDSE campaign: two-photon ionisation and excitation of He and Ar with tRecX — results

For optmat's brief for the two-photon TDSE campaign (unpublished; with its 9 Oct addendum and the decisions in `STATUS.md`).
- Runs on a Slurm cluster on 9–11 Oct 2026; the run 5 checks for addendum 3 (run5b) on 11 Oct. This campaign is separate from round 2: the cluster's `$RUN_ROOT/`.
- Every run exited 0.
- Working log with every job, fix and check: `NOTES.md`.

## What is where

| item | file |
|---|---|
| Run 1 results (He, degenerate) | `run1/summary_run1.md`; every row `run1/results_sae.csv`, `run1/results_hacc.csv` |
| Runs 2–6 results | `summary_runs2to6.md`; every row `run2/results.csv` … `run6/results.csv` |
| Run 5: P(1s2s) against time | the cluster's `run5/hacc/<name>/pop_vs_t.csv` |
| Inputs as run | `run1/…`, `run2/…` … `run6/<model>/<name>/<name>.inp` (tRecX's own copy: `<name>/0000/inpc` on the cluster) |
| Pulse definitions | `make_inputs.py` (run 1), `make_runs.py` (`pulse()`, runs 2–6), `manifest_runs2to6.csv` (every run's colours, FWHM, chirp, peak A, window) |
| Pulse files | the cluster's `pulse/` (A(t)/max\|A\| on a 0.01 au grid; 1.2 GB, not copied); run 1's five files also in `pulse/` here |
| Propagated 1s2p (addendum 4) | `run5c/results.csv`, `run5c/scan_fit.csv`, `run5c/pop_<name>.csv` (1 au grid, with ⟨z⟩ and ⟨p_z⟩) |
| Run 5 checks (addendum 3) | `run5b/results.csv` (every run 5 and run5b run), `run5b/scan_fit.csv`, `run5b/pop_<name>.csv` (populations every 1 au) |
| Extraction | `extract_run1.py`, `extract_runs.py`, `extract_run5b.py`, `extract_run5c.py`; run 4's perturbative amplitudes `run4/pt_two_photon.py` |

**Code.**
- tRecX commit a987d7c (OccDual: bi-orthogonal bound-state populations), GCC 16 build.
- Run 5 used the same commit plus a patch that writes the energies of the projected states
  (`occE_patch_20261010.diff`, uncommitted). haCC has no eigenvalue calculation in this version, so this is
  the only way to get its level energies.

## Conventions used

- **Field and gauge.**
  - The field is defined on the **vector potential**, velocity gauge, linear polarisation along z:
    A(t) = Σᵢ (E0ᵢ/ωᵢ) g(t) sin(ωᵢt).
  - g = exp(−2 ln2 t²/T²), with T the **FWHM of the intensity**. CEP 0 for every colour.
  - E = −dA/dt. The peak field (one colour, transform-limited) equals E0 to ≤ 1e-4; two colours peak at E01 + E02.
  - Two colours share the envelope and the centre time.
- **Intensity.** I0 = ½ ε0 c E0² per colour (E0² × 3.50944758e16 W/cm²).
- **Chirp.**
  - A quadratic spectral phase φ(ω) = ½ φ₂ (ω − ω0)² on A(ω), in the convention A(t) = ∫ A(ω) e^{+iωt} dω with
    A(ω) → A(ω) e^{−iφ}. So **φ₂ > 0 is an up-chirp**, checked from the instantaneous frequency.
  - The spectrum, and so the fluence, is unchanged. **I0 is the transform-limited peak intensity.**
  - 3 fs with ±10 fs² stretches to 9.72 fs (peak field 0.556 E0); 1 fs with ±5 fs² to 13.9 fs (0.269 E0).
- **Window.** [−3T_eff, +3T_eff + 1000 au]; measures at the end of the window unless stated.
- **Energies:** each model's own; **no NIST values were used anywhere**.

| model | level | energy |
|---|---|---|
| One-electron He (Tong–Lin) | Iₚ | 24.5937 eV |
| | 1s2p | 21.056 eV |
| | 1s2s | 20.240 eV (a one-electron 2s, not 1s2s ¹S) |
| haCC He (He19, 1 channel) | Iₚ | 24.541 eV (from its initial energy) |
| | 1s2s | **20.739 eV in the propagation** (run 5 checks, below); 20.636 eV from its Occupation set |
| | 1s2p | **21.2089 eV in the propagation** (addendum 4); the Occupation set gives the same to 0.25 meV |
| One-electron Ar (Tong–Lin) | Iₚ | 15.774 eV |
| | 4s | 11.267 eV |
| | 4p | 12.926 eV |
| | 3d | 13.859 eV |

- The haCC Occupation energies are the states of its bound-state projection, in its production basis (lMax 10,
  rMax 40), relative to its ground state; NIST gives 20.616 eV (1s2s) and 21.218 eV (1s2p).
- **The Occupation set is not the propagated Hamiltonian's spectrum.** The run 5 checks put the 1s2s, as the
  propagation sees it, 102.6 meV above the Occupation value. The Occupation set's ground state (−0.901876 au) also
  differs from the propagated initial state (−0.901867 au). The 1s2p is not shifted (addendum 4): the Occupation 1s2p
  state is the propagated one to 1e-8 au. So the error is in the Occupation set's ¹S excited state. The cause has not
  been found.

## Discretisation

- **One electron:**
  - A2's radial grid: 3 elements of order 14 on [0, 2], 4 of order 12 on [2, 20], ~4.3-bohr elements of order 12 to
    BOX, then 30 polExp functions;
  - Eta 10 (l ≤ 9), m = 0 (and m = 1 for argon); ECS θ 0.15 (He) or 0.3 (Ar) from BOX; cutEnergy 100;
  - explicit RK4 with a fixed step: round 1's rule 2/Emax, capped at 0.03 (He) and 0.015 (Ar);
  - BOX 40 in run 1, BOX 80 in runs 2–6 (run 1 showed the BOX 40 spectra off by up to 8 %).
- **haCC:** He19_Standard, 1 channel (He⁺ 1s), lMax 10, rMax 40. haCC fixes cutEnergy 400 and its own ECS.
- **The ionised fraction:**
  - one electron: P_proj = 1 − Σ_bound OccDual (every state with Re E < 0 of the run's own discretisation), with the
    no-field drift of the same discretisation subtracted;
  - every run: the **iSurff photoelectron spectrum** from 0 to 40.8 eV on 1,200 points;
  - peak integrals over ±0.25 eV around each two-photon peak;
  - haCC: the norm-loss, which equals its spectrum integral to 0.1 %.
  - tRecX's default spectral range (max(12 Uₚ, 3Ω)) came out as ~3e-4 au for these weak fields, which made the first
    run 1 useless. It was rerun with the range set explicitly.
- **Bound populations:** OccDual (one electron; haCC in run 5).

## Convergence (details in the two summaries)

- **Run 1, one electron:**
  - Eta 15: 0.000 %; θ × 2/3: −0.3 %; dt × 2/3 and cutEnergy 150: +1.6 % (after their own drift); BOX 80/120: 0.0 %
    on the projection;
  - the spectrum at BOX ≥ 80 agrees with the projection to ≤ 2 %;
  - P ∝ I² to ≤ 1.5 % over 1e11–1e13.
- **Run 1, haCC:** lMax 15: +0.003 %; rMax 60: −0.23 %; P ∝ I² to ≤ 1.2 %.
- **Run 3, argon at 9 eV, 1e12, m = 0:** Eta 15 0.000 %; BOX 120 +0.18 % (peak); dt × 2/3 0.000 %; cutEnergy 150
  +0.02 %; θ × 2/3 0.000 %.
- **Run 4:** the TDSE 4p populations equal the perturbative prediction to 0.1 % at 1e10 (both m, 30 and 50 fs).
- **Drift** (one-electron no-field norm-loss over the window): He 1.4e-7 (5961 au, BOX 40), He 3.0e-8 (run 6's 1 fs
  window), Ar 2.4e-9.
  - So the one-electron He 1e11 points, and run 6's 1e12 points, rest on the spectrum integral, which is free of
    drift.
  - The drift grows with the window: run 6's chirped pulses have 3.6× longer windows than its no-field run, so use
    the spectrum there.

## Results in brief

**Run 1: helium σ⁽²⁾ (10⁻⁵² cm⁴ s), at 1e12.** σ⁽²⁾ = P / ∫Φ² dt, with Φ = I(t)/ħω and
I(t) = I0 exp(−4 ln2 t²/T²).

| photon energy (eV) | 13 | 14 | 16 | 18 | 20 |
|---|---|---|---|---|---|
| one electron | 7.46 | 7.33 | 6.52 | 6.23 | 9.35 |
| haCC | 11.89 | 12.80 | 10.72 | 9.06 | 10.18 |

**Run 2: helium two-colour.**
- Each colour 1e12, 20 fs, at Ω/2 = 14 and 17 eV (ν = 0.5, 1, 2, 4, 6) and 20 eV (ν = 0.5, 1, and 4 beyond the
  1s2p; ν = 2 dropped, 0.06 eV from the model's 1s2p).
- One run per point plus a one-colour reference at Ω/2 (the per-colour intensity), in both models.
- Reported: the integrals of the cross peak at Ω − Iₚ and of the 2ω₁ and 2ω₂ peaks.
- The ratio P_cross/P_ref includes the two-colour factor 4: haCC 4.00 (ν = 0.5) to 5.14 (ν = 6) at Ω/2 = 14 eV.

**Run 3: argon (one electron, atom = 2 P(m=0) + 4 P(m=1)).**
- Degenerate σ⁽²⁾ from 1.9e-50 (8 eV) to 3.1e-50 cm⁴ s (11 eV), plus the two-colour peaks at Ω/2 = 8.5, 9.5 and
  10.5 eV.
- **11 eV is 0.27 eV below the model's 4s** (11.267 eV): near-resonant in this model. Its m = 0 yield at 1e13 is 15 %
  below I² scaling.
- **At 8 eV** (E_kin 0.23 eV) the spectrum reads 11 % above the projection.

**Run 4: argon 3p → 4p (one electron), 6.463 eV.**
- Lowest-order amplitudes, length gauge:
  - **M₀ = +8.914 au, M₁ = +2.418 au** (M₁/M₀ = +0.2713);
  - R_s = +17.07 au, R_d = +12.09 au, with M₀ = R_s/3 + 4R_d/15 and M₁ = R_d/5.
- **Phase convention:** u = rR real and positive near the nucleus; spherical harmonics with the Condon–Shortley
  phase.
- The TDSE populations confirm |M_m|: TDSE/PT = 0.999 at 1e10 and 0.99 at 1e11, for both m and both durations.
- **Deviation from the agreed plan:** the signs come from this perturbative calculation in the same potential, not
  from an offline projection of ψ(T). tRecX writes ψ(T) and its eigenvectors only in its binary format, with an
  arbitrary eigensolver phase.

**Run 5: helium 1s2s (haCC).**
- 1s2s, 1s3s and 1s3d populations, the ionised fraction, and P(1s2s)(t), for 1, 3 and 10 fs and 3 fs ±10 fs², at
  1e12, 1e13 and 1e14.
- **Superseded by the checks below:** the ~16 % ±10 fs² difference at 1e12 was an end-of-window sampling of a beat;
  the post-pulse means agree to 0.07 %. The 10 fs shortfall against lowest order is haCC's 1s2s lying at 20.739 eV
  in the propagation, not at the 20.636 eV of its Occupation set.
- One haCC state at −0.0778 au relative to He⁺ is never populated and matches no helium singlet: unidentified.
- 1s3d and 1s3p are told apart by the two-photon selection rule: 1s3d is populated and 1s3p is not.

**Run 6: helium broadband (one electron), 14 eV.**
- Yields from the spectrum for 1, 2 and 5 fs and 1 fs ±5 fs², at 1e12 and 1e13.
- The ± chirps give equal yields to 0.01 %.

## Deviations from the brief, all agreed or stated

- **Run 2:** no P1/P2 runs (agreed: the peaks are resolved). The intensity is 1e12 per colour, not 1e11: it is
  perturbative (P ∝ I² to 1.5 % up to 1e13) and above the one-electron drift floor. 20 fs throughout; 10 fs was not
  needed.
- **Run 3:** one electron only (agreed). The optional haCC mMax 1 point (9.5 eV, 6 fs) has not been run (~3 days on
  a whole node); it is available on request.
- **Run 4:** signs from perturbation theory (above).
- **Gauge:** velocity gauge only. A length-gauge rerun was not made. Run 4's perturbative check is in the length
  gauge and agrees with the velocity-gauge TDSE to 0.1 %.

## Addendum 3: the run 5 checks (haCC, He 1s2s)

**Runs.**
- 11 haCC runs (`run5b/`, jobs 479500/479501, 11 Oct; all exit 0, 10–50 min each on 16 ranks), set up as run 5:
  He19, 1 channel, lMax 10, rMax 40, occE build.
- Populations printed every 1 au (run 5 printed every 10 au).
- Extraction: `extract_run5b.py` (sha256 b8b94831…). It also re-reads every run 5 run.

**Measures.**
- P(1s2s) = OccDual of Occupation state 1, averaged after the pulse (t > 2 T_eff, where the field envelope is below 4e-3
  of its peak).
- Lowest order: P = K |S2(Ω)|², with S2(Ω) = ∫ E₊(t)² e^{iΩt} dt computed from each run's own pulse file (dt 0.01 au).
  E₊ ~ e^{−iωt} is the positive-frequency part of the field. This is addendum 3's formula, with K = |M|²/4.

### Check 1: P(1s2s)(t) after the pulse, at 1 au

- **No decay.** The linear trend after the pulse is ≤ 0.3 % per 1000 au in every run5b run.
- **The oscillation is one beat, at 0.7621 au (period 8.245 au, 20.738 eV), in every field run.** It is the propagation's
  1s2s–ground spacing. Its relative rms is 1.4–8.4 %.
- **It is an interference with the projection's floor (check 4).**
  - The OccDual 1s2s projector picks up the ground state with |ε|² = 2.862e-9.
  - So P(t) = |c₁ + ε c₀ e^{iΔt}|², and its relative rms should be √2 |ε| / √P̄.
  - That predicts 1.89 %, 6.56 % and 8.37 % for the 10 fs reference, the −75 meV point and the 3 fs chirp; measured:
    1.89 %, 6.49 % and 8.36 %.
- **Consequence.**
  - Use the post-pulse mean. The cross term averages out, leaving a bias of +|ε|² = 2.9e-9: 0.02 % at 10 fs and 0.9 %
    at 1 fs, 1e12.
  - Any single sample is off by up to ±2√2 × rms.
  - Run 5's end values were single samples. They explain the 16 % chirp difference and most of the 1 fs excess
    (1 fs, 1e12: end 3.645e-7, mean 3.140e-7).

### Check 2: carrier scan, 10 fs, 1e12

| ω − 10.306 (meV) | 2ω (eV) | P(1s2s), mean after the pulse | at the end | TDSE / fitted lowest order |
|---|---|---|---|---|
| −75 | 20.4625 | 1.330e-6 | 1.208e-6 | 0.995 |
| −50 | 20.5125 | 3.754e-6 | 3.550e-6 | 0.997 |
| −25 | 20.5625 | 8.599e-6 | 8.292e-6 | 0.999 |
| 0 | 20.6125 | 1.599e-5 | 1.557e-5 | 1.000 |
| +25 | 20.6625 | 2.412e-5 | 2.361e-5 | 1.000 |
| +50 | 20.7125 | 2.955e-5 | 2.898e-5 | 1.000 |
| +75 | 20.7625 | 2.940e-5 | 2.883e-5 | 1.000 |

- **Fitting P_i = K |S2_i(Ω)|² with K and Ω free gives Ω = 20.7386 eV** (relative residual 1.2e-7) and **K = 3.078 au**.
  - The residual rises to ~5e-4 at ±5 meV, so the position is fixed to a few meV.
  - Fitted Ω is 102.6 meV above the Occupation set's 1s2s (20.6360 eV).
  - The beat in check 1 (20.738 eV) is the same spacing, measured independently.
  - The resonance is at ω = 10.369 eV per photon.
- **4K = 12.31 au².** Addendum 3 quotes Hattig's strength δ = 12.306 au. If δ is |M|² in that formula, haCC agrees to
  0.04 % once its own resonance position is used.

**All run 5 runs against lowest order at the fitted resonance** (post-pulse means; K = 3.078 from the scan):

| pulse | 1e12 | 1e13 | 1e14 |
|---|---|---|---|
| 1 fs | 1.031 | 1.013 | 0.952 |
| 3 fs | 1.009 | 0.983 | 0.748 |
| 10 fs | 1.001 | 0.858 | 0.079 |
| 3 fs, +10 fs² | 1.006 | 0.989 | 0.839 |
| 3 fs, −10 fs² | 1.010 | 0.983 | 0.788 |

- At 1e12, lowest order holds to ≤ 1 % for every duration and chirp, and to 3 % at 1 fs. At 1 fs, the floor accounts
  for 0.9 % of that 3 %.
- Against the Occupation set's 1s2s instead, the 10 fs, 1e12 ratio is 0.539, reproducing addendum 3's 0.52.
- At 1e13 and 1e14 the departures are nonlinear and grow with duration.

### Check 3: the 3 fs ±10 fs² pair, accuracy 1e-9, +3000 au after the pulse

- **Post-pulse means: 8.1659e-7 (+10 fs²) and 8.1716e-7 (−10 fs²), 0.07 % apart.** Lowest order predicts them equal;
  they are equal.
  - The end values (7.733e-7, 8.364e-7) differ by 8 %, from the beat alone.
- Run 5 at the default accuracy (1e-7):
  - post-pulse means 8.159e-7 and 8.192e-7 (0.4 % apart);
  - 10 fs reference at accuracy 1e-9 against 1e-7: 1.5989e-5 against 1.6002e-5 (0.08 %).
  - Accuracy 1e-7 was adequate.
- The ± difference does grow with intensity: 0.6 % at 1e13 and 6.6 % at 1e14 (6.81e-3 vs 6.39e-3). It is real but
  beyond lowest order.

### Check 4: the OccDual floor, no field, the 1 fs window

- **P(1s2s) = 2.862e-9**, equal at the start and end (to 1e-6 relative). Over the window, norm-loss 1.5e-10 and
  P(ground) 0.99999998.
- This run printed only its first and last rows: with no laser, tRecX takes one step. So the floor's time dependence
  was not sampled. Its size is confirmed by the beat amplitudes in check 1.
- It is 0.9 % of the 1 fs, 1e12 population (3.14e-7).

### For optmat

- **Use for haCC's 1s2s:**
  - resonance 2ω = 20.739 eV;
  - P(1s2s) = post-pulse mean − 2.9e-9 (`run5b/results.csv`: `P_mean`, `K_fit`, `ratio_fit`, for run 5 too).
- With these, run 5 is lowest-order to ≤ 1 % at 1e12 for 3 and 10 fs and both chirps, with 4K = 12.31 au².
- **haCC's Occupation energies are not the propagated spectrum for the 1s2s.** The 1s2p was measured in addendum 4
  (below): it is unshifted, so run 2's pole distances stand.

## Addendum 4: haCC's propagated 1s2p

**Runs.**
- 5 haCC runs (`run5c/`, jobs 479626/479627, 11–12 Oct; all exit 0, 1.9–2.0 h each on 16 ranks), set up as run 5b:
  He19, 1 channel, lMax 10, rMax 40, occE build, velocity gauge, accuracy 1e-9.
- Pulse: one colour, 10 fs, 1e9 W/cm², carriers 21.21 eV + {0, ±50, ±100} meV (0 is the main run).
- Window: +3000 au after the pulse; populations on a 1 au grid.
- Extraction: `extract_run5c.py` (sha256 3f24eaba…); `run5c/results.csv` (c604f768…), `run5c/scan_fit.csv`
  (ea64b98b…).

**Change of method.** The brief's P(1s2p) beat cannot be measured.
- The 1s2p projector is ¹P, so by parity it cannot see the ¹S ground state. The no-field floor of P(1s2p) is
  **1.27e-19** (run5b), against 2.862e-9 for the 1s2s.
- A leak beat would be ~1e-8 of P(1s2p). After the pulse, P(1s2p) is flat to a relative rms of 1e-4, with no beat.
- **Instead, the dipole was printed** (`Dipole: length, velocity`, i.e. ⟨z⟩ and ⟨p_z⟩).
  - After the pulse, ⟨z⟩(t) = 2|c₀c₁z₀₁| cos(Δt + φ), with Δ the propagated 1s2p–ground spacing and amplitude
    1.2–1.8e-2 au.
  - Fitted by a cos + b sin + c, after t > 2T.
- The carrier scan was run alongside as the second, independent method.

**Results.**

| measure | value |
|---|---|
| **beat in ⟨z⟩, length and velocity forms, 5 runs** | **Δ = 0.7794133 au = 21.208918 eV**: spread 0.9 µeV over the 10 fits, fit σ 0.03 µeV each |
| **scan fit, P = K\|S₁(Ω)\|²** | **Ω = 21.20892 eV** (0.05 meV grid; relative residual 3.4e-10) |
| Occupation set's 1s2p − ground | 21.20917 eV |
| propagated 1s2p, absolute (E_init + Δ) | −0.12245328 au, against the Occupation 1s2p −0.12245329 au: **equal to 1e-8 au** |
| K | **0.169939 au²** (K = \|z₀₁\|² with E = E₊ + E₋) |
| d_z | **0.41224 au** from the scan, 0.412236 au from the beat amplitude (with P₀, P₁ from OccDual) |
| d_z against NIST's 0.4208 au | −2.0 %; f = 2Δ\|z\|² = 0.2649 against 0.2762 (−4.1 %) |
| P(1s2p), post-pulse mean, at 21.21 eV | 4.688e-4 (end 4.688e-4); ±50 meV 3.86e-4 / 3.76e-4; ±100 meV 2.10e-4 / 1.98e-4 |
| first order | TDSE / K\|S₁\|² = 1.0000 ± 5e-5 at all five carriers |
| norm-loss | 6.9e-10 |
| state 3 | at most 1.3e-17: not populated, no beat |

- **The propagated 1s2p is 21.2089 eV, ±1 µeV statistically.**
  - It is 0.25 meV below the Occupation set's spacing. All of that difference comes from the ground state: the
    Occupation ground state is 9.15e-6 au (0.249 meV) below the propagated initial state.
  - So **the 1s2p is not shifted**, and run 2's pole distances (Occupation 1s2p 21.209 eV) stand.
- **The 1s2s, by contrast, is shifted.**
  - Its propagated position (20.739 eV, addendum 3; beat resolution 1e-4 au, about 3 meV) puts the propagated 1s2s at
    −0.1398 au, against the Occupation set's −0.1435 au.
  - So the Occupation error is in its ¹S excited state, not a common offset. The cause has not been found.
- **Gauge check.** The velocity/length amplitude ratio is 0.82571, against Δ = 0.77941 au for exact eigenstates of one
  Hamiltonian: +5.9 %.
  - The two dipole forms are not equivalent in this 1-channel haCC basis.
  - The frequency is the same in both forms.
  - The d_z above is the length form; the scan's K is gauge-free.
