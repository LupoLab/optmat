# Two-photon TDSE campaign: working notes

## 9 Oct 2026: set-up and run 1

- `make_inputs.py pulse` -> `pulse/gauss20fs_<E>eV_A.dat` (13, 14, 16, 18, 20 eV): A(t)/A0 = exp(-2 ln2 t^2/T^2)
  sin(omega t), T = 20 fs (FWHM of the intensity) = 826.8275 au, CEP 0, on [-3T, 3T] at dt 0.01 au;
  `pulse_check.txt`: 63-97 cycles in the FWHM, max |E|/E0 0.99999, bandwidth 91 meV. Window [-3T, +3T + 1000 au].
  Intensity I0 = (1/2) eps0 c E0^2 = E0^2 x 3.50944758e16 W/cm2 (the brief's convention); A0 = E0/omega; velocity
  gauge.
- `make_inputs.py run1`:
  - sae: Tong-Lin He, BOX 40 (A2's radial grid), Eta 10, ECS 0.15, cutEnergy 100, OccDual 0:70, iSurff;
    15 points (5 energies x 1e11, 1e12, 1e13 W/cm2); convergence at 16 eV, 1e12: Eta 15, BOX 80, dt x 2/3,
    cutEnergy 150, theta 0.1; a no-field run; field-free spectra Eta 1-3 (the code's own energies).
  - hacc: He19_Standard, 1 channel, lMax 10, rMax 40, iSurff channel spectra; 15 points; lMax 15 and rMax 60 at 16 eV,
    1e12; a field-free spectrum (Eigen SmallReal[12]).
- Cluster: `run_scan_pack.slurm` (sae, 24 cores, 2-day limit) and `run_mpi_seq.slurm` (new: haCC runs one after another
  in one 16-rank job, with the abort watchdog), so that the haCC set uses two job slots, not 18.
- Submitted 478286 (sae), 478287/478288 (haCC); all queued on QOSMaxJobsPerUserLimit (15 running jobs per user,
  shared with the other agent's jobs, 12 of 15 at submission).

## 9 Oct: optmat's reply (relayed by John)

See STATUS.md "Plan" for the decisions. Points raised back: the two-colour peak spacing is ν, not ≥ 1 eV; what the
ν = 0 reference is (one colour at Ω/2 with the per-colour intensity, vs two identical colours = 4× intensity);
complex 4p amplitudes need ψ(T) and the eigenvectors written and projected offline (tRecX prints populations only);
the atom combination 2 P(m=0) + 4 P(m=1); haCC bound projection (Occupation) for run 5 to be tested.

## 9 Oct: two problems in the first haCC jobs

- `run_mpi_seq.slurm` bug (mine): mpirun reads stdin and swallowed the rest of the `while read` list, so 478287 ran only
  its first input and ended (30 s). Fixed (`< /dev/null`); resubmitted the 9 unstarted runs as 478302
  (`run1/list_hacc_a2`). 478288 (old script) will stop after its first run; follow-on 478303-ish queued with
  `--dependency=afterany:478288` on `run1/list_hacc_b` (skips .done).
- **haCC cannot compute eigenvalues in this tRecX:** `run_trecx.cpp:798` DEVABORT "at the moment, cannot compute full
  eigenvalues for HaCC" whenever `Eigen: select` is set. So: no haCC field-free spectrum (its Ip still comes from the
  initial energy -<H>), and probably no haCC bound projection (Occupation uses the same full eigensolver) — which
  run 5 (1s2s populations in haCC) needs. To confirm and report.
- The brief's addendum says the two-colour peak lies "2 nu >= 1 eV" from the 2 omega peaks; it is nu (0.5 eV at
  nu = 0.5). Still resolved at 20 fs.
- **haCC + Occupation works** (test `tests/hacc_occ/`, He19 1 channel lMax 4 rMax 30, 950 coefficients, 4 ranks): the
  Occupation path diagonalises directly (not through the Eigen branch that aborts). OccDual{H0+EnergyShift:0} = 0.99997
  under a weak 120.3 nm pulse (correct); Occ{} reads 3.99999 for the ground state (its normalisation is wrong here;
  use OccDual, as in round 2). The state energies are not printed, so the 1s2s index must be identified otherwise
  (by its resonant population in run 5, or by a patch that prints them). Run 5 in haCC is therefore feasible.

## 9 Oct ~14:15: run 1 first look, and a rerun with the spectral range fixed

- First one-electron run 1 (now `run1_v1_spec_kmax0p02/`, 24/24 done in ~10 min each), `extract_run1.py sae`:
  - P_ion (1 - sum of OccDual over the 68 states with Re E < 0) equals the norm-loss to all digits.
  - **Drift:** no-field norm-loss 1.432e-7 over 5961 au. That is the 1e11 signal's size (all 1e11 rows ~1.6e-7), and
    6–9 % of 1e12 (so dt × 2/3 and cutEnergy 150 "move" 1e12 by -6 %: different drift).
  - Drift-corrected, P ∝ I² between 1e12 and 1e13 to 1.5 % (16 eV: (1.638 - 0.143)e-6 × 100 = 1.495e-4 vs 1.517e-4).
  - sigma2 from 1e13 (drift negligible): 7.5e-52 (13 eV), 7.4e-52 (14), 6.6e-52 (16), 6.2e-52 (18), 9.2e-52 (20 eV)
    cm^4 s. Ip 24.5937 eV (Tong-Lin He -0.9038011 au).
  - **The iSurff spectra were useless:** tRecX's default range, max(12 Up, 3 Omega), came out as 2.5e-4 au (k <= 0.022)
    for these weak fields, far below the two-photon peaks (k 0.3–1.1). Spectrum integrals ~1e-15.
- The haCC runs (~95 min each on 16 ranks) had the same default (maxEnergy 3.7e-4): cancelled 478288, 478302, 478303
  (mine) after 1–2 runs.
- Fixed: `Spectrum: radialPoints=1200, maxEnergy=1.5` (40.8 eV, past the three-photon peaks; dE ~ 11 meV at the
  slowest peak). Everything else in the inputs unchanged (diffed). Old outputs kept as `run1_v1_spec_kmax0p02/`
  (cluster and workstation). One-electron rerun 478310 (28 inputs incl. the four eigrun_* spectra). haCC to follow once a
  one-electron spectrum confirms the fix (`run1/list_hacc_noeig`, 17 runs; the haCC Eigen input dropped: haCC
  cannot do Eigen).
- Spectrum fix confirmed (he_16eV_I1e12, rerun): the two-photon peak at E = 7.40 eV (2 omega - Ip = 7.41 eV), integral
  1.5785e-6 (the drift-corrected projection gave 1.495e-6: +5.6 %, to look at with the box/resolution variants);
  three-photon peak 3.1e-11. haCC rerun submitted: 478347, 478348.

## 9 Oct: optmat's second reply (relayed by John)

- Peak spacing is ν (their error); keep 20 fs for run 2; if haCC forces 10 fs, run ν = 0.5 eV at 20 fs anyway.
- ν = 0 reference: one colour at Ω/2, same envelope, per-colour intensity; they apply the factor 4 themselves.
- Agreed: Ar = 2 P(m=0) + 4 P(m=1); run 4 offline projection with the phase convention and the perturbative sign
  check; run 5 through OccDual; 1e11 from the spectrum integral, drift subtracted explicitly at 1e12.
- They ask: say whether the helium level energies quoted come from the one-electron model or from NIST.

## 9 Oct ~15:30: one-electron run 1 (spectral range fixed) and follow-up checks

- `run1/results_sae.csv`. sigma2 from the two-photon peak of the iSurff spectrum (1e11 and 1e12; 1e13 within -1.5 %):
  6.90 (13 eV), 7.99 (14), 6.89 (16), 6.31 (18), 9.55 (20) x 1e-52 cm^4 s; P ∝ I² to <= 1.5 %.
- The drift-corrected projection (1 - bound - drift) differs from the peak integral by -11 % to +9 % (energy-dependent,
  the same at every intensity, so not the drift). At 16 eV the spectrum changes -4.4 % from BOX 40 to 80; the
  projection does not. dt x 2/3 and cutEnergy 150 change the raw projection by -7 % but not the spectrum: their drift
  differs from the base's, so the drift must be subtracted per variant.
- Submitted 478367 (18 cores): BOX 80 and 120 at every energy (1e12), no-field runs for every variant
  (eta15, box80, dt23, cut150, th23, box120), eigrun_box120.
- Energies quoted: one-electron = the model's own spectrum (eig_l*/eigrun_*); haCC = Ip from its initial energy only;
  NIST values are not used anywhere in the runs.

## 10 Oct 23:20: run 1 complete

All 44 one-electron inputs and 17 haCC runs exit 0. `extract_run1.py sae|hacc` -> run1/results_*.csv; summary in
`run1/summary_run1.md`. The one-electron spectrum at BOX 40 is off by up to ±8 %; BOX 80/120 spectra agree with the
drift-subtracted projection to <= 2 %. haCC's norm-loss and spectrum agree to 0.1 %; haCC converged (lMax 15, rMax 60).

## 10 Oct ~23:20: runs 2-6 set up (`make_runs.py`, `manifest_runs2to6.csv`)

- **Pulses** (`make_runs.py pulse()`): A(t) = sum_i (E0_i/omega_i) g(t) sin(omega_i t), CEP 0 per colour, 20 fs unless
  stated; chirp as a quadratic spectral phase on A (phi2 > 0 = up-chirp, checked from the instantaneous frequency;
  3 fs ±10 fs^2 -> 9.72 fs, peak field 0.556 E0; 1 fs ±5 fs^2 -> 13.9 fs, 0.269 E0), I0 the transform-limited peak.
  Generated on the cluster (1.1 GB; not copied to the workstation). A at both file ends ~3e-6 of its peak.
- **Run 2** (He two-colour): Omega/2 = 14 (nu 0.5, 1, 2, 4, 6), 17 (the same), 20 eV (nu 0.5, 1, and 4 beyond the
  model's 1s2p at 21.056 eV; nu 2 dropped: omega1 = 21.0 within 0.06 eV), each colour 1e12 W/cm2, 20 fs; references:
  one colour at Omega/2, 1e12. One electron (BOX 80) and haCC (lMax 10, rMax 40). 16 points per model.
- **Run 3** (Ar, one electron, m = 0 and 1): degenerate 8, 8.5, 9, 10, 11 eV at 1e11/1e12/1e13 (11 eV is 0.27 eV from the
  model's 4s at 11.267 eV: near-resonant in this model); two-colour Omega/2 = 8.5 and 9.5 (nu 0.5-2), 10.5 (nu 0.5 only:
  omega1 must stay 0.3 eV below the model's 4s), 1e12 per colour; references at 9.5, 10.5; convergence at 9 eV 1e12 m = 0
  (Eta 15, BOX 120, dt x 2/3, cutEnergy 150, theta 0.2) with the drift of each; spectra for the labels.
- **Run 4** (Ar 3p -> 4p, one electron): 6.463 eV (half the model's 4p at 12.926 eV), 30 and 50 fs, 1e10 and 1e11, m = 0 and
  1. tRecX gives populations only; the signed amplitudes come from A_0 = R_s/3 + 4 R_d/15, A_1 = R_d/5 with |A_m| from the
  TDSE and the signs of R_s, R_d from lowest-order perturbation theory in the same potential (to write).
- **Run 5** (He 1s2s, haCC): 120.3 nm, 1/3/10 fs and 3 fs ±10 fs^2, 1e12/1e13/1e14, Occupation 0:20 (OccDual); needs the
  occE build to label the states (1s3p and 1s3d nearly degenerate).
- **Run 6** (He broadband, one electron): 14 eV, 1/2/5 fs and 1 fs ±5 fs^2, 1e12 and 1e13, and a no-field run.
- Uploads through the relay host: rsync and scp twice reported success without writing the files; piping through ssh
  (`ssh <relay> 'cat > f' < f`) worked. Verify every upload by md5.
- **Run 4 perturbative check** (`run4/pt_two_photon.py`, Dalgarno-Lewis on a radial grid, Tong-Lin Ar): dr 0.004 / 0.002:
  R_s +17.040 / +17.069, R_d +12.084 / +12.092 au; M_0 +8.902 / +8.914, M_1 +2.417 / +2.418 au; M_1/M_0 +0.2715 / +0.2713.
  Step-extrapolated E_3p, E_4p = -0.579689, -0.104673 au (tRecX's -0.579689, -0.104674). Convention: u = rR real, positive
  near the origin; Condon-Shortley; length gauge. Both radial sums positive: the relative sign of M_0 and M_1 is +.
- **occE build** (`tRecX-occdual-occE-gcc16{,-mpi}`, a987d7c + `occE_patch_20261010.diff`, no commit):
  first version printed the energies with PrintOutput::message, which is silent below verbosity "full" (built
  479250, nothing printed); second version (sha256 of the diff 3aa4093a…) writes `<output>/occupation_energies`
  (master rank): state, Re E, Im E of every Occupation state. Rebuild 479255; test on tests/hacc_occ.
- occE build 479258 OK (patch sha256 bd3c33ee…; binaries 03eabe9b… scalar, 6643177d… MPI). Test (tests/hacc_occ,
  He19 1 channel lMax 4 rMax 30): `occupation_energies` holds -0.901876, -0.143518 (1s2s), -0.122453 (1s2p), -0.060849
  (1s3s), -0.055524 (1s3p/3d) au relative to the He+ threshold; 1s2s excitation 20.64 eV in that small basis.
- Run 5 submitted (479260/479261, TRECX = the occE MPI build): inputs and pulse files checked present.

## 11 Oct: the "silent" upload failures explained

In zsh, `$H:s...` (a host variable followed by a colon and a path beginning with `s`) applies the `:s` parameter
modifier (a substitution) to $H, so rsync/scp got a local path, not host:path, and copied into oddly named local
directories (copies of my own scripts/inputs only, removed). The relay was not dropping files. Write `${H}:path`.

## 11 Oct ~10:00: runs 2-6 extracted (`extract_runs.py`, `summary_runs2to6.md`)

- Fixes on the way: (1) rows keyed by name only (one-electron and haCC share names; haCC overwrote) -> keyed by
  (model, name); (2) run 4 pair spectrum: numpy's positive-frequency part is exp(+i omega t) = E-, and tRecX's Laser
  output is sampled every ~9 au (3 points per cycle) -> the field from the pulse file (dt 0.01) and E+ = conj(...).
- Run 4: TDSE / PT = 0.999 (1e10), 0.99 (1e11) for both m and both durations: the TDSE confirms |M_m| to 0.1 %.
- Run 5: state 5 (-0.05555) is populated by two photons and state 6 (-0.05474) is not: 5 = 1s3d, 6 = 1s3p. State 3
  (-0.07776) is never populated (< 1e-23) and matches no He singlet: unidentified (an artefact of the haCC basis?).
  ±10 fs^2 chirps give 1s2s populations differing by ~16 % at 1e12 (lowest order predicts equal on exact resonance;
  2 x 10.306 = 20.613 eV vs haCC 1s2s 20.636 eV, 23 meV detuned).
- Run 6: the one-electron drift is window-dependent (the chirped windows are 3.6 x longer than the 1 fs no-field run),
  so the projection is biased for the chirped pulses; use the spectrum integral (agrees to 0.2 % at 1e13 for TL pulses).

## 11 Oct ~17:40: optmat's addendum 3 (run 5 checks)

optmat takes runs 1-4 and 6 as they are; run 5 (haCC 1s2s) is 0.46-0.52 of lowest order at 10 fs (1.0 at 3 fs,
1.0-1.2 at 1 fs), as if haCC's two-photon resonance were ~0.12 eV from its reported 1s2s; and the ±10 fs^2 chirps differ
by 16 %. Checks asked, cheapest first:
1. **P(1s2s)(t)** (run5/hacc/*/pop_vs_t.csv, printed every 10 au): **not constant after the pulse.** Spread after the pulse
   (max - min)/mean: 10 fs 1e12 5.3 %, 10 fs 1e13 0.6 %, 3 fs 13 %, 3 fs ±10 fs^2 24 %, 1 fs 38 % (absolute 1e-7 to 9e-7).
   The ±10 fs^2 end values (8.98e-7, 7.75e-7) both lie in a common post-pulse range 7.2-9.1e-7: the 16 % "asymmetry" can
   be this fluctuation sampled at different phases. Note the ground-1s2s beat period is 8.3 au < the 10 au print step
   (aliasing). At 10 fs P(1s2s) peaks during the pulse (1.668e-5) and falls back to 1.56e-5: detuned-like. The
   fluctuation cannot explain the factor 0.5 at 10 fs.
2-4. `make_runs.py run5b` -> run5b/ (manifest_run5b.csv), all haCC as run 5 (occE MPI build), print every 1 au:
   carrier scan at 10 fs 1e12, 10.306 eV + {-75 ... +75 meV} (7 runs); the 10 fs reference with accuracy 1e-9 (haCC's
   default adaptive RK4 accuracy is 1e-7); the 3 fs ±10 fs^2 pair at 1e12 with accuracy 1e-9 and +3000 au after the pulse;
   a no-field run over the 1 fs window (the OccDual floor). Jobs 479500/479501 (16 ranks, in sequence).

## 11 Oct ~22:15: addendum 3 checks 2-4 done and extracted (run5b)

- Jobs 479500 (tp_r5b_a, 17:34-20:43) and 479501 (tp_r5b_b, 17:34-20:18) COMPLETED; all 11 runs exit 0, 9.5-49 min each.
- `extract_run5b.py`, run on the cluster in tdse_two_photon/ (`python extract_run5b.py`, the round-2 runs' Python environment),
  11 Oct 22:12 BST:
  - script sha256 **b8b94831ae4ebf6a8f9f0c14240123fd59e713d165a219d2e3865a833b717df2**;
  - `run5b/results.csv` **3ccb2ec49c5944978bf0d8406d34ceed56bf1d7897c892886666887815e3e84d**;
  - `run5b/scan_fit.csv` **927b0532f991467a3b520b4780479d6ea36ae5447517bef8baeb66742bb9d9e9**;
  - `run5b/pop_<name>.csv` (populations every 1 au).

  All copied to the workstation; the checksums match.
- First pass: the beat period from rfft assumed uniform steps (expec has a few irregular rows: dt 0.48, 2) and gave
  ~7.2 au. Replaced by a non-uniform DFT, which gives 0.7621 au (8.245 au) in every run5b field run.
- **Results** (details in REPORT.md, "Addendum 3"):
  - The scan fit gives haCC's two-photon 1s2s resonance at **20.7386 eV**, +102.6 meV above the Occupation set's
    20.6360 eV. The post-pulse beat is 0.7621 au = 20.738 eV, the same spacing.
  - At that resonance, every 1e12 run 5/5b run is lowest-order with K = 3.078 au (4K = 12.31) to ≤ 1 % (1 fs 3 %).
  - The P(1s2s) beat is an interference with the OccDual floor |eps|^2 = 2.862e-9 (no-field run): predicted rms equals
    the measured rms to ≤ 1 %.
  - ±10 fs^2 post-pulse means agree to 0.07 % (acc 1e-9) and 0.4 % (acc 1e-7). The 16 % in run 5 was end-value
    sampling.
  - Accuracy 1e-9 vs 1e-7: 0.08 %.
- **No-field run:** with no laser, expec has only the first and last rows (tRecX takes one step). The floor is equal at
  both ends.
- **Not explained:** why the Occupation operator's spectrum differs from the propagated one. Its ground state is
  -0.901876 au, against the initial state's -0.901867. The 1s2p (used in run 2's haCC pole distances) is unmeasured.
  Offered: one weak one-photon run near 21.2 eV at a 1 au print step.
- Side observations:
  - the ± chirp norm-loss at 1e12 differs (1.52e-9 vs 1.91e-9); not examined;
  - the run5b titles contain a colon ("run 5 check 3:"); tRecX accepted them, but avoid colons in titles in future.

## 11 Oct ~23:00: optmat's addendum 4 (haCC's propagated 1s2p), run5c set up and started

- **Ask:** the propagated 1s2p–ground spacing, from the beat in P(1s2p) after the pulse.
  - Pulse: one photon at 21.21 eV (58.46 nm), 10 fs, 1e9 W/cm², haCC as run 5b, accuracy 1e-9, +3000 au after the
    pulse, print every 1 au.
  - Also requested: P(1s2p), and the floor if it differs.
  - Optional: a carrier scan at ±50 and ±100 meV fitted to K|S1|², and state 3's beat.
- **Problem with the beat method.** The 1s2p projector is ¹P, so by parity it cannot see the ¹S ground state. The run5b
  no-field P(1s2p) floor is **1.3e-19**, against 2.86e-9 for 1s2s. A ground-state leak beat would have a relative
  amplitude of about 1e-8: not measurable.
- **Changes made:**
  - The scan is run as a primary method, not an optional one.
  - `Dipole: length, velocity` is added: after the pulse, ⟨z⟩(t) rings at the 1s2p–ground spacing with amplitude
    2|c0 c1 z01|, of order 1e-2 au. This works with haCC; see tutorial/201He400.
- **Set-up:**
  - `make_runs.py run5c` (`run5c()`; `hacc_b` gained `dipole=` and `brief=`) → `run5c/`, `manifest_run5c.csv`;
  - five runs, `he_r5c_1s2p_d{0,±50,±100}meV`; d0 is the main run.
  - Estimated P(1s2p) on resonance: 1.5e-3 if |z| = 0.7347, 5e-4 if |z| = 0.42.
- **First submission failed at input parsing** (479624/479625, both cancelled or ended within 1 min):
  - "Dipole: length, velocity" needs a blank line after it, or tRecX reads the next line (Spectrum) as its values.
  - The failed starts were moved, not deleted, to `run5c/failed_parse_20261011/`.
- **Resubmitted as 479626** (tp_r5c_a: d0, d−100, d+100) **and 479627** (tp_r5c_b: d−50, d+50), occE MPI build,
  16 ranks each, 23:00.
- With Dipole on, expec gets a row at every adaptive step (~0.035 au), not every 1 au. Speed is ~38 au per 46 s, so
  ~1.8 h per run.
- Extraction: `extract_run5c.py`. Beat fit on DipLenZ and DipVelZ after t > 2T. It recovers a synthetic test frequency
  to 1e-8 au.

## 12 Oct ~10:15: addendum 4 results (run5c)

- Jobs 479626 and 479627 finished overnight; all 5 runs exit 0, 1.9-2.0 h each.
- **Corrected brief from optmat:** |z| reference 0.4208 au, and K = |d_z|². The runs are unaffected: the extraction
  already used K = |z01|² with E = E+ + E-.
- **First extraction too slow.** It ran on all ~152k rows per run and finished only 1 run in 50 min on the login node.
  Stopped (my PID) and changed to resample onto a 1 au grid, with a coarse-then-fine scan grid. Run time 3m47s.
- Outputs:
  - `extract_run5c.py`: sha256 **3f24eabad044dbe47d6cc85229628023763176b8ae489f1962170a429b1df393**;
  - `run5c/results.csv`: **c604f7681e1e9cc47289b3fb5a14ac6a71335e39129998786983e5fd92d8d586**;
  - `run5c/scan_fit.csv`: **ea64b98bf71eec5c9fe73896d659a3f3f0c0f20cfda916be7cc5f009d221e187**.

  Run on the cluster at 10:12 BST, 12 Oct; copied to the workstation with checksums identical.
- **Results:**
  - ⟨z⟩ beat Δ = 0.7794133 au = **21.208918 eV** (10 fits, spread 0.9 µeV); scan Ω = 21.20892 eV.
  - The Occupation 1s2p equals the propagated 1s2p to 1e-8 au absolute. The 0.25 meV spacing difference is the
    Occupation ground state's offset.
  - K = 0.169939, d_z = 0.41224 au (−2.0 % vs NIST's 0.4208).
  - P(1s2p) = 4.688e-4 at 21.21 eV, first order to 5e-5; state 3 not populated.
  - Velocity/length amplitude 0.8257 vs Δ 0.7794 (+5.9 %): the gauges are not equivalent in this basis.
