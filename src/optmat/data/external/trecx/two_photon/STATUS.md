# Two-photon TDSE campaign (two-photon ionisation and excitation of He and Ar): status

A campaign separate from round 2 (`../ROUND2_STATUS.md`).
- Brief: optmat's specification of the two-photon TDSE campaign (unpublished), received 9 Oct 2026.
- John: "read ... and start them on the cluster"; track separately from round 2.
- Local: `testruns/two_photon/`. Cluster: `$RUN_ROOT/`. Slurm job names `tp_*`.
- Working log: `NOTES.md`.

## Jobs

| job | what | state |
|---|---|---|
| 478286, 478287/88, 478302/03 | run 1, first versions | superseded (spectral range; script bug), kept as `run1_v1_spec_kmax0p02/` |
| 478310, 478367 `tp_run1_sae*` | run 1 one-electron + box/drift checks (44 inputs) | **done** 9 Oct |
| 478347/478348 `tp_run1_hacc_a/b` | run 1 haCC (17 runs) | **done** 10 Oct |
| 479251 `tp_r2346_sae` | runs 2, 3, 4, 6 one-electron (110 inputs, 32 cores, 4-day limit) | submitted 10 Oct 23:20 |
| 479252/479253 `tp_r2_hacc_a/b` | run 2 haCC (16 runs, 8 + 8 in sequence, 16 ranks) | submitted 10 Oct 23:20 |
| 479250, 479255, 479258 `build_occE` | tRecX a987d7c + Occupation energies written to `occupation_energies` | done (479258; the first two: message silent, then a quoting error) |
| 479260/479261 `tp_r5_hacc_a/b` | run 5 haCC (15 runs, Occupation 0:20, occE MPI build) | **done** |

**Runs 2–6 done** (11 Oct, all exit 0; 479251–53 done too): `summary_runs2to6.md`, `run{2..6}/results.csv`.

**Report for optmat: `REPORT.md`** (11 Oct). Open on request: the optional haCC argon mMax 1 point (run 3), a length-gauge rerun.

**Run 1 result:** `run1/summary_run1.md`. sigma2 (1e-52 cm^4 s) at 13/14/16/18/20 eV: one electron 7.46/7.33/6.52/6.23/9.35,
haCC 11.89/12.80/10.72/9.06/10.18 (haCC/one electron 1.09–1.75). Converged to ~2 % (one electron), <= 0.3 % (haCC).
haCC ~45 min per 20 fs run on 16 ranks: affordable for run 2.

## Plan (decisions from optmat's reply, relayed by John 9 Oct)

1. **Run 1** (helium degenerate): as submitted, including 1e13; at 1e11 use the iSurff spectrum integral.
2. **Run 2** (helium two-colour): ONE run per point with both colours (no P1/P2 runs); send the integral of the peak
   at ħ(ω1+ω2) − Ip, and a ν = 0 reference with the same pulse (here: one colour at Ω/2, same envelope, the
   per-colour intensity; stated). haCC if run 1 shows it is affordable, one electron as well if cheap. 10 fs allowed
   if it halves the cost, with the reference at the same pulse.
   - Peak spacing is ν (optmat agreed, 9 Oct). Keep 20 fs; if haCC forces 10 fs, ν = 0.5 eV at 20 fs anyway.
   - ν = 0 reference: one colour at Ω/2, same envelope, per-colour intensity (agreed; optmat applies the factor 4).
3. **Run 3** (argon): one-electron Tong–Lin only, m = 0 and m = 1; the atom = 2 P(m=0) + 4 P(m=1) (stated); report
   the model's 4s and 3d energies. No haCC mMax 0. Optional, lowest priority: one haCC mMax 1 degenerate point,
   9.5 eV, 6 fs (~3 days).
4. **Run 4** (Ar 3p⁵4p): one-electron. Per initial 3p m (0, 1), the complex projections ⟨4p m|ψ(T)⟩ with the phase
   convention, the exact field, and the model's 4p and 4s energies. tRecX prints only populations, so ψ(T) and the
   4p eigenvectors are written and projected offline; cross-check the sign by lowest-order perturbation theory in
   the same potential. A factor-of-2 estimate suffices; a lower-tier number.
5. **Run 5** (He 1s2s dynamics): haCC only (the one-electron model has no singlet 1s2s). Needs a bound projection in
   haCC (Occupation with 1 channel): to test during run 1.
6. **Run 6** (broadband identity): one-electron, against its own run 1 cross section.

## Addendum 3 (11 Oct): run 5 checks

- Check 1 (P(1s2s)(t)) done from the existing runs: not constant after the pulse (5–38 %). Checks 2–4: run5b/, jobs 479500/479501 (11 haCC runs).
- **11 Oct 22:15: all four checks done; reported in REPORT.md "Addendum 3".**
  - haCC's propagated 1s2s resonance is at 20.739 eV, +102.6 meV above its Occupation set's value.
  - With that resonance, run 5 at 1e12 is lowest order to ≤ 1 % (4K = 12.31 au²).
  - The P(t) beat is an interference with the OccDual floor of 2.9e-9.
  - The ± chirp difference vanishes in the post-pulse mean (0.07 %).
  - Open (offered): the propagated haCC 1s2p.

## Addendum 4 (11 Oct): haCC's propagated 1s2p

- **Status:** run5c, 5 haCC runs, jobs 479626 and 479627 (tp_r5c_a, tp_r5c_b), started 23:00 11 Oct.
- **Method:**
  - the ⟨z⟩(t) beat, since the P(1s2p) projector has no ground-state floor by parity (1.3e-19);
  - the scan fit to K|S1|².
- **12 Oct 10:15: done; reported in REPORT.md "Addendum 4".**
  - Propagated 1s2p: 21.20892 eV (beat 21.208918 eV ±1 µeV; scan 21.20892 eV on its 0.05 meV grid). This equals the Occupation 1s2p, so run 2 stands.
  - d_z = 0.41224 au.
