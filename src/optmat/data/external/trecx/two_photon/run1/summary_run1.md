# Run 1 (helium, degenerate two-photon ionisation): summary

20 fs Gaussian (FWHM of the intensity) on A, CEP 0; I0 = (1/2) eps0 c E0^2; sigma2 = P / int Phi^2 dt with
Phi = I(t)/(hbar omega), I(t) = I0 exp(-4 ln2 t^2/T^2). Full rows: results_sae.csv, results_hacc.csv.

**sigma2 (1e-52 cm^4 s), at 1e12 W/cm2:**

| E (eV) | one electron: projection, BOX 40 (drift-subtracted) | one electron: spectrum peak, BOX 40 / 80 / 120 | haCC He19: spectrum peak (rMax 40) | haCC: norm-loss | haCC / one electron |
|---|---|---|---|---|---|
| 13 | 7.456 | 6.904 / 7.448 / 7.595 | 11.890 | 11.887 | 1.595 |
| 14 | 7.328 | 7.988 / 7.400 / 7.305 | 12.804 | 12.805 | 1.747 |
| 16 | 6.524 | 6.891 / 6.590 / 6.599 | 10.723 | 10.724 | 1.644 |
| 18 | 6.228 | 6.307 / 6.264 / 6.189 | 9.063 | 9.068 | 1.455 |
| 20 | 9.346 | 9.540 / 9.285 / 9.404 | 10.175 | 10.182 | 1.089 |

- One electron: the spectrum peak at BOX 40 is off by up to ±8 %; at BOX 80 and 120 it agrees with the
  (box-independent) drift-subtracted projection to <= 2 %. Use the projection (or BOX >= 80 spectra).
- Drift (no-field norm-loss over 5961 au): 1.432e-07 (base), 1.891e-08 (dt x 2/3, cutEnergy 150). At 1e11 the
  one-electron signal (~2e-8) is below it: 1e11 from the spectrum. P ∝ I^2 to <= 1.5 % (1e11–1e13).
- One-electron convergence at 16 eV, 1e12 (projection): Eta 15 0.000 %, theta x 2/3 -0.3 %, dt x 2/3 / cutEnergy 150 +1.6 %
  (after their own drift), BOX 80 / 120 0.0 %.
- haCC: the norm-loss and the spectrum peak agree to <= 0.1 % (1e12, 1e13; 1.7 % at 1e11); lMax 15 +0.003 %, rMax 60 -0.23 %.
  P ∝ I^2 to <= 1.2 %. ~45 min per run on 16 ranks.
- Ip: one electron 24.5937 eV (Tong-Lin, its own spectrum); haCC 24.5410 eV (-E_init). No NIST values used.

