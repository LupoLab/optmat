# pyroraptor helium yields: two electrons, haCC and one electron (J. C. Travers, 2026)

Helium in strong laser pulses, computed by J. C. Travers with pyroraptor, a
time-dependent Schrödinger equation solver (one active electron, both electrons of
helium, and haCC, the hybrid anti-symmetrised coupled channels method of Majety, Zielinski
and Scrinzi, New J. Phys. 17, 063002 (2015)). Unpublished. The source record is
`travers_2026_pyroraptor`; its data files are written from the four result tables here
by `tools/extract/travers_2026_pyroraptor.py` (TC-1092). Tier 3: nothing is fitted to
them, and they rank below published calculations and measurements.

## The runs

The run records (JSON) are in pyroraptor's repository, main at
`f9b35142812833113fcbd8af63554b40f1d95f4d` (16 October 2026), with the reports that
discuss them. The tables here were written from those records by
`records_to_csv.py.txt` (this directory; run with Python 3.14 on 16 October 2026). It
copies every number with full precision and leaves out the records' host and build
strings. `records.sha256` gives the sha256 of each record read, by family and file name.

| Table | Records | What |
|---|---|---|
| `he_two_electron_390nm.csv` | `two_electron/k2b-*.json` (8) | both electrons active; the bound-state analysis of the state at 744, 1244 and 1744 au; boxes R0 = 40, 59, 78, 97 bohr |
| `he_hacc_390nm.csv` | `hacc_390/*.json` (26) | haCC with the He+ 1s channel (`ion_principal` 1, "haCC(1)") or the n <= 2 channels (2, "haCC(5)"), the base, the convergence variations and the neutral-square ladder, 1000 au after the pulse |
| `he_sae_390nm.csv` | `sae_390/*.json` (4) | one electron in Tong and Lin's helium potential, four discretisations, 1000 au after the pulse |
| `he_hacc_800nm.csv` | `hacc_800/h1/*.json`, `h5/*.json` (54) | haCC(1) and haCC(5) in Sato et al.'s six-cycle sin^2 pulse at 800 nm, the base, the variations and the final ones, at the end of the pulse |

The 390 nm pulse is A2's 6 fs Gaussian on the vector potential at 390 nm, as in the
tRecX runs of `travers_2026_trecx` (their part K): A(t) = A0 exp(-2 ln 2 t^2/T^2)
sin(omega t), omega = 0.1168291 au, T = 248.048 au (6 fs, the full width at half maximum
of A^2), max |A| = 0.8734300 au (3.684031e14 W/cm2) and 1.0854638 au (5.68981e14 W/cm2),
then 1000 au without field. The haCC records name the pulse by case (`390`, `390hi`);
the converter checks their max |A| against the intensity to 1e-4. The one-electron
records use the same max |A| (their driver's `K2`, labelled 3.684e14 and 5.690e14). The
800 nm pulse is Sato et al.'s (2016) Eq. 27, E(t) = E0 sin(omega t) sin^2(pi t/(6T)),
at 4, 8 and 14e14 W/cm2 (I = E0^2 in atomic units).

## The measures

- `above_n4` (`ionisation_plus_n5` in the one-electron and 800 nm records): N0 -
  P(ground) - the populations of the bound states with n <= 4, n = round(1/sqrt(-2E))
  with E the state's energy above the He+(1s) threshold. That is the ionisation plus the
  excitation to n >= 5. `above_n5` and `above_n6` likewise cut at n <= 5 and n <= 6.
  The two-electron, haCC and one-electron analyses use the same rule (the runs'
  analysis scripts in pyroraptor's repository).
- `beyond_rc`: 1 - P(both electrons inside 20 bohr) at the end of the pulse, Sato et
  al.'s measure.
- `norm_loss_end`, `ionisation`: as named. Neither is box-converged for two electrons:
  each larger box keeps another Rydberg shell bound.

**The reports compare different cuts.** pyroraptor's reports call haCC's `above_n5` (for
example 7.113e-4 at 3.684e14 W/cm2) "ionisation plus n >= 5" and compare it with the
two-electron and one-electron `above_n4`. By the rule above, `above_n5` is ionisation
plus n >= 6. optmat compares like with like (TC-1092).

## Checked (2026-10-16)

- Every row of every table against its record (`above_n4`, `above_n5`, the norm-loss and
  `beyond_rc`), to every digit: 112 rows.
- The two-electron `above_n4` at 1744 au: 8.7008e-4, 8.4650e-4, 8.4670e-4 and 8.4641e-4
  at R0 = 40, 59, 78 and 97 (3.684e14), and 7.2562e-3, 7.2477e-3, 7.2431e-3 and
  7.2440e-3 (5.69e14). It is box-converged from R0 = 59, as the report states (8.464e-4
  and 7.243e-3 at 78 bohr).
- haCC's `above_n4` and `above_n5` at R0 = 40, 60 and 80 bohr (one ion channel, square
  3) agree to 0.01 %: unlike the two-electron calculation, haCC's bound states are
  converged at 40 bohr.
- haCC(1) and haCC(5) in Sato et al.'s pulse: 2.961e-4, 1.5445e-2 and 0.1536; 3.023e-4,
  1.5825e-2 and 0.1578, the report's values.

## Files

```
fe05ef08b2a25ca9f33b10298af0b1a41c6ef1a3f6320d59b766e748c9fb2160  he_two_electron_390nm.csv
60b3c0df12888762ffe25c8c37827e4d2739fffc414bbb55811388c9d9878259  he_hacc_390nm.csv
f0935bcf97904a4033d9ab3c7b2963a40fdd6d6652c2b53abaa1e2eb813fae7c  he_sae_390nm.csv
8a97719b6d89c526d2234701d54cbf1f970937f680f2143b13200a2fcc29cd82  he_hacc_800nm.csv
3ba3f210ce17f2ba63bf0eb5af1e36ede5322eee253caf2299e78fcdefacb93d  records_to_csv.py.txt
9effd372ac0035b999f0f5b99f56b09bae1f40d0a196ba989b5fb85d5f9302d9  records.sha256
```

## Terms

These are our own calculations, run with pyroraptor (not included), and are
distributed under optmat's MIT licence.
