# Changelog

All notable changes are recorded here, in the style of
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The review register
(`REGISTER.md`) records the changes to data and models in detail. Versions before
0.8.0 were internal and were not published.

## Unreleased

## 0.9.0 (2026-10-24)

Models read only the data they need, so start-up is seconds to tens of seconds faster
on a cold disk, and helium and xenon gain validation data. No recommended value or model
changed. `Source.datasets` is no longer a constructor field (see *Changed*), hence the
minor version.

### Added

- Helium on matched measures (TC-1092): a new source of tier-3 calculations,
  `travers_2026_pyroraptor` (both electrons active, haCC with one and five ion
  channels, and one electron, at 390 nm; haCC in Sato et al.'s 800 nm pulse), and
  `gas.ionisation_validation.helium_matched_ratios`.
- Xenon in Hartmann et al.'s RMT field (TC-1093, TC-1094): one 5p electron in the same
  field and with the same measure (`travers_2026_trecx`, part N), the paper's measured
  neutral abundance (`hartmann_2022`), and `xenon_rmt_ratios` and `xenon_survival`.
- The curated data are now 762 source records: 8481 datasets and 197,804 values, 259
  errata and 780 published formulas.
- `gas.ionisation_rates.ionisation_rate_coverage` and
  `gas.molecular_ionisation_rates.molecular_rate_coverage`: the gases each rate model
  covers, read from the rates file without building a model, so a caller can tell an
  atom from a molecule.
- `gas.models.fitted.fitted_model(name)` and `gas.models.formulas.published_model(name)`:
  one gas model, reading only its own record.
- `Source.dataset(id)` and `load_source(path, defer=True)`: a source's datasets are
  built, and their data files read, when first asked for.
- A catalogue index, `data/catalogue_index.json`, written by
  `tools/catalogue_index.py`: which record holds the datasets and formulas of each
  material and quantity. It is checked against the records on first use; an out-of-date
  index raises the new `CatalogueIndexWarning` and the records are read instead.

### Changed

- The models read only the data they need (no value changed). A model that needed a
  few datasets read the whole catalogue (8,300 files), and a gas's model built every
  gas's: `molecular_rate_model("N2")` read 8,671 files and now reads 35,
  `ionisation_rate_model("He")` 339 and now 34, gas `default_model("He")` 195 and now
  28, `air_model()` 799 and now 19, and the SiO2 models 3,964 and now 92. Datasets are
  looked up by key, each gas's model, rate and parameters are built on demand, searches
  by material, quantity or formula use the index, and a source's data files are read
  only for the datasets used. Cold starts, where every file opened is read from disk
  or scanned, are seconds to tens of seconds shorter.
- `Source.datasets` is a cached property built from the new `Source.entries` (the
  record's dataset tables), no longer a constructor field.

### Fixed

- The manual's helium findings at 390 nm (TC-1092 corrects TC-1050): haCC's helium is
  0.96 to 0.99 of the exact atom, not 0.76 to 0.81, and the atom 1.6 times one
  electron, not 1.7 to 1.8. The earlier ratios compared a norm-loss that counts Rydberg
  excitation with one that does not. Parker et al.'s rate is confirmed as before.

## 0.8.1 (2026-10-16)

Fixes for a downstream propagation code. No value at a positive wavelength, model or
stored fit changed.

### Added

- optmat's DOI, 10.5281/zenodo.23422879 (Zenodo; every version), in `CITATION.cff`
  and the README.
- A `py.typed` marker (PEP 561), so that type checkers of downstream code read optmat's
  annotations.
- `optmat.solid.models.formulas.published_model(name)`: one published formula, reading
  only its source's record.

### Changed

- A model refuses a wavelength at or below zero, and a solid model refuses inf, with
  `ArgumentError`; before, a negative wavelength returned a finite, meaningless value
  (TC-1091). A gas and air keep inf as the static limit; NaN still gives NaN.
- Faster lookups for a propagation code. `solid.models.registry.model` (and so
  `default_model`) builds a published formula from its own record, not from every
  record with formulas (N-BK7, BaF2, MgF2: 0.3 s from about 2 s in a fresh process),
  and `gas.air.air_model` reads only the records with mole fractions, not the whole
  catalogue (about 1.4 s from 4-5 s). No value changes.

## 0.8.0 (2026-10-16)

The first public release, at github.com/LupoLab/optmat. No value, model or stored fit
changed.

### Changed

- Wider support: Python 3.12 and later (was 3.14) and NumPy 2.3 (was 2.5); SciPy stays
  at 1.18. These are the oldest versions that reproduce the stored results: Python 3.12
  for PEP 695 syntax; NumPy 2.3 and SciPy 1.18 because the stored fits and the models'
  identity digests, which hash computed arrays, reproduce bit for bit only from these
  releases. CI tests them.
- The `gas` extra no longer installs pyfluids, which optmat never imported; it calls
  CoolProp directly.
- Coverage is measured with `pytest-cov` (as in solprop), replacing
  `tools/measure_coverage.py`; the gate still requires 100 % statement and branch
  coverage, including ordinary-Python subprocesses, and writes an XML report for CI.
- The register's branch-integration recipe moved from `AGENTS.md` to *Maintaining the
  review register*. Agents create their topic branch and worktree without asking and
  keep the main checkout clean.
- Private development references (task briefs, machine paths) removed from the tRecX
  archive, the register and the source records; wording only.

### Added

- CI (`.github/workflows/ci.yml`): the full gate on Linux and macOS, and the tests on
  Python 3.12 with the lowest supported dependencies; a release workflow publishing to
  PyPI by trusted publishing.
- Every file states its copyright and licence (REUSE): `REUSE.toml`, `LICENSES/`, SPDX
  headers in the code, and `reuse lint` in the gate. Code, source records and fits are
  MIT, prose and figures CC-BY-4.0. The published data the package ships (the values
  transcribed into `data/values/` and the files in `data/external/`) remain their
  publishers', included under a third-party data notice
  (`LICENSES/LicenseRef-optmat-third-party.txt`) so that the package is self-contained;
  each external directory's `PROVENANCE.md` states its publisher's terms.
- `CITATION.cff` and `.zenodo.json`.
- Windows: every text file is read and written as UTF-8, whatever the platform's
  default (Windows decoded the curated records with its code page, which failed on
  their author names and unit signs); `tests/test_encoding.py` keeps it so, and CI
  runs the tests on Windows. A git checkout keeps LF line endings on every platform
  (`.gitattributes`), so the stored configurations' hashes match there too; the
  archives under `data/external/` are checked out byte for byte.

### Fixed

- Citation metadata: Fernley, Taylor and Seaton's page (6457, as reproduced in
  Berkowitz 1997), the co-authors of van der Hart et al. 2007, the pages of van Duijnen
  and Swart 1998 and the first page of Manakov et al. 1986. A source note no longer
  quotes a reference that does not exist.
- The manual: Blow and Wood's α is their Kerr weight; Ammosov, Delone and Krainov define
  n* in their eq. 2; and the static rate is their eq. (1) without its cycle-averaging
  factor.

## 0.7.0 (2026-10-15)

Crystalline silicon's complex refractive index from 31 nm to 143 µm at 20 °C: `fit/Si`,
Franta et al.'s universal dispersion model refitted to Frey et al.'s absolute indices,
Schinke et al.'s absorption edge and Franta et al.'s 20 °C tables, now the recommended
model of `Si` (it had none), with twelve new source records and the manual page
*Crystalline silicon*. Also argon's rounds for OEP 8, RMT's xenon ground-state survival,
a five-section manual with a generated source catalogue, and refit tests that compare
the stored solid, gas and compact fits on their predictions and objective rather than
their parameters. No other recommended value changed; no stored fit changed.

### Added

- A five-section manual: Introduction, User guide, The data, Developers guide, and
  API. The introduction includes installation guidance, a functionality overview,
  the quick start, and an explanation of curation; the user guide groups the
  scientific chapters in reading order, with previous/next navigation.
- A reproducible source catalogue generated during documentation builds from every
  source record, including citations, datasets, formulas, notes and corrections.
- Separate chapters for nonlinear conventions, modelling, support and extrapolation;
  Raman validation; atomic and molecular ionisation validation; and two-photon
  source coverage and bound excitation. Existing page URLs and moved section
  fragments lead to the new locations.
- Dedicated developer instructions for extraction, fitting, documentation and register
  maintenance. Historical solprop comparisons now live in regression fixtures; the
  public validation retains factual comparisons with an identified Luna.jl snapshot.
  The manual, gallery, README and API descriptions use consistent terminology and
  distinguish current defaults from historical validation results.

- The vocabulary quantity `ground_state_survival`: the probability that a pulse leaves
  the atom in its ground state, read at a stated time (TC-1058).
- `hartmann_2022`: RMT's xenon ground-state survival at 759 nm, 1.1 to 2.5e14 W/cm2, from
  the authors' archive, with the pulse read from the fields RMT wrote (a 1/e half-width
  of 186 au: the paper's "FWHM 4.5 fs" is a 5.30 fs intensity FWHM). `hutcheson_2025`
  is recorded without values: its archived populations were computed at about 670 nm,
  not the stated 760 nm (TC-1059).
- Argon's round for OEP 8 (TC-1060 to TC-1068): eight sources.
  - The Ar 4p levels' lifetimes and quenching by 22 gases after two-photon excitation at
    184-189 nm (Sadeghi et al. 2001), confirming Wiese et al.'s radiative widths.
  - Measured and calculated Ar 4p-4s transition probabilities (Asghar et al. 2020;
    Savukov 2018).
  - Calculated photoionisation of every Ar, Kr and Xe np level near threshold (Petrov et
    al. 2011).
  - Kr and Xe lifetimes and transition probabilities (Aymar and Coulombe 1978).
  - Three Stark-effect papers, as notes.

- Argon's records for OEP 8 (TC-1069 to TC-1074):
  - Chang and Setser's Ar 4p lifetimes, quenching and transition probabilities, and
    Nguyen and Sadeghi's state-to-state transfer (both 1978). Together they show that
    Sadeghi et al.'s self-quenching of 4p[1/2]0 and 4p[5/2]2 is low by about 2.
  - Baig et al.'s measured 4p' photoionisation at threshold (2008), 2.8-3.1 times
    Petrov et al.'s calculation.
  - One-photon 4s-4p pressure broadening and shifts in Ar (Sadeghi et al. 2022) and
    Ar-He (Eshel et al. 2016).

  The curated data are now 749 source records: 8357 datasets and 166,428 values, 255
  errata and 771 published formulas.
- **Crystalline silicon's complex index** (TC-1075 to TC-1088): `fit/Si`, the
  recommended model of `Si` from 31 nm to 143 µm at 20 °C (manual page *Crystalline silicon*).
  - Franta et al.'s universal dispersion model of float-zone silicon (the newAD2
    parameter file, now `franta_2017/udm_c_si`, 0.016 high in n in the near infrared
    and 2.6 times high in k at 633 nm), refitted to Frey, Leviton & Madison's absolute
    indices at 1.2-5.5 µm (taken from 295 K with their dn/dT; within 1.4e-4), to
    Schinke et al.'s absorption coefficient at 250-1200 nm, and to Franta et al.'s
    temperature-dependent model at 20 °C (their refractiveindex.info tables) in the
    ultraviolet and the infrared lattice bands.
  - Twelve source records: Frey et al. 2006 (Si and Ge, 30-295 K, with their dn/dT and
    temperature-dependent Sellmeier formulas), Li 1980, Green 2008, Schinke et al. 2015,
    Aspnes & Studna 1983 (Si, Ge, GaAs, GaSb, InP), Jellison 1992, Herzinger et al. 1998,
    Edwards's chapter in Palik's handbook (1985), Philipp & Ehrenreich 1963, and three
    Franta et al. papers (2014, 2018, 2019); franta_2017 gains the newAD2 model and
    tables and the 20 °C, 25 °C and 300 K tables.
  - `optmat.solid.fitting` fits an `absorption_coefficient` as k = alpha lambda /
    (4 pi).
  - The model's k beyond 1.15 µm is an upper bound, not the band-to-band absorption;
    it is for intrinsic float-zone silicon at 20 °C only. The curated data are now 761
    source records: 8467 datasets and 197,754 values, 259 errata and 780 published
    formulas.

### Fixed

- The solid fits' refit tests (`fit/SiO2`, `fit/Si`) compare parameters in units of
  their stored uncertainties, predictions over the model's range and the objective,
  instead of each parameter to 1e-6 relative. Silica's stored fit had never matched a
  refit at that precision: the optimiser stops at a path-dependent point of a flat
  minimum, 0.003 sigma away, at the same chi-squared to 2e-9 (TC-1089). No stored
  value changed.
- The gas dispersion fits' and compact models' refit tests compare predictions inside
  each fit's intervals (to 5 % of the model's band and 1e-4 relative), the full
  objective (to 0.01 in chi-squared) and the compact models' residuals (to 0.01
  sigma), instead of the polarizability to 1e-7 out to 65 nm and the poles to 1e-9.
  Every stored gas fit and compact model had stopped elsewhere on a nearly degenerate
  minimum than a refit of the committed code does, with the same data objective to
  round-off (TC-1090). No stored value changed.

- Numerical regression checks tolerate demonstrated arithmetic differences without
  changing model evaluation or stored fits. Gas Raman golden tests compare canonical
  identity records as well as exact self-consistent hashes; the existing golden
  numerical values are preserved. Nonlinear refit tests check predictions, parameter
  displacement in uncertainty units, dataset scales and the full objective including
  priors. Mutation tests retain sensitivity to substantive regressions (OEP 1, Q32).
- Review of the manual's restructure:
  - The Raman page no longer assigns CH4 and SF6 to the normal-coordinate strengths
    (they use Schrötter & Klöckner's Σ).
  - D2's Q-width ratios are stated with their density dependence.
  - The measured-yield comparison gives the current default's intensity scales and
    ratios, not the static-factor default's.
  - An open question of the ionisation accuracy page and the ratio direction of the
    static-rate table are restored. So are the capillary zero-dispersion sensitivity
    and the bound on the temperature dependence of the index.
  - The CH4 gallery caption's plotted range and the Zhang eq. (8) citation are
    corrected. Cross-references left behind by the page split are retargeted.
  - The moved solprop tables have their column definitions back.

## 0.6.0 (2026-10-12)

The bound two-photon levels of the noble gases and the literature behind OEP 8:
`optmat.gas.two_photon_bound`, the even levels below threshold from NIST with every
curated two-photon strength in one per-level convention and the excited levels'
photoionisation (OEP 6, Phase 2a; OEP 6 is Final). Seven vocabulary quantities for
radiative, collisional, light-shift and amplitude data. About 160 sources on two-photon
dynamics, transition probabilities, lifetimes, quenching, pressure effects, hydrogen's
exact amplitudes and the xenon n2. The maintainer's two-photon tRecX runs, archived,
confirm the helium model and OEP 8's weak-field form. Also the strong-field
dissociation data of OEP 7 (Draft) and helium with two active electrons at 390 nm. No
recommended value changed. The default ionisation-rate model's helium digest changed with
its source record travers_2026_trecx (TC-1050); its values are unchanged.

### Added

- **OEP 6 is Final** (OEP 6, *Two-photon ionisation and absorption of the noble
  gases*): the continuum model `optmat.gas.two_photon` and the bound-level records
  `optmat.gas.two_photon_bound`. Its Phase 2b and Phase 3 pass to OEP 8.

- `travers_2026_trecx`: helium with two active electrons at 390 nm, haCC's helium and
  one electron in the same pulse; `parker_pulse_yield`, `helium_two_electron_ratios`.
  The atom is 1.7 to 1.8 times one electron; Parker et al.'s rate over the pulse is
  within 1 to 5 % of it; haCC's helium is 0.76 to 0.81 of it. The manual's helium table
  labelled the static factor's column `default` since 0.3.0; corrected (TC-1050). The
  curated data are now 734 source records: 8162 datasets and 164,315 values, 255 errata
  and 771 published formulas.

- NIST's near-infrared lines of Ne, Ar, Kr and Xe, held as returned beside the
  resonance lines (TC-983): the transition probabilities among excited levels, for
  two-photon amplitudes, in vacuum wavelengths throughout (re-queried, TC-985).
- The literature of OEP 8 (TC-986 to TC-1010): 24 sources on two-photon dynamics,
  transition probabilities, lifetimes and the xenon n2 near 248 nm. Taken as data: three
  Kr level energies and the ionisation limit of the 86Kr isotope (Brandi et al. 2002),
  the Kr and Xe resonance oscillator strengths (Matthias et al. 1977), and reproduced
  He and Xe alpha and gamma (Manakov et al. 1986). The A values, lifetimes, quenching
  rates and excited-level polarisabilities await vocabulary.
- Five quantities for the bound levels' radiative and collisional data (TC-1011):
  `transition_probability`, `radiative_lifetime`, `level_polarizability`,
  `light_shift_coefficient` and `collisional_quenching_rate_coefficient`. Their data
  (TC-1012 to TC-1020):
  - the A values of Ar 4s-4p and 4s-5p (Wiese et al. 1989), Kr 5s-5p (Dzierżęga et al.
    2000) and Xe 6s-6p (Sabbagh and Sadeghi 1977);
  - the lifetimes of the Xe 6p, 6p' and 7p levels (Bruce et al. 1990), the Kr and Xe
    resonance levels (Matthias et al. 1977) and the Kr 5p levels;
  - Xe quenching by Xe and Ar;
  - the static polarisabilities of the Ar 4s and 4p levels (Werbowy 2026);
  - the Kr 5p[1/2]0 light shift (Brandi et al. 2002).

- The second round of OEP 8's literature (TC-1021 to TC-1041): 22 sources.
  - Two quantities (TC-1040): `pressure_line_coefficient` and
    `two_photon_transition_amplitude` (bound levels). Lifetimes, quenching and
    excited-level photoionisation now also apply to molecules.
  - Hydrogen's exact two-photon amplitudes against the photons' split (Quattropani et
    al. 1982, Drachman et al. 1990).
  - Xe 6p lifetimes and quenching (Bowering et al. 1986, Ku and Setser 1986), pressure
    shifts and broadening (Gornik et al. 1981, Hawkins et al. 1980).
  - J-averaged excited-state photoionisation of Ne, Ar, Kr and Xe (Duzy and Hyman
    1980).
  - The H2 and D2 dissociation energies and the Kr 5p[1/2]0 light shift (Balakrishnan
    et al. 1992, 1994).

- The third round of OEP 8's literature (TC-1042 to TC-1047): the Xe 6p', 7p and 6d
  lifetimes and quenching by Ne and Ar (Inoue et al. 1984), which put the 6p' lifetimes
  well above NIST's; the H2 E,F lifetime and quenching from their original (Kligler and
  Rhodes 1978); the hydrogen 1s-2s amplitude in both gauges (Bassani et al. 1977); He
  1s2p and 1s3p lifetimes and photoionisation (Johansson et al. 2003); two-photon Rabi
  oscillations in excited He (Fushitani et al. 2015).

- The maintainer's tRecX runs of two-photon ionisation and excitation (TC-1048,
  TC-1049), archived with their inputs: the two-photon ionisation cross sections
  of He (two-electron haCC and one-electron) and Ar (one-electron), and the one-electron
  Ar 3p-4p amplitudes, as tier-3 checks. haCC's helium agrees with optmat's model to
  0.90-1.16. The two-colour runs support OEP 8's pair-form rule for the He and Ar
  continua.
  Run 5b (TC-1051) puts haCC's propagated He 1s2s at 20.739 eV: with it, the 1s2s
  populations follow OEP 8's weak-field form to 1 %, and the fitted |M|^2 = 12.31 au^2
  equals Hattig et al.'s strength to 0.04 %.
  Run 5c (TC-1052) finds haCC's propagated 1s2p where its state list puts it (21.2089
  eV), so run 2's pole distances stand; haCC's 1s2p dipole, 2 % below NIST's, is held as
  a check.

  The curated data are now 734 source records: 8151 datasets and 164,296 values, 255
  errata and 771 published formulas.
- Strong-field dissociation data for O2, N2, N2O, NO2, O3 and N2O5 (TC-954 to
  TC-981): eight quantities (dissociation energies, ionic-state populations,
  fragmentation branching ratios, kinetic energy release, dissociation times,
  dissociative recombination coefficients, photodissociation quantum yields and
  transition dipoles); NO2, N2O5, NO, NO2+, NO+, N2++ and O2++ registered by reference
  only; 27 sources curated, among them the density-matrix tables of Yuen and Lin
  (2022, 2023) and Jia et al. (2024), Kleine et al.'s measured N2+ state populations
  (2022), Berkowitz and Eland's N2O+ branching ratios (1977), and Kossyi et al.'s
  recombination coefficients (1992). No dissociation model yet. The curated data are now
  682 source records: 7947 datasets and 163,579 values, 232 errata and 771 published
  formulas.

- **The bound two-photon levels of the noble gases** (OEP 6, Phase 2a; TC-952,
  TC-953): `optmat.gas.two_photon_bound.two_photon_bound_levels(gas)`. It gives every
  even J = 0, 1, 2 level below the ionisation energy from NIST, with its curated
  two-photon strengths in one convention (per level; `sigma0_from_delta`), measured and
  calculated, and its photoionisation cross sections. Helium and xenon are checked,
  neon, argon and krypton unresolved. There is no evaluator (Phase 2b waits on OEP 6
  Q9).
- Five further two-photon sources requested for krypton (TC-946 to TC-950): Ritchie
  1977 and L'Huillier, Jönsson and Wendin 1986 (helium two-photon ionisation), Kulander
  1987 (helium, time-dependent Hartree-Fock) and 1988 (xenon), and Crance and Aymar 1985
  (a note-only record). None holds krypton data: krypton still has no two-photon model
  (OEP 6, Q2). OEP 6 Phase 2 is restricted to curation for now (TC-945). The curated
  data are now 656 source records: 7810 datasets and 162,996 values, 232 errata and 771
  published formulas.

### Changed

- hummer_1998's threshold cross sections are placed at the calculation's own thresholds
  (TC-951, correcting TC-934). gangopadhyay_1986's three-photon table carries a doubt:
  Kulander 1988 reports a conversion-factor error of about 1e6 (TC-949).

## 0.5.0 (2026-10-09)

Two-photon ionisation of the noble gases: `optmat.gas.two_photon`, the generalised
two-photon ionisation cross sections of He, Ne, Ar and Xe (Xe also circular), each the
reference curve the author approved, with thresholds, labelled regions (NaN where
unavailable, never an error on a grid), uncertainty by kind and an identity (OEP 6,
Phase 1; TC-943, TC-944). The two-photon literature curated on the way: four
vocabulary quantities and 66 sources on two-photon ionisation, bound two-photon
excitation and excited-level photoionisation, compared on one footing in the manual
(TC-872 to TC-942). Krypton has no model; no existing value or digest changed.

### Added

- **Two-photon ionisation of the noble gases** (OEP 6, Phase 1; TC-943, TC-944):
  `optmat.gas.two_photon.two_photon_ionisation_model(gas, polarization)`, the generalised
  two-photon ionisation cross section of He, Ne, Ar and Xe (Xe also circular). Each
  model is the reference curve the author approved, with its two-photon thresholds, E1,
  labelled regions (zero below threshold; NaN where unavailable, never an error on a
  grid), its uncertainty kept by kind, `rate`, `loss_coefficient`,
  `pair_cross_section` and an identity. Krypton has no model. New exceptions
  `ArgumentError` and `UnavailableModelError`. No existing model, default or digest
  changes.
- The follow-up two-photon literature (TC-904 to TC-942): the originals the first round
  quoted, and the other calculations it named. Among them are two-photon ionisation of
  He (van der Hart and Bingham 2005; Bachau et al. 1990; Proulx and Shakeshaft 1993),
  Ne and Ar (Moccia et al. 1983; Pan, Gao and Starace 1990; Pindzola and Kelly 1975),
  Ar, Kr and Xe (Fink and Johnson 1990, relativistic, linear and circular) and Xe
  (Gangopadhyay et al. 1986; L'Huillier and Wendin 1987; Bischel et al. 1979; Hodges et
  al. 1981). There are also bound two-photon measurements and calculations (McCann et
  al. 1988; Khambatta et al. 1988; Bokor et al. 1980; Sundholm et al. 1994; and others),
  excited-level photoionisation (Chang and Kim 1982; Pindzola 1981; Hyman 1977; Hummer and
  Storey 1998) and xenon's Rydberg levels (Yoshino and Freeman 1985). Krypton still has
  no vector-drawn or tabulated two-photon ionisation cross section (OEP 6, Q2). The
  manual page's comparison is extended. OEP 6's author decided that the model will not
  raise on propagation grids, and that a time-resolved loss current is a documented gap
  (TC-903). The curated data are now 651 source records: 7785 datasets and 162,533
  values, 232 errata and 771 published formulas.
- Four quantities for two-photon data (TC-872): `two_photon_transition_strength`,
  `two_photon_excitation_cross_section`, `excited_state_photoionization_cross_section`
  and `multiphoton_ionization_polarization_ratio`.
- The two-photon literature of the noble gases (TC-873 to TC-901): the two-photon
  ionisation cross sections of He (Saenz and Lambropoulos 1999; Mihelič and Horvat 2021;
  Nikolopoulos and Lambropoulos 2001; Sato et al. 2011, measured), Ar (Petrov et al.
  2016; Starace and Jiang 1987), Ne to Xe (McGuire 1981, circular), Kr and Xe
  (L'Huillier, Tang and Lambropoulos 1989; L'Huillier and Wendin 1987) and Xe (Lagutin et
  al. 2017, linear and circular; Karamatskou and Santra 2017; McCown et al. 1982,
  measured; Santra and Greene 2004); bound two-photon strengths and cross sections (Hättig
  et al. 1998; Henriksson et al. 2005; Shekhtman et al. 2020; Kröll and Bischel 1990;
  Drag et al. 2021; Dakka et al. 2018) and excited-level photoionisation (with Geohegan
  et al. 1986); and records of the papers that hold nothing to curate. The curated data
  are now 613 source records: 7641 datasets and 155,228 values, 226 errata and 771
  published formulas.
- A manual page, *Two-photon ionisation and absorption: the data*, comparing them on
  one footing, with the thresholds of each gas (TC-902). optmat has no two-photon model
  yet; that is OEP 6's.

### Changed

- The analysis tool behind the third-order manual's tables,
  `tools/fit/nonlinear_support.py`, reads each gas's resonance limits from its support
  record, so `data/gas/nonlinear_support.toml` is the one copy of the selections (OEP 5,
  Q11; TC-871). No number in the tables changed.

## 0.4.0 (2026-10-07)

Where the hyperpolarizability fits are supported: each model's `support`, a public
record of the largest fitted W2, single leg and pair sum and of the gas's resonance
limits, with their sources, for propagation codes that evaluate every four-wave quartet
leg by leg (OEP 5, Final; TC-862 to TC-870). New sources on ionisation and the third
order: many-electron calculations against one electron, R-matrix multiphoton cross
sections, ion yields, saturation intensities, xenon's negative n2 and helium's dynamic
polarizability (TC-852 to TC-861).

### Added

- **OEP 5 is Final** (OEP 5, *A public record of where the hyperpolarizability fits
  are supported*): `HyperpolarizabilityModel.support` and its records (TC-870).

- `HyperpolarizabilityModel.support`, a public record of where each gas's
  hyperpolarizability fit is supported, for propagation codes (OEP 5, Phase 1;
  TC-870). It holds the largest fitted W2, single leg and pair sum, with the fitted
  points that attain them, and the gas's resonance limits E1, E2 for one linear
  polarisation and for any, and the ionisation energy, each with its datasets, kind,
  qualifications and uncertainty; absorption evidence is marked as such and a missing
  limit has its reason. It also gives two derived quantities, `leg_max_in_w2` and
  `band_edge()`. Wavenumbers are in m^-1. The record states facts and carries no rule.
  It has its own digest and leaves every model's digest and evaluated gamma unchanged.
  New classes `HyperpolarizabilitySupport`, `SupportLimit`, `SupportPoint`,
  `SupportUncertainty`, `SupportCandidate` and `AbsorptionEvidence`; the selection
  `data/gas/nonlinear_support.toml`, the stored records
  `data/gas/models/nonlinear_support/` and the tool
  `tools/fit/nonlinear_support_record.py`.

- Where the hyperpolarizability fits are supported, for propagation codes that evaluate
  every four-wave quartet leg by leg (solprop's SEP 6): the legs of every fitted value
  and validator (TC-864); the resonance limits of the 15 gases, a single leg at the
  lowest dipole-allowed level and a pair sum at the lowest two-photon level (TC-865); the
  fits' increase with each squared leg (TC-866); per-leg caps against the data, a Kerr
  carrier and solprop's demonstration quartets (TC-867); and the quantities derived from
  the facts, sqrt(3 W2,max / 4) and min(E1, E2/2), with their use left to the
  propagation code (TC-868). New manual
  sections on *Third-order nonlinearity* and *Third-order accuracy*, and the tool
  `tools/fit/nonlinear_support.py`. The noble gases' NIST ASD levels (TC-862) and Huber &
  Herzberg's states of H2, D2, N2, O2 and CO from the NIST WebBook (TC-863), with a new
  quantity, `electronic_state_energy`. The curated data are now 585 source records: 7441
  datasets and 113,554 values, 209 errata and 771 published formulas.

- `travers_2026_trecx`: haCC's argon in a 20 fs pulse and with mMax 1, and one electron
  at its ionisation energy in both pulses: the atom about 0.8 of one electron, rising by
  3 to 8 % from 6 to 20 fs (TC-861). The curated data are now 585 source records: 7402
  datasets and 113,515 values, 209 errata and 771 published formulas.

- `pathak_2020`, `pathak_2022`, `sato_2016`, `sato_2018`: time-dependent many-electron
  calculations of Ar, Kr and He in short sin^2 pulses, read from the figures (TC-859);
  `travers_2026_trecx`: one electron in the same pulses (TC-858); and
  `published_atom_ratios`: argon 0.87 to 1.10 and krypton 0.91 of one electron, helium
  1.53 to 1.71 (TC-860). The curated data are now 585 source records: 7392 datasets and
  113,505 values, 209 errata and 771 published formulas.

- `mckenna_2004`: McKenna & van der Hart's R-matrix Floquet two- and three-photon
  cross sections of Ne and Ar, read from the vector figures. Argon's at 193 nm is
  matched within 26 % by the default with the sum over channels, and 1.5 to 12 times
  over with its integral, which gives the multiphoton yield a non-integer order
  (TC-856).
- `lai_2017`: Lai et al.'s ion yields of Ne, Ar, Kr and Xe from 0.4 to 3.9 um, linear
  and circular, read from the vector figures; their intensities partly fitted to PPT,
  so a test of the yields' shape (TC-857). The curated data are now 581 source records:
  7358 datasets and 113,427 values, 209 errata and 771 published formulas.

- `uiterwaal_1998`: the saturation intensities of He, Ar and Kr at 248.6 nm, 500 fs,
  from the ponderomotive shift, and their multiphoton cross sections; the quantity
  `saturation_intensity`; `saturation_test`. In this multiphoton regime (Keldysh
  parameter 1.8 to 7) the default's rate is 0.55 to 4.3 times what they imply, the
  static Coulomb factor's 5.7 to 32 (TC-853).
- `lehmberg_1995`: xenon's negative n2 1 nm above its two-photon resonance at 248.4 nm,
  measured and calculated, and CO2's (TC-854).
- `masili_2003`: helium's dynamic polarizability to 0.845 a.u., through the 1s2p
  resonance; the library's lost minus signs restored (TC-855). The curated data are now
  579 source records: 7335 datasets and 102,921 values, 209 errata and 771 published
  formulas.

- `tolstikhina_2014`: the Hartree-Fock weak-field coefficients of He, Ne, Ar, Kr and Xe
  (Table I), and `weak_field_factors`, the atom's weak-field rate over the model's one
  electron: 1.7 to 1.8 for He, Ar, Kr and Xe, the literature's factor 2 less the excess
  of Tong, Zhao & Lin's coefficients; 0.53 for neon, unexplained (TC-852). The curated
  data are now 576 source records: 7322 datasets and 102,878 values, 209 errata and 771
  published formulas.

### Fixed

- *Third-order nonlinearity* said a calculated THG enters the fits with harmonics at
  217 nm and longer; the 650 nm window keeps 694.3 nm and longer, harmonics at 231.4 nm
  (TC-864). The fits are unchanged.
- *Third-order accuracy*'s table of solprop's demonstration quartets showed another
  table's rows (TC-868).
- The support analysis is qualified after OEP 5's review: the band edge min(E1, E2/2)
  is a sufficient bound, an SPM resonance only for one polarisation; the calculations
  beyond the data prove no general lower bound; the polyatomic molecules' absorption
  numbers are evidence, not levels (TC-869).

## 0.3.1 (2026-10-06)

Xenon in Wallace et al.'s pulse (TC-849) and xenon's effective number of electrons 1.5
(TC-850); krypton's 1.35 rests on Wallace et al.'s yields, Santra et al.'s rate being an
independent-electron sum, and the effective number of electrons is a calibration of the
formula, not a count of active electrons (TC-851).

### Fixed

- Santra et al.'s (2006) krypton rate is an independent-electron sum in a
  Herman-Skillman potential, not an all-electron calculation, as 0.3.0's notes and
  krypton's reasons said; krypton's 1.35 electrons rest on Wallace et al.'s yields. The
  manual states that the effective number of electrons calibrates the formula and is not
  a count of active electrons (TC-851).

### Changed

- Xenon's effective number of electrons in the default ionisation rate is 1.5 (was 1,
  provisional), from Wallace et al.'s yields through the one-electron calculation in
  their pulse; its uncertainty factor is 3, for the other absolute measurements, which
  put xenon high against argon at 800 nm (TC-850).

### Added

- `travers_2026_trecx`: one electron of xenon in Wallace et al.'s pulse at phases 0
  and pi/2. Through their test it matches at 1.019, within the calibration, as argon's
  does: the xenon atom ionises like about 1.1 electrons without the spin-orbit
  splitting, about 1.5 in the model's terms. Xenon's effective number of electrons
  stays 1, provisionally (TC-849). The curated data are now 575 source records: 7302
  datasets and 102,858 values, 209 errata and 771 published formulas.

## 0.3.0 (2026-10-04)

The noble gases' default ionisation rate takes the frequency-dependent Coulomb factor,
validated against one-electron calculations of H, He, Ne, Ar, Kr and Xe, and krypton's
effective number of electrons rises to 1.35 (TC-848). Argon's rate at 1e14 W/cm² and
800 nm falls by 24 %; the former default is kept as `default_static`. New one-electron
and all-electron calculations from tRecX and haCC (TC-842 to TC-847).

### Changed

- **The default ionisation rate of the noble gases takes the frequency-dependent
  Coulomb factor** (`coulomb = "frequency"`, c = 5.25; TC-831), as does `one_electron`;
  krypton's effective number of electrons is 1.35 (Santra et al.'s all-electron rate,
  Wallace et al.'s yields); xenon's stays 1, provisionally. The former models are
  `default_static` and `one_electron_static`; the provisional `frequency` and
  `one_electron_frequency` are retired. Argon's rate at 1e14 W/cm2 and 800 nm falls by
  24 %. The manual explains the rate as three parts, each tested against what isolates
  it (TC-848).
- The table of van der Hart's all-electron rates at 390 nm gains the candidate
  `one_electron_frequency`: with it one electron of argon is about 1.15 times the
  all-electron rates below gamma = 2.2, the atom at about 0.87 of one electron, as
  haCC's interim all-electron calculation at 800 nm (TC-842).

### Added

- `travers_2026_trecx`: the ionised fraction of atomic hydrogen at 400 and 800 nm in a
  15 fs pulse, converged to 0.06 %, among the hydrogen comparisons: at 800 nm Zhao et
  al.'s, with the same pulse, are 0.96 times it on average (TC-843). The curated data
  are now 575 source records: 7280 datasets and 102,717 values, 209 errata and 771
  published formulas.
- `travers_2026_trecx`: argon's one-electron ionised fraction at 400 and 1030 nm with a
  30 fs pulse, not fitted by the candidate `frequency` factor, which is 0.87 and 0.95
  times it on average (the static factor 2.0 and 1.18) (TC-844). The curated data are
  now 575 source records: 7284 datasets and 102,752 values, 209 errata and 771
  published formulas.
- `travers_2026_trecx`: one electron of krypton and xenon in A2's 6 fs pulse and their
  static rates: the candidate `frequency` factor is 1.16 (Kr) and 0.94 (Xe) times the
  calculations, krypton's excess its static normalisation (TC-845). The curated data are
  now 575 source records: 7288 datasets and 102,782 values, 209 errata and 771
  published formulas.
- `travers_2026_trecx`: one electron of krypton in Wallace et al.'s pulse at phases 0
  and pi/2; through their test it needs the intensities scaled by 1.045, beyond the
  calibration, where argon's matched: the krypton atom ionises like about 1.35
  electrons, as Santra et al.'s all-electron static rate also says (TC-846). The
  curated data are now 575 source records: 7290 datasets and 102,806 values, 209
  errata and 771 published formulas.
- `travers_2026_trecx`: haCC's all-electron argon in A2's 6 fs pulse and one electron at
  haCC's ionisation energies; `hacc_ratios` in `optmat.gas.ionisation_validation`: the
  argon atom ionises like 0.76 to 0.87 of one m = 0 electron at the same ionisation
  energy (TC-847). The curated data are now 575 source records: 7300 datasets and
  102,834 values, 209 errata and 771 published formulas.

## 0.2.2 (2026-10-02)

solprop SEP 3's findings: the ionisation rates three to four and a half times faster
with the same numbers (TC-837); a smooth derivative of the solids' indices made absolute
(TC-838); D2's viscosity from its correlation (TC-839); the Raman response's status
without the response (TC-840); and fused silica's dispersion near 2 µm explained
(TC-841). Also hydrogen's one-electron calculations (TC-833 to TC-835) and argon's in
Wallace et al.'s pulse (TC-836).

### Added

- One-electron calculations of atomic hydrogen from the figures: `bauer_2018` (Bauer
  and Deng; a 2-10-2 cycle flat top at 400 and 800 nm, and against the frequency at
  two fixed fields, the new quantity `ionization_yield_spectrum`; TC-834) and `le_2016`
  (Le et al.; sin^2 pulses at 400, 800 and 1600 nm; TC-835); Tong and Lin's alpha and
  I_p for hydrogen (`tong_2005`). `hydrogen_model`, `hydrogen_yields`,
  `hydrogen_spectra` and `cycle_pulse` (sin^2 and flat-top pulses for `pulse_yield`)
  in `optmat.gas.ionisation_validation`; a table on the manual page (TC-833). Below
  gamma = 1 the static factor is 1.05 to 1.32 times every set, the candidate
  `frequency` factor 0.92 to 1.03. The curated data are now 575 source records: 7274
  datasets and 102,649 values, 209 errata and 771 published formulas.
- `travers_2026_trecx`: argon's one-electron ionised fraction in Wallace et al.'s 5.5 fs
  pulse at carrier-envelope phases 0 and pi/2 (m = 0 and m = 1), and
  `tdse_intensity_scale_test`, which puts it through Wallace et al.'s test: it matches
  their measured yield at an intensity scale of 1.002 (TC-836). The curated data are
  now 575 source records: 7278 datasets and 102,697 values, 209 errata and 771
  published formulas.
- `raman_response_status(gas)` in `optmat.gas.raman_response`: the status and reason a
  gas's or a mixture's Raman response would have, without building it, or None where
  no constituent is Raman-active (the noble gases). It shares its rule with
  `raman_response` (solprop SEP 3, Q42; TC-840).

### Changed

- The ionisation rates are three to four and a half times faster, with the same
  numbers to the last bit: each rate of the average over a cycle is evaluated alone at
  its own nodes, the nodes once, `logsumexp`'s arithmetic without SciPy's dispatch,
  and the functions of the Keldysh parameter share their transcendentals. Argon's
  default at 800 nm takes 0.2 ms a call and 3.4 us a field, from 0.8 ms and 11 us; N2's
  molecular default 0.4 ms and 9 us. `tools/benchmark_rates.py` times them against
  another tree (solprop SEP 3, Q46; TC-837).
- D2's viscosity is Assael et al.'s correlation (their hydrogen correlation with the
  D2 mass, times the D2 paper's quantum factor), not a cubic spline through their Table
  2. The spline wiggled at 1e-4 (the table's rounding and its degraded 90 K cell), so
  no smooth series could represent it to better than 4e-4. The correlation reproduces
  the table within its rounding (5.1e-4; 1.2e-3 at 90 K), and a Chebyshev series of
  degree 20 represents it to 1e-11. The served values move by up to 5.1e-4 (1.2e-3 at
  90 K); the form `"table"` is removed (solprop SEP 3, Q33; TC-839).
- The silica page states where fit/SiO2's structure near 2 µm comes from: Lithosil
  Q2's OH bands, which the fit keeps unscaled. They give +2.1 % of the GVD at 2 µm over
  Malitson's formula, and the 2.21-2.25 µm group gives 7 % of β3. The absolute
  refractometry can neither confirm nor exclude them. The page also gives the model's
  1e-15 evaluation noise. No values changed (solprop SEP 3, Q37, Q32; TC-841).

### Fixed

- The derivative of a solid's index made absolute (Daimon's CaF2 relative to
  nitrogen, the air-relative formulas) no longer jitters by 1e-8 between neighbouring
  frequencies. The reference gas's refractivity is differenced, not its squared index.
  The derivatives move by up to 2.5e-8, the old error, and the indices are unchanged.
  `FormulaModel.EVALUATION_VERSION` is raised to 2, so the formula models' digests
  change (solprop SEP 3, Q38; TC-838).

## 0.2.1 (2026-10-02)

The equations of state's molar mass, which solprop SEP 3 needs for Gas.Air (TC-832);
argon's one-electron ionised fraction at 800 nm with a 10 fs pulse (TC-830); and a
candidate Coulomb factor that depends on the frequency, as provisional models that are
not the default (TC-831).

### Added

- `travers_2026_trecx`: argon's one-electron ionised fraction at 800 nm with the 10 fs
  pulse of the 400 and 1030 nm sets (m = 0 and m = 1): at gamma = 1.8 the default's
  static Coulomb factor is 1.1, 1.7 and 1.9 times the calculation at 1030, 800 and 400
  nm with one pulse (TC-830). One-electron runs at the all-electron calculation's
  ionisation energies are archived for its comparison. The curated data are now 573
  source records: 7264 datasets and 97,227 values, 209 errata and 771 published
  formulas.
- A candidate Coulomb factor that depends on the frequency, `coulomb = "frequency"`:
  the static factor times (1 + 2 gamma/e)^(-2 n* c s), s = n* omega/2 I_p, c = 5.25
  fitted to one-electron calculations of argon at 400, 800 and 1030 nm; the
  provisional models `frequency` and `one_electron_frequency`. Not the default
  (TC-831).
- `EquationOfState.molar_mass` and `MixtureEquationOfState.molar_mass`, kg/mol: the
  molar mass the backend's equation uses (equal to the transport models'), and a
  mixture's composition mean, for a propagation code's mass terms. Air's is CoolProp's
  pseudo-pure 28.96546 g/mol, not the 28.9586 g/mol Lemmon et al. (2000) print
  (TC-832). No identity digest changes; the golden values add the molar masses.

## 0.2.0 (2026-10-02)

The ionisation rates validated against one-electron calculations at 400, 800 and 1030 nm,
measurements of the plasma's phase shift and Luna.jl and solprop (TC-812 to TC-829); a
rate for circular polarisation; the rates' average over a cycle by quadrature, held
above the maximum field (TC-828), which changes the default rates near and above it.

### Added

- **`travers_2026_trecx`: the 6 fs one-electron yields as the ionised fraction**, with
  convergence runs (`ionised_gauss6_Ar_m0`, `ionised_gauss6_Ar_m1`, `ionised_gauss6_He`;
  TC-812), and the static rates of Ar 3p, m = 0 and m = 1 (`static_Ar_m0`,
  `static_Ar_m1`; TC-813). The validator's 6 fs yields are the ionised fraction.
- **The rates of Luna.jl and solprop against the default** (TC-815): a fixture of both
  codes' rates and pulse yields at 800 nm, their tests on the measured yields, and a
  section of *Strong-field ionisation rates*.
- **Absolute ionisation measured by plasma interferometry**: `nie_2022` (He, Ar, Kr,
  Xe and H2; TC-818) and `nambu_2026` (Ar, Kr, Xe at 800 and 400 nm, linear and
  circular polarisation; TC-817), read from their figures, as the new quantity
  `ionization_phase_shift`; `ionisation_validation.phase_shift_test`,
  `plasma_phase_shift` and `circular_yield` compare a rate with them, and the manual
  states the measured yields as intensity scales (TC-819).
- **The rate in a circularly polarised field**, `IonisationRateModel.circular_rate`
  (with `circular_enhancement` and `circular_saddle_point`): each channel's static
  rate times Barth & Smirnova's (2013) factor for the turning field, summed over the
  shell's orbitals; it fits the measured circular-polarisation yields at the intensity
  scales of the linear ones (TC-820).
- **Static rates of Ne, Ar, Kr, Xe**: `bahl_2016` (one electron, from the publisher's
  data files), `santra_2006` (Kr, spin-orbit resolved), `otobe_2004` (Xe, density
  functional theory); Bahl et al.'s and Santra et al.'s are static references of the
  validator (TC-821).
- **Molecule-to-atom yield ratios and O2's eight-photon cross-section**: `dewitt_2001`,
  `loriot_2008`, `dura_2012` (O2/Xe, N2/Ar; the new quantity `ionization_yield_ratio`,
  with each source's detection convention) and `sharma_2018` (TC-822).
- `lin_2012`: O2/Xe and N2/Ar at 800 and 2000 nm, the measured rise of O2/Xe with the
  wavelength that the default lacks (TC-823). The curated data are now 572 source
  records: 7248 datasets and 96,899 values, 209 errata and 771 published formulas.
- `travers_2026_trecx`: argon's one-electron ionised fraction at 400 nm (10 fs, m = 0
  and m = 1), a reference of the validator: the default's static Coulomb factor is 1.4
  to 3.9 times it at Keldysh parameters 1.4 to 4.2 (TC-824). The curated data are now
  572 source records: 7250 datasets and 96,919 values, 209 errata and 771 published
  formulas.
- `zhao_2016`: one-electron TDSE ionisation probabilities of atomic hydrogen at 600 to
  1200 nm and of He, Ne (400 nm), Ar, Kr, Xe (800 nm); the material `H` (atomic
  hydrogen) for such reference calculations. The static Coulomb factor is 1.75 to 2.2
  times high for krypton and xenon at 800 nm (TC-825), and 1.09 to 1.53 for hydrogen
  at 1200 to 600 nm with Tong and Lin's static correction (TC-827; TC-825 had 2 to 2.6,
  without it).
  The curated data are now 573 source records: 7260 datasets and 97,195 values, 209
  errata and 771 published formulas.
- `travers_2026_trecx`: argon's one-electron ionised fraction at 1030 nm (10 fs, m = 0
  and m = 1): the default's static Coulomb factor is 1.02 to 1.35 times it at Keldysh
  parameters 0.4 to 1.8, against 2.0 at 400 nm with the same pulse at gamma = 1.8
  (TC-826). The curated data are now 573 source records: 7262 datasets and 97,210
  values, 209 errata and 771 published formulas.
- The rates of Luna.jl and solprop through every test at 800 nm (Wallace et al.,
  Nambu et al., Nie et al., Wahlstrand et al., Parker et al.'s helium), a table on the
  manual page beside optmat's default (TC-829).

### Changed

- The one-electron model's agreement with the one-electron calculations is restated:
  1.07 to 1.67 times their ionised fraction for argon at 6 fs, not 0.93 to 1.21 against
  the population beyond 20 bohr (TC-812). The default model is unchanged and, until
  the measured yields agree, follows Wallace et al.'s (TC-814).
- The default ionisation rate's recorded reasons (its `coulomb` choice and argon's
  accuracy) are restated with the new comparisons, which changes the identity digest of
  `rate/default/*`; the golden values are re-recorded, the values unchanged (TC-812).
- The accuracy reasons of the default's Kr and Xe now cite their static references
  (TC-821), which changes the identity digest of `rate/default/*` again; the golden
  values are re-recorded, the values unchanged.
- `zhang_2014`: the excluded row of Table IX is very likely displaced cells; the
  Tong-Lin value is printed in the Muller column (TC-813).
- The rates' average over a cycle is computed by quadrature instead of an expansion in
  the slope and curvature of the rate, which was up to 24 % low near the maximum field,
  made the rate turn over just above it and left rounding noise of 3e-9. Above the
  maximum field a statically corrected rate is now held at its value there. Each
  molecule's maximum field is the largest field its alpha was fitted to or checked
  against (1.61 to 1.96 times the barrier-suppression field; `correction_range`).
  Golden values re-recorded (TC-816, TC-828).
- The `poe` test and coverage tasks run on every core with `pytest-xdist`: the default
  lane and the coverage gate in about 45 s instead of about 4 minutes; the library and
  fits lanes test by test. A plain `pytest` stays serial.

## 0.1.1 (2026-10-01)

OEPs 2 and 3 Final, the lazy imports, and the equation of state's domain fixed at the
melting line (TC-811), which solprop SEP 3 needs.

### Added

- **OEP 3 is Final** (OEP 3, *The molecular gases between their infrared bands*): the
  band rules and windows, several intervals per model, the refitted molecules and the
  settled questions (TC-747 to TC-756, TC-770, TC-797 to TC-809).
- **OEP 2 is Final** (OEP 2, *Compact infrared terms for the molecular gases*): the
  compact models of the molecules and their printed forms (TC-799 to TC-801).
- **Compact models of the molecules** (OEP 2, Phase 1; TC-799): `compact/CO`,
  `compact/CO2`, `compact/CH4`, `compact/N2O` and `compact/SF6`, each its fit's
  electronic terms and covariance with three to five effective infrared poles in place
  of thousands of band-rule poles, within 0.1 of the fit's sigma in the polarizability,
  the GVD and the TOD over its optical intervals. Listed by `registry.models`, never a
  default. `tools/fit/compact.py` re-derives them after a refit.
- `FittedModel.dispersion_sigma_au`, the absolute one-sigma of a dispersion derivative.
- **A mixture's intervals widen across a trace constituent's gaps** (OEP 3, Q13; TC-809):
  where the mixture passes OEP 3's absorption, pole and representation checks at its
  mole fractions, from each molecule's new gap table (`<gas>_gap_bound.csv`,
  `tools/hitran/gap_bound.py`). The full CIPM-2007 air gains 1499-1516 and 1548-1562 nm.
- **Infrared and static validators** (OEP 3, Q19 and Q20; TC-802 to TC-807): Rollefson and
  Havens's CH4 at 1.68-14.8 um, Matsumoto's CO2 at 10.6 um, Tsankova et al.'s static
  CO, and Old et al.'s 1695 and 1817 nm CO2 points, compared at their fitted scale; Legay
  (1955, HCl) and Mitrofanov et al. (2017, air's group delay) recorded without datasets.
  A validator may take a wavelength window. The fits are unchanged. The curated data
  are now 562 source records: 7195 datasets and 95,846 values, 209 errata and 771
  published formulas.
- `tools/hitran/resolution.py`: the fits' representation error inside their optical
  ranges at coarser resolutions (OEP 3, Q18; TC-808).
- **The compact models' printed forms** (OEP 2, Phase 2; TC-801):
  `optmat.gas.models.compact.coefficients` in the polarizability form or as a Sellmeier
  form at a state, with the qualified intervals, the electronic covariance in the printed
  parameters and the Sellmeier form's error against the full model; `format_coefficients`
  prints them to eleven digits. The tables are on *The gas models*, and a gallery example
  compares each compact model with its fit.

### Changed

- Importing the gas and solid models no longer loads scipy.signal, scipy.stats,
  scipy.ndimage, scipy.interpolate or scipy.optimize; they load on first use, which
  takes `optmat.gas.models.registry`'s import from about 1.2 s to 0.4 s (TC-810).
- The first use of a model no longer loads the whole source catalogue (562 records,
  about 2 s): a lookup builds only the registry its name belongs to, for gases and
  solids, and the models read only the records they cite. `silica_bundle()` and SiO2's
  index take 0.4 s instead of 2.1 s from a fresh process, the gas `default_model("Ar")`
  1.2 s instead of 3.0 s (TC-810).

### Fixed

- The equation of state's inverse, `EquationOfState.pressure`, no longer fails at every
  density above the critical temperature where CoolProp's melting line lies below the
  equation's pressure: nitrogen to 284 K, argon to 254 K, methane to 256 K, CO2 to
  328 K, air to 236 K, and helium, neon, hydrogen and deuterium at low temperatures.
  `phase_domain` now ends just below the melting pressure there. The forward law refuses
  states above it as solid, and the domain excludes each equation's triple point. No
  value the forward law or the inverse returned before has changed (TC-811).

## 0.1.0 (2026-10-01)

The first version, tagged locally once OEP 1's phases and OEP 3's refits had merged
(OEP 1, section 6.12).

### Fixed

- Incomplete resolved Raman configurations now fail explicitly instead of silently
  becoming unresolved bands (TC-735). Invalid EUV continuum-constraint steps raise
  `DataError` before division by zero or logarithms of non-positive densities
  (TC-736). Existing curated Raman and EUV outputs are unchanged.
- The dispersion fit gallery now applies the fitter's electric/refractive and infrared
  conventions, wavelength windows, additive corrections and uncertainty options to
  its direct-data comparisons (TC-729). Argon's static dielectric point is shown at
  its actual residual, about 0.41 sigma, instead of the spurious 8.99 sigma; fitted
  error bars use the uncertainties of the fit. Stored models and fits are unchanged.

### Changed

- **Dry air's mixture equation of state** (TC-796; OEP 4, Q19) is the multi-fluid model
  at every temperature, provisional, in place of the fallback that held only at
  290-320 K. It agrees with CIPM-2007's compressibility within 8.6e-6.
- **The declared scale of an uncurated like term** (TC-795; OEP 4, Q17) is 2.5 cm^3/mol
  (was 2), the largest curated |b|. It changes only domain warnings and the mixtures'
  sensitivity.
- **Oxygen's refractivity virial coefficient** (TC-794; OEP 4, Q18): O2 takes Egan et
  al. 2026's published estimate, b_R = 0.210 cm^3/mol at 633 nm with 50 %
  uncertainty, as a provisional like term (it had none). O2's n - 1 at a state rises by
  8.7e-6 at 1 atm and 9.1e-4 at 100 bar; air's optical domain widens from 1300 to 6150
  mol/m^3.
- **Krypton's refractivity virial coefficient** (TC-793; OEP 4, Q18): krypton takes
  Achtermann et al. 1993's optical B_R, 5.96(6) cm^6/mol^2, in place of Orcutt and
  Cole's static 5.6(3). Krypton's n - 1 at a given state rises by 1.4e-6 at 1 atm and
  1.9e-4 at 100 bar; its polarizability moves by 1.2e-7 at 800 nm, within the fit's 78
  ppm.
- **Neon's refractivity virial coefficient** (TC-792, TC-744; OEP 4, Q18): neon takes
  Achtermann et al. 1993's optical B_R, -0.11(2) cm^6/mol^2, in place of Orcutt and
  Cole's static -0.30(10). Neon's n - 1 at a given state rises by 7.9e-6 (relative) at
  1 atm and 7.5e-4 at 100 bar. Its polarizability moves by at most 1.7e-8, and the
  refit's sigma at 800 nm is 8.6 ppm (was 8.7).
- **The molecules refitted with their bands as rules and their windows** (OEP 3,
  Phase 3; TC-797). CO, CO2, CH4 and N2O now:
  - carry every band as a Gaussian rule of HITRAN's lines, and the polar molecules
    their rotational band;
  - hold on several qualified intervals, their optical range split only where the
    pure gas absorbs, plus the windows between the bands.

  CH4 now ends at 1313 nm (was 1400) and N2O at 1950 nm (was 2000). Inside the new
  intervals the polarizability moves by at most 0.05 (CO), 0.26 (CO2), 0.67 (N2O) and
  1.9 (CH4) of the present fits' sigma. CO2's 9.2-10.7 um values (Marchetti) become
  validators, since the model is not qualified there. The static polarizabilities of CO and N2O gain
  their orientation part (+0.65 and +1.45 au).
- Models own their storage (OEP 1, Phase 1): their arrays are copied and read-only,
  and their mappings are read-only `FrozenMapping`s, so changing an array after
  passing it to a model changes nothing, and writing to a model's array raises. A gas
  formula's reference density, the virial ratio and the models silica's Raman kernel
  reads are resolved once per model. `GasRamanResponse` gains `material` and
  `population_floor`.
- The development gate now requires 100% statement and branch coverage in the default
  test lane, with the existing package scope and exclusions. Diagnostic coverage
  commands remain available without thresholds. Agents prefer topic branches and
  separate worktrees, confirming creation or switching before proceeding.
- Branch-only register numbers are explicitly provisional. `poe register-sync`
  previews and applies their renumbering and citation updates before integration,
  and can resolve the subsequent register append conflict without changing entries
  already on `main`.
- Removed duplicate private checks for dataset units and wavelength media already
  rejected by dataset construction. Public runtime guards remain and are tested in
  ordinary Python processes; coverage reports now include those subprocesses.
- The dedicated comparison pages and Luna comparison figure are retained as private
  development notes, outside the published manual and gallery. Numerical reference
  checks remain in the test suite.

### Added

- **OEP 1 is Final** (OEP 1, *The propagation-code interface*): identity, overrides,
  the equation of state, transport, air, the Raman groups, the rate energies, the
  silica bundle, the solid tables, solprop's baseline, the page for propagation codes
  and the release test (TC-745 to TC-789).
- **OEP 4 is Final** (OEP 4, *The refractive index of gas mixtures at finite density*): the
  refractive index of gas mixtures at finite density, its optics, density and example,
  with the author's decisions on its open questions (TC-771 to TC-796).
- **The release test** (OEP 1, Phase 5; TC-789): `poe test-release` builds the wheel
  and sdist, installs each offline into a clean environment, and runs the quickstart,
  the propagation page and the public surface there. The sdist no longer picks up local
  tools' state under `.claude/`.
- **The page *Using optmat from a propagation code*** (OEP 1, Phase 4; TC-788): what a
  propagation code calls, in which units, what to store and where a model holds, with
  snippets the tests run.
- **solprop's baseline at one commit** (OEP 1, Phase 3; TC-787): every solprop fixture
  re-recorded at solprop d577414, with no value changed. New transport and density
  fixtures are compared with optmat's on *Transport of the dilute gases*. solprop
  returns a liquid's density for SF6 at 30 bar, where optmat refuses the state.
- **A helium-argon example and the combined sensitivity** (OEP 4, Phase 3; TC-786):
  the gallery figure of a He/Ar fill from 1 to 60 bar, and a table on *Gas mixtures at
  finite density* of what the local field, the pair terms, the density rule and the
  density method each change in n - 1 and the GVD.
- **The density of gas mixtures** (OEP 4, Phase 2; TC-778 to TC-785):
  - `optmat.gas.mixture_eos.mixture_equation_of_state(composition)` gives a mixture's
    equation of state by one of three methods: CoolProp's multi-fluid model, the
    virial fallback (each species' own equation at its partial density plus measured
    cross second virials), or the partial-density rule. It provides `number_density`,
    `partial_densities` (the optics' keys), `pressure`, `compressibility`,
    `second_virial`, its pairs, method, status and identity.
  - The multi-fluid models were compared with the measured cross, excess and mixture
    virials and the mixture densities, under a rule fixed before the comparison
    (`optmat.gas.mixture_validation`). The He-Ar, Ne-Ar and He-Ne models (Tkaczuk et
    al. 2020) and the air pairs' (Gernert 2013) fail, and their mixtures take the
    fallback: He/Ne validated; He/Ar, Ne/Ar and dry air provisional. Dry air's
    fallback agrees with CIPM-2007's Z within 1.8e-6.
  - Seven source records: Brewer and Vaughn 1969, Martin et al. 1982, Vatter et al.
    1996, Vogl and Hall 1972, Kalfoglou and Miller 1967, Blancett et al. 1970 and
    Kestin et al. 1981. Bose and Cole 1970's pressure virials were added to its
    record.
  - A `thermodynamic` family of six quantities, and new pressure and molar-density
    units.
  - The density sections of *Gas mixtures at finite density*. No existing value
    changed.
- **The optics of gas mixtures at finite density** (OEP 4, Phase 1; TC-771 to TC-775):
  `optmat.gas.mixture.mixture_optics(species)` gives `MixtureOptics`, the
  Lorentz-Lorenz law of the whole mixture with its pair terms. It has the real
  `refractivity` and complex `susceptibility`, both derivatives,
  `refractivity_uncertainty`, every `PairTerm` with its status, the constituents'
  common intervals and the optical domain's `density_bound`. One species is the pure
  law exactly. No unlike-pair term is curated: each is the declared mean of its like
  terms. New quantity `unlike_pair_refractivity_virial_coefficient`. Achtermann et al.
  1991's C_R and D_R are extracted, and Achtermann et al. 1993 and Huot and Bose 1991
  are curated. The new page *Gas mixtures at finite density*. No existing value
  changed.
- **Several qualified intervals per fit** (OEP 3, Phase 2; TC-770): fit configurations
  gain `[model] band_nodes` (the derived band rules, and the rotational system as a term
  class of its own, left out of electric data) and `windows = true` (the derived
  windows as the qualified intervals, with each gap's reason in the range warning).
  `FittedModel.representation_bound_au` gives the bound on the representation error.
  No stored fit uses the keys yet, and every fit is unchanged.
- **Refusal outside a solid table** (OEP 1, Phase 2h; TC-769): `TabulatedModel` has
  `out_of_range`, `"hold"` by default as before. `with_out_of_range("raise")` makes a
  wavelength outside the knots raise `StateError`. The metal tables' digests moved;
  no value changed.
- **The silica nonlinear bundle** (OEP 1, Phase 2g; TC-768):
  `optmat.solid.bundle.silica_bundle(name)`, the consistent triple (n2, f_R, h_R) of
  fused silica for a code in the n0 = 1 convention. Bundles: `default` (Hollenbeck,
  paired), `tensor` (Lin and Agrawal, their f_R) and one per response. Each carries its
  gain check against Stolen and Ippen.
- **The energy per ionisation event** (OEP 1, Phase 2f; TC-767): `channel_rates`,
  `energy_loss_rate`, `energy_per_event` and `energy_kind` on the atomic and molecular
  rate models. The channel rates sum to the rate, and the energy per event is defined
  down to zero field. No rate changes.
- **The Raman groups** (OEP 1, Phase 2e; TC-766): `GasRamanResponse.groups()`, each
  damping group's width law in SI with `decay_rate(N)` down to zero density, and its
  share of the response; `line_id`, a stable identifier of each line (snapshot in
  `tests/data/raman_line_ids.json`); `check_densities(N)`; and the derivation of the two
  tensor parts on the manual page. The response's vibrational `band` is now a short key
  (`v1`, `nu1`, `fermi_1388`, ...), so the Raman digests moved; no value changed.
- **Transport follow-up** (TC-765): four new source records. D2's zero-density table
  is explained: the D2 paper printed the universal collision functional in place of the
  hydrogen correlation (TC-760). SF6's viscosity (0.3 %, 223.555-1000 K) and ideal-gas
  heat capacity are checked against their papers (TC-761). D2 and SF6 are now
  `validated`. O2's heat capacity stops at 1700 K, where CoolProp's ideal part leaves its
  source's stated accuracy.
- **Explicit parameters and overrides** (OEP 1, Phase 2a; TC-764):
  `model.with_parameters(**changes)` on every model of the surface. It rebuilds from
  the class's declared inputs, refuses derived fields, and records the changes and the
  curated digest in the identity, with the status `user`. New constructors:
  `IonisationRateModel.from_parameters`, `MolecularRateModel.from_parameters`,
  `OscillatorModel.from_parameters` (SI poles), `GasRamanResponse.from_lines`,
  `FormulaModel.sellmeier` and `TabulatedModel.from_arrays`. A species may be a custom
  name. The curated builders go through these constructors; no value or digest
  changes. `MolecularRateModel.with_ionization_energy` is now an override, so its
  result has the status `user`.
- **The dilute gases' transport properties** (OEP 1, Phase 2c, second part; TC-760 to
  TC-763): `optmat.gas.transport.transport_model(gas)` for solprop's twelve gases and
  CO2, with the dilute-gas viscosity, the ideal-gas heat capacity (CODATA's R; (5/2) R
  for the noble gases) and the molar mass, their ranges and identity, and Wilke's
  `mixture_viscosity`. Sixteen new source records (twenty with the follow-up) name each
  correlation and equation of state; CoolProp's correlations reproduce their papers' dilute forms to 1e-9. Ne, Kr,
  Xe and N2O take Huber's dilute form (checked against the WebBook's tables, now in
  `data/external/nist_webbook`), D2 Assael et al.'s table (`provisional`), SF6
  CoolProp's unchecked correlation (`provisional`). Two quantities,
  `dilute_viscosity` and `ideal_gas_heat_capacity`, in a new family, `transport`.
- **Named air compositions and the air model** (OEP 1, Phase 2d; TC-759):
  `"dry air (N2, O2, Ar)"` and `"dry air (N2, O2)"`, CIPM-2007 reduced and renormalised;
  `optmat.gas.air.air_model(name)`, the mole-fraction sum of the constituents'
  recommended models with its derivative, intervals and identity, and CoolProp's
  pseudo-pure dry air for its density (`air_equation_of_state`). The full composition
  agrees with Ciddor's standard dry air (a new record, `ciddor_1996`) within 1.1e-8 in
  n from 400 to 1400 nm. The Raman response accepts the two reduced names (`"dry air"`
  is the first) and refuses the full composition.
- **The equation-of-state object** (OEP 1, Phase 2c, first part; TC-758):
  `optmat.gas.eos.equation_of_state(gas)`, with its identity, forward law, speed of
  sound, phase domain and inverse, `pressure(N, T)`, round-tripping to 1e-12; and
  `optmat.gas.eos.pressure`. The real-gas density now refuses a condensed state or one
  beyond the equation, raising `StateError` (TC-731); no curated value changes.
- **The polarizability's frequency derivative** (OEP 1, Phase 2b; TC-757):
  `dpolarizability_domega(wavelength)` on every gas dispersion model, complex, in SI,
  analytic for every class, 0 in the static limit. Fitted models expose their
  reference temperature, the temperature of their infrared terms and the measured
  temperature coefficient (N2, O2, CO2) as metadata in their identity. O3's first
  frequency derivative (`dpolarizability_domega_au`, `dispersion_au` at order 1) is now
  the spline's own, which moves it by at most 4.5e-5 of its maximum.
- **The molecular gases' infrared windows, band rules and bound** (OEP 3, Phase 1;
  TC-747 to TC-756). `tools/hitran/windows.py` derives, from HITRAN's line lists, the
  Gaussian band rules, the folds of the weak in-window bands, the rotational systems,
  the qualified intervals from 12.5 um to each fit's shortest wavelength and the
  representation bound, for CO, CO2, CH4 and N2O, with their provenance; loaders in
  `optmat.gas.infrared` (`band_systems`, `infrared_windows`, `bound_cells`,
  `gauss_rule`, `fold_poles`). The optical ranges are split only where the gas absorbs.
  The fits do not use these tables yet, and no number changes. A new source record,
  Golub and Welsch (1969).
- **Model identity** (OEP 1, Phase 1; TC-745, TC-746). Every model, response and fit of
  the public surface has an `identity`: its name, the quantity it gives, its material and
  sources, its status (`validated`, `provisional`, `untested`, `comparison` or `user`,
  with a reason) and a SHA-256 digest of what it evaluates, stable across processes and
  releases that do not touch it. Dispersion models have qualified `intervals`, one per
  model, with `wavelength_range` their envelope. SF6's molecular rate is `untested`,
  and the CO2 and N2O Raman responses are `provisional`. The new page *Identity, status
  and validity* describes them. The public surface and golden values of every model are
  snapshotted in the tests. No number changes.
- Statement and branch coverage reports for the default, fit-reproduction
  and paper-library test lanes. New contract tests cover
  provenance and errata, model construction, weighted fits and uncertainty,
  Raman and nonlinear response, ionisation and spectral-loader boundaries.
- An About page describing AI-assisted development, paperextract, source verification
  and reproducible fits; clearer project motivation and acknowledgement of
  refractiveindex.info on the front pages.
- Links to LupoLab/optmat and an Edit on GitHub button in the manual.
- **Sources for the accepted OEPs** (TC-737 to TC-744), seven records of notes, no
  datasets: Wilke's mixture-viscosity rule (`wilke_1950`, checked against the page and
  reproduced from its Table I), Buckingham and Pople's dielectric virial
  (`buckingham_1955`), Garberoglio and Harvey's ab initio refractivity virials of He,
  Ne and Ar (`garberoglio_2020`), and the mixture equations of state GERG-2008
  (`kunz_2012`), Lemmon et al.'s air (`lemmon_2000`), Tkaczuk et al.'s noble-gas
  binaries (`tkaczuk_2020`) and Gernert's humid-gas and CCS binaries (`gernert_2013`).
  An open finding: the neon refractivity virial the package uses (Orcutt and Cole's
  static value) is about three times the ab initio and optical values (TC-744).
- **Dispersion sources for the gas models** (TC-704 to TC-724), 21 new source records:
  - helium theory: `puchalski_2020` (the alpha^3 QED static polarizability, 7.5 ppb
    above puchalski_2016) and `pachucki_2019` (the retardation term of the refractive
    index, a -24 ppm omega^2 term that fit/He does not yet include: open);
  - magnetic susceptibilities: `lesiuk_2023b` (Ne and Ar with the relativistic
    correction; the first `magnetic_susceptibility` data) and `lang_2026` (N2, with its
    ab initio Cauchy moments at 303 K, 0.16 % below fit/N2);
  - metrology: `egan_2023` (A_R of Ar and N2 at 633 nm on a piston-gauge scale, with
    its authors' correction `egan_2024_correction` applied as errata) and the review
    `garberoglio_2023`;
  - Rayleigh scattering: `thalman_2014` (N2, Ar, O2; with its 2017 erratum) and
    `wilmouth_2019` (Ar, CO2, SF6, CH4), which resolve the methane disagreement of
    TC-150 (a published formula 14 % high) and show SF6 about 2 % above its fit in
    the UV;
  - the vacuum ultraviolet of Kr and Xe: `kramer_1985` (xenon's negative refractivity
    at 115-117 nm, read from a scanned figure; a check), `laporte_1975`,
    `laporte_1981`, `gerardo_1974`, and the third-harmonic papers `garrett_1986`,
    `cotter_1979`, `junginger_1981`;
  - `chan_1993_erratum`, applied to chan_1992 (it moves only euv/Ar);
  - method only, for compact infrared terms: `voronin_2017`, `mathar_2004`,
    `mathar_2007`, `tang_1972`.
- **Refits with them** (TC-725): Ne's and Ar's relativistic chi_0, N2's first chi_0
  (ab initio), egan_2023 fitted in fit/Ar (-1.3 ppm at 633 nm; Kr and Xe follow), and
  new Rayleigh, theory and metrology validators for He, Ar, N2 and O2. fit/CO2 moved to
  a lower minimum with one oscillator on its 202 nm bound (GVD -0.35 % at 633 nm, +32 %
  at 250 nm; its 343 nm working estimate is now 2 %): open.
- **CO2's two weak bands below 10 eV as lines** (TC-726): Chan et al.'s 1Delta_u and
  1Pi_g strengths as priors, the effective oscillators from 10 eV. The ultraviolet form
  no longer follows a modelling bound; the fit returns to within 0.3 % of its GVD before
  TC-725, and Bideau-Mehu's held-out 181-254 nm values agree (chi-squared 17363 -> 57).
- **A reference temperature for the molecules** (TC-730): `[model]
  reference_temperature_k` and a dataset option `temperature` bring values measured
  elsewhere to it with a measured coefficient. N2, O2, CO2, N2O, H2 and D2 refer to
  293.15 K; Egan's 303 K values of CO2 and O2 are brought there (n - 1 -40 and -21 ppm).
- **Helium's retardation term** (TC-728): Pachucki & Puchalski's -3.73e-5 omega^2 a.u.,
  added like chi_0 to the refractive polarizability (`OscillatorModel.retardation_au`,
  `[model] retardation_au`): fit/He's GVD -24 ppm, and Yang et al.'s measured dispersion
  coefficient agrees (4.8 -> 0.6 sigma). fit/He also takes the Thomas-Reiche-Kuhn sum
  S(0) = 2 as a prior, which keeps its highest oscillator from running away.
- `travers_2026_trecx`: argon's one-electron ionised fraction in Wahlstrand et al.'s
  42 fs pulse, 17 intensities with convergence runs and the bound states subtracted,
  and `long_pulse_yield` to compare a model with it. The `one_electron` model is 1.24
  times it on average, through its resonance structure; the one electron already
  exceeds the measured yield of the atom by 1.3 to 4.1 times (TC-702).
- `travers_2026_trecx`: the 6 fs iSurff yields are marked as lower bounds. The runs keep
  the Coulomb tail at the tSurff surface, which made iSurff 1.4 to 23 % low in the 42 fs
  runs; with the tail switched off it agrees with the ionised fraction (TC-703).
- **Values read from figures** (data rule 8, TC-699): where a source prints no number,
  its figure may be read, from the PDF's vector drawing (the markers' and lines'
  coordinates, the axes calibrated from the drawing's ticks) or from an embedded image,
  by the source's extraction script (`tools/extract/_figure.py`). Such a dataset sets
  `digitised = "vector"` or `"raster"`, names the library figure, states its digitising
  uncertainty, and is checked against every number the source prints; values read from
  an image are never fitted (the loader rejects it). `pypdfium2` is a development
  dependency of the extraction scripts.
- `vanderhart_2007` and `vanderhart_2006`: van der Hart's all-electron R-matrix rates at
  390 nm, the time-dependent ones of Ar (Table I, and its Fig. 4's vector drawing) and
  the Floquet ones of He, Ne and Ar (read from the image of its Fig. 6, which matches
  the vector copy within 4 %); `rmatrix_rates` compares a model with them. Helium's
  Floquet rate agrees with Parker et al.'s; argon's default rate is 1.6 to 2.9 times
  above them near the tunnelling regime and 7 to 12 times in the multiphoton one,
  which is open (TC-700, TC-701).
- `jagau_2018`: coupled-cluster static ionisation rates of H2, CO, CH4 and O2 from the
  supplement of Jagau (2018), with their averages over orientations
  (`orientation_averaged_rates(..., reference="jagau_2018")`). O2's and CH4's correction
  at finite field is now fitted to them (alpha 6.1 and 18.0) instead of taking N2's 14,
  and CH4's rate takes the middle maximum of its Jahn-Teller band, 14.4 eV, where the
  first gave a dependence on the field no alpha reproduced. O2's rate at 1e14 W/cm2 rises
  2.8 times (TC-698).
- **Enhancement proposals.** Major changes to optmat are proposed, reviewed, implemented
  and recorded as optmat Enhancement Proposals (OEPs, `docs/oeps/`). There is a process
  document (OEP 0), a template, an archive and a consistency test
  (`tests/test_proposals.py`). Every OEP is reviewed by at least one human and one
  frontier AI model before it is accepted. OEPs 1 to 4 are accepted: the
  propagation-code interface that solprop's SEP 3 needs, the compact infrared terms,
  the molecular gases between their infrared bands, and the index of gas mixtures at
  finite density that solprop's SEP 11 needs.
- **The curated data:**
  - 557 source records: 7190 datasets and 95,742 values, 209 errata and 770
    published formulas, with the provenance of every value and the evidence of every
    erratum;
  - the extraction scripts that reproduce them from the papers' tables;
  - the review register.
- **Gases:** He, Ne, Ar, Kr, Xe, H2, D2, N2, O2, CO, CO2, CH4, SF6, N2O and O3.
  - A fitted model of the refractive index of each (dipole oscillator-strength models
    with HITRAN's vibrational bands for the infrared-active molecules, and a
    Kramers–Kronig model from the measured absorption for ozone), with covariance and
    leave-one-out robustness records.
  - Real-gas densities, Lorentz–Lorenz with the refractivity virial term, and the
    magnetic term.
  - The complex index of the noble gases in the extreme ultraviolet.
  - The second hyperpolarizability of each gas as a process-aware dispersion model,
    with the conversions from the conventions of the literature.
  - Every published formula of the curated sources as a named model, and a comparison
    with the refractiveindex.info database.
- **Materials:** the records name a material from a registry
  (`optmat.materials`, `data/materials.toml`) in place of a gas species, so that glasses
  and crystals share the core. Solid datasets carry no gas state, and an anisotropic
  crystal's datasets name their axis. The gases' dependency, pyfluids, is the `gas`
  extra (TC-318).
- **Solids (`optmat.solid`), fused silica first:**
  - Franta et al.'s universal dispersion model, implemented from its published
    equations. It reproduces the authors' tables to their last digit, with exact
    derivatives (TC-321).
  - The model interface for solids: complex index, d eps/d omega, group index, GVD.
  - Silica's recommended complex index from 25 nm to 125 um, `fit/SiO2`: Franta et
    al.'s model refitted to the absolute refractometry of the window (Leviton & Frey;
    Malitson, whose indices are relative to air), which Franta's published model misses
    by up to 3.4e-4 in n and 1.5 % in GVD (TC-399).
  - An index relative to room air (Edlen's dry air at any temperature and pressure),
    and Leviton & Frey's temperature-dependent formula of silica, 30-300 K (TC-398).
  - Silica's n2(lambda) with its Raman fraction and Hollenbeck's response (TC-393).
  - A fitting framework for solid dielectric functions (TC-322).
- **The solid-state literature:** 54 new source records (silica's linear index, n2,
  two-photon absorption, ionisation and damage, the Raman response, compilations),
  with data for some sixty solids, and the refractiveindex.info solid entries and
  Franta's newAD2 files as external data (TC-320, TC-325 to TC-392).
- **The documentation:** the manual, with accuracy pages checked against the stored
  fits by the tests, and a gallery of the data and the fits.
- **Independent numerical references:** recorded values check that optmat's named
  formulas, densities, rates and Raman parameters reproduce their reference
  implementations (`tests/test_luna_reference.py`).
- **The manual for publication:**
  - the gallery's figures are SVG (with PNG thumbnails), in sections;
  - the front pages and the quickstart cover the solids, the Raman response and
    ionisation, and a test keeps their counts of records current
    (`tests/test_counts.py`);
  - the snippets of the README and of every page that has them run in the tests;
  - a Read the Docs configuration.
- **No unlicensed third-party files** (TC-727): Babar & Weaver's workbook, McPeak et
  al.'s data files, the NIFS tables, Franta et al.'s newAD2 files and the SCHOTT data
  sheet's text are no longer in the package. They are attached to their papers in the
  curators' paper library, from which the extraction scripts read them. The four newAD2
  parameter sets that formulas use are held in the source records as the formulas'
  terms, and the silica fit starts from `franta_2016/udm_lithosil_q2`.
- `optmat.gas.linewidths`, `optmat.gas.molecular_ionisation_rates` and
  `optmat.gas.molecular_ionisation_validation` load as attributes of `optmat.gas`, and a
  test keeps each subpackage's list of modules complete.
- **The ionisation rates of the noble gases** (`optmat.gas.ionisation_rates`,
  `optmat.gas.ionisation_validation`; `docs/ionisation_rates.md`; TC-652 to TC-667):
  - the rate of He, Ne, Ar, Kr and Xe in a linearly polarised field, cycle averaged
    against the amplitude and instantaneous against the field, in SI, with the
    sources, range and uncertainty of each model;
  - the default model: the rate of Perelomov, Popov and Terent'ev with the static
    Coulomb factor, Tong and Lin's correction of the static rate, Tong, Zhao and Lin's
    coefficients, the fine structure of the ion, and an effective number of electrons
    set from all-electron static rates (1.47 for He, 1.0 for the others);
  - the forms of the literature as named models, for comparison (`ppt_larochelle`,
    `ppt_ilkov`, `ppt_popruzhenko`);
  - an instantaneous rate whose average over a cycle is the cycle-averaged rate (to
    3 %), where the usual division by the factor of the tunnelling limit gives up to
    1.4 times too much at 800 nm and 1.9 at 400 nm;
  - `jagau_2016`: the all-electron static rates of He, Ne and Ar from Jagau's
    supplement (Tables S3 to S5), on which the effective number of electrons of Ne
    and Ar now rests (TC-668, TC-669); `perelomov_1967` and `perelomov_1968`, the
    sources of the static Coulomb factor and of its range (TC-670);
  - `wallace_2016`: the 37 measured Kr+ yields of the supplement's Table I, and the
    models against them (`measured_points_test`; TC-672);
  - **optmat's rate of molecules** (`molecular_rate_model(gas, wavelength)`, the
    model `default`; TC-677 to TC-680): the rate of the noble gases carried to H2,
    D2, N2, O2, CO2, N2O and CH4, with the orientation-averaged structure factors of
    the molecular orbital, the static Coulomb factor, Tong & Lin's form of the
    correction at finite field fitted to all-electron static rates (N2, H2), the
    number of electrons in the ionising shell, the vertical ionisation energy, and
    for H2 and D2 the average over the vibrational ground state; within 3 % (H2) and
    13 % (N2) of the all-electron static rates it was fitted to;
  - CO in the molecular rates, and a test of N2's correction on it: with N2's alpha
    unfitted the default is 1.26 to 1.32 times CO's orientation-averaged
    all-electron static rates (thirteen angles, 0 to 180 degrees); the default is
    given to 1.5 times the barrier-suppression field, below which every molecule's
    rate rises (it stopped rising at 1.7 with alpha 14); Guo et al.'s molecule-to-atom
    comparison (N2/Ar agrees, O2/Xe 4.5 times below); three figures in the manual
    (TC-683);
  - `zhao_2011`: SF6's structure coefficients, read on the page by the maintainer's
    exception (the library's cells of the rotated table are unusable; the script
    checks the reading against the degenerate orbitals and the cells the library
    does hold), and SF6 in the default rate: 80 times N2's at the same intensity, a
    prediction that nothing tests (factor 10; TC-684);
  - `zhao_2011`: CH4's structure coefficients, complex and over three degenerate
    orbitals (`combined_structure_coefficients`; TC-681);
  - **the molecules' rates in the forms of the literature**
    (`optmat.gas.molecular_ionisation_rates`, `docs/molecular_ionisation_rates.md`;
    TC-673, TC-674): the molecular tunnelling rate of Tong, Zhao & Lin averaged over
    orientations, the rate with Talebpour et al.'s effective charges, and the
    hydrogen-like tunnelling rate, each of one electron and as its source gives it.
    They reproduce solprop's three molecular rates to rounding. None is validated and
    there is no default; against the orientation average of all-electron static rates
    (`optmat.gas.molecular_ionisation_validation`) the tunnelling form is 0.9 to 6.6
    times the reference for N2 and twice it for H2;
  - findings: the factor (1 + gamma^2)^(-n*) in the Coulomb term of Ilkov et al.'s
    and Bergé et al.'s rate is derived nowhere and makes the rate too low by 2 at a
    Keldysh parameter of 1; Ilkov et al.'s eq. (8) is misprinted by a factor of 2; the
    factor exp(-2 alpha_+ F) of Luna and solprop is in no paper;
  - validation against two-electron calculations of helium (static, 780 and 390 nm),
    exact one-electron static rates of argon, one-electron TDSE yields of He and Ar,
    and the measured yields of Wallace et al. (matched within their calibration) and
    Wahlstrand et al. (1.7 to 8 times below the model; open);
  - new records: `wallace_2016`, `travers_2026_trecx` (the maintainer's tRecX yields,
    tier 3), the yield fits of `wahlstrand_2018`, the fine-structure levels of the
    ions in `nist_asd_2026`; a new quantity, `level_energy`.
- **The ionisation parameters of the gases** (`optmat.gas.ionisation`, TC-440 to
  TC-444; `docs/ionisation.md`), for the propagation code's rates:
  - the ionisation energies (NIST ASD for the atoms; the primary determinations for
    the molecules, TC-542), adiabatic and vertical;
  - the ionising orbital and its structure coefficients (Tong et al. 2002, Zhao et
    al. 2010);
  - n*, l* and the ADK and PPT coefficient (Ammosov et al.'s eq. (20) corrected);
  - the static polarizabilities of the neutral (the fitted models' electronic static
    limit) and of the cation;
  - the model constants: Tong & Lin's and Zhang et al.'s barrier-suppression
    constants, and Talebpour et al.'s effective charges.

  Validated against the ab initio static-field rates of He, H2 and N2, and Zhang et
  al.'s single-active-electron rates of Ar. Gas-phase
  cations are materials of their own.
- **The Raman lines of the gases** (`optmat.gas.raman`, `optmat.gas.linewidths`,
  TC-442, TC-492, TC-541, TC-542; `docs/raman.md`) of N2, O2, H2, D2, CH4, SF6, CO2
  and N2O (vibrational only): every rotational S-branch and vibrational Q-branch line,
  with its frequency, amplitude (in solprop's convention) and width at a stated
  temperature, density, composition and pump wavelength.
  - The widths are per line, against J and temperature where the data allow, with the
    Dicke narrowing of H2 and D2.
  - The transition polarizabilities are ab initio per J for H2 and D2 (Hunt et al.),
    or come from absolute cross sections. They are dispersed to the pump with Bischel
    & Black's measured law, which follows Rychlewski's ab initio dispersion to 0.8 %.
  - The steady-state gains are validated against the measured gains of H2, D2 and
    CH4.
- **The rotational nonlinear index** of the lines (TC-493): it reproduces the
  published rotational n2 of N2, O2 and D2 to 0-4 %. That settles the factor of two of
  TC-312, which lies in a paper's printed response function.
- **The Raman response functions** (TC-685 to TC-696; `docs/raman_response.md`):
  - of the gases and dry air (`optmat.gas.raman_response`): h(t) and its analytic
    transform, with Hellwarth's tensor form (each line's isotropic and anisotropic
    share, the anisotropic part of the Q lines included), each line's width as
    c + d/n + beta n with its damping group, H2's line shift, the partners of a
    mixture, and the densities between which the lines are independent Lorentzians;
  - of fused silica (`optmat.solid.raman_response`): Hollenbeck & Cantrell's response,
    normalised analytically, with the Kerr fit's Raman fraction (the default), and Lin &
    Agrawal's and Blow & Wood's forms; the size of the kernel in a real-field model,
    K = 2 n0² eps0 c f_R n2, and the gain;
  - validated against measured gains, the decay of O2's revivals and silica's absolute
    gain and nuclear share (`optmat.gas.raman_response_validation`,
    `optmat.solid.raman_response_validation`; `stolen_1973` and `hellwarth_1975` curated),
    and compared with solprop's responses (`tools/solprop_raman_response.py`);
  - the Kerr fit's eta(tau) now uses the analytic transform of the response, and the
    fit takes Stolen & Ippen's measured Raman gain as data in place of the prior on
    f_R (TC-696): f_R = 0.22 +- 0.03 (was 0.15, set by the prior), the long-pulse n2
    2.7 % higher and its electronic part 5.4 % lower at 800 nm.
- **Metals, crystals and more of the literature** (TC-543 to TC-625): 67 source
  records.
  - **Metals** (a new material kind, TC-543): Babar & Weaver's and Werner et al.'s
    permittivities of Ag, Au and Cu, Ordal et al.'s aluminium (51 errata).
  - **Crystals:** the primary indices and formulas of CaF2 (Daimon & Masumura,
    Malitson), MgF2 (Dodge), BaF2 and sapphire (Malitson); silica against temperature
    (Matsuoka, Waxler & Cleek) and in the vacuum ultraviolet (Lamy, Rife & Osantowski,
    Tan et al.).
  - **Solids' nonlinear data:** n2 (Milam & Weber, Weber et al., Adair et al.),
    two-photon absorption (Liu et al.), breakdown (Du, Lenzner, Liu, Mero) and Blow &
    Wood's Raman response.
  - **Gases:** Brion's dipole (e,e) spectra of H2, N2, O2, CO, CO2, N2O and CH4;
    Meath's DOSDs of CH4, CO, CO2 and SF6; Hohm's polarizabilities against
    temperature; static dielectric virials; the CIPM-2007 air; ionisation energies,
    PPT and static-field rates; H2, D2 and O2 lines and widths.
- **More of the metals** (TC-626 to TC-648): 8 new source records and 6 extended.
  - The DESY report SR-74/7 (Hagemann, Gudat & Kunz): Al, Cu, Ag, Au and amorphous
    Al2O3 from the far infrared to x-rays, with 1 - n above 65-122 eV.
  - Weaver et al.'s Physics Data 18-2: Shiles et al.'s aluminium (0.04-1551 eV) and
    Weaver's gold. Babar & Weaver's analysis to 30 eV from Weaver's website.
  - Johnson & Christy's noble-metal films with their errors; Rakic's aluminium;
    Ordal et al.'s infrared compilation (1983) and copper (1985); Leveque et al.
  - McPeak et al.'s film data files, and Stahrenberg et al.'s numbers as
    refractiveindex.info gives them.
- **Sapphire's two rays and N-BK7** (TC-649 to TC-651): Malitson & Dodge's 1972
  Sellmeier pair from Weber's *Handbook of Optical Materials* (2003) and SCHOTT's N-BK7
  data sheet (2023, with TIE-29's conventions), each made absolute; the data sheet's
  internal transmittance (a new quantity, `internal_transmittance`) and dn/dT curated.
- **The window crystals' recommended models** (TC-545, TC-546; `docs/crystals.md`):
  each primary Sellmeier formula made absolute with the gas the index was measured
  against (nitrogen, `index_medium = "nitrogen"`, or room air). A uniaxial crystal's
  recommendation is per axis.
- **The metals** (TC-544, TC-619, TC-646, TC-647; `docs/metals.md`): table models,
  curated datasets joined by range with a reason for each and a monotone interpolant;
  silver, gold and copper from 3-4 nm to 12.4 µm (copper to 56 µm), aluminium from
  0.8 nm to 200 µm. A segment may be a table of 1 - n
  (`refractive_index_decrement`) or Henke et al.'s scattering factors (TC-644, TC-645).
- **The composition of dry air** (`optmat.gas.composition`, TC-548): CIPM-2007.
- **The new sources in the models:** Milam & Weber's silica n2 in the Kerr fit (+0.9 %,
  TC-594); new validators of six gas fits, and five DOSDs' static limits within 0.2 %
  of the fits (TC-606); the ionisation and Raman choices on the originals (TC-547,
  TC-586).
- **Ionisation, Raman and solid-state sources** (TC-494 to TC-540): 36 source
  records.
  - Primary ionisation energies of the molecules. Measured yields, which are figures
    only and calibrated against models.
  - Zhang et al.'s barrier-suppression constants and static rates, Talebpour et al.'s
    effective charges, and the Zhao et al. erratum.
  - The Weber volume (1979): Schrötter & Klöckner's normalised Raman cross sections,
    and Stoicheff's H2 and D2 and Wang's CO2 rotational constants.
  - Hunt et al.'s transition polarizabilities of H2 and D2.
  - Bischel & Black's dispersion of the cross sections.
  - The N2 Q widths of Rosasco et al. and Rahn & Palmer, and the D2 widths of Smyth et
    al. and Keijser et al. (TC-541).
  - Two-photon absorption and n2 of solids: Sheik-Bahae et al. 1991 and 1994,
    Bristow et al., Patankar et al.
- **The ionisation and Raman literature** (TC-445 to TC-491): 43
  source records:
  - ionisation energies, structure coefficients, ion polarizabilities, static-field
    rates;
  - rotational and vibrational constants, Raman cross sections, gains and
    derivatives;
  - linewidths, shifts and their temperature laws.
- **Raman and ionisation for solprop** (TC-439): the scope and the stages of the work
  that lets solprop carry no material parameters. Every Raman and ionisation value
  solprop uses is recorded as a parity target (`tools/solprop_raman_ionisation.py`,
  `tests/data/solprop_{raman,ionisation}.json`), and the missing papers are listed on
  the roadmap.
- **Further gas sources** (TC-400 to TC-438): 31 source records.
  - Hohm's 325-633 nm dispersion in the Kr, Xe, O2, N2O, CO2 and CH4 fits: its shape at
    4e-5, its level free (TC-400).
  - Static dielectric polarizabilities and 633 nm refractivities as validators.
  - The absolute dc Kerr of helium and H2 (Tammer et al. 1992), CH4 and CO2 as
    validators of the hyperpolarizability models (TC-401).
  - Molecular photoabsorption (Brion & Thomson 1984) and dipole sums (Olney et al.
    1997), for the molecules' extreme ultraviolet.
