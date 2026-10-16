# tRecX two-photon ionisation and excitation of He and Ar (J. C. Travers, October 2026)

The two-photon TDSE campaign (register entry TC-1048; optmat's specification of the runs,
"the brief" in the files here, is unpublished, as are its two addenda of 9 October 2026),
computed by J. C. Travers with tRecX (A. Scrinzi,
Comp. Phys. Comm. 270, 108146 (2022), https://doi.org/10.1016/j.cpc.2021.108146).
Unpublished. The source record is `travers_2026_trecx_two_photon`; its three
cross-section files are written from `run1/results_hacc.csv`, `run1/results_sae.csv` and
`run3/results.csv` by `tools/extract/travers_2026_trecx_two_photon.py`, and its two
amplitudes are inline values (run 4, below). Tier 3: validators only, ranked below
published ab initio calculations and measurements; nothing is fitted to them.

Copied on 2026-10-11 from the maintainer's working directory, a local copy of tRecX's
`testruns/two_photon/` directory on a workstation (the runs' own directory; the `.git`
of `testruns/` does not exist, so the files are identified by their sha256 below).

Paths: `$RUN_ROOT` stands for the run directory on the cluster (the campaign's
`tdse_two_photon/` directory, where every run and pulse file lived) and `$TRECX_DATA`
for tRecX's data directory there (which holds `chemical/He19_Standard`). The `.inp`
inputs under `inputs/` are the archived originals, byte for byte, and keep their run
paths, which refer to the cluster run directory.

## The code

- tRecX `a987d7c790e52b6d7f44a5615a665e0abc3e7f4a` (28 September 2026 16:58 +0100,
  "Occupation: add bi-orthogonal populations OccDual{H0:k} = |D_k^T Psi|^2"), the commit
  of the 42 fs argon runs, built with GCC 16 on the cluster. Runs 1-4 and 6 used it as committed.
- **Run 5** used the same commit plus an uncommitted patch that writes the energies of the
  Occupation states to `<output>/occupation_energies` (state, Re E, Im E; master rank):
  `occE_patch_20261010.diff` on the cluster, sha256 beginning `bd3c33ee` (the runs'
  `NOTES.md`; an earlier version, `3aa4093a...`, printed nothing below verbosity "full"
  and was not used), built as job 479258 (`tRecX-occdual-occE-gcc16{,-mpi}`;
  binaries' sha256 beginning `03eabe9b` (scalar) and `6643177d` (MPI)). haCC has no
  eigenvalue calculation in this version (`run_trecx.cpp` aborts when `Eigen: select` is
  set), so the patch is the only source of haCC's level energies. The patch itself is
  not archived here; only the prefixes of its hash are recorded (in `NOTES.md`).
- haCC's helium: the Scrinzi group's `He19_Standard` data set, at
  `$TRECX_DATA/chemical/He19_Standard` on the cluster (the inputs' `Chemical:` line);
  the commit of the `chemical` repository is not recorded by the runs.

## The runs

All on a Slurm cluster, 9-11 October 2026, in `$RUN_ROOT/`; every
run exited 0 (`REPORT.md`, `STATUS.md`). Jobs: 478310 and 478367 (run 1, one electron,
44 inputs), 478347 and 478348 (run 1, haCC, 17 runs), 479251 (runs 2, 3, 4, 6, one
electron, 110 inputs), 479252 and 479253 (run 2, haCC, 16 runs), 479260 and 479261 (run 5,
haCC, 15 runs, the patched build). Superseded: 478286-478288, 478302, 478303 (run 1's
first version, whose spectral range, tRecX's default max(12 Up, 3 Omega), came out near
3e-4 au, far below the two-photon peaks; kept by the runs as `run1_v1_spec_kmax0p02/`, not
archived here).

| Run | What | Model | Pulses | Table |
|---|---|---|---|---|
| 1 | He, degenerate sigma2, 13/14/16/18/20 eV, 1e11/1e12/1e13 W/cm2; convergence at 16 eV, 1e12 | one electron (BOX 40; BOX 80, 120 at 1e12) and haCC | 20 fs Gaussian | `run1/results_sae.csv`, `run1/results_hacc.csv` |
| 2 | He two-colour kernel, Omega/2 = 14, 17, 20 eV, nu 0.5-6 eV, 1e12 per colour, and one-colour references | one electron (BOX 80) and haCC | 20 fs | `run2/results.csv` |
| 3 | Ar degenerate 8-11 eV (1e11/1e12/1e13) and two-colour Omega/2 = 8.5, 9.5, 10.5 eV; m = 0 and 1; convergence at 9 eV | one electron (BOX 80) | 20 fs | `run3/results.csv` |
| 4 | Ar 3p -> 4p at 6.463 eV, 1e10 and 1e11, 30 and 50 fs, m = 0 and 1 | one electron | 30, 50 fs | `run4/results.csv` |
| 5 | He 1s2s at 120.3 nm (10.306 eV), 1, 3, 10 fs and 3 fs +-10 fs2, 1e12/1e13/1e14 | haCC | | `run5/results.csv` |
| 6 | He broadband at 14 eV, 1, 2, 5 fs and 1 fs +-5 fs2, 1e12/1e13 | one electron (BOX 80) | | `run6/results.csv` |

The conventions (field on the vector potential, velocity gauge, linear polarisation along
z, g = exp(-2 ln2 t^2/T^2) with T the full width at half maximum of the intensity, CEP 0
per colour, I0 = (1/2) eps0 c E0^2 = E0^2 x 3.50944758e16 W/cm2 per colour, chirp as a
quadratic spectral phase on A with phi2 > 0 an up-chirp and I0 the transform-limited peak,
window [-3T_eff, +3T_eff + 1000 au]), the models' own energies and the discretisation are
in `REPORT.md`; every run's colours, width, chirp, peak A and window in
`manifest_runs2to6.csv` (runs 2-6) and in the inputs (run 1). No NIST value was used in
the runs.

## Checked (2026-10-11)

- `tools/extract/travers_2026_trecx_two_photon.py` reproduces every sigma2 the runs'
  summaries print, to the printed digits: haCC 11.890, 12.804, 10.723, 9.063, 10.175 and
  one-electron 7.456, 7.328, 6.524, 6.228, 9.346 (1e-52 cm4 s; `run1/summary_run1.md`,
  `REPORT.md`), argon 1.921, 1.916, 2.094, 2.503, 3.052 (1e-50 cm4 s;
  `summary_runs2to6.md`). Its int Phi^2 dt reproduces the table's own
  `sigma2_peak2_cm4s` from `P_peak2` to 1e-12 (haCC).
- `run4/pt_two_photon.py.txt`, run by the curator (Python 3, numpy) on 2026-10-11, prints
  M_0 = +8.914061e+00 and M_1 = +2.418379e+00 au at dr = 0.002 (+8.902447e+00 and
  +2.416898e+00 at 0.004), as `extract_runs.py.txt` holds them (`M_PT`) and `NOTES.md`
  states; `run4/results.csv` gives TDSE/PT 0.9982 to 0.9996 at 1e10 W/cm2.
- The input decks here are the local copies written by `make_inputs.py` and
  `make_runs.py` and uploaded to the cluster; they were not compared with tRecX's own copies (`<name>/0000/inpc` on the cluster). The job lists
  `run*/list_*.txt` name 198 of the 203 decks; the five `run1/sae/eigrun_*` decks (the
  field-free spectra whose energies label the bound states) ran in jobs 478310 and 478367
  from lists not copied. `run1/hacc/eig_l10_r40.inp` is in `run1/list_hacc.txt` but was
  not run: haCC cannot compute eigenvalues in this tRecX, and the haCC jobs ran from a
  list without it (`list_hacc_noeig`, `NOTES.md`).
- The manifest of runs 2-6 has 141 rows, one per (model, name), no duplicate.

## Caveats

- The one-electron norm-loss and projection carry a no-field drift of the discretisation
  (He 1.43e-7 over run 1's 5961 au at BOX 40, 1.9e-8 with step x 2/3 or cutEnergy 150,
  3.0e-8 over run 6's 1 fs window; Ar 2.4e-9): the tables' `P_proj_corr` (runs 2-6) and
  the record's one-electron helium sigma2 subtract the drift of the no-field run of the
  same discretisation; at 1e11 W/cm2 (He) and for run 6's chirped pulses (whose windows are
  3.6 times longer than the no-field run's) only the spectrum is usable.
- The one-electron helium spectra at BOX 40 differ from the projection by -7.4 to +9.0 %
  (energy-dependent; `NOTES.md` gives -11 to +9 % from an earlier reading); at BOX 80 and 120 by at most 1.9 %. The argon spectrum and the
  drift-corrected projection differ by +10.9 % at 8 eV (0.23 eV photoelectrons), -1.2 to
  -3.7 % at 8.5-9.5 eV and +0.5 to +1.1 % at 10-10.5 eV; the runs' summary says
  "elsewhere <= 2 %", which 9 and 9.5 eV exceed.
- Run 4's signs are from lowest-order perturbation theory in the same potential, not from
  a projection of the TDSE wave function (tRecX writes psi(T) only in its binary format,
  with an arbitrary eigensolver phase); the TDSE confirms the moduli only.
- Run 5's haCC state 3 (-0.07776 au from He+ 1s) is never populated and matches no helium
  singlet: unidentified (possibly an artefact of the basis).

## Not archived

- The pulse files: A(t)/max|A| on a 0.01 au grid, written by `make_inputs.py pulse` (run 1)
  and `make_runs.py` (`pulse()`, runs 2-6); 1.2 GB on the cluster in
  `$RUN_ROOT/pulse/`. Run 1's five are also on a workstation in the local copy of
  tRecX's `testruns/two_photon/pulse/` (19,125,957 bytes each), sha256:
  `988144d607573841eca68b27023d9f5273e30bb35574ad423553218d2e969c80` (13 eV),
  `ee655ba5f4526f000e900d1526c9579f499c96115b78e31569b677b411f0c292` (14 eV),
  `3fab19c25d8f4ed1c4f6758a9aa037f043797058058b4b4e8ea2d65a2d19c470` (16 eV),
  `c683cf0ee2d5b1a27b3d0d3d853e9a2be6bdf568f89b69dcb0dc112a77becf03` (18 eV),
  `833b31363081244f7cb4d34e46c71ef5094a8ab6173f897927cdfdfeaa09300e` (20 eV).
  `pulse_check.txt` here gives their cycles, peak field and bandwidth.
- The run directories' outputs (`expec`, `spec_total`/`spec_channels`, `eig`,
  `occupation_energies`, `Laser`, `WARNINGS`, `inpc`): on the cluster only.
- Run 5's P(1s2s) against time: `$RUN_ROOT/run5/hacc/<name>/pop_vs_t.csv` on the cluster.
- Run 5b's populations every 1 au (11 files, 3.5 MB): `run5b/pop_<name>.csv` beside the
  runs on the workstation and on the cluster.
- Run 5c's populations and dipoles every 1 au: `run5c/pop_<name>.csv`, likewise.
- The run-5 patch (above) and the first version of run 1.

## Files

| File | What it is |
|---|---|
| `REPORT.md` | the runs' report for optmat (11 Oct 2026, updated with addendum 3's run 5 checks, 11 Oct, and addendum 4's 1s2p, 12 Oct): conventions, the models' energies, discretisation, convergence, results in brief, deviations from the brief |
| `STATUS.md` | the campaign's job table and the decisions taken with optmat (9 Oct) |
| `NOTES.md` | the working log: every job, fix and check |
| `run1/summary_run1.md` | run 1's summary table and convergence |
| `summary_runs2to6.md` | runs 2-6's summary tables |
| `run1/results_sae.csv`, `run1/results_hacc.csv` | run 1, one row per run (written by `extract_run1.py`) |
| `run2/results.csv` ... `run6/results.csv` | runs 2-6, one row per (model, run) (written by `extract_runs.py`) |
| `manifest_runs2to6.csv` | every run of 2-6: path on the cluster, window, width, chirp, colours, peak A, peak field |
| `pulse_check.txt` | run 1's pulse files: cycles, max|E|/E0, bandwidth, window |
| `make_inputs.py.txt` | writes run 1's pulse files and inputs (the runs' `make_inputs.py`, as text) |
| `make_runs.py.txt` | writes runs 2-6's pulse files, inputs and manifest |
| `extract_run1.py.txt` | reads run 1's outputs into `run1/results_*.csv` |
| `extract_runs.py.txt` | reads runs 2-6's outputs into `run*/results.csv` (with run 4's M_PT) |
| `run4/pt_two_photon.py.txt` | run 4's lowest-order amplitudes by Dalgarno-Lewis on a radial grid |
| `run5b/results.csv` | the run 5 checks of the campaign's addendum 3 (11 Oct 2026, jobs 479500/479501): every run 5 and run 5b haCC run, with the post-pulse mean P(1s2s), the fitted K and resonance and the ratio to lowest order (written by `extract_run5b.py`) |
| `run5b/scan_fit.csv` | the 10 fs carrier scan's fit of K and the resonance (2 omega = 20.7386 eV) |
| `run5b/list_hacc.txt` | run 5b's job list |
| `extract_run5b.py.txt` | reads the run 5 and 5b outputs into `run5b/results.csv` |
| `run5c/results.csv` | the campaign's addendum 4 (11-12 Oct 2026, jobs 479626/479627): haCC's propagated 1s2p from the dipole beat and a five-carrier scan at 1e9 W/cm2, 10 fs (written by `extract_run5c.py`) |
| `run5c/scan_fit.csv` | the scan's fit of K and the resonance (21.20892 eV) |
| `run5c/list_hacc.txt` | run 5c's job list |
| `extract_run5c.py.txt` | reads the run 5c outputs into `run5c/results.csv` |
| `run*/list_*.txt` | the job lists (input paths on the cluster, under `$RUN_ROOT/`) |
| `inputs/run<N>/<model>/<name>.inp` | the 203 input decks as generated (246 kB in all) |

```
e0beac56f64008c91bf898d26f2f54382847a7da9d2fd8266fa182a87a16a746  NOTES.md
c4262f8d3b448a87059c970b98a60d573534eb8714cd981e7c725ce3342edaa6  REPORT.md
ca711a53308628ca45e8c9d91984555f38f2b1270ec4a200fc73364903b7d75c  STATUS.md
4ff8113052477a00fab6d5eefb38792847efba0a80f982efcc5a4f361a994d60  extract_run1.py.txt
3184ec70232b0781ef6b7612fc5277974f6c20cb92f71370561b8ddd1095a022  extract_runs.py.txt
db02745ac575d857b9b0f3c63ead6dbdccc65c8693e9c157e87f0918c7532647  make_inputs.py.txt
486ea4dc4a726a3ab161b17800635e944ab3d34bc977631b22b828ec91a5d99a  make_runs.py.txt
0c01835d7fd98415442096655fc22b9dab36221018260efcf8a12a10c34ef99c  manifest_runs2to6.csv
813c5448c7bb414a4180fe06f6bc0a22f6a0c60f1c52f19fbf275812caf0cd1e  pulse_check.txt
6cd20718339b0d7cc24e6b127ce2085dbdbbbb7a343ea98b8e8c6988467243db  run1/list_hacc.txt
f4064cfd294f89c8570e37cd2c02f45192a9373f6d95ad384742e469425979ce  run1/list_sae.txt
8d059269d66521d538732c48e59f8374f0119a6600404d8a3b6c8d514dd939b5  run1/results_hacc.csv
ab6c7acb2cc949893dac3397780b086d1d14b3a7d1f70f28eb105e57669d0d38  run1/results_sae.csv
7991061c681aa86a00b45096032ef8604d08fe1182b112f1b828ad0f9aae595d  run1/summary_run1.md
91d0681c1f2f863158971e3d7a14992bbd76ceb12e58567fc649c8d6805b08f2  run2/list_hacc.txt
02e1639c6df07115970cb9176bb2367eec1cf69ccb20ff21a0a793fd0ae15504  run2/list_sae.txt
9869f1e6b8b9374f1c3db4754fd6ebdc9b29dca8228472f4bcedea8ff56691e5  run2/results.csv
5b97099d8913efc117aee1533331720248056089b207535a5277b560dd5a7cb1  run3/list_sae.txt
c36272496fd7b354fe7b375963292665f65197dee500a9267cc72be9ffdd8b75  run3/results.csv
cf4ca16641494f873ba0782c77ea2ab0f631e81ea7671df7be0ec0b9c1782033  run4/list_sae.txt
dfe2d09f0e5b54d41e10b26ba7c85d69d6b5f4d261f27b13dcf230f68fc45102  run4/pt_two_photon.py.txt
12548c4ac272461a5c2116f978942d9c9e5c74ceeca6634df8497e01f1a85be9  run4/results.csv
ddda779c5d0242a65ab9b178597d8d58c242b9298f3ffc37f1d995f0d221013c  run5/list_hacc.txt
646aab22e18981f6bbe0c08a4920447512df9afc394d26b2e4485e19def1e850  run5/results.csv
5c447fcaa8733300667c9c59f25f715cca8b51515b6edfe34256d392aa98ac84  run6/list_sae.txt
f3198cc2b7bee969e5cdeec6a5ee6ee14437b30f35f7f24b17848cc204d17813  run6/results.csv
17276df2904c2e350aef3d4f6ee8d6f39a8347c9aa63015173112d3f5674d3b4  summary_runs2to6.md
188ce7ac7d46792ec8243c3150cd9e46553928d5cafa5493835a59b593ecc3d6  extract_run5b.py.txt
bb1ae5cf90d3678724286b490a58df8ea3a524355d11f2a669a126da276e82a6  run5b/list_hacc.txt
3ccb2ec49c5944978bf0d8406d34ceed56bf1d7897c892886666887815e3e84d  run5b/results.csv
927b0532f991467a3b520b4780479d6ea36ae5447517bef8baeb66742bb9d9e9  run5b/scan_fit.csv
b7831fb6979e68be8497c5a70c7b3f019c2a2259a30a2cd2a2bafe88f6094310  extract_run5c.py.txt
98cf76000dff6a14cc392d2a08d95dafc1ed6413bab590d8fde48602ce81b1eb  run5c/list_hacc.txt
c604f7681e1e9cc47289b3fb5a14ac6a71335e39129998786983e5fd92d8d586  run5c/results.csv
ea64b98bf71eec5c9fe73896d659a3f3f0c0f20cfda916be7cc5f009d221e187  run5c/scan_fit.csv
```

The files above were edited for publication on 2026-10-16, paths and wording only (the
cluster's paths to `$RUN_ROOT` and `$TRECX_DATA`, host and machine names, and the name of
the unpublished brief; no number changed), and the hashes above are of the edited files.
The scripts therefore no longer write the archived inputs byte for byte (their comment
lines and paths differ), and the script hashes quoted in `NOTES.md` are of the files as
run. The sha256 of the files as copied on 2026-10-11, before those edits:

```
0e762ba3a43c30f9815c90e54e69d86d4e44bef4bde087d016c327128f78c3d7  NOTES.md
bd1a395fd8b4decd550168c96b28d53fda4e96ec9c1bac321b900a7e9a212a7a  REPORT.md
83616ac5c78d59c84d6a26a5de9685ea6d594e98d40c43a1a8c2a488d5a73676  STATUS.md
1e5100db6f9d50052ed0194a2c4f6952004150aef5a85313aabd406e6cb79d1e  extract_run1.py.txt
9b4f26038788a736229dd339aceec77a946b3ba87e877c5c651aa68506dc592c  extract_runs.py.txt
6d3f8e21922207cff1810661075c7233426ccd7cf0c90d27f7f392efc797327e  make_inputs.py.txt
bac80f54bdfdab7aa5235f370592cb129012ee00a6fa6a14e643eb32c58c55f5  make_runs.py.txt
23459eeb9a8de8997f956c907878ba72f84a0860c1a590924e8d403675756755  manifest_runs2to6.csv
d44cf470779b41090dc69b5fe0e9ec20caae277e86cb93a9c54c9664c89cfdf1  run1/list_hacc.txt
d76060309cab2e6f652e3d6b6188ba2cd0f8cbf3e422fdf45b11155fedec0404  run1/list_sae.txt
3d63fa6d3d8a3a4099f9f6272f39b8aff1a81bb23e4df71178aef6856723adca  run2/list_hacc.txt
5643f7d47f236b71602c9b895ca4a9d27e90e56c43f6140f571dccc597eaca11  run2/list_sae.txt
d94950de0e65652b83bbeb7c5c641e5c7c3cac4794a861fc95a2ebe7d29abf47  run3/list_sae.txt
dea6dab82f4b7c972c38c677ad7a53c2e76306ffcfdd3d048544b0ba1fca8f0e  run4/list_sae.txt
ad7b7a559f7fe13351a7c3bace40d1fc248437e9050855a987c1c36a7176152a  run5/list_hacc.txt
04402045373de4bc238ad3d68beca9eb988be78aa57155f03feff0e6f8976759  run6/list_sae.txt
6171de1cc5fc1795bb4839a722dfcb560120943e9e126165ccdf908c2549fa3f  summary_runs2to6.md
b8b94831ae4ebf6a8f9f0c14240123fd59e713d165a219d2e3865a833b717df2  extract_run5b.py.txt
3222da15e0cd4b406bfea0cb58e0960f13cefc5cea2a32103bea4d240c32b96e  run5b/list_hacc.txt
3f24eabad044dbe47d6cc85229628023763176b8ae489f1962170a429b1df393  extract_run5c.py.txt
29b78e919cd343d94413b1c3c5d00fabcfe6f5f5a3b8dd0a641a394cf1a10194  run5c/list_hacc.txt
```

## sha256 of the inputs

```
226254cc8058995e94252f64850b2845d162a80704650709c2df509493ca1187  inputs/run1/hacc/eig_l10_r40.inp
5b358787c935bef89588e9e66b4b9e062bfa0182b76a41e224b125ab20d86384  inputs/run1/hacc/he_13eV_I1e11.inp
4a02bfbe98ee65c76ac772f5f082bb4ef2bcc92c29f0e8023456e435f10318ef  inputs/run1/hacc/he_13eV_I1e12.inp
cb86bbbbfeff20ac8c22dd3f4f2c4624c538184a9bebc85d4496a58c8951ce10  inputs/run1/hacc/he_13eV_I1e13.inp
f7975962d2c48e29619a559dbcbc58723e9b317a6abff72800ddf2a52537b253  inputs/run1/hacc/he_14eV_I1e11.inp
ee376d4b55ff2c65ccfb36e7870ec9b5372cad411b702f79a0c688d550381877  inputs/run1/hacc/he_14eV_I1e12.inp
e1b81bb14770f96abf4fd16b6b9482ecec01dda1f33ceaa6b554e130d5740640  inputs/run1/hacc/he_14eV_I1e13.inp
5dcb2d661a05ed80c471ffe5675d68f7f3890606b6d4bf5917e6b822c61b037a  inputs/run1/hacc/he_16eV_I1e11.inp
9d9a95a8af06132003ec57315a8370201d0965ba174499eef4e47d51bc24bd83  inputs/run1/hacc/he_16eV_I1e12.inp
ab0cf283edc2eb894690b66bd48bebdd135d3ce06aaa94a4ec15fb87d05ebd30  inputs/run1/hacc/he_16eV_I1e12_l15.inp
c7aef8398a8e4c078c468ab65f4ba056ab7b7a31662a9aaf66e6d8e785807a51  inputs/run1/hacc/he_16eV_I1e12_r60.inp
1436eaba1fb3f3afbce05a9d536a2e2df81a3dba71ce8fddce9ae2b99fc11957  inputs/run1/hacc/he_16eV_I1e13.inp
d96de53c38a3350d4ae288e3b56bcd698bfa6705f5c584bd37feac2823d10045  inputs/run1/hacc/he_18eV_I1e11.inp
643cb1b6db250eef22b1d44c925498f60a2a6777332a0e62233b34c8e5a919ae  inputs/run1/hacc/he_18eV_I1e12.inp
e1407c1e7f50bdd621799cbde4147622bd86a80992528eb12ab49ea37199ad2f  inputs/run1/hacc/he_18eV_I1e13.inp
d6644bf52865ba455b75a29e40491300eeac0c21faf81f54e432b9e36e484cef  inputs/run1/hacc/he_20eV_I1e11.inp
2fa2d7db227a990c32508a117ef0cfa341a3e52d1130b71c370d93acee044800  inputs/run1/hacc/he_20eV_I1e12.inp
3413ef8bae31eecae4ad8a2210aff88dbdbefc778a4ac471870f20454c7a5f77  inputs/run1/hacc/he_20eV_I1e13.inp
611209ec8dda5ca49091c7f1e15539be861c0f466d5498cdb62bf1b6a6bef921  inputs/run1/sae/eig_l1.inp
768ee5a5c5848e0bfe219fcd4fcad1508fb8eac86528f9469e1139543534a4e7  inputs/run1/sae/eig_l2.inp
1825c45cf71ce108a6ecd7b63ec876a6806dbb7350c0e1c35f66997974013d56  inputs/run1/sae/eig_l3.inp
378c545da10a5caa5c04d2b520e643e05de4088908351080963c5207b0f38e8e  inputs/run1/sae/eigrun_base.inp
a881000fc033891e5854e6f2f9f72b4fc958d5a4ab38f311c0877265cce48775  inputs/run1/sae/eigrun_box120.inp
2461d6a293f547ba638a40f2083b76e9e419b800c62d3fbe215f85579f598e35  inputs/run1/sae/eigrun_box80.inp
d2f692ef7a3dcddd9b1f681bd70d9340b87bf68d58e4353b98117be3c9296150  inputs/run1/sae/eigrun_eta15.inp
b13f6c62f176081a37ac90168ddc81291d537a3a1f95e69bd6631d3d6334b3a8  inputs/run1/sae/eigrun_th23.inp
42efdbedf3c2ca478bc26e095fb2128d7963189da002fcae42815d0f2e0edb4b  inputs/run1/sae/he_13eV_I1e11.inp
910451b4f0a077eef8bb6e9716d41a4c05182feba40f22de4f840fd9934e951f  inputs/run1/sae/he_13eV_I1e12.inp
47657a9f249fcb191444e4be010f0199d92d045dc41db2f0c00b2b3311652181  inputs/run1/sae/he_13eV_I1e12_box120.inp
ff7b7ed5709e2abcd8791d6f0469e7e59d9229e0f96f361552019ab6337f8ea5  inputs/run1/sae/he_13eV_I1e12_box80.inp
680af4d61e6d77f9784a770035d529e32a3a642493d3552d5af0a45881dec6c8  inputs/run1/sae/he_13eV_I1e13.inp
3365e60f367c7f6da0792046dfa3f89faa737a2c2e8baead1702111aed475ea9  inputs/run1/sae/he_14eV_I1e11.inp
bfdaffe31c8d8fca612aa313f30ad3d3624f69b1df70d967c497784d4c10b8b4  inputs/run1/sae/he_14eV_I1e12.inp
7177f645c6638876bca0ec77de6b3000659a5e86fc5612d519ebb904dfd8be3c  inputs/run1/sae/he_14eV_I1e12_box120.inp
6852fcc14dfce8a7713f52a3c4c3a1a78d748b4980f32c309dd4e3f78cc393e3  inputs/run1/sae/he_14eV_I1e12_box80.inp
3183f1d83289efd8fb1a1585d5f31181df4a3a0ac46f99050d142db69ac69db6  inputs/run1/sae/he_14eV_I1e13.inp
c283f8a3eb358c0fa3d975d2c8cc889b021616b42d86a8b49eef3ee7feabc789  inputs/run1/sae/he_16eV_I1e11.inp
71060c4e3b40d75971d71e27ece4edbebd3f420a3e091bb4945309973a2f25cd  inputs/run1/sae/he_16eV_I1e12.inp
c92cf238a7a5a68d48527f4ceb49fa4552b6c89c5d334d672a22551bd5b57a12  inputs/run1/sae/he_16eV_I1e12_box120.inp
d260a08c066c93d6b50b2f2e7294883e5109249d0d23d4ee239846b6eb85fb3c  inputs/run1/sae/he_16eV_I1e12_box80.inp
1b58d7853d1e78aad1aca82f64c3cd628856c28b9b511c00a01e38e87f0c5632  inputs/run1/sae/he_16eV_I1e12_cut150.inp
60d4284fe60918420e64c638266083ecbfa7e0b50763552f2f485387707e0d0b  inputs/run1/sae/he_16eV_I1e12_dt23.inp
f457153ac0d89dad2ffcc9b04e0de4c0c7337efaea986c84b92ec209d237c3a1  inputs/run1/sae/he_16eV_I1e12_eta15.inp
aa6ea7c5d7c8cbf05e9adecbb3bef571a858780629edf22b4e40256b8840ce32  inputs/run1/sae/he_16eV_I1e12_th23.inp
88e700a3edae48c0658ab05fde5d088de853b7f3e0e35b3411aedad0158c4787  inputs/run1/sae/he_16eV_I1e13.inp
6755004fb500b83fcca026dd0b0eb77d0eb58ffc5c09654a84baf1379243aa9a  inputs/run1/sae/he_18eV_I1e11.inp
5c4fa925e8b0bc6ebb39a8be8adf57a17bea77d697b76098c029d4a5f0c2d176  inputs/run1/sae/he_18eV_I1e12.inp
c733e2b276f6f7b3e45f053fa95de0abfe007e56c1a9f7bd9cbf01d7351d62fa  inputs/run1/sae/he_18eV_I1e12_box120.inp
8bb7dafa34f4968ab92f93949fab0f1257f0d80cab71c6144c14515c93a37e44  inputs/run1/sae/he_18eV_I1e12_box80.inp
ae2804a2f5016db79cb57cf230f518095ff4ef3c285390bf11bceddaaae7ff54  inputs/run1/sae/he_18eV_I1e13.inp
dca8226ddde6bb8ac388a84583122a23523a206c8e893413b5d68af901046e7b  inputs/run1/sae/he_20eV_I1e11.inp
7ae2045092deafea0305f86dc9fb158db15a71e8a15b3f2a15a0e5e749862931  inputs/run1/sae/he_20eV_I1e12.inp
5632e2c3f1d8fee5c68cad5e18d0064c292f1758cc8d0d7683cf58ee9d11ba64  inputs/run1/sae/he_20eV_I1e12_box120.inp
17fc26fe5069fd354a1021c65fedfc0472b4623789b8c5849a75a606e3c0241f  inputs/run1/sae/he_20eV_I1e12_box80.inp
f0c2b1effae3ecdfd4ab69022545741a78ba38aeed460f3d6a74e50ef81dc4f0  inputs/run1/sae/he_20eV_I1e13.inp
a51c0ba6b91e669584a274374bc3dc283f09d8c01c0bf857a016cbb1390b24f1  inputs/run1/sae/he_nofield.inp
e4f26388c8406a8583a4d66fbf2f8eb6cf800f9bb87b5334bcbad3b552fd3ee4  inputs/run1/sae/he_nofield_box120.inp
bb1e6a8210ac8d87671882c1acc687f8fb0fadb3d1ac5163fc21823e4f88ff5f  inputs/run1/sae/he_nofield_box80.inp
03b60c939de6e924c8bec06a627e5eb21f5d60f8b3355bec59712554eff17f6f  inputs/run1/sae/he_nofield_cut150.inp
e2d170e6e9ad2138dc012aefc0658bada2c19bb235887bf38f1b1ac88e7faf0d  inputs/run1/sae/he_nofield_dt23.inp
4fa10ec2ada44ef643fde93e3cdc683f83e30d568049b96a202a5f95bf167e28  inputs/run1/sae/he_nofield_eta15.inp
f649981417a9bad489cd5433d077d4f6afc21abb1bdcb7465a21d565517f4ec1  inputs/run1/sae/he_nofield_th23.inp
a217ea19a9de270cbd36e9c1a4401da56cec00ad7498399ad638e5159d82e709  inputs/run2/hacc/he_r2_14_nu0.inp
bd6e2d1e62cfc11441e8f5911301ee0ee44a5e758d228c9bfef4de8c54a2a2a7  inputs/run2/hacc/he_r2_14_nu0p5.inp
fdcf0058307ccb02d30423f23fca2e2f559eb25af0dd308a6058a0a935e5bd3c  inputs/run2/hacc/he_r2_14_nu1.inp
e5e6c9aeeb428dda91e00a19313b6464660954e0682b766d2831308dda6a2356  inputs/run2/hacc/he_r2_14_nu2.inp
b8b7ef202e20d672b08b163d27f17d2fe6af2418d5b0ae5376baca1eae226106  inputs/run2/hacc/he_r2_14_nu4.inp
7e6a498f159f2b9e7e0a67c44c6469ef31e6ba23f174714cb4ebc917160c3647  inputs/run2/hacc/he_r2_14_nu6.inp
c512bd5c89f8c6dbd366ed17c91234915d5019ee2f1c178137aec5a27f9b0b9a  inputs/run2/hacc/he_r2_17_nu0.inp
394764f74d3f26e2ea5082df47a566251d4418a4cba84d802d3fb90ca7702964  inputs/run2/hacc/he_r2_17_nu0p5.inp
f15a8a7415488237956b2019b8c5c9408c05516c3e0013ac1c8c6d1e11c8225d  inputs/run2/hacc/he_r2_17_nu1.inp
4ca0f941a853c31315b5bd21d131f9744da8387e0eadcb492ec067d51df28e4a  inputs/run2/hacc/he_r2_17_nu2.inp
42c0045315364e1fa2ea1b7b3bfa0766b82c4b3aa87d4bac1f7e787c8be6df43  inputs/run2/hacc/he_r2_17_nu4.inp
e16e68d1b609e676cf12409bb6b095210db1608aaafcbda974ed516304e1e793  inputs/run2/hacc/he_r2_17_nu6.inp
d13e0ed73a814efc8856a3d991d2b002c739ed6f35df05c6e7dc54d0a8d6f6aa  inputs/run2/hacc/he_r2_20_nu0.inp
b448363c5c4b8191a37bf994df054de9b76fbfaff982435214ca8ab41b4e2938  inputs/run2/hacc/he_r2_20_nu0p5.inp
0f6fa07633f27a72613a411af74fe94951e972347ee513b7f3bb076bb0f3de99  inputs/run2/hacc/he_r2_20_nu1.inp
b45b73fbd56e51c57f0d8dc17ce2e512db0d405bb9aeb314d350cfd20bc9a33b  inputs/run2/hacc/he_r2_20_nu4.inp
b5ce4ce10781c0f082ce6c7b610ec082cefabbb758250d0d974bf22f6d47004e  inputs/run2/sae/he_r2_14_nu0.inp
fe2480464e3732dcc5a7ce3cc5385d0ed38393b0ec94c677679d7087b33b1606  inputs/run2/sae/he_r2_14_nu0p5.inp
5af723d36a8fbefa437c61b8ebc9f01fd3de466f2a982e9720ab45083370c27a  inputs/run2/sae/he_r2_14_nu1.inp
b510fb4ac1c0b28ee453799ae94006e3bd25ea50575247ff84b1365ae285b943  inputs/run2/sae/he_r2_14_nu2.inp
946a26acea9b388bcac3d10f8dc01daf5b47f58824224c0006d097210dd896d0  inputs/run2/sae/he_r2_14_nu4.inp
8fe1bbb50d9dcffe44ffb60823a2f4f6fab43ea2e65bb237c04774410807106d  inputs/run2/sae/he_r2_14_nu6.inp
1cd686ecc61dc8ddfbec36ec3682a254edec8944f124e909244a1ca96fa76283  inputs/run2/sae/he_r2_17_nu0.inp
817e7b05f4cae79401018918426ff12387f18d424f0796796d7d2454116bd8d0  inputs/run2/sae/he_r2_17_nu0p5.inp
5f6e8939d0bf6ec560a3a6cc10343b7b1a774eca3e20e597a22a832941df9d8f  inputs/run2/sae/he_r2_17_nu1.inp
a25af99f2735f3c1603cbd9f60b4b5b5b948007d4d03dad1b53ce58c7e8a686b  inputs/run2/sae/he_r2_17_nu2.inp
40cd5cf33e9c9e17f890956b921f38f357489a9ef31994b2f0a6951e4775d840  inputs/run2/sae/he_r2_17_nu4.inp
affa83f84abcbdcbe13f2a73f272c61d34ae926b5fdc7728cc87dbfeaf73451a  inputs/run2/sae/he_r2_17_nu6.inp
6d496d1d9ce754e0a45c34ff29eb40e7f4bc00178b2598322582bfd0a7b1282d  inputs/run2/sae/he_r2_20_nu0.inp
8cc85cea1c37303202f854858aff105eb52c3402aad95c22a6b2c7055fcba078  inputs/run2/sae/he_r2_20_nu0p5.inp
9580ef20e2a904e9708b5a29007fec259af055d63ab20f369e81c5b5c67c7cf8  inputs/run2/sae/he_r2_20_nu1.inp
310a8689f53b17fbf5dcbfc1b3c3f1b7e88ec73816d6da7a17b21e1c9571b2fb  inputs/run2/sae/he_r2_20_nu4.inp
bda108bd498ad45de1c95778b54d2d071908410926a9de5c4a4213c53d73cca8  inputs/run2/sae/he_r2_nofield.inp
fb4d50e66a505af0c53d088dd6c82f4bab053835bbb33cf40af6ab3610c1a86a  inputs/run3/sae/ar_m0_r3_10eV_I1e11.inp
e2fa1009106ab145d9376a8eff9b127ae500c8bdb634ebeff36c6c4ceabc620d  inputs/run3/sae/ar_m0_r3_10eV_I1e12.inp
b438eba9d705b8096bf7a359d1442aa488ef0cb9723b3c2b933868d9209973e8  inputs/run3/sae/ar_m0_r3_10eV_I1e13.inp
088f6e31a211c06971e178e875e7561e82a9ed003ea35683d7d768e166c198b8  inputs/run3/sae/ar_m0_r3_10p5_nu0p5.inp
634e2349a365d640fcbe10fb091b94cb0792de28f2f1039a1ab2ea33bfc37886  inputs/run3/sae/ar_m0_r3_10p5eV_I1e12.inp
edb38c5f58a69658ece52db89ee9e22866b7411466d5c32817e521dbd4c5d1bc  inputs/run3/sae/ar_m0_r3_11eV_I1e11.inp
80240f37dfcf62a3f341c5d8e3596d7ceb2f0861597bcd6969d94ad5fe1e4b79  inputs/run3/sae/ar_m0_r3_11eV_I1e12.inp
241a40e007396e6c9daf40636b3c8d4dbacab72308ad3ac4a16fb0ecba887a24  inputs/run3/sae/ar_m0_r3_11eV_I1e13.inp
7c4123636207c31bd37d7aedac2f0e12ec7523a4dd86976241578e4b1ac7b61c  inputs/run3/sae/ar_m0_r3_8eV_I1e11.inp
97ef728725bd3a5ef7005b2477f2d4bc5facbaf01d6d0113c60d080158ab36ea  inputs/run3/sae/ar_m0_r3_8eV_I1e12.inp
33705bb5d558f4b474035f0bc9cdde35d9f2311f65916715ee182819dfd364c3  inputs/run3/sae/ar_m0_r3_8eV_I1e13.inp
d5137dad9256d1d7f15f9706df4d6e21b7626cadeabddaf49bcf6982ae2b90cb  inputs/run3/sae/ar_m0_r3_8p5_nu0p5.inp
c918869dd04fc00b4e077751ba0cad84a1ba5b34f3e3cf0b028732ec12566067  inputs/run3/sae/ar_m0_r3_8p5_nu1.inp
c69856a96c198d1a779f6889c7ee04b0e7b3d778788fcd3e311ef7f15f70765a  inputs/run3/sae/ar_m0_r3_8p5_nu1p5.inp
fb63b6210ae85ff7523365a2f66ce1683862db4960c078781acc3bf4e8b7c005  inputs/run3/sae/ar_m0_r3_8p5_nu2.inp
9869d0759a1bb1b4723d6250cee952d2d355953fa9e9eab5f7463836e81d1240  inputs/run3/sae/ar_m0_r3_8p5eV_I1e11.inp
a74e7a9216940d154cd743cf863372a15284307a4918588f09947e42570c0a38  inputs/run3/sae/ar_m0_r3_8p5eV_I1e12.inp
c72bd8136998750ed67d396a845909ef40a70cc3c419f865a4a526fb144fb7dc  inputs/run3/sae/ar_m0_r3_8p5eV_I1e13.inp
78d2c3696d40e5936af3e6be1abef6a2f66a9fc7a366322cec1495b259a3fad9  inputs/run3/sae/ar_m0_r3_9eV_I1e11.inp
9a2bea2cafce0dbae204f97dcf69903d0fde997aaceca093f25e23a2fcb4c337  inputs/run3/sae/ar_m0_r3_9eV_I1e12.inp
4f0d3b12e1f862cdff05d878eb2c6479ded4f1ea787c30e4f523de653d14206f  inputs/run3/sae/ar_m0_r3_9eV_I1e12_box120.inp
e97cf3420ea4c38de2ac895d9f114d0c83f1dfcb08d5b8b5a45a6152e0731db3  inputs/run3/sae/ar_m0_r3_9eV_I1e12_cut150.inp
e1d8b2e52822aec00c2341d921633d10b160966c4adbd743c8228764104cfe7b  inputs/run3/sae/ar_m0_r3_9eV_I1e12_dt23.inp
2a3126a375091bc35b27155af2dbb8e14c98140ffcd9af378ebd1a0dbda7f136  inputs/run3/sae/ar_m0_r3_9eV_I1e12_eta15.inp
28012985d13955d776045f7fb56f3049b36af4bc6c23e0209bba2c52e4927273  inputs/run3/sae/ar_m0_r3_9eV_I1e12_th23.inp
965d9b0873c76120cef3e4fb5d2868d119c6bf84cd9da9f4c778590d5d50ebf4  inputs/run3/sae/ar_m0_r3_9eV_I1e13.inp
94e02445b8bf35ea70bc4d097ea20d805414957b49468e2e69d107fcbcd1dfbc  inputs/run3/sae/ar_m0_r3_9p5_nu0p5.inp
91c775acd5f1eaa6ad619ef216ed04e955d1e2ce759d8a105e54e360f6e14572  inputs/run3/sae/ar_m0_r3_9p5_nu1.inp
cd5dc4585b793d1aac34f5cc0db8b2f9b7c3a424bc4b9ad88708f4471fcff5ad  inputs/run3/sae/ar_m0_r3_9p5_nu1p5.inp
c9f3a0aa93224998af50d44458c06cac093ef6cf07020cff84aac9f1cee34d24  inputs/run3/sae/ar_m0_r3_9p5_nu2.inp
67befca8f264c275efd25bc3dbcdb2504e8306b3542ec5bef0b8e3280a39d2c7  inputs/run3/sae/ar_m0_r3_9p5eV_I1e12.inp
2ff0ea737fa235e77de02987cb07dbec5d9219cdedfbe24dae6d2d131332711b  inputs/run3/sae/ar_m0_r3_nofield.inp
6887fe95d6926129db84ce236308f6aaaf37d2c9c20632931f04bce39e23662f  inputs/run3/sae/ar_m0_r3_nofield_box120.inp
f3454ae567bc4993013fd97d353691935ae3bccca507d0fefef3ae4908aa10fa  inputs/run3/sae/ar_m0_r3_nofield_cut150.inp
b9de6314ef686ae52fd0e58f93a1b3e55873bc14dcec1df476774d90c4b7b659  inputs/run3/sae/ar_m0_r3_nofield_dt23.inp
4c12eb33588afb0803d707459eed73df07883e9c86ec5c44fa8455249e64ed65  inputs/run3/sae/ar_m0_r3_nofield_eta15.inp
da97f17d1af34dbac911a987f25108c4cfb5012221ae0dc4a46c0504de1170fa  inputs/run3/sae/ar_m0_r3_nofield_th23.inp
0b852e08977f77a7a45ed7b1da046eee544ee502f91a6169b1b54d81cf3da57d  inputs/run3/sae/ar_m1_r3_10eV_I1e11.inp
6504a13a6a1dc6900af262661c665567999d0302a9640c6793c3ab23b161def2  inputs/run3/sae/ar_m1_r3_10eV_I1e12.inp
386e8d0caef0d69c8816ade712b64e076401fe505fe380ec26fdfdf2b8119564  inputs/run3/sae/ar_m1_r3_10eV_I1e13.inp
3b0013b7486523ba6e1636e49d2a1a699487fe4c5ed18571e970bff8d785b361  inputs/run3/sae/ar_m1_r3_10p5_nu0p5.inp
4384b774d96c4debdb113792b5acbe8fe67b680c362fd2a8e8fd93a41d933a59  inputs/run3/sae/ar_m1_r3_10p5eV_I1e12.inp
b084dacf02a4c4d589f888c1c05639fe2ff6cc9af9d851c49176619c08973c66  inputs/run3/sae/ar_m1_r3_11eV_I1e11.inp
c6638c1069adf5ce4ad17a5e0d78a03d308c32b43d933b08c39af30503aec671  inputs/run3/sae/ar_m1_r3_11eV_I1e12.inp
69c5388dbc988d55ede64e3c2746317f4195a843f2b728a8dfd1995368df4581  inputs/run3/sae/ar_m1_r3_11eV_I1e13.inp
02cc9383c8ac10c7b4ffdc2cc35496e5e7f7fac143bc9040dc7afac9e2de6da9  inputs/run3/sae/ar_m1_r3_8eV_I1e11.inp
54974bc2214cc0c015801fe2bd7156e40d23c0cde0afd1d073565843da453674  inputs/run3/sae/ar_m1_r3_8eV_I1e12.inp
b2c07b7b86514c9d8784df1d33405610cbd289b80cbfd192ff41ca9983711b15  inputs/run3/sae/ar_m1_r3_8eV_I1e13.inp
d6c1e6172cbe46003ac70106557d08623a2e2832a31b6b46f2942818cc1d13b3  inputs/run3/sae/ar_m1_r3_8p5_nu0p5.inp
760313eef885745ebe93f07ee79d75959a61d40fcb486fbeb1507bc6e2d4a64e  inputs/run3/sae/ar_m1_r3_8p5_nu1.inp
6c421ef9cd744220ce7bf5bc2179f2597cb107cb63c2e9b0f21e36e0566715c2  inputs/run3/sae/ar_m1_r3_8p5_nu1p5.inp
13e2a0cf7aa8d381b066afafb5d1e273c709a728422d00a5281dec29c8b7f213  inputs/run3/sae/ar_m1_r3_8p5_nu2.inp
adf71025dcb691aae185b3e5fbbc020ecd8f3c8c1f6b85a0147195d81cda144c  inputs/run3/sae/ar_m1_r3_8p5eV_I1e11.inp
01d595f47fc6d2166fa5310759c34ccefac9524e2f9d97193cd4bf73d2568ebb  inputs/run3/sae/ar_m1_r3_8p5eV_I1e12.inp
dd6dbca633c8f4074cc8ed96476ae9f7545f78beca8c734e87751f438aac9f3b  inputs/run3/sae/ar_m1_r3_8p5eV_I1e13.inp
6fd13721067768efd63878eb22295343eaf6a4bc277bdf3dfc5d370c1e07edf1  inputs/run3/sae/ar_m1_r3_9eV_I1e11.inp
7c5095fdc26b16102d6785dee8cc248ae31acdc11cca65224dae99bd1e88e336  inputs/run3/sae/ar_m1_r3_9eV_I1e12.inp
754604b4dda38acead2fe52f18b69f951989d4d1fee8ca1452de7cba7c03c287  inputs/run3/sae/ar_m1_r3_9eV_I1e13.inp
8d286ccb6ee76d0b146071e07b9909796f21ccd29fbcb9b58b8384d65af0f619  inputs/run3/sae/ar_m1_r3_9p5_nu0p5.inp
06609199cd69edf8d99843c2287d1589cb7ac8ba4aeb8df03bb6f6e8e73e5c1a  inputs/run3/sae/ar_m1_r3_9p5_nu1.inp
7ce1c7709017ad7eeec6fe5e1d233da41d8bc522118976a830a7fedbcd324809  inputs/run3/sae/ar_m1_r3_9p5_nu1p5.inp
20d96428b9c839effc692c76c9f09af070fb05faf3987519be4ca4cb85d8b939  inputs/run3/sae/ar_m1_r3_9p5_nu2.inp
779e1580b3edbac1110b45dbf5e1459a78ba1760eb1bb8c514f40a3a4ffd44ba  inputs/run3/sae/ar_m1_r3_9p5eV_I1e12.inp
ab8051c45f4c2f4dfdc6e1c4d6e50639817b94adf5f657dd105d59b4c29873ad  inputs/run3/sae/ar_m1_r3_nofield.inp
2e775d67393ef69bb5ac6b134864f1ba7877784061a53126451212c153fa1c92  inputs/run3/sae/eig_ar_m0_b120.inp
884c511397380ade09b5a2a9dab02397f557920fd600f18f19e261d00b12ff8d  inputs/run3/sae/eig_ar_m0_b80.inp
577f294467b3d990be4a092d93dc1cfb6e1c1661d326816576f27d7729a04eb9  inputs/run3/sae/eig_ar_m0_b80_eta15.inp
1ad817d699426370d513f1a6532276715e7d379b9d511c9ea66975841ca9a26e  inputs/run3/sae/eig_ar_m0_b80_th23.inp
29306e3ba7b364a86b8f32a757e2adfe8d28046f83c815400a941122abe89f66  inputs/run3/sae/eig_ar_m1_b120.inp
b8be38a719ecc09005e22f912570deff49dc9d034502237606399b01aa6971b2  inputs/run3/sae/eig_ar_m1_b80.inp
a8a1d91907c6830691c4478950c0390e264fd54af57d932032c84757ddff941d  inputs/run3/sae/eig_ar_m1_b80_eta15.inp
09fe29cf23623ce34bc1ddca65c0ffe1355a542cd4285c45f4e33b30efaa7ad9  inputs/run3/sae/eig_ar_m1_b80_th23.inp
4ec95a6c4d7e6f32b5ca32f4862f8acca3ceb84f74e27c3aaec2f9a327475203  inputs/run4/sae/ar_m0_r4_30fs_I1e10.inp
3e69cb9d3dc87cd1d201d3d2b4d79ad2539007485f6a07427b47962cf25e7fc1  inputs/run4/sae/ar_m0_r4_30fs_I1e11.inp
28f0b81c8358b72bad42cdbd2cb12e9e1a7bc9df58f0a330c4a3754e6d03ba02  inputs/run4/sae/ar_m0_r4_50fs_I1e10.inp
48fbf2a0b9c62ab64d8ef757930c7fae8de5393454a3b62ccc8961abf53fe703  inputs/run4/sae/ar_m0_r4_50fs_I1e11.inp
2f44365e4f994012e82ae6d1ddd7d7535d1ec0a141564fbb8f8f11710b3e6dad  inputs/run4/sae/ar_m1_r4_30fs_I1e10.inp
fb76e10dd411c45fc8cd41263bca465b773655629c4f3ad41520db14d4d50fe5  inputs/run4/sae/ar_m1_r4_30fs_I1e11.inp
17627fd891ce743e26f618e4ef25b6cf429de16fbc40b4d7a899900fd2300a49  inputs/run4/sae/ar_m1_r4_50fs_I1e10.inp
da8721ce00fa1d2064c0c35b7db3273c3274627faf7a57c5ffa3800465ab5240  inputs/run4/sae/ar_m1_r4_50fs_I1e11.inp
884c511397380ade09b5a2a9dab02397f557920fd600f18f19e261d00b12ff8d  inputs/run4/sae/eig_ar_m0_b80.inp
b8be38a719ecc09005e22f912570deff49dc9d034502237606399b01aa6971b2  inputs/run4/sae/eig_ar_m1_b80.inp
eb8b7d92940c895f8833b61d3085019575ebe70bc778584ae6592def4b026e17  inputs/run5/hacc/he_r5_10fs_c0_I1e12.inp
e42f284f94fe7eccb8efd2760b11dada6f107dbd642bd5253373fbc67be35fbd  inputs/run5/hacc/he_r5_10fs_c0_I1e13.inp
f44660f236ff23338282594dd61378938c46283a82f86f2047f23ab68755c5e7  inputs/run5/hacc/he_r5_10fs_c0_I1e14.inp
b4a5b42a774bb8260ab767e388855f49c4d2fe9e141a3df70e617c525b0e7d74  inputs/run5/hacc/he_r5_1fs_c0_I1e12.inp
471ce6405419f3223e9c291fe68f4b6cd012e8d4f277d8470a0fa3b5909e1b4a  inputs/run5/hacc/he_r5_1fs_c0_I1e13.inp
bbba6f4fa154cba24230a1759cffcddf5470d2aea787a81c002de9bea616a97a  inputs/run5/hacc/he_r5_1fs_c0_I1e14.inp
0363b7ecd98ebb2ec935a178270a415d3de17d54dba9f1180df10a0e4c6ae5cf  inputs/run5/hacc/he_r5_3fs_c-10_I1e12.inp
6a8c8ac6ed737b8d71fdd5d6c65c0fd3a83688d4cd246a6d3b2492ad984768c0  inputs/run5/hacc/he_r5_3fs_c-10_I1e13.inp
ae7ff502e67ad46e1314c30e8912669b157facc927b1b2ef83c1dd23a0e8cb43  inputs/run5/hacc/he_r5_3fs_c-10_I1e14.inp
89c3258c26c9c074a9f7b15328675f124b134d78fd03d98d5a004d1138ac53bc  inputs/run5/hacc/he_r5_3fs_c0_I1e12.inp
67b96b53b071d9dd811fa0035312a33d58e1edcd2f96303bdd22d2f92ea4b013  inputs/run5/hacc/he_r5_3fs_c0_I1e13.inp
e8f42864e2ab608434cc17cc2651df0f60acb89dfe21d3e4d6480c74d8166be2  inputs/run5/hacc/he_r5_3fs_c0_I1e14.inp
6f0280cfc823a912c66291da303b7296df58922620a6807a0174f357dff72185  inputs/run5/hacc/he_r5_3fs_c10_I1e12.inp
e11eaf542d518a417ea3e730dcac191045e186ba46ee2bb5e5d58801be4a91a8  inputs/run5/hacc/he_r5_3fs_c10_I1e13.inp
8ce10534f2697dc1b69d324beccac2ad6a0876ed34f467029b90640443db005d  inputs/run5/hacc/he_r5_3fs_c10_I1e14.inp
39ad86b8af088229524565cd9f906502663371b67c37ca2442e0a0ce64ff4bd0  inputs/run6/sae/he_r6_1fs_c-5_I1e12.inp
3b9ff5575b4205e0d8c80d9411e4b0d7b96bbc2b7ae145706aed35c0335b8439  inputs/run6/sae/he_r6_1fs_c-5_I1e13.inp
9cbbfa264ffc88f12c60cdae1ec4582516b2018e121e5160fec7d7bae9b1cb8a  inputs/run6/sae/he_r6_1fs_c0_I1e12.inp
389287b69f417f031441e080a2827dca18a836a442dc866b3d6c77f74a0b1fed  inputs/run6/sae/he_r6_1fs_c0_I1e13.inp
062efcf513aae72a633800ecac45c24db02d56f610ac1782ab0c6bb59d274343  inputs/run6/sae/he_r6_1fs_c5_I1e12.inp
2d438f4f93c9ac226fb4b0ff88a9b5d8be3954a57fdc311b18768cdc664cac8e  inputs/run6/sae/he_r6_1fs_c5_I1e13.inp
bfd0d0307f3a165954c3cbe4fb031396f4b8eaf8be94a55c7ec3e3ef0764e235  inputs/run6/sae/he_r6_2fs_c0_I1e12.inp
e914aba0749ef9c98750cd62a71b38d3dc9d7034d058fa8317b83c0373cc22b2  inputs/run6/sae/he_r6_2fs_c0_I1e13.inp
63f78bceea3084fdebe0c8a600abd3261746905c57415aebb3271c9edd5e41b5  inputs/run6/sae/he_r6_5fs_c0_I1e12.inp
0a54fb132848b23662091d83c56c7e6453cf1a87c5199b065eb144269ae58f44  inputs/run6/sae/he_r6_5fs_c0_I1e13.inp
ae4025499c5efa36fdc876d38dd3b824cd8df8e8b23faa447ab9772901178046  inputs/run6/sae/he_r6_nofield.inp
```
