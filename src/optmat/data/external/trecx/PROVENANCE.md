# tRecX single-active-electron TDSE yields (J. C. Travers, 2026)

Ionisation yields of argon (3p, m = 0 and m = 1) and helium (1s) in a 6 fs pulse at
800 nm against peak intensity, computed by J. C. Travers with tRecX (A. Scrinzi,
Comp. Phys. Comm. 270, 108146 (2022), https://doi.org/10.1016/j.cpc.2021.108146, as the
code's README gives it). Unpublished. The source record is `travers_2026_trecx`; its
data files are written from the two result tables here by
`tools/extract/travers_2026_trecx.py`. These are tier-3 validators (TC-439, TC-444):
nothing is fitted to them, and they rank below published calculations and
measurements.

## The code

- Upstream: `https://gitlab.physik.uni-muenchen.de/AG-Scrinzi/tRecX.git`, branch
  `master`, at `fd6da402331a9fbe4e705c5c87068f2c987622d1` (30 June 2025, "fixed problem
  with PotentialTruncate").
- The working copy the runs were made in is two local commits ahead of it, both
  J. C. Travers's:
  `05f9cdea58863050a4f4af2fd50d17b029b059a8` (9 June 2026 22:36:33 +0100, "fix build on
  mac os") and `96bc605d62737070959a9fc79ba9f4afb90f953a` (Tue Jun 9 23:41:03 2026
  +0100, "arbitrary field loading": `fieldFromFile.cpp/.h`, `pulse.cpp/.h`,
  `run_trecx.cpp`), the head when this was written (2026-09-28).
- **The executable** is dated 9 June 2026 23:31, ten minutes before the second commit.
  The build's own record (`GIT_DIFF`) says "Local modifications from git version:
  05f9cdea ... in files: SRCS/core/CPP/pulse.cpp, SRCS/core/pulse.h,
  SRCS/main/run_trecx.cpp". So the runs used `05f9cdea` with the then uncommitted
  changes of those files, which are presumably those committed as `96bc605d`; that
  they are identical was not checked. The changes add the loading of a field from a
  file; the runs here use tRecX's analytic `gauss` and `cos8` shapes, and the pulse
  each run wrote (its `Laser` file) was checked (below).

## The runs

51 runs, made on 11 June 2026 between 12:35 and 19:42 (the runs' own logs and the file
times), on one machine. Every input file is kept in `inputs/` as run (each equals the
`inpc` copy tRecX wrote into its run directory; checked 2026-09-28):

| Family | Runs | Target | Pulse | Yields |
|---|---|---|---|---|
| `ar_scan_00` ... `_09` | 10 | Ar 3p, m = 0 | `gauss` | chi-loss, norm-loss |
| `ar_m1_00` ... `_09` | 10 | Ar 3p, m = 1 | `gauss` | chi-loss, norm-loss |
| `ar_isurff_00` ... `_06` | 7 | Ar 3p, m = 0 | `cos8` | iSurff (and chi-, norm-loss) |
| `ar_m1isurff_00` ... `_06` | 7 | Ar 3p, m = 1 | `cos8` | iSurff (and chi-, norm-loss) |
| `he_scan_00` ... `_09` | 10 | He 1s | `gauss` | chi-loss, norm-loss |
| `he_isurff_00` ... `_06` | 7 | He 1s | `cos8` | iSurff (and chi-, norm-loss) |

Within a family the argon inputs differ only in the intensity (and the comments that
restate it). The helium inputs differ in the intensity and also in the box, the number
of radial elements that fill it, and the time step:

| Index | I (W/cm2) | BOX (bohr) | Elements in [20, BOX] | fixStep (au) |
|---|---|---|---|---|
| 0 | 1.0000e+14 | 80 | 14 | 0.0283 |
| 1 | 1.5445e+14 | 80 | 14 | 0.0272 |
| 2 | 2.3853e+14 | 80 | 14 | 0.0260 |
| 3 | 3.6840e+14 | 80 | 14 | 0.0246 |
| 4 | 5.6898e+14 | 80 | 14 | 0.0230 |
| 5 | 8.7876e+14 | 80 | 14 | 0.0212 |
| 6 | 1.3572e+15 | 100 | 19 | 0.0193 |
| 7 | 2.0961e+15 | 120 | 23 | 0.0172 |
| 8 | 3.2374e+15 | 150 | 30 | 0.0151 |
| 9 | 5.0000e+15 | 180 | 37 | 0.0129 |

Common to all (read from the inputs and the runs' outputs):

- **Model:** one active electron in a model potential; the other electrons are frozen
  into the potential. Tong and Lin's potentials (J. Phys. B 38, 2593 (2005), Table 1;
  the source record `tong_2005`), written out in the inputs:
  Ar `-(1 + 16.039 exp(-2.007 r) - 25.543 r exp(-4.525 r) + 0.961 exp(-0.443 r))/r`,
  He `-(1 + 1.231 exp(-0.662 r) - 1.325 r exp(-1.236 r) - 0.231 exp(-0.480 r))/r`.
  The parameters are those of `tong_2005` as curated.
- **Initial state and its energy in the model:** Ar 3p, <H0> = -0.579689 au (15.774 eV),
  the same for m = 0 and m = 1; He 1s, <H0> = -0.903801 au (24.594 eV). These are the
  ionisation energies of the model, not the measured ones.
- **Field:** linear polarisation along z, velocity gauge (`iLaserAz[t]<<D/DZ>>`),
  800 nm. The envelope is on the vector potential: the `Laser` files are reproduced by
  A(t) = (E0/omega) g(t) sin(omega t) to 1e-6 of the peak (Gaussian) and 5e-9 (cos8),
  with E0 = sqrt(I / 3.50945e16 W/cm2) au, so that the field is cosine-like, with its
  largest maximum E0 at the peak of the envelope (carrier-envelope phase 0). g(t) is
  exp(-2 ln2 t^2/FWHM^2) (`gauss`) or cos^8(pi t/T) with T = 1333.3 au (`cos8`), both
  with a full width at half maximum **of the intensity envelope g^2** of 2.24844
  optical cycles = 6.0 fs (248.05 au). The field, -dA/dt, is not E0 g(t) cos(omega t)
  (that differs by up to 7 % of E0). The intensity given is that of a continuous wave
  of amplitude E0.
- **Discretisation:** 30 associated Legendre functions in the polar angle (l up to 29
  for m = 0); finite elements in r, orders 14 and 12, on [0, 2], [2, 20] and [20, BOX],
  then an exponentially damped element of order 30 to infinity; infinite-range
  exterior complex scaling from r = BOX, angle 0.3 (Ar) or 0.15 (He); fixed-step RK4;
  `cutEnergy` 60, `operatorThreshold` 1e-9.
- **Propagation:** the Gaussian runs from -3 to +3 FWHM (-744.14 to 744.14 au), where
  the field is 2e-6 of its peak; the cos8 runs over the support of the pulse (-666.65
  to 666.65 au). The propagation ends where the pulse ends: no field-free time follows.

## The three yields

- **chi-loss**, 1 - <chi[0,20]>: the population outside a sphere of 20 bohr at the end
  of the run. It counts ionised electrons and also bound excited (Rydberg) population
  that reaches beyond 20 bohr; slow electrons still inside 20 bohr at the end are not
  counted.
- **norm-loss**, 1 - <norm>: the population absorbed by the complex scaling beyond
  r = BOX by the end of the run. Smaller than the chi-loss by what is still between
  20 bohr and the absorber.
- **iSurff**: the integral of the photoelectron spectrum, int sigma(k) k^2 dk of
  `spec_total` (trapezoid rule over the momenta tRecX wrote), cos8 runs only. It counts
  the continuum alone. The spectrum is cut at the largest momentum written; the
  integrand there is 2.7 % of its maximum for `ar_isurff_00` and `_01`, 0.5 % for
  `ar_isurff_02` and `he_isurff_06`, and below 0.13 % for the others (`all_runs.csv`,
  `spec_tail_frac`).

## The files

| File | What it is | sha256 |
|---|---|---|
| `ar_results.csv` | the yields of the four argon families per intensity, with derived columns (field, Keldysh parameter, shell and m-averaged combinations, ratios) | `8845a8759b4452a87003a43c712529ad2767167976b6f1f423915f9ba545e2c4` |
| `he_results.csv` | the same for helium, with the box and the largest norm during the run | `cbd8ec6f611cab9c54c0c458e8884d8d3fcb93f8ae30fdab59f661ac15636a77` |
| `all_runs.csv` | one row per run: settings read from the input, times, <H0>, the yields from `expec` and from `outf`, the spectrum's cut, the warnings, the norm before the pulse | `41cf430af92e605065048f096eae940cb968da59e731df0a33306dee97de6f14` |
| `extract.py.txt` | first pass: reads each run's input, `expec`, `Laser`, `spec_total` and `WARNINGS`, writes `all_runs.csv` | `5a85cb77322e6be49c94e5cbdff37d5bfc61a9eba9b52d3815fdc69653f4b8db` |
| `extract2.py.txt` | second pass: reads each run's `outf` (the time trace printed during the run), adds its columns to `all_runs.csv` | `e45c34afb094723abaa7d31fdb5ae7e5f482c30e99093d106b87b4dad1a4a888` |
| `tables.py.txt` | writes `ar_results.csv` and `he_results.csv` from `all_runs.csv` (yields from the `outf` columns, `%.4e`) | `18a5007da4585f5709ed511725b5f3bb6081de40bbe321b9b553e99d7c6b719a` |
| `README_ar_intensity_scan.md` | the README of the argon runs' directory, as it is | `29bc2a79746418bfe5ca70faefd37ba8897d37f514f8f0962489fd855bd22e6a` |
| `README_he_intensity_scan.md` | the README of the helium runs' directory, as it is | `35e765b09c6ca78fc2226b647ec20525b268f576ed525dd0a6941ce409aab8d1` |
| `inputs/*.inp` | the 51 inputs as run | below |

Two edits were made for publication, and the sha256 above are of the files as kept here:
the local directory of the runs in the two extraction scripts (`ROOT`) is replaced by a
placeholder, and one comment of the helium README no longer names how the runs were
launched. Nothing else differs from the files as they were used.

The three scripts and the three tables were written on 2026-09-28 (09:36) in
`tRecX/testruns/ppt_review/audit/`, from the run directories
`tRecX/testruns/ar_intensity_scan` and `he_intensity_scan`. They are copied unchanged
apart from the two edits above; the scripts carry the suffix `.txt` so that the
package's linters do not read them as part of optmat (they are kept as evidence). The
run directories themselves (the wavefunction checkpoints, `Laser`,
`expec`, `outf`, `spec_total`) are not archived here.

The READMEs are the working notes of the runs and compare with models that are not
part of this record; their statements on convergence are quoted below and are not
supported by retained runs.

## Checked (2026-09-28)

- Every yield of `ar_results.csv` and `he_results.csv` was recomputed from the run
  directories by an independent script (last row of `expec`; the trapezoid integral of
  `spec_total`) and agrees with the table to the digits written, for all 51 runs but
  one: `he_scan_03`, whose `expec` file holds only its first row. Its yields
  (1.0752e-04, 1.0263e-04) are from the last row of `outf`
  (t = 744.145: norm 0.99989737316998, chi 0.99989247944675), which was read by eye.
  (`he_scan_01` and `he_isurff_01`, `_03` have short `expec` files too, whose last row
  is the final time.)
- The intensity of every row is the one in its runs' inputs.
- The pulse, the potentials and the initial energies, as stated above.

## Caveats

1. **No convergence runs are retained.** Every run has one discretisation. The READMEs
   state, without files to show it: "iSurff is converged at Eta=30, radialPoints=150
   (Eta=40/250 changes it <1 %)"; "cos8 == gauss to <0.1 % on the yield, checked at
   8.8e14"; "5e15 converged independently at BOX=160 (96.5 %) and BOX=180 (96.2 %)";
   a box smaller than the quiver amplitude gave 5.8 % instead of 97 % at 5e15 W/cm2.
   No uncertainty can be stated from convergence.
2. **What the tables themselves show.** The chi-loss of the cos8 runs differs from
   that of the Gaussian runs at the same intensity by at most 1.6 % for yields above
   1e-7 (Ar m = 0: -0.5, -1.0, -1.1, +0.2, +0.04, +0.9 % at indices 1 to 6; Ar m = 1:
   +1.6, -0.7, +0.7, +0.5, +0.4 % at indices 2 to 6; He: -0.2, +0.7, +0.4, +0.2,
   +0.3 % at indices 2 to 6; the README's "<0.1 %" at 8.8e14 is 0.2 % here) and by
   6.5 % at 1.3e-8 (Ar m = 0, index 0). The iSurff yield is 0.68 to
   0.92 of the chi-loss of the same run for Ar m = 0, 0.64 to 0.96 for Ar m = 1
   (indices 1 to 6) and 0.77 to 0.88 for He (indices 2 to 6); the norm-loss of the
   Gaussian runs is 0.71 to 0.99 of their chi-loss for Ar m = 0, 0.81 to 1.00 for Ar
   m = 1 and 0.94 to 0.98 for He (indices 2 to 9). The three measures count different
   things, and their spread (up to 30 to 40 % at low yields, a few percent near
   saturation) is the least uncertainty of "the yield".
3. **The noise floor of the chi- and norm-loss.** Both are 1 minus a number close to
   1. Before the pulse (t = -450 au) the loss is already 1.2e-10 (Ar, Gaussian runs),
   and 4.2e-08, 3.4e-08, 2.7e-08 for He at indices 0, 1, 2, falling to 8e-10 at index 9
   (`all_runs.csv`, `norm_loss_prepulse`). The He yields at indices 0 and 1 (3.1e-07,
   4.8e-07) are of the size of that drift over the run, and the iSurff yields there
   are 1.2e-09 and 8.7e-08 (0.004 and 0.19 of the chi-loss): they are noise. Ar m = 1
   at index 0 (7.9e-10; iSurff 1.1e-10) likewise.
4. **Two unstable runs.** `he_scan_08` (3.2374e15 W/cm2) and `he_scan_09` (5.0e15)
   carry 1973 and 2902 warnings "norm>1, possible instability - try smaller scaling
   angle or smaller operator threshold"; the norm rose to 3.16 and 3.41 during the
   pulse. No other run has that warning (the largest norm is 1.0). Every Gaussian run
   carries the warning "time propagation starts or ends during pulse", which is the
   truncation at 3 FWHM.
5. **The propagation ends with the pulse**, so the chi-loss misses electrons slower
   than about 20 bohr per remaining time, and population still excited is counted by
   chi-loss whatever becomes of it.
6. **One electron.** The yields are those of one electron in one orbital (one m). How
   they combine into the yield of the atom (two electrons in He 1s; six in Ar 3p, two
   with m = 0 and four with |m| = 1) is a model: the columns `shell_*` and `mavg_*` of
   `ar_results.csv` are two such combinations, and are not curated.

## sha256 of the inputs

```
46718e54ac4b346caa20967cfd07bef2c91d6655bce438083a3fc74de3b77177  inputs/ar_isurff_00.inp
d805756a70028a3c3728b7b0150acefe3aa8b26686845fafaa7c3c96486826e5  inputs/ar_isurff_01.inp
5762f3fef4ef06328239225e293e2c9f335f093466fb4a7e59ff577313fc56d9  inputs/ar_isurff_02.inp
4f86d353195c3c082ab9305bb65799087f36d66af4991c465ac24e19d91a42a7  inputs/ar_isurff_03.inp
2d95bd7dbf6f6b4184adc1c65ba64abb70f10f085f3a17966b8abad2ed2b3905  inputs/ar_isurff_04.inp
cdbe01e4749bc4935c395a61d78335bfba25d1dc78f14791ac766bd333b4c300  inputs/ar_isurff_05.inp
4bdbf5c2b98a5301bed260954ccc00502ea7c31b5178b1c39edb02da8e5ef610  inputs/ar_isurff_06.inp
075710952e50b63c0dfd4fdec55c0fd14d062378336496943a9420abdaacf24f  inputs/ar_m1_00.inp
55fd2068294a331c677a58a5175ce00b9c87a1528ee567180d4c0872a8eb5121  inputs/ar_m1_01.inp
6ac28193ca77532e5816254538a7e138eb565666ae3b4ba260c045ccbdb547b3  inputs/ar_m1_02.inp
ad81198b1da2d7d1175857e1054f2626b783cdf0724f977b27e399e1d74ee1ce  inputs/ar_m1_03.inp
33c44e9e723c1d87f2f864e49a785379c16ca6ea9a81f817e923eed64cbba794  inputs/ar_m1_04.inp
652d529ef2129f4d903a506c647f8e23f948c3a71c753d82ad63e45f0de54aca  inputs/ar_m1_05.inp
cc38b9d8e94ae6c6ce58c2db58fe4813d4bd7656d86c2edb7f7efa794ca7a062  inputs/ar_m1_06.inp
dde863a630e38216350e971c55dabf59ca51552b40f8dd0c64efc6672a354b3b  inputs/ar_m1_07.inp
d3ab49c092b86dbc132142441e93eb38ba9eef290cdc18646e02493b45bfde82  inputs/ar_m1_08.inp
b0608ea2f6903a01b0a23950838b3062fb803a62a26f98b97cfa587beb019e23  inputs/ar_m1_09.inp
86029334072a082c43427e3bd9f20acc23cc2f94c4ee4ff3a2429968231f6be2  inputs/ar_m1isurff_00.inp
215f9eb869640d0e8a5477f47e8560a76da134a075ea1964954c13dc2c2d39d0  inputs/ar_m1isurff_01.inp
34ed94312b557cedf3aefd9b46aa7d31617264d901f4d4bced4e2d971a0b60fb  inputs/ar_m1isurff_02.inp
47626e2a85da9a4f7a44dbf6ef230c52adb1b49fdd5588875fed9bc6648cec5f  inputs/ar_m1isurff_03.inp
8ca2e683e62cbc1eaf4b754fe05affe135de587201adbd4f92d9ec77efa85e09  inputs/ar_m1isurff_04.inp
6b3479ff730b095807c37e84ea4bdadd135fe376229af24d726f82a26dffd076  inputs/ar_m1isurff_05.inp
b9fe8a3901aa0ba5e7eee34c4604dd60b27cf8bf6d2d9df39d2e53e9f49e4758  inputs/ar_m1isurff_06.inp
b87853cf769a77f637caf3e64c63845330c01ca3e83ad9161625d938de6fcde7  inputs/ar_scan_00.inp
c9024dff9c90a91bc6a73b55d7e0b9b20861d27057553e6ff9355be5daa6f639  inputs/ar_scan_01.inp
d05285aea37f6e81c618df6c19a796ae9fc547a8a2fad4eb96d7d2d23de3a0c1  inputs/ar_scan_02.inp
051fe7b0d54f5b4f1d30a80d150bf795038bc87d91a81937ca37c5dada4cd343  inputs/ar_scan_03.inp
5f03f4b97dcac3e8570febf6ba31e209a7887014b89d3ab07dc9031178cc2a16  inputs/ar_scan_04.inp
3caba6e1d8fbe58325b7dd51023beb5adc6b8a2d3c0f46fe68e98651a9a146c8  inputs/ar_scan_05.inp
075d77390feb0e2eecc00dce6f1af83784dc725658a470ca2fef2b4ba3ec7329  inputs/ar_scan_06.inp
6a0af512dbd4b49d39e73203790f4fd739d30a094b42ed9570918937164aa7d5  inputs/ar_scan_07.inp
aa4b63c907eb68ec4d2e576c2e576c310e432bc6058074843696eb1bd8629236  inputs/ar_scan_08.inp
3d81192180d216376a6edce2f98b753eff626a323c7ceb188b03f3501f114fa3  inputs/ar_scan_09.inp
7cbb3f4837173ec25b1c22cf7452ddd5f3910d2a68b9c6bed4d918b3bf117229  inputs/he_isurff_00.inp
0466e30c7e388606918639f28ef115fa3bf3e5599a75c1032f74c77de4c24e2a  inputs/he_isurff_01.inp
00df9be78257b0d34783ce0c4ac5a4eca24e2ca681e3ad0d7d7ed49c6fd676fd  inputs/he_isurff_02.inp
1d38307a9f5f2c50c89aa257de805cc538e4ef4769e3dc9acd017d4e0c8d4a43  inputs/he_isurff_03.inp
6335e3c49069d6515c880261a1e71172483392b8257e52449c200f2b8b6061e2  inputs/he_isurff_04.inp
8b6524adaca5214ecf84c5c078e4ca2a80f4b8adf881af76928aeed0d56a18f2  inputs/he_isurff_05.inp
5647fbef8d8c4dfc29fc401ff08a1714c0296db25b55aa3e7377669fd112c857  inputs/he_isurff_06.inp
e1ddca238b0d6d33bb2a5a483f17050b7d9b0250431ebf1521ed33538b68eb68  inputs/he_scan_00.inp
e2db49950104ac61ccd6825b313563ce3de96c8b0870a1e16ca3a56090301327  inputs/he_scan_01.inp
f5e9f4f099e0d4e7fef21b184bd68a584f3b4edf40661a56d1c81060c533690c  inputs/he_scan_02.inp
47e37338727bd173ab4fa8d607e8ab7702a10944811d3f871be001aaef44aba4  inputs/he_scan_03.inp
d2ad92e2607f81ee64aeb915bf7647414bc0bfc82670c6b79d5b3d5f73b133c1  inputs/he_scan_04.inp
cd85e9e5f12b322e0d28916f946f4db078cdf4d603d8e8831cb5f074c6361b1f  inputs/he_scan_05.inp
70dae74fde5ebaf0c31d3ebdcad9d01e1a8fec9829069db651026e6b5a7f3e58  inputs/he_scan_06.inp
63220e110656140fc98542f296b1253b083370dd7f57d4b3221dae557820e19b  inputs/he_scan_07.inp
2f7cd84200f796590903dbe6e1ee30fe3103464791217b4a010a352b2d6631b8  inputs/he_scan_08.inp
d5579c256016a12c694c73fa45c471ea4405fe597be06cc1d85999fdec7aa1ba  inputs/he_scan_09.inp
```

# Argon 3p, 800 nm, 42 fs Gaussian (September 2026)

Ionisation of one active electron (Ar 3p, m = 0, Tong and Lin's potential) by the
pulse of Wahlstrand et al. (PRL 120, 183901 (2018)): 800 nm, a Gaussian of 42 fs full
width at half maximum of the intensity, 17 peak intensities from 4.0e13 to 1.35e14
W/cm2, computed by J. C. Travers with tRecX on 28-29 September 2026. The source
record's datasets `ionised_gauss42_Ar_m0` are written from `ar_42fs/results.csv` by
`tools/extract/travers_2026_trecx.py` (TC-702). Tier 3, as the 6 fs runs above.

## The code

tRecX `a987d7c790e52b6d7f44a5615a665e0abc3e7f4a` (28 September 2026 16:58), one
commit of J. C. Travers's on top of `96bc605d` (above). The commit adds output
columns `OccDual{H0:k}`, |D_k^T psi|^2: the population of field-free eigenstate k in
the bi-orthogonal sense of the complex-scaled Hamiltonian, over the whole space (33
lines in `SRCS/core/projectSubspace.h`, `SRCS/core/CPP/timePropagatorOutput.cpp`,
`SRCS/main/tRecX.cpp`). tRecX's own `Occupation` prints (R_k^H S psi)* c_k over the
unscaled region only, which contains cross terms with the continuum: in a 6 fs test
its n = 6-8 populations drifted by 5-20 % after the pulse and n >= 9 went negative,
while the `OccDual` populations of every shell n <= 8 stayed constant after the pulse
to 1e-9. The norm, the expectation values and iSurff are bit-identical to the
unmodified build. Built as `Scalar` (no MPI) with gcc 13.3 in a clean worktree of that
commit; all 53 runs used that one executable on one workstation, one core per run
(column `commit` of `results.csv`: `a987d7c7` in every row).

## The runs

`ar_42fs/manifest.csv` lists the 53 runs and their settings; `ar_42fs/results.csv`
(one row per run) gives the settings, the initial <H0>, the pulse as run, the
populations by shell at the end of the pulse and 1500 au later (with `OccDual`, and the
`Occupation` values in the columns ending `_occ`), the chi- and norm-loss, iSurff, the
largest norm and the warnings. Families:

| Rows | Runs | What |
|---|---|---|
| `base_*` | 17 | the scan: every intensity, the base discretisation, 1500 au field-free after the pulse |
| `pend_*` | 7 | the first seven intensities, the propagation ending with the pulse |
| `eta40_`, `eta50_`, `nel20_`, `ord14_`, `box120_`, `dt010_`, `th020_`, `cut100_`, `cut150_`, `opthr11_` | 26 | convergence, at 4.0e13, 7.8e13 and 1.35e14 W/cm2: angular functions 40 (50), radial elements 20, order 14, BOX 120, step 0.010, ECS angle 0.2, cutEnergy 100 (150), operator threshold 1e-11 |
| `gauss_7p80e13` | 1 | tRecX's built-in `gauss` instead of the pulse from a file |
| `nofield` | 1 | no field: the numerical drift |

The base discretisation is that of the 6 fs argon runs (30 associated Legendre
functions, radial elements to BOX = 80 bohr, ECS angle 0.3, step 0.015 au, cutEnergy
60, operator threshold 1e-9), with the elements in [20, 80] split at 40 bohr. The
pulse, A(t) = (E0/omega) g(t) sin(omega t), g = exp(-2 ln 2 t^2/T^2), T = 42 fs, read
from a file sampled every ~0.05 au over [-3T, +3T] (g = 3.8e-6 at the ends; tRecX's
own `gauss` has no end, which iSurff refuses), then 1500 au without field. The
pulse file, 9.2 MB and written by the inputs' generator, is not archived: sha256
`43c05821e5f5f42eaf9166b020f45f5b7f5c727713f94e0a8a87bc4db3973107`.

## Checked (against results.csv, 2026-09-29)

- The initial <H0> is -0.579689334 au in every run (spread 8e-12).
- The pulse as run: the full width at half maximum of the intensity envelope is
  41.99977 fs and the peak field E0 to 4e-7 in every run with a field; the built-in
  `gauss` changes the ionised fraction by 3e-6 relative.
- The largest norm is 1.0 in every run. The warnings: only the run with tRecX's own
  `gauss` notes that the propagation starts or ends during the pulse (that Gaussian
  never ends); the three BOX 120 runs warn that an eigenproblem of H0 is satisfied only to
  1e-7, and they agree with BOX 80 to 0.03 %.
- The drift: the norm-loss of the run without a field is 4.8e-9, below every yield
  by a factor of 900 or more.
- The populations of the bound shells are constant from the end of the pulse to
  1500 au later (`pend_*` against `base_*`: the ionised fraction to 1e-13 at 7.8e13).

## The ionised fraction, and its uncertainty

The ionised fraction is 1 - (ground + every bound state up to n = 10) at the end of
the pulse, `ion_proj_le10_pe`. States up to n = 8 are exact in the BOX 80
discretisation and agree with BOX 120 to 1e-8; n = 9 to ~1e-5; n = 10 decays slowly
into the absorber. The sum to n = 10 agrees between BOX 80 and 120 to 0.03 %; the
n = 11 shell holds at most 0.7 % of the ionised population (at 7.8e13 W/cm2, BOX 120),
so the fraction may be high by that much.

Its uncertainty is the largest relative change among the convergence runs at the
nearest convergence intensity (between two, the larger): 2.6 % below 7.8e13 W/cm2
(cutEnergy 150 at 4.0e13), 1.3 % from 7.8e13 to 1.2e14 (the ECS angle at 7.8e13,
through the n = 9 and 10 shells), 0.1 % at 1.35e14. The changes of the angular and
radial discretisation, the step and the threshold are below 1e-5. The cutEnergy
change is one-sided and converging (60 to 100 to 150: +1.97 %, +2.59 % at 4.0e13;
+1.01 %, +1.34 % at 7.8e13; -0.02 % from 60 to 100 at 1.35e14): the yields are low by
about 3 % at 4.0e13 and 1.5 % at 7.8e13. The uncertainty of the ten intensities added
between the first seven was not measured; it is taken from their neighbours.

The excited bound population (n <= 10) is 1.2 to 30 % of the ionised fraction: below
5 % at most intensities, 10 to 30 % at the resonant ones (5.0e13, 5.5e13, 7.8e13,
1.05e14 W/cm2). The yield against intensity is not smooth: its local slope d ln Y /
d ln I swings between 2.4 and 10.7, low just past each channel closing (Ip + Up = n
omega near 4.8, 7.3, 9.9 and 12.5e13 W/cm2); the ten added intensities do not resolve
this structure.

## Not curated

The chi-loss (it counts Rydberg population beyond 20 bohr as ionised: 1 to 29 % more),
the norm-loss (the n = 9 to 11 population leaking into the absorber: up to 12 % more at
the resonances), iSurff (1.4 to 23 % less; see below), the populations by shell, and
the convergence rows.

iSurff falls short (TC-703) because of the Coulomb tail of the untruncated potential at the
tSurff surface (80 bohr), which tSurff assumes free of any potential. Checked with the
same pulse and the potential switched off between 60 and 80 bohr, at 4.0e13, 5.0e13,
7.8e13 and 1.35e14 W/cm2 (the same commit built with GCC 16 on another machine). The
truncated potential has 41 bound states, all of them projected, so its 1 - bound is
complete. iSurff then equals 1 - bound to 0.04 to 0.15 % below 1.35e14 and to 1.7 %
at 1.35e14, and the norm-loss equals it to 0.5 % or better. The spectrum's trapezoid
sum on every second momentum point changes by 0.8 % at 1.35e14, so part of that 1.7 %
is the momentum grid. The check validates projection onto the bound states as a
measure of ionisation. It does not correct the curated fraction: truncation removes
the Rydberg states above n of about 6 and moves the resonances, and it changes the
ionised fraction by -11 to +5 %.

## Files

The inputs as run (106: the runs and their field-free spectra), the scripts that made
the inputs, the table and the report, and the report itself are kept with the runs,
outside this repository:

```
9cc1cd90d094ff67ae4dadb432930c22a329e9e30b326f6ccdd462d675438db3  ar_42fs/results.csv
6a140f38922766a363ae93a53d50bd73df048cc0e1829148b7093a7ddbb24468  ar_42fs/manifest.csv
3c8b76b0a43f63503842015d22614d886d3536f252cf153dc34ced54449bd103  the 106 inputs, concatenated in name order
a96005c0b849ef4b1b3928a3ce0e1fe7dd65185b08beefe46664902d65c065c3  make_inputs.py
78c0d11162a167b2f5d01976fd2a4487766bc621b4bb6f6e23d718f4e112fe8e  extract.py
ea12ac850d938691bb9b2a3d64a82a609562627e89b91ae3d2b4dcfd62dd516c  label_states.py
1422e0e26bfe56c98132c243c556b596239134ed1f801ffb459991c1cf3fdf19  tables.py
491aa24f02adb73372a93a701d3a36ada480408ddd8b90fc430d77abab8e81af  check_trunc.py
f075206318dd4063b4d65fbc1e4f6d3f6aebbe1ec80fab845543c86da06e37e5  check_trunc.csv
0099b099a3ff3f33cf5459e3aeef6d47d4f9fc40338d5bd51a667e0da472ccfc  REPORT.md
```

# The 6 fs runs again: Ar 3p (m = 0, m = 1) and He 1s (September 2026)

The 6 fs scans of the June runs above, repeated by J. C. Travers with tRecX on 29
September 2026 so that the yield is the ionised fraction, with convergence runs. The
source record's datasets `ionised_gauss6_Ar_m0`, `ionised_gauss6_Ar_m1` and
`ionised_gauss6_He` are written from `scan_6fs_v2/results.csv` by
`tools/extract/travers_2026_trecx.py` (TC-812). Tier 3.

## The code

tRecX `a987d7c790e52b6d7f44a5615a665e0abc3e7f4a`, the commit of the 42 fs runs, built
as `Scalar` with GCC 16.2 (OpenBLAS 0.3.34, FFTW 3.3.11, arpack-ng 3.9.1) on a
cluster; that build reproduces the workstation build of the 42 fs runs to about 1e-9
relative. All 122 runs used the one executable, single-core.

## The runs

`scan_6fs_v2/results.csv` holds one row per run: the settings (family, row, intensity,
discretisation), the initial <H0>, the populations by shell at the end of the pulse and
at the end of the run (`OccDual`, and tRecX's `Occupation` in the columns ending
`_occ`), the chi-loss to 20, 40 and BOX bohr, the norm-loss, iSurff (600 momentum
points), the pulse as run, the largest norm, the warnings and the run time.

| Rows | Runs | What |
|---|---|---|
| `ar_m0_NN`, `ar_m1_NN` | 20 | the ten intensities of the June runs, 1.69e13 to 8.44e14 W/cm2 |
| `he_NN` | 8 | rows 0 to 7, 1.0e14 to 2.1e15 W/cm2 (rows 8 and 9 of June were unstable) |
| `<run>_<variation>` | 91 | convergence at Ar m = 0 rows 0, 1, 4, 9; Ar m = 1 rows 0, 2, 4, 9; He rows 0, 4, 7: angular functions 40 (`eta40`; 50 at the highest row), radial elements x 20/14 (`nel20`), order 14 (`ord14`), BOX + 40 (`boxp40`), step x 2/3 (`dt23`), ECS angle x 2/3 (`th23`), operator threshold 1e-11 (`opthr11`), cutEnergy 150 (`cut150`) |
| `<family>_nofield` | 3 | no field: the numerical drift |

The discretisation is that of the June runs (Tong and Lin's potentials, velocity gauge,
30 associated Legendre functions, BOX 80 bohr, 100 and 120 for He rows 6 and 7; ECS
angle 0.3 for Ar, 0.15 for He), with the elements in [20, BOX] split at 40 bohr and
cutEnergy 100. The pulse is the June Gaussian, A(t) = A0 g(t) sin(omega t), 6 fs full
width at half maximum of g^2, read from a file over [-3T, +3T] (the same file and code
as the 42 fs runs), then 1000 au without field. The step is 2/Emax at cutEnergy 100,
capped at 0.015 au (Ar; at the cap except rows 8 and 9) and 0.03 au (He: 0.0129 to
0.0176 au; 0.00785 at row 7).

## Checked (against results.csv and the runs' report, 2026-10-01)

- Every run completed; the initial <H0> is -0.579689334 au (Ar, spread 5e-11) and
  -0.903801106 au (He, spread 1.3e-10).
- The drift: the norm-loss of the run without a field is 1.0e-9 (Ar) and 3.2e-8 (He)
  in 2488 au. A yield within a factor 10 of it is an upper bound only, and is excluded:
  Ar m = 0 row 0, Ar m = 1 rows 0 and 1, He rows 0 and 1.
- The chi-loss to 20 bohr reproduces the June runs' to 0.99 to 1.03 in every row above
  the floor, except He row 2 (0.885), where June's numerical loss was a large part of
  the value: June's He values at rows 0 to 2 are mostly numerical loss.
- The largest norm is 1 at every print step except He row 7 and its variations, where
  it rises to 1.085 (1.04 to 1.12) at the peak of the field. It disappears with a finer
  radial basis, a smaller ECS angle or a larger box, and changes the ionised fraction by
  at most 0.05 % when it does.
- The warnings: a memory advisory for 600 momentum points (every run with a spectrum);
  the eight argon BOX 120 runs, that an eigenproblem of H0 is satisfied to 1e-7 and its
  projector to 1e-9 (as the 42 fs BOX 120 runs; they differ from BOX 80 with the pattern
  of the Rydberg states).
- Helium at the June runs' cutEnergy 60 was not converged in the step (a step 2/3 as
  long changed the yield at 5.7e14 W/cm2 by 67 %); hence cutEnergy 100, converged in
  the step to 0.03 %.

## The ionised fraction, and its uncertainty

The ionised fraction is 1 - (ground + every bound state up to n = 10) at the end of the
pulse, `ion_proj_le10_pe` (n <= 9 at He row 7, where the n = 10 states hold at most
0.05 % of the yield). The n = 11 shell, resolved at BOX 120, holds 0.03 to 1.5 % of
the ionised population. Against the June measures, over the rows above the floor: the
chi-loss to 20 bohr is 1.0 to 1.46 times it (the Rydberg states beyond 20 bohr), the
norm-loss 1.0 to 1.10 times, iSurff 0.82 to 1.02 times (TC-703).

Its uncertainty is the largest relative change among the convergence runs at the
nearest convergence rows below and above, counting only convergence rows above the
floor, and at least the drift over the yield (`six_fs_uncertainty` in the extraction
script reproduces the report's values):

| Family | Rows | Uncertainty | Set by |
|---|---|---|---|
| Ar m = 0 | 1 to 3; 4 to 8; 9 | 1.7 %; 0.6 %; 0.05 % | ECS angle at row 1; BOX at row 4 (both through n = 9, 10); angular functions |
| Ar m = 1 | 2; 3; 4; 5 to 9 | 1.1 %; 0.7 %; 0.4 %; 1.65 % | drift; step at row 2; BOX; angular functions 40 at row 9 (50: 0.6 %) |
| He | 2; 3, 4; 5 to 7 | 0.8 %; 0.4 %; 1.8 % | drift; angular functions 40 at row 4; at row 7 (50: 0.6 %) |

The rows between convergence rows take the larger of their neighbours'; the angular
functions matter only at the highest intensities, so this is probably conservative
there.

## Files

The inputs as run, the run directories and the label spectra are on the cluster with
the runs; the scripts, the tables and the report are kept with the runs, outside this
repository:

```
d3813b030ff9a159909e21d5488e29e2be2d74cd6fb87a74dfe061f9f2ef4df3  scan_6fs_v2/results.csv
a376e28c46c59d4118ce1ef4d39dde1458ca0197755be0bdb316bcd86ac4d483  make_inputs.py
7bfd526389b3c3f5ab8160f7a2a983b0419cb857e62f6a4545c967d510db2f07  extract.py
447d690dccf705566f4550927da46f256ede5e4d54e7ad218fab9222c447b385  label_states.py
e8d03499ed3c673b29126a1b6f74946163a2f84fffc2ee456f3cb50064cb6daf  tables.py
2e9efd8fa3297e12417a90aed6de5b0567ac646c3288885c3db9e7cf5e9eaa49  tables.md
2099ef0eb306697777a3b27951eaceaf8c95e9f084dc98c8f9cfbd9284ac7e45  run_scan_pack.slurm
a5b656958022ec531baa7b9a9f3a04dbfce08517ac2972c98f13915c16866522  REPORT.md
```

# Static rates of Ar 3p in Tong and Lin's potential (September 2026)

The static ionisation rate of one electron in the Ar 3p orbital, m = 0 and m = 1, at
the seven fields of Zhang et al.'s Table IX up to 0.2 au (Phys. Rev. A 90, 043410
(2014); `zhang_2014`), computed by J. C. Travers with tRecX on 29 September 2026. The
datasets `static_Ar_m0` and `static_Ar_m1` are written from `static_rates/results.csv`
by `tools/extract/travers_2026_trecx.py` (TC-813). Tier 3.

## The code and the method

tRecX `a987d7c790e5...` (as above; its population output is not used here), `Scalar`
build on the workstation. The Hamiltonian is H0 + F z (length gauge) with Tong and
Lin's argon potential; m = 0 and m = 1 (`expIm[1]`). The discretisation is that of the
42 fs runs (30 Legendre functions; finite elements of orders 14 and 12 on [0, 2],
[2, 20] and [20, 40] bohr), with exterior complex scaling from 40 bohr at an angle of
0.4 (Zhang et al.'s angle; they scale from the origin). The resonance is traced from
F = 0 (the field-free 3p energy, -0.5796893 au, in every run) in steps of at most
0.004 au, each step the eigenvalue nearest the extrapolated guess by inverse iteration
(ARPACK); the field is the last step, so it is exactly Zhang et al.'s. The rate is
Gamma = -2 Im E. Nothing from Zhang et al. is used but the fields.

## The runs

`static_rates/results.csv`: the 14 base traces (2 orbitals x 7 fields) and 84 traced
variations, with the eigenvalue, the overlap of successive eigenvectors (at least
0.905 in the base traces), Gamma and, for m = 0, Zhang et al.'s value and the ratio.
`static_rates/results_nearest.csv`: 112 direct solves at each field of the varied
discretisations for the eigenvalue nearest the base one: angular functions 40, radial
elements 8, radial order 14, ECS angle 0.3 and 0.6, ECS radius 30, 20 and 60 bohr.

## Checked (2026-10-01)

- Every run exited normally and wrote its trace; the only warnings are 'inverse
  iteration not converged' in two traced runs with the ECS radius at 60 bohr (not used)
  and the by-design note of the direct solves that they solve at one field.
- m = 0 against Zhang et al.'s Tong-Lin column (twice the printed Gamma/2): 0.04 to
  2.4 % at the six fields optmat keeps, each within the uncertainty below except the
  weakest field (0.11 % against 0.016 %). At the excluded field, 0.12177 au, Gamma/2 =
  1.3232e-2, which is the printed Muller cell (1.3229e-2) to 0.02 %: the row's cells
  are very likely displaced (TC-813).
- The uncertainty is the largest change of Gamma among the direct solves, the ECS
  radius at 60 bohr left out: m = 0, 0.016 % at 0.053 au rising to 1.4 % at 0.18 au and
  7.7 % at 0.197 au (the ECS angle 0.3); m = 1, 0.006 to 0.86 % (`static_uncertainty`
  in the extraction script). Moving the ECS radius outwards to 60 bohr loses the broad
  resonance among the discretised continuum states (successive overlaps fall to 0.4 to
  0.75 above 0.1 au) and changes Gamma by up to 108 % at 0.16 au and above; moving it
  inwards to 30 and 20 bohr changes it by at most 1 % below 0.19 au.

## Files

```
2a2a4664541de84e778dbca76e2f99fb31d676d2ff0c48ddba5da0af60b04577  static_rates/results.csv
69cf7f1d98a311af4458e78de7c149e93cc9bc43b2d9e7ac8d9236c9a5f6b2c6  static_rates/results_nearest.csv
635e6e65b9fb53576bee5bcd432ff8fc2e5dc704eee1376673048d4c912e0a1e  the 98 traced inputs, concatenated in name order
ec1fc11e7aebde6e08d4afc52c6402607798ecd5bf88b66b33bd34f3d2a6a999  the 112 direct inputs, concatenated in name order
1e4fd914898ef1aa962cbc4a0f188e805bb20a58f9abed37d9304fc1b507066c  make_inputs.py
fc41dd6314be575c2ada5466e8d5a961f9c98a3556f07ae08e726b3351ef9258  extract.py
4dfb80e0c96406f24d4a1fe2622d213fdf055a28dc15b02129092f490f4dbbb6  tables.md
65c719db6c8b6fed740c05980c53bcb9d7ae26257ccedf478bdebb5a73c55b5d  REPORT.md
```

# Argon 3p at 400 nm (October 2026)

The A2 set-up at 400 nm: Ar 3p, m = 0 and m = 1, a Gaussian of 10 fs (full width at
half maximum of the intensity) on the vector potential, ten peak intensities from 3e13
to 1.5e15 W/cm2, then 1000 au without field. The datasets
`ionised_gauss10_400nm_Ar_m0` and `_m1` are written from
`scan_wavelength/results_400.csv` by `tools/extract/travers_2026_trecx.py` (TC-824).
Tier 3.

## The code and the runs

tRecX `a987d7c790e52b6d7f44a5615a665e0abc3e7f4a`, the GCC 16 cluster build of A2;
single-core runs, 53 in all (m = 0 at ten intensities with convergence at rows 0, 4, 9;
m = 1 at the ten without its own convergence; a run without a field per family), every
one complete. Window [-1240.24, 2240.24] au; fixStep at its cap of 0.015 au throughout.
The convergence list is A2's plus fixStep x 1/2 and cutEnergy 200.

## Checked (against results_400.csv and the runs' report, 2026-10-01)

- Initial <H0> -0.5796893340 au in every run; the pulse's FWHM 10.0000017 fs and peak
  field to 5e-11.
- The drift is 1.39e-9 in 3480 au; every yield is at least 250 times above 10 times it.
- The largest norm is 1 at every print step except the three cutEnergy 200 runs (1 +
  3e-10).
- The uncertainty, by the rule of the A2 section above (`six_fs_uncertainty` with this
  table): 0.32 % to row 3 (the ECS angle), 0.17 % to row 8, 0.046 % at row 9; m = 1
  takes m = 0's at the same row. The time step and cutEnergy are converged (0.001 %
  and, one-sided, 0.22 %).

## Files

```
73cb9601ac3fc491d59f5d8edbf6a3a371303ffdde402524af20af7d80db0f27  scan_wavelength/results_400.csv
78173376ac259a92080e9fe0cdd3ae8ee23c8f73e80509617ac28d9d42272ee3  make_inputs.py
4ebbde09bb5cffaa6539a22b2e68cd9a5e5d2ea7047ac5b4d64e040b2da65b17  extract.py
447d690dccf705566f4550927da46f256ede5e4d54e7ad218fab9222c447b385  label_states.py
aacdf5659fb27e9211fb30a6d69091feeb37c737f3d5afda2b91e6954674274b  tables.py
5b1b3972b6ac1583336d693d3f35b1c955c4983968c441e811ee8cb33d4cd025  tables_400.md
2099ef0eb306697777a3b27951eaceaf8c95e9f084dc98c8f9cfbd9284ac7e45  run_scan_pack.slurm
06a08d8056de7152098b0c1c42c966de6a2412b5c3259aea05f65fad179195a2  REPORT.md
```

# Argon 3p at 1030 nm (October 2026)

The 400 nm set-up at 1030 nm: Ar 3p, m = 0 and m = 1, the same 10 fs Gaussian (2.9
cycles), ten peak intensities from 1e13 to 5e14 W/cm2, then 1000 au without field.
BOX 80 at rows 0 to 7 and 120 at rows 8 and 9, where the quiver radius (49 and 61
bohr) is more than half of 80. The datasets `ionised_gauss10_1030nm_Ar_m0` and `_m1`
are written from `scan_wavelength/results_1030.csv` by
`tools/extract/travers_2026_trecx.py` (TC-826). Tier 3.

## The code and the runs

tRecX `a987d7c790e52b6d7f44a5615a665e0abc3e7f4a`, the build of the 400 nm runs; 59
single-core runs (m = 0 at ten intensities with convergence at rows 1, 2, 5 and 9; m = 1
at the ten without its own convergence; a run without a field per family), every one
complete. fixStep at its cap of 0.015 au except row 9 (0.0141). The convergence list
is the 400 nm one without fixStep x 1/2 and cutEnergy 200, plus the ECS angle x 1/2
(0.15) and, at row 9, BOX 160 and Eta 50; state labels for BOX 160 were made for it.

## Checked (against results_1030.csv and the runs' report, 2026-10-02)

- Initial <H0> -0.5796893340 au (BOX 80), -0.5796893339 (BOX 120 and 160); the pulse's
  FWHM 10.0000038 fs and peak field to 6e-11.
- The drift is 1.39e-9 in 3480 au. Rows 0 and 1 of m = 0 and 0 to 2 of m = 1 are
  within a factor 10 of it, and are excluded (upper bounds only).
- The uncertainty, by the rule of the A2 section above (`six_fs_uncertainty` with this
  table): 4.2 % at row 2 (2.4 times the floor; time step and cutEnergy partly drift),
  3.3 % at rows 3 and 4 (from row 2), 0.94 % to row 8 (the ECS angle x 1/2, through n
  = 9 and 10 states left unlabelled), 0.22 % at row 9; m = 1 takes m = 0's at the same
  row, and 7.1 % (its drift) at row 3.

## Files

```
c8b2fc4c1fe3aa51c8cdd6d2d0b5fbff344b63679be3eca69ee321e921b274c3  scan_wavelength/results_1030.csv
c90da110c49a18ab99e92f79d336ece5fbb12e8832f0c3e0ddb1a0be5b0809a5  make_inputs.py
4ebbde09bb5cffaa6539a22b2e68cd9a5e5d2ea7047ac5b4d64e040b2da65b17  extract.py
447d690dccf705566f4550927da46f256ede5e4d54e7ad218fab9222c447b385  label_states.py
c0555bf3959c935557b30c995fcf4cb5db926b542aa625747f575133558c0508  tables.py
ce513012bb043500bc0798654740fb70b5d3db3b688d0202fa1e1674dcfd97f8  tables_1030.md
2099ef0eb306697777a3b27951eaceaf8c95e9f084dc98c8f9cfbd9284ac7e45  run_scan_pack.slurm
5ef944e486596f1f152b69103a40ab442294a4d4c8befba85e296f059a4f79d0  REPORT.md
```

# Argon 3p at 800 nm with the 10 fs pulse (October 2026)

The 400 and 1030 nm set-up at 800 nm: Ar 3p, m = 0 and m = 1, the same 10 fs Gaussian
(3.75 cycles), ten peak intensities from 1.5e13 to 8e14 W/cm2, then 1000 au without
field; BOX 80, and 120 at row 9. The datasets `ionised_gauss10_800nm_Ar_m0` and `_m1`
are written from `scan_wavelength/results_800.csv` by
`tools/extract/travers_2026_trecx.py` (TC-830). Tier 3.

## The code and the runs

tRecX `a987d7c790e52b6d7f44a5615a665e0abc3e7f4a`, the build of the 400 nm runs; 67
runs in the table (m = 0 at ten intensities with convergence at rows 0, 1, 5 and 9 with
the 400 nm list plus the ECS angle x 1/2; m = 1 at the ten; a run without a field per
family), every one complete. fixStep 0.015 au (0.0142 at row 9).

## Checked (against results_800.csv and the runs' report, 2026-10-03)

- Initial <H0> -0.5796893340 au; the pulse's FWHM 10.0000037 fs and peak field to
  1e-11.
- The drift is 1.39e-9 in 3480 au; row 0 of m = 0 and rows 0 and 1 of m = 1 are within
  a factor 10 of it and excluded.
- The uncertainty by `six_fs_uncertainty`: 2.3 % to row 4 (the ECS angle x 1/2, through
  n = 9 and 10), 0.79 % to row 8, 0.096 % at row 9; m = 1 takes m = 0's.
- The BOX 160 run at row 9 warned "norm>1, possible instability": the norm of the
  unscaled region exceeded 1 at three print steps near the peak of the field (1.14 to
  1.56) and was 1 before and after; its end-of-pulse yield agrees with BOX 120 to
  0.028 %. Kept as the box check, with this caveat.

## Files

```
e8706fdf8407aaafd33532b7873adc3db23a71b1ef16354e753963737b09e5b1  scan_wavelength/results_800.csv
c90da110c49a18ab99e92f79d336ece5fbb12e8832f0c3e0ddb1a0be5b0809a5  make_inputs.py
4ebbde09bb5cffaa6539a22b2e68cd9a5e5d2ea7047ac5b4d64e040b2da65b17  extract.py
c0555bf3959c935557b30c995fcf4cb5db926b542aa625747f575133558c0508  tables.py
c901d7ea253d47c9b9b484cb7304c108680ebdcd19975dfb9e3ac103de6f9d9d  tables_800.md
09428b2dbcc94d690718bd80a76b32ffd05d3df1031e6370003c7297181dcb63  REPORT.md
```

# One electron at haCC's ionisation energies (October 2026)

Archived in TC-830; curated as datasets in TC-847 (`normloss_ipmatched_*`,
`ionised_ipmatched_*`), with haCC's runs below: Tong and Lin's argon potential with its
short-range terms scaled by beta (Zhao et al. 2016, Eq. 3) so that the 3p energy is
haCC's configuration-interaction ionisation energy, 0.560275 au (beta = 0.9930278467,
`results_ci.csv`), or the effective energy of haCC's 3-channel mMax 1 set-up, 0.563940
au (beta = 0.9943617816, `results_c3m1.csv`); A2's 6 fs pulse at haCC's three
intensities (4.03e13, 9.61e13, 2.29e14 W/cm2), m = 0 and m = 1, the ionised fraction
and the norm-loss 1000 au after the pulse, with A2's convergence at 9.61e13 (0.62 % for
the ionised fraction; the norm-loss -1.6 % at BOX 120). They are the like-for-like
reference of the all-electron haCC runs.

## Files

```
6730b12c776ae987e0c26425d9fda9d4accb4ba01c8a85bedbaea1fb47989ead  ip_matched/results_c3m1.csv
2f490e3c2e8095c304775700e40817a9c1e4b8c38b2f82dceb8c39bf383a0c95  ip_matched/results_ci.csv
24abe9def220d0c980920c6b4b684d56e08711aae3b52f50962c62e636b767a8  beta_search.py
1a1b00784f1ce983f341a47fa39214b1a849ba86a2d3f1346d4482c15dcd0076  beta_result_ci.txt
959199a76afded505c5f04aba248900d62bd62b1f076e49e39a53bbf0cacfc77  beta_result_c3m1.txt
9b1b683bea2e62b7b63c1fec041585f8bb6c3af798eb112b1e2f3ef40a4fd278  make_inputs.py
4ebbde09bb5cffaa6539a22b2e68cd9a5e5d2ea7047ac5b4d64e040b2da65b17  extract.py
c0555bf3959c935557b30c995fcf4cb5db926b542aa625747f575133558c0508  tables.py
aec931b3c4c1ba3dc2e726e6af4b91aab65716ef3ecb7eaf74e24328313029c0  tables_ci.md
a6a1398ff49917bf56938941df3ac477bda8dc2d2f1914806c9653bb1339a86f  tables_c3m1.md
09c31776511935c48f2b3d4a7fc0fe42a365bf88a4d82a23229ec167f00aec68  REPORT.md
```

# Argon 3p in Wallace et al.'s pulse (October 2026)

A2's set-up with Wallace et al.'s (2016) pulse: Ar 3p, m = 0 and m = 1, a Gaussian of
5.5 fs (full width at half maximum of the intensity, 2.06 cycles) on the vector
potential at 800 nm, carrier-envelope phases 0 and pi/2, twelve peak intensities from
5e13 to 6e14 W/cm2, then 1000 au without field; BOX 80, and 120 at rows 9 to 11. The
datasets `ionised_wallace_cep0_Ar_m0`, `_m1` and `ionised_wallace_cep90_Ar_m0`, `_m1`
are written from `scan_wallace/results_cep0.csv` and `results_cep90.csv` by
`tools/extract/travers_2026_trecx.py` (TC-836). Tier 3.

## The code and the runs

tRecX `a987d7c790e52b6d7f44a5615a665e0abc3e7f4a`, the build of the 400 nm runs; 60 runs
at phase 0 (m = 0 with convergence at rows 0, 6 and 11 with the 400 nm list plus the
ECS angle x 1/2 and, at row 11, Eta 50; m = 1; a run without a field per family) and 26
at phase pi/2 (no convergence runs), every one complete. Window [-682.13, 1682.13] au.

## Checked (against the tables and the runs' report, 2026-10-03)

- Initial <H0> -0.5796893340 au; the pulse's intensity FWHM 5.500002 fs.
- The drift is 9.5e-10 in 2364 au; every yield is far above 10 times it.
- The largest norm is 1 (+2e-10 in the cutEnergy 200 runs).
- The uncertainty by `six_fs_uncertainty` (phase 0) and `wallace_uncertainty`: 1.5 % to
  row 5 (the ECS angle x 1/2, through n = 9 and 10), 0.66 % to row 10, 0.33 % at row 11;
  phase pi/2 and m = 1 take phase 0's m = 0 at the same row.

## Files

```
e55524b29cbae31ac669a56da8c9875085448bf3c9df8b0622ac15530772e929  scan_wallace/results_cep0.csv
d256699367123c27ec4dc620fe7b183104ad51736fe907e65df30f41440364ff  scan_wallace/results_cep90.csv
d262087dffb8267a9bf0ee46008356085e393916f4a2b2f434a00474436f0d77  make_inputs.py
b4ed1f5f64ac9afcc56acd6eb1009e7fe4091224f5408eaf14f1022b976d1626  extract.py
6de6efdb3f82d585d92429a53ed370130d88efdf81eb93f409723a43c1ea6ec0  tables.py
d1f34b01bf47f8bf4abfc429faea14852206b44239a08f07c32c674b08b30ac0  tables_cep0.md
d8383e27a739a4ec3cfbf7748943cdadc58215d3b969232e1054cf76e4063bce  tables_cep90.md
de8a59da23940854b77c20f1d03c6b9c51a511502c076fc86debf476e8f54425  REPORT.md
```

# Atomic hydrogen, 400 and 800 nm, 15 fs (October 2026)

H 1s in the bare Coulomb potential, a Gaussian of 15 fs (full width at half maximum of
the intensity) on the vector potential at 400 and 800 nm, carrier-envelope phase 0,
ten peak intensities from 5e13 to 3e14 W/cm2, then 1000 au without field; BOX 80, ECS
angle 0.15. The datasets `ionised_gauss15_400nm_H` and `_800nm_H` are written from
`scan_hydrogen/results_400.csv` and `results_800.csv` by
`tools/extract/travers_2026_trecx.py` (TC-843). Tier 3.

## The code and the runs

tRecX `a987d7c790e52b6d7f44a5615a665e0abc3e7f4a`, the build of the argon runs; 42 runs
per wavelength (ten intensities, convergence at rows 0, 4 and 9 with the argon list
less the ECS angle x 1/2, the base angle being 0.15 already; a run without a field),
every one complete. The step is set by the spectral cut (0.016 to 0.019 au).

## Checked (against the tables and the runs' report, 2026-10-03)

- Initial <H0> -0.5000000 au; the pulse's intensity FWHM 15.00002 fs (400 nm) and
  15.000004 fs (800 nm), peak field to 1e-7.
- The drift is 2.6e-9 (400 nm) and 2.1e-9 (800 nm) in 4721 au; every yield is far
  above 10 times it.
- The uncertainty by `six_fs_uncertainty`: at most 0.058 %.
- Hydrogen's levels are degenerate in l, so the states are labelled by n alone, which
  is what 1 - bound(n <= 10) needs.

## Files

```
6461780326b00737e3651e4ecf20afc91e8c3bcfaea00affd59b875ec5fd747a  scan_hydrogen/results_400.csv
429b9f8d00bc3d10136342b067b99948a7f6230db19d474dd8ad992ef275b941  scan_hydrogen/results_800.csv
92e83e1b5e8106464a743256a3a3d7720538d5891b46d817c47c740b5fcc8db1  make_inputs.py
8ef3504671b962da351c0a0fb163a731498c212e500b7483445873622345fa4a  extract.py
69cd0666edc8b810067644e12883fb57f5f80df4e82f1991b0aec6aee6070b11  tables.py
7ba306fd394d6dfd6e879bc1f1a7d04f989cf270b8b4f8cce49894a6a8939a43  tables_h400.md
4d63d3f5094b72d32e7f2b5a26406362ddbef6a4bdf9adb37043d456c01a72f0  tables_h800.md
94e3171a681741660ccf36a90f03c0a5632b4f6f65a7aa9be9abee903ad9729d  REPORT.md
```

# Argon 3p at 400 and 1030 nm with a 30 fs pulse (October 2026)

The 400 and 1030 nm set-ups and intensity grids of the 10 fs runs with a 30 fs Gaussian
(full width at half maximum of the intensity) on the vector potential, then 1000 au
without field; window [-3720.72, 4720.72] au. The datasets
`ionised_gauss30_400nm_Ar_m0`, `_m1` and `ionised_gauss30_1030nm_Ar_m0`, `_m1` are
written from `scan_wavelength/results_400_30fs.csv` and `results_1030_30fs.csv` by
`tools/extract/travers_2026_trecx.py` (TC-844). Tier 3.

## The code and the runs

tRecX `a987d7c790e52b6d7f44a5615a665e0abc3e7f4a`; 53 runs at 400 nm and 59 at 1030 nm
with the convergence rows and lists of the 10 fs runs, every one complete (32 min to 9
h each).

## Checked (against the tables and the runs' report, 2026-10-03)

- Initial <H0> -0.5796893340 au; the pulse's intensity FWHM 29.9988 fs (400 nm) and
  29.99999 fs (1030 nm), peak field to 5e-11.
- The drift is 3.4e-9 in 8441 au. At 1030 nm rows 0 and 1 of m = 0 and 0 to 2 of m = 1
  are within a factor 10 of it, and excluded.
- The uncertainty by `six_fs_uncertainty`: 0.34 % to row 3, 0.18 % to row 8, 0.15 % at
  row 9 (400 nm; cutEnergy, one-sided); 4.1 % at row 2, 3.9 % at rows 3 and 4, 0.42 %
  to row 8, 0.14 % at row 9 (1030 nm); m = 1 takes m = 0's.

## Files

```
bd468d3f4b2dc90538aec7768acb8d47d92761d55891f053c696e1fdd3cff31a  scan_wavelength/results_1030_30fs.csv
1fb0ceee6b40ac7859f1a34a39fc967cbcbd5a68c1200532034c9d4d4d2c8155  scan_wavelength/results_400_30fs.csv
3fe028cd5133b357a0d99b343545805e345affa8ef9b4c47f4e3c11553bc06ae  make_inputs.py
4ebbde09bb5cffaa6539a22b2e68cd9a5e5d2ea7047ac5b4d64e040b2da65b17  extract.py
9a1545efc51b461fb75197e76b75acf391065258706c6d057a23bd9ba7e7ffb2  tables.py
e6e4dff85c0d9cfff13112c9fc98b34150441026ccebb678983389b4dffe9fd6  tables_400_30fs.md
a3bc623b0032ced4a5f59398d30d29f974babd9f8e7e154ab206981f3231a5c4  tables_1030_30fs.md
63ec305f4cc292aa85c42b46548de5750f9bbdc82b36c49745ca0a47874909b3  REPORT.md
```

# Krypton and xenon, 800 nm, 6 fs, and their static rates (October 2026)

Kr 4p and Xe 5p (m = 0, one active electron) in zhang_2014's potentials (Kr: Green,
Sellin and Zachor's form, Table II; Xe: Tong and Lin's form, Table I), their valence
energies not adjusted (Kr 14.0370 eV, +0.27 % against the measured; Xe 12.1040 eV,
-0.21 %). A2's 6 fs pulse file and ten intensities; a radial grid refined over A2's
inside 20 bohr (6 elements of order 16 on [0, 2], 8 of order 16 on [2, 20]); the step
capped by the deep 1s (0.00136 au Xe, 0.00346 au Kr). Static rates by part C's trace
of the complex-scaled resonance (BOX 40, angle 0.4) at 0.03 to 0.10 au. The datasets
`ionised_gauss6_Kr_m0`, `_Xe_m0`, `static_Kr_m0`, `static_Xe_m0` are written from
`scan_kr_xe/results.csv` and `static_results.csv` by
`tools/extract/travers_2026_trecx.py` (TC-845). Tier 3.

## Checked (against the tables and the runs' report, 2026-10-03)

- The potentials' parameters against zhang_2014's table images; every run's initial
  <H0> at the field-free value to 5e-11.
- The runs without a field lose nothing to 12 digits: every row is above the floor.
- The uncertainty by `six_fs_uncertainty` (with A2's radial grid as a variation): 1.3 %
  (Xe) and 1.5 % (Kr) to row 3 (cutEnergy 150 at row 0), 0.31 % and 0.46 % to row 8,
  0.11 % at row 9; the static rates' by `kr_xe_static_uncertainty`: at most 0.3 %,
  except 1.2 % for Xe at 0.10 au (over the barrier; the ECS angle).
- All 78 scan runs, 78 spectra, 44 label spectra and 90 static traces complete, exit 0;
  the only warning tRecX's memory advisory.

## Files

```
118c76616dfd9e1f1ee3989a39c25cfaddce44160c614b4115762dcb9c23fb75  scan_kr_xe/results.csv
8efd8ebb750e2c4a70b83d2b87512e3a9ea4892c1d510093b82a30e72d7a34aa  scan_kr_xe/static_results.csv
efa354a2b65d27434d439d6e0632a306087269ab45192d9486b006bfe1d684d8  make_inputs.py
e89c441bfd2c5cfba7a62d272aaf0699b13201a1ea7d893b1dd50bbdcf5f5b26  make_static.py
79eb69638dec801b4ad3b657e77077969011fa21b25703fafe395d057aae0f0d  extract.py
5a11dc7aaca23f5b13276aa6fafe66b80ddcc400b7b7299553190b6b006a3c7a  extract_static.py
3a8b78f2f689cf42c341905bf53d87637151c06b505e4c6e4680ea8c1dbac4ee  tables.py
dca029b2ea178e19833115d0cd5f7ee4a89af420bff9add6608be4f89d4e09f5  spectra.py
7833a51f3958d4c50dd85f26ac5b9e8ae1b5d5d094a60a73e066aa998e7ff7fd  tables_g.md
c7ffae61cac763ae469816148f00f9f8a5e885e0550d9a0070d773aea50c5b07  static_tables.md
f1dfa195dcd922e6dcb4a3e7d86dbf41e57ab0ad9b56883686163217c151078a  REPORT.md
```

# Krypton 4p in Wallace et al.'s pulse (October 2026)

Kr 4p, m = 0, in part G's potential (zhang_2014 Table II) with its refined radial grid
and step cap (0.00346 au), the argon Wallace-pulse set-up (5.5 fs, phases 0 and pi/2,
twelve intensities from 5e13 to 6e14 W/cm2); BOX 80, and 120 at rows 9 to 11. The
datasets `ionised_wallace_cep0_Kr_m0` and `_cep90_Kr_m0` are written from
`scan_wallace/results_kr_cep0.csv` and `results_kr_cep90.csv` by
`tools/extract/travers_2026_trecx.py` (TC-846). Tier 3.

## Checked (against the tables and the runs' report, 2026-10-04)

- Initial <H0> -0.5158494127 au (-0.5158455 on A2's grid); no drift to 12 digits.
- The largest norm 1.000000; 108 runs (54 and their spectra), all exit 0.
- The uncertainty by `wallace_uncertainty` with these tables: 0.70 % at row 0, 0.73 %
  to row 10, 0.27 % at row 11; phase pi/2 takes phase 0's.
- Warnings: one GMRES tolerance warning (`kr_m0_11_th23`, its yield within 0.012 % of
  the base) and the known eigen-decomposition check at BOX 120 (`kr_m0_11_a2grid`).

## Files

```
e6790bb62a93c404d988ae68bb58456d222fb6a89fb1e5a9a40a70a7d9819207  scan_wallace/results_kr_cep0.csv
8d934dcbb7913a1345996b9f9a4ab62e8164b64788a9e52317970cdf7411a033  scan_wallace/results_kr_cep90.csv
d262087dffb8267a9bf0ee46008356085e393916f4a2b2f434a00474436f0d77  make_inputs.py
3482ca16d48d7f87f1c69213669fc39311e0e728f1788c47be74f5b56f025356  extract.py
3a8b78f2f689cf42c341905bf53d87637151c06b505e4c6e4680ea8c1dbac4ee  tables.py
1cce7089d321ee63da8ccf7d42883ce4eb5ed2102a4a2f62cf113a66e5ea78d1  tables_kr_cep0.md
e2d8b9bf5d1d3c06cc1ea5968071669f182fbeb3cf4238af540168912e389110  tables_kr_cep90.md
dbd6518cddee1418459db060c11c59511ea2cbc36b1a79404afa31f0604d6eca  REPORT.md
```

# The argon atom in haCC, 800 nm, 6 fs (October 2026)

tRecX's hybrid anti-symmetrised coupled channels (haCC) for argon with all electrons:
the Scrinzi group's `Ar_qz` data set (repository `chemical`, commit
88c17cfe12c150199bcfad92d6909267edebea21, 2025-02-01, "added Ar_qz, quadruple zeta
basis"; one neutral and four ion configuration-interaction states; its 3p ionisation
energy 15.25 eV, 3.2 % below the measured; no spin-orbit coupling). The code's history
shows the data set added the same day as a "preliminary" XUV tutorial (tRecX commit
71d8b6be, `tutorial/203ArHaCC`), so it is not validated by its author for strong
infrared fields. A2's 6 fs pulse at 800 nm at 4.03e13, 9.61e13 and 2.29e14 W/cm2; the
yield is the norm-loss 1000 au after the pulse (haCC cannot project on the neutral's
bound states; its argon iSurff integral is 250 to 290 times the norm-loss in every run,
a property of the data set, and is not used). The datasets
`hacc_normloss_gauss6_Ar_c4` and `_c3` are written from `ar_hacc/results.csv` by
`tools/extract/travers_2026_trecx.py` (TC-847). Tier 3.

## The code and the runs

tRecX `a987d7c790e52b6d7f44a5615a665e0abc3e7f4a`, MPI build (16 processes; the decisive
3-channel mMax 1 run on a whole node). The table holds every run of E1 (helium) and
E2 (argon), including the MPI checks whose one- and four-channel runs at 9.61e13 W/cm2
are complete over A2's window and are the curated ones at that intensity.

## Checked (against the table and the runs' report, 2026-10-04)

- Every curated run complete over A2's window (t_last 1744.14 au); the initial energies
  give the set-ups' effective ionisation energies (4 channels 0.5631718 au, 3 channels
  0.5626988, 1 channel 0.5608551).
- The uncertainty by `hacc_uncertainty`: the one-channel variations at the same
  intensity (accuracy, rMax 120, radial order 19 and 12, lMax 30 and 40) and the
  set-up's own lMax 30 run: 1.9 %, 7.3 % and 1.2 % (lMax 20 is high; radial order 24 is
  a failure of haCC and is left out).
- **Not in this table:** the 3-channel mMax 1 run (`e2_9p61e13_c3m1`) finished after
  the table was written (4 Oct); the runs' report gives its norm-loss, 2.2203e-3 at
  +1000 au. It is to be added when the table is re-extracted.

## Files

```
47f1177c15fe1fd3a6b4ef4e54672b470c5560970328c3867e9f7307e9775e3e  ar_hacc/results.csv
c33c3f0636ce75bfe8affcdb2bc5b4470691e2ac35e5fc96b9f636fedce38bc3  E2/extract.py
67426923ebc242fb65a0ab083d336e0dbef4619bfa9ec621c92ce282883ea984  E2/make_inputs.py
b20b18319fc411a04adf8b5a34755b02bd7cd2d055c1bc31ed416a2e152d4d7f  REPORT.md
f4caa391c9b8f10d794312a56b097cbf0be22e115b226fd8af08aa3317ca792e  chemical_Ar_qz_He19_sha256.txt
```

# Xenon 5p in Wallace et al.'s pulse (October 2026)

Xe 5p, m = 0, in part G's potential (zhang_2014 Table I) with its refined radial grid
and step cap (0.00136 au), the argon Wallace-pulse set-up (5.5 fs, phases 0 and pi/2,
twelve intensities from 5e13 to 6e14 W/cm2); BOX 80, and 120 at rows 9 to 11. The
datasets `ionised_wallace_cep0_Xe_m0` and `_cep90_Xe_m0` are written from
`scan_wallace/results_xe_cep0.csv` and `results_xe_cep90.csv` (copies of the runs'
`runs_xe_cep0/results.csv` and `runs_xe_cep90/results.csv`) by
`tools/extract/travers_2026_trecx.py` (TC-849). Tier 3.

## Checked (against the tables and the runs' report, 2026-10-06)

- Initial <H0> -0.4448128882 au (-0.4447651734 on A2's grid); no drift to 12 digits.
- The largest norm 1 + 2e-11; 108 runs (54 and their spectra), all exit 0.
- 23 phase-0 runs (rows 0 to 6, the convergence runs at rows 0 and 6, three at row 11)
  were stopped by the job's 2-day limit 35 to 99.7 % through and continued from their
  checkpoints (`run_scan_pack_restart.slurm`); the column `restarted` marks them. The
  norm continues across the restart (xe_m0_00: 0.986203, 0.986189, 0.986178 from
  t = 820.03 au). Six of them lack the expectation row at the kill time, all after the
  pulse (792 to 986 au); the measures are read at the end of the pulse and at +1000 au.
  The runs' `extract.py` was changed for them (the initial energy from the first
  pass's output, the wall times summed); the phase-0 file was extracted with it.
- The uncertainty by `wallace_uncertainty` with these tables: 0.84 % at row 0, 1.2 %
  to row 10 (cutEnergy 150 at row 6), 0.03 % at row 11; phase pi/2 takes phase 0's.
- The Keldysh parameter is 1.42 to 0.41 with the potential's I_p (the report's prose
  gives 1.41 at row 0, its table 1.42).
- Warnings: one GMRES tolerance warning (`xe_m0_11_th23`, its yield within 0.015 % of
  the base) and the known eigen-decomposition check at BOX 120 (`xe_m0_11_a2grid`).

## Files

```
8230ec5b51e99d723fc116bf434a3c710856596c378d73045a82d886627d95c9  scan_wallace/results_xe_cep0.csv
48d592fd59bd48ba56feae79de59da54e0eeb9f4d1968faa502fe69650116569  scan_wallace/results_xe_cep90.csv
d262087dffb8267a9bf0ee46008356085e393916f4a2b2f434a00474436f0d77  make_inputs.py
be9128de20078786347a10c811057406a2ee49017cabbc9d071d5c93b3fe9530  extract.py
3a8b78f2f689cf42c341905bf53d87637151c06b505e4c6e4680ea8c1dbac4ee  tables.py
b0ee978be1dafaef2a1ff1a74392f933fcd44ff17fbb8adf9a0e5658c325c9cf  run_scan_pack_restart.slurm
8b1ccd46ff92a2b44feabb1b4fe1bd3ff98b110d76abe50557e36ad80a38fc22  tables_xe_cep0.md
27a35370f7d0b7db8421163219cedba6020d06d3b794c69a167a5e937de2f987  tables_xe_cep90.md
3c4dfae73b457b5dcbb3dbf58796f5ebc52ee77870c927c6eebf3ae3fba38e94  REPORT.md
```

# Part M: one electron in the sin^2 pulses of published many-electron calculations (October 2026)

Ar 3p (Tong and Lin's potential, and Pathak et al. 2020's Muller-Kooiman form), Kr 4p
(part G's) and He 1s, m = 0 and m = 1, in E(t) = E0 sin(omega t) sin^2(pi t / tau): argon
at 1200 nm, 3 cycles, 1-4e14 W/cm2 (M1, Pathak et al. 2020); krypton at 800 nm, 3
cycles, 1.6-2.5e14 (M2, Pathak et al. 2022); argon at 800 nm, 3 cycles, 3-6e14 (M3, Sato
et al. 2018); helium at 800 nm, 6 cycles, 4-14e14 (M4, Sato et al. 2016). The datasets
`outside20_sin2_*` and `ionised_sin2_*` are written from `scan_sin2/results_*.csv`
(copies of the runs' `runs_<part>/results.csv`) by `tools/extract/travers_2026_trecx.py`
(TC-858). Tier 3.

## Checked (against the tables and the runs' report, 2026-10-06)

- The pulses: the field as run within 2.9e-9 E0 of the formula, A(tau) <= 2.2e-10 A0;
  peak |E| / E0 = 0.936525 (3 cycles) and 0.983190 (6 cycles).
- The measure 1 - <chi[0,20]>, tRecX's characteristic function of r < 20 au on the
  unscaled region (an element boundary of A2's grid at 20 au): BOX + 40 changes it by at
  most 0.001 %, where the norm-loss at tau moves by 16 to 48 %.
- Convergence at the middle intensity of each part: at most 0.36 % (argon, krypton),
  1.2 % (helium, angular functions 30 to 50); the uncertainties recomputed by
  `six_fs_uncertainty` (m = 1 takes the same atom's m = 0 changes).
- Muller-Kooiman's 3p at -15.820 eV, as Pathak et al. state; drift below 1.2e-8; the
  largest norm 1.00000000; the known warnings only (the eigen-decomposition check at
  BOX >= 120; one GMRES tolerance warning in a convergence run, within 0.003 %).

## Files

```
9b4a62a52617b8671e46119463c141aeb6c7ef73720dba6896daf7f9ad3f6bdf  scan_sin2/results_m1_ar.csv
d966b6cb1b47f41dc9e1592dc0f42b639b10efcdc4aa37d32b0c121d38780536  scan_sin2/results_m1_mk.csv
e0463f7cf085fd6aaa13c8d9a5829e9b77f6be7ce9ad593a2f24a46c5217d3de  scan_sin2/results_m2_kr.csv
7450aa4703a654da62cf9bb66392fa702fae87b15acbe582db660442f64680a9  scan_sin2/results_m3_ar.csv
176b92990aab940d28f393449e117881a4f59991df28d5b3e26217eb673cdbef  scan_sin2/results_m4_he.csv
1d414a7adf95fdfdd166cbf1d6057124faf06446372c0440038c72f0dee4198d  make_inputs.py
7b10eeb93fb72fef49865dcf65ba8f652130145fd749ad7a75d838cff4428d77  extract.py
99a39c90bbdf742b383f9d11323d562e7d1795e96571880f6bca94aeb1065167  tables.py
1fdad234c7d1e018d0a0f448b440e597c92a3157716de1c7f9bc65aa4c1fe455  summary.py
d1d13876d8aa94caf4c03fb79c197407baf879f998341207c642c9870f4662ca  label_states.py
60eb8e18983b67fc6804d23abf64f2dca56d4d2765c8f8ab4e14029639c8d742  summary.md
9afcce75daf54bdc4f3ee2686f4246b1fe1dc196071f4e08dfe3c37b888659c3  REPORT.md
4b181d44cd290800d5cf7da443f1eddbee4b7e4eab35a5ed71b6161dde792f07  pulse_check_m1_ar.txt
0cc39154abd6e1df52298f5fd3062dfaf59cf8f0158f1392a635c0202c525e54  pulse_check_m1_mk.txt
64205b5d0da028429ab26835952ffb9b2c327648c64ab4d049413ebc425ee781  pulse_check_m2_kr.txt
118613f45b3809828d0f3f37cb4d2a61d2b58acf40c1eacb793fdc110b2c800c  pulse_check_m3_ar.txt
35310d60893378d76e3036ddaa67d32f502994e12c7a9a3ff3cce9ce58cf1e8f  pulse_check_m4_he.txt
```

# The argon atom at 20 fs, and haCC's three-channel mMax 1 (E3; October 2026)

`ar_hacc/results.csv` refreshed from the runs' `E2/results.csv` of 6 October 2026
(`E2/extract.py`, sha256 below; the pulse window now read from each input, every older
row identical to the last digit): `e2_9p61e13_c3m1` complete (norm-loss at +1000 au
2.2203e-3, 74.7 h), and the E3 runs `e3_9p61e13_c3_20fs` (8.4373e-3) and
`e3_9p61e13_c3l30_20fs` (8.1776e-3): three ion channels, mMax 0, lMax 20 and 30, a 20 fs
Gaussian on the vector potential at 800 nm (7.50 cycles in the full width at half
maximum of the intensity, phase 0, window [-3T, +3T + 1000 au]), 9.605215e13 W/cm2.
`ip_matched/results_c3_6fs.csv` and `results_c3_20fs.csv`: one electron in Tong and
Lin's argon scaled by beta = 0.9939109898 to the set-up's 0.5626988 au, m = 0 and m = 1,
in A2's 6 fs pulse and the 20 fs one (the runs' `E3/sae/runs_c3_*/results.csv`).
`hacc_normloss_gauss20_Ar_c3`, `hacc_normloss_gauss6_Ar_c3m1`, `normloss_` and
`ionised_ipmatched_c3_*`, `_c3_20fs_*` (TC-861). Tier 3.

## Checked (against the tables and the runs' report, 2026-10-06)

- haCC's initial energies -0.56269883987679 and -0.56269883987676 au, as E2's
  three-channel runs; the largest norm 1; exit 0; the inverse-overlap warning at 5.7e-9
  and 4.4e-9, as in every argon run.
- The pulse file: full width 19.998 fs on its grid, peak |E| / E0 0.999998.
- One electron at 20 fs: every convergence variation within 0.25 % (ionised fraction)
  and 0.42 % (norm-loss); drift 2.0e-9 over the window. The 6 fs set has no
  convergence runs of its own and takes those of the c3m1 set (the same pulse and
  discretisation).
- The `.git` of `testruns/` does not exist: the extraction is identified by the
  sha256 of `E2/extract.py`.

## Files

```
180986cccda64e627e899592c85a1809f3e130521fb846b175a1866c9d694264  ar_hacc/results.csv
0321b716ff7dfadaf64cd311ed9917481a68265894c05328a02ba6bad20d9a1e  ip_matched/results_c3_6fs.csv
7e7bb97fa55bb60e66045249e0b3935e276a3a6954799b7d3b3d21f317f94962  ip_matched/results_c3_20fs.csv
ae655b0ba45b614711c79825cad061bfd716efb8131e1f35d15db2711e95470f  E2/extract.py
457db7a6ab2c9091d02b721de3f139313c6008297ffd32e4366bbd8c92c352eb  E3/make_inputs.py
b29942fa2023de552507d17d7960afb774504c53f42fe9da02f1585d1623c521  E3/pulse_check.txt
23f32e79a20ff9e4714a6b87bc67ce1ae31810cd57af08e4f8fbe25ba041b3d0  E3/sae/make_inputs.py
4ebbde09bb5cffaa6539a22b2e68cd9a5e5d2ea7047ac5b4d64e040b2da65b17  E3/sae/extract.py
e09d9b5f0f48d08e62fda70247e6051cbc32bb146c9d5a3ec323ac96f86da7ff  E3/sae/beta_result.txt
386fa0ee2e37aff9ec2f7873024cb1a1a96abd92a2b96f39bd88c961dce49de3  E3/sae/tables_c3_20fs.md
37f76a4fbeda68a0f42381813244447d04949ab7a83c6974481556a222735330  E3/sae/tables_c3_6fs.md
b114c64989b196f47e3237650f77d9945bfcfbac143a7393caf428207df222e3  REPORT.md
```

# Helium with two active electrons at 390 nm (part K, step 2; October 2026)

Helium with both electrons active in tRecX (MPI build a987d7c, 16 ranks; total M = 0,
|m| <= 1, l <= 15 per electron with Lshape[3;5], BOX 40; and the variations), haCC's
helium (He19_Standard, one ion channel, lMax 40) and one electron of Tong and Lin's
helium at three ionisation energies, in A2's 6 fs Gaussian at 390 nm (4.61 cycles), at
5.68981e14 and 3.684031e14 W/cm2; the norm-loss 1000 au after the pulse. The tables
are copies of the runs' `runs2/results.csv` (`extract_runs2.py`, run on the cluster on 11
October 2026), `hacc2_results.csv` and `sae390_{base,m5,tl}/runs_*/results.csv`. The
datasets `twoelectron_normloss_gauss6_390nm_He*`, `hacc_normloss_gauss6_390nm_He`,
`normloss_` and `ionised_ipmatched_390nm_*_He` are written from them by
`tools/extract/travers_2026_trecx.py` (TC-1050). Tier 3.

## Checked (against the tables and the runs' report, 2026-10-11)

- Two electrons: every run exit 0; the initial energy -2.9018542200 au (the |m| <= 2
  basis -2.9031255440, BOX 60 -2.9018540962); norms never above 1; drift 3.9e-10 per au;
  the only warning the chi[0,20] quadrature note (no element boundary at 20 au; the
  norm-loss is unaffected). BOX 60's values from its own log (its queued continuation
  rewrote the expectation file without propagating). The values agree with the report's
  section 8 to the last digit.
- haCC: initial energy -0.90186662 au; the 'Large spectral constraint of 10 %' warning
  (haCC fixes cutEnergy 400). The five-channel mMax 1 run (lMax 20) finished its
  propagation (in two jobs, the second from the checkpoint) but was cancelled in its
  iSurff pass. Its row was refreshed from its expectation file on 11 October 18:10
  (`refresh_hacc2_c5m1.py`; the previous table dbded4a6… kept with the runs): done
  True, the norm-loss 6.5847e-3 at +1000 au, +2.08 % on one channel at lMax 20, as the
  report; the initial energy agrees with the expectation file to 1e-12. The other
  rows are unchanged; the refresh adds a `note` column, empty but for this row.
- One electron: beta 0.9969113371 and 0.9989292176 for the two-electron bases'
  ionisation energies, 1 for Tong and Lin's; A2's convergence list at 5.69e14 for the
  base (norm-loss BOX 80 -4.95 %, the ionised fraction <= 1.7 %).

## Files

```
04f8136d632e8e7157ca55118b3570dd020fd8a1ad74e108d9db81a014659a1c  he_two_electron/results_2e.csv
9bea326a18ddce4aa76b9147e8b5e55a30f14b5bf71a376f1ad4d94181422e52  he_two_electron/results_hacc.csv
c581d3a1fc956b7c399db4b5f987d3f70afadb1f2a538acdee4002c3ce3e847c  he_two_electron/results_sae_base.csv
9b6b9274826c2781b443f606f5c45418741c41b8779c3a7ca1cd537d50474bac  he_two_electron/results_sae_m5.csv
f5a2a1c9712bfc394983946a16acf8dcaca49d6a715bd1a48cead5c843a41aab  he_two_electron/results_sae_tl.csv
b1526ff8ac0d3cbdaecb4ef77b8f878b67ce2ac98e7139038375c32c7a74f411  extract_runs2.py
8630104f31523035998a4d99c8495f0775a015c5d4e39450e76cb566b90fd8f0  step2.py
3e16cb8a00bb8c10c25574f445609750787a11e1806dfda6f77beb704b612e77  pulse_check.txt
f7bab4e823f649d6d762711256e1e69206abc31d1db3b3fff10d088ed4d93240  PARAMETERS_two_electron.md
67f7251ac7065c3db3f05b6bf3c9ae5d9476ec3f1b7532c309cc9e1800eb57f1  REPORT.md
```

# Part N: one electron of xenon in Hartmann et al.'s RMT pulse (October 2026)

Xenon 5p (m = 0 and m = 1) in part G's potential (Tong and Lin's form with Zhang, Lan and
Lu's parameters; part G's refined inner grid, Eta 30, BOX 80, irECS angle 0.3 from BOX,
cutEnergy 100, fixStep 0.00136 au), its short-range terms scaled by beta =
0.9997362280 so that the 5p energy is RMT's -0.444000 au (Hartmann et al. 2022's
`GS_Energy_desired`; |dE| 4.5e-13), and unscaled (-0.4448129 au) at three intensities.
The pulse is the one RMT wrote (`hartmann_2022`): E(t) = E0 exp(-((t - 900)/186)^2)
cos(0.06 (t - 900)), t in au, with E0 = 0.05336 sqrt(I / 1e14 W/cm2) au, as A(t) =
-int_0^t E dt' on [0, 3000] au from a file (velocity gauge), then 1000 au without field.
tRecX a987d7c (OccDual build), one core per run; 41 runs and their 41 field-free
spectra, made on 12 to 15 October 2026. `results.csv` and `tables.md` are the runs'
own tables (the runs' `extract.py`, kept with the runs), `manifest.csv` the run list,
`inputs/` every input as run (each run and its `_eig` spectrum). The datasets
`ground_survival_rmtpulse_*_Xe` are written from `results.csv` by
`tools/extract/travers_2026_trecx.py` (TC-1093). Tier 3.

## Checked (against the tables and the runs' report, 2026-10-16)

- Every run exit 0; initial energy -0.44400000013 au (scaled) and -0.4448128882 au
  (unscaled), the 5p state at index 11 (m = 0) and 6 (m = 1); OccDual 1 in the 5p at
  t = 0.
- The pulse file: -dA/dt peaks at 0.99999994 E0 at t = 900.00 au and reproduces the
  formula to 6.1e-8 E0; A(3000) = -1.2e-10 E0; the intensity's full width at half
  maximum 219.00 au (the runs' `pulse_check.txt`, kept with the runs).
- The ground-state population at 1500 au equals that at 3000 au to 1e-5 relative in
  every run. Convergence at 1.5 and 2.5e14 W/cm2, m = 0 (Eta 40 and 50, 20/14 radial
  elements, order 14, BOX 120, the step x 2/3, irECS angle 0.2, operatorThreshold 1e-11,
  A2's inner grid, cutEnergy 150): the depletion 1 - P changes by at most 0.43 %
  (cutEnergy 150 at 1.5e14). Without field over [0, 4000] au, P = 1 to every printed
  digit and the norm-loss <= 2e-12. The one warning ('Large number of 600 radial
  spectral points') concerns the iSurff spectrum's memory.
- The values agree with the report's tables to the last printed digit.

## Files

```
84d3a368bb69da48e20f952e34cdcd071d634382f01cf50ec9d7bcb70dc093dc  xe_rmt/results.csv
b5310e3da7981c55ff65f497f77e2fbd048c1491a1fafc4fde266f5f7b4e71a2  xe_rmt/manifest.csv
1de17b10adbf602b2ddc9551ad4bef474fbe034a7d9018eb374c67c136fb04e4  xe_rmt/tables.md
3cff4fa4daa7c5d1c143fd4d0dd644799edcd80616d79401c9f512f62ae44942  extract.py
a5f3de325efe57d53cc7569bd697a39ea3fb9a73d8fb0025512ce6dca8ae0982  make_inputs.py
f3245c35017435ffb6c6efe426fea8d7536613220a94033fa4dfe74a53231690  beta_search_xe.py
6c1582b15b91ace5a48390ae99ddb96990c41fc7b1b7e2ceb27c8d5dec814fc2  pulse_check.txt
05b9ebeb6f705bf0f26b4bdde83f47a1cd232c6c8ad78cf5dbc32ebde6c8a57e  REPORT.md
```

`results.csv`, `tables.md` and `inputs/` are as the runs left them. In `manifest.csv`
the runs' directory is replaced by the placeholder `$RUN_ROOT`; as run, its sha256 was
`a71cc3db8a0c60294c351d8a9158544a731fe5e9edf8a02a87bcdeb26b2e5e8d`. The last five files
(the runs' scripts, pulse check and report) are kept with the runs, outside this
repository. The inputs name the pulse file by its path on the cluster.

## Terms

These are our own calculations, run with the tRecX code (not included), and are
distributed under optmat's MIT licence.
