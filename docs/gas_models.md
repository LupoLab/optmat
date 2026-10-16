# The gas models

Read {doc}`models` for selection and {doc}`fitting` for construction. This chapter
records the supporting evidence and coefficients for each gas; {doc}`accuracy`
interprets the uncertainty. Start with {doc}`quickstart` for a first calculation,
and use {doc}`/the-data/sources` to inspect the cited records.

Every gas has one fitted model of its refractive index (`fit/<gas>`), recommended by
default, and each published formula of the curated sources beside it for comparison
({doc}`models`). The fit and its model are described on {doc}`fitting`. This page says
what each model rests on, how far the formulas in common use are from it, and what the
curation found about the data along the way.

## What each model rests on

The *backbone* of a model is the data that set its level and its dispersion: ab initio
theory, precision metrology or a dataset of stated accuracy. The rest either constrain
the shape without their own level ("level free": the dataset's scale is fitted) or are
validators, compared with the fit but not fitted. The range is the vacuum-wavelength
range over which the model is declared valid. The uncertainty is the fit's statistical
one-sigma uncertainty of $n - 1$ (the level) at 800 nm; {doc}`accuracy` gives the
uncertainties of the dispersion, which matter more for propagation, with the working
estimates.

<!-- gas-models:start -->
| Gas | What the model rests on | Range (nm) | $n - 1$, $1\sigma$ at 800 nm | Reduced $\chi^2$ |
|---|---|---|---|---|
| He | Puchalski et al. 2016's ab initio dispersion coefficients with every relativistic and QED correction; Pachucki & Puchalski's retardation term; Schmidt et al. 2007's microwave static polarizability. The experiments are tested against it, not needed by it. | 60-static | 0.1 ppm | 0.58 |
| Ne | Lesiuk et al. 2020's ab initio moments; Gaiser & Fellmuth's dielectric-constant gas thermometry; the dispersion of Cuthbertson (level free), Bideau-Mehu and Shelton's He/Ne ratio. | 76-static | 8.6 ppm | 2.63 |
| Ar | Lesiuk & Jeziorski 2023's ab initio moments; dielectric-constant gas thermometry; Egan et al.'s refractometry; the shapes of Peck & Fisher (level free), Bideau-Mehu and Gordon. | 108-static | 2 ppm | 0.50 |
| Kr | Dielectric static values (Orcutt & Cole, Huot & Bose); Birch and Hohm & Kerl at 633 nm; Burns et al. at 458-633 nm; Smith et al. 1976 in the ultraviolet; Bideau-Mehu; Shelton's Ar/Kr ratio; Hohm's 325-633 nm dispersion. | 125-static | 78 ppm | 1.23 |
| Xe | Huot & Bose's dielectric static value; Burns et al.; Hohm & Kerl at 633 nm; Abjean et al. 1971 and Bideau-Mehu to 147 nm; Shelton's Ar/Xe ratio; Hohm's 325-633 nm dispersion. | 148-static | 1.7e+02 ppm | 0.96 |
| H2 | Raj et al. 2018's ab initio average of normal hydrogen (0.1 % allowed); the shapes of Peck & Hung and Smith et al.; Shelton's He/H2 ratio. | 112-static | 1.3e+02 ppm | 0.22 |
| D2 | Raj et al. 2018's ab initio average of normal deuterium; Shelton's D2/H2 and H2/D2 ratios. | 112-static | 1.4e+02 ppm | 0.14 |
| N2 | Peck & Khanna; the 633 nm metrology (Birch, Stone, Achtermann, Egan); static values; Griesmann's vacuum ultraviolet (level free); Shelton's He/N2 ratio. | 145-static | 15 ppm | 1.16 |
| O2 | Egan et al. 2026 at 633 and 1542 nm; Birch; Smith et al. 1976 above the Schumann-Runge bands; Hohm's 325-633 nm dispersion; Shelton's N2/O2 ratio. | 195-static | 5.4 ppm | 0.23 |
| CO | Smith et al.'s 167-288 nm interferometry; the Cuthbertsons' visible dispersion; Watson's green line; Zhao et al.'s 1566 nm refractometry; HITRAN's bands. | 167-2288 + 2 intervals to 12500 | 42 ppm | 0.37 |
| CO2 | Old, Gentili & Peck; the Cuthbertsons and Watson & Ramaswamy (level free); the 633 nm metrology; Egan et al. 2026; Hohm's dispersion; Shelton's ratios; Marchetti's 10 um values as validators; HITRAN's bands; Chan et al.'s two weak bands below 10 eV as lines. | 250-1400 + 7 intervals to 8941 | 8.1 ppm | 4.01 |
| CH4 | The Cuthbertsons' and Watson & Ramaswamy's dispersion (level free); Achtermann; Burns et al.; Zhao et al.'s 1566 nm value; Shelton's 14 ratios; Hohm's 325-633 nm dispersion (level free); HITRAN's bands. | 240-1123 + 4 intervals to 12088 | 44 ppm | 0.68 |
| SF6 | Swift; Vukovic et al.; Achtermann; Burns et al.; Shelton's 8 ratios; HITRAN's bands. | 230-5000 | 3.3e+02 ppm | 4.07 |
| N2O | Hohm's 325-633 nm dispersion and the Cuthbertsons' 1914 dispersion (both level free); Birch; Egan et al.; HITRAN's bands; the Thomas-Reiche-Kuhn sum as a prior. | 320-1516 + 3 intervals to 7048 | 58 ppm | 3.01 |
| O3 | The measured absorption spectrum by Kramers-Kronig (Serdyuchenko et al. 2014 at 293 K and three others); one fitted strength above 12 eV; Swift et al., the Cuthbertsons' level and the ab initio static value with Naves et al.'s vibrational correction. | 185-3800 | 0.96 % | 0.61 |
<!-- gas-models:end -->

**Range:**
- **The atoms** extend to the static limit. An atom has no infrared resonance, so
  between the near infrared and zero frequency its polarizability is fixed by its
  dipole sums $S(-2)$ and $S(-4)$, which the backbones determine. The mid-infrared
  points of Rusch (argon, 6.76 µm) and Marchetti (helium, 10.57 µm) agree within their
  few-per-mille uncertainties.
- **The infrared-active molecules** (CO, CO2, CH4, N2O) hold on several intervals
  (OEP 3; {doc}`fitting`). Their bands are Gaussian rules of HITRAN's lines.
  - Their optical range is split only where the pure gas absorbs. CH4 now ends at
    1313 nm and has gaps at 1123-1137, 1149-1177 and 1187.4-1187.9 nm; N2O has gaps at
    1516-1518 nm and from 1950 nm.
  - Beyond the optical range they hold again in the windows between the bands, where
    the bands' representation meets the tolerance. CO2's last window is 8.4-8.9 µm: its
    laser bands at 9.4 and 10.4 µm are kept, so 8.9-12.5 µm, where Marchetti's values
    lie, is not a qualified interval (TC-797).
  - The range warning names each gap's reason. The column gives the first interval,
    and how many follow and where the last ends.
- **Hydrogen and deuterium** rest on theory; the experiments' levels scatter by about
  1e-3 about it (TC-95).
- **Ozone** is modelled differently. Its measured absorption enters by Kramers–Kronig,
  so the model is complex, and its uncertainty is set by the structure of the
  absorption spectrum rather than by the fit (see {doc}`accuracy`).

**Temperature:** the molecules' polarizabilities refer to 293.15 K (N2, O2, CO2, N2O,
H2, D2). Egan et al.'s 303 K values that set the levels of CO2 and O2 are brought there
with their own measured coefficients, which moves the models' $n - 1$ by -40 ppm (CO2)
and -21 ppm (O2). CH4 and CO rest
on 297.5 K refractometry and SF6 mostly on 373 K values, whose temperature
coefficients are not curated. The atoms need no reference temperature (see the
conventions page, TC-730).

## The molecules' qualified intervals

CO, CO2, CH4 and N2O carry their infrared bands as Gaussian rules of HITRAN's lines, and
are qualified on several intervals (OEP 3; {doc}`fitting`, {doc}`data`). For each
interval:
- the largest smoothed absorption of the pure gas at 101325 Pa and 296 K in it, at a
  resolution of 1 cm⁻¹;
- the largest single-line Lorentz peak, so that a user can judge their own path length,
  pressure and bandwidth;
- the reason for the gap that follows it, as the derivation labels the gap over its
  whole span.

Inside a fit's optical range only absorption splits it (TC-749). There the gap begins
where the absorption reaches 0.1 m⁻¹, whatever the gap's label. Beyond the range, the
representation of the bands decides (`representation`; `kept band` where a band is kept
as poles rather than folded), and `width` marks a window narrower than 50 cm⁻¹.

The models are the background at that resolution, not the resolved structure of the
weak lines. A picosecond pulse near a weak band sees structure the model does not hold.
The representation bound is a deterministic bound, given beside the statistical sigma
and not combined with it (`FittedModel.representation_bound_au`).

<!-- molecule-intervals:start -->
| Gas | Interval (µm) | Largest $a_\sigma$ (m⁻¹) | Largest line peak (m⁻¹) | The gap that follows |
|---|---|---:|---:|---|
| CO | 0.1673-2.288 | 0.023 | 0.27 | representation: kept band |
| CO | 2.565-4.34 | 2.6e-09 | 4.7e-08 | absorption |
| CO | 5.191-12.5 | 3.5e-05 | 0.00032 | — |
| CO2 | 0.2504-1.4 | 0.0086 | 0.019 | representation |
| CO2 | 2.302-2.38 | 1.2e-06 | 1.8e-06 | representation |
| CO2 | 2.381-2.475 | 1.3e-05 | 3.5e-05 | representation |
| CO2 | 2.541-2.59 | 0.00024 | 0.0003 | representation |
| CO2 | 3.203-3.35 | 3.5e-05 | 0.00012 | representation |
| CO2 | 3.357-3.593 | 0.00011 | 0.00023 | representation |
| CO2 | 5.549-6.946 | 1e-05 | 9.5e-06 | representation: kept band |
| CO2 | 8.435-8.941 | 6.7e-07 | 1.6e-06 | — |
| CH4 | 0.2401-1.123 | 0.09 | 0.21 | width |
| CH4 | 1.137-1.149 | 0.1 | 0.39 | absorption |
| CH4 | 1.177-1.187 | 0.099 | 0.2 | absorption |
| CH4 | 1.188-1.313 | 0.097 | 0.13 | representation |
| CH4 | 11.12-12.09 | 2e-06 | 5e-06 | — |
| N2O | 0.3202-1.516 | 0.098 | 0.028 | absorption |
| N2O | 1.518-1.95 | 0.098 | 0.18 | representation |
| N2O | 4.814-4.956 | 0.0013 | 0.003 | representation |
| N2O | 6.752-7.048 | 0.0056 | 0.0098 | — |
<!-- molecule-intervals:end -->

Pure CO2 is not qualified at its laser lines: its 9.4 and 10.4 µm bands are kept, and
the last window ends at 8.94 µm.

The intervals are what [OEP 3](oeps/archive/oep-0003-infrared-windows.md) proposed and
implemented: its theory, decisions and phase reports are the record.

## The compact coefficients

Each molecule's compact model ({doc}`fitting`, *Compact infrared terms*; OEP 2) as
printed coefficients, from `optmat.gas.models.compact.coefficients` and
`format_coefficients`, to eleven significant digits: the rounding then changes the
polarizability, $D_2$ and $D_3$ by at most 0.002 of the fit's sigma. The polarizability
form goes through the package's density conversion (Lorentz-Lorenz, the real-gas
density and the refractivity virial; {doc}`conventions`), which keeps the stated
tolerance at any density. The Sellmeier form is $n - 1$ in the dilute limit at
101325 Pa and 293.15 K, the temperature the molecular polarizabilities refer to, with
its error against the full model stated: about $10^{-4}$ of $n - 1$, mostly the
Lorentz-Lorenz term it drops. Each table holds only in its intervals; between and
beyond them the gas absorbs, and the compact poles are effective ones, not band
centres.

<!-- compact-coefficients:start -->
`compact/CO`, the polarizability form; qualified at 167.3-2000.0 nm. $\alpha = \sum_k f_k / (E_k^2 - \omega^2)$ in atomic units, $E_k = hc / \lambda_k$.

| Term | Kind | $f_k$ | $\lambda_k$ (µm) |
|---|---|---:|---:|
| effective 2 | effective | 4.2653378498e+00 | 7.3643744199e-02 |
| effective 1 | effective | 1.8512815886e-01 | 1.4715780853e-01 |
| compact 3 | infrared | 8.2279160419e-08 | 2.3394546499e+00 |
| compact 2 | infrared | 1.6094654462e-08 | 2.7915109723e+00 |
| compact 1 | infrared | 1.1641060935e-05 | 4.6649107831e+00 |
| rotational 0.0 | rotational | 2.0923184033e-08 | 2.5412095195e+02 |

`compact/CO`, the Sellmeier form at 101325 Pa and 293.15 K ($N$ = 2.5045083675e+25 m$^{-3}$); qualified at 167.3-2000.0 nm. $n - 1 = \sum_k B_k \lambda^2 / (\lambda^2 - C_k)$; against the full model, at most 8.7e-05 relative in $n - 1$ and 0.00021 of the largest $\beta_2$.

| Term | Kind | $B_k$ | $C_k$ (µm²) |
|---|---|---:|---:|
| effective 2 | effective | 2.5983614722e-04 | 5.4234010597e-03 |
| effective 1 | effective | 4.5031199478e-05 | 2.1655420612e-02 |
| compact 3 | infrared | 5.0581711831e-09 | 5.4730480591e+00 |
| compact 2 | infrared | 1.4087526484e-09 | 7.7925335087e+00 |
| compact 1 | infrared | 2.8454676784e-06 | 2.1761392615e+01 |
| rotational 0.0 | rotational | 1.5176902074e-05 | 6.4577458222e+04 |

`compact/CO2`, the polarizability form; qualified at 250.4-1400.0 nm. $\alpha = \sum_k f_k / (E_k^2 - \omega^2)$ in atomic units, $E_k = hc / \lambda_k$.

| Term | Kind | $f_k$ | $\lambda_k$ (µm) |
|---|---|---:|---:|
| effective 2 | effective | 6.4907086934e+00 | 6.4889749029e-02 |
| effective 1 | effective | 7.7236988164e-01 | 1.0667109688e-01 |
| CO2 1Pi_g (1pi_g -> 3s sigma_g), about 9.30 eV | line | 7.6027031505e-03 | 1.3332000000e-01 |
| CO2 1Delta_u (1pi_g -> pi*), about 8.38 eV | line | 5.9884165655e-03 | 1.4795000000e-01 |
| compact 5 | infrared | 2.1372528019e-09 | 1.4319392695e+00 |
| compact 4 | infrared | 8.2792836701e-10 | 1.4334495406e+00 |
| compact 3 | infrared | 8.0538606278e-07 | 2.2958630257e+00 |
| compact 2 | infrared | 1.1273466434e-04 | 4.1495064626e+00 |
| compact 1 | infrared | 1.4322012500e-05 | 1.3101376123e+01 |
| rotational 0.0 | rotational | 2.1672907678e-16 | 6.5691891555e+02 |

`compact/CO2`, the Sellmeier form at 101325 Pa and 293.15 K ($N$ = 2.5169069720e+25 m$^{-3}$); qualified at 250.4-1400.0 nm. $n - 1 = \sum_k B_k \lambda^2 / (\lambda^2 - C_k)$; against the full model, at most 8.8e-05 relative in $n - 1$ and 0.00017 of the largest $\beta_2$.

| Term | Kind | $B_k$ | $C_k$ (µm²) |
|---|---|---:|---:|
| effective 2 | effective | 3.0850584292e-04 | 4.2106795291e-03 |
| effective 1 | effective | 9.9205999884e-05 | 1.1378722910e-02 |
| CO2 1Pi_g (1pi_g -> 3s sigma_g), about 9.30 eV | line | 1.5253788847e-06 | 1.7774222401e-02 |
| CO2 1Delta_u (1pi_g -> pi*), about 8.38 eV | line | 1.4796569170e-06 | 2.1889202501e-02 |
| compact 5 | infrared | 4.9467980443e-11 | 2.0504500716e+00 |
| compact 4 | infrared | 1.9203335482e-11 | 2.0547775854e+00 |
| compact 3 | infrared | 4.7919816272e-08 | 5.2709870326e+00 |
| compact 2 | infrared | 2.1911366006e-05 | 1.7218403883e+01 |
| compact 1 | infrared | 2.7749615169e-05 | 1.7164605631e+02 |
| rotational 0.0 | rotational | 1.0557468691e-12 | 4.3154246161e+05 |

`compact/CH4`, the polarizability form; qualified at 240.1-1070.0 nm. $\alpha = \sum_k f_k / (E_k^2 - \omega^2)$ in atomic units, $E_k = hc / \lambda_k$.

| Term | Kind | $f_k$ | $\lambda_k$ (µm) |
|---|---|---:|---:|
| effective 2 | effective | 4.1267565855e+03 | 1.5360712677e-03 |
| effective 1 | effective | 2.5306220540e+00 | 1.0154430067e-01 |
| compact 5 | infrared | 6.8182074197e-09 | 1.1249545641e+00 |
| compact 4 | infrared | 2.1158851236e-08 | 1.1576252434e+00 |
| compact 3 | infrared | 7.4914852172e-08 | 1.3729279715e+00 |
| compact 2 | infrared | 1.0832835868e-05 | 2.8784429688e+00 |
| compact 1 | infrared | 9.2912506578e-06 | 6.6302255555e+00 |
| rotational 0.0 | rotational | 4.8380075997e-11 | 1.0171149383e+02 |

`compact/CH4`, the Sellmeier form at 101325 Pa and 293.15 K ($N$ = 2.5081503359e+25 m$^{-3}$); qualified at 240.1-1070.0 nm. $n - 1 = \sum_k B_k \lambda^2 / (\lambda^2 - C_k)$; against the full model, at most 0.00012 relative in $n - 1$ and 0.00021 of the largest $\beta_2$.

| Term | Kind | $B_k$ | $C_k$ (µm²) |
|---|---|---:|---:|
| effective 2 | effective | 1.0953100472e-04 | 2.3595149395e-06 |
| effective 1 | effective | 2.9352417863e-04 | 1.0311244999e-02 |
| compact 5 | infrared | 9.7061200818e-11 | 1.2655227714e+00 |
| compact 4 | infrared | 3.1895803645e-10 | 1.3400962041e+00 |
| compact 3 | infrared | 1.5884330706e-09 | 1.8849312149e+00 |
| compact 2 | infrared | 1.0096315379e-06 | 8.2854339244e+00 |
| compact 1 | infrared | 4.5944789867e-06 | 4.3959890917e+01 |
| rotational 0.0 | rotational | 5.6300480167e-09 | 1.0345227977e+04 |

`compact/N2O`, the polarizability form; qualified at 320.2-1500.0 nm. $\alpha = \sum_k f_k / (E_k^2 - \omega^2)$ in atomic units, $E_k = hc / \lambda_k$.

| Term | Kind | $f_k$ | $\lambda_k$ (µm) |
|---|---|---:|---:|
| effective 2 | effective | 2.0888813178e+01 | 3.3472221509e-02 |
| effective 1 | effective | 1.1188454137e+00 | 1.2507109255e-01 |
| compact 5 | infrared | 7.5438533852e-10 | 1.5159040405e+00 |
| compact 4 | infrared | 3.1374048071e-09 | 1.5174692101e+00 |
| compact 3 | infrared | 4.6537983265e-07 | 2.2819960154e+00 |
| compact 2 | infrared | 6.3167428291e-05 | 4.2594630608e+00 |
| compact 1 | infrared | 1.8062801923e-05 | 8.5288772204e+00 |
| rotational 0.0 | rotational | 1.0534861604e-08 | 5.3412336540e+02 |

`compact/N2O`, the Sellmeier form at 101325 Pa and 293.15 K ($N$ = 2.5179981372e+25 m$^{-3}$); qualified at 320.2-1500.0 nm. $n - 1 = \sum_k B_k \lambda^2 / (\lambda^2 - C_k)$; against the full model, at most 8.3e-05 relative in $n - 1$ and 0.00018 of the largest $\beta_2$.

| Term | Kind | $B_k$ | $C_k$ (µm²) |
|---|---|---:|---:|
| effective 2 | effective | 2.6429574227e-04 | 1.1203896127e-03 |
| effective 1 | effective | 1.9764748866e-04 | 1.5642778191e-02 |
| compact 5 | infrared | 1.9576902942e-11 | 2.2979650600e+00 |
| compact 4 | infrared | 8.1586379407e-11 | 2.3027128035e+00 |
| compact 3 | infrared | 2.7368100121e-08 | 5.2075058145e+00 |
| compact 2 | infrared | 1.2942263447e-05 | 1.8143025567e+01 |
| compact 1 | infrared | 1.4838027503e-05 | 7.2741746641e+01 |
| rotational 0.0 | rotational | 3.3940584367e-05 | 2.8528776946e+05 |

`compact/SF6`, the polarizability form; qualified at 230.0-5000.0 nm. $\alpha = \sum_k f_k / (E_k^2 - \omega^2)$ in atomic units, $E_k = hc / \lambda_k$.

| Term | Kind | $f_k$ | $\lambda_k$ (µm) |
|---|---|---:|---:|
| effective 2 | effective | 1.7687344964e+01 | 5.9307116398e-02 |
| effective 1 | effective | 5.8154314433e-03 | 1.8225339665e-01 |
| compact (the parent's) 3 | infrared | 3.0068868641e-06 | 6.3390105109e+00 |
| compact (the parent's) 2 | infrared | 2.2162382890e-04 | 1.0581503679e+01 |
| compact (the parent's) 1 | infrared | 1.5200997054e-05 | 1.6244212056e+01 |

`compact/SF6`, the Sellmeier form at 101325 Pa and 293.15 K ($N$ = 2.5338300116e+25 m$^{-3}$); qualified at 230.0-5000.0 nm. $n - 1 = \sum_k B_k \lambda^2 / (\lambda^2 - C_k)$; against the full model, at most 0.00021 relative in $n - 1$ and 0.00026 of the largest $\beta_2$.

| Term | Kind | $B_k$ | $C_k$ (µm²) |
|---|---|---:|---:|
| effective 2 | effective | 7.0697761300e-04 | 3.5173340554e-03 |
| effective 1 | effective | 2.1951419072e-06 | 3.3216300590e-02 |
| compact (the parent's) 3 | infrared | 1.3730597305e-06 | 4.0183054258e+01 |
| compact (the parent's) 2 | infrared | 2.8199449285e-04 | 1.1196822012e+02 |
| compact (the parent's) 1 | infrared | 4.5582566469e-05 | 2.6387442531e+02 |
<!-- compact-coefficients:end -->

## Compared with the formulas in common use

- **Helium:** the largest difference for propagation is its dispersion. Several published or legacy helium
  formulas give a group-velocity dispersion in the visible and near
  infrared 5–14 % away from the ab initio value. The fit, Stone & Stejskal's formula and
  Puchalski et al.'s own series agree to four digits.
- **Argon and krypton:** the common formulas are good in the visible, but their
  vacuum-ultraviolet extensions (unpublished fits distributed with propagation codes)
  err by 5–15 % near the resonances, where resonant dispersive waves are phase-matched.
- **CH4, SF6 and N2O:** the published formulas (He et al. 2021) have no infrared term.
  - They miss Shelton's near-infrared phase-matching ratios by 3 % (CH4) and 28 % (SF6)
    at 1319 nm.
  - They overstate SF6's group-velocity dispersion at 800 nm by 55 %.

The {doc}`accuracy` page gives the statistical and empirical uncertainties. The
{ref}`dispersion gallery <sphx_glr_auto_gallery_accuracy_plot_gvd_comparison.py>`
compares named formulas and fits over their stated ranges.

## What the curation found

- **A reference state misread.** Ermolov et al.'s helium formula is a fit at 1 atm,
  although the paper says 1 bar. At 1 bar it lies 1.3 % above every source (TC-5).
- **Hohm's absolute levels.** His 1994 polarizabilities (Mol. Phys. 81) lie 0.04-0.47 %
  above the fits (Kr +0.20 %, Xe +0.47 %), and his own 1990 values of the same gases are
  0.1-0.3 % lower. His shapes agree: argon's with theory to about 4e-5 rms (7e-5 at
  325 nm). The fits take his dispersion with the level free (TC-400).
- **The classical absolute levels are off:**
  - argon: Peck & Fisher +0.035 %, Zhang et al. +0.026 % (TC-73);
  - the light gases: the old refractometry is high by 0.1–0.3 %, consistent with traces
    of heavier gases (TC-90, TC-93).

  The fits take these datasets' shapes, not their levels.
- **Mansfield & Peck's helium dispersion** departs from theory by 1e-4 (TC-89).
- **Neon's ab initio dispersion** is 2.5 sigma below every experiment (TC-92). This is
  the one open physics question of the noble-gas models.
- **The magnetic term and the virial term matter at the 1e-5 level**, and bring the
  metrology into agreement: Stone's helium refractivity at 633 nm closes to 0.8 ppm, and
  Egan et al.'s argon value to 20 ppm of Lesiuk's theory.
- **Misprints**, each corrected with its evidence (107 errata and 13 corrected formula
  coefficients in all):
  - Bideau-Mehu's neon and xenon formulas;
  - Börzsönyi's xenon $C_1$ (a factor of 1000);
  - Mansfield & Peck's theoretical formulas;
  - many table values.

  Bideau-Mehu's krypton formula does not reproduce its own table and is kept as
  printed.
- **Tensions in the molecular fits.** The SF6, CO2 and N2O fits carry reduced
  chi-squared of 4.1, 3.3 and 3.0:
  - SF6's two sets of Shelton ratios prefer HITRAN's band strengths times 1.036 (TC-146);
  - N2O's two visible dispersions, Hohm's and the Cuthbertsons', differ in shape by about
    1e-4, and the fit follows Hohm's (TC-178).

## Cross-checks

- **Shelton's phase-matching ratios.** These density ratios tie the near-infrared
  dispersion of pairs of gases together to 0.02–0.3 %. The fits use them, with the
  partner gas held at its recommended model, and every fitted ratio lies within 1.9
  sigma.
- **Argon, the one inconsistency.** Shelton's set is self-consistent, but the argon
  model's He/Ar ratio is 0.3 % below his (TC-96). Fitting his N2/Ar ratio would pull
  nitrogen onto that difference, so it is a validator.
- **Börzsönyi et al.'s group-delay dispersions at 800 nm.** The formulas in common use
  were fitted to these. They scatter about the fits by −12 % (neon) to +21 % (helium,
  within its 16 % uncertainty), and by +10 % for nitrogen (12 sigma); they are
  validators (TC-97).

The gallery page *Cross-checks* shows both.

## Robustness

Each fit is refitted without each of its datasets (gallery page *Robustness*):
- **The backbones** move their models by many sigma when dropped, by design:
  - helium's ab initio moments;
  - neon's and argon's dielectric values and moments;
  - hydrogen's and deuterium's theory.
- **Every other dataset** (100 in all): 80 move their model by less than 2 sigma
  anywhere in range and 93 by less than 4. The seven above 4 sigma are each the only
  data over part of a model's range, so the shift appears where the model would
  otherwise be extrapolated:
  - the ultraviolet of CO (Smith et al. 1976, 63 sigma), N2O (Hohm's 325 nm point, 34),
    O2 and Kr (Smith et al. 1976, 13 and 5.4) and CH4 (Shelton's 1986 harmonics, 5.8);
  - CO's near infrared (Zhao et al.'s 1566 nm value, 7.7);
  - N2O's level (Egan et al. 2019, 5.6).
- **Hydrogen's level** once rested on one 633 nm metrology value. It is now on theory,
  and three ab initio calculations agree to 0.15 % (TC-95).
- **Oxygen's level**, which rested on Birch alone, is now fixed as well by Egan et al.'s
  2026 refractometry, in agreement (TC-147).
- **Near-degenerate minima.** A far effective oscillator enters only through $f/E^2$,
  so the fits' cost has near-degenerate minima. The ones a wider search finds differ
  from the stored fits by at most 0.44 sigma anywhere in range. CO2's two effective
  oscillators were the exception: with the lower one free down to 202 nm, the minimum
  followed its bound and the ultraviolet dispersion moved by up to 30 % at 250 nm
  (TC-725). Its two weak bands below 10 eV are now lines with priors from the measured
  photoabsorption, and the effective oscillators start at 10 eV (TC-726).

The compact coefficients are what
[OEP 2](oeps/archive/oep-0002-compact-infrared.md) proposed and implemented: its theory,
decisions and phase reports are the record.
