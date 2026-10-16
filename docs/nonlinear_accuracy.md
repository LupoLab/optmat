# Third-order accuracy

Read {doc}`/nonlinear` and {doc}`/user-guide/nonlinear-support` first. This chapter
compares statistical uncertainty, dataset sensitivity and independent measurements.
{doc}`/user-guide/nonlinear-extrapolation` treats evaluation beyond fitted data, and
the {ref}`uncertainty gallery <sphx_glr_auto_gallery_nonlinear_plot_nonlinear_accuracy.py>`
plots the distinction between statistical and working estimates.

## What the numbers mean

The Kerr gamma sets the nonlinear index, $n_2 = N\gamma_\text{Kerr}/(8\varepsilon_0^2 c)$,
and so everything that scales with it. A relative error $\varepsilon$ in the Kerr gamma is
an error $\varepsilon$ in:
- the nonlinear phase, the B-integral and the spectral broadening by self-phase
  modulation;
- the critical power of self-focusing (which scales as $1/n_2$);
- the square of the soliton order ($N^2 \propto n_2$), so $\varepsilon/2$ in the order;
- the nonlinear phase mismatch of every phase-matched process driven by the pulse
  (dispersive waves, modulational instability, four-wave mixing).

For the third harmonic the generated field scales with the THG gamma, so its power scales
as $\gamma_\text{THG}^2$: an error $\varepsilon$ in gamma is $2\varepsilon$ in the
harmonic's power.

As with the linear index, a level error in gamma is equivalent to a density error, since
every third-order effect scales as $N\gamma$. The linear index is known to 1e-4 or better,
far below the uncertainty of a gas fill. Gamma is not: the best absolute scale (the ESHG
scale, below) is good to about 1 %, and older values in use are off by 3-30 %. Here the
level matters as much as the shape.

The shape is the dependence on the wavelength and on the process. Between a fundamental
of 1800 and 400 nm the Kerr gamma rises by 12 % (He) to 86 % (Xe). At 800 nm the THG gamma
exceeds the Kerr one by 6 % (He) to 39 % (Xe). A single constant coefficient cannot give both.

## Three measures and the working estimate

Each value below comes with three measures of its uncertainty:

- **Statistical:** the one-sigma uncertainty propagated from the fit's covariance, scaled
  by the reduced chi-squared where it exceeds one
  ({meth}`~optmat.gas.nonlinear.HyperpolarizabilityModel.gamma_sigma_au`). It assumes the
  dispersion form and the stated uncertainties.
- **Leave-one-out:** the largest relative change when one fitted dataset is left out and
  the gas refitted (`tools/fit/nonlinear_robustness.py`, records in
  `data/gas/models/nonlinear_robustness/`). A calculation's processes leave together (they
  share one scale). A refit with fewer values than parameters is not counted: there the
  priors alone decide and the change measures nothing. An exactly determined refit is
  counted.
- **Floor:** a systematic uncertainty that neither of the other two contains, stated per
  gas:
  - **He, 0.1 %:** the allowance added to each of Bishop & Pipin's calculated values.
    Shelton & Rice judge the calculation better than 0.1 %, and it is nonrelativistic.
  - **Ne, Ar, Kr, Xe, H2, D2, N2, O2, 1 %:** the ESHG scale. Every one of Shelton's
    values rests on Bishop & Pipin's helium through measured ratios. Shelton states the
    ratios to about 0.4 % (0.1-0.7 %) and the curves of Ar, H2 and N2 to about 0.2 %.
    He then adds that helium pair interactions "may increase chi^(3) ... by 1%". All his
    values "should be adjusted upward by about 1%" (Ne half as much), "the effective
    accuracy of our results is only about 1%, even for the atoms and diatomic molecules"
    (Shelton 1990, p. 13). The adjustment is not applied.
  - **CO2, CH4, SF6, 2 %:** the same scale, and the fit is to the total ESHG gamma, whose
    vibrational part Shelton could not separate for these molecules.
  - **CO, 5 %:** the level rests on one ESHG value, Ward & Miller's 1730 +- 50 au (2.9 %,
    electronic plus vibrational).
  - **O3, 12 %:** the allowance for the method dependence of Naves et al.'s CCSD(T)
    values.
  - **N2O, none beyond the statistical:** the n2's stated uncertainty and the priors on
    its borrowed dispersion are already in the statistical figure.

The **working estimate** is the largest of the three, rounded up to the next of 0.1,
0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20 %. It is mechanical, and the test
recomputes it; the judgement is in the floors and in the notes below.

## At 800 nm

The Kerr gamma (SPM, XPM) and the THG gamma at a fundamental of 800 nm. "Largest when
left out" names the dataset, or the calculation (its scale group), whose removal moves
the Kerr gamma most.

<!-- nonlinear-accuracy-800:start -->
| Gas | Kerr: statistical | Kerr: leave-one-out | Largest when left out | THG: statistical | THG: leave-one-out | Floor | Kerr estimate | THG estimate |
|---|---|---|---|---|---|---|---|---|
| He | 0.015 % | 0.012 % | `bishop_1989/omega_thg_zzzz` | 0.036 % | 0.011 % | 0.1 % | 0.1 % | 0.1 % |
| Ne | 0.23 % | 0.19 % | `shelton_1992/gamma_present` | 0.23 % | 0.25 % | 1 % | 1 % | 1 % |
| Ar | 0.15 % | 0.31 % | `shelton_1990/table3_Ar_literature` | 0.16 % | 0.4 % | 1 % | 1 % | 1 % |
| Kr | 0.58 % | 3.1 % | `shelton_1990/table3_Kr_present` | 0.86 % | 0.19 % | 1 % | 4 % | 1 % |
| Xe | 0.39 % | 1.4 % | `shelton_1990/table3_Xe_present` | 0.71 % | 3.1 % | 1 % | 1.5 % | 4 % |
| H2 | 0.1 % | 0.15 % | `bishop_1991` | 0.13 % | 0.14 % | 1 % | 1 % | 1 % |
| D2 | 0.085 % | 0.1 % | `bishop_1991` | 0.12 % | 0.054 % | 1 % | 1 % | 1 % |
| N2 | 0.17 % | 0.47 % | `hattig_ccsd` | 0.14 % | 0.68 % | 1 % | 1 % | 1 % |
| O2 | 0.24 % | 1.6 % | `shelton_1990/table3_O2_present` | 0.39 % | 0.21 % | 1 % | 2 % | 1 % |
| CO2 | 0.51 % | 2.6 % | `shelton_1990/table3_CO2_present` | 0.36 % | 3.1 % | 2 % | 3 % | 4 % |
| CH4 | 0.3 % | 0.74 % | `shelton_1990/table3_CH4_present` | 0.28 % | 0.092 % | 2 % | 2 % | 2 % |
| SF6 | 0.31 % | 3.2 % | `shelton_1990/table3_SF6_present` | 0.27 % | 0.34 % | 2 % | 4 % | 2 % |
| CO | 2.8 % | 10 % | `ward_1979/chi3_CO` | 2.8 % | 10 % | 5 % | 10 % | 10 % |
| N2O | 13 % | none | none | 14 % | none | 0 % | 15 % | 15 % |
| O3 | 10 % | 0.082 % | `naves_2011/gamma_static_ccsd_t_expt` | n/a | n/a | 12 % | 15 % | n/a (resonant) |
<!-- nonlinear-accuracy-800:end -->

## At 400 and 1800 nm

The Kerr gamma at the two ends of the range most used (n/a where outside the model's
declared range):

<!-- nonlinear-accuracy-bands:start -->
| Gas | 400 nm: statistical | leave-one-out | estimate | 1800 nm: statistical | leave-one-out | estimate |
|---|---|---|---|---|---|---|
| He | 0.12 % | 0.046 % | 0.2 % | 0.017 % | 0.0031 % | 0.1 % |
| Ne | 0.24 % | 0.23 % | 1 % | 0.27 % | 0.63 % | 1 % |
| Ar | 0.44 % | 0.7 % | 1 % | 0.23 % | 0.37 % | 1 % |
| Kr | 3.3 % | 0.003 % | 4 % | 1 % | 6.3 % | 7 % |
| Xe | 2.9 % | 6.8 % | 7 % | 0.42 % | 2.1 % | 3 % |
| H2 | 0.54 % | 0.4 % | 1 % | 0.11 % | 0.11 % | 1 % |
| D2 | 0.53 % | 0.4 % | 1 % | 0.093 % | 0.17 % | 1 % |
| N2 | 0.41 % | 0.93 % | 1 % | 0.33 % | 0.41 % | 1 % |
| O2 | 1.5 % | 0.23 % | 1.5 % | 0.5 % | 3.7 % | 4 % |
| CO2 | 1.2 % | 7.1 % | 10 % | 1.1 % | 5.9 % | 7 % |
| CH4 | 1.2 % | 0.026 % | 2 % | 0.64 % | 1.7 % | 2 % |
| SF6 | 0.78 % | 0.25 % | 2 % | 0.63 % | 6.7 % | 7 % |
| CO | 3 % | 10 % | 10 % | 3 % | 10 % | 10 % |
| N2O | 16 % | none | 20 % | 14 % | none | 15 % |
| O3 | n/a | n/a | n/a | 7.3 % | 0.1 % | 15 % |
<!-- nonlinear-accuracy-bands:end -->

At 400 nm the Kerr gamma is at the edge of its data. The largest fitted $W_2$ corresponds
to a Kerr fundamental of 374-401 nm for the gases with ESHG data (263 nm for helium's
calculation). For N2O and O3 it corresponds to 664 and 695 nm, so N2O's 400 nm value
rests on its borrowed dispersion. The single legs of the data reach much further than
the Kerr diagonal: to the 229-244 nm harmonics of the visible ESHG and calculated THG
(152 nm for helium; {doc}`/user-guide/nonlinear-support`).
At 1800 nm the Kerr gamma is within a few per cent of the static limit, fixed by
Shelton's 1064 and 1319 nm values. Leaving those two out, and extrapolating the visible
data, moves it by 3.7 % (O2), 5.9 % (CO2), 6.3 % (Kr) and 6.7 % (SF6). This is why
those estimates are large.

<span id="beyond-the-fitted-data-bounds-on-the-legs"></span>
<span id="on-the-data-and-within-the-fitted-w2"></span>
<span id="against-a-carrier"></span>
<span id="the-polynomial-beyond-the-data"></span>

See {doc}`/user-guide/nonlinear-extrapolation`.

## Gas by gas

- **He:** on theory. The fit reproduces Bishop & Pipin's four processes with residuals of
  about 0.015 % rms.
  Their tabulated values are in the reduced-mass atomic unit of 4He, and are converted to
  standard atomic units (+0.096 %) since 2026-09-26, as Shelton & Rice do.
- **Ne:** Shelton's 1990 and 1992 ESHG values agree. Hättig & Jørgensen's CCSD process
  dependence enters at a fitted scale of 1.030 (their values 3 % above the fit).
- **Ar:** Shelton's 19 rescaled literature values, his two infrared values and Hättig &
  Jørgensen (scale 1.011) agree (chi-squared 24 for 36).
- **Kr:** only Shelton's six values (two infrared, four rescaled literature), which
  disagree: chi-squared 20 for 2 degrees of freedom. The statistical uncertainty is
  scaled up by 3.2, and leaving out the infrared pair moves the Kerr gamma by 3.1 % at
  800 nm and 6.3 % at 1800 nm. There is no calculated process dependence (B' is held at
  its prior, 0).
- **Xe:** Tarazkar et al.'s CCSD calculation (scale 1.009) gives the process dependence.
  Leaving out the rescaled literature moves the Kerr gamma at 400 nm by 6.8 % and the THG
  at 800 nm by 3.1 %.
- **H2, D2:** Bishop et al.'s 1991 electronic calculation agrees with Shelton's
  electronic ESHG (fitted scales 0.997 and 0.995). The vibrational and rotational part
  is subtracted with Shelton's own tabulated values.
- **N2:** Hättig & Jørgensen (scale 0.994). The vibrational and rotational part is
  subtracted with Shelton's formula (Table V).
- **O2:** Shelton only, with B' held at 0. The infrared pair carries the long-wavelength
  side, hence 4 % at 1800 nm.
- **CO2:** Sekino & Bartlett's CCSD(T) (scale 0.940: their values 6 % below the fit).
  Leaving out the rescaled literature moves the Kerr gamma at 400 nm by 7.1 %.
- **CH4:** Shelton only; his 1319 nm point is excluded, as Shelton excluded it. Hättig &
  Jørgensen's CCSD lies 11-15 % below and is a validator; fitting its shape was tried and
  rejected (TC-314).
- **SF6:** Shelton only, with B' held at 0.
- **CO:** the level is Ward & Miller's ESHG value. Sekino & Bartlett's CCSD(T) gives the
  process dependence at a scale of 1.113, so their values are 11 % above the measurement.
  Without the measurement the level would be theirs, 10 % higher.
- **N2O:** one value, Wahlstrand et al.'s 2012 n2 at 800 nm (electronic plus
  vibrational), with CO2's A and B as priors of 50 % and 100 %. Nothing else of any kind
  is in the library.
- **O3:** Naves et al.'s CCSD(T) static and Kerr values, 700 nm and longer; the third
  harmonic of 800 nm falls in the Hartley band and is not given.

## What the independent measurements say

The fits use the ESHG scale and, for the process dependence, calculations. Every other
measurement is a validator, compared with the fit at its own process and wavelengths.
Each ratio below is the measured gamma over the fitted model, averaged over the
dataset's values; conversions follow {doc}`/user-guide/nonlinear-conventions`.

<!-- nonlinear-validators:start -->
| Validator | Source | Data / model (mean over the values) |
|---|---|---|
| Wahlstrand et al. 2012 (PRA), SSSI n2, 800 nm | `wahlstrand_2012` | Ar 0.89; N2 0.87; O2 1.03 |
| Wahlstrand et al. 2012 (PRL), SSSI n2, 800 nm | `wahlstrand_2012b` | He 0.84; Ne 0.94; Ar 0.92; Kr 0.92; Xe 0.89 |
| Wahlstrand et al. 2015, SSSI electronic n2, 800 nm | `wahlstrand_2015` | H2 1.01; D2 0.93 |
| Zahedpour et al. 2015, SSSI n2, 400-2400 nm | `zahedpour_2015` | Ar 0.95; N2 (n2_N2) 0.96; N2 (n2_adjusted_N2) 0.96; O2 0.95 |
| Zahedpour et al. 2019, SSSI n2, 3-10.5 um | `zahedpour_2019` | Ar 0.88; N2 0.99; O2 1.02 |
| Wang et al. 2013, SPM in a hollow fibre | `wang_2013` | Ne 1.02; Ar (n2_Ar_800) 1.05; Ar (n2_Ar_1800) 0.68 |
| Nibbering et al. 1997, 'instantaneous' n2, 800 nm | `nibbering_1997` | Ar 1.39; Xe 1.30; N2 2.91; O2 6.08; SF6 1.33 |
| Börzsönyi et al. 2010, n2 (header corrected), 800 nm | `borzsonyi_2010` | Ne 1.98; Ar 1.88; Xe 1.32; N2 0.83 |
| Loriot et al. 2009 (with erratum), n2, 800 nm | `loriot_2009` | Ar 0.97; N2 1.36; O2 1.86 |
| Geints et al. 2025, self-focusing, 800 nm | `geints_2025` | CO2 0.97 |
| Pigeon et al. 2018, four-wave mixing, 10.6 um | `pigeon_2018` | Kr 1.34; Xe 0.88 |
| Lehmeier et al. 1985, THG, 1055 nm | `lehmeier_1985` | He 0.96; Ne 0.69; Ar 0.77; Kr 0.90; Xe 0.94; N2 0.88 |
| Ward & New 1969, THG, 694 nm | `ward_1969` | Ne 0.86; Ar 0.93; Kr 1.13; Xe 0.91; H2 0.97; N2 1.03; CO2 1.11 |
| Ward & Elliott 1984, THG, 694 nm | `ward_1984` | CH4 0.92; SF6 0.98 |
| Ward & Miller 1979, ESHG, 694 nm | `ward_1979` | H2 0.95; N2 0.97; O2 0.97; CO2 0.97; SF6 0.99; CO 1.00 |
| Shelton & Buckingham 1982, ESHG, 515 nm | `shelton_1982` | He 0.90 |
| Lundeen et al. 1983, CARS, 532 nm | `lundeen_1983` | He 0.83; Ar 1.12; Kr 1.20; D2 1.01; N2 1.15; O2 1.39; CO2 1.29; CH4 1.16; SF6 1.31; CO 1.04 |
| Hahn & Lee 1995, CARS, 532 nm | `hahn_1995` | CO2: not compared (no Stokes wavelength); CH4: not compared (no Stokes wavelength); CO: not compared (no Stokes wavelength) |
| Buckingham & Dunmur 1968, dc Kerr, 633 nm | `buckingham_1968` | He 1.21; Ne 0.92; Ar 0.95; Kr 0.99; Xe 1.03; SF6 1.62 |
| Shelton & Rice 1994, dc Kerr (review), 515-633 nm | `shelton_1994` | He 1.11 |
| Tammer et al. 1992, absolute dc Kerr (field calibrated), 633 nm | `tammer_1992` | He 1.00; H2 1.03 |
| Buckingham & Orr 1969, dc Kerr, 633 nm | `buckingham_1969` | CH4 1.03 |
| Gentle et al. 1989, dc Kerr, 633 nm | `gentle_1989` | CO2 1.63 |
| Hättig & Jørgensen 1998, CCSD Kerr and THG | `hattig_1998` | CH4 (ccsd_kerr_CH4) 0.89; CH4 (ccsd_thg_CH4) 0.85 |
| Maroulis 1999, CCSD(T) static | `maroulis_1999` | O3 0.97 |
| Maroulis 1994, static (correlated) | `maroulis_1994` | O3 1.20 |
<!-- nonlinear-validators:end -->

- **The Maryland n2** (single-shot supercontinuum spectral interferometry, SSSI) mostly
  lie below the ESHG scale: 0.84-1.03 at 800 nm, 0.95-0.96 at 400-2400 nm (2015) and
  0.88-1.02 at 3-10.5 µm (2019). Their stated uncertainties are 7-15 %.
  - Their absolute scale is a nitrogen flow whose index corresponds to 290.9 K at 1 atm,
    and the 2012 and 2015 values are reduced at that density (0.7 % above 293.15 K;
    TC-316). Zahedpour et al.'s scale comes through the anisotropy and an unstated
    density, and stays at 293.15 K.
  - Shelton's suggested +1 % would widen the gap, not close it.
  - The 5-15 % offset is the largest open question for the Kerr level. The fits keep the
    ESHG scale, which is tied to an ab initio helium.
- **Other n2:**
  - Wang et al. agree at 800 nm (1.05 for Ar, 1.02 for Ne) but not at 1800 nm (0.68).
    Their 1800 nm point is inconsistent with the Maryland group's flat infrared argon
    (TC-313, open).
  - Loriot et al.'s argon agrees (0.97). Their molecules include the rotational response
    (1.36 for N2, 1.86 for O2), as do Nibbering et al.'s "instantaneous" values (N2 2.9,
    O2 6.1).
  - Börzsönyi et al.'s values are 0.8-2 times the fits' after correcting their table's
    header, which reads 10^-19 where 10^-20 fits their text and their own literature
    column (an erratum on circumstantial evidence; TC-316). As printed they were 8-20
    times.
- **Third-harmonic generation:**
  - Lehmeier et al.'s THG lies 4-31 % below the fitted model (Ne 0.69, Ar 0.77, N2 0.88).
  - Ward & New's and Ward & Elliott's THG agree within about 10-15 %.
- **An independent ESHG scale:** Ward & Miller's ESHG at 694 nm agrees within 5 %
  (0.95-1.00).
- **CARS:** Lundeen et al.'s absolute values scatter from 0.83 to 1.39, with their stated
  13 % scale. Hahn & Lee's per-amagat values need a Stokes wavelength that is not
  recorded.
- **dc Kerr:**
  - **Tammer et al.'s** field-calibrated absolute measurement (1992; the field from the
    Stark effect of OCS) gives helium 1.00 +- 0.02 (TC-401). Helium's high values in
    the older dc Kerr data (Buckingham & Dunmur 1.21; Shelton & Rice's review average
    1.11) are those measurements', not the model's.
  - **Tammer's H2** (1.03 +- 0.03, electronic plus vibrational, the orientation term
    subtracted) is the only absolute dc Kerr level independent of both the ESHG scale
    and helium theory.
  - **Buckingham & Orr's CH4** (1969) agrees (1.03, +-5 %).
  - **Buckingham & Dunmur's atoms** agree within 8 %.
  - **Molecules with a vibrational part:** SF6's 1.62 presumably contains the
    vibrational contribution that a static field brings, which the model does not
    have. CO2's 1.63 +- 0.42 (Gentle et al.; the orientation term extrapolated away,
    their uncertainty being the fit's scatter) is 1.5 sigma high.
- **Theory:**
  - Maroulis's two static O3 values bracket the model (0.97 and 1.20), which is the
    method dependence the 12 % allows for.
  - Hättig & Jørgensen's CH4 lies 11-15 % below.

<span id="against-solprop"></span>


## Compared with the state of the art

- **The reference data are Shelton's 1990 ESHG scale and his dispersion fits** (Table
  IV: gamma_0, A, B against $W_2$ for He to Xe, H2, D2, N2, O2), on Bishop & Pipin's
  helium. The fits here are built on them and agree within 0.4 % in gamma_0 (TC-314).
- **What is added:**
  - the process dependence at fourth order (B') from ab initio calculations, each at a
    fitted scale;
  - the conversion of every other source convention into one currency, so that n2,
    THG, CARS and dc Kerr can be compared with the same curve;
  - a covariance and a leave-one-out record per gas;
  - CO, N2O and O3, for which no dispersion fit existed.
- **The independent absolute n2 of the Maryland group** are the modern counterpart.
  Their stated 7-15 % covers most of the offset from the ESHG scale; the ESHG scale is
  the more precise.
- **Ab initio:** coupled-cluster values for Ne, Ar, Kr, H2 and N2 lie within a few per
  cent of the ESHG scale (TC-313). Helium's (Bishop & Pipin) is better than any
  measurement.

The fits have not been reviewed outside this project.

## What the uncertainties leave out

- **The truncation of the dispersion form:** outside the fitted $W_2$ ({doc}`/user-guide/nonlinear-support`'s
  table gives each gas's shortest Kerr and THG fundamental; the model warns there),
  beyond the fitted legs (see {doc}`/user-guide/nonlinear-extrapolation`),
  and near resonance (O3; the resonance limits of every gas are in
  {doc}`/user-guide/nonlinear-support`).
- **Shelton's possible +1 % helium pair-interaction bias:** it is covered by the floor,
  not applied.
- **The vibrational and rotational parts:** the models are electronic, except the
  unseparated ESHG vibrational part of CO2, CH4 and SF6 and whatever the N2O n2
  contains. In propagation the Raman response is separate.
- **Local-field factors:** $n_2$ is given in the dilute limit. The Lorentz local-field
  factor raises the nonlinear polarization by about $8(n-1)/3$, or 7e-4 for argon at
  1 bar and 800 nm, growing with density.
- **The unresolved disagreements:**
  - the Maryland n2 against the ESHG scale;
  - Wang et al.'s 1800 nm argon;
  - Börzsönyi et al.'s values, even corrected;
  - krypton's internal disagreement.

## What would improve it

The papers listed on the {doc}`roadmap`, most useful first:
- an absolute dc Kerr of helium (Tammer et al. 1992) and the older absolute dc Kerr of
  the molecules;
- Shelton's 1987 vibrational gamma of N2 and O2;
- the Maryland calibration (Wahlstrand et al. 2016);
- any gamma at all for N2O;
- the supplementary material behind the O3 calculations.

<span id="solprops-demonstration-quartets"></span>

See the [representative four-wave-mixing quartets](user-guide/nonlinear-extrapolation.md#representative-four-wave-mixing-quartets).
