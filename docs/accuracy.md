(accuracy-and-the-state-of-the-art)=
# Gas dispersion: accuracy and limitations

This chapter explains statistical uncertainty, sensitivity to dataset selection and
working estimates for gas dispersion. Read {doc}`models` and {doc}`fitting` first,
or return to the {doc}`quickstart` for a calculation. The tables are checked against
stored results; the {ref}`uncertainty gallery <sphx_glr_auto_gallery_accuracy_plot_dispersion_uncertainty.py>`
shows their wavelength dependence. Source records are in {doc}`/the-data/sources`.

## What the numbers mean

A model can be wrong in two ways, and they matter very differently.

- **The level**: a relative error $\varepsilon$ that scales $n - 1$ at every wavelength.
  It scales every dispersion coefficient $\beta_k$ by the same $\varepsilon$, so it is
  equivalent to a density error of $\varepsilon$. One ppm of density is 0.29 mK at 293 K,
  or 0.1 Pa at 1 bar.
- **The shape**: an error in how $n - 1$ changes with wavelength. It changes the group
  delay, the group-velocity dispersion (GVD, $\beta_2$), the third-order dispersion (TOD,
  $\beta_3$) and every phase-matching condition, and near a resonance it can be large.

For scale, at 800 nm, 1 bar and 20 °C (from the fits):

| Gas | $n - 1$ | $k(n-1)$ | $\beta_2$ |
|---|---|---|---|
| He | 3.19e-5 | 251 rad/m | 0.99 fs²/m |
| Ar | 2.57e-4 | 2019 rad/m | 18.4 fs²/m |
| N2 | 2.72e-4 | 2139 rad/m | 19.8 fs²/m |
| Xe | 6.36e-4 | 4997 rad/m | 87.9 fs²/m |

So a level error of 100 ppm is a phase error of 0.2 rad per metre of argon at 1 bar and
a GVD error of 0.01 %. Who needs what:

- **Refractive-index gas metrology** (optical pressure and temperature standards,
  dielectric-constant gas thermometry) works at the ppm level and below. Its references are
  ab initio polarizabilities (Puchalski et al. 2016 for helium; Lesiuk et al. for neon and
  argon) and the metrology these are checked against (Schmidt et al. 2007; Gaiser &
  Fellmuth; Stone & Stejskal 2004; Egan et al. 2019). The helium and argon models here
  reproduce those references (0.1 ppm and 3 ppm at 800 nm); they add nothing at this level,
  and such work should cite the references directly.
- **Pulse propagation** in gas-filled fibres and cells depends on the shape: GVD and TOD,
  zero-dispersion and phase-matching wavelengths. Here the level hardly matters: a gas fill
  is rarely known to better than a few 0.1 % in density (an order-of-magnitude statement
  about pressure and temperature measurement, not a curated number), which is thousands of
  ppm. What matters is the dispersion to about 0.1-1 %, and the index near resonances.
  In a 75 µm-radius capillary a 1 % change of the gas's GVD moves the zero-dispersion
  wavelength by about 1 nm (helium) to 2 nm (argon) at 1 bar; the
  {ref}`capillary example <sphx_glr_auto_gallery_accuracy_plot_capillary_zdw.py>`
  shows how gas dispersion and pressure set the zero-dispersion wavelength.
- **Phase matching near resonances** (resonant dispersive waves and harmonics in the
  vacuum ultraviolet) depends on $n - 1$ close to the first resonance, where the published
  formulas err by 5-15 % and the fits' $n - 1$ is known to 0.01 % (helium at 65 nm), 0.1 %
  (argon at 120 nm), 0.3-0.4 % (krypton at 140 nm, xenon at 150 nm) and 1.5 % (neon at
  120 nm; the page The gas models, TC-74, TC-91, TC-94).
- **The extreme ultraviolet** (harmonic generation, absorption) needs the complex index
  to a few per cent; see the page on the extreme ultraviolet (TC-99).

## How well the fits know the dispersion

Each fitted model has a parameter covariance, and
{meth}`~optmat.gas.models.fitted.FittedModel.dispersion_sigma_relative` propagates
it into the relative uncertainty of $\beta_k$ (in the dilute limit it does not depend on
the state). That **statistical** uncertainty is a lower bound: it assumes the model form
(two to four effective oscillators) and treats the stated uncertainties as independent. A
second measure is **empirical**: the largest change of $\beta_2$ when one dataset is left
out of the fit (`tools/fit/robustness.py`). Agreement between the two supports the statistical estimate but does not test
all systematic errors or alternative model forms. Where one dataset moves the GVD by many statistical sigma, that dataset alone
fixes the near-infrared dispersion and the other data, extrapolated by the model, disagree
with it.

At 800 nm (TC-100):

<!-- accuracy-table:start -->
| Gas  |  $\beta_2$, statistical $1\sigma$  |  Largest leave-one-out shift of $\beta_2$  |  Dataset left out  |  Working estimate |
| --- | --- | --- | --- | --- |
| He  |  0.00040 %  |  0.62 %  |  puchalski_2016/cauchy_total  |  < 0.01 % |
| Ne  |  0.21 %  |  0.34 %  |  cuthbertson_1932/experiment_Ne  |  0.3 % |
| Ar  |  0.016 %  |  0.010 %  |  lesiuk_2023/cauchy_Ar  |  0.05 % |
| Kr  |  0.11 %  |  0.30 %  |  hohm_1994b/alpha_Kr  |  0.5 % |
| Xe  |  0.090 %  |  0.13 %  |  hohm_1994b/alpha_Xe  |  0.4 % |
| N2  |  0.012 %  |  0.17 %  |  peck_1966/experiment  |  0.2 % |
| H2  |  0.027 %  |  0.069 %  |  peck_1977/experiment  |  0.1 % |
| D2  |  0.016 %  |  0.029 %  |  shelton_1990/density_H2_over_D2  |  0.1 % |
| O2  |  0.055 %  |  0.040 %  |  smith_1976/experiment_O2  |  0.1 % |
| CH4  |  0.076 %  |  0.078 %  |  hohm_1993/alpha_CH4  |  0.2 % |
| SF6  |  0.27 %  |  1.3 %  |  shelton_1985_sf6/density_N2_over_SF6  |  2 % |
| N2O  |  0.66 %  |  1.7 %  |  cuthbertson_1914/n2o_dispersion  |  3 % |
| CO2  |  0.054 %  |  0.21 %  |  old_1971/experiment  |  0.2 % |
| O3  |  0.52 %  |  1.1 %  |  cuthbertson_1914/ozone_adopted  |  50 % |
| CO  |  0.21 %  |  0.29 %  |  cuthbertson_1920/experiment_CO  |  0.3 % |
<!-- accuracy-table:end -->

In the ultraviolet (343 nm) and the infrared (1800 nm), the same three measures ("n/a"
where the wavelength is outside the model's declared range):

<!-- accuracy-bands:start -->
| Gas | 343 nm: statistical | 343 nm: leave-one-out | 343 nm: working | 1800 nm: statistical | 1800 nm: leave-one-out | 1800 nm: working |
|---|---|---|---|---|---|---|
| He | 0.00037 % | 0.77 % | < 0.01 % | 0.00040 % | 0.60 % | < 0.01 % |
| Ne | 0.56 % | 2.0 % | 1 % | 0.22 % | 0.68 % | 0.5 % |
| Ar | 0.025 % | 0.14 % | 0.1 % | 0.017 % | 0.013 % | 0.05 % |
| Kr | 0.15 % | 0.77 % | 1 % | 0.15 % | 0.30 % | 0.5 % |
| Xe | 0.16 % | 0.36 % | 1 % | 0.12 % | 0.20 % | 0.5 % |
| N2 | 0.16 % | 0.081 % | 0.2 % | 0.034 % | 0.19 % | 0.2 % |
| H2 | 0.12 % | 0.18 % | 0.2 % | 0.042 % | 0.10 % | 0.1 % |
| D2 | 0.022 % | 0.32 % | 0.3 % | 0.021 % | 0.052 % | 0.1 % |
| O2 | 0.17 % | 0.42 % | 1 % | 0.079 % | 0.055 % | 0.1 % |
| CH4 | 0.12 % | 0.71 % | 1 % | n/a | n/a | n/a |
| SF6 | 0.22 % | 2.0 % | 2 % | 0.075 % | 0.41 % | 1 % |
| N2O | 1.5 % | 27 % | 5 % | 0.58 % | 3.1 % | 5 % |
| CO2 | 0.35 % | 0.47 % | 1 % | n/a | n/a | n/a |
| O3 | 0.012 % | 0.025 % | 50 % | 0.46 % | 0.93 % | 1 % |
| CO | 0.19 % | 3.8 % | 4 % | 0.51 % | 1.4 % | 1.5 % |
<!-- accuracy-bands:end -->

In the ultraviolet each model leans on one dataset or its theory backbone (dropping it
moves the GVD at 343 nm by 1-4 %): helium's and argon's ab initio moments, xenon's
Bideau-Mehu (its 0.1 % point at 170.2 nm), oxygen's Smith et al., methane's and SF6's
Shelton ratios (whose harmonics reach 230-340 nm), CO2's Shelton ratios, and Hohm's 325 nm
point for N2O (without it N2O's range would stop at 450 nm). Helium's and argon's working estimates stay at their
theory's accuracy; the others are set near the leave-one-out shift. In the infrared the
estimates follow 800 nm, except neon and N2O (visible data only).

The working estimate is a judgement, gas by gas:

- **He**: the theory (Puchalski's dispersion coefficients with relativistic and QED
  corrections) is the backbone. Without it the experiments alone give a GVD 0.7 %
  different; with it they fit at a reduced chi-squared of 0.58. Independent evaluations of
  the theory (Stone & Stejskal's formula, Puchalski's own series) agree with the fit's GVD
  to 1e-4. The refractive model adds Pachucki & Puchalski's retardation term, -24 ppm of
  the GVD (6 statistical sigma); with it Yang et al.'s measured two-colour dispersion
  coefficient agrees to 0.6 of its sigma, without it 4.8 (TC-705, TC-728).
- **Ne**: the leave-one-out shifts are about 1.5 statistical sigma; the ab initio
  dispersion is 2.5 sigma below every experiment (TC-92), an open question.
- **Ar**: statistical and empirical agree; the limit is Lesiuk & Jeziorski's $S(-4)$
  (0.05 %). Shelton's He/Ar ratio is 0.3 % from the model at 1.9 sigma (TC-96), open.
- **Kr, Xe**: the near-infrared dispersion rests on Shelton's ratios to argon; without them
  the visible and ultraviolet data, extrapolated, gave a GVD 0.95 % (Kr) and 0.35 % (Xe)
  lower. Adding the ratios moved the fits by 8 and 3 statistical sigma, which is how the
  statistical figure was found to be too small (TC-96). Hohm's 325-633 nm dispersion
  (1994), fitted for its shape with his level free (TC-400), now agrees with the ratios:
  dropping it moves the GVD at 800 nm by 0.30 % (Kr) and 0.13 % (Xe), the largest
  leave-one-out shifts there.
- **N2**: Peck & Khanna fix the near-infrared dispersion. Shelton's He/N2 ratio lies
  0.3-0.4 % from the fit (1-2 sigma), and his N2/Ar ratio -0.1 %.
- **H2**: on theory since 2026-09-25 (TC-95): Raj et al.'s ab initio average (182 nm to
  the static limit) with a 0.1 % allowance, the experiments' levels free or within their
  stated systematics. Peck & Hung's shape agrees; Shelton's He/H2 ratio to 0.02 %.
- **D2**: Raj et al.'s theory and Shelton's ratios to hydrogen, which agree within their
  errors at the 0.1 % allowance.
- **O2**: Egan et al.'s 2026 refractometry at 633 and 1542 nm (7e-6) fixes both the level,
  in agreement with Birch, and the near-infrared dispersion; with Hohm's 325-633 nm
  dispersion (its shape at 4e-5, TC-400), no dataset moves the GVD at 800 nm by more than
  0.04 % (TC-147).
- **CH4**: Shelton's 14 ratios (477-1319 nm), two visible dispersion sets, Achtermann at
  633 nm, Zhao et al.'s 2026 value at 1566 nm and Hohm's 325-633 nm dispersion (1993)
  agree (reduced chi-squared 0.69); the working estimate includes nitrogen's 0.2 %.
  Since the refit with its bands as rules (TC-797) the optical range ends at 1313 nm,
  so Shelton's 1319 nm ratio lies just beyond it. With
  Hohm's ultraviolet shape the largest leave-one-out shift at 343 nm fell from 2.1 % to
  0.71 % (Shelton's 1986 ratios), hence 1 % there (TC-400). The static dielectric value
  of Bose et al. (1972) agrees with the model's, fixed infrared bands and all, to 59 ppm.
- **SF6**: the infrared bands (HITRAN, fixed) give anomalous dispersion that cancels the
  electronic dispersion near 1.2 µm, where the gas's own GVD crosses zero (the relative
  uncertainty diverges there). The two Shelton ratio sets prefer HITRAN's strengths times
  1.036 +- 0.009 (TC-146), and dropping either set moves the GVD at 800 nm by up to 1.3 % (1.8 % before the refit of
  2026-09-29, TC-725):
  2 %.
- **CO2**: Old, Gentili & Peck's 480-1600 nm refractivity (7 ppm) sets the shape; the
  HITRAN bands are needed (their strength, scanned, is what the data want), now as
  Gaussian rules of their lines (OEP 3). The optical range stops at 1.4 µm because weak
  combination bands lie from there to 2.3 µm. The model holds again in seven windows
  from 2.30 to 8.94 µm; its laser bands are kept, so it is not qualified at 8.9-12.5 µm,
  where Marchetti's values lie: they are validators, 0.76 % above the model (TC-797).
  Reduced chi-squared 4.0 (3.4 before the refit; TC-148, TC-400, TC-797). The two weak bands below 10 eV (1Delta_u, 1Pi_g; Chan et al.) are
  lines with priors and the effective oscillators start at 10 eV: with a free lower
  oscillator the ultraviolet shape was set by its bound (TC-725, TC-726). The largest
  leave-one-out shift of the GVD at 343 nm is 0.47 % (without Shelton's 1986 ratios),
  so the working estimate there is 1 %. Bideau-Mehu et al.'s 181-254 nm values scatter by
  0.3 % about the model and are a validator (TC-178).
- **N2O**: the dispersion rests on Hohm's 325-633 nm values (their shape at 4e-5, the
  level free; TC-400) as well as the Cuthbertsons' 1914 (480-671 nm), and on the infrared
  bands. The two visible sets differ in shape by about 1e-4 (reduced chi-squared 3.0),
  and the fit follows Hohm's (TC-178). Dropping Hohm's set moves the GVD 1.6 % at 800 nm
  and 27 % at 343 nm, where his 325 nm point is the only datum. Working estimates: 3 % at 800 nm, 5 % at 343 and
  1800 nm.
- **CO**: Smith et al.'s 167-288 nm interferometry, the Cuthbertsons' visible
  dispersion, Watson's green line and Zhao et al.'s 1566 nm refractometry (38 ppm) agree
  (reduced chi-squared 0.37). The ultraviolet rests on Smith alone: without it the GVD at
  343 nm moves 3.9 %, hence 4 %; at 800 nm 0.3 %; at 1800 nm, beyond the data but for
  Zhao's point, 1.5 %. Nothing is measured between 290 and 480 nm.
- **O3**: a different kind of model. The measured absorption (Serdyuchenko et
  al. 2014, 293 K, smoothed to 50 cm^-1) enters by Kramers-Kronig. Only the strength at
  12-20.6 eV is fitted, to Swift et al., the Cuthbertsons' adopted level and the ab
  initio static value plus Naves et al.'s vibrational correction. These agree (reduced
  chi-squared 0.6), and the fitted strength (2.55) lies near Kumar & Thakkar's B1.
  - **n - 1:** known to about 2 %. The largest leave-one-out shift is 1.9 %. Sharipov's
    small vibrational correction instead of Naves's would lower it by 2.5 %.
  - **GVD at 343 nm (Huggins band) and 800 nm (Chappuis-Wulf edge):** set by the band
    structure of the measured spectrum. Building the spectrum from Bogumil's or Brion's
    independent data changes it by 16-52 %, the 243 K spectrum by 23-34 %, and 25 or
    100 cm^-1 smoothing by 6-40 %. Hence the 50 % working estimate. The leave-one-out
    shifts (1 %) and the statistical figure do not measure this.
  - **GVD at 1800 nm:** clear of the bands; every such change is below 0.5 %.
  - **The Cuthbertsons' visible dispersion:** 20 % steeper than the spectrum allows (a
    validator, TC-218).

A gas fitted to a ratio inherits its partner model's uncertainty, which its covariance
does not contain (the partner is held fixed); the working estimates include it. The TOD is
known about as well as the GVD relatively (the gallery page shows both). Near the short
end of each range the uncertainties grow: see the gallery page.

## Compared with the state of the art

| Gas | Reference work in the library | Relationship to the optmat model | Coverage and differences | Limits |
|---|---|---|---|---|
| He | Puchalski et al. 2016 ab initio $\alpha(\omega)$ with QED; Schmidt et al. 2007 | Equal (reproduces Puchalski's $\alpha$ to 0.1 ppm) | A usable form down to 60 nm; the density conversion | None at the level of the theory |
| Ne | Lesiuk et al. 2020 ab initio moments; Gaiser & Fellmuth 2018 DCGT | Equal in level; dispersion carries TC-92, with the photoabsorption now on the other side of theory (TC-126) | The vacuum ultraviolet to 76 nm, as a physical pole model | No data below 140 nm; TC-92 |
| Ar | Lesiuk & Jeziorski 2023 ab initio; DCGT; Egan et al. 2019 | Equal (built on them; 3 ppm at 800 nm) | 108-140 nm, where Cauchy series diverge and the Kingston fit errs by up to 15 % | Sparse vacuum-ultraviolet data (Chashchina not in the library) |
| Kr, Xe | Bideau-Mehu 1976 and 1981; Abjean et al. 1971; Smith et al. 1976; Burns et al. 1986; Kumar & Thakkar 2010 distributions; dielectric static values | On par: built on much the same data. Kumar & Thakkar's pseudospectra lie +0.18 % (Kr) and +0.33 % (Xe) above the fits in level, +0.34 % and +0.72 % in GVD at 800 nm (TC-127) | Consistent level (Börzsönyi's is 0.3-0.9 % high); poles at the real lines | No converged ab initio dispersion (Pawłowski et al.'s CC3 krypton moments are not, TC-175); the near infrared rests on Shelton's ratios |
| N2 | Peck & Khanna 1966; Griesmann & Burnett 1999; the 633 nm metrology | On par in the visible and near infrared | One model from 145 nm to the static limit | Refers to 293.15 K; the temperature dependence (1.18 ppm/K) is not carried by the model (TC-730); chi_0 from ab initio theory alone (TC-708) |
| H2, D2 | Raj et al. 2018, Rychlewski 1983 and Bishop & Cheung 1980 theory; Peck & Hung 1977 | Equal: built on Raj et al.'s ab initio average | One consistent model from 112 nm, experiments reconciled | The experiments' levels scatter by 1e-3 about theory (TC-95); temperature dependence (0.03 % between 273 and 296 K) not modelled |
| O2 | Smith et al. 1976; Křen 2011 on Zhang 2008; Egan et al. 2026 | On par, with Egan 2026 at its core | The static limit | From 195 nm only |
| CO2 | Old, Gentili & Peck 1971; the 633 nm metrology; Marchetti 2007 (10 µm) as validators | On par: built on them, with HITRAN's bands | One model from 250 nm to 1.4 µm and in seven windows from 2.30 to 8.94 µm | Ultraviolet from Hohm's 325 nm point and Shelton's harmonics; Bideau-Mehu 1973 (181-254 nm) scatters by 0.3 %, a validator (TC-178) |
| CH4, SF6, N2O | He et al. 2021 formulas (Rayleigh-derived, no infrared term) | Fitted with HITRAN's infrared bands (TC-104, TC-146) | Near-infrared dispersion and phase matching: the formulas miss Shelton's ratios by 3-28 % and SF6's GVD by 55 % | N2O's dispersion rests on visible data and Hohm's 325 nm point (TC-178); the SF6 infrared strength question (TC-146) |
| EUV (He, Ne, Ar, Kr, Xe) | Henke et al. 1993 (the CXRO tables); Chan et al., Samson et al. absorption | Complementary | Near the absorption thresholds, where Henke's tables are not reliable | The approximations of TC-99 and TC-126 |

What is new is the combination, not any single number: one physically constrained model
per gas from the vacuum ultraviolet to the static limit (poles at the real resonances, so
extrapolation is controlled), the reference state handled correctly (real-gas density,
Lorentz-Lorenz, the refractivity virial and the magnetic term; the published formulas lack
these, which is how Ermolov's 1 bar and 1 atm were confused), a stated uncertainty with
its covariance, and the provenance of every input. The fits have not been reviewed
outside this project.

## What the uncertainties leave out

- The model form: two to four effective oscillators. The empirical column above measures
  its effect where a dataset constrains the near infrared on its own.
- Correlated systematics: a theory allowance (hydrogen's 0.1 %) enters as independent per
  point, and dataset scales cover only the systematics the sources state.
- The partner of a ratio, held fixed.
- The temperature dependence of hydrogen's rotationally averaged polarizability.
- Unresolved disagreements between sources: TC-92 (neon, decided on theory), TC-96 (argon
  against Shelton), TC-146 (SF6's infrared strength).
- The molecules' line structure near their weak near-infrared bands. Inside the optical
  ranges the infrared terms are Gaussian rules, and the difference from the line-by-line
  sum is reported, not enforced (OEP 3, TC-749). Against the sum smoothed to 1 cm$^{-1}$,
  the group-velocity dispersion stays within 1 % of the larger of its value and its
  electronic part up to about 900 nm for every molecule. Beyond, the fine structure of the
  overtone bands exceeds that. A pulse much broader than the smoothing sees less of it:
  smoothed to 100 cm$^{-1}$, the 1 % holds to 1.47 µm for CO, 1.18 µm for CO2 and
  1.14 µm for N2O. For CH4 (940 nm) it does not move, because there the fit's own poles
  lie at the edges of its narrow absorbing gaps (TC-808).

## What would improve it

The papers listed on the {doc}`roadmap` under *Papers and data we are looking for*:
the vacuum-ultraviolet refractometry of the 1960s-70s, converged ab initio dispersion
and magnetic susceptibilities for krypton and xenon, and more dispersion and infrared
data for the polyatomic gases.
