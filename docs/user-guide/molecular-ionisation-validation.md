# Molecular ionisation: validation and limitations

Read {doc}`/molecular_ionisation_rates` for the working example, rate and model choices. Here the comparisons distinguish all-electron static rates, vibrational effects and measured yields. The callable validation functions are in {mod}`optmat.gas.molecular_ionisation_validation`.

<span id="validation"></span>

## Against all-electron static rates

Majety & Scrinzi (2015, haCC) give static rates of N2 at seven angles between the
field and the molecular axis, of CO at thirteen (0 to 180°: CO is polar, and its rate
differs by up to 2.1 between the two directions along the axis) and of H2 at two.
Jagau (2018, coupled cluster in complex-scaled Gaussian bases; its supplement's
tables) gives O2 at nine angles (at 0.06 and 0.10 au), CO at thirteen (with 0° the
reverse of Majety's), H2 at two, and CH4 in three orientations, which Lebedev's
14-point rule averages as 0.3, 0.3 and 0.4 (along a CH bond towards H and towards
C, and bisecting two). All are for fixed nuclei, with all electrons. Their average
over the orientations
({func}`~optmat.gas.molecular_ionisation_validation.orientation_average`) is what the
rate of an unaligned molecule is to reproduce in the static limit. The models at fixed
nuclei over that average, and the anisotropy of the reference:

<!-- molecular-rate-static:start -->
| Molecule | Reference |  | 0.025 au | 0.03 au | 0.04 au | 0.05 au | 0.06 au | 0.07 au | 0.08 au | 0.09 au | 0.1 au | 0.12 au | 0.14 au |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N2 | majety_2015 | reference: rate at the first angle over its average |  | 2.57 | 2.47 | 2.36 | 2.21 | 2.02 | 1.83 | 1.67 | 1.53 | 1.33 | 1.21 |
| N2 | majety_2015 | `default` |  | 0.99 | 0.99 | 1.03 | 1.07 | 1.11 | 1.13 | 1.13 | 1.10 | 0.98 | 0.81 |
| N2 | majety_2015 | `moadk_zhao_2010` |  | 0.90 | 1.10 | 1.39 | 1.78 | 2.24 | 2.78 | 3.39 | 4.04 | 5.37 | 6.64 |
| N2 | majety_2015 | `moadk_zhao_2017` |  | 1.44 | 1.76 | 2.23 | 2.85 | 3.60 | 4.46 | 5.44 | 6.47 | 8.60 | 10.64 |
| N2 | majety_2015 | `moadk_tong_2002` |  | 0.50 | 0.61 | 0.78 | 0.99 | 1.25 | 1.55 | 1.89 | 2.25 | 2.99 | 3.69 |
| N2 | majety_2015 | `adk_hydrogenic` |  | 0.50 | 0.60 | 0.76 | 0.97 | 1.22 | 1.51 | 1.83 | 2.17 | 2.86 | 3.51 |
| CO | majety_2015 | reference: rate at the first angle over its average | 2.69 | 1.42 | 1.36 | 1.24 | 1.12 | 1.02 | 0.92 | 0.83 |  |  |  |
| CO | majety_2015 | `default` | 2.25 | 1.26 | 1.29 | 1.29 | 1.31 | 1.32 | 1.31 | 1.26 |  |  |  |
| CO | majety_2015 | `moadk_zhao_2010` | 2.16 | 1.38 | 1.83 | 2.38 | 3.12 | 4.08 | 5.24 | 6.53 |  |  |  |
| CO | majety_2015 | `moadk_zhao_2017` | 3.48 | 2.22 | 2.96 | 3.85 | 5.07 | 6.63 | 8.52 | 10.63 |  |  |  |
| CO | majety_2015 | `moadk_tong_2002` | 0.67 | 0.42 | 0.56 | 0.73 | 0.95 | 1.25 | 1.59 | 1.98 |  |  |  |
| CO | majety_2015 | `adk_hydrogenic` | 1.01 | 0.64 | 0.85 | 1.10 | 1.44 | 1.87 | 2.39 | 2.96 |  |  |  |
| CO | jagau_2018 | reference: rate at the first angle over its average |  |  |  |  | 2.05 |  |  | 1.87 |  |  |  |
| CO | jagau_2018 | `default` |  |  |  |  | 1.54 |  |  | 1.14 |  |  |  |
| CO | jagau_2018 | `moadk_zhao_2010` |  |  |  |  | 3.67 |  |  | 5.95 |  |  |  |
| CO | jagau_2018 | `moadk_zhao_2017` |  |  |  |  | 5.94 |  |  | 9.69 |  |  |  |
| CO | jagau_2018 | `moadk_tong_2002` |  |  |  |  | 1.12 |  |  | 1.81 |  |  |  |
| CO | jagau_2018 | `adk_hydrogenic` |  |  |  |  | 1.69 |  |  | 2.70 |  |  |  |
| H2 | majety_2015 | reference: rate at the first angle over its average |  |  |  |  | 1.19 |  | 1.17 |  | 1.15 | 1.13 |  |
| H2 | majety_2015 | `default` |  |  |  |  | 1.00 |  | 0.97 |  | 1.00 | 1.03 |  |
| H2 | majety_2015 | `moadk_zhao_2010` |  |  |  |  | 2.05 |  | 1.87 |  | 1.98 | 2.19 |  |
| H2 | majety_2015 | `moadk_zhao_2017` |  |  |  |  | 1.91 |  | 1.75 |  | 1.85 | 2.04 |  |
| H2 | majety_2015 | `moadk_tong_2002` |  |  |  |  | 4.06 |  | 3.71 |  | 3.92 | 4.33 |  |
| H2 | majety_2015 | `adk_hydrogenic` |  |  |  |  | 2.98 |  | 2.72 |  | 2.88 | 3.18 |  |
| H2 | jagau_2018 | reference: rate at the first angle over its average |  |  |  |  | 1.52 |  | 1.20 |  | 1.14 | 1.13 | 1.12 |
| H2 | jagau_2018 | `default` |  |  |  |  | 1.30 |  | 0.99 |  | 1.01 | 1.04 | 1.04 |
| H2 | jagau_2018 | `moadk_zhao_2010` |  |  |  |  | 2.67 |  | 1.91 |  | 1.99 | 2.20 | 2.42 |
| H2 | jagau_2018 | `moadk_zhao_2017` |  |  |  |  | 2.49 |  | 1.78 |  | 1.87 | 2.06 | 2.26 |
| H2 | jagau_2018 | `moadk_tong_2002` |  |  |  |  | 5.28 |  | 3.78 |  | 3.95 | 4.36 | 4.79 |
| H2 | jagau_2018 | `adk_hydrogenic` |  |  |  |  | 3.88 |  | 2.77 |  | 2.90 | 3.20 | 3.51 |
| O2 | jagau_2018 | reference: rate at the first angle over its average |  |  |  |  | 0.50 |  |  |  | 0.86 |  |  |
| O2 | jagau_2018 | `default` |  |  |  |  | 1.10 |  |  |  | 0.95 |  |  |
| O2 | jagau_2018 | `moadk_zhao_2010` |  |  |  |  | 1.95 |  |  |  | 2.81 |  |  |
| O2 | jagau_2018 | `moadk_zhao_2017` |  |  |  |  | 2.96 |  |  |  | 4.25 |  |  |
| O2 | jagau_2018 | `moadk_tong_2002` |  |  |  |  | 2.77 |  |  |  | 3.98 |  |  |
| O2 | jagau_2018 | `adk_hydrogenic` |  |  |  |  | 20.11 |  |  |  | 25.56 |  |  |
| CH4 | jagau_2018 | reference: rate at the first angle over its average |  |  |  |  | 0.84 |  | 0.89 |  | 0.90 | 0.94 |  |
| CH4 | jagau_2018 | `default` |  |  |  |  | 1.27 |  | 1.25 |  | 1.03 | 0.74 |  |
| CH4 | jagau_2018 | `adk_hydrogenic` |  |  |  |  | 5.10 |  | 5.12 |  | 5.40 | 5.53 |  |
<!-- molecular-rate-static:end -->

- **The default is within 3 % of H2's references and 13 % of N2's to 0.12 au**
  (1.5 times N2's barrier-suppression field), within 10 % of O2's and 27 % of CH4's;
  at 0.14 au it is 19 % low for N2. Its constants were fitted to these references:
  this shows that one constant per molecule is enough, not that the rate is right
  elsewhere.
- **The two references agree where both exist.** For H2 the default is within 4 %
  of both from 0.08 au (Jagau's 0.06 au rate, 2e-5 au, is near his method's floor).
  For CO, Majety's average is 18 % above Jagau's at 0.06 au and 13 % at 0.09 au, as
  Jagau reports; the default is 1.54 and 1.14 times Jagau's.
- **With N2's $\alpha$ O2's rate was 3 to 8 times below Jagau's** at 0.06 to 0.10 au:
  the correction outran it. Jagau places O2's barrier suppression at about 0.12 au,
  not the 0.05 au of the atomic formula, and beyond 0.10 au his rate keeps rising
  as a tunnelling rate would, where the default's range has ended (TC-698).
- **CO tests that N2's constant carries to another molecule.** It takes N2's $\alpha$
  without a fit, and the default is 1.26 to 1.32 times CO's reference from 0.03 to
  0.09 au (0.45 to 1.36 times CO's barrier-suppression field): the dependence on the
  field is right and the magnitude 30 % high, within CO's factor of 1.5. Fitted to CO,
  $\alpha$ would absorb the offset as a slope. At 0.025 au every model jumps against
  the reference, by 1.6 to 1.8, and the reference's anisotropy from 1.4 to 2.7: that
  point is taken as the reference's, not the models' (TC-683).
- **The forms have no correction at finite field**, and fail as the uncorrected
  tunnelling rate of an atom does: for N2 the `moadk` form with Zhao et al.'s
  coefficients of 2010 is within 10 % of the reference at 0.03 to 0.04 au and 2.8
  times too high at the barrier-suppression field, 0.08 au (5e14 W/cm² is 0.12 au).
  With their coefficients of 2017 it is 1.6 times higher still.
- **For H2 the form with Zhao et al.'s coefficients is twice too high at every
  field**, and with Tong et al.'s four times. The reference is for fixed nuclei at
  1.4 bohr, where the ionisation energy is 0.6 to 0.9 eV above the adiabatic one that
  the form is evaluated with; with 16.06 eV, the vertical energy of Awasthi et al.'s
  (2008) configuration-interaction calculation, the form is 0.8 to 1.3 times the
  reference (TC-676).
- **The anisotropy of N2 falls with the field**: the rate along the axis is 2.6 times
  the average at 0.03 au and 1.2 times at 0.14 au. The tunnelling form has one
  anisotropy at all fields.
- For H2 only two angles exist; the average takes the rate as $w_\perp +
  (w_\parallel - w_\perp)\cos^2\theta$.

No reference of this kind exists for CO2, N2O or SF6.

```{figure} /auto_gallery/ionisation/images/sphx_glr_plot_molecular_ionisation_001.*
:target: ../auto_gallery/ionisation/plot_molecular_ionisation.html
:alt: The molecular rates of the default and the forms, and against the all-electron static rates

The default and the forms for N2, CO, O2 and H2 at 800 nm; the all-electron static
rates of N2 by angle with their average over orientations; the default at fixed nuclei
and the static forms over that average for N2, CO and H2.
```

**The vibrational factor of H2** against Saenz's rates averaged over the vibrational
ground state, over his rates at 1.4 bohr (the field along the axis):

<!-- molecular-rate-vibration:start -->
| F (au) | Saenz | exp(c/F^2) | Model / Saenz |
|---|---|---|---|
| 0.06 | 1.562 | 1.581 | 1.013 |
| 0.08 | 1.320 | 1.294 | 0.981 |
| 0.1 | 1.158 | 1.179 | 1.019 |
| 0.12 | 1.157 | 1.121 | 0.969 |
<!-- molecular-rate-vibration:end -->

## Against measured yields

Wahlstrand et al. (2018) measured the absolute yield of N2 and O2 as of the noble
gases, with a 42 fs pulse at 800 nm and the intensity calibrated through the Kerr
effect (13 %). For the atoms the rate with the static Coulomb factor (`default_static`,
the atoms' default until TC-848) is 1.7 to 3.9 times above that measurement, and 1.0
within calibration against Wallace et al.'s ({doc}`ionisation-measured`, TC-660). The
last column divides the default's ratio for the molecule by that for an atom of the
same experiment at the same intensity, with the atom's static factor as the
molecule's, which removes the intensity scale to first order:

<!-- molecular-rate-wahlstrand:start -->
| Molecule | I (TW/cm2) | Measured yield | `default` | `moadk_zhao_2010` | `ppt_talebpour` | `ppt_talebpour_ilkov` | `default`, over the atom's `default_static` |
|---|---|---|---|---|---|---|---|
| N2 | 75 | 3.40e-04 | 8.79 | 0.92 | 2.98 | 1.23 | 2.67 (Ar) |
| N2 | 89 | 9.90e-04 | 7.75 | 1.35 | 2.87 | 1.29 | 2.10 (Ar) |
| N2 | 105 | 2.88e-03 | 6.47 | 1.76 | 2.65 | 1.29 | 1.65 (Ar) |
| N2 | 125 | 8.40e-03 | 5.09 | 2.05 | 2.32 | 1.23 | 1.29 (Ar) |
| O2 | 40 | 3.43e-04 | 5.02 | 0.33 | 1.40 | 0.67 | 0.65 (Xe) |
| O2 | 55 | 1.28e-03 | 6.16 | 1.02 | 2.07 | 1.11 | 0.81 (Xe) |
| O2 | 76 | 4.76e-03 | 6.56 | 2.19 | 2.70 | 1.63 | 1.19 (Xe) |
| O2 | 105 | 1.77e-02 | 5.91 | 3.40 | 3.04 | 2.04 | 2.39 (Xe) |
<!-- molecular-rate-wahlstrand:end -->

- **N2 is 5 to 9 times above the measurement**, and 1.3 to 2.7 times above it relative
  to argon. Part of the latter is argon's rate being low at these fields: argon's
  all-electron static rates are 1.2 to 1.6 times its model's below 1.7e14 W/cm²
  (TC-669).
- **O2 is 5.0 to 6.5 times above**, as N2 is, since its $\alpha$ was fitted to Jagau's
  rates (with N2's it was 2.2 to 2.9). Relative to xenon it is 0.65 to 2.4, which
  brackets the measurement (0.34 to 0.90 with N2's $\alpha$; TC-698).
- The rate's dependence on the frequency is tested for no molecule by a calculation.
  For the atoms it rests on one-electron and helium calculations from γ = 0.2 to 2.8
  (TC-653); at these intensities γ is 1.1 to 1.9. The molecules' default keeps the
  static Coulomb factor; the atoms' default changed to the frequency-dependent one in
  TC-848, and whether the molecules follow is open.

**Guo et al. (1998)** measured N2⁺ beside Ar⁺ and O2⁺ beside Xe⁺ with a 30 fs pulse at
800 nm. Their yields are in figures only and are not taken, but the text states two
comparisons: N2 behaves like argon, and "the O2⁺ LP ion yield is roughly ten times
below Xe⁺ in the unsaturated range of the data". A ratio of yields at one intensity
does not depend on the intensity scale, which Guo et al. set with ADK's rate of argon.
The default's yields of a 30 fs pulse:

| I (W/cm²) | 2e13 | 3e13 | 5e13 | 7e13 | 1e14 | 1.5e14 |
|---|---|---|---|---|---|---|
| N2 / Ar | 1.41 | 1.28 | 1.13 | 1.04 | 0.94 | 0.84 |
| O2 / Xe | 0.035 | 0.039 | 0.045 | 0.055 | 0.089 | 0.235 |

- **N2 against argon agrees.**
- **O2 against xenon is 2 to 3 times lower than stated** (4.5 times with N2's
  $\alpha$), and changes little with the intensity
  where xenon is not saturated, so that no error of the intensity scale moves it.
  Wahlstrand et al.'s O2 against xenon points the same way, by 1.1 to 3. Two
  measurements put the default's O2 too low relative to xenon. Whether O2 is low or
  xenon high is not decided by a ratio: against Wahlstrand et al.'s absolute yields
  xenon is 4 to 8 times high and O2 5.0 to 6.5 (TC-683, TC-698).

**DeWitt et al. (2001), Loriot et al. (2008) and Durá et al. (2012)** measured the same
pairs, and their ratios are curated (`ionization_yield_ratio`, TC-822). They differ in
how the ion detector's efficiency is treated, which moves O2/Xe by a factor of about 2:
Guo et al. and Durá et al. correct for it ($\propto Z/\sqrt{M}$), DeWitt et al. do not,
and Loriot et al.'s measurement is optical.

- **O2 against xenon:** the default is 0.23 to 0.70 of DeWitt et al.'s uncorrected
  ratio (0.47 to 1.4 with Guo et al.'s correction), and three to seven times below
  Durá et al.'s at 3e13 W/cm², whose ratio rises from 0.115 at 0.8 µm to 0.27 at 1.4
  µm where the default's does not change. The absolute measurements put the default's
  xenon high against argon by 2 to 3 (TC-819), which would close most of the gap at
  800 nm but not the dependence on wavelength.
- **N2 against argon:** the measurements disagree, DeWitt et al. 1.05 to 2.03, Loriot et
  al. 0.72 to 0.89, Durá et al. 0.71 to 0.89; the default, 0.84 to 1.39, lies between.

```{figure} /auto_gallery/ionisation/images/sphx_glr_plot_molecular_ionisation_measured_001.*
:target: ../auto_gallery/ionisation/plot_molecular_ionisation_measured.html
:alt: The models against Wahlstrand et al.'s measured yields of N2 and O2

Wahlstrand et al.'s yields of N2 and O2 with the band of their intensity's uncertainty,
the default's and the forms'; the default over the measurement for N2, O2, Ar, Kr
and Xe.
```

## Open questions

- **The dependence on the frequency** of the molecules' rate is tested by no
  calculation; a time-dependent calculation of N2 or H2 at 800 nm would test it.
- **CO2, N2O and SF6** take N2's correction untested; no all-electron static rate
  of them is published. For CO, the one molecule where the transfer can be tested, it
  holds with a constant offset of 1.3 (TC-683).
- **The measured yields of Wahlstrand et al.** lie below the default for the
  molecules as for the atoms (TC-660, TC-680).
- **O2 against xenon**: four measurements put the default's O2 too low relative to
  xenon, by 1.1 to 7 depending on the detection convention, and one (Durá et al.)
  measures a rise with wavelength the default lacks (TC-683, TC-822).
- **SF6's rate**, 80 times N2's at the same intensity on Zhao et al.'s coefficients,
  is tested by nothing: a measured SF6⁺ yield beside an atom's or N2's would decide it
  (TC-684).
- **Structure coefficients.** Those of Zhao et al. (2017) differ from the same
  group's of 2010 by 60 % in the rate for N2 and 50 % for O2, without a new
  comparison with a calculation or a measurement (TC-675). N2O and CH4 have one
  tabulation each (Zhao et al. 2017, 2011).
- **CO2, N2O and SF6** have no calculated reference of any kind; O2 and CH4 have one
  method's (CCSD, without triples, which raise CO's rate by 20 %).
