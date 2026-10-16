# Historical comparison regression tables

These tables were moved from the manual during its documentation review. They
retain numerical regression expectations against the named snapshots in the
`solprop_*.json` fixtures. They do not describe the current downstream package.
The tests check the solprop comparison columns against those fixtures; the other
columns of the accuracy table repeat the manual's and are checked there.

<!-- accuracy-table:start -->
| Gas | $\beta_2$, statistical $1\sigma$ | Largest leave-one-out shift of $\beta_2$ | Dataset left out | Working estimate | solprop's formula against the fit |
|---|---|---|---|---|---|
| He | 0.00040 % | 0.62 % | puchalski_2016/cauchy_total | < 0.01 % | -13.65 % |
| Ne | 0.21 % | 0.34 % | cuthbertson_1932/experiment_Ne | 0.3 % | -0.45 % |
| Ar | 0.016 % | 0.010 % | lesiuk_2023/cauchy_Ar | 0.05 % | -1.30 % |
| Kr | 0.11 % | 0.30 % | hohm_1994b/alpha_Kr | 0.5 % | -1.60 % |
| Xe | 0.090 % | 0.13 % | hohm_1994b/alpha_Xe | 0.4 % | +1.13 % |
| N2 | 0.012 % | 0.17 % | peck_1966/experiment | 0.2 % | -0.87 % |
| H2 | 0.027 % | 0.069 % | peck_1977/experiment | 0.1 % | -0.04 % |
| D2 | 0.016 % | 0.029 % | shelton_1990/density_H2_over_D2 | 0.1 % | +2.88 % |
| O2 | 0.055 % | 0.040 % | smith_1976/experiment_O2 | 0.1 % | +0.51 % |
| CH4 | 0.076 % | 0.078 % | hohm_1993/alpha_CH4 | 0.2 % | +4.35 % |
| SF6 | 0.27 % | 1.3 % | shelton_1985_sf6/density_N2_over_SF6 | 2 % | +54.82 % |
| N2O | 0.66 % | 1.7 % | cuthbertson_1914/n2o_dispersion | 3 % | +11.22 % |
| CO2 | 0.054 % | 0.21 % | old_1971/experiment | 0.2 % | n/a |
| O3 | 0.52 % | 1.1 % | cuthbertson_1914/ozone_adopted | 50 % | n/a |
| CO | 0.21 % | 0.29 % | cuthbertson_1920/experiment_CO | 0.3 % | n/a |
<!-- accuracy-table:end -->

## Transport

<!-- transport-solprop:start -->
solprop's values themselves, recorded at solprop d577414
(`tests/data/solprop_transport.json`), against optmat's: the largest relative
difference in the viscosity over solprop's range, and in the heat-capacity ratio
$\gamma = c_p/(c_p - R)$. The molar masses agree exactly.

| Gas | solprop's range, K | Viscosity | $\gamma$ |
|---|---|---:|---:|
| He | 200-600 | −4.9e-4 at 200 K | 2.3e-7 |
| Ne | 200-600 | −2.9e-4 at 200 K | 2.0e-7 |
| Ar | 200-600 | −1.6e-3 at 200 K | 8.0e-6 |
| Kr | 200-600 | −3.3e-3 at 200 K | 1.0e-5 |
| Xe | 200-600 | −5.4e-3 at 200 K | 2.9e-5 |
| N2 | 200-600 | −1.7e-3 at 200 K | 4.0e-6 |
| H2 | 200-600 | −2.3e-4 at 200 K | 7.2e-7 |
| D2 | 200-600 | −6.8e-3 at 600 K | 1.8e-7 |
| O2 | 200-600 | −1.6e-3 at 200 K | 5.3e-6 |
| CH4 | 200-600 | −1.2e-3 at 250 K | 4.5e-6 |
| SF6 | 225-600 | +4.7e-3 at 225 K | 3.0e-6 |
| N2O | 200-500 | −5.2e-3 at 200 K | 1.1e-5 |

solprop's heat-capacity ratio is CoolProp's at 1e-3 bar, not the ideal-gas part, hence
the 1e-5 differences for the heavier gases.
<!-- transport-solprop:end -->



<!-- nonlinear-solprop:start -->
| Gas | solprop's source | solprop gamma_T (au) | / Kerr, 400 nm | / Kerr, 800 nm | / Kerr, 1800 nm | / THG, 800 nm |
|---|---|---:|---:|---:|---:|---:|
| He | lehmeier | 43.53 | 0.897 | 0.981 | 1.003 | 0.925 |
| Ne | lehmeier | 78.35 | 0.633 | 0.706 | 0.727 | 0.657 |
| Ar | lehmeier | 1021 | 0.631 | 0.812 | 0.862 | 0.690 |
| Kr | lehmeier | 2777 | 0.683 | 0.971 | 1.049 | 0.773 |
| Xe | lehmeier | 8132 | 0.626 | 1.047 | 1.164 | 0.755 |
| H2 | shelton | 686.5 | 0.697 | 0.921 | 0.986 | 0.768 |
| D2 | shelton | 686.5 | 0.720 | 0.947 | 1.013 | 0.792 |
| N2 | lehmeier | 917.6 | 0.735 | 0.931 | 0.984 | 0.800 |
| O2 | zahedpour | 975.4 | 0.647 | 0.932 | 0.999 | 0.739 |
| CO2 | none | n/a | n/a | n/a | n/a | n/a |
| CH4 | shelton | 2689 | 0.628 | 0.933 | 1.018 | 0.721 |
| SF6 | shelton | 1404 | 0.781 | 0.945 | 0.976 | 0.841 |
| CO | none | n/a | n/a | n/a | n/a | n/a |
| N2O | wahlstrand | 2062 | 0.782 | 1.051 | 1.125 | 0.867 |
| O3 | none | n/a | n/a | n/a | n/a | n/a |
<!-- nonlinear-solprop:end -->


Raman lines, solprop (`tests/data/solprop_raman.json`, TC-439) against the curated
values. The vibrational amplitude compares solprop's
$(4\pi\varepsilon_0\,\partial\alpha/\partial q)^2 / (2\mu\Omega\varepsilon_0)$ with the
curated $|\langle 0|\bar\alpha|1\rangle|^2/(\hbar\varepsilon_0)$.

<!-- raman-solprop:start -->
| Gas | B: solprop / optmat (cm-1) | D: solprop / optmat (cm-1) | Delta alpha: solprop / optmat - 1 | nu: solprop / optmat (cm-1) | Vibrational A: solprop / optmat - 1 |
|---|---|---|---|---|---|
| N2 | 1.9900 / 1.9896 | 5.740e-06 / 5.760e-06 | +0.0% | 2330.00 / 2329.92 | +0.4% |
| O2 | 1.4377 / 1.4377 | 4.839e-06 / 4.839e-06 | +0.0% | 1556.39 / 1556.39 | -1.4% |
| H2 | 58.9000 / 59.3392 | 5.000e-02 / 4.599e-02 | -4.5% | 4155.10 / 4160.20 | +13.1% |
| D2 | 29.3000 / 29.9105 | 2.100e-02 / 1.134e-02 | +0.8% | 2987.00 / 2993.50 | +31.1% |
| CH4 |  |  |  | 2914.00 / 2916.70 | -10.4% |
| SF6 |  |  |  | 775.00 / 775.00 | -11.0% |
<!-- raman-solprop:end -->


Raman responses, solprop (`tests/data/solprop_raman_response.json`, solprop 2dc84d0,
295 K, 1 amagat) against optmat's, part by part. The early amplitude is the
least-squares scale of optmat's response on solprop's over the first 100 fs, before any
line has dephased; the late one compares the two responses' largest magnitudes near
20 ps.

<!-- response-solprop:start -->
| Gas | Part | Early amplitude: optmat / solprop | At 20 ps: optmat / solprop |
|---|---|---|---|
| N2 | rotation | 1.000 | 0.958 |
| N2 | vibration | 1.033 | 5.47 |
| O2 | rotation | 1.000 | 1.02 |
| O2 | vibration | 1.081 | 0.155 |
| H2 | rotation | 1.101 | 1.24 |
| H2 | vibration | 0.897 | 0.688 |
| D2 | rotation | 0.981 | 1.21 |
| D2 | vibration | 0.772 | 0.374 |
| CH4 | vibration | 1.116 | 1.07 |
| SF6 | vibration | 1.124 | 1.14 |
<!-- response-solprop:end -->


<!-- ionisation-rate-legacy:start -->
| Gas | I (W/cm2) | gamma | rate: solprop / Luna | optmat / Luna | optmat / solprop | 30 fs yield: Luna | solprop | optmat |
|---|---|---|---|---|---|---|---|---|
| He | 5.00e+13 | 2.03 | 0.89 | 1.62 | 1.82 | 2.8e-12 | 2.5e-12 | 4.7e-12 |
| He | 1.00e+14 | 1.43 | 0.84 | 1.24 | 1.47 | 8.5e-09 | 7.2e-09 | 1.1e-08 |
| He | 2.00e+14 | 1.01 | 0.79 | 0.99 | 1.26 | 7.4e-06 | 5.9e-06 | 7.6e-06 |
| He | 5.00e+14 | 0.64 | 0.68 | 0.79 | 1.16 | 0.0073 | 0.0051 | 0.006 |
| Ne | 5.00e+13 | 1.90 | 0.81 | 0.50 | 0.62 | 3e-09 | 2.5e-09 | 1.6e-09 |
| Ne | 1.00e+14 | 1.34 | 0.74 | 0.38 | 0.52 | 2.8e-06 | 2.1e-06 | 1.1e-06 |
| Ne | 2.00e+14 | 0.95 | 0.65 | 0.30 | 0.46 | 0.00084 | 0.00056 | 0.00026 |
| Ne | 5.00e+14 | 0.60 | 0.51 | 0.24 | 0.47 | 0.22 | 0.13 | 0.061 |
| Ar | 5.00e+13 | 1.62 | 0.62 | 1.13 | 1.81 | 9.5e-05 | 6.1e-05 | 0.00011 |
| Ar | 1.00e+14 | 1.15 | 0.51 | 0.97 | 1.89 | 0.0081 | 0.0044 | 0.0081 |
| Ar | 2.00e+14 | 0.81 | 0.39 | 0.92 | 2.36 | 0.24 | 0.11 | 0.22 |
| Ar | 5.00e+14 | 0.51 | 0.23 | 1.03 | 4.56 | 1 | 0.83 | 1 |
| Kr | 5.00e+13 | 1.53 | 0.53 | 0.67 | 1.26 | 0.0039 | 0.0022 | 0.0029 |
| Kr | 1.00e+14 | 1.08 | 0.41 | 0.41 | 1.01 | 0.21 | 0.096 | 0.1 |
| Kr | 2.00e+14 | 0.77 | 0.28 | 0.25 | 0.88 | 1 | 0.86 | 0.84 |
| Kr | 5.00e+14 | 0.48 | 0.13 | 0.11 | 0.82 | 1 | 1 | 1 |
| Xe | 5.00e+13 | 1.42 | 0.40 | 0.64 | 1.58 | 0.09 | 0.04 | 0.068 |
| Xe | 1.00e+14 | 1.01 | 0.28 | 0.36 | 1.29 | 0.94 | 0.59 | 0.71 |
| Xe | 2.00e+14 | 0.71 | 0.16 | 0.19 | 1.16 | 1 | 1 | 1 |
| Xe | 5.00e+14 | 0.45 |  |  | 0.85 |  | 1 | 1 |
<!-- ionisation-rate-legacy:end -->


<!-- ionisation-rate-legacy-tests:start -->
| Test | Gas | Luna | solprop | optmat |
|---|---|---|---|---|
| Wallace et al. 2016: intensity factor (calibration 1.3 %) | Ar | 1.05 | 0.94 | 1.03 |
| Wallace et al. 2016: intensity factor (calibration 1.5 %) | Kr | 0.95 | 1.06 | 1.02 |
| Wallace et al. 2016: intensity factor (calibration 2.5 %) | Xe | 0.90 | 1.13 | 1.03 |
| `nambu_2026/fig6_Ar_LP`: intensity factor | Ar | 1.08 | 1.19 | 1.06 |
| `nambu_2026/fig6_Kr_LP`: intensity factor | Kr | 0.85 | 0.99 | 0.97 |
| `nambu_2026/fig6_Xe_LP`: intensity factor | Xe | 0.88 | 1.03 | 0.92 |
| `nie_2022/fig3_He`: intensity factor | He | 0.90 | 0.93 | 0.91 |
| `nie_2022/fig4_Ar_LP`: intensity factor | Ar | 0.98 | 1.12 | 0.98 |
| `nie_2022/fig4_Kr_LP`: intensity factor | Kr | 0.80 | 0.92 | 0.89 |
| `nie_2022/fig4_Xe_LP`: intensity factor | Xe | 0.75 | 0.89 | 0.80 |
| Wahlstrand et al. 2018: rate over measured | Ar | 0.97 to 3.10 | 0.64 to 1.60 | 1.20 to 3.06 |
| Wahlstrand et al. 2018: rate over measured | Kr | 1.44 to 7.31 | 0.88 to 3.53 | 1.36 to 3.80 |
| Wahlstrand et al. 2018: rate over measured | Xe | 1.81 to 9.37 | 1.09 to 4.18 | 2.90 to 6.75 |
| Parker et al. 2007 (780 nm, two electrons, to 1e15 W/cm2): rate over theirs, the codes at 800 nm | He | 0.59 to 1.44 | 0.50 to 0.88 | 0.76 to 1.05 (at 780 nm: 0.85 to 1.06) |
<!-- ionisation-rate-legacy-tests:end -->


<!-- molecular-rate-solprop:start -->
| Molecule | I_p, optmat (eV) | I_p, solprop (eV) | Difference (meV) | Rate ratio at 2e10 V/m | at 4e10 V/m |
|---|---|---|---|---|---|
| H2 | 15.4258 | 15.4261 | -0.3 | 1.001 | 1.000 |
| D2 | 15.4665 | 15.4670 | -0.5 | 1.001 | 1.001 |
| N2 | 15.5807 | 15.5812 | -0.5 | 1.001 | 1.001 |
| O2 | 12.0696 | 12.0697 | -0.1 | 1.000 | 1.000 |
| CH4 | 12.6100 | 12.6152 | -5.2 | 1.010 | 1.005 |
| SF6 | 15.3500 | 15.3200 | +30.0 | 0.938 | 0.968 |
| N2O | 12.8860 | 12.8982 | -12.2 | 1.024 | 1.012 |
<!-- molecular-rate-solprop:end -->


<!-- ionisation-solprop:start -->
| Gas | solprop I_p (eV) | - optmat (meV) | alpha: solprop / optmat - 1 | alpha+: solprop / optmat - 1 | solprop (n, l, Z) |
|---|---|---|---|---|---|
| He | 24.5882 | +0.8 | -4.6% | -0.0% | (1, 0, 1) |
| Ne | 21.5650 | +0.5 | -10.7% | -5.3% | (2, 1, 1) |
| Ar | 15.7608 | +1.2 | -2.8% | -0.4% | (3, 1, 1) |
| Kr | 13.9921 | -7.5 |  |  | (4, 1, 1) |
| Xe | 12.1308 | +1.0 |  |  | (5, 1, 1) |
| H2 | 15.4261 | +0.3 |  |  |  |
| D2 | 15.4670 | +0.5 |  |  |  |
| N2 | 15.5812 | +0.5 |  |  | (2, 0, 0.9) |
| O2 | 12.0697 | +0.1 |  |  | (2, 0, 0.53) |
| CH4 | 12.6152 | +5.2 |  |  |  |
| SF6 | 15.3200 | -30.0 |  |  |  |
| N2O | 12.8982 | +12.2 |  |  |  |
<!-- ionisation-solprop:end -->


<!-- crystals:start -->
| Material | Axis | Model | Range (nm) | Sample | Index relative to | n at 500 nm | n at 1 µm | solprop, at 1 µm |
|---|---|---|---|---|---|---|---|---|
| CaF2 |  | `daimon_2002/sellmeier_20C` | 138-2326 | 20 °C | nitrogen, 1013.25 hPa | 1.43690 | 1.42931 | −3.9e-4 |
| MgF2 | o | `dodge_1984/sellmeier_o` | 203-7042 | 19 °C | room air (inferred) | 1.38016 | 1.37396 | −3.7e-4 |
| MgF2 | e | `dodge_1984/sellmeier_e` | 203-7042 | 19 °C | room air (inferred) | 1.39207 | 1.38557 | none |
| BaF2 |  | `malitson_1964/sellmeier` | 265-10349 | 25 °C | room air (inferred) | 1.47824 | 1.46895 | −3.9e-4 |
| Al2O3 | o | `weber_2003/sellmeier_o` | 220-5001 | 20 °C | room air (inferred) | 1.77475 | 1.75616 | −4.8e-4 |
| Al2O3 | e | `weber_2003/sellmeier_e` | 220-5001 | 20 °C | room air (inferred) | 1.76655 | 1.74828 | none |
| N-BK7 |  | `schott_2023/sellmeier` | 313-2326 | 20 °C, catalogue | air, 1013.3 hPa | 1.52184 | 1.50791 | −4.1e-4 |
<!-- crystals:end -->
