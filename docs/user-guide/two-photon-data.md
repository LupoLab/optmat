# Two-photon ionisation: source coverage and comparisons

Start with {doc}`/two_photon` for the cross-section convention and model interface. This chapter describes the recorded evidence and its limits; source details are in the {doc}`/the-data/sources` catalogue.

## The curated two-photon ionisation cross sections

Every curve read from a figure is `use = "check"` (data rule 8): those read from a
vector drawing by their script, those read from a scanned image (*raster*) by theirs.
The tables list every two-photon ionisation dataset of each gas, with the photon
energies it covers inside the window $[I_p/2, E_1)$. The comparisons give the total
cross section, the sum of the partials where a source draws only those, at photon
energies off resonance. Each value is followed by its ratio to the first row.

## Helium

<!-- two-photon-coverage-He:start -->
| Dataset | Form | Use | Polarisation | Channel | Photon energies in the window (eV) |
|---|---|---|---|---|---|
| `bachau_1990/fig1a_W2_L2` | raster | check | linear | total | 12.339 to 14.962 |
| `bachau_1990/tbl1_sigma2_L2` | table | check | linear | total | 13.606 |
| `bachau_1990/tbl1_sigma2_qd` | table | check | linear | total | 13.606 |
| `kulander_1987/tbl1b_sigma2_tdhf` | table | check | linear | total | 13.600 to 17.700 |
| `lhuillier_1986/fig11_one_electron` | raster | check | linear | total | 12.406 to 21.143 |
| `lhuillier_1986/fig11_rpa` | raster | check | linear | total | 12.406 to 21.206 |
| `mihelic_2021/fig3a_sigma2_L0` | vector | check | linear | total | 12.300 to 21.192 |
| `mihelic_2021/fig3a_sigma2_L2` | vector | check | linear | total | 12.300 to 21.204 |
| `mihelic_2021/fig3b_sigma2_L0` | vector | check | linear | total | outside |
| `mihelic_2021/fig3b_sigma2_L2` | vector | check | linear | total | outside |
| `mihelic_2021/fig5_sigma2_L0_1ses` | vector | check | linear | 1s eps s | outside |
| `mihelic_2021/fig5_sigma2_L0_2pep` | vector | check | linear | 2p eps p | outside |
| `mihelic_2021/fig5_sigma2_L0_2ses` | vector | check | linear | 2s eps s | outside |
| `mihelic_2021/fig5_sigma2_L0` | vector | check | linear | total | outside |
| `mihelic_2021/fig5_sigma2_L2_1sed` | vector | check | linear | 1s eps d | outside |
| `mihelic_2021/fig5_sigma2_L2_2pef` | vector | check | linear | 2p eps f | outside |
| `mihelic_2021/fig5_sigma2_L2_2pep` | vector | check | linear | 2p eps p | outside |
| `mihelic_2021/fig5_sigma2_L2_2sed` | vector | check | linear | 2s eps d | outside |
| `mihelic_2021/fig5_sigma2_L2` | vector | check | linear | total | outside |
| `nikolopoulos_2001/fig3_sigma2_D` | vector | check | linear | total | 12.748 to 21.200 |
| `nikolopoulos_2001/fig4_sigma2_D` | vector | check | linear | total | 12.350 to 21.067 |
| `nikolopoulos_2001/fig4_sigma2_S` | vector | check | linear | total | 12.350 to 21.167 |
| `nikolopoulos_2001/fig5_sigma2_double` | vector | check | linear | total | outside |
| `nikolopoulos_2001/fig5_sigma2_single` | vector | check | linear | total | outside |
| `nikolopoulos_2001/text_sigma2_single_45eV` | text | check | linear | total | outside |
| `proulx_1993/tbl1_rate2` | table | check | linear | total | 14.327 to 20.449 |
| `proulx_1993/tbl2_sigma2` | table | check | linear | total | outside |
| `ritchie_1977/tbl1_E_nonorth` | table | check | linear | total | 12.408 to 21.198 |
| `ritchie_1977/tbl1_E` | table | check | linear | total | 12.408 to 21.198 |
| `ritchie_1977/tbl1_NE_nonorth` | table | check | linear | total | 12.408 to 21.198 |
| `ritchie_1977/tbl1_NE` | table | check | linear | total | 12.408 to 21.198 |
| `ritchie_1977/tbl2_E_nonorth` | table | check | linear | total | 12.408 to 20.762 |
| `ritchie_1977/tbl2_E` | table | check | linear | total | 12.408 to 20.762 |
| `ritchie_1977/tbl2_NE_nonorth` | table | check | linear | total | 12.408 to 20.762 |
| `ritchie_1977/tbl2_NE` | table | check | linear | total | 12.408 to 20.762 |
| `ritchie_1977/text_length_1000A` | text | check | linear | total | 12.408 |
| `saenz_1999/fig1a_sigma2_1S` | vector | check | linear | total | 12.651 to 21.214 |
| `saenz_1999/fig1b_sigma2_1D` | vector | check | linear | total | 12.650 to 21.217 |
| `saenz_1999/fig4_sigma2_1D_length` | vector | check | linear | total | 12.650 to 21.212 |
| `saenz_1999/fig4_sigma2_1D_velocity` | vector | check | linear | total | 12.650 to 21.217 |
| `saenz_1999/text_sigma2_1D_complex_scaling` | text | check | linear | total | 13.500 |
| `saenz_1999/text_sigma2_1D_length` | text | check | linear | total | 13.500 |
| `saenz_1999/text_sigma2_1D_velocity` | text | check | linear | total | 13.500 |
| `saenz_1999/text_sigma2_1S_complex_scaling` | text | check | linear | total | 13.500 |
| `saenz_1999/text_sigma2_1S_length` | text | check | linear | total | 13.500 |
| `saenz_1999/text_sigma2_1S_velocity` | text | check | linear | total | 13.500 |
| `sato_2011/tbl2_lopt` | table | check | linear | total | 20.193 |
| `sato_2011/tbl2_obs` | table | check | linear | total | 20.193 |
| `sato_2011/tbl2_rmf` | table | check | linear | total | 20.193 |
| `sato_2011/text_1e10` | text | check | linear | total | 20.193 |
| `travers_2026_trecx_two_photon/sigma2_hacc_He` | table | check | linear | total | 13.000 to 20.000 |
| `travers_2026_trecx_two_photon/sigma2_sae_He` | table | check | linear | total | 13.000 to 20.000 |
| `vanderhart_2005/fig2_sigma2_rmf` | vector | check | linear | total | 12.324 to 21.215 |
| `vanderhart_2005/fig2_sigma2_saenz` | vector | check | linear | total | 12.651 to 21.215 |
<!-- two-photon-coverage-He:end -->

<!-- two-photon-comparison-He:start -->
| $\sigma^{(2)}$ (cm⁴ s) at | 13 eV | 14 eV | 16 eV | 18 eV | 20 eV |
|---|---|---|---|---|---|
| Saenz & Lambropoulos 1999, ¹S + ¹D (velocity) | 1.48e-51 | 1.43e-51 | 1.09e-51 | 8.24e-52 | 7.75e-52 |
| Mihelič & Horvat 2021, L = 0 + 2 | 1.32e-51 (0.90) | 1.27e-51 (0.89) | 1.04e-51 (0.96) | 8.55e-52 (1.04) | 8.74e-52 (1.13) |
| Nikolopoulos & Lambropoulos 2001, S + D (Fig. 4) | 1.49e-51 (1.01) | 1.31e-51 (0.92) | 9.83e-52 (0.90) | 6.97e-52 (0.85) | 6.98e-52 (0.90) |
| van der Hart & Bingham 2005, R-matrix Floquet (1e13 W/cm²) | 1.37e-51 (0.93) | 1.31e-51 (0.92) | 1.07e-51 (0.98) | 8.51e-52 (1.03) | 8.28e-52 (1.07) |
<!-- two-photon-comparison-He:end -->

The three calculations agree within about 15 % from threshold to 20 eV. Mihelič and
Horvat's own curves and Saenz and Lambropoulos's differ the most near threshold; there
Saenz's scaled axis matters (TC-873, TC-874). Sato et al.'s FEL measurement (TC-876)
has one point in the window, 2.6(6)×10⁻⁵² cm⁴ s at 61.4 nm (20.19 eV, 1e13 W cm⁻²),
three times below the calculations. Its other points lie above $E_1$, one on the 1s2p
resonance. Sato's LOPT column is Nikolopoulos and Lambropoulos's D wave alone (their
Fig. 3), not the S + D total (TC-875, TC-876). Van der Hart and Bingham's R-matrix
Floquet calculation at 1e13 W cm⁻² (TC-904), which Sato quotes, is 0.92 to 1.07 of
Saenz's. Bachau et al.'s and Proulx and Shakeshaft's tables (TC-905, TC-906) lie above
the window.

## Neon

<!-- two-photon-coverage-Ne:start -->
| Dataset | Form | Use | Polarisation | Channel | Photon energies in the window (eV) |
|---|---|---|---|---|---|
| `mcguire_1981/fig1_sigma2_Ne` | raster | check | circular | total | 11.069 to 16.553 |
| `mckenna_2004/fig1_sigma2_Ne` | vector | check | linear | total | 10.800 to 16.667 |
| `mckenna_2004/fig2_sigma2_Ne` | vector | check | linear | total | outside |
| `moccia_1983/fig1_sigma2_L_Ne` | raster | check | linear | total | 11.976 to 15.777 |
| `moccia_1983/fig1_sigma2_V_Ne` | raster | check | linear | total | 11.157 to 16.193 |
<!-- two-photon-coverage-Ne:end -->

<!-- two-photon-comparison-Ne:start -->
| $\sigma^{(2)}$ (cm⁴ s) at | 11 eV | 12 eV | 14 eV | 16 eV |
|---|---|---|---|---|
| McKenna & van der Hart 2004 | 1.82e-51 | 1.71e-51 | 1.56e-51 | 1.6e-51 |
| Moccia, Rahman & Rizzo 1983, velocity (raster) | — | 1.05e-51 (0.62) | 5.69e-52 (0.37) | 2.21e-52 (0.14) |
| Moccia, Rahman & Rizzo 1983, length (raster) | — | 3.28e-52 (0.19) | 1.23e-52 (0.08) | 3.27e-53 (0.02) |
| McGuire 1981, circular (raster) | — | 1.01e-51 (0.59) | 1.16e-51 (0.75) | 1.53e-51 (0.96) |
<!-- two-photon-comparison-Ne:end -->

One calculation from a vector drawing, McKenna and van der Hart's. The only
independent one is Moccia, Rahman and Rizzo's perturbative calculation (raster;
TC-910). It is strongly gauge dependent: its velocity form is 0.62 to 0.14 of McKenna's
from 12 to 16 eV, its length form 0.19 to 0.02. McGuire's is for circular polarisation,
read from a scan with a piecewise axis.

## Argon

<!-- two-photon-coverage-Ar:start -->
| Dataset | Form | Use | Polarisation | Channel | Photon energies in the window (eV) |
|---|---|---|---|---|---|
| `fink_1990/fig3_circular` | raster | check | circular | total | 7.997 to 11.578 |
| `fink_1990/fig3_linear` | raster | check | linear | total | 8.006 to 11.581 |
| `mcguire_1981/fig2_sigma2_Ar` | raster | check | circular | total | 8.919 to 11.373 |
| `mckenna_2004/fig5_sigma2_Ar` | vector | check | linear | total | 7.912 to 11.623 |
| `mckenna_2004/fig7_sigma2_Ar` | vector | check | linear | total | outside |
| `moccia_1983/fig2_sigma2_L_Ar` | raster | check | linear | total | 8.027 to 11.562 |
| `moccia_1983/fig2_sigma2_V_Ar` | raster | check | linear | total | 7.894 to 11.562 |
| `pan_1990b/fig1b_all_L_p1S_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1S | 8.101 to 10.899 |
| `pan_1990b/fig1b_all_V_p1S_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1S | 8.199 to 11.002 |
| `pan_1990b/fig1b_noshakeup_L_p1S_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1S | 8.301 to 11.299 |
| `pan_1990b/fig1b_noshakeup_V_p1S_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1S | 8.098 to 9.299 |
| `pan_1990b/fig2_all_L_p1D_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1D | 8.102 to 9.601 |
| `pan_1990b/fig3_all_L_f1D_Ar` | raster | check | linear | 3p5(2P) epsilon-f 1D | 8.199 to 11.600 |
| `pan_1990b/fig4a_hf_L_Ar` | raster | check | linear | total (3p5 epsilon-p 1S + epsilon-p 1D + epsilon-f 1D) | 8.101 to 11.198 |
| `pan_1990b/fig4a_hf_V_Ar` | raster | check | linear | total (3p5 epsilon-p 1S + epsilon-p 1D + epsilon-f 1D) | 8.099 to 11.602 |
| `pan_1990b/fig4b_noshakeup_L_Ar` | raster | check | linear | total (3p5 epsilon-p 1S + epsilon-p 1D + epsilon-f 1D) | 8.098 to 11.301 |
| `pan_1990b/fig4b_noshakeup_V_Ar` | raster | check | linear | total (3p5 epsilon-p 1S + epsilon-p 1D + epsilon-f 1D) | 8.101 to 11.500 |
| `pan_1990b/fig4c_all_L_Ar` | raster | check | linear | total (3p5 epsilon-p 1S + epsilon-p 1D + epsilon-f 1D) | 8.401 to 11.098 |
| `pan_1990b/fig4c_all_V_Ar` | raster | check | linear | total (3p5 epsilon-p 1S + epsilon-p 1D + epsilon-f 1D) | 8.098 to 11.098 |
| `pan_1990b/fig5_all_L_circular_Ar` | raster | check | circular | total (3p5 epsilon-p 1D + epsilon-f 1D; no 1S for circular) | 8.099 to 10.002 |
| `pan_1990b/fig5_all_V_circular_Ar` | raster | check | circular | total (3p5 epsilon-p 1D + epsilon-f 1D; no 1S for circular) | 9.701 to 11.199 |
| `pan_1990b/fig5_hf_L_circular_Ar` | raster | check | circular | total (3p5 epsilon-p 1D + epsilon-f 1D; no 1S for circular) | 8.299 to 11.500 |
| `pan_1990b/fig5_hf_V_circular_Ar` | raster | check | circular | total (3p5 epsilon-p 1D + epsilon-f 1D; no 1S for circular) | 9.601 to 10.002 |
| `petrov_2016/fig10_corr_L` | vector | check | linear | total | 8.000 to 11.617 |
| `petrov_2016/fig10_corr_V` | vector | check | linear | total | 8.000 to 11.617 |
| `petrov_2016/fig10_corrcp_L` | vector | check | linear | total | 8.000 to 11.619 |
| `petrov_2016/fig10_corrcp_V` | vector | check | linear | total | 8.000 to 11.623 |
| `petrov_2016/fig6_corr_L_ep1S` | vector | check | linear | 3p5 ep (1S) | 8.000 to 11.619 |
| `petrov_2016/fig6_corr_V_ep1S` | vector | check | linear | 3p5 ep (1S) | 8.000 to 11.621 |
| `petrov_2016/fig6_corrcp_L_ep1S` | vector | check | linear | 3p5 ep (1S) | 8.000 to 11.430 |
| `petrov_2016/fig6_corrcp_V_ep1S` | vector | check | linear | 3p5 ep (1S) | 8.000 to 11.466 |
| `petrov_2016/fig7_corr_L_ep1D` | vector | check | linear | 3p5 ep (1D) | 8.000 to 11.623 |
| `petrov_2016/fig7_corrcp_L_ep1D` | vector | check | linear | 3p5 ep (1D) | 8.000 to 11.319 |
| `petrov_2016/fig8_corr_L_ef1D` | vector | check | linear | 3p5 ef (1D) | 8.000 to 11.617 |
| `petrov_2016/fig8_corrcp_L_ef1D` | vector | check | linear | 3p5 ef (1D) | 8.000 to 11.619 |
| `pindzola_1975/fig2_L_p1D_Ar` | raster | check | linear | 3p5 kp 1D | 8.000 to 9.001 |
| `pindzola_1975/fig2_V_p1D_Ar` | raster | check | linear | 3p5 kp 1D | 8.001 to 9.001 |
| `pindzola_1975/fig3_L_p1S_Ar` | raster | check | linear | 3p5 kp 1S | 7.998 to 11.002 |
| `pindzola_1975/fig3_V_p1S_Ar` | raster | check | linear | 3p5 kp 1S | 7.997 to 10.500 |
| `pindzola_1975/fig4_L_f1D_Ar` | raster | check | linear | 3p5 kf 1D | 8.001 to 11.502 |
| `pindzola_1975/fig4_V_f1D_Ar` | raster | check | linear | 3p5 kf 1D | 7.998 to 11.502 |
| `pindzola_1975/fig6_LGI_f1D_Ar` | raster | check | linear | 3p5 kf 1D | 9.498 to 10.998 |
| `pindzola_1975/fig6_VGI_f1D_Ar` | raster | check | linear | 3p5 kf 1D | 9.498 to 11.502 |
| `starace_1987/fig2_esi_f1D_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-f 1D | 8.498 to 11.501 |
| `starace_1987/fig2_gsc_f1D_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-f 1D | 8.498 to 11.599 |
| `starace_1987/fig2_hf_f1D_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-f 1D | 8.499 to 11.602 |
| `starace_1987/fig2_iic_f1D_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-f 1D | 8.601 to 11.498 |
| `starace_1987/fig3_esi_p1D_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1D | 8.399 to 9.702 |
| `starace_1987/fig3_hf_p1D_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1D | 8.202 to 10.398 |
| `starace_1987/fig3_iic_p1D_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1D | 8.401 to 10.301 |
| `starace_1987/fig4_esi_p1S_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1S | 8.398 to 11.101 |
| `starace_1987/fig4_gsc_p1S_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1S | 8.402 to 10.700 |
| `starace_1987/fig4_hf_p1S_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1S | 8.202 to 10.702 |
| `starace_1987/fig4_iic_p1S_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-p 1S | 8.402 to 10.798 |
| `starace_1987/fig5_esi_shifted_f1D_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-f 1D | 8.401 to 11.402 |
| `starace_1987/fig5_iic_shifted_f1D_linear_Ar` | raster | check | linear | 3p5(2P) epsilon-f 1D | 8.799 to 9.902 |
| `starace_1987/fig6_esi_circular_Ar` | raster | check | circular | total | 8.100 to 11.100 |
| `starace_1987/fig6_iic_circular_Ar` | raster | check | circular | total | 8.100 to 11.003 |
| `starace_1987/fig7_esi_linear_Ar` | raster | check | linear | total | 8.098 to 9.400 |
| `starace_1987/fig7_iic_linear_Ar` | raster | check | linear | total | 8.399 to 10.799 |
| `travers_2026_trecx_two_photon/sigma2_sae_Ar` | table | check | linear | total | 8.000 to 11.000 |
<!-- two-photon-coverage-Ar:end -->

<!-- two-photon-comparison-Ar:start -->
| $\sigma^{(2)}$ (cm⁴ s) at | 8 eV | 9 eV | 10 eV | 11 eV |
|---|---|---|---|---|
| McKenna & van der Hart 2004 | 1.27e-50 | 1.24e-50 | 1.17e-50 | 1.13e-50 |
| Petrov et al. 2016, corr + CP, length | 1.15e-50 (0.91) | 1.14e-50 (0.92) | 1.07e-50 (0.92) | 1.18e-50 (1.05) |
| Petrov et al. 2016, corr + CP, velocity | 1.26e-50 (0.99) | 1.19e-50 (0.95) | 1.09e-50 (0.93) | 1.06e-50 (0.94) |
| Petrov et al. 2016, corr, length | 1.02e-50 (0.81) | 9.83e-51 (0.79) | 8.83e-51 (0.75) | 8.02e-51 (0.71) |
| Starace & Jiang 1987, IIC, length (raster) | — | 1.3e-50 (1.04) | 1.32e-50 (1.13) | 1.43e-50 (1.27) |
| Starace & Jiang 1987, ESI, length (raster) | — | 1.04e-50 (0.84) | 9.26e-51 (0.79) | 8.14e-51 (0.72) |
| Pan, Gao & Starace 1990, length (raster) | — | 9.21e-51 (0.74) | 8.09e-51 (0.69) | 7.06e-51 (0.63) |
| Pan, Gao & Starace 1990, velocity (raster) | — | 8.77e-51 (0.71) | 8.31e-51 (0.71) | 7.53e-51 (0.67) |
| Pindzola & Kelly 1975, second order, length (raster) | — | 1.17e-50 (0.94) | 1.19e-50 (1.01) | 1.15e-50 (1.02) |
| Moccia, Rahman & Rizzo 1983, velocity (raster) | 7.77e-51 (0.61) | 8.57e-51 (0.69) | 9.04e-51 (0.77) | 9.16e-51 (0.81) |
| Moccia, Rahman & Rizzo 1983, length (raster) | — | 1.09e-50 (0.87) | 9.34e-51 (0.80) | 8.49e-51 (0.75) |
| Fink & Johnson 1990, linear (raster) | — | 1.03e-50 (0.83) | 8.77e-51 (0.75) | 7.42e-51 (0.66) |
| McGuire 1981, circular (raster) | — | 1.53e-50 (1.23) | 1.72e-50 (1.47) | 2.1e-50 (1.87) |
<!-- two-photon-comparison-Ar:end -->

The two all-electron calculations agree within 10 % from 8 to 11 eV: Petrov et al. with
core polarisation (both gauges), and McKenna and van der Hart. Without core polarisation,
Petrov et al.'s is 0.71 to 0.81 of it (TC-877). Starace and Jiang's two
transition-matrix approximations (raster) bracket them, 0.72 to 1.27 (TC-878).
McGuire's circular cross section is 1.2 to 1.9 times the linear calculations, more than
his own ratio $R$ allows (at most 3/2; TC-879).

Three effects together are therefore about ±10 % in argon below 11 eV:
- independent calculation;
- gauge;
- core polarisation.

Without core polarisation, or in the simpler approximations, the spread is about ±30 %.
The older many-body calculations, all raster, lie mostly below the two all-electron
ones:
- Pan, Gao and Starace (0.63 to 0.74; TC-908);
- Moccia, Rahman and Rizzo (0.61 to 0.87; TC-910);
- Fink and Johnson's relativistic calculation (0.66 to 0.83; TC-912);
- Pindzola and Kelly's second order (0.94 to 1.02; TC-909).

No measured argon two-photon ionisation cross section was found. Xenakis et al. observed
argon two-photon ionisation at 82.9 nm but measured no cross section (TC-911).

## Krypton

<!-- two-photon-coverage-Kr:start -->
| Dataset | Form | Use | Polarisation | Channel | Photon energies in the window (eV) |
|---|---|---|---|---|---|
| `fink_1990/fig4_circular` | raster | check | circular | total | 7.347 to 9.943 |
| `fink_1990/fig4_linear` | raster | check | linear | total | 7.616 to 9.946 |
| `lhuillier_1989/fig10_total_Kr` | raster | check | linear | total | 7.470 to 8.648 |
| `lhuillier_1989/fig6_J0_Kr` | raster | check | linear | J = 0 final states (partial) | 7.130 to 7.208 |
| `lhuillier_1989/fig6_total_Kr` | raster | check | linear | total (J = 0 and J = 2 final states) | 7.130 to 7.221 |
| `mcguire_1981/fig3_sigma2_Kr` | raster | check | circular | total | 7.749 to 9.918 |
<!-- two-photon-coverage-Kr:end -->

<!-- two-photon-comparison-Kr:start -->
| $\sigma^{(2)}$ (cm⁴ s) at | 7.2 eV | 7.5 eV | 8 eV | 8.5 eV | 9.5 eV |
|---|---|---|---|---|---|
| L'Huillier, Tang & Lambropoulos 1989, MQDT (raster) | 3.11e-50 | — | — | — | — |
| L'Huillier, Tang & Lambropoulos 1989, above ²P₁/₂ (raster) | — | 1.14e-50 | 3.13e-51 | 1.83e-51 | — |
| Fink & Johnson 1990, linear (raster) | — | — | 1.9e-50 | 1.68e-50 | 1.24e-50 |
| Fink & Johnson 1990, circular (raster) | — | 1.56e-50 | 1.58e-50 | 1.49e-50 | 1.21e-50 |
| McGuire 1981, circular (raster) | — | — | 4.25e-50 | 4.86e-50 | 5.63e-50 |
<!-- two-photon-comparison-Kr:end -->

**Krypton has no vector-drawn or tabulated two-photon ionisation cross section and no
measurement.** Its only curves are read from scans:
- L'Huillier, Tang and Lambropoulos's MQDT, between the thresholds and above ²P₁/₂;
- Fink and Johnson's relativistic calculation, linear and circular, from the ²P₁/₂
  threshold (TC-912);
- McGuire's, circular polarisation.

Above ²P₁/₂ the MQDT and relativistic calculations differ by up to a factor of 6 at
8 eV. McGuire's circular values are 2.7 to 4.7 times Fink and Johnson's from 8 to 9.5
eV.

The data do not support a model of krypton.

## Xenon

<!-- two-photon-coverage-Xe:start -->
| Dataset | Form | Use | Polarisation | Channel | Photon energies in the window (eV) |
|---|---|---|---|---|---|
| `bischel_1979/alpha_Xe` | table | check | not stated | total | 6.424 |
| `fink_1990/fig5_circular` | raster | check | circular | total | 6.803 to 8.378 |
| `fink_1990/fig5_linear` | raster | check | linear | total | 6.803 to 8.381 |
| `fink_1990/fig6_circular_p12` | raster | check | circular | Xe+ 5p5 2P1/2 core | 6.803 to 7.894 |
| `fink_1990/fig6_circular_p32` | raster | check | circular | Xe+ 5p5 2P3/2 core | 6.966 to 8.003 |
| `fink_1990/fig6_linear_p12` | raster | check | linear | Xe+ 5p5 2P1/2 core | 6.966 to 7.891 |
| `fink_1990/fig6_linear_p32` | raster | check | linear | Xe+ 5p5 2P3/2 core | 6.963 to 8.000 |
| `gangopadhyay_1986/fig1` | raster | check | linear | Xe+ 2P3/2 (the only open core between the thresholds) | 6.103 to 6.598 |
| `gangopadhyay_1986/tab1_jc12` | table | check | linear | Xe+ 2P1/2 (residual core J = 1/2) | 7.000 to 8.435 |
| `gangopadhyay_1986/tab1_jc32` | table | check | linear | Xe+ 2P3/2 (residual core J = 3/2) | 7.000 to 8.435 |
| `hodges_1981/alpha_Xe` | text | check | not stated | total | 6.424 |
| `karamatskou_2017/fig1_interchannel` | vector | check | linear | total | outside |
| `karamatskou_2017/fig1_intrachannel` | vector | check | linear | total | outside |
| `karamatskou_2017/fig2_interchannel` | vector | check | linear | total | outside |
| `karamatskou_2017/fig3_interchannel_shifted` | vector | check | linear | total | outside |
| `kligler_1978/sigma_over_I_Xe` | table | check | not stated | total | 6.424 |
| `kulander_1988/text_sigma2_Xe` | text | check | linear | total | 6.424 |
| `lagutin_2017/fig3a_hf` | vector | check | linear | total | 7.037 to 8.417 |
| `lagutin_2017/fig3a_pf` | vector | check | linear | total | 7.020 to 8.405 |
| `lagutin_2017/fig3a_pfcp` | vector | check | linear | total | 7.075 to 8.435 |
| `lagutin_2017/fig3b_hf` | vector | check | linear | total | outside |
| `lagutin_2017/fig3b_pf` | vector | check | linear | total | outside |
| `lagutin_2017/fig3b_pfcp` | vector | check | linear | total | outside |
| `lagutin_2017/fig4_cipfcp_length` | vector | check | linear | total | 6.722 to 8.436 |
| `lagutin_2017/fig4_cipfcp_velocity` | vector | check | linear | total | 6.721 to 8.436 |
| `lagutin_2017/fig4_pf_length` | vector | check | linear | total | 6.722 to 8.426 |
| `lagutin_2017/fig4_pf_velocity` | vector | check | linear | total | 6.722 to 8.393 |
| `lagutin_2017/fig6_present_w` | vector | check | linear | total | 6.722 to 8.436 |
| `lagutin_2017/fig6_present_wo` | vector | check | linear | total | 6.722 to 8.436 |
| `lagutin_2017/fig7_cipf_length` | vector | check | linear | 5p3/2^5 ep3/2 (J = 0) | 6.803 to 8.435 |
| `lagutin_2017/fig7_cipf_velocity` | vector | check | linear | 5p3/2^5 ep3/2 (J = 0) | 6.823 to 8.435 |
| `lagutin_2017/fig7_noshakeup_length` | vector | check | linear | 5p3/2^5 ep3/2 (J = 0) | 6.815 to 8.437 |
| `lagutin_2017/fig7_noshakeup_velocity` | vector | check | linear | 5p3/2^5 ep3/2 (J = 0) | 6.804 to 8.428 |
| `lagutin_2017/fig8_cipfcp_circular` | vector | check | circular | total | 6.721 to 8.436 |
| `lagutin_2017/fig8_cipfcp_linear` | vector | check | linear | total | 6.721 to 8.436 |
| `lhuillier_1987/fig6_hf_Xe` | raster | check | linear | 5p shell (5p -> eps p, eps f) | 6.286 to 8.434 |
| `lhuillier_1987/fig6_ldrpa_Xe` | raster | check | linear | 5p shell (5p -> eps p, eps f) | 6.313 to 8.436 |
| `lhuillier_1987/fig6_rpae_Xe` | raster | check | linear | 5p shell (5p -> eps p, eps f) | 6.286 to 8.436 |
| `lhuillier_1987/fig7_hf_Xe` | raster | check | linear | 5p shell (5p -> eps p, eps f) | outside |
| `lhuillier_1987/fig7_ldrpa_Xe` | raster | check | linear | 5p shell (5p -> eps p, eps f) | outside |
| `lhuillier_1987/fig7_rpae_Xe` | raster | check | linear | 5p shell (5p -> eps p, eps f) | outside |
| `lhuillier_1987b/fig2_lda` | raster | check | linear | 5p shell | 6.346 to 8.002 |
| `lhuillier_1987b/fig2_ldrpa_first` | raster | check | linear | 5p shell | 6.089 to 8.348 |
| `lhuillier_1987b/fig2_ldrpa` | raster | check | linear | 5p shell | 6.102 to 8.414 |
| `lhuillier_1987b/fig3_rpae_first` | raster | check | linear | 5p shell | 6.316 to 8.436 |
| `lhuillier_1987b/fig4_ldrpa` | raster | check | linear | 5p shell | 6.074 to 8.411 |
| `lhuillier_1987b/fig4_rpae` | raster | check | linear | 5p shell | 6.289 to 8.411 |
| `lhuillier_1989/fig10_total_Xe` | raster | check | linear | total | 6.850 to 8.028 |
| `lhuillier_1989/fig6_J0_Xe` | raster | check | linear | J = 0 final states (partial) | 6.386 to 6.558 |
| `lhuillier_1989/fig6_total_Xe` | raster | check | linear | total (J = 0 and J = 2 final states) | 6.386 to 6.570 |
| `mccann_1988/sigma2_Xe_118_355` | text | check | linear | total | outside |
| `mccown_1982/alpha_Xe` | table | primary | natural | total | 6.424 |
| `mccown_1982/bischel_Xe` | table | check | not stated | total | 6.424 |
| `mccown_1982/hodges_Xe` | table | check | not stated | total | 6.424 |
| `mccown_1982/mcguire_Xe` | table | check | circular | total | 6.424 |
| `mcguire_1981/fig4_sigma2_Xe` | raster | check | circular | total | 6.390 to 8.211 |
| `santra_2004/sigma2_Xe` | text | check | linear | total | 6.420 |
<!-- two-photon-coverage-Xe:end -->

<!-- two-photon-comparison-Xe:start -->
| $\sigma^{(2)}$ (cm⁴ s) at | 6.42 eV | 7 eV | 7.5 eV | 8 eV |
|---|---|---|---|---|
| Lagutin et al. 2017, CIPFCP, linear | — | 4.27e-50 | 3.7e-50 | 3.44e-50 |
| Lagutin et al. 2017, CIPFCP, circular | — | 4.46e-50 (1.05) | 3.99e-50 (1.08) | 3.15e-50 (0.92) |
| Karamatskou & Santra 2017, TDCIS (Fig. 1) | — | — | — | — |
| Fink & Johnson 1990, linear (raster) | — | 3.59e-50 (0.84) | 2.77e-50 (0.75) | 3.41e-50 (0.99) |
| Gangopadhyay et al. 1986, MQDT, Table I (sum of cores) | — | 1.98e-50 (0.46) | 1.38e-50 (0.37) | 2.36e-50 (0.69) |
| Gangopadhyay et al. 1986, MQDT, between the thresholds (raster) | 4.77e-51 | — | — | — |
| L'Huillier & Wendin 1987 (J. Phys. B), RPAE (raster) | 2.61e-50 | 2.02e-50 (0.47) | 1.35e-50 (0.36) | 9.35e-51 (0.27) |
| L'Huillier & Wendin 1987, RPAE (raster) | 2.59e-50 | 2.05e-50 (0.48) | 1.46e-50 (0.39) | 9.39e-51 (0.27) |
| L'Huillier, Tang & Lambropoulos 1989, MQDT (raster) | 5.71e-50 | — | — | — |
| L'Huillier, Tang & Lambropoulos 1989, above ²P₁/₂ (raster) | — | 2.51e-50 (0.59) | 1.71e-50 (0.46) | 2.09e-50 (0.61) |
| McCown, Ediger & Eden 1982, measured | 4.12e-50 | — | — | — |
| Bischel et al. 1979, estimate | 1.03e-52 | — | — | — |
| Hodges, Lee & Moseley 1981, lower limit | 5.15e-52 | — | — | — |
| Santra & Greene 2004 | 6.2e-50 | — | — | — |
| McGuire 1981, circular (raster) | 8.35e-50 | 1.19e-49 (2.80) | 1.29e-49 (3.50) | 1.37e-49 (3.99) |
<!-- two-photon-comparison-Xe:end -->

Above the upper threshold the calculations disagree by factors of 2 to 4:
- **Lagutin et al.'s** relativistic, correlated calculation (TC-880) begins at the
  ²P₁/₂ threshold (6.72 eV) and has nothing between the thresholds.
- **L'Huillier and Wendin's RPAE and L'Huillier, Tang and Lambropoulos's MQDT** (raster)
  are 0.3 to 0.6 times it.
- **Fink and Johnson's** relativistic calculation (raster; TC-912) is 0.75 to 0.99 of
  it, and Lagutin et al. redraw it.
- **Gangopadhyay et al.'s MQDT Table I** (TC-914) is 0.37 to 0.69 of it. Every later
  drawing of "Gangopadhyay et al." is this table, above the ²P₁/₂ threshold.
- **Karamatskou and Santra's TDCIS** (TC-881) starts above the window, at 9 eV. There
  it is 20 to 200 times Lagutin et al.'s off resonance.
- **The only measurement**, McCown et al.'s at 193 nm (4.1×10⁻⁵⁰ cm⁴ s; TC-884), lies
  between the thresholds, where only the raster MQDT curve and Santra and Greene's
  single-threshold calculation reach. It is drawn by L'Huillier, Tang and Lambropoulos
  1.39 times higher than its printed value, and by Lagutin et al. and L'Huillier and
  Wendin at it (TC-882, TC-883). Gangopadhyay et al.'s raster curve between the
  thresholds (TC-914) reads 0.11 of it at the line centre. Bischel et al.'s estimate and
  Hodges et al.'s lower limit at 193 nm, both with McCown's definition, are 400 and 80
  times below it (TC-917, TC-918).

**Circular polarisation:** in the window, Lagutin et al.'s circular cross section is
0.92 to 1.08 times their linear one; McGuire's circular values are 2.8 to 4 times
Lagutin's linear ones.


See {doc}`/oeps/archive/oep-0006-two-photon-ionisation` for the completed proposal and its phase reports.
