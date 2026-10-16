# Quantity vocabulary

Use this reference when interpreting a dataset or adding a source record. Read
{doc}`/data` for the record structure and {doc}`/conventions` for the package’s
physical conventions. The tables below describe stored quantities, rather than the
parameters of the model interfaces in the {doc}`/api/index`.

## Quantities

Every quantity belongs to a **family**, which is also the directory of its files:
`refractive` (the real linear response and the dispersion it produces), `absorption`
(cross sections and the oscillator-strength distribution), `third_order`, `raman` (the
anisotropic and vibrational response, and the lines), `ionisation`, `damage` (laser
damage and breakdown of solids), `composition` (the mole fractions of a gas mixture
such as dry air; TC-543), `transport` (the dilute gas's viscosity and the ideal gas's
heat capacity, OEP 1) and `thermodynamic` (the virial coefficients and densities of the
gases and their mixtures, OEP 4). A quantity of a pair of gases is filed under the first
of the pair in sorted order, with the other as `parameters.partner`; a mole fraction is
always the material's. A metal's optical constants are `refractive_index` and
`extinction_coefficient`, or `permittivity_real` and `permittivity_imaginary` where the
source tabulates the permittivity. The vocabulary is {mod}`optmat.vocabulary`; the table
below is checked against it (`tests/test_datasets.py`). "Needs `state`" applies to a
gas; a solid has no state.

A solid's datasets also state, in `parameters`, whenever the source does (TC-319):

- `temperature_k`, the temperature of the sample;
- `form`: `bulk`, `film` or `fibre` (a film's deposition and thickness go in `notes`);
- `grade`: the product or melt as the source names it (`"Corning 7940"`,
  `"Suprasil 1"`, `"Lithosil Q2"`);
- `axis`, for an anisotropic crystal: one of the material's axes (see *Materials*);
- `polarization`, `technique` and `pulse_fs` for a nonlinear or damage measurement;
- `index_medium` for a refractive index (and its temperature coefficient): `vacuum`
  for an absolute index, `air` for one relative to the room air of the measurement
  (with `air_temperature_k` and `air_pressure_pa` where stated), or `not stated`
  (TC-398);
- `assumed = true` for a value a source takes as a model input (a band gap, a
  collision time, an effective mass in a simulation) rather than measures or derives.
  Such a value is `use = "check"`, with `reproduced_from` where the source cites one.

`form`, `site` (`bulk` or `surface`), `index_medium` and an absorption coefficient's
`base` (`"e"` or `"10"`) take only the values listed.

A gas's datasets may also state, in free form:

- `initial_state`: the level a transition leaves, as the source labels it (for
  example `"X 1Sigma_g+ (v = 0, N = 1), ortho-H2"` for an ionisation energy);
- `ortho_para`: the ortho:para ratio of a hydrogen or deuterium sample (`"3:1"`);
- `temperature_range_k`: a temperature range, where the source gives no single value
  (`"296.0-298.5"`). A model that needs one number takes it from the choices file
  (`data/gas/raman.toml`), which gives its reason.

<!-- quantities:start -->
| Quantity | Family | Meaning | `unit` | Abscissae | Needs `state` | Required `parameters` | Extra columns |
|---|---|---|---|---|---|---|---|
| `refractivity` | refractive | n - 1 at the dataset's state | `1` | `wavelength` | yes |  |  |
| `polarizability` | refractive | per-particle dipole polarizability; parameters.temperature_k for a temperature-dependent (rotationally averaged) value | `SI`, `angstrom3`, `au`, `cm3` | `wavelength` | no |  |  |
| `molar_refractivity` | refractive | A_R = N_A alpha / (3 eps0); parameters.temperature_k where stated | `cm3/mol` | `wavelength` | no |  |  |
| `refractivity_per_density` | refractive | (n - 1) / N at the dataset's state | `cm3` | `wavelength` | yes |  |  |
| `refractivity_ratio` | refractive | (n - 1) of the gas over (n - 1) of parameters.reference_material at the same density | `1` | `wavelength` | no | `reference_material` |  |
| `refractivity_difference` | refractive | n(lambda) - n(parameters.reference_wavelength_nm) at the dataset's state | `1` | `wavelength` | yes | `reference_wavelength_nm` |  |
| `partial_refractivity` | refractive | one stated contribution to n - 1 at the dataset's state (parameters.contribution: what it is and how it was obtained) | `1` | `wavelength` | yes | `contribution` |  |
| `polarizability_correction` | refractive | a correction to the electronic polarizability (parameters.part, one of POLARIZABILITY_CORRECTIONS), never a total | `SI`, `angstrom3`, `au`, `cm3` | `wavelength` | no | `part` |  |
| `magnetic_susceptibility` | refractive | the diamagnetic (or paramagnetic) susceptibility: per particle (au, the magnetizability chi_0) or molar (cm3/mol, m3/mol: A_mu or chi_M) | `au`, `cm3/mol`, `m3/mol` | `none`, `wavelength` | no |  |  |
| `cauchy_moment` | refractive | dipole-oscillator-strength sum S(k), the order k as the abscissa | `au` | `order` | no |  |  |
| `phase_match_density_ratio` | refractive | density ratio rho_gas / rho_reference at which two gases phase-match a harmonic (parameters.reference_material, .harmonic) | `1` | `wavelength` | no | `harmonic`, `reference_material` |  |
| `phase_match_density` | refractive | density at which a harmonic is phase-matched against a stated mismatch (parameters.harmonic, .phase_mismatch_per_m) | as printed | `wavelength` | no | `harmonic`, `phase_mismatch_per_m` |  |
| `phase_mismatch_per_density` | refractive | wave-vector mismatch of a harmonic per particle, Delta k / N (parameters.harmonic) | as printed | `wavelength` | no |  |  |
| `group_delay` | refractive | group delay per unit length at the dataset's state | `fs/m` | `wavelength` | yes |  |  |
| `group_delay_dispersion` | refractive | GDD per unit length at the dataset's state | `fs2/m` | `wavelength` | yes |  |  |
| `third_order_dispersion` | refractive | TOD per unit length at the dataset's state | `fs3/m` | `wavelength` | yes |  |  |
| `dispersion_kappa` | refractive | kappa = -lambda^2 d/dlambda[(n - 1)/lambda] at the dataset's state | `1` | `wavelength` | yes |  |  |
| `verdet_constant` | refractive | Verdet constant at the dataset's state (unit as printed) | as printed | `wavelength` | yes |  |  |
| `verdet_constant_ratio` | refractive | Verdet constant over its value at parameters.reference_wavelength_nm (the medium of the reference wavelength as the dataset's) | `1` | `wavelength` | no | `reference_wavelength_nm` |  |
| `rayleigh_cross_section` | refractive | Rayleigh scattering cross section per particle (parameters.king_factor as used) | `cm2` | `wavelength` | no |  |  |
| `rayleigh_cross_section_ratio` | refractive | Rayleigh cross section of the gas over that of parameters.reference_material | `1` | `wavelength` | no | `reference_material` |  |
| `refractivity_virial_coefficient` | refractive | refractivity (or dielectric) virial coefficient of parameters.order (2 for B_R, the default; 3 for C_R; 4 for D_R); parameters.temperature_k | `cm12/mol4`, `cm6/mol2`, `cm9/mol3` | `wavelength` | no |  |  |
| `unlike_pair_refractivity_virial_coefficient` | refractive | the refractivity (or dielectric) virial coefficient B_st of the dataset's material with parameters.partner, the first in sorted order (the cell inf is static); parameters.temperature_k, .contribution (induced, or total including orientation) | `cm6/mol2` | `wavelength` | no | `contribution`, `partner`, `temperature_k` |  |
| `refractive_index` | refractive | n, the real part of a solid's complex index n + ik (parameters.axis for an anisotropic crystal; .temperature_k, .grade, .form where stated) | `1` | `wavelength` | no |  |  |
| `permittivity_real` | refractive | epsilon_1, the real part of a solid's relative permittivity epsilon_1 + i epsilon_2 = (n + ik)^2 (parameters as for refractive_index) | `1` | `wavelength` | no |  |  |
| `refractive_index_decrement` | refractive | delta = 1 - n, the decrement of the real part of a solid's complex index below 1, as x-ray tables print it (parameters as for refractive_index) | `1` | `wavelength` | no |  |  |
| `refractive_index_deviation` | refractive | a solid's measured n minus the value of the formula parameters.reference (a formula id, source/formula) at the same wavelength | `1` | `wavelength` | no | `reference` |  |
| `thermo_optic_coefficient` | refractive | dn/dT of a solid (parameters.temperature_k) | `1/K` | `wavelength` | no |  |  |
| `absorption_cross_section` | absorption | total photoabsorption cross section per particle, below or above the ionisation threshold (parameters.temperature_k) | `Mb`, `cm2` | `wavelength` | no |  |  |
| `photoionization_cross_section` | absorption | photoionisation cross section per particle, above the ionisation threshold | `Mb`, `cm2` | `wavelength` | no |  |  |
| `excited_state_photoionization_cross_section` | absorption | one-photon photoionisation cross section of an excited level of the atom or molecule (parameters.initial_level: its configuration, term and J, or the source's label), against the photon's wavelength | `Mb`, `cm2` | `wavelength` | no | `initial_level` |  |
| `photodissociation_quantum_yield` | absorption | number of a stated product (or of a channel's events) formed per photon absorbed, against the photon's wavelength (parameters.product; .temperature_k where stated) | `1` | `wavelength` | no | `product` |  |
| `transition_dipole_moment` | absorption | electronic transition dipole moment between two states of the dataset's species (parameters.states: the pair, the row label; .geometry: the nuclear geometry, e.g. a bond length, where the source gives more than one; .component: parallel or perpendicular to the axis, where stated) | `D`, `au` | `label`, `none` | no | `states` |  |
| `two_photon_transition_strength` | absorption | rotationally averaged two-photon transition strength delta^TP of the bound transition from the ground state to the labelled final level, in atomic units (parameters.polarization: linear_parallel or circular; .frequencies: degenerate, or the second photon's energy; .rotational_average: the rotational average's F, G and H as the source states them) | `au` | `label` | no | `frequencies`, `polarization`, `rotational_average` |  |
| `two_photon_transition_amplitude` | absorption | signed second-order amplitude of the bound two-photon transition from the ground state to a stated final level, both photons polarised along z and both orderings summed, against the first photon's energy (the second's is the level energy less it), as the source defines it (parameters.final_level; .gauge; .polarization; .definition: the source's equation; .m_zz_per_unit: the zz element of the two-photon kernel M_zz, in atomic units, per unit of the printed amplitude; .method). Bound final levels only | `1`, `au` | `label`, `wavelength` | no | `definition`, `final_level`, `gauge`, `m_zz_per_unit`, `method`, `polarization` |  |
| `two_photon_excitation_cross_section` | absorption | frequency-integrated generalised cross section sigma_0^(2) of the two-photon excitation of a bound level from the ground state, W = sigma_0^(2) g Phi^2 with Phi the photon flux and g the line profile normalised to one (parameters.frequency_variable: angular or ordinary, the variable g is normalised in; .polarization; .photon_statistics: corrected, or the G2(0) of the measurement); a label names the final level, a wavelength is the exciting photon's | `cm4`, `m4` | `label`, `wavelength` | no | `frequency_variable`, `photon_statistics`, `polarization` |  |
| `integrated_cross_section` | absorption | absorption cross section integrated from the abscissa to the x_upper column (unit as printed: cm2 nm, cm2 cm-1) | as printed | `wavelength` | no |  | `x_upper` |
| `polarizability_imaginary` | absorption | imaginary part of the polarizability above threshold | `SI`, `angstrom3`, `au`, `cm3` | `wavelength` | no |  |  |
| `oscillator_strength` | absorption | absorption oscillator strength of the transition at the abscissa, or of the labelled transition | `1` | `label`, `wavelength` | no |  |  |
| `transition_probability` | absorption | Einstein coefficient A_ki of spontaneous emission from the upper level k to the lower level i of the transition at the abscissa, or of the labelled transition (parameters.notation: the level notation of the labels, e.g. paschen or racah) | `1/s` | `label`, `wavelength` | no | `notation` |  |
| `radiative_lifetime` | absorption | radiative lifetime of the labelled level of the atom or molecule, the inverse of its total rate of spontaneous decay in the absence of collisions (parameters.method: how it was obtained, e.g. measured and extrapolated to zero pressure, semiempirical or calculated) | `fs`, `ns`, `ps`, `s` | `label` | no | `method` |  |
| `level_polarizability` | absorption | dipole polarizability of a stated level of the atom (parameters.configuration, .term, .j), at the abscissa's wavelength or static (inf): the level shifts by -(1/2) alpha F^2 in a static field F and by -(1/4) alpha F^2 in an oscillating field of amplitude F (parameters.component: scalar, vector or tensor, with .definition the source's shift formula for a vector or tensor part; .complex_part: real or imaginary; .method). Not the gas's polarizability, which is polarizability | `SI`, `angstrom3`, `au`, `cm3`, `kHz/(kV/cm)2` | `wavelength` | no | `complex_part`, `component`, `configuration`, `j`, `method`, `term` |  |
| `light_shift_coefficient` | absorption | shift of the frequency of the multiphoton transition from the ground level to a stated final level per unit laser intensity, in the linear regime, positive to the blue, against the exciting photon's wavelength (parameters.final_level; .photons; .intensity_definition: how the source defines the intensity; .polarization where stated) | `Hz/(W/cm2)` | `wavelength` | no | `final_level`, `intensity_definition`, `photons` |  |
| `pressure_line_coefficient` | absorption | collisional (pressure) broadening or shift of a bound line of the atom, per unit perturber pressure or number density, in the impact regime (parameters.effect: broadening, or shift, positive to the blue; .transition: the line in the source's notation; .photons: 1 or 2; .frequency_scale: transition, the line's own frequency, laser, the exciting photon's, or not stated; .width: fwhm, hwhm, not stated, or none for a shift; .frequency_variable: ordinary or angular; .perturber: self or a gas; .temperature_k where stated). A two-photon line's broadening is the decay of the ground-to-final coherence: its pure dephasing is what remains after the radiative and quenching parts | `GHz/amagat`, `MHz/mbar`, `MHz/torr`, `cm-1/amagat`, `cm3/s` | `label`, `temperature` | no | `effect`, `frequency_scale`, `frequency_variable`, `perturber`, `photons`, `transition`, `width` |  |
| `collisional_quenching_rate_coefficient` | absorption | rate coefficient of the collisional removal of the labelled excited level of the atom or molecule by a perturber, from 1/tau = 1/tau_rad + k2 [M] + k3 [M]^2 (parameters.perturber: a gas, or self; .order: 2, bimolecular, or 3, termolecular; .channel: total, or the product level of a state-to-state coefficient; .temperature_k where stated). A loss of population, not the pure dephasing of the level | `cm3/s`, `cm6/s`, `m3/s`, `m6/s` | `label` | no | `channel`, `order`, `perturber` |  |
| `interval_oscillator_strength` | absorption | oscillator strength integrated from the abscissa to the x_upper column | `1` | `wavelength` | no |  | `x_upper` |
| `partial_cauchy_moment` | absorption | S(parameters.order) of the oscillator strength from the abscissa to the x_upper column | `au` | `wavelength` | no | `order` | `x_upper` |
| `differential_oscillator_strength` | absorption | df/dE of the continuum at the photon energy | `1/eV`, `1/hartree`, `1/rydberg` | `wavelength` | no |  |  |
| `extinction_coefficient` | absorption | k, the imaginary part of a solid's complex index n + ik (k > 0 absorbs; parameters as for refractive_index) | `1` | `wavelength` | no |  |  |
| `permittivity_imaginary` | absorption | epsilon_2, the imaginary part of a solid's relative permittivity (epsilon_2 > 0 absorbs; parameters as for refractive_index) | `1` | `wavelength` | no |  |  |
| `internal_transmittance` | absorption | the fraction of the light entering a sample that reaches its far surface, reflection losses excluded, as a glass catalogue prints it (parameters.thickness_mm: the path length) | `1` | `wavelength` | no | `thickness_mm` |  |
| `absorption_coefficient` | absorption | attenuation coefficient of the intensity per length in a solid (parameters.base: e for I = I0 exp(-a z), 10 for I = I0 10^(-a z)) | `cm-1`, `m-1` | `wavelength` | no | `base` |  |
| `band_gap` | absorption | band-gap energy of a solid (parameters.kind: optical, tauc, urbach, photoconductivity, exciton, mobility, effective (fitted to a model), assumed for a value a model takes as input, or not stated) | `eV` | `label`, `none` | no | `kind` |  |
| `hyperpolarizability` | third_order | second hyperpolarizability gamma per molecule (parameters.process, .convention, .part) | `SI`, `au`, `esu` | `wavelength` | no | `convention`, `part`, `process` |  |
| `hyperpolarizability_ratio` | third_order | gamma of the gas over gamma of parameters.reference_material (parameters.process, .part) | `1` | `wavelength` | no | `part`, `process`, `reference_material` |  |
| `nonlinear_index` | third_order | n2, n = n0 + n2 I, at the dataset's state for a gas (parameters.part, .coefficient; for a solid also .pulse_fs, .technique, .polarization where stated) | `cm2/W`, `m2/W` | `wavelength` | yes | `coefficient`, `part` |  |
| `nonlinear_index_field` | third_order | n2 per squared field, n = n0 + n2 <E^2> (the esu n2 of the solid-state literature; parameters.definition: the source's field average, e.g. '<E^2>' of the real field, or the squared amplitude; .part, .coefficient) | `esu`, `m2/V2` | `wavelength` | no | `coefficient`, `definition`, `part` |  |
| `nonlinear_index_ratio` | third_order | n2 of the gas over n2 of parameters.reference_material at the same density (parameters.part) | `1` | `wavelength` | no | `part`, `reference_material` |  |
| `nonlinear_susceptibility` | third_order | third-order susceptibility chi(3), per volume at the dataset's state or per unit density as printed (parameters.process, .convention, .part) | as printed | `wavelength` | no | `convention`, `part`, `process` |  |
| `kerr_constant` | third_order | dc Kerr constant at the dataset's state, as printed (parameters.temperature_k) | as printed | `wavelength` | yes |  |  |
| `molar_kerr_constant` | third_order | the zero-density molar Kerr constant A_K (parameters.temperature_k; against temperature, parameters.wavelength_nm) | as printed | `wavelength`, `temperature` | no |  |  |
| `kerr_virial_coefficient` | third_order | Kerr virial coefficient, B_K (parameters.order 2, the default) or C_K (order 3) (parameters.temperature_k; against temperature, parameters.wavelength_nm) | as printed | `wavelength`, `temperature` | no |  |  |
| `hyperpolarizability_density_coefficient` | third_order | the density (pair) coefficient of chi(3) per molecule, chi(3)/N = gamma-term [1 + b N + ...] (parameters.process, .part; parameters.pair for a mixture) | as printed | `wavelength` | no | `part`, `process` |  |
| `two_photon_absorption_coefficient` | third_order | two-photon absorption coefficient beta of a solid, dI/dz = -alpha I - beta I^2 (parameters.pulse_fs, .technique where stated) | `cm/GW`, `cm/MW`, `cm/W`, `m/W` | `wavelength` | no |  |  |
| `multiphoton_absorption_coefficient` | third_order | N-photon absorption coefficient of a solid, dI/dz = -beta_N I^N (parameters.photons; unit as printed) | as printed | `wavelength` | no | `photons` |  |
| `critical_power` | third_order | critical power for self-focusing (parameters.definition: the formula the source uses, or measured) | `MW`, `W` | `wavelength` | no | `definition` |  |
| `polarizability_anisotropy` | raman | polarizability anisotropy Delta alpha per molecule | `SI`, `angstrom3`, `au`, `cm3` | `wavelength` | no |  |  |
| `polarizability_derivative` | raman | derivative of the polarizability with respect to a normal coordinate (parameters.component: mean or anisotropy; .mode; .coordinate where it is not the mass-weighted normal coordinate: a bond length or a symmetry coordinate, as the source defines it; unit as printed) | as printed | `wavelength` | no | `component`, `mode` |  |
| `depolarization_ratio` | raman | depolarization ratio of scattered light (parameters.kind: rayleigh, raman, rotational_raman; .polarization: linear or natural) | `1` | `wavelength` | no | `kind`, `polarization` |  |
| `raman_shift` | raman | Raman shift or band origin of the labelled transition (parameters.transition when not per row) | `cm-1` | `label`, `none` | no |  |  |
| `rotational_constant` | raman | spectroscopic constant of the rotational ladder (parameters.constant: B, D, H, alpha_e, ...; .vibrational_state_v) | `GHz`, `MHz`, `cm-1` | `label`, `none` | no | `constant` |  |
| `vibrational_term_value` | raman | vibrational energy G(v) of the level labelled by the row (the vibrational quantum number v), from parameters.reference: v0 (the ground level) or minimum (the potential minimum); .electronic_state | `cm-1` | `label`, `none` | no | `reference` |  |
| `raman_scattering_activity` | raman | Raman scattering activity of a normal mode, the invariant combination g (45 a'^2 + 7 g'^2) of the derivatives of the mean (a') and anisotropic (g') polarizability with respect to the mass-weighted normal coordinate, g the mode's degeneracy (parameters.mode; unit as printed, e.g. angstrom4/amu) | as printed | `label`, `none`, `wavelength` | no | `mode` |  |
| `raman_cross_section_ratio` | raman | Raman scattering cross section of a transition relative to that of parameters.reference_material's reference transition (.reference_transition, usually the N2 Q branch) at the same excitation (parameters.transition; .geometry: how the scattering was observed) | `1` | `label`, `none`, `wavelength` | no | `reference_material`, `transition` |  |
| `raman_cross_section` | raman | differential Raman scattering cross section of a transition (parameters.transition; unit as printed, e.g. cm2/sr) | as printed | `wavelength` | no | `transition` |  |
| `raman_linewidth` | raman | linewidth of a Raman transition, or its pressure-broadening coefficient (parameters.transition, .width: hwhm, fwhm or not stated; .perturber for a coefficient: self, a gas id or air; .term for one term of a width law, e.g. the Dicke coefficient A of A/rho; .temperature_k where stated; unit as printed) | as printed | `density`, `label`, `none`, `pressure`, `temperature` | no | `transition`, `width` |  |
| `raman_line_shift` | raman | pressure (collisional) shift of a Raman line, or its coefficient (parameters.transition, .perturber, .temperature_k where stated; unit as printed) | as printed | `density`, `label`, `none`, `pressure`, `temperature` | no | `perturber`, `transition` |  |
| `broadening_temperature_exponent` | raman | exponent N of a broadening coefficient's temperature law gamma(T) = gamma(T0) (T0 / T)^N (parameters.transition, .perturber, .reference_temperature_k) | `1` | `label`, `none` | no | `perturber`, `transition` |  |
| `raman_transition_polarizability` | raman | matrix element of a polarizability invariant between the two states of a Raman transition, alpha_if (parameters.transition, .invariant: mean or anisotropy, or the component parallel or perpendicular to the molecular axis) | `SI`, `angstrom3`, `au`, `cm3` | `label`, `none`, `wavelength` | no | `invariant`, `transition` |  |
| `raman_gain` | raman | steady-state stimulated Raman gain coefficient of a transition at the dataset's state (parameters.transition, .pump_nm when not the abscissa; unit as printed, e.g. cm/GW) | as printed | `density`, `label`, `none`, `pressure`, `wavelength` | yes | `transition` |  |
| `dephasing_time` | raman | coherence (T2) time of a Raman transition (parameters.transition; unit as printed) | as printed | `density`, `label`, `none`, `pressure`, `temperature` | no | `transition` |  |
| `raman_fraction` | raman | fraction f_R of the delayed (Raman) response in a solid's nonlinear response (parameters.definition: how the source defines it) | `1` | `label`, `none` | no | `definition` |  |
| `raman_mode` | raman | one mode of a multi-mode Raman response of a solid: its amplitude, with the columns shift (cm-1) and the gaussian_width and lorentzian_width (as printed, parameters.width: fwhm or hwhm) | as printed | `label`, `none` | no | `width` | `gaussian_width`, `lorentzian_width`, `shift` |
| `raman_response_time` | raman | a time constant of a model Raman response (parameters.constant: tau1, tau2, ...; .model) | `fs`, `ns`, `ps`, `s` | `label`, `none` | no | `constant`, `model` |  |
| `ionization_energy` | ionisation | ionisation energy (parameters.kind: adiabatic or vertical; .final_state) | `cm-1`, `eV`, `hartree` | `label`, `none` | no | `kind` |  |
| `level_energy` | ionisation | energy of a level of an atom or ion above its ground level (parameters.configuration, .term, .j) | `cm-1`, `eV` | `label`, `none` | no | `configuration`, `j`, `term` |  |
| `electronic_state_energy` | ionisation | energy of an excited electronic state of a molecule above its ground state (parameters.state; .kind: Te, the potential minimum, T0, the lowest vibrational level, or onset, the long-wavelength onset of a band; .one_photon and .two_photon from the ground state: allowed, forbidden or spin_forbidden, where stated) | `cm-1`, `eV` | `label`, `none` | no | `kind`, `state` |  |
| `appearance_energy` | ionisation | appearance energy of a fragment ion (parameters.fragment) | `cm-1`, `eV` | `label`, `none` | no | `fragment` |  |
| `dissociation_energy` | ionisation | energy of a dissociation limit of the dataset's species (neutral or ion) into the stated products (parameters.products, with their states where stated; .kind: D0 from the lowest level, De from the potential minimum, or threshold: the limit measured from the ground state of parameters.reference_state, e.g. the neutral molecule for a dication's limit) | `cm-1`, `eV`, `hartree`, `kJ/mol` | `label`, `none` | no | `kind`, `products` |  |
| `ionic_state_population` | ionisation | population of an electronic state of the molecular ion that a pulse leaves, as a fraction of the ions or of the molecules (parameters.state; .normalisation: ions or molecules; .wavelength_nm, .pulse_fs, .polarization: the pulse; .intensity_w_cm2 where one number; .method: measured, or the calculation), against peak intensity or as one number | `1` | `intensity`, `label`, `none` | no | `normalisation`, `polarization`, `pulse_fs`, `state`, `wavelength_nm` |  |
| `fragmentation_branching_ratio` | ionisation | yield of one fragmentation channel of the dataset's species over that of a reference (parameters.channel: the products; .reference: the parent ion, all ions, all channels, or another channel; .excitation: photoionisation at the abscissa's photon energy, or a strong-field pulse with .wavelength_nm, .pulse_fs, .polarization; .detection_correction: how the detector's efficiency was treated), against intensity, the photon's wavelength or energy, or as one number | `1` | `intensity`, `label`, `none`, `wavelength` | no | `channel`, `detection_correction`, `excitation`, `reference` |  |
| `kinetic_energy_release` | ionisation | total kinetic energy the fragments of a channel carry away (parameters.channel; .kind: peak, mean, asymptotic or calculated; .state: the parent state where assigned; the pulse where measured) | `eV` | `label`, `none` | no | `channel`, `kind` |  |
| `dissociation_time` | ionisation | time constant of a dissociation of the dataset's species (parameters.channel; .kind: lifetime of the dissociating state, appearance time of a fragment, or calculated; a collisional time states its .pressure and partner) | `fs`, `ns`, `ps`, `s` | `label`, `none`, `pressure` | no | `channel`, `kind` |  |
| `dissociative_recombination_rate_coefficient` | ionisation | rate coefficient of the dissociative recombination of the dataset's molecular ion with an electron, k in dn_e/dt = -k n_e n_ion (parameters.products, or all channels; .law: value, the coefficient at the abscissa's electron temperature, or power_law, the value at .reference_temperature_k scaled by (T_ref/T_e)^.exponent; .temperature_range_k where stated) | `cm3/s`, `m3/s` | `label`, `none`, `temperature` | no | `law`, `products` |  |
| `ionization_rate` | ionisation | field-ionisation rate against intensity or field (parameters.wavelength_nm, inf for a static field; .polarization: linear, circular or static; unit as printed) | as printed | `field`, `intensity` | no | `polarization`, `wavelength_nm` |  |
| `ionization_yield` | ionisation | ion yield or probability against intensity (parameters.wavelength_nm, .pulse_fs, .polarization; .pulse_envelope and .cep where stated; unit as printed) | as printed | `intensity` | no | `polarization`, `pulse_fs`, `wavelength_nm` |  |
| `ground_state_survival` | ionisation | probability that a pulse leaves the atom in its initial (ground) state, against the peak intensity, read at a stated time; its complement counts excitation as well as ionisation (parameters.wavelength_nm, .pulse_fs, .polarization, .read_at: when the population is read; .pulse_envelope and .cep where stated) | `1` | `intensity` | no | `polarization`, `pulse_fs`, `read_at`, `wavelength_nm` |  |
| `ionization_yield_spectrum` | ionisation | ion yield or probability against the photon energy or wavelength at a fixed peak intensity (parameters.intensity_w_cm2, .pulse_cycles: the pulse's length in optical cycles, which fixes its duration at each wavelength, .polarization; .pulse_envelope and .cep where stated; unit as printed) | as printed | `wavelength` | no | `intensity_w_cm2`, `polarization`, `pulse_cycles` |  |
| `ionization_yield_ratio` | ionisation | ratio of the material's ionisation yield, as the source defines it (singly charged ions, or all channels, or electrons), to a companion species' measured with it, against intensity, wavelength or as one number (parameters.companion; .wavelength_nm, .pulse_fs, .polarization where one value; .detection_correction: how the ion detector's efficiency was treated, e.g. none, Z/sqrt(M), optical) | `1` | `intensity`, `label`, `wavelength` | no | `companion`, `detection_correction`, `polarization` |  |
| `ionization_phase_shift` | ionisation | phase shift of a probe beam from the electrons a pump pulse leaves, averaged over the probe's profile and normalised to a gas pressure, against the pump's peak intensity (parameters.wavelength_nm, .pulse_fs, .polarization: the pump; .probe_wavelength_nm; .pressure_atm: the pressure normalised to; the beams' sizes where stated) | `rad` | `intensity` | no | `polarization`, `pressure_atm`, `probe_wavelength_nm`, `pulse_fs`, `wavelength_nm` |  |
| `tunnelling_structure_coefficient` | ionisation | asymptotic coefficient of the ionising orbital in a tunnelling theory, one per partial wave l, the row label (parameters.model: moadk, adk or wfat; .orbital; .m; .normalisation: how the coefficient is defined) | `1` | `label`, `none` | no | `m`, `model`, `normalisation`, `orbital` |  |
| `ionization_model_parameter` | ionisation | a constant fitted within a named ionisation-rate model, not a property of the material alone (parameters.model, .parameter: its symbol; unit as printed) | as printed | `label`, `none` | no | `model`, `parameter` |  |
| `saturation_intensity` | ionisation | peak intensity at which a pulse saturates the ionisation of an atom, as the source defines it (parameters.definition; .photons, .pulse_fs, .pulse_shape: the pulse) | `W/cm2`, `W/m2` | `wavelength` | no | `definition`, `pulse_fs` |  |
| `multiphoton_ionization_cross_section` | ionisation | generalised N-photon ionisation cross section (parameters.photons; unit as printed, cm^2N s^(N-1)) | as printed | `wavelength` | no | `photons` |  |
| `multiphoton_ionization_polarization_ratio` | ionisation | ratio of the generalised N-photon ionisation cross section for circular polarisation to that for linear (parameters.photons) | `1` | `wavelength` | no | `photons` |  |
| `effective_mass` | ionisation | effective mass of a solid's carriers, in electron masses (parameters.carrier: electron, hole, reduced or not stated) | `m_e` | `label`, `none` | no | `carrier` |  |
| `electron_collision_time` | ionisation | momentum-relaxation (collision) time of a solid's conduction electrons in a Drude model | `fs`, `ns`, `ps`, `s` | `label`, `none` | no |  |  |
| `electron_trapping_time` | ionisation | decay (trapping or recombination) time of the conduction-electron density | `fs`, `ns`, `ps`, `s` | `label`, `none` | no |  |  |
| `avalanche_coefficient` | ionisation | impact-ionisation (avalanche) coefficient of a solid (unit as printed: cm2/J, 1/s, ...; parameters.definition) | as printed | `label`, `none`, `wavelength` | no | `definition` |  |
| `atom_density` | ionisation | density of the particles available for ionisation in a solid (parameters.counts: molecules, atoms or valence_electrons) | `cm-3`, `m-3` | `label`, `none` | no | `counts` |  |
| `damage_threshold_fluence` | damage | laser-induced damage or breakdown threshold fluence (parameters.site: bulk or surface; .criterion; .shots; the pulse length as the abscissa (duration, with .wavelength_nm) or as .pulse_fs) | `J/cm2`, `J/m2` | `duration`, `wavelength` | no | `site` |  |
| `damage_threshold_intensity` | damage | laser-induced damage or breakdown threshold intensity (parameters as for damage_threshold_fluence) | `W/cm2`, `W/m2` | `duration`, `wavelength` | no | `site` |  |
| `damage_threshold_energy` | damage | pulse energy at the damage or breakdown threshold for a stated focusing (parameters as for damage_threshold_fluence, and .numerical_aperture) | `J`, `nJ`, `uJ` | `duration`, `wavelength` | no | `site` |  |
| `mole_fraction` | composition | amount-of-substance fraction of the material (a gas) in a mixture (parameters.mixture: the mixture, as the source defines it, such as 'dry air') | `1`, `mol/mol`, `ppm` | `none` | no | `mixture` |  |
| `dilute_viscosity` | transport | the shear viscosity of the gas in the limit of zero density, a function of temperature alone | `Pa s` | `temperature` | no |  |  |
| `ideal_gas_heat_capacity` | transport | the molar isobaric heat capacity of the ideal gas, c_p0(T) | `J/(mol K)` | `temperature` | no |  |  |
| `second_virial_coefficient` | thermodynamic | the second (pressure, density) virial coefficient B(T) of the pure gas, Z = 1 + B rho + ... | `cm3/mol` | `temperature` | no |  |  |
| `cross_second_virial_coefficient` | thermodynamic | the unlike-pair second virial coefficient B_12(T) of the dataset's material with parameters.partner, B_mix = sum_st x_s x_t B_st | `cm3/mol` | `temperature` | no | `partner` |  |
| `excess_second_virial_coefficient` | thermodynamic | the excess second virial coefficient E = B_12 - (B_11 + B_22)/2 of the dataset's material with parameters.partner | `cm3/mol` | `temperature` | no | `partner` |  |
| `mixture_second_virial_coefficient` | thermodynamic | the second virial coefficient B_mix(T) of a binary mixture of the dataset's material, at mole fraction parameters.mole_fraction, with parameters.partner | `cm3/mol` | `temperature` | no | `mole_fraction`, `partner` |  |
| `mixture_compressibility_factor` | thermodynamic | Z = p / (rho R T) of a binary mixture of the dataset's material, at mole fraction parameters.mole_fraction, with parameters.partner; parameters.temperature_k | `1` | `density`, `pressure` | no | `mole_fraction`, `partner`, `temperature_k` |  |
| `mixture_density` | thermodynamic | the molar (or mass) density of a binary mixture of the dataset's material, at mole fraction parameters.mole_fraction, with parameters.partner; parameters.temperature_k | `g/cm3`, `kg/m3`, `mol/L`, `mol/cm3`, `mol/m3` | `pressure` | no | `mole_fraction`, `partner`, `temperature_k` |  |
<!-- quantities:end -->

The first four refractive quantities each determine the polarizability at their
frequency ({func}`optmat.gas.reduction.polarizability`). The others constrain the dispersion indirectly;
the fits give them forward models. A `polarizability_correction` (a zero-point or pure
vibrational correction) is never read as a total.

The abscissae and their units:

| Abscissa | Units | Cells |
|---|---|---|
| `wavelength` | `m`, `nm`, `um`, `angstrom`, `cm-1`, `um-1`, `eV`, `hartree`, `rydberg` | numbers; `inf` for the static limit; `not stated` |
| `order` | `1` | the order of a moment |
| `temperature` | `K`, `C` | numbers |
| `pressure` | `Pa`, `kPa`, `MPa`, `bar`, `atm`, `torr`, `psia` (lbf/in^2) | numbers |
| `density` | `m-3`, `cm-3`, `amagat`, `mol/m3`, `mol/cm3`, `mol/L` | numbers |
| `intensity` | `W/m2`, `W/cm2` | numbers |
| `duration` | `s`, `ns`, `ps`, `fs` | numbers (a pulse length) |
| `field` | `V/m`, `V/cm`, `au` | numbers (a static or peak field; `au` is E_h / (e a0)) |
| `label` | none | names of the rows (a transition, a state) |
| `none` | none | empty: the value stands alone, described by its parameters |

`Dataset.vacuum_wavelength()` converts a wavelength abscissa, `Dataset.abscissa_values()`
any numeric one to SI, and `Dataset.column_values(name)` an extra column.

## Third-order quantities

The third-order values keep the source's own quantity, convention and conditions; nothing
is converted on entry. Each states:
- **`process`:**
  - `static`;
  - `kerr`, the intensity-dependent index gamma(-w; w, w, -w);
  - `eshg`, gamma(-2w; w, w, 0);
  - `dc_kerr`, gamma(-w; w, 0, 0);
  - `thg`, gamma(-3w; w, w, w);
  - `cars`, gamma(-w_as; w_p, w_p, -w_s);
  - `other`.
  The abscissa is the fundamental (input) wavelength, or the pump for `cars`,
  with the Stokes wavelength in `parameters` or `notes`.
- **`convention`:** one of the named definitions of `optmat.vocabulary.CONVENTIONS`:
  - `T`, the Taylor-series (Buckingham-Orr) convention of Shelton & Rice (p = ... + gamma
    E^3/6);
  - `ward_1969` (Ward's per-molecule chi in esu);
  - `lehmeier_1985` (P = 4 eps0 chi E E E);
  - `chi_cars_per_amagat`;
  - `rosasco_1985_hellwarth`;
  - `hellwarth_sigma`;
  - `maker_terhune_c`;
  - `bjorklund_1975`.

  A new definition is added there, with its relation to gamma_T where known.
- **`part`:** `electronic`, `vibrational`, `rotational`, `vibrational_rotational`
  (Shelton's gamma^vR), `electronic_vibrational`, `total`, `instantaneous` (a
  pulse-dependent mixture, as Nibbering's), or `unknown`.
- **`coefficient`** (for n2): `self` (n = n0 + n2 I of one beam) or `probe` (the index a
  weak probe sees, which is twice the self coefficient for the same polarization).
- **Measurement conditions:** the pulse duration and the probe wavelength go in
  `parameters` (`pulse_fs`, `probe_nm`) where the value depends on them.

A published dispersion fit, such as gamma = gamma0 [1 + A nu_L^2 + B nu_L^4] with nu_L^2
the sum of the squares of the four wavenumbers, is a `[[formulas]]` entry without a
`form`. Its coefficients are strings as printed, and its variable and fitted process are
stated.
