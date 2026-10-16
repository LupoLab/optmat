# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The vocabulary of the curated data: quantities, abscissae, parts, conventions.

Every dataset states one *quantity* from :data:`QUANTITIES`. The quantity fixes:

- the *family* of the data, which is also the directory of its files
  (``data/values/<family>/<material>/``);
- the units its values may carry (after the dataset's scale);
- the *abscissae* its values may be tabulated against (:data:`ABSCISSAE`): a
  wavelength, the order of a moment, a temperature, a pressure, a density, an
  intensity, an electric field, a free label, or none;
- whether it needs the state (pressure, temperature) it refers to;
- the parameters it must state;
- the extra per-row columns its files must carry (the upper edge of an interval).

Third-order values also state their process (:data:`PROCESSES`), the contributions they
contain (:data:`PARTS`) and, for a hyperpolarizability or susceptibility, the convention
of its definition (:data:`CONVENTIONS`). The documentation page *The curated data*
lists everything here.
"""

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Final

__all__ = [
    "ABSCISSAE",
    "CONVENTIONS",
    "FAMILIES",
    "PARTS",
    "POLARIZABILITY_CORRECTIONS",
    "POLARIZABILITY_QUANTITIES",
    "PROCESSES",
    "QUANTITIES",
    "Quantity",
]

FAMILIES: Final = (
    "refractive",
    "absorption",
    "third_order",
    "raman",
    "ionisation",
    "damage",
    "composition",
    "transport",
    "thermodynamic",
)
"""The families of quantities, each a directory of ``data/values/``."""

_WAVELENGTH_UNITS: Final = frozenset(
    {"m", "nm", "um", "angstrom", "cm-1", "um-1", "eV", "hartree", "rydberg"}
)

ABSCISSAE: Final = MappingProxyType(
    {
        "wavelength": _WAVELENGTH_UNITS,
        "order": frozenset({"1"}),
        "temperature": frozenset({"K", "C"}),
        "pressure": frozenset({"Pa", "kPa", "MPa", "bar", "atm", "torr", "psia"}),
        "density": frozenset({"m-3", "cm-3", "amagat", "mol/m3", "mol/cm3", "mol/L"}),
        "intensity": frozenset({"W/m2", "W/cm2"}),
        "field": frozenset({"V/m", "V/cm", "au"}),
        "duration": frozenset({"s", "ns", "ps", "fs"}),
        "label": frozenset({""}),
        "none": frozenset({""}),
    }
)
"""The kinds of abscissa, each with its units.

``wavelength`` is a wavelength, a vacuum (or air) wavenumber or a photon energy; the
cell ``inf`` is the static limit and ``not stated`` a wavelength the source does not
give. ``order`` is the order of a moment. ``duration`` is a pulse duration (a
laser-damage threshold against pulse length). ``field`` is the amplitude of a static
(or peak) electric field; ``au`` is the atomic unit of field, E_h / (e a_0). A
``label`` abscissa is not a number: its cells name the rows (a transition, a state).
``none`` has empty cells: the value stands alone (an ionisation energy), described by
its parameters.
"""

_POLARIZABILITY_UNITS: Final = frozenset({"au", "angstrom3", "cm3", "SI"})


@dataclass(frozen=True)
class Quantity:
    """What a dataset's values are, and what a dataset of them must state.

    Attributes
    ----------
    meaning : str
        The definition, with its parameters.
    family : str
        One of :data:`FAMILIES`.
    units : frozenset of str or None
        The units the values may carry after ``scale``; None for any unit (stated as
        printed, the quantity's meaning says how to read it).
    abscissae : frozenset of str
        The abscissae the values may be tabulated against.
    needs_state : bool
        Whether the values refer to a state (they depend on density).
    required : tuple of str
        Parameters every dataset of the quantity must state.
    columns : tuple of str
        Extra per-row columns every data file of the quantity must carry.
    """

    meaning: str
    family: str
    units: frozenset[str] | None = None
    abscissae: frozenset[str] = field(default=frozenset({"wavelength"}))
    needs_state: bool = False
    required: tuple[str, ...] = ()
    columns: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Check the family and the abscissae.

        Raises
        ------
        ValueError
            If either is unknown.
        """
        if self.family not in FAMILIES:
            msg = f"unknown family {self.family!r}"
            raise ValueError(msg)
        if not self.abscissae <= set(ABSCISSAE):
            msg = f"unknown abscissae {sorted(self.abscissae - set(ABSCISSAE))}"
            raise ValueError(msg)


_ONE: Final = frozenset({"1"})
_LINES: Final = frozenset({"label", "none"})
_CONDITIONS: Final = frozenset({"label", "none", "temperature", "pressure", "density"})
_TIMES: Final = frozenset({"fs", "ps", "ns", "s"})
_SPECTRAL_OR_NONE: Final = frozenset({"wavelength", "label", "none"})

QUANTITIES: Final = MappingProxyType(
    {
        # The linear (real) response and the dispersion it produces.
        "refractivity": Quantity(
            "n - 1 at the dataset's state", "refractive", _ONE, needs_state=True
        ),
        "polarizability": Quantity(
            "per-particle dipole polarizability; parameters.temperature_k for a "
            "temperature-dependent (rotationally averaged) value",
            "refractive",
            _POLARIZABILITY_UNITS,
        ),
        "molar_refractivity": Quantity(
            "A_R = N_A alpha / (3 eps0); parameters.temperature_k where stated",
            "refractive",
            frozenset({"cm3/mol"}),
        ),
        "refractivity_per_density": Quantity(
            "(n - 1) / N at the dataset's state",
            "refractive",
            frozenset({"cm3"}),
            needs_state=True,
        ),
        "refractivity_ratio": Quantity(
            "(n - 1) of the gas over (n - 1) of parameters.reference_material at the "
            "same density",
            "refractive",
            _ONE,
            required=("reference_material",),
        ),
        "refractivity_difference": Quantity(
            "n(lambda) - n(parameters.reference_wavelength_nm) at the dataset's state",
            "refractive",
            _ONE,
            needs_state=True,
            required=("reference_wavelength_nm",),
        ),
        "partial_refractivity": Quantity(
            "one stated contribution to n - 1 at the dataset's state "
            "(parameters.contribution: what it is and how it was obtained)",
            "refractive",
            _ONE,
            needs_state=True,
            required=("contribution",),
        ),
        "polarizability_correction": Quantity(
            "a correction to the electronic polarizability (parameters.part, one of "
            "POLARIZABILITY_CORRECTIONS), never a total",
            "refractive",
            _POLARIZABILITY_UNITS,
            required=("part",),
        ),
        "magnetic_susceptibility": Quantity(
            "the diamagnetic (or paramagnetic) susceptibility: per particle (au, the "
            "magnetizability chi_0) or molar (cm3/mol, m3/mol: A_mu or chi_M)",
            "refractive",
            frozenset({"au", "cm3/mol", "m3/mol"}),
            abscissae=frozenset({"wavelength", "none"}),
        ),
        "cauchy_moment": Quantity(
            "dipole-oscillator-strength sum S(k), the order k as the abscissa",
            "refractive",
            frozenset({"au"}),
            abscissae=frozenset({"order"}),
        ),
        "phase_match_density_ratio": Quantity(
            "density ratio rho_gas / rho_reference at which two gases phase-match a "
            "harmonic (parameters.reference_material, .harmonic)",
            "refractive",
            _ONE,
            required=("reference_material", "harmonic"),
        ),
        "phase_match_density": Quantity(
            "density at which a harmonic is phase-matched against a stated mismatch "
            "(parameters.harmonic, .phase_mismatch_per_m)",
            "refractive",
            abscissae=frozenset({"wavelength"}),
            required=("harmonic", "phase_mismatch_per_m"),
        ),
        "phase_mismatch_per_density": Quantity(
            "wave-vector mismatch of a harmonic per particle, Delta k / N "
            "(parameters.harmonic)",
            "refractive",
        ),
        "group_delay": Quantity(
            "group delay per unit length at the dataset's state",
            "refractive",
            frozenset({"fs/m"}),
            needs_state=True,
        ),
        "group_delay_dispersion": Quantity(
            "GDD per unit length at the dataset's state",
            "refractive",
            frozenset({"fs2/m"}),
            needs_state=True,
        ),
        "third_order_dispersion": Quantity(
            "TOD per unit length at the dataset's state",
            "refractive",
            frozenset({"fs3/m"}),
            needs_state=True,
        ),
        "dispersion_kappa": Quantity(
            "kappa = -lambda^2 d/dlambda[(n - 1)/lambda] at the dataset's state",
            "refractive",
            _ONE,
            needs_state=True,
        ),
        "verdet_constant": Quantity(
            "Verdet constant at the dataset's state (unit as printed)",
            "refractive",
            needs_state=True,
        ),
        "verdet_constant_ratio": Quantity(
            "Verdet constant over its value at parameters.reference_wavelength_nm "
            "(the medium of the reference wavelength as the dataset's)",
            "refractive",
            _ONE,
            required=("reference_wavelength_nm",),
        ),
        "rayleigh_cross_section": Quantity(
            "Rayleigh scattering cross section per particle (parameters.king_factor "
            "as used)",
            "refractive",
            frozenset({"cm2"}),
        ),
        "rayleigh_cross_section_ratio": Quantity(
            "Rayleigh cross section of the gas over that of "
            "parameters.reference_material",
            "refractive",
            _ONE,
            required=("reference_material",),
        ),
        "refractivity_virial_coefficient": Quantity(
            "refractivity (or dielectric) virial coefficient of parameters.order "
            "(2 for B_R, the default; 3 for C_R; 4 for D_R); "
            "parameters.temperature_k",
            "refractive",
            frozenset({"cm6/mol2", "cm9/mol3", "cm12/mol4"}),
        ),
        # OEP 4, section 5: a pair is not a material; filed under the first species
        # of the pair in sorted order, the other as the partner.
        "unlike_pair_refractivity_virial_coefficient": Quantity(
            "the refractivity (or dielectric) virial coefficient B_st of the "
            "dataset's material with parameters.partner, the first in sorted order "
            "(the cell inf is static); parameters.temperature_k, .contribution "
            "(induced, or total including orientation)",
            "refractive",
            frozenset({"cm6/mol2"}),
            required=("partner", "temperature_k", "contribution"),
        ),
        # The complex index of a solid, n + ik (TC-319).
        "refractive_index": Quantity(
            "n, the real part of a solid's complex index n + ik (parameters.axis for "
            "an anisotropic crystal; .temperature_k, .grade, .form where stated)",
            "refractive",
            _ONE,
        ),
        "permittivity_real": Quantity(
            "epsilon_1, the real part of a solid's relative permittivity epsilon_1 + "
            "i epsilon_2 = (n + ik)^2 (parameters as for refractive_index)",
            "refractive",
            _ONE,
        ),
        "refractive_index_decrement": Quantity(
            "delta = 1 - n, the decrement of the real part of a solid's complex index "
            "below 1, as x-ray tables print it (parameters as for refractive_index)",
            "refractive",
            _ONE,
        ),
        "refractive_index_deviation": Quantity(
            "a solid's measured n minus the value of the formula parameters.reference "
            "(a formula id, source/formula) at the same wavelength",
            "refractive",
            _ONE,
            required=("reference",),
        ),
        "thermo_optic_coefficient": Quantity(
            "dn/dT of a solid (parameters.temperature_k)",
            "refractive",
            frozenset({"1/K"}),
        ),
        # Absorption and the oscillator-strength distribution.
        "absorption_cross_section": Quantity(
            "total photoabsorption cross section per particle, below or above the "
            "ionisation threshold (parameters.temperature_k)",
            "absorption",
            frozenset({"cm2", "Mb"}),
        ),
        "photoionization_cross_section": Quantity(
            "photoionisation cross section per particle, above the ionisation "
            "threshold",
            "absorption",
            frozenset({"cm2", "Mb"}),
        ),
        "excited_state_photoionization_cross_section": Quantity(
            "one-photon photoionisation cross section of an excited level of the atom "
            "or molecule "
            "(parameters.initial_level: its configuration, term and J, or the "
            "source's label), against the photon's wavelength",
            "absorption",
            frozenset({"cm2", "Mb"}),
            required=("initial_level",),
        ),
        "photodissociation_quantum_yield": Quantity(
            "number of a stated product (or of a channel's events) formed per photon "
            "absorbed, against the photon's wavelength (parameters.product; "
            ".temperature_k where stated)",
            "absorption",
            _ONE,
            required=("product",),
        ),
        "transition_dipole_moment": Quantity(
            "electronic transition dipole moment between two states of the dataset's "
            "species (parameters.states: the pair, the row label; .geometry: the "
            "nuclear geometry, e.g. a bond length, where the source gives more than "
            "one; .component: parallel or perpendicular to the axis, where stated)",
            "absorption",
            frozenset({"au", "D"}),
            abscissae=_LINES,
            required=("states",),
        ),
        "two_photon_transition_strength": Quantity(
            "rotationally averaged two-photon transition strength delta^TP of the "
            "bound transition from the ground state to the labelled final level, in "
            "atomic units (parameters.polarization: linear_parallel or circular; "
            ".frequencies: degenerate, or the second photon's energy; "
            ".rotational_average: the rotational average's F, G and H as the source "
            "states them)",
            "absorption",
            frozenset({"au"}),
            abscissae=frozenset({"label"}),
            required=("polarization", "frequencies", "rotational_average"),
        ),
        "two_photon_transition_amplitude": Quantity(
            "signed second-order amplitude of the bound two-photon transition from the "
            "ground state to a stated final level, both photons polarised along z and "
            "both orderings summed, against the first photon's energy (the second's is "
            "the level energy less it), as the source defines it (parameters."
            "final_level; .gauge; .polarization; .definition: the source's equation; "
            ".m_zz_per_unit: the zz element of the two-photon kernel M_zz, in atomic "
            "units, per unit of the printed amplitude; .method). Bound final levels "
            "only",
            "absorption",
            frozenset({"1", "au"}),
            abscissae=frozenset({"wavelength", "label"}),
            required=(
                "final_level",
                "gauge",
                "polarization",
                "definition",
                "m_zz_per_unit",
                "method",
            ),
        ),
        "two_photon_excitation_cross_section": Quantity(
            "frequency-integrated generalised cross section sigma_0^(2) of the "
            "two-photon excitation of a bound level from the ground state, W = "
            "sigma_0^(2) g Phi^2 with Phi the photon flux and g the line profile "
            "normalised to one (parameters.frequency_variable: angular or ordinary, "
            "the variable g is normalised in; .polarization; .photon_statistics: "
            "corrected, or the G2(0) of the measurement); a label names the final "
            "level, a wavelength is the exciting photon's",
            "absorption",
            frozenset({"cm4", "m4"}),
            abscissae=frozenset({"label", "wavelength"}),
            required=("frequency_variable", "polarization", "photon_statistics"),
        ),
        "integrated_cross_section": Quantity(
            "absorption cross section integrated from the abscissa to the x_upper "
            "column (unit as printed: cm2 nm, cm2 cm-1)",
            "absorption",
            columns=("x_upper",),
        ),
        "polarizability_imaginary": Quantity(
            "imaginary part of the polarizability above threshold",
            "absorption",
            _POLARIZABILITY_UNITS,
        ),
        "oscillator_strength": Quantity(
            "absorption oscillator strength of the transition at the abscissa, or of "
            "the labelled transition",
            "absorption",
            _ONE,
            abscissae=frozenset({"wavelength", "label"}),
        ),
        "transition_probability": Quantity(
            "Einstein coefficient A_ki of spontaneous emission from the upper level k "
            "to the lower level i of the transition at the abscissa, or of the "
            "labelled transition (parameters.notation: the level notation of the "
            "labels, e.g. paschen or racah)",
            "absorption",
            frozenset({"1/s"}),
            abscissae=frozenset({"wavelength", "label"}),
            required=("notation",),
        ),
        "radiative_lifetime": Quantity(
            "radiative lifetime of the labelled level of the atom or molecule, the "
            "inverse of its "
            "total rate of spontaneous decay in the absence of collisions "
            "(parameters.method: how it was obtained, e.g. measured and extrapolated "
            "to zero pressure, semiempirical or calculated)",
            "absorption",
            _TIMES,
            abscissae=frozenset({"label"}),
            required=("method",),
        ),
        "level_polarizability": Quantity(
            "dipole polarizability of a stated level of the atom (parameters."
            "configuration, .term, .j), at the abscissa's wavelength or static (inf): "
            "the level shifts by -(1/2) alpha F^2 in a static field F and by -(1/4) "
            "alpha F^2 in an oscillating field of amplitude F (parameters.component: "
            "scalar, vector or tensor, with .definition the source's shift formula for "
            "a vector or tensor part; .complex_part: real or imaginary; .method). Not "
            "the "
            "gas's polarizability, which is polarizability",
            "absorption",
            _POLARIZABILITY_UNITS | {"kHz/(kV/cm)2"},
            required=(
                "configuration",
                "term",
                "j",
                "component",
                "complex_part",
                "method",
            ),
        ),
        "light_shift_coefficient": Quantity(
            "shift of the frequency of the multiphoton transition from the ground "
            "level to a stated final level per unit laser intensity, in the linear "
            "regime, positive to the blue, against the exciting photon's wavelength "
            "(parameters.final_level; .photons; .intensity_definition: how the source "
            "defines the intensity; .polarization where stated)",
            "absorption",
            frozenset({"Hz/(W/cm2)"}),
            required=("final_level", "photons", "intensity_definition"),
        ),
        "pressure_line_coefficient": Quantity(
            "collisional (pressure) broadening or shift of a bound line of the atom, "
            "per unit perturber pressure or number density, in the impact regime "
            "(parameters.effect: broadening, or shift, positive to the blue; "
            ".transition: the line in the source's notation; .photons: 1 or 2; "
            ".frequency_scale: transition, the line's own frequency, laser, the "
            "exciting photon's, or not stated; .width: fwhm, hwhm, not stated, or none "
            "for a shift; .frequency_variable: ordinary or angular; .perturber: self "
            "or a gas; .temperature_k where stated). A two-photon line's broadening is "
            "the decay of the ground-to-final coherence: its pure dephasing is what "
            "remains after the radiative and quenching parts",
            "absorption",
            frozenset({"MHz/torr", "MHz/mbar", "GHz/amagat", "cm-1/amagat", "cm3/s"}),
            abscissae=frozenset({"label", "temperature"}),
            required=(
                "effect",
                "transition",
                "photons",
                "frequency_scale",
                "width",
                "frequency_variable",
                "perturber",
            ),
        ),
        "collisional_quenching_rate_coefficient": Quantity(
            "rate coefficient of the collisional removal of the labelled excited level "
            "of the atom or molecule by a perturber, from 1/tau = 1/tau_rad + k2 [M] "
            "+ k3 [M]^2 "
            "(parameters.perturber: a gas, or self; .order: 2, bimolecular, or 3, "
            "termolecular; .channel: total, or the product level of a state-to-state "
            "coefficient; .temperature_k where stated). A loss of population, not the "
            "pure dephasing of the level",
            "absorption",
            frozenset({"cm3/s", "m3/s", "cm6/s", "m6/s"}),
            abscissae=frozenset({"label"}),
            required=("perturber", "order", "channel"),
        ),
        "interval_oscillator_strength": Quantity(
            "oscillator strength integrated from the abscissa to the x_upper column",
            "absorption",
            _ONE,
            columns=("x_upper",),
        ),
        "partial_cauchy_moment": Quantity(
            "S(parameters.order) of the oscillator strength from the abscissa to the "
            "x_upper column",
            "absorption",
            frozenset({"au"}),
            required=("order",),
            columns=("x_upper",),
        ),
        "differential_oscillator_strength": Quantity(
            "df/dE of the continuum at the photon energy",
            "absorption",
            frozenset({"1/eV", "1/hartree", "1/rydberg"}),
        ),
        "extinction_coefficient": Quantity(
            "k, the imaginary part of a solid's complex index n + ik (k > 0 absorbs; "
            "parameters as for refractive_index)",
            "absorption",
            _ONE,
        ),
        "permittivity_imaginary": Quantity(
            "epsilon_2, the imaginary part of a solid's relative permittivity "
            "(epsilon_2 > 0 absorbs; parameters as for refractive_index)",
            "absorption",
            _ONE,
        ),
        "internal_transmittance": Quantity(
            "the fraction of the light entering a sample that reaches its far "
            "surface, reflection losses excluded, as a glass catalogue prints it "
            "(parameters.thickness_mm: the path length)",
            "absorption",
            _ONE,
            required=("thickness_mm",),
        ),
        "absorption_coefficient": Quantity(
            "attenuation coefficient of the intensity per length in a solid "
            "(parameters.base: e for I = I0 exp(-a z), 10 for I = I0 10^(-a z))",
            "absorption",
            frozenset({"cm-1", "m-1"}),
            required=("base",),
        ),
        "band_gap": Quantity(
            "band-gap energy of a solid (parameters.kind: optical, tauc, urbach, "
            "photoconductivity, exciton, mobility, effective (fitted to a model), "
            "assumed for a value a model takes as input, or not stated)",
            "absorption",
            frozenset({"eV"}),
            abscissae=_LINES,
            required=("kind",),
        ),
        # The third-order response.
        "hyperpolarizability": Quantity(
            "second hyperpolarizability gamma per molecule (parameters.process, "
            ".convention, .part)",
            "third_order",
            frozenset({"au", "esu", "SI"}),
            required=("process", "convention", "part"),
        ),
        "hyperpolarizability_ratio": Quantity(
            "gamma of the gas over gamma of parameters.reference_material "
            "(parameters.process, .part)",
            "third_order",
            _ONE,
            required=("reference_material", "process", "part"),
        ),
        "nonlinear_index": Quantity(
            "n2, n = n0 + n2 I, at the dataset's state for a gas (parameters.part, "
            ".coefficient; for a solid also .pulse_fs, .technique, .polarization "
            "where stated)",
            "third_order",
            frozenset({"cm2/W", "m2/W"}),
            needs_state=True,
            required=("part", "coefficient"),
        ),
        "nonlinear_index_field": Quantity(
            "n2 per squared field, n = n0 + n2 <E^2> (the esu n2 of the solid-state "
            "literature; parameters.definition: the source's field average, e.g. "
            "'<E^2>' of the real field, or the squared amplitude; .part, "
            ".coefficient)",
            "third_order",
            frozenset({"esu", "m2/V2"}),
            required=("definition", "part", "coefficient"),
        ),
        "nonlinear_index_ratio": Quantity(
            "n2 of the gas over n2 of parameters.reference_material at the same "
            "density (parameters.part)",
            "third_order",
            _ONE,
            required=("reference_material", "part"),
        ),
        "nonlinear_susceptibility": Quantity(
            "third-order susceptibility chi(3), per volume at the dataset's state or "
            "per unit density as printed (parameters.process, .convention, .part)",
            "third_order",
            required=("process", "convention", "part"),
        ),
        "kerr_constant": Quantity(
            "dc Kerr constant at the dataset's state, as printed "
            "(parameters.temperature_k)",
            "third_order",
            needs_state=True,
        ),
        "molar_kerr_constant": Quantity(
            "the zero-density molar Kerr constant A_K (parameters.temperature_k; "
            "against temperature, parameters.wavelength_nm)",
            "third_order",
            abscissae=frozenset({"wavelength", "temperature"}),
        ),
        "kerr_virial_coefficient": Quantity(
            "Kerr virial coefficient, B_K (parameters.order 2, the default) or C_K "
            "(order 3) (parameters.temperature_k; against temperature, "
            "parameters.wavelength_nm)",
            "third_order",
            abscissae=frozenset({"wavelength", "temperature"}),
        ),
        "hyperpolarizability_density_coefficient": Quantity(
            "the density (pair) coefficient of chi(3) per molecule, "
            "chi(3)/N = gamma-term [1 + b N + ...] (parameters.process, .part; "
            "parameters.pair for a mixture)",
            "third_order",
            required=("process", "part"),
        ),
        "two_photon_absorption_coefficient": Quantity(
            "two-photon absorption coefficient beta of a solid, dI/dz = -alpha I - "
            "beta I^2 (parameters.pulse_fs, .technique where stated)",
            "third_order",
            frozenset({"m/W", "cm/W", "cm/MW", "cm/GW"}),
        ),
        "multiphoton_absorption_coefficient": Quantity(
            "N-photon absorption coefficient of a solid, dI/dz = -beta_N I^N "
            "(parameters.photons; unit as printed)",
            "third_order",
            required=("photons",),
        ),
        "critical_power": Quantity(
            "critical power for self-focusing (parameters.definition: the formula "
            "the source uses, or measured)",
            "third_order",
            frozenset({"W", "MW"}),
            required=("definition",),
        ),
        # Raman: the anisotropic and vibrational response.
        "polarizability_anisotropy": Quantity(
            "polarizability anisotropy Delta alpha per molecule",
            "raman",
            _POLARIZABILITY_UNITS,
        ),
        "polarizability_derivative": Quantity(
            "derivative of the polarizability with respect to a normal coordinate "
            "(parameters.component: mean or anisotropy; .mode; .coordinate where it is "
            "not the mass-weighted normal coordinate: a bond length or a symmetry "
            "coordinate, as the source defines it; unit as printed)",
            "raman",
            required=("component", "mode"),
        ),
        "depolarization_ratio": Quantity(
            "depolarization ratio of scattered light (parameters.kind: rayleigh, "
            "raman, rotational_raman; .polarization: linear or natural)",
            "raman",
            _ONE,
            required=("kind", "polarization"),
        ),
        "raman_shift": Quantity(
            "Raman shift or band origin of the labelled transition "
            "(parameters.transition when not per row)",
            "raman",
            frozenset({"cm-1"}),
            abscissae=_LINES,
        ),
        "rotational_constant": Quantity(
            "spectroscopic constant of the rotational ladder (parameters.constant: B, "
            "D, H, alpha_e, ...; .vibrational_state_v)",
            "raman",
            frozenset({"cm-1", "MHz", "GHz"}),
            abscissae=_LINES,
            required=("constant",),
        ),
        "vibrational_term_value": Quantity(
            "vibrational energy G(v) of the level labelled by the row (the vibrational "
            "quantum number v), from parameters.reference: v0 (the ground level) or "
            "minimum (the potential minimum); .electronic_state",
            "raman",
            frozenset({"cm-1"}),
            abscissae=_LINES,
            required=("reference",),
        ),
        "raman_scattering_activity": Quantity(
            "Raman scattering activity of a normal mode, the invariant combination "
            "g (45 a'^2 + 7 g'^2) of the derivatives of the mean (a') and anisotropic "
            "(g') polarizability with respect to the mass-weighted normal coordinate, "
            "g the mode's degeneracy (parameters.mode; unit as printed, e.g. "
            "angstrom4/amu)",
            "raman",
            abscissae=_SPECTRAL_OR_NONE,
            required=("mode",),
        ),
        "raman_cross_section_ratio": Quantity(
            "Raman scattering cross section of a transition relative to that of "
            "parameters.reference_material's reference transition "
            "(.reference_transition, usually the N2 Q branch) at the same excitation "
            "(parameters.transition; .geometry: how the scattering was observed)",
            "raman",
            _ONE,
            abscissae=_SPECTRAL_OR_NONE,
            required=("transition", "reference_material"),
        ),
        "raman_cross_section": Quantity(
            "differential Raman scattering cross section of a transition "
            "(parameters.transition; unit as printed, e.g. cm2/sr)",
            "raman",
            required=("transition",),
        ),
        "raman_linewidth": Quantity(
            "linewidth of a Raman transition, or its pressure-broadening coefficient "
            "(parameters.transition, .width: hwhm, fwhm or not stated; .perturber for "
            "a coefficient: self, a gas id or air; .term for one term of a width law, "
            "e.g. the Dicke coefficient A of A/rho; .temperature_k where stated; unit "
            "as printed)",
            "raman",
            abscissae=_CONDITIONS,
            required=("transition", "width"),
        ),
        "raman_line_shift": Quantity(
            "pressure (collisional) shift of a Raman line, or its coefficient "
            "(parameters.transition, .perturber, .temperature_k where stated; unit as "
            "printed)",
            "raman",
            abscissae=_CONDITIONS,
            required=("transition", "perturber"),
        ),
        "broadening_temperature_exponent": Quantity(
            "exponent N of a broadening coefficient's temperature law gamma(T) = "
            "gamma(T0) (T0 / T)^N (parameters.transition, .perturber, "
            ".reference_temperature_k)",
            "raman",
            _ONE,
            abscissae=_LINES,
            required=("transition", "perturber"),
        ),
        "raman_transition_polarizability": Quantity(
            "matrix element of a polarizability invariant between the two states of a "
            "Raman transition, alpha_if (parameters.transition, .invariant: mean or "
            "anisotropy, or the component parallel or perpendicular to the molecular "
            "axis)",
            "raman",
            _POLARIZABILITY_UNITS,
            abscissae=_SPECTRAL_OR_NONE,
            required=("transition", "invariant"),
        ),
        "raman_gain": Quantity(
            "steady-state stimulated Raman gain coefficient of a transition at the "
            "dataset's state (parameters.transition, .pump_nm when not the abscissa; "
            "unit as printed, e.g. cm/GW)",
            "raman",
            abscissae=frozenset({"wavelength", "pressure", "density", "label", "none"}),
            needs_state=True,
            required=("transition",),
        ),
        "dephasing_time": Quantity(
            "coherence (T2) time of a Raman transition (parameters.transition; unit "
            "as printed)",
            "raman",
            abscissae=_CONDITIONS,
            required=("transition",),
        ),
        "raman_fraction": Quantity(
            "fraction f_R of the delayed (Raman) response in a solid's nonlinear "
            "response (parameters.definition: how the source defines it)",
            "raman",
            _ONE,
            abscissae=_LINES,
            required=("definition",),
        ),
        "raman_mode": Quantity(
            "one mode of a multi-mode Raman response of a solid: its amplitude, with "
            "the columns shift (cm-1) and the gaussian_width and lorentzian_width "
            "(as printed, parameters.width: fwhm or hwhm)",
            "raman",
            abscissae=_LINES,
            required=("width",),
            columns=("shift", "gaussian_width", "lorentzian_width"),
        ),
        "raman_response_time": Quantity(
            "a time constant of a model Raman response (parameters.constant: tau1, "
            "tau2, ...; .model)",
            "raman",
            _TIMES,
            abscissae=_LINES,
            required=("constant", "model"),
        ),
        # Ionisation.
        "ionization_energy": Quantity(
            "ionisation energy (parameters.kind: adiabatic or vertical; .final_state)",
            "ionisation",
            frozenset({"eV", "cm-1", "hartree"}),
            abscissae=_LINES,
            required=("kind",),
        ),
        "level_energy": Quantity(
            "energy of a level of an atom or ion above its ground level "
            "(parameters.configuration, .term, .j)",
            "ionisation",
            frozenset({"cm-1", "eV"}),
            abscissae=_LINES,
            required=("configuration", "term", "j"),
        ),
        "electronic_state_energy": Quantity(
            "energy of an excited electronic state of a molecule above its ground "
            "state (parameters.state; .kind: Te, the potential minimum, T0, the "
            "lowest vibrational level, or onset, the long-wavelength onset of a band; "
            ".one_photon and .two_photon from the ground state: allowed, forbidden or "
            "spin_forbidden, where stated)",
            "ionisation",
            frozenset({"cm-1", "eV"}),
            abscissae=_LINES,
            required=("state", "kind"),
        ),
        "appearance_energy": Quantity(
            "appearance energy of a fragment ion (parameters.fragment)",
            "ionisation",
            frozenset({"eV", "cm-1"}),
            abscissae=_LINES,
            required=("fragment",),
        ),
        "dissociation_energy": Quantity(
            "energy of a dissociation limit of the dataset's species (neutral or ion) "
            "into the stated products (parameters.products, with their states where "
            "stated; .kind: D0 from the lowest level, De from the potential minimum, "
            "or threshold: the limit measured from the ground state of parameters."
            "reference_state, e.g. the neutral molecule for a dication's limit)",
            "ionisation",
            frozenset({"eV", "cm-1", "hartree", "kJ/mol"}),
            abscissae=_LINES,
            required=("products", "kind"),
        ),
        "ionic_state_population": Quantity(
            "population of an electronic state of the molecular ion that a pulse "
            "leaves, as a fraction of the ions or of the molecules (parameters.state; "
            ".normalisation: ions or molecules; .wavelength_nm, .pulse_fs, "
            ".polarization: the pulse; .intensity_w_cm2 where one number; .method: "
            "measured, or the calculation), against peak intensity or as one number",
            "ionisation",
            _ONE,
            abscissae=frozenset({"intensity", "label", "none"}),
            required=(
                "state",
                "normalisation",
                "wavelength_nm",
                "pulse_fs",
                "polarization",
            ),
        ),
        "fragmentation_branching_ratio": Quantity(
            "yield of one fragmentation channel of the dataset's species over that of "
            "a reference (parameters.channel: the products; .reference: the parent "
            "ion, all ions, all channels, or another channel; .excitation: "
            "photoionisation at the abscissa's photon energy, or a strong-field pulse "
            "with .wavelength_nm, .pulse_fs, .polarization; .detection_correction: how "
            "the detector's efficiency was treated), against intensity, the photon's "
            "wavelength or energy, or as one number",
            "ionisation",
            _ONE,
            abscissae=frozenset({"intensity", "wavelength", "label", "none"}),
            required=("channel", "reference", "excitation", "detection_correction"),
        ),
        "kinetic_energy_release": Quantity(
            "total kinetic energy the fragments of a channel carry away "
            "(parameters.channel; .kind: peak, mean, asymptotic or calculated; "
            ".state: the parent state where assigned; the pulse where measured)",
            "ionisation",
            frozenset({"eV"}),
            abscissae=_LINES,
            required=("channel", "kind"),
        ),
        "dissociation_time": Quantity(
            "time constant of a dissociation of the dataset's species (parameters."
            "channel; .kind: lifetime of the dissociating state, appearance time of a "
            "fragment, or calculated; a collisional time states its .pressure and "
            "partner)",
            "ionisation",
            _TIMES,
            abscissae=frozenset({"label", "none", "pressure"}),
            required=("channel", "kind"),
        ),
        "dissociative_recombination_rate_coefficient": Quantity(
            "rate coefficient of the dissociative recombination of the dataset's "
            "molecular ion with an electron, k in dn_e/dt = -k n_e n_ion "
            "(parameters.products, or all channels; .law: value, the coefficient at "
            "the abscissa's electron temperature, or power_law, the value at "
            ".reference_temperature_k scaled by (T_ref/T_e)^.exponent; .temperature_"
            "range_k where stated)",
            "ionisation",
            frozenset({"cm3/s", "m3/s"}),
            abscissae=frozenset({"temperature", "label", "none"}),
            required=("products", "law"),
        ),
        "ionization_rate": Quantity(
            "field-ionisation rate against intensity or field (parameters."
            "wavelength_nm, inf for a static field; .polarization: linear, circular "
            "or static; unit as printed)",
            "ionisation",
            abscissae=frozenset({"intensity", "field"}),
            required=("wavelength_nm", "polarization"),
        ),
        "ionization_yield": Quantity(
            "ion yield or probability against intensity (parameters.wavelength_nm, "
            ".pulse_fs, .polarization; .pulse_envelope and .cep where stated; unit as "
            "printed)",
            "ionisation",
            abscissae=frozenset({"intensity"}),
            required=("wavelength_nm", "pulse_fs", "polarization"),
        ),
        "ground_state_survival": Quantity(
            "probability that a pulse leaves the atom in its initial (ground) state, "
            "against the peak intensity, read at a stated time; its complement counts "
            "excitation as well as ionisation (parameters.wavelength_nm, .pulse_fs, "
            ".polarization, .read_at: when the population is read; .pulse_envelope "
            "and .cep where stated)",
            "ionisation",
            _ONE,
            abscissae=frozenset({"intensity"}),
            required=("wavelength_nm", "pulse_fs", "polarization", "read_at"),
        ),
        "ionization_yield_spectrum": Quantity(
            "ion yield or probability against the photon energy or wavelength at a "
            "fixed peak intensity (parameters.intensity_w_cm2, .pulse_cycles: the "
            "pulse's length in optical cycles, which fixes its duration at each "
            "wavelength, .polarization; .pulse_envelope and .cep where stated; unit as "
            "printed)",
            "ionisation",
            abscissae=frozenset({"wavelength"}),
            required=("intensity_w_cm2", "pulse_cycles", "polarization"),
        ),
        "ionization_yield_ratio": Quantity(
            "ratio of the material's ionisation yield, as the source defines it "
            "(singly charged ions, or all channels, or electrons), to a companion "
            "species' measured with it, against intensity, wavelength or as one "
            "number (parameters.companion; .wavelength_nm, .pulse_fs, .polarization "
            "where one value; .detection_correction: how the ion detector's "
            "efficiency was treated, e.g. none, Z/sqrt(M), optical)",
            "ionisation",
            _ONE,
            abscissae=frozenset({"intensity", "wavelength", "label"}),
            required=("companion", "polarization", "detection_correction"),
        ),
        "ionization_phase_shift": Quantity(
            "phase shift of a probe beam from the electrons a pump pulse leaves, "
            "averaged over the probe's profile and normalised to a gas pressure, "
            "against the pump's peak intensity (parameters.wavelength_nm, .pulse_fs, "
            ".polarization: the pump; .probe_wavelength_nm; .pressure_atm: the "
            "pressure normalised to; the beams' sizes where stated)",
            "ionisation",
            frozenset({"rad"}),
            abscissae=frozenset({"intensity"}),
            required=(
                "wavelength_nm",
                "pulse_fs",
                "polarization",
                "probe_wavelength_nm",
                "pressure_atm",
            ),
        ),
        "tunnelling_structure_coefficient": Quantity(
            "asymptotic coefficient of the ionising orbital in a tunnelling theory, "
            "one per partial wave l, the row label (parameters.model: moadk, adk or "
            "wfat; .orbital; .m; .normalisation: how the coefficient is defined)",
            "ionisation",
            _ONE,
            abscissae=_LINES,
            required=("model", "orbital", "m", "normalisation"),
        ),
        "ionization_model_parameter": Quantity(
            "a constant fitted within a named ionisation-rate model, not a property "
            "of the material alone (parameters.model, .parameter: its symbol; unit as "
            "printed)",
            "ionisation",
            abscissae=_LINES,
            required=("model", "parameter"),
        ),
        "saturation_intensity": Quantity(
            "peak intensity at which a pulse saturates the ionisation of an atom, as "
            "the source defines it (parameters.definition; .photons, .pulse_fs, "
            ".pulse_shape: the pulse)",
            "ionisation",
            frozenset({"W/cm2", "W/m2"}),
            abscissae=frozenset({"wavelength"}),
            required=("definition", "pulse_fs"),
        ),
        "multiphoton_ionization_cross_section": Quantity(
            "generalised N-photon ionisation cross section (parameters.photons; unit "
            "as printed, cm^2N s^(N-1))",
            "ionisation",
            required=("photons",),
        ),
        "multiphoton_ionization_polarization_ratio": Quantity(
            "ratio of the generalised N-photon ionisation cross section for circular "
            "polarisation to that for linear (parameters.photons)",
            "ionisation",
            _ONE,
            required=("photons",),
        ),
        "effective_mass": Quantity(
            "effective mass of a solid's carriers, in electron masses "
            "(parameters.carrier: electron, hole, reduced or not stated)",
            "ionisation",
            frozenset({"m_e"}),
            abscissae=_LINES,
            required=("carrier",),
        ),
        "electron_collision_time": Quantity(
            "momentum-relaxation (collision) time of a solid's conduction electrons "
            "in a Drude model",
            "ionisation",
            _TIMES,
            abscissae=_LINES,
        ),
        "electron_trapping_time": Quantity(
            "decay (trapping or recombination) time of the conduction-electron density",
            "ionisation",
            _TIMES,
            abscissae=_LINES,
        ),
        "avalanche_coefficient": Quantity(
            "impact-ionisation (avalanche) coefficient of a solid (unit as printed: "
            "cm2/J, 1/s, ...; parameters.definition)",
            "ionisation",
            abscissae=_SPECTRAL_OR_NONE,
            required=("definition",),
        ),
        "atom_density": Quantity(
            "density of the particles available for ionisation in a solid "
            "(parameters.counts: molecules, atoms or valence_electrons)",
            "ionisation",
            frozenset({"m-3", "cm-3"}),
            abscissae=_LINES,
            required=("counts",),
        ),
        # Laser damage and breakdown of solids.
        "damage_threshold_fluence": Quantity(
            "laser-induced damage or breakdown threshold fluence (parameters.site: "
            "bulk or surface; .criterion; .shots; the pulse length as the abscissa "
            "(duration, with .wavelength_nm) or as .pulse_fs)",
            "damage",
            frozenset({"J/cm2", "J/m2"}),
            abscissae=frozenset({"wavelength", "duration"}),
            required=("site",),
        ),
        "damage_threshold_intensity": Quantity(
            "laser-induced damage or breakdown threshold intensity (parameters as for "
            "damage_threshold_fluence)",
            "damage",
            frozenset({"W/cm2", "W/m2"}),
            abscissae=frozenset({"wavelength", "duration"}),
            required=("site",),
        ),
        "damage_threshold_energy": Quantity(
            "pulse energy at the damage or breakdown threshold for a stated focusing "
            "(parameters as for damage_threshold_fluence, and .numerical_aperture)",
            "damage",
            frozenset({"J", "uJ", "nJ"}),
            abscissae=frozenset({"wavelength", "duration"}),
            required=("site",),
        ),
        # The composition of a gas mixture (TC-543).
        "mole_fraction": Quantity(
            "amount-of-substance fraction of the material (a gas) in a mixture "
            "(parameters.mixture: the mixture, as the source defines it, such as "
            "'dry air')",
            "composition",
            frozenset({"1", "ppm", "mol/mol"}),
            abscissae=frozenset({"none"}),
            required=("mixture",),
        ),
        # The dilute gases' transport and ideal-gas properties (OEP 1, section 4.3).
        "dilute_viscosity": Quantity(
            "the shear viscosity of the gas in the limit of zero density, a function "
            "of temperature alone",
            "transport",
            frozenset({"Pa s"}),
            abscissae=frozenset({"temperature"}),
        ),
        "ideal_gas_heat_capacity": Quantity(
            "the molar isobaric heat capacity of the ideal gas, c_p0(T)",
            "transport",
            frozenset({"J/(mol K)"}),
            abscissae=frozenset({"temperature"}),
        ),
        # The density of gas mixtures (OEP 4, sections 4.2 and 5). A pair is filed
        # under the first species in sorted order, with the other as the partner;
        # mole_fraction is always the material's.
        "second_virial_coefficient": Quantity(
            "the second (pressure, density) virial coefficient B(T) of the pure gas, "
            "Z = 1 + B rho + ...",
            "thermodynamic",
            frozenset({"cm3/mol"}),
            abscissae=frozenset({"temperature"}),
        ),
        "cross_second_virial_coefficient": Quantity(
            "the unlike-pair second virial coefficient B_12(T) of the dataset's "
            "material with parameters.partner, B_mix = sum_st x_s x_t B_st",
            "thermodynamic",
            frozenset({"cm3/mol"}),
            abscissae=frozenset({"temperature"}),
            required=("partner",),
        ),
        "excess_second_virial_coefficient": Quantity(
            "the excess second virial coefficient E = B_12 - (B_11 + B_22)/2 of the "
            "dataset's material with parameters.partner",
            "thermodynamic",
            frozenset({"cm3/mol"}),
            abscissae=frozenset({"temperature"}),
            required=("partner",),
        ),
        "mixture_second_virial_coefficient": Quantity(
            "the second virial coefficient B_mix(T) of a binary mixture of the "
            "dataset's material, at mole fraction parameters.mole_fraction, with "
            "parameters.partner",
            "thermodynamic",
            frozenset({"cm3/mol"}),
            abscissae=frozenset({"temperature"}),
            required=("partner", "mole_fraction"),
        ),
        "mixture_compressibility_factor": Quantity(
            "Z = p / (rho R T) of a binary mixture of the dataset's material, at mole "
            "fraction parameters.mole_fraction, with parameters.partner; "
            "parameters.temperature_k",
            "thermodynamic",
            _ONE,
            abscissae=frozenset({"pressure", "density"}),
            required=("partner", "mole_fraction", "temperature_k"),
        ),
        "mixture_density": Quantity(
            "the molar (or mass) density of a binary mixture of the dataset's "
            "material, at mole fraction parameters.mole_fraction, with "
            "parameters.partner; parameters.temperature_k",
            "thermodynamic",
            frozenset({"mol/m3", "mol/cm3", "mol/L", "kg/m3", "g/cm3"}),
            abscissae=frozenset({"pressure"}),
            required=("partner", "mole_fraction", "temperature_k"),
        ),
    }
)
"""Every quantity a dataset may hold."""

POLARIZABILITY_QUANTITIES: Final = frozenset(
    {"refractivity", "polarizability", "molar_refractivity", "refractivity_per_density"}
)
"""Quantities that each determine the polarizability at their frequency."""

PROCESSES: Final = frozenset(
    {"static", "kerr", "eshg", "dc_kerr", "thg", "cars", "other"}
)
"""Third-order processes: static, the intensity-dependent index gamma(-w; w, w, -w)
("kerr"), gamma(-2w; w, w, 0) ("eshg"), gamma(-w; w, 0, 0) ("dc_kerr"), gamma(-3w; w,
w, w) ("thg"), coherent anti-Stokes Raman gamma(-w_as; w_p, w_p, -w_s) ("cars")."""

PARTS: Final = frozenset(
    {
        "electronic",
        "vibrational",
        "rotational",
        "vibrational_rotational",
        "electronic_vibrational",
        "total",
        "instantaneous",
        "unknown",
    }
)
"""Which contributions a third-order value contains. ``vibrational_rotational`` is
Shelton's gamma^vR; ``instantaneous`` a pulse-dependent mixture (Nibbering's)."""

POLARIZABILITY_CORRECTIONS: Final = frozenset(
    {"zero_point_vibrational", "pure_vibrational", "vibrational", "thermal"}
)
"""The parts of a ``polarizability_correction``: the zero-point vibrational average,
the pure vibrational polarizability, their sum, and a thermal (rovibrational)
correction at parameters.temperature_k."""

CONVENTIONS: Final = MappingProxyType(
    {
        "T": "Taylor series (Buckingham-Orr; Shelton & Rice 1994): p = ... + gamma "
        "E^3 / 6",
        "ward_1969": "Ward's per-molecule chi in esu (Ward & New 1969, Ward & Miller "
        "1979, Ward & Elliott 1984, Lundeen et al. 1983): gamma_T = 6 chi",
        "lehmeier_1985": "Lehmeier et al.'s gamma, P = 4 eps0 chi E E E: gamma_T = 24 "
        "gamma",
        "chi_cars_per_amagat": "CARS chi(3) per amagat (Farrow et al., Hahn & Lee): "
        "gamma_T = 2 chi / N0",
        "rosasco_1985_hellwarth": "Rosasco & Hurst's 3 sigma / 24 per amagat, "
        "Hellwarth's sigma = N gamma_T / 3",
        "hellwarth_sigma": "Hellwarth's electronic sigma of P = (1/2) sigma E E^2 "
        "(Gaussian esu)",
        "maker_terhune_c": "Maker & Terhune's c tensor (Gaussian esu), as used by "
        "Hellwarth 1990",
        "bjorklund_1975": "Bjorklund's chi(3) for harmonic generation in esu (Mahon "
        "et al. 1979)",
    }
)
"""The definitions a third-order value may be printed in (``parameters.convention``)."""
