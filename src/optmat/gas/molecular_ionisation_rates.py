# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Strong-field ionisation rates of molecules.

**The rate** (the model ``default``; TC-677 to TC-680) carries the rate of the noble
gases (:mod:`optmat.gas.ionisation_rates`) to molecules. For an unaligned molecule, at
the fields and frequencies of a laser pulse,

.. math::

   W(F, \omega) = N_\text{eff}\, V(F) \sum_{m'} \langle B^2(m') \rangle\,
       \frac{\kappa^{1 - 2n^*}}{2^{|m'|}\,|m'|!\, I_p}\;
       \frac{w^\text{sr}_{m'}(F, \omega)}{|C_{\kappa l}|^2 f(l, m')}\;
       Q(F)\, S(F):

- :math:`\langle B^2(m') \rangle`: the structure factors of the ionising orbital of
  Tong, Zhao and Lin (2002), averaged over the orientations of the molecule
  (:func:`isotropic_structure_factors`), which take the place of an atom's
  :math:`|C_{\kappa l}|^2 f(l, m)`: for one partial wave the two are the same;
- :math:`w^\text{sr}_{m'}`: the short-range rate of Perelomov, Popov and Terent'ev
  (1966, eq. 54) for the projection :math:`m'`, with its dependence on the
  frequency, and :math:`Q = (2F_0/F)^{2n^*}` Perelomov and Popov's static Coulomb
  factor;
- :math:`S = \exp[-\alpha F/(I_p \kappa^3)]`: Tong and Lin's (2005) form of the
  correction at finite field, with :math:`\alpha` fitted here to the average over
  orientations of all-electron static rates (N2, H2);
- :math:`N_\text{eff}`: the number of electrons in the ionising shell, the
  weak-field limit of the many-electron theory (Tolstikhin, Madsen and Morishita
  2014);
- :math:`V`: the average of the rate over the vibrational ground state, where a
  reference gives it (H2), 1 otherwise.

The ionisation energy is the **vertical** one, and the rate that of fixed nuclei at
the neutral's geometry times :math:`V`. The instantaneous and the cycle-averaged rate
are related as for the noble gases (:data:`optmat.gas.ionisation_rates.Averaging`).

**The forms of the literature.** Three rates are in use for molecules in
pulse-propagation codes. Each is given as its source gives it, in SI, so that results
computed with it can be reproduced and the forms compared with one another and with
the references (TC-673):

``moadk``
    The molecular tunnelling rate of Tong, Zhao and Lin (2002). The ionising orbital
    is written by its asymptotic tail, :math:`\Psi^m \to \sum_l C_l\, r^{Z/\kappa -
    1} e^{-\kappa r}\, Y_{lm}` (their eqs. 1-2), and for a molecule at the Euler
    angles :math:`\mathbf R` to a static field :math:`F` (eqs. 8-9)

    .. math::

       w(F, \mathbf R) = \sum_{m'} \frac{B^2(m')}{2^{|m'|}\, |m'|!}\,
           \frac{1}{\kappa^{2Z/\kappa - 1}}
           \left(\frac{2\kappa^3}{F}\right)^{2Z/\kappa - |m'| - 1}
           e^{-2\kappa^3/3F},
       \qquad
       B(m') = \sum_l C_l\, D^l_{m'm}(\mathbf R)\, Q(l, m'),

    with :math:`Q(l, m) = (-1)^m [(2l+1)(l+|m|)!/2(l-|m|)!]^{1/2}`. Averaged over the
    orientations of the molecule the cross terms between partial waves vanish,
    :math:`\langle B^2(m') \rangle = \sum_l C_l^2\, Q^2(l, m')/(2l+1)`
    (:func:`isotropic_structure_factors`).
``ppt``
    The rate of Perelomov, Popov and Terent'ev for an electron with :math:`l = m =
    0` and an **effective charge** :math:`Z_\text{eff}` fitted to measured ion yields
    (Talebpour et al. 1998, 1999): :math:`n^* = Z_\text{eff}/\kappa`, the Hartree
    coefficient with :math:`l^* = n^* - 1`, and a Coulomb factor
    (:func:`optmat.gas.ionisation_rates.coulomb_factor`): the static one, within
    which the charges were fitted (``ppt_talebpour``), or Ilkov et al.'s
    (``ppt_talebpour_ilkov``), within which the charges were not fitted.
``adk``
    The tunnelling rate of a hydrogen-like s electron built on the ionisation energy
    alone: :math:`Z = 1`, :math:`l^* = 0`, no angular factor.

Each is the rate of **one electron**. To apply an occupancy, set
:attr:`MolecularRateModel.electrons` consistently with the coefficients.
The rates of the ``moadk`` and ``adk`` forms are static: their cycle-averaged rate is
the static rate times the factor of the tunnelling limit, :math:`(3F/\pi F_0)^{1/2}`.
The ``ppt`` form is cycle averaged, and its instantaneous rate is that divided by the
same factor.

None of these is validated as the rate of a molecule. Against the orientation
average of all-electron static rates the ``moadk`` form is twice too high for H2, and
for N2 from 0.9 to 6.6 times as the field rises to twice its barrier-suppression
field (:mod:`optmat.gas.molecular_ionisation_validation`).

References
----------
- X. M. Tong, Z. X. Zhao and C. D. Lin, "Theory of molecular tunneling ionization",
  Phys. Rev. A 66, 033402 (2002).
- S.-F. Zhao, C. Jin, A.-T. Le, T. F. Jiang and C. D. Lin, "Determination of structure
  parameters in strong-field tunneling ionization theory of molecules", Phys. Rev. A
  81, 033423 (2010).
- A. Talebpour, J. Yang and S. L. Chin, "Semi-empirical model for the rate of
  tunnel ionization of N2 and O2 molecule in an intense Ti:sapphire laser pulse",
  Opt. Commun. 163, 29 (1999).
- A. Talebpour, S. Larochelle and S. L. Chin, "Suppressed tunnelling ionization of
  the D2 molecule in an intense Ti:sapphire laser pulse", J. Phys. B 31, L49 (1998).
- M. V. Ammosov, N. B. Delone and V. P. Krainov, "Tunnel ionization of complex
  atoms and of atomic ions in an alternating electromagnetic field", Sov. Phys. JETP
  64, 1191 (1986).
- A. M. Perelomov, V. S. Popov and M. V. Terent'ev, "Ionization of atoms in an
  alternating electric field", Sov. Phys. JETP 23, 924 (1966).
- A. M. Perelomov and V. S. Popov, "Ionization of atoms in an alternating electrical
  field. III", Sov. Phys. JETP 25, 336 (1967).
- X. M. Tong and C. D. Lin, "Empirical formula for static field ionization rates of
  atoms and molecules by lasers in the barrier-suppression regime", J. Phys. B 38,
  2593 (2005).
- O. I. Tolstikhin, L. B. Madsen and T. Morishita, "Weak-field asymptotic theory of
  tunneling ionization in many-electron atomic and molecular systems", Phys. Rev. A
  89, 013421 (2014).
- V. P. Majety and A. Scrinzi, "Static field ionization rates for multi-electron
  atoms and small molecules", J. Phys. B 48, 245603 (2015).
- A. Saenz, "Behavior of molecular hydrogen exposed to strong dc, ac, or low-frequency
  laser fields. I. Bond softening and enhanced ionization", Phys. Rev. A 66, 063407
  (2002).
"""

import functools
import math
import tomllib
import warnings
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import Any, ClassVar, Final, Literal, cast, get_args

import numpy as np
from jaxtyping import Float
from numpy.typing import ArrayLike

from optmat import _catalogue
from optmat.constants import AU_ANGULAR_FREQUENCY, AU_ENERGY, AU_FIELD, HBAR, C
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas._tunnelling import isotropic_structure_factors
from optmat.gas.gases import Gas, as_gas
from optmat.gas.ionisation import (
    _energy,  # pyright: ignore[reportPrivateUsage]
    _row,  # pyright: ignore[reportPrivateUsage]
    hartree_coefficient_squared,
    ionisation_parameters,
)
from optmat.gas.ionisation_rates import (
    Averaging,
    ChannelSum,
    CoulombFactor,
    _log_coulomb_factor,
    _log_cycle_averages,  # pyright: ignore[reportPrivateUsage] - one average for both
    _log_short_range_rate,
    _log_short_range_strength,
    _log_sum_exp,  # pyright: ignore[reportPrivateUsage] - scipy's arithmetic, faster
    _species,  # pyright: ignore[reportPrivateUsage] - one rule for both rate kinds
    cycle_average_factor,
    keldysh_parameter,
    ppt_exponent_function,
)
from optmat.identity import Identified, Status, parse_status, source_ids

__all__ = [
    "EVALUATION_VERSION",
    "MOLECULAR_RATES_FILE",
    "MolecularCorrection",
    "MolecularForm",
    "MolecularRateModel",
    "isotropic_structure_factors",
    "moadk_rate_au",
    "molecular_rate_coverage",
    "molecular_rate_model",
    "molecular_rate_models",
]

EVALUATION_VERSION: Final = 2
"""The version of this module's evaluation code, in the identity of every model that
evaluates through it; raised whenever the code changes its output (OEP 1, section 4.8).
"""

MOLECULAR_RATES_FILE: Final = GAS_DATA_DIRECTORY / "molecular_ionisation_rates.toml"
"""The named molecular rate models, each with its reasons."""

MolecularForm = Literal["molecular", "moadk", "ppt", "adk"]
"""The form of a molecular rate: see the module's description."""

MolecularCorrection = Literal["tong-lin", "none"]
"""The correction of the static rate at finite field."""

_VIBRATIONAL_TAIL: Final = 2.0
"""Standard deviations of the ground state to which the vibrational factor follows it.

Within two standard deviations lies 95 % of a Gaussian ground state; beyond, the
anharmonic potential and the curvature of the ionisation energy in the bond length
make the Gaussian cumulant an overestimate. For H2 in a static field the switch is at
0.029 au (3e13 W/cm2), where the rate is 1e-9 of its value at 0.1 au (TC-679).
"""

_MINIMUM_WAVELENGTH: Final = 350e-9
"""Shortest wavelength, m, at which the rate is given without a warning.

The frequency dependence rests on the comparisons of the noble gases, made at 390 nm
and above (TC-653).
"""

_MAXIMUM_FIELD: Final = 1.5
"""The range of a correction, in barrier-suppression fields, where none is given.

The models of ``molecular_ionisation_rates.toml`` give each molecule's
``correction_range``: the largest field its alpha was fitted to or checked against
(TC-828). Within it the rate rises with the field; with N2's alpha of 14 it rises to
2.2 times the barrier-suppression field.
"""


def _log_moadk_rate(
    field_au: Float[np.ndarray, " *shape"],
    ionization_energy_au: float,
    structure_factors: Mapping[int, float],
    charge: float,
) -> Float[np.ndarray, " *shape"]:
    """Return the logarithm of :func:`moadk_rate_au`.

    Parameters
    ----------
    field_au : numpy.ndarray
        Static field, atomic units, > 0.
    ionization_energy_au : float
        Ionisation energy, hartree.
    structure_factors : mapping of int to float
        ``B^2(m')`` by ``m'``.
    charge : float
        Charge of the ion.

    Returns
    -------
    numpy.ndarray
        The logarithm of the rate in atomic units.
    """
    kappa = math.sqrt(2 * ionization_energy_au)
    power = 2 * charge / kappa
    ratio = np.log(2 * kappa**3 / field_au)
    terms = [
        math.log(factor)
        - abs(m) * math.log(2)
        - math.lgamma(abs(m) + 1)
        + (1 - power) * math.log(kappa)
        + (power - abs(m) - 1) * ratio
        # A partial wave printed as 0 (H2's l = 4) gives projections of no weight.
        for m, factor in structure_factors.items()
        if factor > 0
    ]
    return np.asarray(_log_sum_exp(np.stack(terms)) - 2 * kappa**3 / (3 * field_au))


def _log_moadk_term(
    field_au: Float[np.ndarray, " *shape"],
    ionization_energy_au: float,
    factor: float,
    m: int,
    charge: float,
) -> Float[np.ndarray, " *shape"]:
    """Return the logarithm of one projection's term of :func:`moadk_rate_au`.

    Parameters
    ----------
    field_au : numpy.ndarray
        Static field, atomic units, > 0.
    ionization_energy_au : float
        Ionisation energy, hartree.
    factor : float
        ``B^2(m')``, > 0.
    m : int
        The projection ``m'``.
    charge : float
        Charge of the ion.

    Returns
    -------
    numpy.ndarray
        The logarithm of the term in atomic units.
    """
    kappa = math.sqrt(2 * ionization_energy_au)
    power = 2 * charge / kappa
    return np.asarray(
        math.log(factor)
        - abs(m) * math.log(2)
        - math.lgamma(abs(m) + 1)
        + (1 - power) * math.log(kappa)
        + (power - abs(m) - 1) * np.log(2 * kappa**3 / field_au)
        - 2 * kappa**3 / (3 * field_au)
    )


def moadk_rate_au(
    field_au: Float[np.ndarray, " *shape"],
    ionization_energy_au: float,
    structure_factors: Mapping[int, float],
    charge: float = 1.0,
) -> Float[np.ndarray, " *shape"]:
    r"""Return the static molecular tunnelling rate of one electron.

    Tong, Zhao and Lin (2002), eqs. (7) and (9):

    .. math::

       w(F) = \sum_{m'} \frac{B^2(m')}{2^{|m'|}\, |m'|!}\,
           \frac{1}{\kappa^{2Z/\kappa - 1}}
           \left(\frac{2\kappa^3}{F}\right)^{2Z/\kappa - |m'| - 1}
           e^{-2\kappa^3/3F}.

    Parameters
    ----------
    field_au : numpy.ndarray
        Static field, atomic units, > 0.
    ionization_energy_au : float
        Ionisation energy of the orbital, hartree, > 0.
    structure_factors : mapping of int to float
        :math:`B^2(m')` by :math:`m'`, every :math:`m'` of either sign:
        :func:`isotropic_structure_factors` for an unaligned molecule.
    charge : float, optional
        Charge :math:`Z` of the ion left behind.

    Returns
    -------
    numpy.ndarray
        Rate, atomic units, with the shape of ``field_au``.

    Examples
    --------
    Hydrogen 1s (:math:`C_0 = 2`) gives :math:`(4/F) e^{-2/3F}`:

    >>> import numpy as np
    >>> factors = isotropic_structure_factors({0: 2.0})
    >>> w = moadk_rate_au(np.array(0.05), 0.5, factors)
    >>> f"{float(w) / (4 / 0.05 * np.exp(-2 / 0.15)):.12f}"
    '1.000000000000'
    """
    return np.asarray(
        np.exp(
            _log_moadk_rate(
                np.asarray(field_au, dtype=float),
                ionization_energy_au,
                structure_factors,
                charge,
            )
        )
    )


@dataclass(frozen=True)
class MolecularRateModel(Identified):
    """A named form of the ionisation rate of a molecule, at one wavelength.

    Built by :func:`molecular_rate_model` from the curated parameters, or by
    :meth:`from_parameters` from a caller's; :meth:`with_parameters
    <optmat.identity.Identified.with_parameters>` rebuilds it with any of
    :attr:`INPUTS` changed. The orientation average ``structure_factors`` and the
    Hartree coefficient are derived, and an override naming them is refused.

    Attributes
    ----------
    gas : Gas or str
        The molecule, or a custom species' name.
    wavelength : float
        Vacuum wavelength of the field, m; ``math.inf`` for a static field. The
        ``molecular`` and ``ppt`` forms depend on it.
    model : str
        The model's name in ``molecular_ionisation_rates.toml``.
    form : {"molecular", "moadk", "ppt", "adk"}
        The form of the rate.
    ionization_energy : float
        The ionisation energy the rate is evaluated with, J.
    ionization_energy_kind : str
        ``adiabatic``, ``vertical`` or ``not stated``.
    charge : float
        The charge seen by the electron: 1, or the fitted effective charge of the
        ``ppt`` form.
    structure_factors : dict of int to float
        ``<B^2(m')>`` of the ``molecular`` and ``moadk`` forms, averaged over
        orientations; empty for the other forms.
    coefficient_squared : float
        The Hartree coefficient of the ``ppt`` and ``adk`` forms; NaN otherwise.
    coulomb : {"static", "ilkov", "popruzhenko"}
        The Coulomb factor of the ``molecular`` and ``ppt`` forms.
    electrons : float
        The factor on the rate of the one electron: 1 as the sources give the rates;
        the number of electrons in the ionising shell in the ``default`` model.
    static_correction : {"tong-lin", "none"}
        The correction of the static rate at finite field.
    tong_lin_alpha : float
        Its constant; NaN without a correction.
    correction_range : float
        The largest field, in barrier-suppression fields, at which the correction is
        used: above it the rate is held at its value there (:attr:`maximum_field`).
    vibrational_constant : float
        :math:`c` of the vibrational factor :math:`V = e^{c/F^2}`, atomic units of
        field squared; 0 for fixed nuclei.
    averaging : {"consistent", "tunnelling"}
        How the cycle-averaged and the instantaneous rates are related
        (:data:`optmat.gas.ionisation_rates.Averaging`).
    channel_sum : {"integral", "sum"}
        How the sum over the above-threshold channels is taken.
    adiabatic_ionization_energy : float
        The adiabatic ionisation energy, J, where it differs from the energy of the
        rate: the energy of the molecule's ion in its ground state. NaN where not set.
    uncertainty_factor : float
        The factor within which the rate is expected to be right; NaN for the forms
        of the literature.
    sources : dict of str to str
        The dataset or paper of each part.
    reasons : dict of str to str
        Why each choice was made.
    """

    gas: Gas | str
    wavelength: float
    model: str
    form: MolecularForm
    ionization_energy: float
    ionization_energy_kind: str
    charge: float = 1.0
    structure_factors: Mapping[int, float] = field(default_factory=dict)
    coefficient_squared: float = math.nan
    coulomb: CoulombFactor = "ilkov"
    electrons: float = 1.0
    static_correction: MolecularCorrection = "none"
    tong_lin_alpha: float = math.nan
    correction_range: float = _MAXIMUM_FIELD
    vibrational_constant: float = 0.0
    averaging: Averaging = "tunnelling"
    channel_sum: ChannelSum = "integral"
    adiabatic_ionization_energy: float = math.nan
    uncertainty_factor: float = math.nan
    sources: Mapping[str, str] = field(default_factory=dict, repr=False)
    reasons: Mapping[str, str] = field(default_factory=dict, repr=False)

    IDENTITY_KIND: ClassVar[str] = "ionisation_rate"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = (
        "optmat.gas.ionisation_rates",
        "optmat.gas.molecular_ionisation_rates",
    )
    INPUTS: ClassVar[tuple[str, ...]] = (
        "adiabatic_ionization_energy",
        "averaging",
        "channel_sum",
        "charge",
        "correction_range",
        "coulomb",
        "electrons",
        "form",
        "ionization_energy",
        "ionization_energy_kind",
        "static_correction",
        "structure_coefficients",
        "tong_lin_alpha",
        "uncertainty_factor",
        "vibrational_constant",
    )
    DERIVED: ClassVar[tuple[str, ...]] = ("coefficient_squared", "structure_factors")

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the model's name, gas and source records.

        Returns
        -------
        tuple
            ``("molecular_rate/<model>/<gas>@<wavelength>", gas, sources)``, the
            wavelength in nm or ``static``.
        """
        at = (
            "static"
            if math.isinf(self.wavelength)
            else f"{self.wavelength * 1e9:.6g}nm"
        )
        name = f"molecular_rate/{self.model}/{self.gas}@{at}"
        return name, str(self.gas), source_ids(*self.sources.values())

    def __post_init__(self) -> None:
        """Check the wavelength, the energy, the factor and the correction.

        Raises
        ------
        ValueError
            If one of the first four is not positive, the correction has no constant
            or no positive range, or the vibrational constant is negative.
        """
        super().__post_init__()
        for name in ("wavelength", "ionization_energy", "charge", "electrons"):
            if not getattr(self, name) > 0:
                msg = f"{name} must be positive, got {getattr(self, name)!r}"
                raise ValueError(msg)
        if self.static_correction == "tong-lin" and not self.correction_range > 0:
            msg = (
                f"the correction's range must be positive, got "
                f"{self.correction_range!r}"
            )
            raise ValueError(msg)
        if self.static_correction == "tong-lin" and not self.tong_lin_alpha > 0:
            msg = f"Tong and Lin's alpha must be positive, got {self.tong_lin_alpha!r}"
            raise ValueError(msg)
        if self.vibrational_constant < 0:
            msg = "the vibrational constant must not be negative"
            raise ValueError(msg)

    @property
    def ionization_energy_au(self) -> float:
        """The ionisation energy in hartree.

        Returns
        -------
        float
            ``I_p / E_h``.
        """
        return self.ionization_energy / AU_ENERGY

    @property
    def omega_au(self) -> float:
        """The angular frequency of the field in atomic units.

        Returns
        -------
        float
            ``hbar omega / E_h``; 0 for a static field.
        """
        if math.isinf(self.wavelength):
            return 0.0
        return 2 * math.pi * HBAR * C / (self.wavelength * AU_ENERGY)

    @property
    def barrier_suppression_field(self) -> float:
        """The field at which the barrier of a Coulomb potential is suppressed, V/m.

        Returns
        -------
        float
            ``I_p^2 / 4Z`` in atomic units, with ``Z = 1``: beyond it no tunnelling
            rate holds, and none of these forms has a correction.
        """
        return self.ionization_energy_au**2 / 4 * AU_FIELD

    @property
    def maximum_field(self) -> float:
        """The largest field at which the rate is given without a warning, V/m.

        Returns
        -------
        float
            With Tong and Lin's correction, :attr:`correction_range` times the
            barrier-suppression field: the largest field its alpha was fitted to or
            checked against (TC-828). The barrier-suppression field without a
            correction.
        """
        if self.static_correction == "tong-lin":
            return self.correction_range * self.barrier_suppression_field
        return self.barrier_suppression_field

    @property
    def citation(self) -> str:
        """The papers of the form in force.

        Returns
        -------
        str
            One line naming them.
        """
        structure = self.sources.get("structure")
        energy = self.sources.get("ionization_energy")
        forms = {
            "molecular": "Perelomov, Popov & Terent'ev, Sov. Phys. JETP 23, 924 "
            "(1966), eq. (54), with the structure factors of Tong, Zhao & Lin, Phys. "
            "Rev. A 66, 033402 (2002), averaged over orientations (coefficients: "
            f"{structure}), "
            f"the {self.coulomb} Coulomb factor, Tong & Lin's correction (alpha "
            f"{self.tong_lin_alpha:g}, optmat's fit), {self.electrons:g} electrons"
            + (
                f", vibrational factor exp({self.vibrational_constant:g}/F^2)"
                if self.vibrational_constant
                else ""
            ),
            "moadk": "Tong, Zhao & Lin, Phys. Rev. A 66, 033402 (2002), eqs. (8)-(9), "
            f"averaged over orientations; coefficients: {structure}",
            "ppt": "Perelomov, Popov & Terent'ev, Sov. Phys. JETP 23, 924 (1966), eq. "
            f"(54), l = m = 0, Coulomb factor {self.coulomb}; effective charge "
            f"{self.charge:g}: {self.sources.get('charge')}",
            "adk": "Ammosov, Delone & Krainov, Sov. Phys. JETP 64, 1191 (1986), for a "
            "hydrogen-like s electron (Z = 1, l* = 0)",
        }
        return (
            f"{forms[self.form]}; ionisation energy "
            f"({self.ionization_energy_kind}): {energy}; "
            f"electrons {self.electrons:g}"
        )

    def keldysh(self, field_strength: ArrayLike) -> Float[np.ndarray, " *shape"]:
        """Return the Keldysh parameter at a field.

        Parameters
        ----------
        field_strength : array_like
            Field, V/m, > 0.

        Returns
        -------
        numpy.ndarray
            ``omega kappa / F``; 0 for a static field.
        """
        strength = np.abs(np.asarray(field_strength, dtype=float)) / AU_FIELD
        return keldysh_parameter(strength, self.omega_au, self.ionization_energy_au)

    def _log_corrections(
        self, field_au: Float[np.ndarray, " n"], omega_au: float
    ) -> Float[np.ndarray, " n"]:
        r"""Return the logarithm of the static correction and the vibrational factor.

        **The vibrational factor** is the average of the rate over the vibrational
        ground state, over the rate at the neutral's bond length. The bond length
        enters the rate through the ionisation energy, and the rate depends on the
        energy chiefly through its exponent, :math:`(2F_0/3F)\, g(\gamma)`
        (:func:`optmat.gas.ionisation_rates.ppt_exponent_function`), whose
        derivative is :math:`-(2\kappa/F)\, h(\gamma)`, :math:`h = g + \gamma g'/3`.
        For a Gaussian spread :math:`\sigma_I` of the energy over the ground state,

        .. math::

           \ln V = \frac{\sigma_I^2}{2} \left(\frac{2\kappa}{F}\right)^2 h(\gamma)^2
                 = \frac{c}{F^2}\, h(\gamma)^2 = \frac{x^2}{2},
           \qquad x = \frac{\sqrt{2c}\, h(\gamma)}{F},

        exact for a Gaussian and a linear dependence of the energy on the bond
        length. :math:`h \to 1` in the tunnelling limit, where :math:`c` is fitted,
        and :math:`h \approx \ln(2\gamma)/\gamma` in the multiphoton limit, where
        the rate depends on the energy only through the number of photons.
        :math:`x` is the number of standard deviations out in the ground state at
        which the average is dominated; beyond :data:`_VIBRATIONAL_TAIL` of them
        neither the Gaussian nor the linear dependence holds, and the factor is
        continued as the contribution of that bond length alone, :math:`\ln V = kx -
        k^2/2`: continuous, with a continuous slope, so that the rate still falls to
        zero with the field (TC-679).

        Parameters
        ----------
        field_au : numpy.ndarray
            Field, atomic units, > 0.
        omega_au : float
            Angular frequency, atomic units; 0 for the static rate.

        Returns
        -------
        numpy.ndarray
            ``ln S + ln V``.
        """
        energy = self.ionization_energy_au
        log = np.zeros_like(field_au)
        if self.vibrational_constant > 0:
            if omega_au == 0:
                weight = np.ones_like(field_au)
            else:
                gamma = keldysh_parameter(field_au, omega_au, energy)
                # g'(gamma) by a central difference: g is smooth on the scale of gamma.
                step = 1e-4 * gamma
                slope = (
                    ppt_exponent_function(gamma + step)
                    - ppt_exponent_function(gamma - step)
                ) / (2 * step)
                weight = ppt_exponent_function(gamma) + gamma * slope / 3
            x = math.sqrt(2 * self.vibrational_constant) * weight / field_au
            tail = _VIBRATIONAL_TAIL
            log = np.where(x <= tail, x**2 / 2, tail * x - tail**2 / 2)
        if self.static_correction == "tong-lin":
            # Tong and Lin (2005), eq. (2), with Z = 1.
            log = log - self.tong_lin_alpha * field_au / (energy * (2 * energy) ** 1.5)
        return np.asarray(log)

    def _log_terms(
        self, field_au: Float[np.ndarray, " n"]
    ) -> tuple[Float[np.ndarray, " n"], Float[np.ndarray, " n"]]:
        r"""Return the logarithms of the static rate and of the rate with its frequency.

        As :meth:`optmat.gas.ionisation_rates.IonisationRateModel._log_terms`, with a
        term for each projection :math:`m'` of the orbital on the field, whose
        :math:`\langle B^2(m') \rangle \kappa^{1-2n^*}/(2^{|m'|} |m'|!\, I_p)` takes
        the place of an atom's :math:`|C_{\kappa l}|^2 f(l, m)`.

        Parameters
        ----------
        field_au : numpy.ndarray
            Field, atomic units, > 0.

        Returns
        -------
        tuple of numpy.ndarray
            The logarithms of the static rate and of the cycle-averaged rate over its
            tunnelling limit's factor, both with the corrections, one electron.
        """
        return self._log_term(0, field_au), self._log_term(1, field_au)

    def _log_term(
        self, which: int, field_au: Float[np.ndarray, " n"]
    ) -> Float[np.ndarray, " n"]:
        r"""Return one of the two logarithms of :meth:`_log_terms`, evaluating it alone.

        The projections :math:`\pm m'` of an orientation average have the same
        :math:`B^2(m')`, and a term depends on :math:`|m'|` only, so each distinct term
        is evaluated once and enters the sum once for each projection.

        Parameters
        ----------
        which : int
            0, the static rate, or 1, the rate with its frequency.
        field_au : numpy.ndarray
            Field, atomic units, > 0.

        Returns
        -------
        numpy.ndarray
            The logarithm.
        """
        energy = self.ionization_energy_au
        omega = 0.0 if which == 0 else self.omega_au
        kappa = math.sqrt(2 * energy)
        power = 2 * self.charge / kappa
        terms: list[np.ndarray] = []
        evaluated: dict[tuple[int, float], np.ndarray] = {}
        for m, factor in self.structure_factors.items():
            # A partial wave printed as 0 (H2's l = 4) gives projections of no weight.
            if factor <= 0:
                continue
            key = (abs(m), factor)
            if key not in evaluated:
                if omega == 0:
                    evaluated[key] = _log_moadk_term(
                        field_au, energy, factor, m, self.charge
                    )
                else:
                    log_strength = (
                        math.log(factor)
                        + (1 - power) * math.log(kappa)
                        - abs(m) * math.log(2)
                        - math.lgamma(abs(m) + 1)
                        - math.log(energy)
                    )
                    evaluated[key] = (
                        _log_short_range_strength(
                            field_au, omega, energy, log_strength, m, self.channel_sum
                        )
                        + _log_coulomb_factor(
                            field_au, omega, energy, self.coulomb, self.charge
                        )
                        - np.log(cycle_average_factor(field_au, energy))
                    )
            terms.append(evaluated[key])
        return np.asarray(
            _log_sum_exp(np.stack(terms)) + self._log_corrections(field_au, omega)
        )

    def _log_molecular_rate(
        self, field_au: Float[np.ndarray, " n"], *, cycle_averaged: bool
    ) -> Float[np.ndarray, " n"]:
        r"""Return the logarithm of the rate of the ``molecular`` form.

        The cycle-averaged rate is :math:`N_\text{eff} B a_0` and the instantaneous
        rate :math:`N_\text{eff} B a_0/a`, with :math:`a` and :math:`a_0` the averages
        over a cycle of the rate and of its static limit (TC-662).

        Parameters
        ----------
        field_au : numpy.ndarray
            Field, atomic units, > 0.
        cycle_averaged : bool
            The rate against the amplitude, or the instantaneous rate.

        Returns
        -------
        numpy.ndarray
            The logarithm of the rate.
        """
        centre = self._log_terms(field_au)
        log_static, log_enhanced = centre
        log_electrons = math.log(self.electrons)
        if self.averaging == "tunnelling":
            log_rate = log_enhanced + log_electrons
            if cycle_averaged:
                log_rate = log_rate + np.log(
                    cycle_average_factor(field_au, self.ionization_energy_au)
                )
            return np.asarray(log_rate)
        if not cycle_averaged and self.omega_au == 0:
            return np.asarray(log_static + log_electrons)
        log_static_average, log_average = _log_cycle_averages(
            self._log_terms, field_au, centre=centre, log_term=self._log_term
        )
        log_rate = log_enhanced + log_electrons + log_static_average
        if not cycle_averaged:
            log_rate = log_rate - log_average
        return np.asarray(log_rate)

    def _log_rate_au(
        self, field_au: Float[np.ndarray, " n"], *, cycle_averaged: bool
    ) -> Float[np.ndarray, " n"]:
        """Return the logarithm of the rate in atomic units.

        Parameters
        ----------
        field_au : numpy.ndarray
            Field, atomic units, > 0.
        cycle_averaged : bool
            The rate against the amplitude, or the instantaneous rate.

        Returns
        -------
        numpy.ndarray
            The logarithm of the rate.
        """
        if self.form == "molecular":
            return self._log_molecular_rate(field_au, cycle_averaged=cycle_averaged)
        energy = self.ionization_energy_au
        # log (3F / pi F_0)^(1/2): the cycle average of the tunnelling limit.
        average = 0.5 * np.log(3 * field_au / (math.pi * (2 * energy) ** 1.5))
        if self.form == "ppt":
            log_rate = _log_short_range_rate(
                field_au,
                self.omega_au,
                energy,
                self.coefficient_squared,
                0,
                0,
                "integral",
            ) + _log_coulomb_factor(
                field_au, self.omega_au, energy, self.coulomb, self.charge
            )
            if not cycle_averaged:
                log_rate = log_rate - average
        else:
            if self.form == "moadk":
                log_rate = _log_moadk_rate(
                    field_au, energy, self.structure_factors, self.charge
                )
            else:
                power = 2 * self.charge / math.sqrt(2 * energy) - 1
                ratio = 2 * (2 * energy) ** 1.5 / field_au
                log_rate = (
                    math.log(self.coefficient_squared * energy)
                    + power * np.log(ratio)
                    - ratio / 3
                )
            if cycle_averaged:
                log_rate = log_rate + average
        return np.asarray(log_rate + math.log(self.electrons))

    def _rate(
        self, field_strength: ArrayLike, *, cycle_averaged: bool
    ) -> Float[np.ndarray, " *shape"]:
        """Return a rate in SI, zero at zero field, warning above the range.

        Above :attr:`maximum_field` a rate with Tong and Lin's correction is held at its
        value there (TC-816, TC-828); one without keeps its formula.

        Parameters
        ----------
        field_strength : array_like
            Field, V/m.
        cycle_averaged : bool
            The rate against the amplitude, or the instantaneous rate.

        Returns
        -------
        numpy.ndarray
            Rate per molecule, 1/s.
        """
        strength = np.asarray(np.abs(np.asarray(field_strength, dtype=float)))
        if np.any(strength > self.maximum_field):
            warnings.warn(
                f"fields above {self.maximum_field:.3g} V/m are beyond the range of "
                f"the {self.form} form of {self.gas}'s rate"
                + (
                    f" ({self.correction_range:g} times the barrier-suppression field)"
                    if self.static_correction == "tong-lin"
                    else " (the barrier-suppression field: it has no correction)"
                )
                + (
                    "; the rate is held at its value there"
                    if self.static_correction == "tong-lin"
                    else ""
                ),
                ModelRangeWarning,
                stacklevel=3,
            )
            if self.static_correction == "tong-lin":
                strength = np.asarray(np.minimum(strength, self.maximum_field))
        if self.form == "molecular" and self.wavelength < _MINIMUM_WAVELENGTH:
            warnings.warn(
                f"the rate's dependence on the frequency is not tested below "
                f"{_MINIMUM_WAVELENGTH * 1e9:.0f} nm",
                ModelRangeWarning,
                stacklevel=3,
            )
        if self.form == "ppt" and self.omega_au == 0:
            msg = "the ppt form needs a finite wavelength"
            raise ValueError(msg)
        flat = strength.ravel() / AU_FIELD
        rate = np.zeros_like(flat)
        positive = flat > 0
        if np.any(positive):
            with np.errstate(under="ignore"):
                rate[positive] = np.exp(
                    self._log_rate_au(flat[positive], cycle_averaged=cycle_averaged)
                )
        return (rate * AU_ANGULAR_FREQUENCY).reshape(strength.shape)

    def cycle_averaged_rate(self, amplitude: ArrayLike) -> Float[np.ndarray, " *shape"]:
        r"""Return the rate averaged over a cycle of a field of constant amplitude.

        For the ``moadk`` and ``adk`` forms, the static rate at the amplitude times
        :math:`(3F/\pi F_0)^{1/2}`.

        Parameters
        ----------
        amplitude : array_like
            Amplitude of the field, V/m; its sign is ignored.

        Returns
        -------
        numpy.ndarray
            Rate per molecule, 1/s, with the shape of ``amplitude``; 0 at zero field.

        Raises
        ------
        ValueError
            If the form is ``ppt`` and the wavelength is infinite.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.

        Examples
        --------
        >>> model = molecular_rate_model("N2", 800e-9, "moadk_zhao_2010")
        >>> f"{float(model.cycle_averaged_rate(2.74e10)):.3e}"
        '2.077e+11'
        """
        return self._rate(amplitude, cycle_averaged=True)

    def instantaneous_rate(
        self, field_strength: ArrayLike
    ) -> Float[np.ndarray, " *shape"]:
        r"""Return the rate at the instantaneous value of a carrier-resolved field.

        For the ``ppt`` form, the cycle-averaged rate at an amplitude equal to the
        instantaneous field, divided by :math:`(3F/\pi F_0)^{1/2}`.

        Parameters
        ----------
        field_strength : array_like
            The field, V/m; its sign is ignored.

        Returns
        -------
        numpy.ndarray
            Rate per molecule, 1/s, with the shape of ``field_strength``; 0 at zero
            field.

        Raises
        ------
        ValueError
            If the form is ``ppt`` and the wavelength is infinite.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.

        Examples
        --------
        >>> model = molecular_rate_model("N2", 800e-9, "moadk_zhao_2010")
        >>> f"{float(model.instantaneous_rate(2.74e10)):.3e}"
        '1.019e+12'
        """
        return self._rate(field_strength, cycle_averaged=False)

    def channel_rates(self, field_strength: ArrayLike) -> Float[np.ndarray, " c ..."]:
        """Return the rate of each channel: one, the instantaneous rate.

        For a linearly polarised field (quasi-static along an elliptical one, which
        optmat does not validate).

        Parameters
        ----------
        field_strength : array_like
            The field, V/m; its sign is ignored.

        Returns
        -------
        numpy.ndarray
            Rates per molecule, 1/s, shape ``(1,) + shape(field)``.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.
        """
        return self.instantaneous_rate(field_strength)[None, ...]

    def energy_loss_rate(
        self, field_strength: ArrayLike
    ) -> Float[np.ndarray, " *shape"]:
        """Return the energy the ionisation takes from the field, ``I_p Gamma``.

        For a linearly polarised field.

        Parameters
        ----------
        field_strength : array_like
            The field, V/m; its sign is ignored.

        Returns
        -------
        numpy.ndarray
            The energy-loss rate, W per molecule, the shape of the field.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.
        """
        return self.ionization_energy * self.instantaneous_rate(field_strength)

    def energy_per_event(
        self, field_strength: ArrayLike
    ) -> Float[np.ndarray, " *shape"]:
        """Return the ionisation energy of an event: the one the exponent uses.

        Parameters
        ----------
        field_strength : array_like
            The field, V/m.

        Returns
        -------
        numpy.ndarray
            :attr:`ionization_energy`, J, the shape of the field.
        """
        return np.full(np.shape(field_strength), self.ionization_energy)

    @property
    def energy_kind(self) -> str:
        """The kind of the ionisation energy the exponent uses.

        Returns
        -------
        str
            ``vertical``, ``adiabatic`` or ``not stated``, as the curated dataset
            states it (:attr:`ionization_energy_kind`).
        """
        return self.ionization_energy_kind

    def at_wavelength(self, wavelength: float) -> "MolecularRateModel":
        """Return the same model at another wavelength.

        Parameters
        ----------
        wavelength : float
            Vacuum wavelength, m; ``math.inf`` for a static field.

        Returns
        -------
        MolecularRateModel
            The model at that wavelength.
        """
        return replace(self, wavelength=wavelength)

    @classmethod
    def from_parameters(  # noqa: PLR0913 - the model's parameters
        cls,
        species: Gas | str,
        wavelength: float,
        *,
        form: MolecularForm,
        ionization_energy: float,
        ionization_energy_kind: str = "not stated",
        structure_coefficients: Mapping[int, float] | None = None,
        structure_factors: Mapping[int, float] | None = None,
        charge: float = 1.0,
        electrons: float = 1.0,
        coulomb: CoulombFactor = "ilkov",
        static_correction: MolecularCorrection = "none",
        tong_lin_alpha: float = math.nan,
        correction_range: float = _MAXIMUM_FIELD,
        vibrational_constant: float = 0.0,
        averaging: Averaging = "tunnelling",
        channel_sum: ChannelSum = "integral",
        adiabatic_ionization_energy: float = math.nan,
        uncertainty_factor: float = math.nan,
        citation: str = "",
        model: str = "user",
        sources: Mapping[str, str] | None = None,
        reasons: Mapping[str, str] | None = None,
    ) -> "MolecularRateModel":
        """Return a molecular rate model from explicit parameters.

        The curated models are built through it. The orientation average of the
        ``molecular`` and ``moadk`` forms is computed from the structure
        coefficients, and the Hartree coefficient of the ``ppt`` and ``adk`` forms
        from the energy and the charge.

        Parameters
        ----------
        species : Gas or str
            A curated molecule, or a custom species' name, the identity's material.
        wavelength : float
            Vacuum wavelength, m; ``math.inf`` for a static field.
        form : {"molecular", "moadk", "ppt", "adk"}
            The form of the rate.
        ionization_energy : float
            J.
        ionization_energy_kind : str, optional
            ``adiabatic``, ``vertical`` or ``not stated``.
        structure_coefficients : Mapping or None, optional
            ``C_l`` by partial wave (``sqrt(sum_m |C_lm|^2)``), in Tong et al.'s
            normalisation, for the ``molecular`` and ``moadk`` forms.
        structure_factors : Mapping or None, optional
            The orientation average itself, as a curated model stores it, in place
            of the coefficients.
        charge : float, optional
            The charge seen by the electron (the ``ppt`` form's effective charge).
        electrons : float, optional
            The factor on the one electron's rate.
        coulomb : str, optional
            The Coulomb factor.
        static_correction : {"tong-lin", "none"}, optional
            The correction at finite field, with ``tong_lin_alpha``.
        tong_lin_alpha : float, optional
            Its constant.
        correction_range : float, optional
            The largest field at which it is used, in barrier-suppression fields, 1.5
            by default; above it the rate is held.
        vibrational_constant : float, optional
            The vibrational factor's constant, atomic units of field squared.
        averaging : str, optional
            ``consistent`` or ``tunnelling``.
        channel_sum : str, optional
            ``integral`` or ``sum``.
        adiabatic_ionization_energy : float, optional
            J, where it differs from the rate's.
        uncertainty_factor : float, optional
            A stated uncertainty factor; NaN (none) by default.
        citation : str, optional
            Where the parameters come from.
        model : str, optional
            The model's name, ``"user"`` by default.
        sources : Mapping or None, optional
            The source of each part.
        reasons : Mapping or None, optional
            Why each choice was made.

        Returns
        -------
        MolecularRateModel
            The model, status ``user``.

        Raises
        ------
        ValueError
            If the structure is missing, given twice, or given to a form without
            one, or as the model's checks do.

        Examples
        --------
        >>> from optmat.constants import E_CHARGE
        >>> model = MolecularRateModel.from_parameters(
        ...     "N2*", 800e-9, form="moadk", ionization_energy=15.58 * E_CHARGE,
        ...     structure_coefficients={0: 2.02, 2: 0.78, 4: 0.04})
        >>> model.identity.material, sorted(model.structure_factors)
        ('N2*', [-4, -3, -2, -1, 0, 1, 2, 3, 4])
        """
        structured = form in {"molecular", "moadk"}
        given = (structure_coefficients is not None) + (structure_factors is not None)
        if structured != (given == 1) or given > 1:
            msg = (
                f"the {form} form takes the structure coefficients (or factors) once"
                if structured
                else f"the {form} form takes no structure coefficients"
            )
            raise ValueError(msg)
        factors = (
            isotropic_structure_factors(structure_coefficients)
            if structure_coefficients is not None
            else dict(structure_factors or {})
        )
        return cls(
            gas=_species(species),
            wavelength=wavelength,
            model=model,
            form=form,
            ionization_energy=ionization_energy,
            ionization_energy_kind=ionization_energy_kind,
            charge=charge,
            structure_factors=factors,
            coefficient_squared=_coefficient(
                form, ionization_energy / AU_ENERGY, charge
            ),
            coulomb=coulomb,
            electrons=electrons,
            static_correction=static_correction,
            tong_lin_alpha=tong_lin_alpha,
            correction_range=correction_range,
            vibrational_constant=vibrational_constant,
            averaging=averaging,
            channel_sum=channel_sum,
            adiabatic_ionization_energy=adiabatic_ionization_energy,
            uncertainty_factor=uncertainty_factor,
            sources=({"citation": citation} if citation else {}) | dict(sources or {}),
            reasons=dict(reasons or {}),
        )

    def _override(self, changes: Mapping[str, Any]) -> "MolecularRateModel":
        """Return the model with some inputs changed, through :meth:`from_parameters`.

        Parameters
        ----------
        changes : Mapping
            Inputs to new values.

        Returns
        -------
        MolecularRateModel
            The rebuilt model; its structure factors kept unless the coefficients
            change.
        """
        inputs: dict[str, Any] = {
            name: getattr(self, name)
            for name in self.INPUTS
            if name != "structure_coefficients"
        } | dict(changes)
        if "structure_coefficients" not in changes and self.structure_factors:
            inputs["structure_factors"] = self.structure_factors
        return type(self).from_parameters(
            self.gas,
            self.wavelength,
            model=self.model,
            sources=self.sources,
            reasons=self.reasons,
            **inputs,
        )

    def _stated_uncertainty(self) -> Mapping[str, Any]:
        """Return the stated uncertainty an override drops.

        Returns
        -------
        Mapping
            ``uncertainty_factor`` to NaN.
        """
        return {"uncertainty_factor": math.nan}

    def with_ionization_energy(
        self, ionization_energy: float, kind: str
    ) -> "MolecularRateModel":
        """Return the same model evaluated with another ionisation energy.

        An override (:meth:`~optmat.identity.Identified.with_parameters`): the
        Hartree coefficient of the ``ppt`` and ``adk`` forms is recomputed; the
        structure coefficients and the effective charge are not, though each was
        determined with an energy of its own.

        Parameters
        ----------
        ionization_energy : float
            The energy, J, > 0.
        kind : str
            ``adiabatic``, ``vertical`` or a description of where it is from.

        Returns
        -------
        MolecularRateModel
            The model with that energy, status ``user``.
        """
        return self.with_parameters(
            ionization_energy=ionization_energy, ionization_energy_kind=kind
        )


def _coefficient(form: str, ionization_energy_au: float, charge: float) -> float:
    """Return the Hartree coefficient of a form, NaN for the molecular form.

    Parameters
    ----------
    form : str
        ``moadk``, ``ppt`` or ``adk``.
    ionization_energy_au : float
        Ionisation energy, hartree.
    charge : float
        The charge seen by the electron.

    Returns
    -------
    float
        ``|C_{n* l*}|^2`` with ``l* = n* - 1`` (``ppt``: Ilkov et al.'s choice for a
        ground state) or ``l* = 0`` (``adk``: the hydrogen-like form).
    """
    n_star = charge / math.sqrt(2 * ionization_energy_au)
    if form == "ppt":
        return hartree_coefficient_squared(n_star, n_star - 1)
    if form == "adk":
        return hartree_coefficient_squared(n_star, 0.0)
    return math.nan


def _text(table: object) -> dict[str, str]:
    """Return a table of reasons as strings.

    Parameters
    ----------
    table : object
        A mapping of ``molecular_ionisation_rates.toml``, or anything else.

    Returns
    -------
    dict of str to str
        The table's entries, or an empty dict.
    """
    if not isinstance(table, Mapping):
        return {}
    return {str(k): str(v) for k, v in dict(table).items()}  # pyright: ignore[reportUnknownArgumentType, reportUnknownVariableType]


def _build_molecular(
    name: str, model: Mapping[str, object], gas: Gas, entry: Mapping[str, object]
) -> MolecularRateModel:
    """Assemble the ``molecular`` form of one molecule from its entry.

    Parameters
    ----------
    name : str
        The model's name.
    model : mapping
        Its table of ``molecular_ionisation_rates.toml``.
    gas : Gas
        The molecule.
    entry : mapping
        The molecule's table within the model.

    Returns
    -------
    MolecularRateModel
        The model at 800 nm.

    Raises
    ------
    DataError
        If a named tabulation or energy does not exist.
    """
    parameters = ionisation_parameters(gas)
    tabulation = str(entry.get("structure", model["structure"]))
    if tabulation not in parameters.structure:
        msg = f"{name}: no {tabulation} structure coefficients for {gas.value}"
        raise DataError(msg)
    structure = parameters.structure[tabulation]
    reference = str(entry["ionization_energy"])
    by_key = _catalogue.by_key()
    dataset, row = _row(by_key, reference)
    energy, _ = _energy(dataset, row)
    kind = str(dataset.parameters.get("kind", "not stated"))
    sources = {
        "structure": structure.source,
        "ionization_energy": reference,
        "adiabatic_ionization_energy": parameters.ionization_energy_source,
        "electrons": f"molecular_ionisation_rates.toml ({name})",
        "static_correction": f"molecular_ionisation_rates.toml ({name}): optmat's fit",
    }
    correction = str(model.get("static_correction", "tong-lin"))
    averaging = str(model.get("averaging", "consistent"))
    channel_sum = str(model.get("channel_sum", "integral"))
    if (
        correction not in get_args(MolecularCorrection)
        or averaging not in get_args(Averaging)
        or channel_sum not in get_args(ChannelSum)
    ):
        msg = f"{name}: unknown correction, averaging or channel sum"
        raise DataError(msg)
    return MolecularRateModel.from_parameters(
        gas,
        800e-9,
        model=name,
        form="molecular",
        ionization_energy=energy,
        ionization_energy_kind=kind,
        # The average over orientations needs only sum_m |C_lm|^2 per l.
        structure_coefficients={
            ell: math.sqrt(q) for ell, q in structure.partial_wave_squares.items()
        },
        coulomb=cast("CoulombFactor", str(model.get("coulomb", "static"))),
        electrons=float(str(entry["electrons"])),
        static_correction=cast("MolecularCorrection", correction),
        tong_lin_alpha=float(str(entry.get("tong_lin_alpha", "nan"))),
        correction_range=float(str(entry.get("correction_range", _MAXIMUM_FIELD))),
        vibrational_constant=float(str(entry.get("vibrational_constant", 0.0))),
        averaging=cast("Averaging", averaging),
        channel_sum=cast("ChannelSum", channel_sum),
        adiabatic_ionization_energy=parameters.ionization_energy,
        uncertainty_factor=float(str(entry.get("uncertainty_factor", "nan"))),
        sources=sources,
        reasons={**_text(model.get("reasons")), **_text(entry.get("reasons"))},
    )


def _build(name: str, model: Mapping[str, object], gas: Gas) -> MolecularRateModel:
    """Assemble one model of one molecule.

    Parameters
    ----------
    name : str
        The model's name.
    model : mapping
        Its table of ``molecular_ionisation_rates.toml``.
    gas : Gas
        The molecule.

    Returns
    -------
    MolecularRateModel
        The model at 800 nm.

    Raises
    ------
    DataError
        If the form is unknown or the molecule lacks what it needs.
    """
    parameters = ionisation_parameters(gas)
    if parameters.orbital_l is not None:
        msg = f"{gas.value} is an atom: see optmat.gas.ionisation_rates"
        raise DataError(msg)
    form = str(model["form"])
    coulomb = str(model.get("coulomb", "ilkov"))
    if form not in get_args(MolecularForm) or coulomb not in get_args(CoulombFactor):
        msg = f"{name}: unknown form {form!r} or Coulomb factor {coulomb!r}"
        raise DataError(msg)
    sources = {"ionization_energy": parameters.ionization_energy_source}
    coefficients: Mapping[int, float] | None = None
    charge = 1.0
    if form == "moadk":
        tabulation = str(model["structure"])
        if tabulation not in parameters.structure:
            msg = f"{name}: no {tabulation} structure coefficients for {gas.value}"
            raise DataError(msg)
        structure = parameters.structure[tabulation]
        coefficients = structure.coefficients
        sources["structure"] = structure.source
    elif form == "ppt":
        if parameters.effective_charge is None:
            msg = f"{name}: no effective charge is curated for {gas.value}"
            raise DataError(msg)
        charge = parameters.effective_charge
        sources["charge"] = parameters.effective_charge_source
    reasons = model.get("reasons", {})
    return MolecularRateModel.from_parameters(
        gas,
        800e-9,
        model=name,
        form=cast("MolecularForm", form),
        ionization_energy=parameters.ionization_energy,
        ionization_energy_kind=parameters.ionization_energy_kind,
        charge=charge,
        structure_coefficients=coefficients,
        coulomb=cast("CoulombFactor", coulomb),
        sources=sources,
        reasons={str(k): str(v) for k, v in dict(reasons).items()}
        if isinstance(reasons, Mapping)
        else {},
    )


@functools.cache
def _models_table() -> dict[str, Mapping[str, Any]]:
    """Read the models of ``molecular_ionisation_rates.toml``.

    Returns
    -------
    dict of str to mapping
        Each model's table, by name.
    """
    with MOLECULAR_RATES_FILE.open("rb") as handle:
        models = tomllib.load(handle)["models"]
    return {str(name): model for name, model in models.items()}


def _covered(name: str) -> dict[Gas, str]:
    """Return the molecules a model covers.

    Parameters
    ----------
    name : str
        A model of ``molecular_ionisation_rates.toml``.

    Returns
    -------
    dict of Gas to str
        Each molecule and its name in the file, in the file's order.
    """
    return {as_gas(gas): str(gas) for gas in _models_table()[name]["gases"]}


@functools.cache
def _rate_model(name: str, gas: Gas) -> MolecularRateModel:
    """Build one model of one molecule, reading only the datasets it cites.

    Parameters
    ----------
    name : str
        A model of ``molecular_ionisation_rates.toml``.
    gas : Gas
        A molecule it covers.

    Returns
    -------
    MolecularRateModel
        The model at 800 nm, with its status.
    """
    model = _models_table()[name]
    # optmat's rate is validated where the data do not say otherwise; the forms of
    # the literature are kept for comparison (TC-673).
    default = "validated" if name == "default" else "comparison"
    if model["form"] == "molecular":
        listed = _covered(name)[gas]
        entry = model["gases"][listed]
        return replace(
            _build_molecular(name, model, gas, entry),
            **_status(entry, f"{name}/{listed}", default),
        )
    return replace(_build(name, model, gas), **_status(model, name, default))


def molecular_rate_coverage() -> dict[str, tuple[Gas, ...]]:
    """Return the gases each molecular rate model covers, building none of them.

    Returns
    -------
    dict of str to tuple of Gas
        Each model of :func:`molecular_rate_models` and its gases, in file order,
        read from ``molecular_ionisation_rates.toml`` alone; the ``default`` model's
        are the molecules optmat gives an ionisation rate for.

    Examples
    --------
    >>> [gas.value for gas in molecular_rate_coverage()["ppt_talebpour"]]
    ['D2', 'N2', 'O2']
    """
    return {name: tuple(_covered(name)) for name in _models_table()}


@functools.cache
def molecular_rate_models() -> dict[str, dict[Gas, MolecularRateModel]]:
    """Return every named molecular rate model, at 800 nm.

    Returns
    -------
    dict of str to dict of Gas to MolecularRateModel
        The models by name and molecule.

    Examples
    --------
    >>> sorted(molecular_rate_models())  # doctest: +NORMALIZE_WHITESPACE
    ['adk_hydrogenic', 'default', 'moadk_tong_2002', 'moadk_zhao_2010',
     'moadk_zhao_2017', 'ppt_talebpour', 'ppt_talebpour_ilkov']
    """
    return {
        name: {gas: _rate_model(name, gas) for gas in _covered(name)}
        for name in _models_table()
    }


def _status(entry: Mapping[str, object], where: str, default: Status) -> dict[str, str]:
    """Return the status fields a model's entry gives it.

    Parameters
    ----------
    entry : Mapping
        The model's or molecule's entry, with optional ``status`` and
        ``status_reason``.
    where : str
        The entry, for error messages.
    default : str
        The status where the entry declares none.

    Returns
    -------
    dict of str to str
        ``status`` and ``status_reason``.
    """
    status, reason = parse_status(entry, where, default=default)
    return {"status": status, "status_reason": reason}


def molecular_rate_model(
    gas: Gas | str, wavelength: float, model: str = "default"
) -> MolecularRateModel:
    """Return the ionisation rate of a molecule, or a named form of it.

    ``default`` is optmat's rate (TC-677 to TC-680); the others are the forms of the
    literature (TC-673).

    Parameters
    ----------
    gas : Gas or str
        The molecule.
    wavelength : float
        Vacuum wavelength of the field, m, > 0; ``math.inf`` for a static field
        (not for the ``ppt`` form).
    model : str, optional
        A name of :func:`molecular_rate_models`: ``default``, ``moadk_zhao_2010``,
        ``moadk_zhao_2017``, ``moadk_tong_2002``, ``ppt_talebpour``,
        ``ppt_talebpour_ilkov`` or ``adk_hydrogenic``.

    Returns
    -------
    MolecularRateModel
        The model.

    Raises
    ------
    DataError
        If the model is unknown or is not given for the molecule.

    Examples
    --------
    >>> model = molecular_rate_model("O2", 800e-9, "ppt_talebpour")
    >>> model.charge, model.form
    (0.53, 'ppt')
    >>> nitrogen = molecular_rate_model("N2", 800e-9)
    >>> nitrogen.form, nitrogen.electrons, nitrogen.ionization_energy_kind
    ('molecular', 2.0, 'vertical')
    """
    models = _models_table()
    if model not in models:
        msg = f"unknown molecular rate model {model!r}; the models are {sorted(models)}"
        raise DataError(msg)
    species = as_gas(gas)
    covered = _covered(model)
    if species not in covered:
        given = ", ".join(sorted(g.value for g in covered))
        msg = f"the model {model} is not given for {species.value}; it is for {given}"
        raise DataError(msg)
    return _rate_model(model, species).at_wavelength(wavelength)
