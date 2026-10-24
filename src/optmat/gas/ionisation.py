# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The parameters of the strong-field ionisation rates of the gases (TC-439).

Tunnelling and multiphoton ionisation rates (ADK, PPT, molecular ADK) are evaluated
from a handful of properties of the atom or molecule. This module supplies them,
curated, with their sources and uncertainties. The rates of the noble gases are
:mod:`optmat.gas.ionisation_rates` (TC-652); the molecules' stay with the propagation
code:

- the first **ionisation energy** :math:`I_p` (adiabatic or vertical, as chosen and
  stated per gas);
- the **ionising orbital**: its label, orbital angular momentum :math:`l` (an atom's
  valence shell) and magnetic quantum number :math:`m` along the molecular axis;
- its **structure coefficients** :math:`C_l`, the asymptotic expansion of the orbital
  in partial waves (Tong, Zhao & Lin 2002, eqs. 1-2), from each curated tabulation;
- the **static polarizabilities** of the neutral and of its cation, for the Stark shift
  of the ionisation energy and the dipole correction of the PPT rate. The neutral's is
  the static limit of the electronic part of the fitted dispersion model here (the
  infrared, nuclear part does not follow an optical field), with the fit's
  uncertainty;
- the **barrier-suppression constant** of Tong & Lin (2005) where tabulated.

From :math:`I_p` follow the effective principal quantum number
:math:`n^* = Z / \kappa`, :math:`\kappa = \sqrt{2 I_p}` (atomic units, :math:`Z` the
charge of the ion), and the asymptotic coefficient of the ADK and PPT rates in the
Hartree approximation,

.. math::

   |C_{n^* l^*}|^2
   = \frac{2^{2 n^*}}{n^* \, \Gamma(n^* + l^* + 1) \, \Gamma(n^* - l^*)},
   \qquad l^* = n^* - 1,

(Ammosov, Delone & Krainov 1986, eqs. 17 and 19, :math:`l^* = n_0^* - 1` with
:math:`n_0^*` the ground state's), and its Stirling form (their eq. 20, which is
misprinted; see :func:`adk_coefficient_squared`).

Every choice (which ionisation energy, which tabulation) is made in
``data/gas/ionisation.toml`` with its reason, and every value points to the curated
dataset it comes from.
"""

import functools
import math
import tomllib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace

import numpy as np

from optmat import _catalogue
from optmat.constants import AU_ENERGY, AU_POLARIZABILITY, E_CHARGE, HBAR, C
from optmat.datasets import Dataset
from optmat.exceptions import DataError
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.gases import Gas, as_gas
from optmat.gas.reduction import polarizability as dataset_polarizability

__all__ = [
    "IONISATION_FILE",
    "IonisationParameters",
    "StructureCoefficients",
    "adk_coefficient_squared",
    "combined_structure_coefficients",
    "electronic_static_polarizability",
    "hartree_coefficient_squared",
    "ionisation_parameters",
    "structure_coefficients",
]

IONISATION_FILE = GAS_DATA_DIRECTORY / "ionisation.toml"
"""The per-gas choices of the ionisation parameters, each with its reason."""

_ENERGY_UNITS: dict[str, float] = {
    "eV": E_CHARGE,
    "hartree": AU_ENERGY,
    # A wavenumber E / (h c) in cm^-1: h c in J cm.
    "cm-1": 2 * math.pi * HBAR * C * 100,
}


def hartree_coefficient_squared(n_star: float, l_star: float) -> float:
    """Return the asymptotic coefficient ``|C_{n* l*}|^2`` of the ADK and PPT rates.

    Ammosov, Delone & Krainov (1986), eq. (19): the Coulomb asymptotic wave function
    with the effective quantum numbers in place of ``n`` and ``l``.

    Parameters
    ----------
    n_star : float
        Effective principal quantum number ``Z / sqrt(2 I_p)``, > 0.
    l_star : float
        Effective orbital quantum number, ``n* - 1`` in ADK's choice.

    Returns
    -------
    float
        ``2^(2 n*) / (n* Gamma(n* + l* + 1) Gamma(n* - l*))``.

    Examples
    --------
    Hydrogen (``n* = 1``, ``l* = 0``) gives 4, the square of ``C_0 = 2``:

    >>> hartree_coefficient_squared(1.0, 0.0)
    4.0
    """
    return 2 ** (2 * n_star) / (
        n_star * math.gamma(n_star + l_star + 1) * math.gamma(n_star - l_star)
    )


def adk_coefficient_squared(n_star: float, l_star: float) -> float:
    """Return the Stirling form of ``|C_{n* l*}|^2`` (Ammosov's eq. 20, corrected).

    Stirling's formula for the two Gamma functions of eq. (19) gives
    ``(4 e^2 / (n*^2 - l*^2))^n* ((n* - l*) / (n* + l*))^(l* + 1/2) / (2 pi n*)``.
    Ammosov, Delone & Krainov print eq. (20) with ``4 e`` for ``4 e^2`` and the ratio
    inverted, ``(n* + l*) / (n* - l*)``: as printed it contradicts both their eq. (19)
    and, at ``l* = 0``, their eq. (11) (TC-441).

    Parameters
    ----------
    n_star : float
        Effective principal quantum number, > 0.
    l_star : float
        Effective orbital quantum number, ``|l*| < n*``.

    Returns
    -------
    float
        The corrected eq. (20).

    Examples
    --------
    At ``l* = 0`` it is their eq. (11) squared, ``(2 e / n*)^(2 n*) / (2 pi n*)``:

    >>> f"{adk_coefficient_squared(1.0, 0.0):.4f}"
    '4.7040'
    """
    return (
        (4 * math.e**2 / (n_star**2 - l_star**2)) ** n_star
        * ((n_star - l_star) / (n_star + l_star)) ** (l_star + 0.5)
        / (2 * math.pi * n_star)
    )


@dataclass(frozen=True)
class StructureCoefficients:
    """The asymptotic coefficients of an ionising orbital, from one tabulation.

    Attributes
    ----------
    source : str
        The curated dataset's key.
    model : str
        The theory the coefficients belong to (``moadk``, ``adk``, ``wfat``).
    orbital : str
        The orbital, as the source names it.
    m : int
        The magnetic quantum number along the molecular axis (an atom's ``m``); 0
        for a combination of several (a non-linear molecule).
    coefficients : dict of int to float
        ``C_l`` by partial wave ``l``, as printed (signs kept). For a combination of
        datasets (a non-linear molecule), the root of :attr:`partial_wave_squares`:
        what the average over orientations needs, not the coefficients of an aligned
        molecule.
    normalisation : str
        The source's definition of the coefficients.
    projections : dict of (int, int) to complex
        ``C_lm`` by ``(l, m)`` and orbital index, as printed: the real part the value,
        the imaginary part the dataset's ``imaginary`` column where it has one. Keys
        are ``(orbital, l, m)``, the orbitals numbered in the order given.
    """

    source: str
    model: str
    orbital: str
    m: int
    coefficients: Mapping[int, float]
    normalisation: str
    projections: Mapping[tuple[int, int, int], complex] = field(default_factory=dict)

    @property
    def partial_wave_squares(self) -> dict[int, float]:
        r"""The sum over ``m`` of :math:`|C_{lm}|^2` per ``l``, averaged over orbitals.

        By the orthogonality of the rotation matrices, the average over orientations
        of the structure factor depends on the coefficients only through these sums
        (:func:`optmat.gas.molecular_ionisation_rates.isotropic_structure_factors`).
        The orbitals of a degenerate shell are rotations of one another and have the
        same sums; they are averaged, and their number is part of the count of
        electrons.

        Returns
        -------
        dict of int to float
            :math:`\langle \sum_m |C_{lm}|^2 \rangle` by ``l``.
        """
        orbitals = {key[0] for key in self.projections}
        squares: dict[int, float] = {}
        for (_, l, _), value in self.projections.items():  # noqa: E741 - l is the partial wave
            squares[l] = squares.get(l, 0.0) + abs(value) ** 2 / len(orbitals)
        return dict(sorted(squares.items()))


def structure_coefficients(dataset: Dataset) -> StructureCoefficients:
    """Read one tabulation of structure coefficients.

    Parameters
    ----------
    dataset : Dataset
        A ``tunnelling_structure_coefficient`` dataset, one row per partial wave.

    Returns
    -------
    StructureCoefficients
        The coefficients by ``l``.

    Raises
    ------
    DataError
        If the dataset holds another quantity or a row is not an integer ``l``.
    """
    if dataset.quantity != "tunnelling_structure_coefficient":
        msg = f"{dataset.key} is a {dataset.quantity}, not structure coefficients"
        raise DataError(msg)
    try:
        partial_waves = [int(label) for label in dataset.abscissae]
    except ValueError as error:
        msg = f"{dataset.key}: every row label must be a partial wave l"
        raise DataError(msg) from error
    values = dataset.scaled_values()
    imaginary = (
        dataset.column_values("imaginary")
        if "imaginary" in dataset.columns
        else np.zeros(len(values))
    )
    parameters = dataset.parameters
    m = int(str(parameters["m"]))
    return StructureCoefficients(
        source=dataset.key,
        model=str(parameters["model"]),
        orbital=str(parameters["orbital"]),
        m=m,
        coefficients=dict(zip(partial_waves, map(float, values), strict=True)),
        normalisation=str(parameters["normalisation"]),
        projections={
            (0, l, m): complex(float(re), float(im))
            for l, re, im in zip(partial_waves, values, imaginary, strict=True)  # noqa: E741 - l is the partial wave
        },
    )


def combined_structure_coefficients(
    tabulation: Sequence[Dataset],
) -> StructureCoefficients:
    """Read the structure coefficients of an orbital given in several datasets.

    A non-linear molecule's orbital has coefficients at several ``m`` for each ``l``,
    and a degenerate shell several orbitals: one dataset per orbital and ``m``
    (the parameter ``orbital`` naming the orbital). The result keeps every
    coefficient in :attr:`StructureCoefficients.projections`; its ``coefficients``
    are the roots of :attr:`StructureCoefficients.partial_wave_squares`.

    Parameters
    ----------
    tabulation : sequence of Dataset
        ``tunnelling_structure_coefficient`` datasets of one source.

    Returns
    -------
    StructureCoefficients
        The combination.

    Raises
    ------
    DataError
        If a dataset holds another quantity, or two hold the same orbital and ``m``.
    """
    parts = [structure_coefficients(dataset) for dataset in tabulation]
    orbitals = list(dict.fromkeys(part.orbital for part in parts))
    projections: dict[tuple[int, int, int], complex] = {}
    for part in parts:
        index = orbitals.index(part.orbital)
        for (_, l, m), value in part.projections.items():  # noqa: E741 - l is the partial wave
            if (index, l, m) in projections:
                msg = f"{part.source}: (l, m) = ({l}, {m}) of {part.orbital} twice"
                raise DataError(msg)
            projections[index, l, m] = value
    combined = StructureCoefficients(
        source=", ".join(part.source for part in parts),
        model=parts[0].model,
        orbital=" / ".join(orbitals),
        m=0,
        coefficients={},
        normalisation=parts[0].normalisation,
        projections=projections,
    )
    return replace(
        combined,
        coefficients={
            l: math.sqrt(square)
            for l, square in combined.partial_wave_squares.items()  # noqa: E741 - l is the partial wave
        },
    )


def electronic_static_polarizability(gas: Gas | str) -> tuple[float, float, str]:
    """Return the static electronic polarizability of a gas from its fitted model.

    The fitted model's terms are the electronic oscillator strengths, then the fixed
    infrared (vibrational) bands and the magnetic term. The electronic static limit is
    the sum of the fitted terms' ``f_i / E_i^2`` alone: in an optical field the nuclei
    do not follow, so a Stark shift sees the electronic polarizability (TC-441).

    Parameters
    ----------
    gas : Gas or str
        The gas.

    Returns
    -------
    tuple of (float, float, str)
        The polarizability and its one-sigma uncertainty (Birge-scaled), C m^2/V, and
        the model's name.
    """
    # Loaded here: the fitted models bring the whole dispersion chain.
    from optmat.gas.models.fitted import FittedModel, fitted_model  # noqa: PLC0415

    model: FittedModel = fitted_model(f"fit/{as_gas(gas).value}")
    n_fitted = len(model.strengths) - model.n_fixed
    strengths = np.asarray(model.strengths[:n_fitted])
    energies = np.asarray(model.energies[:n_fitted])
    alpha_au = float(np.sum(strengths / energies**2))
    sigma_au = float(model.polarizability_sigma_au(np.zeros(1))[0])
    return alpha_au * AU_POLARIZABILITY, sigma_au * AU_POLARIZABILITY, model.name


@dataclass(frozen=True)
class IonisationParameters:
    """The curated parameters of a gas's strong-field ionisation rates.

    Attributes
    ----------
    gas : Gas
        The gas.
    ionization_energy : float
        First ionisation energy, J.
    ionization_energy_sigma : float
        Its one-sigma uncertainty, J; NaN where the source states none.
    ionization_energy_kind : str
        ``adiabatic`` or ``vertical`` (the two coincide for an atom).
    ionization_energy_source : str
        The dataset key.
    orbital : str
        The ionising orbital.
    orbital_l : int or None
        The orbital angular momentum of an atom's valence shell (the ``l`` of the ADK
        and PPT angular factor); None for a molecule.
    orbital_m : int
        The magnetic quantum number of the orbital (along a molecule's axis).
    structure : dict of str to StructureCoefficients
        The curated tabulations of the orbital's coefficients, by name.
    default_structure : str or None
        The name of the recommended tabulation, or None where none is published.
    polarizability : float
        Static electronic polarizability of the neutral, C m^2/V.
    polarizability_sigma : float
        Its one-sigma uncertainty, C m^2/V.
    polarizability_source : str
        The model it is the static limit of.
    vertical_ionization_energy : float or None
        The vertical ionisation energy of a molecule, J, where it differs and is
        curated: the energy at the neutral's geometry, which a tunnelling rate in the
        fixed-nuclei approximation sees (TC-444).
    vertical_ionization_energy_source : str
        The dataset key, or an empty string.
    ion_polarizability : float or None
        Static polarizability of the cation, C m^2/V, where curated.
    ion_polarizability_source : str
        The dataset key, or an empty string.
    barrier_suppression_alpha : float or None
        Tong & Lin's (2005) alpha, where tabulated.
    barrier_suppression_source : str
        The dataset key, or an empty string.
    barrier_suppression_zhang : tuple of float or None
        Zhang, Lan & Lu's (2014) a1, a2, a3: the rate is ``W_ADK exp(a1 x^2 + a2 x +
        a3)`` with ``x = F / F_b`` (the sign of their eq. 8 corrected, TC-542).
    barrier_suppression_zhang_source : str
        The dataset key, or an empty string.
    effective_charge : float or None
        Talebpour, Yang & Chin's (1999) effective charge of N2 or O2: a constant fitted
        within their PPT rate with l = m = 0, not a property of the molecule (for O2
        fitted with I_p = 12.55 eV, TC-542).
    effective_charge_source : str
        The dataset key, or an empty string.
    reasons : dict of str to str
        Why each choice was made (from ``ionisation.toml``).
    """

    gas: Gas
    ionization_energy: float
    ionization_energy_sigma: float
    ionization_energy_kind: str
    ionization_energy_source: str
    orbital: str
    orbital_l: int | None
    orbital_m: int
    structure: Mapping[str, StructureCoefficients]
    default_structure: str | None
    polarizability: float
    polarizability_sigma: float
    polarizability_source: str
    vertical_ionization_energy: float | None = None
    vertical_ionization_energy_source: str = ""
    ion_polarizability: float | None = None
    ion_polarizability_source: str = ""
    barrier_suppression_alpha: float | None = None
    barrier_suppression_source: str = ""
    barrier_suppression_zhang: tuple[float, float, float] | None = None
    barrier_suppression_zhang_source: str = ""
    effective_charge: float | None = None
    effective_charge_source: str = ""
    reasons: Mapping[str, str] = field(default_factory=dict)

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
    def kappa(self) -> float:
        """``kappa = sqrt(2 I_p)`` in atomic units.

        Returns
        -------
        float
            The decay constant of the orbital's asymptotic tail, 1/a0.
        """
        return math.sqrt(2 * self.ionization_energy_au)

    @property
    def n_star(self) -> float:
        """The effective principal quantum number ``Z / kappa`` of a singly charged ion.

        Returns
        -------
        float
            ``1 / sqrt(2 I_p)`` with ``I_p`` in hartree.
        """
        return 1.0 / self.kappa

    @property
    def l_star(self) -> float:
        """The effective orbital quantum number of ADK, ``n* - 1`` (Ammosov, eq. 17).

        Returns
        -------
        float
            ``n* - 1``.
        """
        return self.n_star - 1.0

    @property
    def hartree_coefficient_squared(self) -> float:
        """The PPT and ADK ``|C_{n* l*}|^2`` in the Hartree approximation.

        Returns
        -------
        float
            :func:`hartree_coefficient_squared` at ``n*`` and ``l* = n* - 1``.
        """
        return hartree_coefficient_squared(self.n_star, self.l_star)

    @property
    def polarizability_difference(self) -> float | None:
        """The neutral's polarizability less the cation's, C m^2/V.

        In a static field ``F`` each state's energy falls by ``alpha F^2 / 2``, so the
        ionisation energy rises by ``(alpha_0 - alpha_+) F^2 / 2``: the Stark shift of
        the PPT rate.

        Returns
        -------
        float or None
            ``alpha_0 - alpha_+``, or None where the cation's is not curated.
        """
        if self.ion_polarizability is None:
            return None
        return self.polarizability - self.ion_polarizability


def _dataset(by_key: Mapping[str, Dataset], key: str) -> Dataset:
    """Return a whole dataset named in ``ionisation.toml``.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    key : str
        ``source/dataset``.

    Returns
    -------
    Dataset
        The dataset.

    Raises
    ------
    DataError
        If it does not exist.
    """
    if key not in by_key:
        msg = f"ionisation.toml names an unknown dataset {key!r}"
        raise DataError(msg)
    return by_key[key]


def _row(by_key: Mapping[str, Dataset], reference: str) -> tuple[Dataset, int]:
    """Return a dataset and a row index from ``key`` or ``key#row``.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    reference : str
        ``source/dataset`` or ``source/dataset#row``.

    Returns
    -------
    tuple of (Dataset, int)
        The dataset and the index of the row (0 when it has one row).

    Raises
    ------
    DataError
        If the dataset or the row does not exist, or a row is needed and not given.
    """
    key, _, row = reference.partition("#")
    if key not in by_key:
        msg = f"ionisation.toml names an unknown dataset {key!r}"
        raise DataError(msg)
    dataset = by_key[key]
    if not row:
        if len(dataset) != 1:
            msg = f"{key} has {len(dataset)} rows: name one as {key}#<row>"
            raise DataError(msg)
        return dataset, 0
    if row not in dataset.rows:
        msg = f"{key} has no row {row!r}"
        raise DataError(msg)
    return dataset, dataset.rows.index(row)


def _energy(dataset: Dataset, index: int) -> tuple[float, float]:
    """Return an ionisation energy and its uncertainty in J.

    Parameters
    ----------
    dataset : Dataset
        An ``ionization_energy`` dataset.
    index : int
        The row.

    Returns
    -------
    tuple of (float, float)
        The energy and its one-sigma uncertainty (NaN when unstated), J.

    Raises
    ------
    DataError
        If the dataset is not an ionisation energy.
    """
    if dataset.quantity != "ionization_energy":
        msg = f"{dataset.key} is a {dataset.quantity}, not an ionization_energy"
        raise DataError(msg)
    factor = _ENERGY_UNITS[dataset.unit]
    value = float(dataset.scaled_values()[index]) * factor
    sigma = float(dataset.standard_uncertainty()[index]) * factor
    return value, sigma


def _strings(entry: Mapping[str, object], name: str) -> dict[str, str]:
    """Return a sub-table of an entry as strings.

    Parameters
    ----------
    entry : mapping
        A gas's table in ``ionisation.toml``.
    name : str
        The sub-table.

    Returns
    -------
    dict of str to str
        Its keys and values, or an empty dict when absent.

    Raises
    ------
    DataError
        If it is not a table.
    """
    table = entry.get(name, {})
    if not isinstance(table, Mapping):
        msg = f"ionisation.toml: {name} must be a table"
        raise DataError(msg)
    return {str(k): str(v) for k, v in table.items()}  # pyright: ignore[reportUnknownVariableType]


def _structures(
    entry: Mapping[str, object], by_key: Mapping[str, Dataset]
) -> dict[str, StructureCoefficients]:
    """Read a gas's tabulations of structure coefficients.

    Parameters
    ----------
    entry : mapping
        Its table in ``ionisation.toml``.
    by_key : mapping of str to Dataset
        The curated datasets.

    Returns
    -------
    dict of str to StructureCoefficients
        The tabulations by name: one dataset, or a list of them combined
        (:func:`combined_structure_coefficients`).

    Raises
    ------
    DataError
        If ``structure`` is not a table, or a named dataset does not exist.
    """
    tables = entry.get("structure", {})
    if not isinstance(tables, Mapping):
        msg = "ionisation.toml: structure must be a table"
        raise DataError(msg)
    structure: dict[str, StructureCoefficients] = {}
    for name, key in tables.items():  # pyright: ignore[reportUnknownVariableType]
        if isinstance(key, list):
            structure[str(name)] = combined_structure_coefficients(
                [_dataset(by_key, str(k)) for k in key]  # pyright: ignore[reportUnknownVariableType, reportUnknownArgumentType]
            )
        else:
            structure[str(name)] = structure_coefficients(_dataset(by_key, str(key)))
    return structure


def _parameters(
    gas: Gas, entry: Mapping[str, object], by_key: Mapping[str, Dataset]
) -> IonisationParameters:
    """Build one gas's parameters from its entry in ``ionisation.toml``.

    Parameters
    ----------
    gas : Gas
        The gas.
    entry : mapping
        Its table in ``ionisation.toml``.
    by_key : mapping of str to Dataset
        The curated datasets.

    Returns
    -------
    IonisationParameters
        The parameters.

    Raises
    ------
    DataError
        If a named dataset does not exist, holds another material or quantity, or the
        default tabulation is not among the tabulations.
    """
    reasons = _strings(entry, "reasons")
    energy_dataset, energy_row = _row(by_key, str(entry["ionization_energy"]))
    if energy_dataset.material != gas.value:
        msg = f"{energy_dataset.key} is {energy_dataset.material}'s, not {gas.value}'s"
        raise DataError(msg)
    energy, energy_sigma = _energy(energy_dataset, energy_row)
    vertical = None
    vertical_source = ""
    if "vertical_ionization_energy" in entry:
        vertical_source = str(entry["vertical_ionization_energy"])
        vertical_dataset, vertical_row = _row(by_key, vertical_source)
        if (
            vertical_dataset.material != gas.value
            or vertical_dataset.parameters.get("kind") != "vertical"
        ):
            msg = f"{vertical_dataset.key} is not a vertical energy of {gas.value}"
            raise DataError(msg)
        vertical = _energy(vertical_dataset, vertical_row)[0]
    structure = _structures(entry, by_key)
    default = entry.get("default_structure")
    if default is not None and default not in structure:
        msg = f"{gas.value}: default_structure {default!r} is not a tabulation"
        raise DataError(msg)
    alpha, alpha_sigma, alpha_source = electronic_static_polarizability(gas)
    ion_alpha = None
    ion_source = ""
    if "ion_polarizability" in entry:
        ion_source = str(entry["ion_polarizability"])
        ion_dataset, ion_row = _row(by_key, ion_source)
        if ion_dataset.material != f"{gas.value}+":
            msg = f"{ion_dataset.key} is not {gas.value}+'s polarizability"
            raise DataError(msg)
        ion_alpha = float(dataset_polarizability(ion_dataset)[ion_row])
    bsi = None
    bsi_source = ""
    if "barrier_suppression_alpha" in entry:
        bsi_source = str(entry["barrier_suppression_alpha"])
        bsi_dataset, bsi_row = _row(by_key, bsi_source)
        bsi = float(bsi_dataset.scaled_values()[bsi_row])
    zhang = None
    zhang_source = ""
    if "barrier_suppression_zhang" in entry:
        zhang_source = str(entry["barrier_suppression_zhang"])
        values = _dataset(by_key, zhang_source).scaled_values()
        zhang = (float(values[0]), float(values[1]), float(values[2]))
    charge = None
    charge_source = ""
    if "effective_charge" in entry:
        charge_source = str(entry["effective_charge"])
        charge_dataset, charge_row = _row(by_key, charge_source)
        charge = float(charge_dataset.scaled_values()[charge_row])
    orbital_l = entry.get("orbital_l")
    return IonisationParameters(
        gas=gas,
        ionization_energy=energy,
        ionization_energy_sigma=energy_sigma,
        ionization_energy_kind=str(energy_dataset.parameters["kind"]),
        ionization_energy_source=str(entry["ionization_energy"]),
        orbital=str(entry["orbital"]),
        orbital_l=None if orbital_l is None else int(str(orbital_l)),
        orbital_m=int(str(entry.get("orbital_m", 0))),
        structure=structure,
        default_structure=None if default is None else str(default),
        polarizability=alpha,
        polarizability_sigma=alpha_sigma,
        polarizability_source=alpha_source,
        vertical_ionization_energy=vertical,
        vertical_ionization_energy_source=vertical_source,
        ion_polarizability=ion_alpha,
        ion_polarizability_source=ion_source,
        barrier_suppression_alpha=bsi,
        barrier_suppression_source=bsi_source,
        barrier_suppression_zhang=zhang,
        barrier_suppression_zhang_source=zhang_source,
        effective_charge=charge,
        effective_charge_source=charge_source,
        reasons=reasons,
    )


@functools.cache
def _entries() -> dict[Gas, Mapping[str, object]]:
    """Read ``ionisation.toml``.

    Returns
    -------
    dict of Gas to mapping
        Each gas's entry.
    """
    with IONISATION_FILE.open("rb") as handle:
        table = tomllib.load(handle)
    return {as_gas(name): entry for name, entry in table.items()}


@functools.cache
def _gas_parameters(species: Gas) -> IonisationParameters:
    """Build one gas's parameters, reading only the datasets it cites.

    Parameters
    ----------
    species : Gas
        A gas of ``ionisation.toml``.

    Returns
    -------
    IonisationParameters
        Its parameters.
    """
    return _parameters(species, _entries()[species], _catalogue.by_key())


def ionisation_parameters(gas: Gas | str) -> IonisationParameters:
    """Return the curated ionisation parameters of a gas.

    Parameters
    ----------
    gas : Gas or str
        The gas.

    Returns
    -------
    IonisationParameters
        Its parameters, each with its source.

    Raises
    ------
    DataError
        If the gas has no ionisation parameters.
    """
    species = as_gas(gas)
    if species not in _entries():
        msg = f"no ionisation parameters for {species.value}"
        raise DataError(msg)
    return _gas_parameters(species)
