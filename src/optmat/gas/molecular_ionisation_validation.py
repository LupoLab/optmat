# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The molecular rates against their references.

The references that test a rate of an **unaligned** molecule are rates resolved in the
angle between the field and the molecular axis, averaged over that angle. Majety and
Scrinzi (2015) give all-electron static rates of N2 at seven angles, of CO at thirteen
and of H2 at two; Jagau (2018) gives coupled-cluster rates of O2 at nine angles, of CO
at thirteen, of H2 at two and of CH4 in three orientations.
:func:`orientation_averaged_rates` averages them over the orientations and compares the
models of :mod:`optmat.gas.molecular_ionisation_rates` with the averages. The
references are for fixed nuclei at the equilibrium distance, and the models are
compared at fixed nuclei.

:func:`vibrational_factors` compares the vibrational factor of H2 with Saenz's (2002)
rates averaged over the vibrational ground state. :func:`fitted_constants` recomputes
the constants that the default model's correction and vibrational factor were fitted
to, so that a test can hold the curated values to their data.
"""

import math
import re
import warnings
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import Final

import numpy as np
from jaxtyping import Float
from scipy import interpolate

from optmat.constants import AU_ANGULAR_FREQUENCY, AU_FIELD
from optmat.datasets import Dataset, datasets
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas.gases import Gas, as_gas
from optmat.gas.molecular_ionisation_rates import molecular_rate_model

__all__ = [
    "ANGLE_RESOLVED",
    "FITTED_REFERENCES",
    "JAGAU_2018",
    "FittedConstants",
    "OrientationAverage",
    "VibrationalFactors",
    "fitted_constants",
    "orientation_average",
    "orientation_averaged_rates",
    "vibrational_factors",
]

ANGLE_RESOLVED: Final[Mapping[str, tuple[str, str]]] = {
    "N2": ("majety_2015/tableA1_N2_{angle}deg", "Majety & Scrinzi 2015, haCC(5)"),
    "H2": ("majety_2015/table{angle}_hacc5_H2", "Majety & Scrinzi 2015, haCC(5)"),
    "CO": ("majety_2015/tableA2_CO_{angle}deg", "Majety & Scrinzi 2015, haCC(6)"),
}
"""The angle-resolved static rates of each molecule: the datasets' keys and a label.

For N2 the angles are 0 to 90 degrees in steps of 15, for CO, which is polar, 0 to 180
(0 degrees: the field from the O atom to the C atom). For H2 only the parallel and the
perpendicular rates exist (``2_parallel`` and ``3_perpendicular``).
"""

JAGAU_2018: Final[Mapping[str, tuple[str, str]]] = {
    "H2": (r"jagau_2018/tableS3_fci_H2_(\d+)deg", "Jagau 2018, full CI"),
    "CO": (r"jagau_2018/tableS6_ccsd_CO_(\d+)deg", "Jagau 2018, CCSD"),
    "O2": (r"jagau_2018/tableS(?:9|10)_ccsd_O2_(\d+)deg", "Jagau 2018, CCSD"),
    "CH4": (r"jagau_2018/tableS8_ccsd_CH4_([ABC])", "Jagau 2018, CCSD"),
}
"""Jagau's (2018) correlated static rates of each molecule: a pattern of the keys, and a
label.

O2's angles are 0 to 90 degrees, at every field for 0 and 90 (Table S9) and at 0.06 and
0.10 au for those between (Table S10): the average exists at 0.06 and 0.10 au. CO's are
0 to 180 degrees with 0 degrees the field from C to O, the reverse of Majety and
Scrinzi's. CH4's orientations are A (along a CH bond, towards H), B (along a CH bond,
towards C) and C (bisecting two CH bonds), averaged with the weights 0.3, 0.3 and 0.4
of Lebedev's 14-point rule.
"""

REFERENCES: Final[Mapping[str, Mapping[str, tuple[str, str]]]] = {
    "majety_2015": ANGLE_RESOLVED,
    "jagau_2018": JAGAU_2018,
}

FITTED_REFERENCES: Final[Mapping[str, str]] = {
    "H2": "majety_2015",
    "N2": "majety_2015",
    "CO": "majety_2015",
    "O2": "jagau_2018",
    "CH4": "jagau_2018",
}
"""The reference each molecule's alpha is fitted to (or, for CO, tested against)."""

_H2_TABLES: Final = {0.0: "2_parallel", 90.0: "3_perpendicular"}
_QUARTER_TURN: Final = 90.0
_HALF_TURN: Final = 180.0
_TETRAHEDRAL: Final = {
    "A": (0.0, 0.3),
    "C": (math.degrees(math.acos(1 / math.sqrt(3))), 0.4),
    "B": (180.0, 0.3),
}
"""CH4's orientations: the angle between the field and one CH bond, and the weight.

A points along a CH bond towards H, B along a CH bond towards C, and C bisects two CH
bonds. With the hydrogens at four corners of a cube, the eight corner directions are
the four A and the four B, and the six face normals the C. Lebedev's 14-point rule on
the sphere, weight 1/15 on each face normal and 3/40 on each corner, is exact for
spherical harmonics to degree 5: the average over orientations is 0.3 A + 0.3 B +
0.4 C.
"""


def orientation_average(
    angles: Float[np.ndarray, " n_angle"], rates: Float[np.ndarray, " n_angle *shape"]
) -> Float[np.ndarray, " *shape"]:
    r"""Return the average of a rate over the orientations of a linear molecule.

    For a molecule without a permanent dipole the rate is even about 90 degrees, and
    the average over the sphere is :math:`\int_0^{\pi/2} w(\theta) \sin\theta\,
    d\theta`. The rate is a smooth function of :math:`\cos^2\theta`. For two angles
    it is taken as :math:`w_\perp + (w_\parallel - w_\perp) \cos^2\theta`, whose
    average is :math:`(w_\parallel + 2 w_\perp)/3`. For more, the logarithm of the
    rate is interpolated in :math:`\cos^2\theta` by a shape-preserving cubic
    (which does not oscillate between the angles as a polynomial through all of them
    can) and the rate integrated by Gauss-Legendre quadrature. For a polar molecule,
    whose angles run to 180 degrees, the average is :math:`\frac12 \int_0^\pi
    w(\theta) \sin\theta\, d\theta`, and the logarithm is interpolated in
    :math:`\cos\theta`.

    Parameters
    ----------
    angles : numpy.ndarray
        Angles between the field and the molecular axis, degrees, in [0, 90],
        including 0 and 90, or in [0, 180], including 0 and 180.
    rates : numpy.ndarray
        The rates at those angles, first axis the angle, > 0.

    Returns
    -------
    numpy.ndarray
        The average, with the shape of the remaining axes.

    Raises
    ------
    ValueError
        If the angles do not span 0 to 90 or 0 to 180 degrees.

    Examples
    --------
    The average of :math:`1 + \cos^2\theta` over the sphere is 4/3:

    >>> import numpy as np
    >>> angles = np.arange(0.0, 91.0, 15.0)
    >>> rates = 1 + np.cos(np.radians(angles)) ** 2
    >>> f"{float(orientation_average(angles, rates)):.3f}"
    '1.333'
    """
    angles = np.asarray(angles, dtype=float)
    rates = np.asarray(rates, dtype=float)
    order = np.argsort(angles)
    angles, rates = angles[order], rates[order]
    polar = angles[-1] == _HALF_TURN
    if angles[0] != 0.0 or angles[-1] not in (_QUARTER_TURN, _HALF_TURN):
        msg = "the angles must span 0 to 90 or 0 to 180 degrees"
        raise ValueError(msg)
    nodes, weights = np.polynomial.legendre.leggauss(48)
    flat = rates.reshape(len(angles), -1)
    if polar:
        # The average over the sphere is half the integral over u = cos(theta) from
        # -1 to 1, and the rate a smooth function of u.
        x = np.cos(np.radians(angles))
        u = nodes
        weights = 0.5 * weights
    else:
        # The average over the sphere is the integral over u = cos(theta) from 0 to
        # 1, and the rate a smooth function of u^2.
        x = np.cos(np.radians(angles)) ** 2
        u = (0.5 * (nodes + 1)) ** 2
        weights = 0.5 * weights
        two = 2
        if len(angles) == two:
            return np.asarray((flat[0] + 2 * flat[1]) / 3).reshape(rates.shape[1:])
    average = np.empty(flat.shape[1])
    ascending = np.argsort(x)
    for k in range(flat.shape[1]):
        logarithm = interpolate.PchipInterpolator(
            x[ascending], np.log(flat[ascending, k])
        )
        average[k] = float(np.sum(weights * np.exp(logarithm(u))))
    return average.reshape(rates.shape[1:])


@dataclass(frozen=True)
class OrientationAverage:
    """A molecular rate against the orientation average of angle-resolved rates.

    Attributes
    ----------
    gas : Gas
        The molecule.
    model : str
        The molecular rate model.
    label : str
        The reference.
    angles : numpy.ndarray
        The angles of the reference, degrees.
    field : numpy.ndarray
        The static fields, V/m.
    resolved : numpy.ndarray
        The reference's rates, atomic units, first axis the angle.
    expected : numpy.ndarray
        Their average over orientations, atomic units.
    predicted : numpy.ndarray
        The model's static rate, atomic units.
    """

    gas: Gas
    model: str
    label: str
    angles: Float[np.ndarray, " n_angle"]
    field: Float[np.ndarray, " n"]
    resolved: Float[np.ndarray, " n_angle n"]
    expected: Float[np.ndarray, " n"]
    predicted: Float[np.ndarray, " n"]

    @property
    def ratio(self) -> Float[np.ndarray, " n"]:
        """The model over the reference.

        Returns
        -------
        numpy.ndarray
            ``predicted / expected``.
        """
        return self.predicted / self.expected

    @property
    def anisotropy(self) -> Float[np.ndarray, " n"]:
        """The reference's rate along the axis over its average.

        Returns
        -------
        numpy.ndarray
            ``resolved[0] / expected``.
        """
        return self.resolved[0] / self.expected


def _usable(dataset: Dataset) -> tuple[list[str], Float[np.ndarray, " n"]]:
    """Return the fields (as printed) and rates of a dataset without excluded rows.

    Parameters
    ----------
    dataset : Dataset
        A static ionisation rate against the field in atomic units.

    Returns
    -------
    tuple
        The printed fields and the rates, atomic units.
    """
    excluded = set(dataset.excluded_rows)
    rows = [str(row) for row in dataset.rows]
    keep = [k for k, row in enumerate(rows) if row not in excluded]
    values = np.asarray(dataset.scaled_values(), dtype=float)
    return [rows[k] for k in keep], values[keep]


def _found(reference: str, gas: str, by_key: Mapping[str, Dataset]) -> dict[float, str]:
    """Return the keys of a reference's angle-resolved datasets of a gas, by angle.

    Parameters
    ----------
    reference : str
        A key of the references (``majety_2015``, ``jagau_2018``).
    gas : str
        The molecule.
    by_key : mapping of str to Dataset
        The curated datasets.

    Returns
    -------
    dict of float to str
        Dataset keys by the angle between the field and the axis, degrees (for CH4,
        between the field and one CH bond).
    """
    pattern, _ = REFERENCES[reference][gas]
    if reference == "majety_2015" and gas == "H2":
        return {angle: pattern.format(angle=name) for angle, name in _H2_TABLES.items()}
    if reference == "majety_2015":
        pattern = re.escape(pattern).replace(r"\{angle\}", r"(\d+)")
    matcher = re.compile(pattern + "$")
    found = {}
    for key in by_key:
        match = matcher.match(key)
        if match is not None:
            label = match.group(1)
            angle = _TETRAHEDRAL[label][0] if label in _TETRAHEDRAL else float(label)
            found[angle] = key
    return found


def orientation_averaged_rates(
    gas: Gas | str, model: str = "default", reference: str = "majety_2015"
) -> OrientationAverage:
    """Compare a molecular rate with the orientation average of resolved rates.

    Parameters
    ----------
    gas : Gas or str
        N2, CO or H2 (``majety_2015``); O2, CO, H2 or CH4 (``jagau_2018``).
    model : str, optional
        A model of :func:`optmat.gas.molecular_ionisation_rates.molecular_rate_model`
        with a static limit (not the ``ppt`` forms). It is evaluated at fixed nuclei:
        without its vibrational factor.
    reference : str, optional
        ``majety_2015``, Majety and Scrinzi's haCC rates, or ``jagau_2018``, Jagau's
        coupled-cluster rates.

    Returns
    -------
    OrientationAverage
        The comparison, at the fields the reference gives at every angle.

    Raises
    ------
    DataError
        If no angle-resolved rate of the reference is curated for the molecule.
    """
    species = as_gas(gas)
    if reference not in REFERENCES or species.value not in REFERENCES[reference]:
        msg = (
            f"no angle-resolved static rates of {reference} are curated for "
            f"{species.value}"
        )
        raise DataError(msg)
    label = REFERENCES[reference][species.value][1]
    by_key = {d.key: d for d in datasets()}
    found = _found(reference, species.value, by_key)
    angles = np.array(sorted(found))
    tables = [_usable(by_key[found[angle]]) for angle in angles]
    common = [row for row in tables[0][0] if all(row in rows for rows, _ in tables)]
    resolved = np.array(
        [[values[rows.index(row)] for row in common] for rows, values in tables]
    )
    field = np.array([float(row) for row in common]) * AU_FIELD
    # The references are for fixed nuclei: the vibrational factor is left out.
    rate = replace(
        molecular_rate_model(species, math.inf, model), vibrational_constant=0.0
    )
    with warnings.catch_warnings():
        # The references run to and beyond the barrier-suppression field, where the
        # comparison shows how the uncorrected forms fail.
        warnings.simplefilter("ignore", ModelRangeWarning)
        predicted = rate.instantaneous_rate(field) / AU_ANGULAR_FREQUENCY
    if species.value == "CH4":
        weights = np.array([weight for _, weight in sorted(_TETRAHEDRAL.values())])
        expected = weights @ resolved
    else:
        expected = orientation_average(angles, resolved)
    return OrientationAverage(
        species, model, label, angles, field, resolved, expected, predicted
    )


@dataclass(frozen=True)
class VibrationalFactors:
    """The vibrational factor of a model against Saenz's averaged rates.

    Attributes
    ----------
    field : numpy.ndarray
        The static fields, V/m.
    expected : numpy.ndarray
        Saenz's rate averaged over the vibrational ground state (v = 0) over his rate
        at 1.4 bohr, the field along the axis.
    predicted : numpy.ndarray
        The model's factor :math:`e^{c/F^2}`.
    """

    field: Float[np.ndarray, " n"]
    expected: Float[np.ndarray, " n"]
    predicted: Float[np.ndarray, " n"]

    @property
    def ratio(self) -> Float[np.ndarray, " n"]:
        """The model over the reference.

        Returns
        -------
        numpy.ndarray
            ``predicted / expected``.
        """
        return self.predicted / self.expected


def _saenz_ratio() -> tuple[Float[np.ndarray, " n"], Float[np.ndarray, " n"]]:
    """Return H2's vibrationally averaged rate over the fixed-nuclei rate.

    Returns
    -------
    tuple of numpy.ndarray
        The fields, atomic units, and the ratio: saenz_2002 Table I row (a), v = 0,
        over Saenz's fixed-nuclei rates at 1.4 bohr as majety_2015 reprints them.
    """
    by_key = {d.key: d for d in datasets()}
    averaged = by_key["saenz_2002/tableI_a_v0"]
    fixed = by_key["majety_2015/table2_parallel_saenz_2000_H2"]
    fixed_rows, fixed_values = _usable(fixed)
    averaged_rows, averaged_values = _usable(averaged)
    fields = []
    ratios = []
    for row, value in zip(fixed_rows, fixed_values, strict=True):
        matches = [
            k for k, other in enumerate(averaged_rows) if float(other) == float(row)
        ]
        if matches:
            fields.append(float(row))
            ratios.append(averaged_values[matches[0]] / value)
    return np.array(fields), np.array(ratios)


def vibrational_factors(model: str = "default") -> VibrationalFactors:
    """Compare H2's vibrational factor with Saenz's averaged rates.

    Parameters
    ----------
    model : str, optional
        A model with a vibrational factor for H2.

    Returns
    -------
    VibrationalFactors
        The factor of the model and of the reference, at the fields of both
        (0.06 to 0.12 au).

    Examples
    --------
    >>> comparison = vibrational_factors()
    >>> [f"{x:.2f}" for x in comparison.expected]
    ['1.56', '1.32', '1.16', '1.16']
    """
    field, expected = _saenz_ratio()
    constant = molecular_rate_model("H2", math.inf, model).vibrational_constant
    return VibrationalFactors(field * AU_FIELD, expected, np.exp(constant / field**2))


@dataclass(frozen=True)
class FittedConstants:
    """The constants of the default model, as fitted to their references.

    Attributes
    ----------
    alpha : dict of str to float
        Tong and Lin's alpha of each molecule of :data:`FITTED_REFERENCES` by least
        squares in the logarithm, with the model's electron count, against the
        orientation-averaged all-electron static rates. H2's, N2's, O2's and CH4's
        are the curated ones; CO takes N2's, and its own is what the test of that
        transfer would ask for.
    free_electrons : dict of str to float
        The electron count when it is fitted too.
    free_alpha : dict of str to float
        Alpha when the count is fitted too.
    vibrational_constant : float
        :math:`c` of H2's vibrational factor, by least squares in its logarithm.
    """

    alpha: Mapping[str, float]
    free_electrons: Mapping[str, float]
    free_alpha: Mapping[str, float]
    vibrational_constant: float


def fitted_constants(model: str = "default") -> FittedConstants:
    r"""Recompute the default model's fitted constants from their references.

    With :math:`y = \ln(W_\text{reference}/W_0)`, :math:`W_0` the model at fixed
    nuclei without its correction and its count, and :math:`x = F/(I_p\kappa^3)`:
    :math:`y = \ln N - \alpha x`, :math:`N` fixed at the model's count, gives
    :math:`\alpha`; with :math:`N` free, both. :math:`\ln V = c/F^2` against Saenz's
    ratios gives :math:`c`.

    Parameters
    ----------
    model : str, optional
        The model whose electron count, energy and coefficients are used.

    Returns
    -------
    FittedConstants
        The constants.
    """
    alpha: dict[str, float] = {}
    free_electrons: dict[str, float] = {}
    free_alpha: dict[str, float] = {}
    for gas, reference in FITTED_REFERENCES.items():
        comparison = orientation_averaged_rates(gas, model, reference)
        rate = molecular_rate_model(gas, math.inf, model)
        bare = replace(
            rate,
            static_correction="none",
            electrons=1.0,
            vibrational_constant=0.0,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ModelRangeWarning)
            predicted = bare.instantaneous_rate(comparison.field) / AU_ANGULAR_FREQUENCY
        y = np.log(comparison.expected / predicted)
        energy = rate.ionization_energy_au
        x = comparison.field / AU_FIELD / (energy * (2 * energy) ** 1.5)
        count = math.log(rate.electrons)
        alpha[gas] = float(np.sum((count - y) * x) / np.sum(x * x))
        design = np.stack([np.ones_like(x), -x], axis=1)
        (log_count, slope), *_ = np.linalg.lstsq(design, y, rcond=None)
        free_electrons[gas] = math.exp(float(log_count))
        free_alpha[gas] = float(slope)
    field, ratio = _saenz_ratio()
    constant = float(np.sum(np.log(ratio) / field**2) / np.sum(1 / field**4))
    return FittedConstants(alpha, free_electrons, free_alpha, constant)
