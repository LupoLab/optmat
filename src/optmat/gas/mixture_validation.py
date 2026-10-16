# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The multi-fluid mixture models against the measured data (OEP 4, section 8; TC-778).

Every curated measured virial coefficient of a pair (cross, excess or of a mixture) and
every mixture compressibility factor or density between 250 and 350 K, up to the
mixture's optical-domain density, is compared with CoolProp's multi-fluid model of the
pair. A point passes within twice the combined uncertainty; a pair qualifies when at
least 95 % of its counted points pass in each kind, with at least three counted points.
The rule was fixed before the comparison (TC-778); ``data/gas/mixture_pairs.toml``
records its result per pair, and a test keeps the two in step.
"""

import functools
import math
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Final

from optmat.conditions import State
from optmat.constants import K_B, N_A
from optmat.datasets import Dataset, datasets
from optmat.exceptions import ModelRangeWarning
from optmat.gas.mixture import mixture_optics
from optmat.gas.mixture_eos import (
    TRUNCATION_BOUND,
    MixtureEquationOfState,
    _binary,  # pyright: ignore[reportPrivateUsage]
    _pair_records,  # pyright: ignore[reportPrivateUsage]
    mixture_equation_of_state,
)

__all__ = [
    "MODEL_DENSITY_UNCERTAINTY",
    "TEMPERATURES",
    "Comparison",
    "PairResult",
    "comparisons",
    "effective_cross_virial",
    "fallback_comparisons",
    "fallback_results",
    "pair_results",
]

TEMPERATURES: Final = (250.0, 350.0)
"""The temperatures of the comparison, K (OEP 4, section 7.4)."""

MODEL_DENSITY_UNCERTAINTY: Final = {
    ("N2", "O2"): 0.001,
    ("Ar", "N2"): 0.005,
    ("Ar", "O2"): 0.010,
}
"""The multi-fluid models' own stated relative density uncertainty, read as standard
(Gernert's gas-phase statements, gernert_2013); Tkaczuk et al. state none."""

_VIRIALS: Final = (
    "cross_second_virial_coefficient",
    "excess_second_virial_coefficient",
    "mixture_second_virial_coefficient",
)
_DENSITIES: Final = ("mixture_compressibility_factor", "mixture_density")


@dataclass(frozen=True)
class Comparison:
    """One measured point against the multi-fluid model.

    Attributes
    ----------
    dataset : str
        The dataset, ``<source>/<id>``.
    pair : tuple of str
        The pair, sorted.
    kind : str
        ``virial`` or ``density``.
    temperature : float
        K.
    abscissa : float
        The temperature (a virial), or the pressure (Pa) or density (m^-3) of a point.
    measured, model : float
        The values, in the dataset's unit.
    sigma : float
        The combined standard uncertainty; NaN where the measurement states none.
    counted : bool
        Whether the point counts for the pair's qualification.
    passed : bool
        Whether ``|model - measured| <= 2 sigma``.
    """

    dataset: str
    pair: tuple[str, str]
    kind: str
    temperature: float
    abscissa: float
    measured: float
    model: float
    sigma: float
    counted: bool
    passed: bool


@dataclass(frozen=True)
class PairResult:
    """A pair's comparison under TC-778's rule.

    Attributes
    ----------
    pair : tuple of str
        The pair, sorted.
    counted : dict
        Per kind, the number of counted points.
    passed : dict
        Per kind, the number that pass.
    qualified : bool
        At least 95 % pass in each kind with counted points, and at least three
        counted in all.
    """

    pair: tuple[str, str]
    counted: Mapping[str, int]
    passed: Mapping[str, int]
    qualified: bool


def _quiet_eos(
    composition: Mapping[str, float], method: str = "multi-fluid"
) -> MixtureEquationOfState:
    """Return a composition's equation by one method, its warning silenced.

    Parameters
    ----------
    composition : Mapping
        Species to mole fraction.
    method : str, optional
        The density method; the multi-fluid model by default.

    Returns
    -------
    MixtureEquationOfState
        The equation, qualified or not.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        return mixture_equation_of_state(composition, method=method)


@functools.cache
def effective_cross_virial(
    species: str, partner: str, temperature: float, fraction: float = 0.5
) -> float:
    """Return a multi-fluid model's effective cross second virial coefficient.

    Parameters
    ----------
    species, partner : str
        The pair.
    temperature : float
        K.
    fraction : float, optional
        The mole fraction of ``species``; 1/2 by default.

    Returns
    -------
    float
        ``[B_mix(x) - x^2 B_11 - (1 - x)^2 B_22] / [2 x (1 - x)]``, m^3/mol.

    Examples
    --------
    >>> f"{effective_cross_virial('Ar', 'N2', 300.0) * 1e6:.2f}"
    '-10.03'
    """
    mixture = _quiet_eos({species: fraction, partner: 1 - fraction})
    pure = [_quiet_eos({g: 1.0}).second_virial(temperature) for g in (species, partner)]
    b_mix = mixture.second_virial(temperature)
    like = fraction**2 * pure[0] + (1 - fraction) ** 2 * pure[1]
    return (b_mix - like) / (2 * fraction * (1 - fraction))


def _pair(dataset: Dataset) -> tuple[str, str]:
    """Return a dataset's pair, sorted.

    Parameters
    ----------
    dataset : Dataset
        A pair's dataset.

    Returns
    -------
    tuple of str
        The material and its partner, sorted.
    """
    first, second = sorted((str(dataset.material), str(dataset.parameters["partner"])))
    return first, second


def _virial_points(dataset: Dataset) -> list[Comparison]:
    """Compare a virial dataset's points in :data:`TEMPERATURES`.

    Parameters
    ----------
    dataset : Dataset
        A cross, excess or mixture second virial coefficient.

    Returns
    -------
    list of Comparison
        One per point in range.
    """
    material, partner = str(dataset.material), str(dataset.parameters["partner"])
    pair = _pair(dataset)
    fraction = float(str(dataset.parameters.get("mole_fraction", "0.5")))
    sigmas = dataset.standard_uncertainty()
    out = []
    for temperature, measured, sigma in zip(
        dataset.abscissa_values(), dataset.scaled_values(), sigmas, strict=True
    ):
        if not TEMPERATURES[0] <= temperature <= TEMPERATURES[1]:
            continue
        cross = effective_cross_virial(material, partner, float(temperature), fraction)
        pure = [
            _quiet_eos({g: 1.0}).second_virial(float(temperature))
            for g in (material, partner)
        ]
        if dataset.quantity == "cross_second_virial_coefficient":
            model = cross
        elif dataset.quantity == "excess_second_virial_coefficient":
            model = cross - (pure[0] + pure[1]) / 2
        else:
            model = _quiet_eos(
                {material: fraction, partner: 1 - fraction}
            ).second_virial(float(temperature))
        model *= 1e6  # m^3/mol -> cm^3/mol
        counted = dataset.use == "primary" and not math.isnan(sigma)
        out.append(
            Comparison(
                dataset=dataset.key,
                pair=pair,
                kind="virial",
                temperature=float(temperature),
                abscissa=float(temperature),
                measured=float(measured),
                model=model,
                sigma=float(sigma),
                counted=counted,
                passed=bool(counted and abs(model - measured) <= 2 * sigma),
            )
        )
    return out


def _density_points(dataset: Dataset, method: str = "multi-fluid") -> list[Comparison]:
    """Compare a mixture density or Z dataset's points in range.

    Parameters
    ----------
    dataset : Dataset
        A mixture compressibility factor or density.
    method : str, optional
        The density method compared: the multi-fluid model by default, or the
        fallback (``virial``), whose points are also kept within its truncation
        bound and its measured temperatures.

    Returns
    -------
    list of Comparison
        One per point in :data:`TEMPERATURES` within the optical-domain bound.
    """
    material, partner = str(dataset.material), str(dataset.parameters["partner"])
    pair = _pair(dataset)
    temperature = float(str(dataset.parameters["temperature_k"]))
    if not TEMPERATURES[0] <= temperature <= TEMPERATURES[1]:
        return []
    fraction = float(str(dataset.parameters["mole_fraction"]))
    composition = {material: fraction, partner: 1 - fraction}
    eos = _quiet_eos(composition, method)
    low, high = eos.temperature_range
    if not low <= temperature <= high:
        return []
    bound = mixture_optics([material, partner]).density_bound(composition)
    model_relative = MODEL_DENSITY_UNCERTAINTY.get(pair, 0.0)
    out = []
    for abscissa, measured, sigma in zip(
        dataset.abscissa_values(),
        dataset.scaled_values(),
        dataset.standard_uncertainty(),
        strict=True,
    ):
        with warnings.catch_warnings():
            # Beyond the truncation bound the fallback warns; those points are left
            # out below.
            warnings.simplefilter("ignore", ModelRangeWarning)
            if dataset.abscissa == "pressure":
                state = State(pressure=float(abscissa), temperature=temperature)
                density = eos.number_density(state)
            else:
                density = float(abscissa)
                state = State(
                    pressure=float(eos.pressure(density, temperature)),
                    temperature=temperature,
                )
        if density > bound:
            continue
        if method == "virial" and (
            abs(density / N_A * eos._cross_sum(temperature)) > TRUNCATION_BOUND  # pyright: ignore[reportPrivateUsage]
        ):
            continue
        if dataset.quantity == "mixture_compressibility_factor":
            model = state.pressure / (density * K_B * temperature)
        else:
            model = _in_unit(density, dataset.unit, composition)
        combined = math.hypot(float(sigma), model_relative * float(measured))
        counted = dataset.use == "primary" and not math.isnan(sigma)
        out.append(
            Comparison(
                dataset=dataset.key,
                pair=pair,
                kind="density",
                temperature=temperature,
                abscissa=float(abscissa),
                measured=float(measured),
                model=float(model),
                sigma=combined,
                counted=counted,
                passed=bool(counted and abs(model - measured) <= 2 * combined),
            )
        )
    return out


_MOLAR_MASS: Final = {"He": 4.002602e-3, "Ne": 20.1797e-3, "Ar": 39.948e-3}


def _in_unit(density: float, unit: str, composition: Mapping[str, float]) -> float:
    """Return a number density in a dataset's density unit.

    Parameters
    ----------
    density : float
        Number density, m^-3.
    unit : str
        The dataset's unit.
    composition : Mapping
        Species to mole fraction, for a mass density.

    Returns
    -------
    float
        In ``unit``.
    """
    molar = density / N_A
    scale = {"mol/m3": 1.0, "mol/cm3": 1e-6, "mol/L": 1e-3}
    if unit in scale:
        return molar * scale[unit]
    mass = molar * sum(x * _MOLAR_MASS[s] for s, x in composition.items())
    return mass if unit == "kg/m3" else mass * 1e-3


@functools.cache
def comparisons() -> tuple[Comparison, ...]:
    """Return every point of the comparison, under TC-778's rule.

    Returns
    -------
    tuple of Comparison
        In the datasets' order; pairs without a fitted multi-fluid model are left
        out.
    """
    out: list[Comparison] = []
    for dataset in datasets():
        if dataset.quantity not in _VIRIALS + _DENSITIES:
            continue
        # Only a pair with a fitted multi-fluid model has one to compare.
        if _binary(*_pair(dataset))[1] != "fitted":
            continue
        if dataset.quantity in _VIRIALS:
            out.extend(_virial_points(dataset))
        else:
            out.extend(_density_points(dataset))
    return tuple(out)


@functools.cache
def fallback_comparisons() -> tuple[Comparison, ...]:
    """Return the fallback's points against the measured mixture densities.

    The same rule as :func:`comparisons`, for the pairs with a cross virial, within
    the fallback's truncation bound and measured temperatures.

    Returns
    -------
    tuple of Comparison
        In the datasets' order.
    """
    records = _pair_records()
    out: list[Comparison] = []
    for dataset in datasets():
        if dataset.quantity not in _DENSITIES:
            continue
        if not records.get(_pair(dataset), {}).get("cross_virial"):
            continue
        out.extend(_density_points(dataset, "virial"))
    return tuple(out)


def _results(points: tuple[Comparison, ...]) -> dict[tuple[str, str], PairResult]:
    """Return each pair's qualification from its points.

    Parameters
    ----------
    points : tuple of Comparison
        The compared points.

    Returns
    -------
    dict
        Pair (sorted) to its result.
    """
    by_pair: dict[tuple[str, str], list[Comparison]] = {}
    for point in points:
        by_pair.setdefault(point.pair, []).append(point)
    out = {}
    for pair, group in sorted(by_pair.items()):
        kinds = ("virial", "density")
        counted = {k: sum(p.counted for p in group if p.kind == k) for k in kinds}
        passed = {k: sum(p.passed for p in group if p.kind == k) for k in kinds}
        enough = sum(counted.values()) >= 3
        fractions = [passed[k] / counted[k] for k in counted if counted[k]]
        qualified = enough and all(f >= 0.95 for f in fractions)
        out[pair] = PairResult(pair, counted, passed, qualified)
    return out


def pair_results() -> dict[tuple[str, str], PairResult]:
    """Return each compared pair's multi-fluid qualification under TC-778's rule.

    Returns
    -------
    dict
        Pair (sorted) to its result.
    """
    return _results(comparisons())


def fallback_results() -> dict[tuple[str, str], PairResult]:
    """Return each pair's fallback qualification against the measured densities.

    Returns
    -------
    dict
        Pair (sorted) to its result; pairs without density data are absent.
    """
    return _results(fallback_comparisons())
