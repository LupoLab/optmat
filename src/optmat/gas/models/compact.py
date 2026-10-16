# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Compact infrared terms for the molecular fits (OEP 2).

A molecular fit (``fit/<gas>``) carries its vibrational bands as fixed poles: thousands
for CO2 and N2O. Its compact model, ``compact/<gas>``, keeps the fit's electronic terms,
covariance and rotational system unchanged and replaces the vibrational poles by a few
effective ones, derived from the stored fit by a constrained least-squares fit to their
contribution over the compact model's intervals. Two moments of the bands are kept
exactly: the optical-side limit ``M_0 = sum f`` and the static contribution
``M_-2 = sum f / E^2``.

The intervals are the parent's qualified intervals inside its optical range
(``range_nm``), optionally ending at a configured wavelength (``longest_nm``); the
infrared windows beyond are left to the parent (OEP 2, Part II). Every compact pole is
bounded to an absorbing region, between two intervals or outside all of them, so none
can lie inside an interval. The fit is weighted by the parent's absolute statistical
one-sigma of ``alpha``, ``D_2`` and ``D_3`` (``D_m = d^m (omega alpha) / d omega^m``),
and the reduction is accepted when every residual on a four-times finer check grid is
below the configured fraction of that sigma.

The reduction's settings are ``data/gas/compact/<gas>.toml``; ``tools/fit/compact.py``
derives and writes ``data/gas/models/compact/<gas>.json``, which :func:`compact_models`
loads. A compact model is never a default.
"""

import dataclasses
import functools
import hashlib
import itertools
import json
import math
import tomllib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar, Final

import numpy as np
from jaxtyping import Float

from optmat import _validity
from optmat.conditions import State
from optmat.constants import AU_ANGULAR_FREQUENCY, AU_POLARIZABILITY, EPS0
from optmat.constants import C as C_LIGHT
from optmat.exceptions import DataError, OptmatError
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.eos import number_density
from optmat.gas.gases import Gas, as_gas
from optmat.gas.infrared import gauss_rule
from optmat.gas.models._defaults import defaults
from optmat.gas.models.fitted import FittedModel, fitted_models
from optmat.identity import FrozenMapping

__all__ = [
    "COMPACT_CONFIG_DIRECTORY",
    "COMPACT_MODEL_DIRECTORY",
    "CompactCoefficients",
    "CompactConfig",
    "CompactModel",
    "Reduction",
    "Region",
    "absorbing_regions",
    "coefficients",
    "compact_intervals",
    "compact_models",
    "compact_record",
    "derive",
    "format_coefficients",
    "load_compact_config",
    "parent_digest",
    "pole_dispersion_au",
    "reproduction",
    "residuals_sigma",
    "statistical_sigma_au",
    "vibrational_terms",
]

COMPACT_CONFIG_DIRECTORY: Final = GAS_DATA_DIRECTORY / "compact"
"""The reduction configurations, ``<gas>.toml``."""

COMPACT_MODEL_DIRECTORY: Final = GAS_DATA_DIRECTORY / "models" / "compact"
"""The stored compact models, ``<gas>.json``."""

_HARTREE_NM: Final = 45.56335252767
"""hc / E_h, nm: a photon energy in hartree is this over its vacuum wavelength in nm."""

_ORDERS: Final = (0, 2, 3)
"""The derivatives the reduction fits: alpha, D_2 and D_3."""

_ORDER_NAMES: Final = {0: "alpha", 2: "D2", 3: "D3"}

_HARTREE_CM: Final = 219474.6313632
"""E_h / hc, cm^-1."""

_MOMENT_TOLERANCE: Final = 1e-10
"""The relative error the inner solve allows on M_0 and M_-2 (OEP 2, section 7)."""


@dataclass(frozen=True)
class CompactConfig:
    """The settings of one reduction (``data/gas/compact/<gas>.toml``).

    Attributes
    ----------
    gas : Gas
        The species.
    description : str
        What the reduction gives.
    parent : str
        The parent fit, ``fit/<gas>``.
    poles : tuple of int
        The number of compact poles in each absorbing region, in ascending energy
        (:func:`absorbing_regions`). A region given at least as many poles as it holds
        keeps the parent's poles unchanged; a region given none is represented by the
        others.
    longest_nm : float or None
        Where the compact intervals end, nm; None for the end of the optical range.
    grid_points : int
        Fitted points per interval, geometric in wavelength.
    check_factor : int
        The check grid has this many times as many points.
    tolerance_sigma : float
        The largest residual allowed, in units of the parent's sigma.
    random_starts : int
        Log-uniform random starting points, beyond the fixed ones.
    seed : int
        The seed of the random starts.
    digest : str
        SHA-256 of the configuration file.
    """

    gas: Gas
    description: str
    parent: str
    poles: tuple[int, ...]
    longest_nm: float | None = None
    grid_points: int = 400
    check_factor: int = 4
    tolerance_sigma: float = 0.1
    random_starts: int = 0
    seed: int = 1
    digest: str = ""


def load_compact_config(name: str | Path) -> CompactConfig:
    """Load a reduction configuration.

    Parameters
    ----------
    name : str or pathlib.Path
        The gas (``data/gas/compact/<gas>.toml``) or a path.

    Returns
    -------
    CompactConfig
        The settings.

    Raises
    ------
    DataError
        If the file has unknown keys or a negative pole count.
    """
    path = Path(name)
    if path.suffix != ".toml":
        path = COMPACT_CONFIG_DIRECTORY / f"{name}.toml"
    raw = path.read_bytes()
    record = tomllib.loads(raw.decode("utf-8"))
    allowed = {f.name for f in dataclasses.fields(CompactConfig)} - {"digest"}
    unknown = set(record) - allowed
    if unknown:
        msg = f"{path.name}: unknown keys {sorted(unknown)}"
        raise DataError(msg)
    poles = tuple(int(k) for k in record["poles"])
    if any(k < 0 for k in poles):
        msg = f"{path.name}: a region cannot have {min(poles)} poles"
        raise DataError(msg)
    longest = record.get("longest_nm")
    return CompactConfig(
        gas=as_gas(record["gas"]),
        description=str(record["description"]),
        parent=str(record["parent"]),
        poles=poles,
        longest_nm=None if longest is None else float(longest),
        grid_points=int(record.get("grid_points", 400)),
        check_factor=int(record.get("check_factor", 4)),
        tolerance_sigma=float(record.get("tolerance_sigma", 0.1)),
        random_starts=int(record.get("random_starts", 0)),
        seed=int(record.get("seed", 1)),
        digest=hashlib.sha256(raw).hexdigest(),
    )


def pole_dispersion_au(
    strengths: Float[np.ndarray, " k"],
    energies: Float[np.ndarray, " k"],
    omega: Float[np.ndarray, " n"],
    order: int,
) -> Float[np.ndarray, "k n"]:
    """Return each pole's ``D_m``: ``alpha`` for order 0, else ``d^m (omega alpha)``.

    With ``omega f / (E^2 - omega^2) = (f / 2) [1 / (E - omega) - 1 / (E + omega)]``,
    ``D_m = (f / 2) m! [(E - omega)^-(m+1) - (-1)^m (E + omega)^-(m+1)]`` for
    ``m >= 1``, on either side of the pole.

    Parameters
    ----------
    strengths : numpy.ndarray
        Oscillator strengths.
    energies : numpy.ndarray
        Energies, hartree.
    omega : numpy.ndarray
        Photon energies, hartree.
    order : int
        0, or the derivative's order ``m >= 1``.

    Returns
    -------
    numpy.ndarray
        One row per pole, atomic units.
    """
    f = np.asarray(strengths, dtype=float)[:, None]
    e = np.asarray(energies, dtype=float)[:, None]
    w = np.asarray(omega, dtype=float)[None, :]
    if order == 0:
        return f / (e**2 - w**2)
    m = order
    return (
        f
        * 0.5
        * math.factorial(m)
        * ((e - w) ** -(m + 1) - (-1) ** m * (e + w) ** -(m + 1))
    )


def vibrational_terms(
    parent: FittedModel,
) -> tuple[Float[np.ndarray, " b"], Float[np.ndarray, " b"]]:
    """Return the parent's vibrational poles: its fixed terms after the rotational.

    Parameters
    ----------
    parent : FittedModel
        The fit.

    Returns
    -------
    tuple of numpy.ndarray
        Strengths and energies (hartree).
    """
    start = len(parent.strengths) - parent.n_fixed + parent.n_rotational
    return np.array(parent.strengths[start:]), np.array(parent.energies[start:])


def _cut(a: float, b: float, low: float, high: float) -> tuple[float, float]:
    """Return an interval cut to a range.

    An end the range does not cut is the interval's own float, not a round trip through
    nm, so that the compact intervals lie inside the parent's exactly.

    Parameters
    ----------
    a : float
        The interval's start, m.
    b : float
        Its end, m.
    low : float
        The range's start, nm.
    high : float
        Its end, nm.

    Returns
    -------
    tuple of float
        The cut interval, m; empty (end before start) if they do not overlap.
    """
    start = a if a * 1e9 >= low * (1 - 1e-12) else low * 1e-9
    end = b if b * 1e9 <= high * (1 + 1e-12) else high * 1e-9
    return start, end


def compact_intervals(
    parent: FittedModel, longest_nm: float | None = None
) -> _validity.Intervals:
    """Return the compact model's intervals: the parent's inside its optical range.

    Parameters
    ----------
    parent : FittedModel
        The fit.
    longest_nm : float or None, optional
        Where the intervals end instead, nm (inside the optical range).

    Returns
    -------
    tuple of tuple of float
        ``(low, high)`` vacuum wavelengths, m, ascending.

    Raises
    ------
    DataError
        If no interval is left.
    """
    low, high = (float(x) for x in parent.record["range_nm"])  # pyright: ignore[reportGeneralTypeIssues]
    if longest_nm is not None:
        high = min(high, longest_nm)

    cut = (_cut(a, b, low, high) for a, b in parent.intervals)
    kept = tuple((a, b) for a, b in cut if b > a)
    if not kept:
        msg = f"{parent.name}: no qualified interval below {high} nm"
        raise DataError(msg)
    return kept


@dataclass(frozen=True)
class Region:
    """An absorbing region: frequencies outside every compact interval, with poles.

    Attributes
    ----------
    low : float
        The region's lowest parent pole, hartree.
    high : float
        Its highest, hartree.
    strengths : numpy.ndarray
        The parent's poles in it: strengths.
    energies : numpy.ndarray
        Their energies, hartree, ascending.
    """

    low: float
    high: float
    strengths: Float[np.ndarray, " r"] = field(repr=False)
    energies: Float[np.ndarray, " r"] = field(repr=False)


def absorbing_regions(
    parent: FittedModel, intervals: _validity.Intervals
) -> tuple[Region, ...]:
    """Return the absorbing regions of the parent's vibrational poles, by energy.

    The frequencies outside every interval divide into the region below the lowest
    interval's frequencies, those between two adjacent intervals and the region above
    the highest (where a fold puts a share of the top band system, OEP 3). Only the
    regions with parent poles count.

    Parameters
    ----------
    parent : FittedModel
        The fit.
    intervals : tuple of tuple of float
        The compact intervals, m.

    Returns
    -------
    tuple of Region
        In ascending energy.

    Raises
    ------
    DataError
        If a parent pole lies inside an interval.
    """
    strengths, energies = vibrational_terms(parent)
    bands = sorted(
        (_HARTREE_NM / (b * 1e9), _HARTREE_NM / (a * 1e9)) for a, b in intervals
    )
    edges = [0.0, *itertools.chain.from_iterable(bands), math.inf]
    regions = []
    for low, high in itertools.pairwise(edges):
        inside = (energies > low) & (energies < high)
        if (low, high) in bands:
            if np.any(inside):
                msg = (
                    f"{parent.name}: a vibrational pole lies inside a compact interval"
                )
                raise DataError(msg)
            continue
        if np.any(inside):
            order = np.argsort(energies[inside])
            e, f = energies[inside][order], strengths[inside][order]
            regions.append(Region(float(e[0]), float(e[-1]), f, e))
    return tuple(regions)


def _grid(intervals: _validity.Intervals, points: int) -> Float[np.ndarray, " n"]:
    """Return photon energies geometric in wavelength over every interval.

    Parameters
    ----------
    intervals : tuple of tuple of float
        The intervals, m.
    points : int
        Points per interval, ends included.

    Returns
    -------
    numpy.ndarray
        Photon energies, hartree.
    """
    return np.concatenate(
        [_HARTREE_NM / np.geomspace(a * 1e9, b * 1e9, points) for a, b in intervals]
    )


def statistical_sigma_au(
    parent: FittedModel, omega: Float[np.ndarray, " n"]
) -> dict[int, Float[np.ndarray, " n"]]:
    """Return the parent's absolute one-sigma of ``alpha``, ``D_2`` and ``D_3``.

    The statistical uncertainty from the covariance, Birge-scaled (OEP 2, section 4).

    Parameters
    ----------
    parent : FittedModel
        The fit.
    omega : numpy.ndarray
        Photon energies, hartree.

    Returns
    -------
    dict of int to numpy.ndarray
        By order (0, 2, 3), atomic units.
    """
    return {
        0: parent.polarizability_sigma_au(omega),
        2: parent.dispersion_sigma_au(omega, 2),
        3: parent.dispersion_sigma_au(omega, 3),
    }


def residuals_sigma(
    parent: FittedModel,
    intervals: _validity.Intervals,
    strengths: Sequence[float],
    energies: Sequence[float],
    points: int,
) -> dict[str, Float[np.ndarray, " n"]]:
    """Return the compact poles' residuals against the parent's, in its sigma.

    Parameters
    ----------
    parent : FittedModel
        The fit.
    intervals : tuple of tuple of float
        The compact intervals, m.
    strengths : sequence of float
        The compact vibrational poles: strengths.
    energies : sequence of float
        Their energies, hartree.
    points : int
        Points per interval.

    Returns
    -------
    dict of str to numpy.ndarray
        ``alpha``, ``D2`` and ``D3``: (compact - parent) / sigma at each point.
    """
    omega = _grid(intervals, points)
    sigma = statistical_sigma_au(parent, omega)
    f_parent, e_parent = vibrational_terms(parent)
    f, e = np.asarray(strengths, dtype=float), np.asarray(energies, dtype=float)
    return {
        _ORDER_NAMES[m]: (
            pole_dispersion_au(f, e, omega, m).sum(axis=0)
            - pole_dispersion_au(f_parent, e_parent, omega, m).sum(axis=0)
        )
        / sigma[m]
        for m in _ORDERS
    }


@dataclass(frozen=True)
class Reduction:
    """A derived set of compact vibrational poles.

    Attributes
    ----------
    strengths : tuple of float
        The poles' strengths, ascending in energy.
    energies : tuple of float
        Their energies, hartree.
    region : tuple of int
        The absorbing region of each pole.
    kept : tuple of bool
        Whether each pole is the parent's own (a region kept unchanged).
    worst : dict of str to float
        The largest ``|residual| / sigma`` on the check grid, per quantity.
    moment_errors : tuple of float
        The relative errors of ``M_0`` and ``M_-2``.
    start : str
        The starting point the chosen solution came from.
    passed : bool
        Whether every residual is below the tolerance.
    """

    strengths: tuple[float, ...]
    energies: tuple[float, ...]
    region: tuple[int, ...]
    kept: tuple[bool, ...]
    worst: Mapping[str, float]
    moment_errors: tuple[float, float]
    start: str
    passed: bool


def _inner(
    design: Float[np.ndarray, "n k"],
    target: Float[np.ndarray, " n"],
    constraint: Float[np.ndarray, "2 k"],
) -> tuple[float, Float[np.ndarray, " k"]] | None:
    """Solve the equality-constrained non-negative least squares exactly.

    Minimise ``|A x - b|`` subject to ``C x = (1, 1)`` and ``x >= 0`` by enumerating the
    supports; on each, the null-space method (OEP 2, section 7).

    Parameters
    ----------
    design : numpy.ndarray
        ``A``, one column per pole.
    target : numpy.ndarray
        ``b``.
    constraint : numpy.ndarray
        ``C``: the scaled ``M_0`` and ``M_-2`` rows.

    Returns
    -------
    tuple or None
        The cost and the solution, or None if no support is feasible.
    """
    k = design.shape[1]
    ones = np.ones(2)
    best: tuple[float, np.ndarray] | None = None
    for size in range(1, k + 1):
        for support in itertools.combinations(range(k), size):
            columns = list(support)
            a, c = design[:, columns], constraint[:, columns]
            particular = np.linalg.lstsq(c, ones, rcond=None)[0]
            rank = np.linalg.matrix_rank(c)
            null = np.linalg.svd(c)[2][rank:].T
            x = particular
            if null.shape[1]:
                step = np.linalg.lstsq(a @ null, target - a @ particular, rcond=None)[0]
                x = particular + null @ step
            if np.any(x < 0) or np.max(np.abs(c @ x - ones)) > _MOMENT_TOLERANCE:
                continue
            cost = float(np.sum((a @ x - target) ** 2))
            if best is None or cost < best[0]:
                solution = np.zeros(k)
                solution[columns] = x
                best = (cost, solution)
    return best


@dataclass(frozen=True)
class _Problem:
    """The reduction's fixed quantities on one grid."""

    omega: Float[np.ndarray, " n"]
    sigma: Mapping[int, Float[np.ndarray, " n"]]
    target: Float[np.ndarray, " m"]
    m0: float
    m_minus2: float

    def design(self, energies: Float[np.ndarray, " k"]) -> Float[np.ndarray, "m k"]:
        """Return the weighted ``D_m`` of unit-strength poles, times ``M_0``.

        Parameters
        ----------
        energies : numpy.ndarray
            The poles' energies, hartree.

        Returns
        -------
        numpy.ndarray
            One column per pole.
        """
        unit = np.ones_like(energies)
        rows = [
            pole_dispersion_au(unit, energies, self.omega, m) / self.sigma[m]
            for m in _ORDERS
        ]
        return np.concatenate(rows, axis=1).T * self.m0

    def constraint(self, energies: Float[np.ndarray, " k"]) -> Float[np.ndarray, "2 k"]:
        """Return the ``M_0`` and ``M_-2`` rows for strengths in units of ``M_0``.

        Parameters
        ----------
        energies : numpy.ndarray
            The poles' energies, hartree.

        Returns
        -------
        numpy.ndarray
            Both rows equal one at a solution.
        """
        return np.vstack(
            [np.ones_like(energies), self.m0 / self.m_minus2 / energies**2]
        )


def _problem(
    parent: FittedModel,
    intervals: _validity.Intervals,
    points: int,
    kept: tuple[Float[np.ndarray, " j"], Float[np.ndarray, " j"]],
) -> _Problem:
    """Return the reduction's fixed quantities on a grid.

    Parameters
    ----------
    parent : FittedModel
        The fit.
    intervals : tuple of tuple of float
        The compact intervals, m.
    points : int
        Points per interval.
    kept : tuple of numpy.ndarray
        The strengths and energies of the parent's poles kept unchanged.

    Returns
    -------
    _Problem
        The grid, the sigmas, and the target and moments the searched poles take on:
        the parent's vibrational poles less the kept ones.
    """
    omega = _grid(intervals, points)
    sigma = statistical_sigma_au(parent, omega)
    f, e = vibrational_terms(parent)
    g, h = kept
    target = np.concatenate(
        [
            (
                pole_dispersion_au(f, e, omega, m).sum(axis=0)
                - pole_dispersion_au(g, h, omega, m).sum(axis=0)
            )
            / sigma[m]
            for m in _ORDERS
        ]
    )
    m0 = float(f.sum() - g.sum())
    m_minus2 = float(np.sum(f / e**2) - np.sum(g / h**2))
    return _Problem(omega, sigma, target, m0, m_minus2)


def _per_region(
    regions: Sequence[Region], counts: Sequence[int], rule: Any
) -> Float[np.ndarray, " k"]:
    """Return a start's energies, region by region.

    Parameters
    ----------
    regions : sequence of Region
        The searched regions.
    counts : sequence of int
        Their pole counts.
    rule : callable
        ``rule(region, count)``, the region's starting energies.

    Returns
    -------
    numpy.ndarray
        The energies, hartree, ascending within each region.
    """
    return np.concatenate(
        [np.sort(rule(r, k)) for r, k in zip(regions, counts, strict=True)]
    )


def _spread(region: Region, count: int) -> Float[np.ndarray, " k"]:
    """Return energies spread geometrically inside a region.

    Parameters
    ----------
    region : Region
        The region.
    count : int
        How many.

    Returns
    -------
    numpy.ndarray
        Energies, hartree, strictly between its lowest and highest pole.
    """
    return np.geomspace(region.low, region.high, count + 2)[1:-1]


def _ends(region: Region, count: int) -> Float[np.ndarray, " k"]:
    """Return energies spread geometrically from a region's lowest pole to its highest.

    Parameters
    ----------
    region : Region
        The region.
    count : int
        How many.

    Returns
    -------
    numpy.ndarray
        Energies, hartree.
    """
    return np.geomspace(region.low, region.high, count)


def _gauss_start(region: Region, count: int) -> Float[np.ndarray, " k"]:
    """Return the nodes of the Gaussian rule of a region's poles (OEP 2, section 4).

    Parameters
    ----------
    region : Region
        The region.
    count : int
        How many nodes.

    Returns
    -------
    numpy.ndarray
        Energies, hartree.
    """
    return gauss_rule(region.strengths, region.energies, count)[1]


def _strongest(region: Region, count: int) -> Float[np.ndarray, " k"]:
    """Return the energies of a region's strongest poles, at least 5 % apart.

    Parameters
    ----------
    region : Region
        The region.
    count : int
        How many; spread energies fill in when too few are 5 % apart.

    Returns
    -------
    numpy.ndarray
        Energies, hartree.
    """
    chosen: list[float] = []
    for index in np.argsort(-region.strengths):
        energy = float(region.energies[index])
        if all(abs(energy / other - 1) >= 0.05 for other in chosen):
            chosen.append(energy)
        if len(chosen) == count:
            break
    return np.concatenate([chosen, _spread(region, count)[len(chosen) :]])


def _starts(
    regions: Sequence[Region], counts: Sequence[int], random: int, seed: int
) -> list[tuple[str, Float[np.ndarray, " k"]]]:
    """Return the labelled starting energies of the searched poles (OEP 2, section 7).

    Parameters
    ----------
    regions : sequence of Region
        The searched regions.
    counts : sequence of int
        Their pole counts.
    random : int
        How many log-uniform random starts.
    seed : int
        Their seed.

    Returns
    -------
    list of tuple
        ``(label, energies)``, ascending within each region.
    """
    starts = [
        ("gauss", _per_region(regions, counts, _gauss_start)),
        ("strongest", _per_region(regions, counts, _strongest)),
        ("spread", _per_region(regions, counts, _spread)),
        ("spread to the ends", _per_region(regions, counts, _ends)),
    ]
    generator = np.random.default_rng(seed)
    for n in range(random):
        energies = _per_region(
            regions,
            counts,
            lambda r, k: np.exp(
                generator.uniform(math.log(r.low), math.log(r.high), k)
            ),
        )
        starts.append((f"random {n + 1}", energies))
    return starts


def _bounds(
    regions: Sequence[Region], counts: Sequence[int], mean: float
) -> tuple[Float[np.ndarray, " k"], Float[np.ndarray, " k"]]:
    """Return the searched poles' log-energy bounds, feasible for both moments.

    A solution with ``F >= 0`` keeps ``M_0`` and ``M_-2`` exactly if and only if some
    pole lies at or below ``mean = (M_0 / M_-2)^(1/2)`` and some at or above it (OEP 2,
    section 7): the lowest searched pole is bounded below it and the highest above it.

    Parameters
    ----------
    regions : sequence of Region
        The searched regions, at least one.
    counts : sequence of int
        Their pole counts.
    mean : float
        ``(M_0 / M_-2)^(1/2)`` of the poles they take on, hartree.

    Returns
    -------
    tuple of numpy.ndarray
        Lower and upper bounds on ``log E``.

    Raises
    ------
    DataError
        If the bounds cannot straddle the mean.
    """
    pairs = zip(regions, counts, strict=True)
    low = np.array([math.log(r.low) for r, k in pairs for _ in range(k)], dtype=float)
    pairs = zip(regions, counts, strict=True)
    high = np.array([math.log(r.high) for r, k in pairs for _ in range(k)], dtype=float)
    high[0] = min(high[0], math.log(mean))
    low[-1] = max(low[-1], math.log(mean))
    if not np.all(low <= high):
        msg = "the compact poles cannot keep both moments: give the regions more poles"
        raise DataError(msg)
    return low, high


def _worst(
    problem: _Problem, energies: Float[np.ndarray, " k"], x: Float[np.ndarray, " k"]
) -> dict[str, float]:
    """Return the largest weighted residual of each quantity.

    Parameters
    ----------
    problem : _Problem
        The check grid.
    energies : numpy.ndarray
        The searched poles' energies, hartree.
    x : numpy.ndarray
        Their strengths, in units of ``M_0``.

    Returns
    -------
    dict of str to float
        ``|residual| / sigma`` for ``alpha``, ``D2`` and ``D3``.
    """
    residual = problem.design(energies) @ x - problem.target
    n = problem.omega.size
    return {
        _ORDER_NAMES[m]: float(np.max(np.abs(residual[j * n : (j + 1) * n])))
        for j, m in enumerate(_ORDERS)
    }


def _search(
    fit: _Problem,
    check: _Problem,
    searched: Sequence[tuple[Region, int]],
    config: CompactConfig,
) -> list[tuple[float, str, np.ndarray, np.ndarray, dict[str, float]]]:
    """Search the poles' energies from every start (OEP 2, section 7).

    Parameters
    ----------
    fit : _Problem
        The fitted grid.
    check : _Problem
        The check grid.
    searched : sequence of tuple
        The searched regions and their pole counts.
    config : CompactConfig
        The settings: the random starts and their seed.

    Returns
    -------
    list of tuple
        Per start: the cost, its label, the energies (hartree), the strengths and the
        largest residuals on the check grid.
    """
    candidates = []
    s_regions = [r for r, _ in searched]
    s_counts = [k for _, k in searched]
    low, high = _bounds(s_regions, s_counts, math.sqrt(fit.m0 / fit.m_minus2))
    free = low < high

    def energies_of(x: Float[np.ndarray, " p"]) -> Float[np.ndarray, " k"]:
        """Return the searched poles' energies from the free log energies.

        Parameters
        ----------
        x : numpy.ndarray
            The free poles' log energies.

        Returns
        -------
        numpy.ndarray
            Every searched pole's energy, hartree.
        """
        log_e = low.copy()
        log_e[free] = x
        return np.exp(log_e)

    def solve(
        energies: Float[np.ndarray, " k"],
    ) -> tuple[float, Float[np.ndarray, " k"]]:
        """Return the exact strengths for some energies.

        Parameters
        ----------
        energies : numpy.ndarray
            The poles' energies, hartree, within the bounds.

        Returns
        -------
        tuple
            The cost and the strengths in units of ``M_0``.
        """
        # Feasible for every energies within the bounds (_bounds).
        solution = _inner(fit.design(energies), fit.target, fit.constraint(energies))
        assert solution is not None
        return solution

    def residual(x: Float[np.ndarray, " p"]) -> Float[np.ndarray, " m"]:
        """Return the weighted residual at the exact strengths.

        Parameters
        ----------
        x : numpy.ndarray
            The free poles' log energies.

        Returns
        -------
        numpy.ndarray
            The residual vector on the fitted grid.
        """
        energies = energies_of(x)
        return fit.design(energies) @ solve(energies)[1] - fit.target

    from scipy import optimize  # noqa: PLC0415 - needed only to derive a model

    inset = 1e-9 * (high - low)
    for label, start in _starts(s_regions, s_counts, config.random_starts, config.seed):
        x0 = np.clip(np.log(start), low + inset, high - inset)[free]
        if x0.size:
            bounds = (low[free], high[free])
            x0 = optimize.least_squares(residual, x0, bounds=bounds, max_nfev=600).x
        energies = energies_of(x0)
        cost, x = solve(energies)
        candidates.append(
            (cost, label, energies, x * fit.m0, _worst(check, energies, x))
        )
        if not x0.size:
            break  # every searched pole is fixed by its bounds
    return candidates


def derive(parent: FittedModel, config: CompactConfig) -> Reduction:
    """Derive the compact vibrational poles of a fit.

    A region given at least as many poles as it holds keeps the parent's; the others'
    poles are searched, taking on the rest of the target and of both moments.

    Parameters
    ----------
    parent : FittedModel
        The fit.
    config : CompactConfig
        The reduction's settings.

    Returns
    -------
    Reduction
        The chosen poles: of the candidates from every start, the lowest cost among
        those meeting the tolerance on the check grid, else the lowest cost.

    Raises
    ------
    DataError
        If the pole counts do not match the regions, or the poles cannot keep both
        moments.
    """
    intervals = compact_intervals(parent, config.longest_nm)
    regions = absorbing_regions(parent, intervals)
    if len(config.poles) != len(regions):
        found = ", ".join(
            f"{r.low * _HARTREE_CM:.1f}-{r.high * _HARTREE_CM:.1f} cm-1"
            for r in regions
        )
        msg = (
            f"{parent.name}: {len(config.poles)} pole counts for {len(regions)} "
            f"regions ({found})"
        )
        raise DataError(msg)
    keep = [
        k >= len(np.unique(r.energies))
        for r, k in zip(regions, config.poles, strict=True)
    ]
    held = [r for r, k in zip(regions, keep, strict=True) if k]
    kept = (
        np.concatenate([np.zeros(0), *(r.strengths for r in held)]),
        np.concatenate([np.zeros(0), *(r.energies for r in held)]),
    )
    kept_region = [i for i, r in enumerate(regions) if keep[i] for _ in r.energies]
    searched = [
        (i, r, k)
        for i, (r, k) in enumerate(zip(regions, config.poles, strict=True))
        if k and not keep[i]
    ]
    fit = _problem(parent, intervals, config.grid_points, kept)
    check = _problem(parent, intervals, config.grid_points * config.check_factor, kept)
    total_f, total_e = vibrational_terms(parent)
    candidates = []
    if not searched:
        # Every region kept, or the rest left to no pole at all.
        if abs(fit.m0) > 1e-12 * total_f.sum():
            msg = "the compact poles cannot keep both moments: give the regions poles"
            raise DataError(msg)
        worst = _worst(check, np.zeros(0), np.zeros(0))
        candidates.append((0.0, "the parent's poles", np.zeros(0), np.zeros(0), worst))
    else:
        candidates = _search(fit, check, [(r, k) for _, r, k in searched], config)
    tolerance = config.tolerance_sigma
    _, label, energies, strengths, worst = min(
        candidates, key=lambda c: (max(c[4].values()) >= tolerance, c[0])
    )
    s_index = [i for i, _, k in searched for _ in range(k)]
    all_f = np.concatenate([kept[0], strengths])
    all_e = np.concatenate([kept[1], energies])
    regions_of = np.array([*kept_region, *s_index], dtype=int)
    is_kept = np.arange(all_e.size) < kept[1].size
    order = np.argsort(all_e, kind="stable")
    order = order[all_f[order] > 0]
    f, e = all_f[order], all_e[order]
    return Reduction(
        strengths=tuple(float(x) for x in f),
        energies=tuple(float(x) for x in e),
        region=tuple(int(x) for x in regions_of[order]),
        kept=tuple(bool(x) for x in is_kept[order]),
        worst=worst,
        moment_errors=(
            float(f.sum() / total_f.sum() - 1),
            float(np.sum(f / e**2) / np.sum(total_f / total_e**2) - 1),
        ),
        start=label,
        passed=max(worst.values()) < tolerance,
    )


def parent_digest(parent: str) -> str:
    """Return the SHA-256 of a stored parent fit.

    Parameters
    ----------
    parent : str
        ``fit/<gas>``.

    Returns
    -------
    str
        The digest of ``data/gas/models/fitted/<gas>.json``.
    """
    name = parent.removeprefix("fit/")
    path = GAS_DATA_DIRECTORY / "models" / "fitted" / f"{name}.json"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact_record(
    parent: FittedModel, config: CompactConfig, reduction: Reduction
) -> dict[str, object]:
    """Return a reduction as the stored compact model's record.

    Parameters
    ----------
    parent : FittedModel
        The fit.
    config : CompactConfig
        The settings.
    reduction : Reduction
        The derived poles.

    Returns
    -------
    dict
        JSON-serialisable: the parent's identity, the settings' digest, the intervals,
        the regions, the poles, the moments and the achieved residuals.
    """
    intervals = compact_intervals(parent, config.longest_nm)
    regions = absorbing_regions(parent, intervals)
    # A gap between two compact intervals is one of the parent's: the intervals are
    # the parent's own, cut only at the ends.
    by_start = {
        low: reason
        for (_, low), reason in zip(parent.intervals, parent.gap_reasons, strict=False)
    }
    reasons = [by_start[end] for (_, end), _ in itertools.pairwise(intervals)]
    f, e = vibrational_terms(parent)
    return {
        "name": config.gas.value,
        "gas": config.gas.value,
        "parent": config.parent,
        "parent_sha256": parent_digest(config.parent),
        "parent_config_sha256": str(parent.record["config_sha256"]),
        "config_sha256": config.digest,
        "longest_nm": config.longest_nm,
        "intervals_nm": [[a * 1e9, b * 1e9] for a, b in intervals],
        "gap_reasons": reasons,
        "regions_hartree": [[r.low, r.high] for r in regions],
        "poles": list(config.poles),
        "strengths": list(reduction.strengths),
        "energies_hartree": list(reduction.energies),
        "region": list(reduction.region),
        "kept": list(reduction.kept),
        "parent_moments": {"M0": float(f.sum()), "M-2": float(np.sum(f / e**2))},
        "moment_relative_errors": list(reduction.moment_errors),
        "tolerance_sigma": config.tolerance_sigma,
        "check_points": config.grid_points * config.check_factor,
        "residuals_sigma": dict(reduction.worst),
        "passed": reduction.passed,
        "start": reduction.start,
        "n_rotational": parent.n_rotational,
        "reference_temperature_k": parent.reference_temperature,
    }


@dataclass(frozen=True)
class CompactModel(FittedModel):
    """A fit with its vibrational poles replaced by a few compact ones (OEP 2).

    The electronic terms, covariance and rotational system are the parent's; the
    trailing fixed terms are the compact poles. Its intervals are the parent's optical
    ones: beyond them, the range warning names the parent, which is qualified in its
    infrared windows.
    """

    EVALUATION_VERSION: ClassVar[int] = 1

    @property
    def parent(self) -> str:
        """The parent fit's name.

        Returns
        -------
        str
            ``fit/<gas>``.
        """
        return str(self.record["parent"])

    @property
    def residuals_sigma(self) -> dict[str, float]:
        """The reduction's largest residuals, in the parent's sigma.

        Returns
        -------
        dict of str to float
            ``alpha``, ``D2`` and ``D3`` on the check grid.
        """
        stored = self.record["residuals_sigma"]
        assert isinstance(stored, Mapping)
        return {str(k): float(v) for k, v in stored.items()}  # pyright: ignore[reportUnknownArgumentType, reportUnknownVariableType]

    def _gap_note(self, wavelength: Float[np.ndarray, "*shape"]) -> str:
        """Return the gap's reason, and the parent where it is qualified.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.

        Returns
        -------
        str
            The gap's reason, then ``"; <parent> is qualified there"`` for any
            wavelength in the parent's intervals.
        """
        note = super()._gap_note(wavelength)
        values = np.ravel(np.asarray(wavelength, dtype=float))
        parent = fitted_models()[self.parent]
        theirs = ~_validity.outside(values, parent.intervals)
        if np.any(theirs & _validity.outside(values, self.intervals)):
            note += f"; {self.parent} is qualified there"
        return note

    def _sources(self) -> tuple[str, ...]:
        """Return the parent's sources.

        Returns
        -------
        tuple of str
            Sorted source ids.
        """
        return fitted_models()[self.parent]._sources()


@functools.cache
def compact_models() -> dict[str, CompactModel]:
    """Return every stored compact model, by name.

    Returns
    -------
    dict of str to CompactModel
        Models named ``compact/<gas>``.
    """
    models = {}
    for path in sorted(COMPACT_MODEL_DIRECTORY.glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        parent = fitted_models()[record["parent"]]
        n_terms = len(parent.strengths) - parent.n_fixed + parent.n_rotational
        labels = [
            "compact (the parent's)" if kept else "compact" for kept in record["kept"]
        ]
        labels = [f"{label} {k + 1}" for k, label in enumerate(labels)]
        # The parent's own floats: the stored nm are for reading.
        intervals = compact_intervals(parent, record["longest_nm"])
        gas = as_gas(record["gas"])
        name = f"compact/{record['name']}"
        default = defaults().get(gas, (record["parent"], ""))[0]
        models[name] = CompactModel(
            gas=gas,
            name=name,
            status="comparison",
            status_reason=(
                f"derived from {record['parent']} for print (OEP 2); never a default, "
                f"the default of {gas.value} is {default}"
            ),
            citation=f"derived here from {record['parent']} ({record['name']}.toml, "
            f"sha256 {record['config_sha256'][:12]})",
            wavelength_range=(intervals[0][0], intervals[-1][1]),
            intervals=intervals,
            kind="compact",
            strengths=(*parent.strengths[:n_terms], *record["strengths"]),
            energies=(*parent.energies[:n_terms], *record["energies_hartree"]),
            labels=(*parent.labels[:n_terms], *labels),
            offset_au=parent.offset_au,
            retardation_au=parent.retardation_au,
            absorption=parent.absorption,
            n_lines=parent.n_lines,
            n_fixed=parent.n_rotational + len(record["strengths"]),
            covariance=np.array(parent.covariance),
            reduced_chi2=parent.reduced_chi2,
            record=record,
        )
    return models


def _inside(
    omega: Float[np.ndarray, " n"], intervals: _validity.Intervals
) -> np.ndarray:
    """Return where photon energies lie in an interval or above the first's frequencies.

    Parameters
    ----------
    omega : numpy.ndarray
        Photon energies, hartree.
    intervals : tuple of tuple of float
        The intervals, m.

    Returns
    -------
    numpy.ndarray
        True in an interval or at a shorter wavelength than every interval.
    """
    wavelength = _HARTREE_NM / np.asarray(omega) * 1e-9
    return ~_validity.outside(wavelength, intervals) | (wavelength < intervals[0][0])


def reproduction(model: CompactModel) -> dict[str, dict[str, float]]:
    """Re-evaluate the parent's datasets and validators with the compact poles.

    Each fitted dataset and validator goes through the fit's own observation code at
    its own wavelengths (harmonics included), with the parent's fitted scales and
    partner models fixed (OEP 2, sections 5 and 7).

    Parameters
    ----------
    model : CompactModel
        The compact model.

    Returns
    -------
    dict
        Per dataset: ``n`` (the values in the compact intervals or at shorter
        wavelengths), ``max_sigma`` (their largest change in their own sigma),
        ``chi2_change``, and ``beyond_max_sigma`` for the values outside (reported,
        not bounded). A static value (zero frequency) has ``static_relative``.
    """
    # Imported here: the fitting code loads every dataset.
    from optmat.gas import fitting  # noqa: PLC0415

    parent = fitted_models()[model.parent]
    config = fitting.load_config(str(parent.record["name"]))
    start = len(model.strengths) - model.n_fixed + model.n_rotational
    compact = dataclasses.replace(
        config,
        infrared=tuple(
            (label, f, e)
            for label, f, e in zip(
                model.labels[start:],
                model.strengths[start:],
                model.energies[start:],
                strict=True,
            )
        ),
    )
    n_fitted = len(parent.strengths) - parent.n_fixed
    strengths = np.asarray(parent.strengths)[:n_fitted]
    energies = np.asarray(parent.energies)[:n_fitted]
    scales = parent.record["scales"]
    assert isinstance(scales, Mapping)
    found: dict[str, fitting.Observation | str] = {
        o.key: o for o in fitting.observations(config)
    }
    for key, obs in fitting._validator_observations(config).items():  # pyright: ignore[reportPrivateUsage]
        # A dataset fitted in one window and compared in another is reported twice.
        found[f"{key} (validator)" if key in found else key] = obs
    report = {}
    for key, obs in found.items():
        if isinstance(obs, str) or obs.orders is not None:
            continue  # no comparison, or electronic sums the infrared does not enter
        dataset = key.removesuffix(" (validator)")
        scale = float(scales.get(dataset, 1.0))  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
        before = scale * fitting._model_values(strengths, energies, obs, config)  # pyright: ignore[reportPrivateUsage]
        after = scale * fitting._model_values(strengths, energies, obs, compact)  # pyright: ignore[reportPrivateUsage]
        entry = {}
        static = obs.omega == 0
        if np.any(static):
            ratio = after[static] / before[static] - 1
            entry["static_relative"] = float(np.max(np.abs(ratio)))
        if np.any(~static):
            omega = np.where(static, 1.0, obs.omega)  # a static value is not dynamic
            inside = _inside(omega, model.intervals) & ~static
            if obs.harmonic:
                inside &= _inside(obs.harmonic * omega, model.intervals)
            stated = obs.sigma > 1e-200
            change = (after - before) / np.where(stated, obs.sigma, np.inf)
            residual = (before - obs.alpha) / np.where(stated, obs.sigma, np.inf)
            entry["n"] = float(np.count_nonzero(inside))
            entry["max_sigma"] = float(np.max(np.abs(change[inside]), initial=0.0))
            entry["chi2_change"] = float(
                np.sum((2 * residual * change + change**2)[inside])
            )
            beyond = ~inside & ~static
            if np.any(beyond):
                entry["beyond_max_sigma"] = float(np.max(np.abs(change[beyond])))
        report[key] = entry
    return report


_KINDS: Final = ("line", "effective", "rotational", "infrared")
"""The kinds of a printed term, in the order of a compact model's terms."""

_DIGITS: Final = 11
"""Significant digits of the printed coefficients: the rounding then changes alpha, D_2
and D_3 by at most 0.002 of the parent's sigma (CO2's poles near 1.43 um need the
eleventh digit; OEP 2 found ten enough on the earlier parents)."""


@dataclass(frozen=True)
class CompactCoefficients:
    """A compact model's coefficients, for print (OEP 2, section 6).

    Terms are in descending energy (ascending wavelength): the electronic ones, then the
    infrared ones and the rotational term. In the polarizability form
    ``alpha(lambda) = chi0 + (e^2 / m_e) sum_k f_k / (omega_k^2 - omega^2)``, with
    ``omega_k = 2 pi c / lambda_k``, and the index follows through the package's density
    conversion ({doc}`conventions`). In the Sellmeier form, at one state,
    ``n - 1 = constant + sum_k B_k lambda^2 / (lambda^2 - C_k)``: the dilute limit,
    without Lorentz-Lorenz or the refractivity virial, whose error against the full
    parent is stated.

    Attributes
    ----------
    model : str
        ``compact/<gas>``.
    parent : str
        ``fit/<gas>``.
    model_digest : str
        The compact model's identity digest.
    parent_digest : str
        The parent's identity digest.
    form : str
        ``polarizability`` or ``sellmeier``.
    labels : tuple of str
        The terms' labels.
    kinds : tuple of str
        ``line``, ``effective``, ``infrared`` or ``rotational``.
    wavelengths : tuple of float
        ``lambda_k``, vacuum, m.
    strengths : tuple of float
        ``f_k``, dimensionless (the polarizability form; empty otherwise).
    chi0 : float
        The constant term of the polarizability, C m^2/V (zero for these molecules).
    B : tuple of float
        The Sellmeier strengths, dimensionless (the Sellmeier form; empty otherwise).
    C : tuple of float
        The Sellmeier resonances ``lambda_k^2``, m^2 (the Sellmeier form).
    constant : float
        The Sellmeier form's constant, ``N chi0 / (2 epsilon_0)``.
    pressure : float or None
        The Sellmeier form's state: pressure, Pa.
    temperature : float or None
        Temperature, K.
    number_density : float or None
        The real-gas number density at that state, m^-3 (:mod:`optmat.gas.eos`).
    intervals : tuple of tuple of float
        Where the coefficients are qualified, vacuum wavelengths, m; between them the
        gas absorbs.
    residuals_sigma : Mapping
        The reduction's largest residuals, in the parent's sigma.
    stated_error : Mapping
        The Sellmeier form's error against the full parent at its state over the
        intervals: ``n_minus_1``, the largest relative difference in ``n - 1``, and
        ``beta2``, the largest difference in the group-velocity dispersion relative to
        its largest value (a relative difference is unbounded at a zero of ``beta2``).
        Empty for the polarizability form.
    covariance : numpy.ndarray
        The covariance of the electronic terms' printed parameters,
        :attr:`covariance_parameters`, Birge-scaled as the parent's.
    covariance_parameters : tuple of str
        What each row of the covariance is: ``f <label>`` and ``lambda <label>`` (m) in
        the polarizability form, ``B <label>`` and ``C <label>`` (m^2) in the Sellmeier
        form.
    """

    model: str
    parent: str
    model_digest: str
    parent_digest: str
    form: str
    labels: tuple[str, ...]
    kinds: tuple[str, ...]
    wavelengths: tuple[float, ...]
    strengths: tuple[float, ...]
    chi0: float
    B: tuple[float, ...]
    C: tuple[float, ...]
    constant: float
    pressure: float | None
    temperature: float | None
    number_density: float | None
    intervals: tuple[tuple[float, float], ...]
    residuals_sigma: Mapping[str, float]
    stated_error: Mapping[str, float]
    covariance: Float[np.ndarray, "q q"] = field(repr=False)
    covariance_parameters: tuple[str, ...] = field(repr=False)

    def __post_init__(self) -> None:
        """Make the covariance read-only and the mappings immutable."""
        covariance = np.array(self.covariance, dtype=float)
        covariance.setflags(write=False)
        object.__setattr__(self, "covariance", covariance)
        object.__setattr__(self, "residuals_sigma", FrozenMapping(self.residuals_sigma))
        object.__setattr__(self, "stated_error", FrozenMapping(self.stated_error))


def _beta2(
    refractivity: Any, wavelength: Float[np.ndarray, " n"]
) -> Float[np.ndarray, " n"]:
    """Return the group-velocity dispersion by a centred second difference.

    Parameters
    ----------
    refractivity : callable
        ``n - 1`` of the vacuum wavelength, m.
    wavelength : numpy.ndarray
        Vacuum wavelengths, m.

    Returns
    -------
    numpy.ndarray
        ``beta_2``, s^2/m.
    """
    omega = 2 * np.pi * C_LIGHT / wavelength
    step = omega * 1e-3
    # The vacuum part of beta = omega n / c is linear in omega and drops out of the
    # second difference; leaving it out avoids cancelling it in floating point.
    beta = [
        (omega + j * step)
        / C_LIGHT
        * refractivity(2 * np.pi * C_LIGHT / (omega + j * step))
        for j in (-1, 0, 1)
    ]
    return (beta[0] - 2 * beta[1] + beta[2]) / step**2


def _printed_covariance(
    model: CompactModel, electronic: Sequence[int], scale: float | None
) -> tuple[Float[np.ndarray, "q q"], tuple[str, ...]]:
    """Return the electronic terms' covariance in their printed parameters.

    The parent's covariance is over the fitted strengths and the effective terms'
    energies (Birge-scaled here). The printed ones are ``f`` and ``lambda`` (m), with
    ``d lambda / d E = -lambda / E``; or, in the Sellmeier form, ``B = scale f / E^2``
    and ``C = lambda^2``, with ``dB/df = B/f``, ``dB/dE = -2B/E`` and
    ``dC/dE = -2C/E``.

    Parameters
    ----------
    model : CompactModel
        The compact model.
    electronic : sequence of int
        The electronic terms' indices, in printed order.
    scale : float or None
        ``N e^2 / (2 epsilon_0 m_e)`` in atomic units for the Sellmeier form; None for
        the polarizability form.

    Returns
    -------
    tuple
        The covariance, and what each of its rows is.
    """
    covariance = np.asarray(model.covariance)
    if model.reduced_chi2 > 1:
        covariance = covariance * model.reduced_chi2
    n_fitted = len(model.strengths) - model.n_fixed
    effective = [k for k in electronic if k >= model.n_lines]
    energy_row = {k: n_fitted + k - model.n_lines for k in effective}
    jacobian = np.zeros((len(electronic) + len(effective), covariance.shape[0]))
    names = ("f", "lambda") if scale is None else ("B", "C")
    for row, k in enumerate(electronic):
        f, e = model.strengths[k], model.energies[k]
        printed = 1.0 if scale is None else scale / e**2
        jacobian[row, k] = printed
        if k in energy_row and scale is not None:
            jacobian[row, energy_row[k]] = -2 * printed * f / e
    for j, k in enumerate(effective):
        e = model.energies[k]
        lam = 2 * np.pi * C_LIGHT / (e * AU_ANGULAR_FREQUENCY)
        printed = lam if scale is None else lam**2
        factor = 1.0 if scale is None else 2.0
        jacobian[len(electronic) + j, energy_row[k]] = -factor * printed / e
    parameters = tuple(
        [f"{names[0]} {model.labels[k]}" for k in electronic]
        + [f"{names[1]} {model.labels[k]}" for k in effective]
    )
    return jacobian @ covariance @ jacobian.T, parameters


def _stated_error(
    model: CompactModel,
    state: State,
    b: Float[np.ndarray, " k"],
    c: Float[np.ndarray, " k"],
    constant: float,
) -> dict[str, float]:
    """Return the Sellmeier form's error against the full parent at its state.

    Parameters
    ----------
    model : CompactModel
        The compact model.
    state : State
        The form's state.
    b : numpy.ndarray
        Its strengths.
    c : numpy.ndarray
        Its resonances, m^2.
    constant : float
        Its constant.

    Returns
    -------
    dict of str to float
        ``n_minus_1``, the largest relative difference in ``n - 1``, and ``beta2``, the
        largest difference in ``beta_2`` over the largest ``|beta_2|``, on the check
        grid over every interval.
    """
    parent = fitted_models()[model.parent]

    def sellmeier(x: Float[np.ndarray, " n"]) -> Float[np.ndarray, " n"]:
        """Return the Sellmeier form's ``n - 1``.

        Parameters
        ----------
        x : numpy.ndarray
            Vacuum wavelengths, m.

        Returns
        -------
        numpy.ndarray
            ``n - 1``.
        """
        x2 = np.asarray(x)[..., None] ** 2
        return constant + np.sum(b * x2 / (x2 - c), axis=-1)

    def full(x: Float[np.ndarray, " n"]) -> Float[np.ndarray, " n"]:
        """Return the full parent's ``n - 1`` at the state.

        Parameters
        ----------
        x : numpy.ndarray
            Vacuum wavelengths, m.

        Returns
        -------
        numpy.ndarray
            ``n - 1``, without range warnings (the stencil reaches past the ends).
        """
        return parent.refractivity(x, state, warn=False)

    points = int(str(model.record["check_points"]))
    grid = np.concatenate([np.geomspace(a, z, points) for a, z in model.intervals])
    exact = _beta2(full, grid)
    return {
        "n_minus_1": float(np.max(np.abs(sellmeier(grid) / full(grid) - 1))),
        "beta2": float(np.max(np.abs(_beta2(sellmeier, grid) - exact)))
        / float(np.max(np.abs(exact))),
    }


def coefficients(
    gas: Gas | str,
    form: str = "polarizability",
    *,
    pressure: float | None = None,
    temperature: float | None = None,
) -> CompactCoefficients:
    """Return a compact model's coefficients for print (OEP 2, section 6).

    Parameters
    ----------
    gas : Gas or str
        CO, CO2, CH4, N2O or SF6.
    form : str, optional
        ``polarizability`` (the poles, with the package's density conversion) or
        ``sellmeier`` (``n - 1`` in the dilute limit at one state).
    pressure : float, optional
        The Sellmeier form's pressure, Pa (within the equation of state's domain).
    temperature : float, optional
        Its temperature, K.

    Returns
    -------
    CompactCoefficients
        The terms in descending energy, the qualified intervals, the reduction's
        residuals and, for the Sellmeier form, its state and stated error.

    Raises
    ------
    OptmatError
        For a gas without a compact model (an atom, H2, D2, N2, O2, O3).
    ValueError
        For another form, or the Sellmeier form without a state.

    Examples
    --------
    >>> table = coefficients("CO2")
    >>> table.kinds[:4], table.kinds[-1]
    (('effective', 'effective', 'line', 'line'), 'rotational')
    >>> air = coefficients("CO2", "sellmeier", pressure=101325.0, temperature=293.15)
    >>> f"{air.number_density:.4e}"
    '2.5169e+25'
    """
    name = f"compact/{as_gas(gas).value}"
    model = compact_models().get(name)
    if model is None:
        msg = (
            f"{as_gas(gas).value} has no compact model: it has no fixed infrared terms"
        )
        raise OptmatError(msg)
    if form not in ("polarizability", "sellmeier"):
        msg = f"form is polarizability or sellmeier, not {form!r}"
        raise ValueError(msg)
    if form == "sellmeier" and (pressure is None or temperature is None):
        msg = "the Sellmeier form needs a pressure and a temperature"
        raise ValueError(msg)
    energies = np.asarray(model.energies)
    order = np.argsort(-energies, kind="stable")
    n_fitted = len(model.strengths) - model.n_fixed
    kinds = np.array(
        ["line"] * model.n_lines
        + ["effective"] * (n_fitted - model.n_lines)
        + ["rotational"] * model.n_rotational
        + ["infrared"] * (model.n_fixed - model.n_rotational)
    )[order]
    f, e = np.asarray(model.strengths)[order], energies[order]
    wavelength = 2 * np.pi * C_LIGHT / (e * AU_ANGULAR_FREQUENCY)
    chi0 = model.offset_au * AU_POLARIZABILITY
    electronic = [int(k) for k in order if k < n_fitted]
    state = density = scale = None
    b = c = np.zeros(0)
    stated: dict[str, float] = {}
    constant = 0.0
    if pressure is not None and temperature is not None and form == "sellmeier":
        state = State(pressure, temperature)
        density = float(number_density(model.gas, state))
        scale = density * AU_POLARIZABILITY / (2 * EPS0)
        b, c, constant = scale * f / e**2, wavelength**2, density * chi0 / (2 * EPS0)
        stated = _stated_error(model, state, b, c, constant)
    covariance, parameters = _printed_covariance(model, electronic, scale)
    return CompactCoefficients(
        model=name,
        parent=model.parent,
        model_digest=model.identity.digest,
        parent_digest=fitted_models()[model.parent].identity.digest,
        form=form,
        labels=tuple(model.labels[k] for k in order),
        kinds=tuple(str(x) for x in kinds),
        wavelengths=tuple(float(x) for x in wavelength),
        strengths=tuple(float(x) for x in f) if state is None else (),
        chi0=chi0,
        B=tuple(float(x) for x in b),
        C=tuple(float(x) for x in c),
        constant=constant,
        pressure=None if state is None else pressure,
        temperature=None if state is None else temperature,
        number_density=density,
        intervals=model.intervals,
        residuals_sigma=model.residuals_sigma,
        stated_error=stated,
        covariance=covariance,
        covariance_parameters=parameters,
    )


def format_coefficients(table: CompactCoefficients) -> str:
    """Return a compact model's coefficients as a Markdown table, to 11 digits.

    Eleven significant digits keep the rounding within a hundredth of the parent's
    sigma (OEP 2, section 7).

    Parameters
    ----------
    table : CompactCoefficients
        The coefficients.

    Returns
    -------
    str
        A caption line and the table, in descending energy; wavelengths in um.
    """
    intervals = ", ".join(f"{a * 1e9:.1f}-{b * 1e9:.1f}" for a, b in table.intervals)
    digits = f".{_DIGITS - 1}e"
    if table.form == "polarizability":
        caption = (
            f"`{table.model}`, the polarizability form; qualified at {intervals} nm. "
            "$\\alpha = \\sum_k f_k / (E_k^2 - \\omega^2)$ in atomic units, "
            "$E_k = hc / \\lambda_k$."
        )
        header = "| Term | Kind | $f_k$ | $\\lambda_k$ (µm) |"
        rows = [
            f"| {label} | {kind} | {f:{digits}} | {lam * 1e6:{digits}} |"
            for label, kind, f, lam in zip(
                table.labels,
                table.kinds,
                table.strengths,
                table.wavelengths,
                strict=True,
            )
        ]
    else:
        caption = (
            f"`{table.model}`, the Sellmeier form at {table.pressure:g} Pa and "
            f"{table.temperature:g} K "
            f"($N$ = {table.number_density:{digits}} m$^{{-3}}$); "
            f"qualified at {intervals} nm. $n - 1 = \\sum_k B_k \\lambda^2 / "
            "(\\lambda^2 - C_k)$; against the full model, at most "
            f"{table.stated_error['n_minus_1']:.2g} relative in $n - 1$ and "
            f"{table.stated_error['beta2']:.2g} of the largest $\\beta_2$."
        )
        header = "| Term | Kind | $B_k$ | $C_k$ (µm²) |"
        rows = [
            f"| {label} | {kind} | {b:{digits}} | {c * 1e12:{digits}} |"
            for label, kind, b, c in zip(
                table.labels, table.kinds, table.B, table.C, strict=True
            )
        ]
    return "\n".join([caption, "", header, "|---|---|---:|---:|", *rows])
