# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The derivation of the two-photon ionisation models (OEP 6, Phase 1).

A model (:class:`~optmat.gas.two_photon.TwoPhotonIonisationModel`) is derived from:

- ``data/gas/two_photon.toml``: per gas, the NIST datasets of its thresholds and of
  :math:`E_1`; per polarisation, the approved reference (TC-943), its validators, its
  uncertainty and the width of the resonance wing below :math:`E_1`;
- the curated datasets those name.

The rules, all stated in the manual (*Two-photon ionisation and absorption*):

- **Thresholds:** half the neutral's ionisation energy, and half of it plus the ion's
  :math:`{}^2P_{1/2}` level; :math:`E_1` is the level, one photon.
- **The reference's vertices:** each dataset's drawn vertices, in photon angular
  frequency and m^4 s. Where a vector drawing's staircase puts several vertices at one
  frequency (the drawing's quantisation), they are one point at their geometric mean.
  Partial cross sections are each interpolated, linearly in log sigma, on the union of
  their vertices over their common coverage and summed. Only vertices from the lowest
  threshold up to :math:`E_1` are kept.
- **Regions:**
  - ``below_threshold`` up to the lowest threshold;
  - ``unavailable`` from there to the reference's first vertex, and from its last
    vertex (or :math:`E_1`, whichever is lower) up;
  - ``unsupported_detail``, for an LS-coupled reference, between the ion's two
    thresholds, and the resonance wing, ``wing_ev`` below :math:`E_1`;
  - ``supported`` elsewhere on the reference.
- **Spreads:** each validator interpolated the same way at the reference's vertices in
  the supported region that it covers, as a ratio to the reference.

``tools/fit/two_photon_record.py`` writes the stored models, and a test checks they are
fresh.
"""

import math
from collections.abc import Mapping, Sequence
from typing import Any, Final

import numpy as np

from optmat.constants import HBAR
from optmat.datasets import Dataset, datasets
from optmat.exceptions import DataError
from optmat.gas.gases import as_gas
from optmat.gas.ionisation import _ENERGY_UNITS
from optmat.gas.two_photon import (
    DERIVATION_VERSION,
    TwoPhotonIonisationModel,
    TwoPhotonReference,
    TwoPhotonRegion,
    TwoPhotonSpread,
    TwoPhotonThreshold,
    TwoPhotonUncertainty,
    _selection,
)
from optmat.identity import parse_status, sha256_digest

__all__ = ["derive", "summed_curve", "vertices"]

_TO_M4S: Final = {"cm4 s": 1e-8}
"""The cross-section units a reference may carry, to m^4 s."""
_EV: Final = _ENERGY_UNITS["eV"]


def _datasets(keys: Sequence[str]) -> list[Dataset]:
    """Return the curated datasets with these keys.

    Parameters
    ----------
    keys : sequence of str
        ``source/dataset`` keys.

    Returns
    -------
    list of Dataset
        The datasets, in the order given.

    Raises
    ------
    DataError
        If one does not exist.
    """
    by_key = {d.key: d for d in datasets()}
    missing = [k for k in keys if k not in by_key]
    if missing:
        msg = f"two_photon.toml names unknown datasets {missing}"
        raise DataError(msg)
    return [by_key[k] for k in keys]


def vertices(dataset: Dataset) -> tuple[np.ndarray, np.ndarray]:
    """Return a cross-section dataset's vertices: angular frequency and m^4 s.

    Vertices that share a frequency (a vector drawing's vertical step) become one, at
    the geometric mean of their values.

    Parameters
    ----------
    dataset : Dataset
        A two-photon ``multiphoton_ionization_cross_section`` dataset in cm4 s.

    Returns
    -------
    tuple of numpy.ndarray
        Strictly increasing angular frequencies, rad/s, and the cross sections, m^4 s.

    Raises
    ------
    DataError
        If it is not a two-photon cross section in a known unit, or a value is not
        positive.
    """
    if (
        dataset.quantity != "multiphoton_ionization_cross_section"
        or str(dataset.parameters.get("photons")) != "2"
        or dataset.unit not in _TO_M4S
    ):
        msg = f"{dataset.key} is not a two-photon cross section in {sorted(_TO_M4S)}"
        raise DataError(msg)
    omega = dataset.angular_frequency()
    sigma = dataset.scaled_values() * _TO_M4S[dataset.unit]
    if np.any(sigma <= 0):
        msg = f"{dataset.key}: a cross section is not positive"
        raise DataError(msg)
    order = np.argsort(omega, kind="stable")
    omega, log_sigma = omega[order], np.log(sigma[order])
    unique, start = np.unique(omega, return_index=True)
    sums = np.add.reduceat(log_sigma, start)
    counts = np.diff(np.append(start, len(omega)))
    return unique, np.exp(sums / counts)


def _interpolate(omega: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return a curve interpolated linearly in log y; NaN outside its vertices.

    Parameters
    ----------
    omega : numpy.ndarray
        Where, rad/s.
    x : numpy.ndarray
        The curve's frequencies, increasing.
    y : numpy.ndarray
        Its positive values.

    Returns
    -------
    numpy.ndarray
        The values, NaN outside ``[x[0], x[-1]]``.
    """
    out = np.exp(np.interp(omega, x, np.log(y)))
    out[(omega < x[0]) | (omega > x[-1])] = np.nan
    return out


def summed_curve(keys: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    """Return a calculation's total: its datasets summed on their common coverage.

    Parameters
    ----------
    keys : sequence of str
        The datasets: one, or the partial cross sections of one calculation.

    Returns
    -------
    tuple of numpy.ndarray
        Angular frequencies, rad/s, and the summed cross section, m^4 s.

    Raises
    ------
    DataError
        If a dataset is unknown or not a cross section, or the partials do not
        overlap.
    """
    curves = [vertices(d) for d in _datasets(keys)]
    if len(curves) == 1:
        return curves[0]
    low = max(x[0] for x, _ in curves)
    high = min(x[-1] for x, _ in curves)
    grid = np.unique(np.concatenate([x for x, _ in curves]))
    grid = grid[(grid >= low) & (grid <= high)]
    if len(grid) < 2:
        msg = f"{list(keys)}: the partial cross sections do not overlap"
        raise DataError(msg)
    return grid, np.sum([_interpolate(grid, x, y) for x, y in curves], axis=0)


def _level(key: str) -> float:
    """Return a NIST level or ionisation energy, J.

    Parameters
    ----------
    key : str
        A one-row ``level_energy`` or ``ionization_energy`` dataset.

    Returns
    -------
    float
        Its energy, J.
    """
    (dataset,) = _datasets([key])
    return float(dataset.scaled_values()[0]) * _ENERGY_UNITS[dataset.unit]


def thresholds(gas: str) -> tuple[tuple[TwoPhotonThreshold, ...], TwoPhotonThreshold]:
    """Return a gas's two-photon thresholds and its first one-photon resonance.

    Parameters
    ----------
    gas : str
        The gas.

    Returns
    -------
    tuple
        The thresholds (lowest first) and :math:`E_1`.
    """
    entry = _selection()[gas]
    keys = entry["thresholds"]
    ionisation = _level(keys[0])
    ion = f"{gas}+"
    names = ["He+ 2S1/2"] if gas == "He" else [f"{ion} 2P3/2", f"{ion} 2P1/2"]
    energies = [ionisation, *(ionisation + _level(k) for k in keys[1:])]
    chain = [[keys[0]], *([keys[0], k] for k in keys[1:])]
    found = tuple(
        TwoPhotonThreshold(name, energy / (2 * HBAR), 2, energy, tuple(datasets_))
        for name, energy, datasets_ in zip(names, energies, chain, strict=True)
    )
    e1 = _level(entry["e1"])
    return found, TwoPhotonThreshold("E1", e1 / HBAR, 1, e1, (entry["e1"],))


def _regions(
    first: float,
    last: float,
    found: Sequence[TwoPhotonThreshold],
    e1: TwoPhotonThreshold,
    choice: Mapping[str, Any],
) -> tuple[TwoPhotonRegion, ...]:
    """Return the regions of a model.

    Parameters
    ----------
    first, last : float
        The reference's first and last vertices, rad/s.
    found : sequence of TwoPhotonThreshold
        The thresholds.
    e1 : TwoPhotonThreshold
        The first one-photon resonance.
    choice : Mapping
        The polarisation's selection.

    Returns
    -------
    tuple of TwoPhotonRegion
        Bands tiling [0, inf).
    """
    lower = found[0].omega
    end = min(np.nextafter(last, math.inf), e1.omega)
    wing = e1.omega - float(choice["wing_ev"]) * _EV / HBAR
    detail: list[tuple[float, float, str]] = []
    if not choice["resolves_fine_structure"] and len(found) > 1:
        detail.append(
            (
                max(first, lower),
                found[1].omega,
                "An LS-coupled reference has one threshold: between the ion's two "
                "fine-structure thresholds its values do not resolve the channels "
                "or their autoionising series.",
            )
        )
    detail.append(
        (
            max(wing, first),
            end,
            f"Within {choice['wing_ev']} eV of E1: the wing of the first one-photon "
            "resonance, which the calculations place at their own energies, not "
            "NIST's.",
        )
    )
    regions = [
        TwoPhotonRegion(
            "below_threshold",
            0.0,
            lower,
            "Below the lowest "
            "two-photon threshold: no direct continuum ionisation. Bound two-photon "
            "excitation and (2+1) ionisation are not modelled (OEP 6, Phase 2).",
        )
    ]
    if first > lower:
        regions.append(
            TwoPhotonRegion(
                "unavailable",
                lower,
                first,
                "Above the threshold but below the reference's first vertex: no value.",
            )
        )
    cursor = first
    for start, stop, reason in sorted(detail):
        if start > cursor:
            regions.append(
                TwoPhotonRegion("supported", cursor, start, "The reference's values.")
            )
        regions.append(TwoPhotonRegion("unsupported_detail", start, stop, reason))
        cursor = stop
    if end < e1.omega:
        regions.append(
            TwoPhotonRegion(
                "unavailable",
                end,
                e1.omega,
                "Above the reference's last vertex, below E1: no value.",
            )
        )
    regions.append(
        TwoPhotonRegion(
            "unavailable",
            e1.omega,
            math.inf,
            "At and above E1, the first one-photon resonance: not qualified by this "
            "model (OEP 6).",
        )
    )
    return tuple(regions)


def _spreads(
    omega: np.ndarray,
    sigma: np.ndarray,
    regions: Sequence[TwoPhotonRegion],
    validators: Sequence[Mapping[str, Any]],
) -> tuple[TwoPhotonSpread, ...]:
    """Return each validator against the reference over the supported region.

    Parameters
    ----------
    omega, sigma : numpy.ndarray
        The reference.
    regions : sequence of TwoPhotonRegion
        Its regions.
    validators : sequence of Mapping
        The selection's validators.

    Returns
    -------
    tuple of TwoPhotonSpread
        One per validator.

    Raises
    ------
    DataError
        If a validator does not overlap the supported region.
    """
    supported = np.zeros(omega.shape, dtype=bool)
    for r in regions:
        if r.label == "supported":
            supported |= (omega >= r.lower) & (omega < r.upper)
    spreads = []
    for v in validators:
        x, y = summed_curve(v["datasets"])
        ratio = _interpolate(omega, x, y) / sigma
        used = supported & np.isfinite(ratio)
        if not np.any(used):
            msg = f"{v['label']}: no overlap with the supported region"
            raise DataError(msg)
        spreads.append(
            TwoPhotonSpread(
                label=str(v["label"]),
                datasets=tuple(v["datasets"]),
                form=str(v["form"]),
                role=str(v["role"]),
                lower=float(omega[used][0]),
                upper=float(omega[used][-1]),
                points=int(np.count_nonzero(used)),
                lower_ratio=float(ratio[used].min()),
                upper_ratio=float(ratio[used].max()),
            )
        )
    return tuple(spreads)


def _digest(dataset: Dataset) -> str:
    """Return a dataset's digest: what it holds and how it is to be read.

    Parameters
    ----------
    dataset : Dataset
        The dataset.

    Returns
    -------
    str
        SHA-256, 64 hexadecimal digits.
    """
    return sha256_digest(
        [
            dataset.key,
            dataset.unit,
            dataset.scale,
            dataset.abscissa_unit,
            dataset.wavelength_medium,
            list(dataset.abscissae),
            list(dataset.values),
        ]
    )


def derive(gas: str, polarization: str) -> TwoPhotonIonisationModel:
    """Derive a gas's two-photon ionisation model from its selection and data.

    Parameters
    ----------
    gas : str
        The gas.
    polarization : str
        ``"linear"`` or ``"circular"``.

    Returns
    -------
    TwoPhotonIonisationModel
        The model, as it is to be stored.

    Raises
    ------
    DataError
        If the selection names no reference for it, or its data are unusable.
    """
    species = as_gas(gas)
    entry = _selection().get(species.value, {})
    if polarization not in entry:
        msg = f"two_photon.toml has no {polarization} reference for {species.value}"
        raise DataError(msg)
    choice = entry[polarization]
    found, e1 = thresholds(species.value)
    x, y = summed_curve(choice["reference"])
    keep = (x >= found[0].omega) & (x <= e1.omega)
    omega, sigma = x[keep], y[keep]
    regions = _regions(float(omega[0]), float(omega[-1]), found, e1, choice)
    spreads = _spreads(omega, sigma, regions, choice["validators"])
    source = choice["source_uncertainty"]
    keys = {
        *choice["reference"],
        *(k for v in choice["validators"] for k in v["datasets"]),
        *(k for t in (*found, e1) for k in t.datasets),
    }
    status, reason = parse_status(
        choice, f"{species.value}.{polarization}", default="validated"
    )
    return TwoPhotonIonisationModel(
        gas=species,
        polarization=polarization,
        derivation_version=DERIVATION_VERSION,
        omega=omega,
        sigma=sigma,
        thresholds=found,
        e1=e1,
        regions=regions,
        reference=TwoPhotonReference(
            datasets=tuple(choice["reference"]),
            gauge=str(choice["gauge"]),
            method=str(choice["method"]),
            intensity_w_cm2=(
                float(choice["intensity_w_cm2"])
                if "intensity_w_cm2" in choice
                else None
            ),
            reason=" ".join(str(choice["reason"]).split()),
            registers=tuple(choice["registers"]),
        ),
        uncertainty=TwoPhotonUncertainty(
            source=None if isinstance(source, str) else float(source),
            source_level=str(choice["source_level"]),
            source_provenance=str(choice["source_provenance"]),
            digitising=float(choice["digitising"]),
            digitising_provenance=str(choice["digitising_provenance"]),
            spreads=spreads,
        ),
        input_digests={"selection": sha256_digest(entry)}
        | {d.key: _digest(d) for d in _datasets(sorted(keys))},
        status=status,
        status_reason=reason,
    )
