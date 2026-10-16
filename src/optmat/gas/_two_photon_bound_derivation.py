# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The derivation of the bound two-photon level records (OEP 6, Phase 2a).

A record (:class:`~optmat.gas.two_photon_bound.TwoPhotonBoundLevels`) is derived from:

- the NIST level list held as returned,
  ``data/external/nist_asd/<gas>_I_levels.tsv`` (``nist_asd_2026``): every even level
  (its term not starred) with :math:`J = 0, 1, 2` below the ionisation energy
  (``nist_asd_2026/ie_<gas>_I``). For helium, singlets only: two photons from the
  1 1S ground level do not reach a triplet;
- ``data/gas/two_photon_bound.toml``: each curated strength and photoionisation cross
  section, with the level it belongs to (TC-952).

Strengths are put in one convention (TC-952):

- ``per_level``: :math:`\delta^{TP}` as printed;
- ``one_component``: a J = 2 level's single tensor component times 4/3;
- ``sigma0``: a frequency-integrated cross section as printed, in m^4, and
  :math:`\delta^{TP}` derived from it.

Each assigned strength's :math:`\sigma_0^{(2)}` is evaluated at its level's NIST photon
frequency. ``tools/fit/two_photon_bound_record.py`` writes the records, and a test
checks that they are fresh.
"""

import csv
import functools
import tomllib
from collections.abc import Mapping
from typing import Any, Final

import numpy as np

from optmat.constants import HBAR
from optmat.datasets import DATA_DIRECTORY, Dataset, datasets, sources
from optmat.exceptions import DataError
from optmat.gas.gases import as_gas
from optmat.gas.ionisation import _ENERGY_UNITS
from optmat.gas.two_photon_bound import (
    _SELECTION,
    BOUND_DERIVATION_VERSION,
    SIGMA0_PER_DELTA,
    BoundLevel,
    BoundPhotoionisation,
    BoundStrength,
    TwoPhotonBoundLevels,
    sigma0_from_delta,
)
from optmat.identity import sha256_digest

__all__ = ["derive_bound", "nist_levels"]

_LEVELS: Final = DATA_DIRECTORY / "external" / "nist_asd"
_CM: Final = _ENERGY_UNITS["cm-1"]
_TO_M4: Final = {"cm4": 1e-8, "m4": 1.0}
_TO_M2: Final = {"cm2": 1e-4, "Mb": 1e-22}
_PRINTED: Final = "[]()?+x "
"""NIST's marks around a level: brackets for derived values, ? for doubt, +x for
an unknown offset; stripped to read the number, kept in ``printed``."""


@functools.cache
def _selection() -> dict[str, Any]:
    """Return the selection file, parsed.

    Returns
    -------
    dict
        Its tables, by gas.
    """
    with _SELECTION.open("rb") as f:
        return tomllib.load(f)


@functools.cache
def _by_key() -> dict[str, Dataset]:
    """Return the curated datasets, by key.

    Returns
    -------
    dict of str to Dataset
        Every dataset.
    """
    return {d.key: d for d in datasets()}


def _dataset(key: str) -> Dataset:
    """Return a curated dataset.

    Parameters
    ----------
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
    if key not in _by_key():
        msg = f"two_photon_bound.toml names an unknown dataset {key!r}"
        raise DataError(msg)
    return _by_key()[key]


def _ionisation(gas: str) -> float:
    """Return a gas's ionisation energy, J.

    Parameters
    ----------
    gas : str
        The gas.

    Returns
    -------
    float
        From ``nist_asd_2026/ie_<gas>_I``.
    """
    d = _dataset(f"nist_asd_2026/ie_{gas}_I")
    return float(d.scaled_values()[0]) * _ENERGY_UNITS[d.unit]


def nist_levels(gas: str) -> tuple[BoundLevel, ...]:
    """Return a gas's even J = 0, 1, 2 levels below its ionisation energy.

    Parameters
    ----------
    gas : str
        The gas.

    Returns
    -------
    tuple of BoundLevel
        By energy (helium's singlets only).
    """
    limit = _ionisation(gas)
    path = _LEVELS / f"{gas}_I_levels.tsv"
    found = []
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            term, j = row["Term"], row["J"]
            printed = row["Level (cm-1)"]
            number = printed.strip(_PRINTED)
            if "*" in term or j not in {"0", "1", "2"} or not number:
                continue
            if gas == "He" and not term.startswith("1"):
                continue
            energy = float(number) * _CM
            if 0 < energy < limit:
                found.append(
                    BoundLevel(
                        configuration=row["Configuration"],
                        term=term,
                        j=int(j),
                        energy=energy,
                        printed=printed,
                        reach="other" if j == "1" else "linear",
                    )
                )
    return tuple(sorted(found, key=lambda level: level.energy))


def _row(dataset: Dataset, row: str) -> float:
    """Return one row's value, after scale.

    Parameters
    ----------
    dataset : Dataset
        The dataset.
    row : str
        The row id.

    Returns
    -------
    float
        Its value in the dataset's unit.

    Raises
    ------
    DataError
        If the row does not exist.
    """
    if row not in dataset.rows:
        msg = f"{dataset.key} has no row {row!r}"
        raise DataError(msg)
    return float(dataset.scaled_values()[dataset.rows.index(row)])


def _strength(
    entry: Mapping[str, Any], levels: Mapping[str, BoundLevel]
) -> BoundStrength:
    """Return one strength in the canonical convention.

    Parameters
    ----------
    entry : Mapping
        Its selection entry.
    levels : Mapping of str to BoundLevel
        The gas's levels, by key.

    Returns
    -------
    BoundStrength
        The strength.

    Raises
    ------
    DataError
        For an unknown level or convention, or a unit that is not the convention's.
    """
    dataset = _dataset(entry["dataset"])
    key = entry.get("level", "")
    if key and key not in levels:
        msg = f"{dataset.key}: no listed level {key!r}"
        raise DataError(msg)
    value = _row(dataset, entry["row"])
    convention = entry["convention"]
    omega = levels[key].energy / (2 * HBAR) if key else None
    if convention in {"per_level", "one_component"} and dataset.unit == "au":
        delta = value * (4 / 3 if convention == "one_component" else 1.0)
        sigma0 = None if omega is None else sigma0_from_delta(delta, omega)
    elif convention == "sigma0" and dataset.unit in _TO_M4:
        sigma0 = value * _TO_M4[dataset.unit]
        delta = None if omega is None else sigma0 / sigma0_from_delta(1.0, omega)
    else:
        msg = f"{dataset.key}: convention {convention!r} with unit {dataset.unit!r}"
        raise DataError(msg)
    return BoundStrength(
        level=key,
        unassigned=str(entry.get("unassigned", "")),
        dataset=dataset.key,
        row=entry["row"],
        source=dataset.source,
        kind="measured"
        if sources()[dataset.source].kind == "measurement"
        else "calculated",
        delta=delta,
        sigma0=sigma0,
        convention=convention,
        derived=convention != "per_level",
        reason=" ".join(str(entry.get("reason", "")).split()),
    )


def _photoionisation(
    entry: Mapping[str, Any], levels: Mapping[str, BoundLevel]
) -> BoundPhotoionisation:
    """Return one photoionisation cross section in SI.

    Parameters
    ----------
    entry : Mapping
        Its selection entry.
    levels : Mapping of str to BoundLevel
        The gas's levels, by key.

    Returns
    -------
    BoundPhotoionisation
        The cross sections, rad/s and m^2.

    Raises
    ------
    DataError
        For an unknown level, or a dataset that is not such a cross section.
    """
    dataset = _dataset(entry["dataset"])
    key = entry.get("level", "")
    if key and key not in levels:
        msg = f"{dataset.key}: no listed level {key!r}"
        raise DataError(msg)
    if (
        dataset.quantity != "excited_state_photoionization_cross_section"
        or dataset.unit not in _TO_M2
    ):
        msg = f"{dataset.key} is not an excited level's photoionisation cross section"
        raise DataError(msg)
    omega = dataset.angular_frequency()
    sigma = dataset.scaled_values() * _TO_M2[dataset.unit]
    return BoundPhotoionisation(
        level=key,
        unassigned=str(entry.get("unassigned", "")),
        dataset=dataset.key,
        source=dataset.source,
        omega=tuple(None if np.isnan(w) else float(w) for w in omega),
        sigma=tuple(float(s) for s in sigma),
        use=dataset.use,
    )


def _digest(dataset: Dataset) -> str:
    """Return a dataset's digest.

    Parameters
    ----------
    dataset : Dataset
        The dataset.

    Returns
    -------
    str
        SHA-256 of its key, unit, scale, abscissae and values.
    """
    return sha256_digest(
        [
            dataset.key,
            dataset.unit,
            dataset.scale,
            dataset.abscissa_unit,
            list(dataset.abscissae),
            list(dataset.values),
        ]
    )


def derive_bound(gas: str) -> TwoPhotonBoundLevels:
    """Derive a gas's bound two-photon level record.

    Parameters
    ----------
    gas : str
        He, Ne, Ar, Kr or Xe.

    Returns
    -------
    TwoPhotonBoundLevels
        The record, as it is to be stored.
    """
    species = as_gas(gas)
    entry = _selection()[species.value]
    levels = nist_levels(species.value)
    by_key = {level.key: level for level in levels}
    strengths = tuple(_strength(e, by_key) for e in entry.get("strengths", []))
    photo = tuple(_photoionisation(e, by_key) for e in entry.get("photoionisation", []))
    keys = sorted(
        {s.dataset for s in (*strengths, *photo)} | {f"nist_asd_2026/ie_{gas}_I"}
    )
    path = _LEVELS / f"{species.value}_I_levels.tsv"
    return TwoPhotonBoundLevels(
        gas=species,
        derivation_version=BOUND_DERIVATION_VERSION,
        ionisation=_ionisation(species.value),
        levels=levels,
        strengths=strengths,
        photoionisation=photo,
        assessment=str(entry["assessment"]),
        assessment_reason=" ".join(str(entry["assessment_reason"]).split()),
        input_digests={
            "selection": sha256_digest(entry),
            "levels": sha256_digest(path.read_text(encoding="utf-8")),
            "sigma0_per_delta": sha256_digest(SIGMA0_PER_DELTA),
        }
        | {k: _digest(_dataset(k)) for k in keys},
        status="validated",
    )
