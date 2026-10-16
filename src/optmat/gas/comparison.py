# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The refractiveindex.info database, as a comparison source only.

A pinned snapshot of the database's gas entries (CC0; see
``data/external/refractiveindex_info/PROVENANCE.md``) is converted to
``data/external/refractiveindex_info.json`` by
``tools/convert_refractiveindex_info.py``.
This module reads that file. Its formula entries become models (kind "comparison"),
so the gallery can compare them with the curated data and the fits; its tabulated
entries are returned as arrays. Nothing here enters a fit: every value the database
takes from a paper is curated from the paper itself.

The formula types used by the gas entries, from the database's documentation (``lambda``
in um): type 2, ``n^2 - 1 = C1 + sum C_{2i} lambda^2 / (lambda^2 - C_{2i+1})``; type 5,
``n = C1 + sum C_{2i} lambda^{C_{2i+1}}``; type 6, ``n - 1 = C1 + sum C_{2i} /
(C_{2i+1} - lambda^-2)``.
"""

import functools
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, ClassVar, Final

import numpy as np
from jaxtyping import Complex, Float

from optmat.conditions import ONE_ATM_0C, State
from optmat.constants import C
from optmat.datasets import DATA_DIRECTORY
from optmat.exceptions import DataError
from optmat.gas.eos import number_density
from optmat.gas.gases import Gas, as_gas
from optmat.gas.lorentz_lorenz import (
    dpolarizability_dsusceptibility,
    polarizability_from_refractivity,
    polarizability_from_susceptibility,
)
from optmat.gas.models.base import Model

__all__ = [
    "COMMIT",
    "DatabaseModel",
    "TabulatedEntry",
    "database_models",
    "tabulated_entries",
]

COMMIT: Final = "c5c2f188e848453def5970e347399d653df2ffc2"
"""The refractiveindex.info-database commit the snapshot was taken from."""

_FORMULA_TYPES: Final = frozenset({"formula 2", "formula 5", "formula 6"})


@dataclass(frozen=True)
class DatabaseModel(Model):
    """A refractiveindex.info formula entry as a polarizability model.

    Attributes
    ----------
    formula_type : str
        "formula 2", "formula 5" or "formula 6".
    coefficients : tuple of float
        ``C1, C2, ...`` as stored.
    state : State
        The entry's stated conditions (0 degrees C and 1 atm where it states none).

    Its overrides are the coefficients and the intervals.
    """

    formula_type: str = ""
    coefficients: tuple[float, ...] = field(default=(), repr=False)
    state: State = ONE_ATM_0C

    EVALUATION_VERSION: ClassVar[int] = 1
    INPUTS: ClassVar[tuple[str, ...]] = ("coefficients", "intervals")

    def _polarizability(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Float[np.ndarray, "*shape"]:
        """Evaluate the entry's formula and convert to a polarizability.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            Polarizability, C m^2/V.
        """
        lam = wavelength * 1e6
        c = self.coefficients
        density = self.reference_density
        pairs = list(zip(c[1::2], c[2::2], strict=False))
        if self.formula_type == "formula 2":
            chi = c[0] + sum(b * lam**2 / (lam**2 - p) for b, p in pairs)
            return polarizability_from_susceptibility(
                chi, density, self.refractivity_virial
            )
        if self.formula_type == "formula 5":
            dn = c[0] - 1.0 + sum(b * lam**p for b, p in pairs)
        else:  # formula 6
            dn = c[0] + sum(b / (p - lam**-2.0) for b, p in pairs)
        return polarizability_from_refractivity(dn, density, self.refractivity_virial)

    def _dcomplex_polarizability_domega(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the entry's analytic d alpha / d omega.

        Each formula is differentiated in ``u = 1/lambda`` (um^-1), proportional to
        ``omega``, and reduced through inverse Lorentz-Lorenz's Jacobian at the
        entry's density. A form-2 entry is NaN in the static limit, as its value is.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            C m^2 V^-1 s (real, as complex).
        """
        u = 1e-6 / wavelength
        c = self.coefficients
        pairs = list(zip(c[1::2], c[2::2], strict=False))
        if self.formula_type == "formula 2":
            with np.errstate(divide="ignore", invalid="ignore"):
                lam = 1 / u
                chi = c[0] + sum(b * lam**2 / (lam**2 - p) for b, p in pairs)
            dchi = sum(b * 2 * p * u / (1 - p * u**2) ** 2 for b, p in pairs)
            dchi = np.where(np.isfinite(chi), dchi, np.nan)
        else:
            if self.formula_type == "formula 5":
                dn = c[0] - 1.0 + sum(b * u**-p for b, p in pairs)
                ddn = sum(-p * b * u ** (-p - 1) for b, p in pairs)
            else:  # formula 6
                dn = c[0] + sum(b / (p - u**2) for b, p in pairs)
                ddn = sum(b * 2 * u / (p - u**2) ** 2 for b, p in pairs)
            chi = dn * (dn + 2)
            dchi = 2 * (1 + dn) * ddn
        jacobian = dpolarizability_dsusceptibility(
            chi, self.reference_density, self.refractivity_virial
        )
        return np.asarray(jacobian * dchi * 1e-6 / (2 * math.pi * C), dtype=complex)

    @functools.cached_property
    def reference_density(self) -> float:
        """The number density of the entry's stated conditions, resolved once.

        Returns
        -------
        float
            The density, m^-3.
        """
        return number_density(self.gas, self.state)

    def _identity_extras(self) -> Mapping[str, Any]:
        """Return the values resolved once: the reference density and the virial.

        Returns
        -------
        Mapping
            ``refractivity_virial`` and ``reference_density``.
        """
        return {
            **super()._identity_extras(),
            "reference_density": self.reference_density,
        }


@dataclass(frozen=True)
class TabulatedEntry:
    """A refractiveindex.info tabulated entry.

    Attributes
    ----------
    gas : Gas
        The species.
    name : str
        ``refractiveindex.info/<gas>/<entry>``.
    references : str
        The entry's references, as stored (HTML).
    state : State
        The entry's stated conditions (0 degrees C and 1 atm where it states none).
    wavelength : numpy.ndarray
        Wavelengths, m (the database's wavelengths are vacuum wavelengths).
    refractive_index : numpy.ndarray
        ``n``.
    """

    gas: Gas
    name: str
    references: str
    state: State
    wavelength: Float[np.ndarray, " n"] = field(repr=False)
    refractive_index: Float[np.ndarray, " n"] = field(repr=False)


def _state(conditions: Mapping[str, object]) -> State:
    """Return an entry's state from its conditions.

    Parameters
    ----------
    conditions : mapping
        The entry's CONDITIONS (temperature in K, pressure in Pa), possibly empty.

    Returns
    -------
    State
        The state, completed with 0 degrees C and 1 atm where a value is missing.
    """
    temperature = conditions.get("temperature", ONE_ATM_0C.temperature)
    pressure = conditions.get("pressure", ONE_ATM_0C.pressure)
    return State(float(pressure), float(temperature))  # pyright: ignore[reportArgumentType]


@functools.cache
def _records() -> tuple[Mapping[str, object], ...]:
    """Return the snapshot's records.

    Returns
    -------
    tuple of mapping
        One record per ``nk`` entry.
    """
    path = DATA_DIRECTORY / "external" / "refractiveindex_info.json"
    return tuple(json.loads(path.read_text(encoding="utf-8")))


def _blocks(record: Mapping[str, object]) -> Sequence[Mapping[str, object]]:
    """Return a record's data blocks.

    Parameters
    ----------
    record : mapping
        A snapshot record.

    Returns
    -------
    sequence of mapping
        Its DATA blocks.

    Raises
    ------
    DataError
        If the record has not exactly one data block.
    """
    blocks = record["data"]
    if not isinstance(blocks, Sequence) or len(blocks) != 1:
        msg = f"refractiveindex.info {record['file']}: expected one data block"
        raise DataError(msg)
    return blocks  # pyright: ignore[reportReturnType]


@functools.cache
def database_models() -> dict[str, DatabaseModel]:
    """Return the snapshot's formula entries as models, by name.

    Returns
    -------
    dict of str to DatabaseModel
        Models named ``refractiveindex.info/<gas>/<entry>``.

    Examples
    --------
    >>> model = database_models()["refractiveindex.info/Ar/Bideau-Mehu"]
    >>> f"{float(model.refractivity(253.7e-9)) * 1e6:.1f}"
    '303.2'
    """
    models = {}
    for record in _records():
        (block,) = _blocks(record)
        if block["type"] not in _FORMULA_TYPES:
            continue
        gas = as_gas(str(record["gas"]))
        name = f"refractiveindex.info/{gas.value}/{record['entry']}"
        low, high = (float(x) * 1e-6 for x in block["wavelength_range_um"])  # pyright: ignore[reportGeneralTypeIssues]
        models[name] = DatabaseModel(
            gas=gas,
            name=name,
            citation=str(record["references"]),
            wavelength_range=(low, high),
            kind="comparison",
            status="comparison",
            status_reason="a refractiveindex.info entry, kept for comparison",
            formula_type=str(block["type"]),
            coefficients=tuple(float(x) for x in block["coefficients"]),  # pyright: ignore[reportGeneralTypeIssues]
            state=_state(record["conditions"]),  # pyright: ignore[reportArgumentType]
        )
    return models


@functools.cache
def tabulated_entries() -> tuple[TabulatedEntry, ...]:
    """Return the snapshot's tabulated entries.

    Returns
    -------
    tuple of TabulatedEntry
        The entries, each with its wavelengths and refractive indices.

    Examples
    --------
    >>> {entry.gas.value for entry in tabulated_entries()} >= {"H2", "Kr"}
    True
    """
    entries = []
    for record in _records():
        (block,) = _blocks(record)
        if block["type"] != "tabulated n":
            continue
        rows = np.array(block["rows"], dtype=float)
        gas = as_gas(str(record["gas"]))
        entries.append(
            TabulatedEntry(
                gas=gas,
                name=f"refractiveindex.info/{gas.value}/{record['entry']}",
                references=str(record["references"]),
                state=_state(record["conditions"]),  # pyright: ignore[reportArgumentType]
                wavelength=rows[:, 0] * 1e-6,
                refractive_index=rows[:, 1],
            )
        )
    return tuple(entries)
