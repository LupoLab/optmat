# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Fused silica's Raman responses against measurements (TC-693, TC-694).

As :mod:`optmat.gas.raman_response_validation`, for :mod:`optmat.solid.raman_response`:

- :func:`silica_gains`: the absolute Raman gain of fused quartz (Stolen & Ippen 1973),
  for each model with its own f_R and optmat's long-pulse n2;
- :func:`silica_low_frequency_gain`: the shape of the gain below 25 cm^-1 (Stolen et
  al. 1989's Eq. (10));
- :func:`silica_raman_fractions`: the nuclear share of the long-pulse n2 from
  spontaneous Raman spectra (Hellwarth, Cherlow & Yang 1975), against each model's
  f_R.
"""

import math
from collections.abc import Mapping

import numpy as np

from optmat._catalogue import records as sources_of
from optmat.constants import C
from optmat.exceptions import DataError
from optmat.gas.raman_response_validation import ResponseComparison, dataset, quiet
from optmat.solid.raman_response import silica_raman_responses

__all__ = ["silica_gains", "silica_low_frequency_gain", "silica_raman_fractions"]

_WAVENUMBER = 2 * math.pi * C * 100
"""Rad/s per cm^-1."""


def silica_low_frequency_gain() -> list[ResponseComparison]:
    """Return silica's gain below 25 cm^-1 against Stolen et al. 1989's Eq. (10).

    The gain normalised to its maximum (near 440 cm^-1),
    ``a Omega + b Omega^3`` with Omega in cm^-1 (``stolen_1989/low_frequency_gain``).

    Returns
    -------
    list of ResponseComparison
        Per model and shift (dimensionless).
    """
    (formula,) = (
        f
        for f in sources_of()["stolen_1989"].formulas
        if f["id"] == "low_frequency_gain"
    )
    coefficients = formula["coefficients"]
    if not isinstance(coefficients, Mapping):
        msg = "stolen_1989/low_frequency_gain has no coefficients"
        raise DataError(msg)
    a = float(str(coefficients["a"]))  # pyright: ignore[reportUnknownArgumentType]
    b = float(str(coefficients["b"]))  # pyright: ignore[reportUnknownArgumentType]
    grid = np.linspace(1.0, 1300.0, 13000) * _WAVENUMBER
    out = []
    for name, model in silica_raman_responses().items():
        peak = float(np.max(np.imag(model.spectrum(grid))))
        for shift in (5.0, 10.0, 20.0):
            out.append(
                ResponseComparison(
                    material="SiO2",
                    reference="stolen_1989/low_frequency_gain",
                    label=f"{shift:g} cm-1",
                    quantity="gain / peak gain",
                    unit="1",
                    expected=a * shift + b * shift**3,
                    sigma=math.nan,
                    predicted=float(np.imag(model.spectrum(shift * _WAVENUMBER)))
                    / peak,
                    model=name,
                )
            )
    return out


def silica_gains() -> list[ResponseComparison]:
    """Return fused quartz's measured Raman gain against each silica model's.

    Stolen & Ippen 1973: pump 526.0 nm, Stokes shift 330 cm^-1; each model with its
    own f_R and optmat's long-pulse n2 at the pump (the Kerr fit, TC-393).

    Returns
    -------
    list of ResponseComparison
        Per model (cm/W).
    """
    record = dataset("stolen_1973/gain_amplifier")
    pump = float(str(record.parameters["pump_nm"])) * 1e-9
    shift = float(str(record.parameters["shift_cm1"]))
    out = []
    for name, model in silica_raman_responses().items():
        out.append(
            ResponseComparison(
                material="SiO2",
                reference=record.key,
                label=f"{pump * 1e9:.1f} nm, {shift:g} cm-1",
                quantity="gain",
                unit="cm/W",
                expected=float(record.scaled_values()[0]),
                sigma=float(record.standard_uncertainty()[0]),
                predicted=float(quiet(model.gain, pump, shift * _WAVENUMBER)) * 100,
                model=name,
            )
        )
    return out


def silica_raman_fractions() -> list[ResponseComparison]:
    """Return the nuclear share of silica's long-pulse n2 against each model's f_R.

    Returns
    -------
    list of ResponseComparison
        Per model (dimensionless).
    """
    record = dataset("hellwarth_1975/nuclear_fraction_SiO2")
    return [
        ResponseComparison(
            material="SiO2",
            reference=record.key,
            label="fused quartz, linear polarisation",
            quantity="Raman share of n2",
            unit="1",
            expected=float(record.scaled_values()[0]),
            sigma=float(record.standard_uncertainty()[0]),
            predicted=model.raman_fraction,
            model=name,
        )
        for name, model in silica_raman_responses().items()
    ]
