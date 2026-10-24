# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The gases' Raman responses against measurements (TC-690).

Each function compares one kind of measurement with the responses of
:mod:`optmat.gas.raman_response` and returns :class:`ResponseComparison` records, which
the manual page *The Raman response* tabulates (fused silica's are in
:mod:`optmat.solid.raman_response_validation`). Only published measurements are
references, and none of them entered the parameters of the response it tests.

- :func:`gas_gains`: steady-state Stokes gains (H2, D2, CH4), from the imaginary part
  of the response's transform at the line, so that the neighbouring lines and the
  anisotropic share count;
- :func:`rotational_dephasing`: the decay of the rotational revivals (O2).
"""

import math
import warnings
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from optmat import _catalogue
from optmat.constants import EPS0, K_B, C
from optmat.datasets import Dataset
from optmat.exceptions import ModelRangeWarning
from optmat.gas.linewidths import AMAGAT
from optmat.gas.raman import raman_parameters
from optmat.gas.raman_response import GasRamanResponse, raman_response

__all__ = [
    "GAIN_MEASUREMENTS",
    "ResponseComparison",
    "dataset",
    "gas_gains",
    "quiet",
    "rotational_dephasing",
]

_WAVENUMBER = 2 * math.pi * C * 100
"""Rad/s per cm^-1."""

GAIN_MEASUREMENTS = (
    ("ottusch_1988/gain_H2", "H2", 1, 532.0, 60.0),
    ("bischel_1986b/gain_H2", "H2", 1, None, 20.0),
    ("ottusch_1988/gain_D2", "D2", 2, 532.0, 60.0),
    ("ottusch_1988/gain_CH4", "CH4", None, 532.0, 115.0),
)
"""The measured gains (TC-477, TC-478, as on the page *Raman lines*): the dataset, the
gas, the line's J (None for an unresolved band), the pump wavelength (nm; None where the
dataset's abscissa is the pump) and the density (amagat). Bischel & Dyer's are their
high-density limit, evaluated at 20 amagat, where it holds to 1 % (their Eq. (8))."""


@dataclass(frozen=True)
class ResponseComparison:
    """A measured quantity against a response's prediction.

    Attributes
    ----------
    material : str
        The gas, mixture or solid.
    reference : str
        The curated dataset (``source/dataset``) or formula.
    label : str
        The condition compared (a line, a pump, a shift).
    quantity : str
        What is compared.
    unit : str
        The unit of ``expected`` and ``predicted``.
    expected : float
        The measured (or source's) value.
    sigma : float
        Its stated standard uncertainty (NaN where none is stated).
    predicted : float
        The response's value.
    model : str
        The response's model.
    """

    material: str
    reference: str
    label: str
    quantity: str
    unit: str
    expected: float
    sigma: float
    predicted: float
    model: str = "default"

    @property
    def ratio(self) -> float:
        """The prediction over the measurement.

        Returns
        -------
        float
            ``predicted / expected``.
        """
        return self.predicted / self.expected


def quiet[**P, R](function: Callable[P, R], *args: P.args, **kwargs: P.kwargs) -> R:
    """Call a function with the range warnings silenced (the validators' conditions).

    Parameters
    ----------
    function : callable
        The function.
    *args, **kwargs
        Its arguments.

    Returns
    -------
    R
        Its value.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        return function(*args, **kwargs)


def dataset(key: str) -> Dataset:
    """Return a curated dataset.

    Parameters
    ----------
    key : str
        ``source/dataset``.

    Returns
    -------
    Dataset
        The dataset.
    """
    return _catalogue.dataset(key)


def _response_gain(
    response: GasRamanResponse,
    pump_nm: float,
    amagat: float,
    line: np.ndarray,
) -> float:
    """Return the steady-state gain at a line's centre from the response's transform.

    ``g = omega_s N Im H(Omega) / (c^2 eps0)``, pump and Stokes polarised alike (the
    gas's index taken as 1): :func:`optmat.gas.raman.stokes_gain` for one isolated
    line, with every other line's tail added.

    Parameters
    ----------
    response : GasRamanResponse
        The response, built at the pump wavelength.
    pump_nm : float
        The pump's vacuum wavelength, nm.
    amagat : float
        The density, amagat.
    line : numpy.ndarray
        A mask selecting the line (its strongest member is taken).

    Returns
    -------
    float
        The gain, cm/GW.
    """
    density = amagat * AMAGAT
    omega = quiet(response.frequencies, density)
    index = int(np.flatnonzero(line)[np.argmax(response.amplitude[line])])
    shift = float(omega[index])
    im = float(np.imag(quiet(response.spectrum, shift, density)))
    omega_s = 2 * math.pi * C / (pump_nm * 1e-9) - shift
    return omega_s * density * im / (C**2 * EPS0) * 1e11  # m/W to cm/GW


def gas_gains() -> list[ResponseComparison]:
    """Return the measured Stokes gains against the responses'.

    Returns
    -------
    list of ResponseComparison
        One per measured gain (cm/GW).
    """
    out = []
    for key, gas, j, pump, amagat in GAIN_MEASUREMENTS:
        record = dataset(key)
        factor = 1e9 if record.unit == "cm/W" else 1.0
        values = record.scaled_values() * factor
        sigmas = record.standard_uncertainty() * factor
        pumps = [pump] if pump is not None else list(record.abscissa_values())
        for value, sigma, pump_nm in zip(values, sigmas, pumps, strict=False):
            response = raman_response(gas, 295.0, wavelength=pump_nm * 1e-9)
            line = np.array(response.branch) == "Q"
            if j is not None:
                line &= response.j == j
            line_name = "Q" if j is None else f"Q({j})"
            out.append(
                ResponseComparison(
                    material=gas,
                    reference=key,
                    label=f"{line_name}, {pump_nm:g} nm, {amagat:g} amagat",
                    quantity="gain",
                    unit="cm/GW",
                    expected=float(value),
                    sigma=float(sigma),
                    predicted=_response_gain(response, pump_nm, amagat, line),
                )
            )
    return out


def _revival_decay(
    response: GasRamanResponse, density: float, b: float, periods: float
) -> float:
    """Return the decay time of the rotational revivals, as Chen et al. fit it.

    The peak of |h| at each revival, every ``2 pi / Delta Omega`` (the spacing of the
    two lowest S lines: 4B, or 8B where only one parity of levels exists, as for O2),
    up to ``periods`` times the full period ``1/(2Bc)``, fitted to an exponential
    (Chen et al. 2007, Fig. 5: 'an exponential fit to the declining revival peaks').

    Parameters
    ----------
    response : GasRamanResponse
        The response.
    density : float
        Number density, m^-3.
    b : float
        The rotational constant, cm^-1.
    periods : float
        The fitted span, in full revival periods.

    Returns
    -------
    float
        The decay time, s.
    """
    rotation = np.array(response.band) == "rotation"
    lowest = np.sort(response.frequency[rotation])[:2]
    beat = 2 * math.pi / float(lowest[1] - lowest[0])
    period = 1 / (2 * b * C * 100)
    times, peaks = [], []
    for k in range(1, int(periods * period / beat + 1e-9) + 1):
        t = np.linspace(k * beat - 0.1 * beat, k * beat + 0.1 * beat, 4001)
        h = quiet(response.response, t, density, "anisotropic")
        times.append(k * beat)
        peaks.append(float(np.max(np.abs(h))))
    slope = np.polyfit(times, np.log(peaks), 1)[0]
    return -1 / slope


def rotational_dephasing() -> list[ResponseComparison]:
    """Return the measured decay of the rotational revivals against the responses'.

    Returns
    -------
    list of ResponseComparison
        Chen et al. 2007's O2 at 5.1 atm (ps).
    """
    record = dataset("chen_2007/T2_O2")
    pressure = float(record.abscissa_values()[0])  # Pa (the loader's SI)
    temperature = 295.0  # "room temperature"
    density = pressure / (K_B * temperature)
    response = raman_response("O2", temperature)
    rotational = raman_parameters("O2").rotational
    assert rotational is not None
    b = rotational.constants.b
    return [
        ResponseComparison(
            material="O2",
            reference=record.key,
            label=f"S branch, {pressure / 101325:.3g} atm, revivals to 1.25 periods",
            quantity="revival decay time",
            unit="ps",
            expected=float(record.scaled_values()[0]),
            sigma=math.nan,
            predicted=_revival_decay(response, density, b, 1.25) * 1e12,
        )
    ]
