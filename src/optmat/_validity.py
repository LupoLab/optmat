# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The qualified intervals of a model: where it may be evaluated without a warning.

A model's source supports it over one or more disjoint intervals of vacuum wavelength
(OEP 1, section 4.1; OEP 3). ``wavelength_range`` is their envelope. Every model has
one interval, its range, unless its source qualifies several.
"""

import itertools

import numpy as np
from jaxtyping import Bool, Float

type Intervals = tuple[tuple[float, float], ...]


def qualified_intervals(
    intervals: Intervals, envelope: tuple[float, float], name: str
) -> Intervals:
    """Return a model's intervals, checked; its envelope where none are given.

    Parameters
    ----------
    intervals : tuple of tuple of float
        ``(low, high)`` vacuum wavelengths, m; empty for the envelope alone.
    envelope : tuple of float
        The model's ``wavelength_range``, m.
    name : str
        The model's name, for error messages.

    Returns
    -------
    tuple of tuple of float
        The intervals as floats.

    Raises
    ------
    ValueError
        If an interval is empty or reversed, the intervals overlap, touch or are out
        of order, or their envelope is not ``envelope``.
    """
    checked = tuple((float(low), float(high)) for low, high in intervals or (envelope,))
    if any(not low < high for low, high in checked):
        msg = f"{name}: every interval must have low < high, got {checked}"
        raise ValueError(msg)
    if any(a[1] >= b[0] for a, b in itertools.pairwise(checked)):
        msg = f"{name}: the intervals must be disjoint and ascending, got {checked}"
        raise ValueError(msg)
    if (checked[0][0], checked[-1][1]) != (float(envelope[0]), float(envelope[1])):
        msg = f"{name}: the intervals {checked} do not span its range {envelope}"
        raise ValueError(msg)
    return checked


def outside(
    wavelength: Float[np.ndarray, "*shape"], intervals: Intervals
) -> Bool[np.ndarray, "*shape"]:
    """Return where wavelengths lie outside every interval.

    Parameters
    ----------
    wavelength : numpy.ndarray
        Vacuum wavelength, m.
    intervals : tuple of tuple of float
        The qualified intervals, m.

    Returns
    -------
    numpy.ndarray
        True outside every interval. A NaN wavelength is not outside, so it does not
        warn.
    """
    result = np.ones(np.shape(wavelength), dtype=bool)
    for low, high in intervals:
        result &= (wavelength < low) | (wavelength > high)
    return result


def describe(intervals: Intervals) -> str:
    """Return the intervals as text, in nm.

    Parameters
    ----------
    intervals : tuple of tuple of float
        The qualified intervals, m.

    Returns
    -------
    str
        Such as ``"167-2000 nm"``, or ``"250-1940, 2210-2590 nm"`` for several.
    """
    parts = ", ".join(f"{low * 1e9:.4g}-{high * 1e9:.4g}" for low, high in intervals)
    return f"{parts} nm"
