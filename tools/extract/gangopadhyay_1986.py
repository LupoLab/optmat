# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Gangopadhyay, Tang, Lambropoulos and Shakeshaft (1986): MQDT GCS of Xe.

The paper's generalised cross sections (Eq. 1: W = sigma F^N, F the photon flux) are
printed in two tables, above the 2P1/2 threshold, and drawn in two figures, between
the two thresholds. Read here: the two tables and Fig. 1 (two-photon); Fig. 2, the
three-photon GCS between the thresholds, is not read (the record's notes). The PDF is
a scan: every page is one 1-bit image, 5060 x 6898 pixels (about 610 dpi), read here
from its native pixels (``_figure.image``). PDF page 1 is journal p. 2998.

* **Table I** (p. 3001, library table ``tbl_84c11fc5ced5``): two-photon GCS in cm^4 s
  with the residual core in J = 3/2 and J = 1/2, against the total energy (eV), 14 to
  16.871 eV.
* **Table II** (p. 3002, ``tbl_8ff4b5716399``): three-photon GCS in cm^6 s^2, the same
  two cores, 14 to 14.506 eV.
* **Fig. 1** (p. 3000, library figure ``fig_5c0c5705aedd``): the two-photon GCS inside
  the autoionising region, linear in units of 1e-50 cm^4 s (ticks 0 to 6, the top one
  labelled '6 x 10^-50'), against an unlabelled energy axis, 12.2 to 13.4: the total
  energy of the two photons (the text gives the autoionising region as '12.13 to
  13.45 eV', and Fig. 3's caption: 'The energies E represent the total energy above
  the ground state'). One solid line, 12.20 to 13.20 eV.

Every table cell is printed with its own power of ten (the J = 3/2 and J = 1/2
columns change from 10^-50 to 10^-46, 10^-76 to 10^-72 down the rows): the value is
the mantissa times that power, and ``printed`` keeps the cell. The data files' ``x``
is the **photon** energy in eV, the printed total energy divided by the number of
photons (2 or 3; the record's decision); the row is the printed total energy.

How Fig. 1 is read (no frame; two axis lines):

* **the axis lines**: the left (vertical) and bottom lines are fitted to the image
  (the centre of the dark run nearest the expected position in every other row or
  column, outliers beyond 1.5 pixels rejected); a pixel (column, row) is mapped to
  *corrected* coordinates, the column it would have on the bottom line and the row
  it would have on the left one, through the lines' slopes (the scan is turned by
  0.24 degree, the left line leaning 0.0043 column per row).
* **the ticks**: the x ticks hang below the bottom line, the y ticks stand left of
  the left line; each is found as a column (row) dark over the whole stretch from
  ``TICK_INSET[0]`` to ``TICK_INSET[1]`` pixels outside the line, and must lie where
  ``X_TICKS`` and ``Y_TICKS`` write it, within ``TICK_TOLERANCE``. The 12.2 tick is
  the left line's continuation below the corner (checked: dark there); the 0 tick is
  the bottom line's continuation left of the corner (found as a tick).
* **the scale**: the tick spacing is not uniform (x: 259 pixels per 0.2 eV at the
  left, 253 at the right; y: 257 pixels per unit at the bottom, 274.5 at the top, 6.8 %
  more), a distortion of the drawing or the print that a single linear scale would
  leave as residuals of up to 3.5 (x) and 8.9 (y) pixels at the ticks (reported). The
  scale is therefore piecewise linear between neighbouring ticks (``Scale``), exact
  at every tick, extended beyond the end ticks with the end intervals.
* **the curve**: the axis lines and the ticks are cleared (``CLEAR`` pixels of
  corrected distance from each line), then every connected piece of ink (8-connected)
  lying wholly inside a box of ``LETTERING`` (the 'J = 0', 'J = 2' labels and their
  arrows; each box must hold some) and every speck smaller than ``SPECK`` pixels are
  erased. The curve is sampled on a grid of total energies (``GRID``); at a sample's
  column the ink must form a single run no taller than ``RUN`` pixels (a steep flank
  or a resonance's peak gives a taller run, an arrowhead touching the line a second
  one, and is not read), whose centre is the value. A value is kept only if a
  neighbouring grid point was read and agrees with it (``SLOPE``, ``SLACK``).

Every value's ``printed`` cell gives its pixel (column, row) in the page image and the
total energy at it. One pixel is 0.00078 eV of total energy and 0.0037 to 0.0039 in
units of 1e-50 cm^4 s; the line is 3 to 6 pixels high where it is flat.

**Checks** (:func:`checks`; each raises if missed): every tick; the arrows of Fig. 1
that point at the features whose energies Fig. 3 prints, 'E = 12.602' (a) and 'E =
12.95' (b), must point there within ``ENERGY`` (the curve's dip and bump there are
reported); Table I as other records' figures draw it, against the sum of its two
columns at the photon energy: lagutin_2017's Fig. 8 squares (vector, read from its
drawing with that script's calibration; within ``DRAWN``), lhuillier_1987's Fig. 6
filled circles (raster, found and placed by that script; within ``CIRCLE_ENERGY``
and ``CIRCLE_VALUE``) and lhuillier_1989's Fig. 10(a) crosses (raster, placed with
that script's calibration; within ``CROSS_ENERGY`` and ``CROSS_VALUE``). Reported:
the drawn peak tops, the single-line tick residuals, the value at twice the 193 nm
photon energy against mccown_1982, and lhuillier_1989's solid line where its fifth
cross lies on it.
"""

import csv
import re
from collections.abc import Sequence
from dataclasses import dataclass
from decimal import Decimal

import lagutin_2017
import lhuillier_1987
import lhuillier_1989
import numpy as np
from _figure import image
from _library import DATA, Point, table, write
from scipy import ndimage

DIRECTORY = "Gangopadhyay_1986_TheoryAutoionizationXe"
SOURCE = "gangopadhyay_1986"
PAGE = 3
"""Fig. 1 is on PDF page 3 (journal p. 3000)."""
SHAPE = (6898, 5060)
DARK = 128

TABLES = {
    "tbl_84c11fc5ced5": (2, ("tab1_jc32", "tab1_jc12")),
    "tbl_8ff4b5716399": (3, ("tab2_jc32", "tab2_jc12")),
}
"""Each table: the number of photons, and the datasets of its two GCS columns."""
FIG1 = "fig1"


def _file(key: str) -> str:
    return f"Xe/{SOURCE}__{key}.csv"


# A cell such as '0.96 \times 10^{-50}': the mantissa and the power of ten as printed.
_CELL = re.compile(r"^\s*([0-9.]+)\s*\\times\s*10\^\{(-?\d+)\}\s*$")


def _cell_value(cell: str) -> str:
    """Return the number a cell prints, mantissa times its own power of ten."""
    match = _CELL.match(cell)
    if match is None:
        msg = f"not a number times a power of ten: {cell!r}"
        raise ValueError(msg)
    return f"{match.group(1)}e{match.group(2)}"


def _photon_energy(total: str, photons: int) -> str:
    """Return the photon energy (eV) of a printed total energy, exactly."""
    value = Decimal(total) / photons
    return f"{value:.6f}".rstrip("0").rstrip(".")


def _tables() -> dict[str, list[Point]]:
    """Return the two tables' columns, against the photon energy."""
    out: dict[str, list[Point]] = {}
    for table_id, (photons, keys) in TABLES.items():
        rows = table(DIRECTORY, table_id)
        header, body = rows[0], rows[1:]
        if not header[0].startswith("Total energy (eV)") or len(header) != 3:
            msg = f"{table_id}: an unexpected header {header}"
            raise ValueError(msg)
        for column, key in enumerate(keys, start=1):
            out[_file(key)] = [
                Point(
                    row[0].strip(),
                    _photon_energy(row[0].strip(), photons),
                    _cell_value(row[column]),
                    "",
                    f"{row[0].strip()} | {row[column].strip()}",
                )
                for row in body
            ]
    return out


# ---------------------------------------------------------------------------- Fig. 1

LEFT = 2922
"""The left axis line's approximate column."""
BOTTOM = 6077
"""The bottom axis line's approximate row."""
TOP = 4250
"""The top of the region read (the left line runs up to row 4180)."""
RIGHT = 4620
"""The bottom line's approximate right end."""
TICK_INSET = (8, 20)
"""A tick is dark over this whole stretch outside its line, pixels."""
TICK_TOLERANCE = 2.0
"""A tick lies within this of where it is written, corrected pixels."""
X_TICKS = (
    (3180.0, 12.4),
    (3438.5, 12.6),
    (3695.5, 12.8),
    (3949.5, 13.0),
    (4203.0, 13.2),
    (4456.0, 13.4),
)
"""The x ticks below the bottom line: (corrected column, total energy in eV)."""
Y_TICKS = (
    (6076.0, 0.0),
    (5819.0, 1.0),
    (5556.5, 2.0),
    (5294.0, 3.0),
    (5024.0, 4.0),
    (4753.0, 5.0),
    (4478.0, 6.0),
)
"""The y ticks left of the left line: (corrected row, value in 1e-50 cm^4 s)."""
CORNER_TICK = 12.2
"""The total energy (eV) of the left line, whose continuation below is its tick."""
CLEAR = 6
"""Ink within this many corrected pixels of an axis line is cleared."""
LETTERING = (
    (4285, 3795, 4348, 3918),  # 'J = 2' over the 12.99 eV peak
    (5246, 4077, 5310, 4205),  # 'J = 0' over the 13.12 eV peak
    (5530, 3428, 5603, 3579),  # 'J = 2' over the 12.75 eV peak
    (5650, 3563, 5712, 3628),  # its arrow
    (5572, 3784, 5634, 3903),  # 'J = 2' over 12.95 eV
    (5692, 3850, 5902, 3890),  # its long arrow
    (5722, 3676, 5785, 3797),  # 'J = 2' over 12.89 eV
    (5830, 3768, 5906, 3816),  # its arrow
    (5753, 3378, 5814, 3500),  # 'J = 2' over 12.60 eV
    (5862, 3428, 5950, 3458),  # its arrow
    (5655, 2929, 5716, 3044),  # 'J = 2' over 12.25 eV
    (5732, 2944, 5754, 2957),  # the mark under it
    (5898, 2950, 5958, 2989),  # its arrow
    (5736, 3042, 5798, 3174),  # 'J = 0' over the 12.29 eV peak
)
"""Boxes (top, left, bottom, right; image pixels) of lettering and arrows."""
SPECK = 8
"""Pieces of ink smaller than this (pixels) are dust."""
RUN = 12
"""The tallest run (pixels) a value is read from: the line is 3-6 pixels high."""
SLOPE = 1.5
SLACK = 4.0
"""Neighbouring values agree if their rows differ by at most SLOPE pixels per pixel of
column plus SLACK pixels."""


def _grid(first: float, last: float, step: float) -> tuple[float, ...]:
    """Return the grid from first to last (inclusive) in steps, rounded to 4 places."""
    count = round((last - first) / step)
    return tuple(round(first + step * k, 4) for k in range(count + 1))


GRID = _grid(12.202, 13.198, 0.002)
"""The total energies (eV) where the curve is sampled."""
ENERGY = 0.006
"""Fig. 1's arrows lie at Fig. 3's printed resonance energies within this (eV, total):
Fig. 3's last digit (12.95) and eight pixels."""


@dataclass(frozen=True)
class Line:
    """A fitted axis line: position = slope * along + intercept, in pixels."""

    slope: float
    intercept: float

    def at(self, along: float) -> float:
        """Return the line's position at a row (vertical) or column (horizontal)."""
        return self.slope * along + self.intercept


def _spans(mask: np.ndarray) -> list[tuple[int, int]]:
    """Return the first and last index of each run of True."""
    index = np.flatnonzero(mask)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [(int(group[0]), int(group[-1])) for group in groups]


def _fit_line(dark: np.ndarray, vertical: bool, near: int, span: range) -> Line:
    """Fit an axis line near a column (vertical) or row, over a span of rows or columns.

    The centre of the dark run (at most 9 pixels wide) nearest the expected position
    within 15 pixels, in every other row or column; outliers beyond 1.5 pixels (a tick,
    lettering or the curve meeting the line) are rejected, three times.
    """
    points = []
    for along in span[::2]:
        cut = slice(near - 15, near + 16)
        line = dark[along, cut] if vertical else dark[cut, along]
        runs = [(a + b) / 2 + near - 15 for a, b in _spans(line) if b - a < 9]
        if runs:
            points.append((along, min(runs, key=lambda c: abs(c - near))))
    data = np.array(points)
    for _ in range(3):
        slope, intercept = np.polyfit(data[:, 0], data[:, 1], 1)
        data = data[np.abs(data[:, 1] - (slope * data[:, 0] + intercept)) < 1.5]
    if len(data) < 0.8 * len(span[::2]):
        msg = f"an axis line near {near}: only {len(data)} points on it"
        raise ValueError(msg)
    slope, intercept = np.polyfit(data[:, 0], data[:, 1], 1)
    return Line(float(slope), float(intercept))


@dataclass(frozen=True)
class Axes:
    """The two axis lines, and the corrected coordinates they define."""

    left: Line
    bottom: Line

    @property
    def corner(self) -> tuple[float, float]:
        """The corner (column, row), where the two lines meet."""
        left, bottom = self.left, self.bottom
        column = (left.slope * bottom.intercept + left.intercept) / (
            1 - left.slope * bottom.slope
        )
        return column, bottom.at(column)

    def corrected(self, column: float, row: float) -> tuple[float, float]:
        """Return a pixel's corrected (column, row): sheared onto the axes."""
        c0, r0 = self.corner
        return (
            column - self.left.slope * (row - r0),
            row - self.bottom.slope * (column - c0),
        )

    def column(self, corrected: float, row: float) -> float:
        """Return the image column of a corrected column at an image row."""
        return corrected + self.left.slope * (row - self.corner[1])


@dataclass(frozen=True)
class Scale:
    """A piecewise-linear scale through tick positions (exact at every tick)."""

    positions: tuple[float, ...]
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        """Require the positions to be monotonic with the values."""
        steps = np.diff(self.positions) / np.diff(self.values)
        if not (np.all(steps > 0) or np.all(steps < 0)):
            msg = "the ticks are not in order"
            raise ValueError(msg)

    def _sorted(self) -> tuple[np.ndarray, np.ndarray]:
        order = np.argsort(self.positions)
        return (
            np.asarray(self.positions, dtype=float)[order],
            np.asarray(self.values, dtype=float)[order],
        )

    def value(self, position: float) -> float:
        """Return the value at a position, extended linearly beyond the end ticks."""
        p, v = self._sorted()
        k = int(np.clip(np.searchsorted(p, position) - 1, 0, len(p) - 2))
        return float(v[k] + (v[k + 1] - v[k]) * (position - p[k]) / (p[k + 1] - p[k]))

    def position(self, value: float) -> float:
        """Return the position of a value (the inverse of :meth:`value`)."""
        order = np.argsort(self.values)
        v = np.asarray(self.values, dtype=float)[order]
        p = np.asarray(self.positions, dtype=float)[order]
        k = int(np.clip(np.searchsorted(v, value) - 1, 0, len(v) - 2))
        return float(p[k] + (p[k + 1] - p[k]) * (value - v[k]) / (v[k + 1] - v[k]))

    def linear_residual(self) -> float:
        """Return the worst tick residual (pixels) of one straight-line fit."""
        slope, intercept = np.polyfit(self.values, self.positions, 1)
        fitted = slope * np.asarray(self.values) + intercept
        return float(np.max(np.abs(fitted - np.asarray(self.positions))))


def _match(found: Sequence[float], expected: Sequence[float], what: str) -> None:
    """Raise unless the found positions are the expected ones, each within tolerance."""
    if len(found) != len(expected) or any(
        abs(f - e) > TICK_TOLERANCE
        for f, e in zip(sorted(found), sorted(expected), strict=True)
    ):
        found_ = [round(f, 1) for f in sorted(found)]
        msg = f"{what}: found {found_}, expected {sorted(expected)}"
        raise ValueError(msg)


def _ticks(dark: np.ndarray, axes: Axes) -> tuple[list[float], list[float]]:
    """Return the corrected positions of the ticks found below and left of the axes."""
    inner, outer = TICK_INSET
    columns = range(LEFT + 15, RIGHT)
    found = []
    for column in columns:
        row = round(axes.bottom.at(column))
        found.append(bool(dark[row + inner : row + outer, column].all()))
    bottom = []
    for a, b in _spans(np.array(found)):
        column = (a + b) / 2 + columns.start
        bottom.append(axes.corrected(column, axes.bottom.at(column))[0])
    rows = range(TOP, BOTTOM + 15)
    found = []
    for row in rows:
        column = round(axes.left.at(row))
        found.append(bool(dark[row, column - outer : column - inner].all()))
    left = []
    for a, b in _spans(np.array(found)):
        row = (a + b) / 2 + rows.start
        left.append(axes.corrected(axes.left.at(row), row)[1])
    return bottom, left


@dataclass(frozen=True)
class Calibration:
    """Fig. 1's axis lines and scales: total energy (eV) and 1e-50 cm^4 s."""

    axes: Axes
    x: Scale
    y: Scale


def calibration(dark: np.ndarray) -> Calibration:
    """Fit Fig. 1's axis lines, check every tick and return the scales."""
    axes = Axes(
        _fit_line(dark, True, LEFT, range(TOP, BOTTOM - 20)),
        _fit_line(dark, False, BOTTOM, range(LEFT + 20, RIGHT)),
    )
    bottom, left = _ticks(dark, axes)
    _match(bottom, [p for p, _ in X_TICKS], "Fig. 1: x ticks")
    _match(left, [p for p, _ in Y_TICKS], "Fig. 1: y ticks")
    c0, r0 = axes.corner
    # The 12.2 tick: the left line continues below the corner (broken in the scan).
    below = dark[round(r0) + 4 : round(r0) + 28, round(c0) - 2 : round(c0) + 3]
    if below.any(axis=1).sum() < 15:
        msg = "Fig. 1: no 12.2 tick below the corner"
        raise ValueError(msg)
    x = Scale((c0, *(p for p, _ in X_TICKS)), (CORNER_TICK, *(v for _, v in X_TICKS)))
    y = Scale(tuple(p for p, _ in Y_TICKS), tuple(v for _, v in Y_TICKS))
    return Calibration(axes, x, y)


@dataclass(frozen=True)
class Ink:
    """Fig. 1's curve: the ink of the box ``rows`` x ``columns``, all else erased."""

    rows: range
    columns: range
    mask: np.ndarray


def _ink(dark: np.ndarray, axes: Axes) -> Ink:
    """Clear the axes, erase the lettering, arrows and dust; return the curve's ink."""
    rows, columns = range(TOP, BOTTOM + 1), range(LEFT, RIGHT + 1)
    box = dark[rows.start : rows.stop, columns.start : columns.stop].copy()
    r, c = np.meshgrid(
        np.arange(rows.start, rows.stop),
        np.arange(columns.start, columns.stop),
        indexing="ij",
    )
    c0, r0 = axes.corner
    corrected_c = c - axes.left.slope * (r - r0)
    corrected_r = r - axes.bottom.slope * (c - c0)
    box &= (corrected_c > c0 + CLEAR) & (corrected_r < r0 - CLEAR)
    labels, _ = ndimage.label(box, structure=np.ones((3, 3), dtype=int))
    used = set()
    for index, where in enumerate(ndimage.find_objects(labels), start=1):
        if where is None:
            continue
        top, bottom = where[0].start + rows.start, where[0].stop + rows.start
        left, right = where[1].start + columns.start, where[1].stop + columns.start
        piece = labels[where] == index
        inside = [
            k
            for k, (t, w, b, e) in enumerate(LETTERING)
            if t <= top and bottom <= b and w <= left and right <= e
        ]
        if inside or piece.sum() < SPECK:
            box[where] &= ~piece
            used.update(inside)
    if len(used) != len(LETTERING):
        missing = [LETTERING[k] for k in range(len(LETTERING)) if k not in used]
        msg = f"Fig. 1: no lettering wholly inside {missing}"
        raise ValueError(msg)
    return Ink(rows, columns, box)


def _run(ink: Ink, column: int) -> float | None:
    """Return the centre row of the curve's single short run in an image column."""
    c = column - ink.columns.start
    if not 0 <= c < ink.mask.shape[1]:
        return None
    spans = _spans(ink.mask[:, c])
    if len(spans) != 1 or spans[0][1] - spans[0][0] + 1 > RUN:
        return None
    return (spans[0][0] + spans[0][1]) / 2 + ink.rows.start


@dataclass(frozen=True)
class Sample:
    """One value read from Fig. 1: where in the image and what it is."""

    grid: float
    column: int
    row: float
    total: float
    """The total energy at the pixel, eV."""
    sigma: float
    """In 1e-50 cm^4 s."""


def _sample(ink: Ink, cal: Calibration, energy: float, guess: float) -> Sample | None:
    """Read the curve at a total energy, or None where it cannot be read.

    The column is that of the energy at the curve's row (found by iteration from a
    guess of the row).
    """
    corrected = cal.x.position(energy)
    column = round(cal.axes.column(corrected, guess))
    for _ in range(3):
        row = _run(ink, column)
        if row is None:
            return None
        new = round(cal.axes.column(corrected, row))
        if new == column:
            break
        column = new
    else:
        return None
    c, r = cal.axes.corrected(column, row)
    return Sample(energy, column, row, cal.x.value(c), cal.y.value(r))


def _consistent(samples: list[Sample]) -> list[Sample]:
    """Return the samples that agree with a sample at a neighbouring grid point."""
    index = {energy: k for k, energy in enumerate(GRID)}
    by_index = {index[s.grid]: s for s in samples}

    def agrees(a: Sample, b: Sample) -> bool:
        return abs(a.row - b.row) <= SLOPE * abs(a.column - b.column) + SLACK

    return [
        sample
        for k, sample in sorted(by_index.items())
        if any(agrees(sample, by_index[j]) for j in (k - 1, k + 1) if j in by_index)
    ]


@dataclass(frozen=True)
class Reading:
    """Fig. 1 as read: the calibration, the curve's ink and the samples."""

    calibration: Calibration
    ink: Ink
    samples: list[Sample]


def reading() -> Reading:
    """Return Fig. 1's reading, the ticks checked."""
    grey = image(DIRECTORY, PAGE)
    if grey.shape != SHAPE:
        msg = f"page {PAGE}: an image of {grey.shape}, not {SHAPE}"
        raise ValueError(msg)
    dark = grey < DARK
    cal = calibration(dark)
    ink = _ink(dark, cal.axes)
    samples = []
    guess = float(BOTTOM - 100)
    for energy in GRID:
        sample = _sample(ink, cal, energy, guess)
        if sample is not None:
            samples.append(sample)
            guess = sample.row
    return Reading(cal, ink, _consistent(samples))


def _extreme(read: Reading, low: float, high: float, lowest: bool) -> float:
    """Return the total energy (eV) of the curve's lowest or highest ink in a window.

    The ink's bottom (``lowest``) or top edge over the columns of the window; the
    energy is that of the middle of the columns sharing the extreme row.
    """
    cal, ink = read.calibration, read.ink
    first = round(cal.x.position(low)) - ink.columns.start
    last = round(cal.x.position(high)) - ink.columns.start
    edges = []
    for c in range(first, last + 1):
        rows = np.flatnonzero(ink.mask[:, c])
        if rows.size:
            edges.append((c, int(rows.max() if lowest else rows.min())))
    extreme = max(r for _, r in edges) if lowest else min(r for _, r in edges)
    at = [c for c, r in edges if r == extreme]
    column = (at[0] + at[-1]) / 2 + ink.columns.start
    row = extreme + ink.rows.start
    return cal.x.value(cal.axes.corrected(column, row)[0])


def _arrow_tip(dark: np.ndarray, cal: Calibration, box: tuple[int, ...]) -> float:
    """Return the total energy (eV) an arrow points at: its lowest ink's column.

    The arrow is the only ink in its box (one of ``LETTERING``); its tip is the
    middle of the dark pixels of its lowest row.
    """
    top, left, bottom, right = box
    rows, cols = np.nonzero(dark[top:bottom, left:right])
    tip = rows.max()
    column = float(cols[rows == tip].mean()) + left
    return cal.x.value(cal.axes.corrected(column, float(tip + top))[0])


def _peak_top(read: Reading, low: float, high: float) -> tuple[float, float]:
    """Return the total energy (eV) and value (1e-50 cm^4 s) of a drawn peak's top."""
    cal, ink = read.calibration, read.ink
    first = round(cal.x.position(low)) - ink.columns.start
    last = round(cal.x.position(high)) - ink.columns.start
    window = ink.mask[:, first : last + 1]
    rows, cols = np.nonzero(window)
    k = int(np.argmin(rows))
    c, r = cal.axes.corrected(
        cols[k] + first + ink.columns.start, rows[k] + ink.rows.start
    )
    return cal.x.value(c), cal.y.value(r)


# The resonance energies Fig. 3 prints (total energy, eV), the arrows of Fig. 1 that
# mark them (their boxes in LETTERING) and the windows of the curve's feature there:
# (a) 'E = 12.602' (labelled J = 0 there; Fig. 1's label says J = 2), where the curve
# drops into a dip; (b) 'E = 12.95' (J = 2), the top of a bump.
FIG3 = {
    "12.602": ((5862, 3428, 5950, 3458), (12.58, 12.63, True)),
    "12.95": ((5692, 3850, 5902, 3890), (12.93, 12.97, False)),
}
# The drawn peaks (total energy windows, eV), reported.
PEAKS = {
    "J = 0, 12.29": (12.27, 12.31),
    "J = 2, 12.75": (12.73, 12.77),
    "J = 2, 12.99": (12.97, 13.01),
    "J = 0, 13.12": (13.10, 13.14),
}
HC_EV_NM = 1239.841984
"""h c in eV nm (CODATA 2018, exact)."""
MCCOWN_SIGMA = 4.11e-50
"""mccown_1982's alpha = 4e-32 cm^4 W^-1 at 193 nm as sigma = alpha hbar omega (that
record's notes), cm^4 s."""


LAGUTIN = "fig8a"
"""lagutin_2017's plot that draws Table I as squares (Fig. 8, top)."""
SQUARE = (3.0, 5.0)
"""Its squares are closed four-sided paths this wide and high, pt."""
DRAWN = 0.01
"""lagutin_2017's squares agree with the sums of Table I within this (relative): its
drawing's quantum (0.13 pt on 30 pt per decade) and Table I's two figures."""


def lagutin_squares() -> list[tuple[float, float]]:
    """Return lagutin_2017's Fig. 8 squares: photon energy (eV), value (cm^4 s).

    Read from its vector drawing with that script's calibration: the closed,
    unfilled four-sided paths of the square's size inside the plot, outside its
    legend's box.
    """
    panel = lagutin_2017.PANELS[LAGUTIN]
    box = lagutin_2017.frame(LAGUTIN)
    x_axis, y_axis = lagutin_2017.axes(LAGUTIN)
    out = []
    for stroke in lagutin_2017._strokes(panel.page):
        v = stroke.vertices
        width, height = np.ptp(v[:, 0]), np.ptp(v[:, 1])
        if (
            stroke.filled
            or len(v) != 5
            or not SQUARE[0] < width < SQUARE[1]
            or not SQUARE[0] < height < SQUARE[1]
            or not np.all(lagutin_2017._inside(v, box))
            or any(np.all(lagutin_2017._inside(v, b)) for b in panel.boxes)
        ):
            continue
        centre = v[:4].mean(axis=0)
        out.append(
            (
                float(x_axis(centre[0])),
                float(y_axis(centre[1])) * lagutin_2017.UNIT,
            )
        )
    return sorted(out)


def _check_lagutin(out: dict[str, object], failures: list[str]) -> None:
    """Compare lagutin_2017's squares with the sums of Table I's two columns."""
    tables = _tables()
    jc32, jc12 = tables[_file("tab1_jc32")], tables[_file("tab1_jc12")]
    sums = {
        float(a.x): float(a.value) + float(b.value)
        for a, b in zip(jc32, jc12, strict=True)
    }
    squares = lagutin_squares()
    if len(squares) != 6:
        failures.append(f"lagutin_2017's Fig. 8: {len(squares)} squares, not 6")
        return
    for energy, value in squares:
        photon = min(sums, key=lambda e: abs(e - energy))
        ratio = value / sums[photon]
        out[f"lagutin_2017 square at {energy:.4f} eV over Table I sum"] = round(
            ratio, 4
        )
        if abs(energy - photon) > 0.005 or abs(ratio - 1) > DRAWN:
            failures.append(f"lagutin_2017's square at {energy:.4f} eV: {ratio:.4f}")


RYDBERG_EV = 13.605693122994
"""R_inf h c in eV (CODATA 2018)."""
CM1_PER_EV = 8065.543937
"""1 eV in cm^-1 (CODATA 2018)."""
CIRCLE_FILL = 0.96
"""lhuillier_1987's Fig. 6: a filled marker whose box is less dark than this is a
circle (Gangopadhyay et al.'s; 0.92-0.955), a fuller one a square (McGuire's,
0.97-0.995)."""
CIRCLE_ENERGY = 0.003
CIRCLE_VALUE = 0.06
"""lhuillier_1987's circles lie within this of Table I's photon energies (Ry) and sums
(relative): hand-placed markers, 19-24 pixels across, on a scan (one pixel 1.3 %)."""
CROSS_SIZE = (20, 70)
"""lhuillier_1989's Fig. 10(a) crosses: connected pieces this many pixels high and
wide."""
CROSS_ENERGY = 0.5
CROSS_VALUE = 0.08
"""lhuillier_1989's crosses lie within this of Table I's total energies (1e3 cm^-1)
and sums (relative): hand-placed, 36-38 pixels across, on a log axis whose ticks fit
within 6.5 pixels (that script's tolerance)."""


def _sums() -> dict[str, tuple[float, float]]:
    """Return Table I's rows: total energy (as printed) -> (eV, sum of its columns)."""
    tables = _tables()
    jc32, jc12 = tables[_file("tab1_jc32")], tables[_file("tab1_jc12")]
    return {
        a.row: (float(a.row), float(a.value) + float(b.value))
        for a, b in zip(jc32, jc12, strict=True)
    }


def _check_lhuillier_1987(out: dict[str, object], failures: list[str]) -> None:
    """Compare lhuillier_1987's Fig. 6 circles with the sums of Table I.

    The filled markers are found by that script (an opening of its marker box) and
    placed with its calibration (photon energy in Ry, 1e-49 cm^4 s); the circles are
    those less dark than :data:`CIRCLE_FILL` over their box.
    """
    panel = next(p for p in lhuillier_1987.PANELS if p.markers)
    dark = image(lhuillier_1987.DIRECTORY, lhuillier_1987.PAGE) < lhuillier_1987.DARK
    frame, _, x_axis, y_axis = lhuillier_1987._axes(dark, panel)
    circles = []
    for top, left, bottom, right in lhuillier_1987._markers(dark, panel):
        if dark[top:bottom, left:right].mean() < CIRCLE_FILL:
            c, r = frame.corrected((left + right - 1) / 2, (top + bottom - 1) / 2)
            circles.append((float(x_axis(c)), float(y_axis(r)) * 1e-49))
    sums = _sums()
    if len(circles) != len(sums) - 1:
        failures.append(f"lhuillier_1987's Fig. 6: {len(circles)} circles, not 6")
        return
    for ry, value in sorted(circles):
        row = min(sums, key=lambda k: abs(sums[k][0] / 2 / RYDBERG_EV - ry))
        energy, total = sums[row]
        ratio = value / total
        out[f"lhuillier_1987 circle at {ry:.4f} Ry over Table I {row} eV"] = round(
            ratio, 3
        )
        if (
            abs(energy / 2 / RYDBERG_EV - ry) > CIRCLE_ENERGY
            or abs(ratio - 1) > CIRCLE_VALUE
        ):
            failures.append(f"lhuillier_1987's circle at {ry:.4f} Ry: {ratio:.3f}")


def _check_lhuillier_1989(out: dict[str, object], failures: list[str]) -> None:
    """Compare lhuillier_1989's Fig. 10(a) crosses with the sums of Table I.

    The crosses are the connected pieces of ink of a cross's size inside that
    script's Xe frame that lie near one of Table I's total energies (the letters 'Xe'
    and '(a)' lie elsewhere), placed with that script's calibration (total energy in
    1e3 cm^-1, cm^4 s). The cross at 129e3 cm^-1 lies on the solid line (that record's
    notes) and is not separated from it.
    """
    panel = next(p for p in lhuillier_1989.PANELS if p.figure == "10" and p.gas == "Xe")
    dark = image(lhuillier_1989.DIRECTORY, panel.page) < lhuillier_1989.DARK
    frame, _, x_axis, y_axis = lhuillier_1989._axes(dark, panel)
    left, right, top, bottom = panel.frame
    inner = dark[top + 8 : bottom - 8, left + 8 : right - 8]
    labels, _ = ndimage.label(inner, structure=np.ones((3, 3), dtype=int))
    sums = {row: (e * CM1_PER_EV / 1e3, s) for row, (e, s) in _sums().items()}
    found = []
    for index, where in enumerate(ndimage.find_objects(labels), start=1):
        if where is None:
            continue
        high = where[0].stop - where[0].start
        wide = where[1].stop - where[1].start
        if not (CROSS_SIZE[0] < high < CROSS_SIZE[1]) or not (
            CROSS_SIZE[0] < wide < CROSS_SIZE[1]
        ):
            continue
        rows, cols = np.nonzero(labels[where] == index)
        c, r = frame.corrected(
            float(cols.mean()) + where[1].start + left + 8,
            float(rows.mean()) + where[0].start + top + 8,
        )
        energy = float(x_axis(c))
        row = min(sums, key=lambda k: abs(sums[k][0] - energy))
        if abs(sums[row][0] - energy) < 2 * CROSS_ENERGY:
            found.append((row, energy, float(y_axis(r))))
    if len(found) != 4:
        failures.append(f"lhuillier_1989's Fig. 10(a): {len(found)} crosses, not 4")
        return
    for row, energy, value in sorted(found):
        ratio = value / sums[row][1]
        out[f"lhuillier_1989 cross at {energy:.2f}e3 cm-1 over Table I {row} eV"] = (
            round(ratio, 3)
        )
        if abs(sums[row][0] - energy) > CROSS_ENERGY or abs(ratio - 1) > CROSS_VALUE:
            failures.append(f"lhuillier_1989's cross at {energy:.2f}: {ratio:.3f}")
    # The fifth cross (16.001 eV) lies on lhuillier_1989's solid line: report the
    # line there, interpolated in its committed file.
    path = DATA / "values" / "ionisation" / "Xe" / "lhuillier_1989__fig10_total_Xe.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        line = [(float(r["row"]), float(r["value"])) for r in csv.DictReader(handle)]
    at, total = sums["16.001"]
    x, y = zip(*line, strict=True)
    out["lhuillier_1989 line at Table I 16.001 eV over its sum"] = round(
        float(np.interp(at, x, y)) / total, 3
    )


def checks() -> dict[str, object]:
    """Return the checks of Fig. 1 against Fig. 3's printed energies, and reports.

    Raises if a printed energy is missed by more than :data:`ENERGY`.
    """
    read = reading()
    dark = image(DIRECTORY, PAGE) < DARK
    out: dict[str, object] = {}
    failures = []
    for printed, (box, (low, high, lowest)) in FIG3.items():
        tip = _arrow_tip(dark, read.calibration, box)
        feature = _extreme(read, low, high, lowest)
        out[f"Fig. 3 'E = {printed}': arrow tip, curve's feature (eV)"] = (
            round(tip, 4),
            round(feature, 4),
        )
        if abs(tip - float(printed)) > ENERGY:
            failures.append(f"Fig. 1's arrow at {tip:.4f} eV, Fig. 3 {printed}")
    for name, window in PEAKS.items():
        energy, value = _peak_top(read, *window)
        out[f"drawn peak top {name} (eV, 1e-50 cm4 s)"] = (
            round(energy, 4),
            round(value, 3),
        )
    out["linear-fit tick residuals x, y (pixels)"] = (
        round(read.calibration.x.linear_residual(), 1),
        round(read.calibration.y.linear_residual(), 1),
    )
    total_193 = 2 * HC_EV_NM / 193.0
    near = [s for s in read.samples if abs(s.total - total_193) < 0.004]
    out["at 2 x 193 nm (eV total; 1e-50 cm4 s; over mccown_1982)"] = (
        round(total_193, 4),
        [round(s.sigma, 3) for s in near],
        [round(s.sigma * 1e-50 / MCCOWN_SIGMA, 3) for s in near],
    )
    _check_lagutin(out, failures)
    _check_lhuillier_1987(out, failures)
    _check_lhuillier_1989(out, failures)
    if failures:
        raise ValueError("; ".join(failures))
    return out


def _fig1() -> list[Point]:
    """Return Fig. 1's values against the photon energy (half the total), eV."""
    return [
        Point(
            f"{s.grid:.3f}",
            f"{s.total / 2:.5f}",
            f"{s.sigma:.4f}",
            "",
            f"pixel ({s.column}, {s.row:.1f}); total energy {s.total:.4f} eV",
        )
        for s in reading().samples
    ]


def extract() -> dict[str, list[Point]]:
    """Return the tables' columns and Fig. 1's curve, after the checks."""
    checks()
    return {**_tables(), _file(FIG1): _fig1()}


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
