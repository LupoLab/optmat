# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract L'Huillier and Wendin (1987): RPAE two-photon cross sections of Xe.

The two-photon ionisation cross sections of the Xe 5p shell are in figures only (data
rule 8, TC-699). The PDF is a scan: every page is one 1-bit image, 5169 x 7129 pixels
(600 dpi), read here from its native pixels (``_figure.image``).

* **Fig. 6** (p. 4753, library figure ``fig_1027dbe08272``): sigma_2 in 1e-49 cm^4 s
  on a logarithmic axis (1e-2 to 1e5) against the photon energy, 0.46 to 0.86 Ry, the
  discrete-resonance region from the two-photon HF threshold: Hartree-Fock (solid
  line), RPAE with nonlinear screening (dashed) and LDRPA with nonlinear screening
  (dotted). McGuire's squares, Gangopadhyay et al.'s circles and McCown et al.'s cross
  are not taken.
* **Fig. 7** (p. 4753, library figure ``fig_39a6b612098f``): the same three
  calculations above the one-photon threshold (above-threshold ionisation), 0.98 to
  1.40 Ry, 1e-2 to 1e2 in 1e-49 cm^4 s.

How each panel is read:

* **the frame**: its four lines are fitted to the image (the centre of the dark run
  near the expected position in every other row or column, outliers beyond 1.5 pixels
  rejected); a pixel (column, row) is mapped to *corrected* coordinates, the column it
  would have on the bottom frame line and the row it would have on the left one,
  through the mean slopes of the frame's lines (the scan is turned by up to 0.2
  degree).
* **the ticks**: the inward tick marks are found along the bottom and left frame
  lines (a row or column dark from 5 to 10 pixels inside the line) and must be where
  ``x_ticks`` and ``y_ticks`` write them, within ``TICK_TOLERANCE`` pixels; dark
  stretches there that are not ticks (a curve meeting the frame) are written in
  ``not_x_ticks`` and ``not_y_ticks`` and must be found too. Labelled frame lines are
  end values (Fig. 6: the bottom, 1e-2; Fig. 7: the bottom, 1e-2, the top, 1e2, and
  the right, 1.4 Ry). The axes are calibrated from those positions
  (``_figure.Axis``), every tick checked on the fitted scale.
* **the curves**: the frame lines and tick marks are erased, then the lettering (boxes
  written in ``labels``; and McCown et al.'s cross, ``cross``) and the filled markers
  of Fig. 6 (pieces that survive an opening by a ``MARKER`` x ``MARKER`` square, found
  in ``marker_box``; their count is checked) with a margin of ``MARKER_MARGIN``
  pixels. The remaining ink is split
  into connected pieces (8-connected): pieces at least ``SOLID`` pixels long or high are
  the **solid** line, pieces no larger than ``DOT`` pixels the **dotted** line, the rest
  the **dashed** line. Each curve is sampled on a grid of photon energies (``grid``);
  at a sample's column the curve's ink must form a single run no taller than ``RUN``
  pixels (a steep flank, a resonance's peak, or two curves' pieces in one column give a
  taller run or two, and are not read), whose centre is the value. A gap between dashes
  or dots is bridged by the linear interpolation between the nearest readable columns
  within ``GAP`` pixels, if no ink of that curve lies between. A value is kept only if
  the value at a neighbouring grid point was read and agrees with it (``SLOPE``,
  ``SLACK``); each curve is read only within its ``ranges``.

Every value's ``printed`` cell gives its pixel (column, row) in the page image. One
pixel is 0.00024 Ry and 1.33 % in Fig. 6, 0.00026 Ry and 1.30 % in Fig. 7; the lines
are 4-7 pixels wide.

Check (``_check_mccown``): McCown et al.'s cross in Fig. 6 lies at the photon energy of
193 nm (0.472 Ry, the text's '0.47 Ry'); its height is compared with McCown et al.'s
printed coefficient, 4e-32 cm^4 W^-1, as sigma = alpha hbar omega.
"""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from _figure import Axis, image
from _library import Point, write
from scipy import ndimage

DIRECTORY = "LHuillier_1987_LinearNonlinearScreening"
SHAPE = (7129, 5169)
DARK = 128
PAGE = 7
# How far, pixels, inside a frame line a tick is looked for (a row or column dark over
# this whole stretch: the minor ticks are 10-12 pixels long), how far from where it is
# written it may be found, and how much of the frame and of each tick (across it) is
# erased before the curves are read.
TICK_INSET = (5, 10)
TICK_TOLERANCE = 2.0
FRAME_CLEAR = 5
TICK_CLEAR = 5
# The ticks' largest residual from the fitted scales, pixels.
AXIS_TOLERANCE = 2.0
# A piece of ink this long or high (pixels) is part of the solid line; a piece no
# larger than DOT is a dot; dashes lie between.
SOLID = 150
DOT = 10
# The filled markers (squares and circles, 20-26 pixels) survive an opening by a
# MARKER-pixel square, the lines (at most 7 pixels wide away from the resonances) do
# not; a marker is erased with MARKER_MARGIN pixels around it.
MARKER = 13
MARKER_MARGIN = 4
# The tallest run (pixels) a value is read from: the lines are 4-7 pixels wide; a
# taller run is a steep stretch (more than about 1 pixel per pixel), a peak, or two
# pieces together.
RUN = 12
# How far (pixels) either side a sample in a gap of the dashes or dots looks for the
# line.
GAP = 25
CURVES = ("hf", "rpae", "ldrpa")
"""The solid, dashed and dotted lines."""


@dataclass(frozen=True)
class Panel:
    """One panel: where it is, its ticks and how its curves are read.

    Tick positions are corrected pixel coordinates (the module's docstring): a column
    for the x axis, a row for the y axis. ``where`` is ``left``, ``right``, ``bottom``
    or ``top`` for a frame line that carries a labelled end value, or ``tick``.
    """

    figure: str
    frame: tuple[int, int, int, int]
    """The frame lines' approximate columns (left, right) and rows (top, bottom)."""
    x_ticks: tuple[tuple[float, float, str], ...]
    """(position, photon energy in Ry, where)."""
    y_ticks: tuple[tuple[float, float, str], ...]
    """(position, sigma_2 in 1e-49 cm^4 s, where)."""
    ends: tuple[float, float]
    """The photon energies (Ry) the text gives for the frame's left and right lines."""
    grid: tuple[float, ...]
    """The photon energies (Ry) where the curves are sampled."""
    files: dict[str, str]
    """The data file of each curve read (``hf``, ``rpae``, ``ldrpa``)."""
    not_x_ticks: tuple[float, ...] = ()
    """Corrected columns of dark stretches at the bottom line that are not ticks."""
    not_y_ticks: tuple[float, ...] = ()
    """Corrected rows of dark stretches at the left line that are not ticks."""
    labels: tuple[tuple[int, int, int, int], ...] = ()
    """Boxes of lettering inside the frame (top, left, bottom, right; image pixels)."""
    marker_box: tuple[int, int, int, int] | None = None
    """Where the filled markers are (top, left, bottom, right; image pixels)."""
    markers: int = 0
    """How many filled markers lie there."""
    ranges: dict[str, tuple[float, float]] | None = None
    """Where each curve is read (Ry), if not over the whole grid."""
    cross: tuple[int, int, int, int] | None = None
    """McCown et al.'s cross (top, left, bottom, right; image pixels), erased."""


def _grid(first: float, last: float, step: float) -> tuple[float, ...]:
    """Return the grid from first to last (inclusive) in steps, rounded to 4 places."""
    count = round((last - first) / step)
    return tuple(round(first + step * k, 4) for k in range(count + 1))


PANELS = (
    Panel(
        figure="6",
        frame=(764, 2405, 4496, 5703),
        # The bottom ticks every 0.05 Ry, labelled 0.5 to 0.8; the frame lines are
        # not labelled (the text: 0.46 to 0.86 Ry, checked).
        x_ticks=(
            (928.0, 0.50, "tick"),
            (1133.5, 0.55, "tick"),
            (1338.5, 0.60, "tick"),
            (1544.5, 0.65, "tick"),
            (1750.0, 0.70, "tick"),
            (1954.0, 0.75, "tick"),
            (2158.5, 0.80, "tick"),
            (2364.5, 0.85, "tick"),
        ),
        # A tick every decade, labelled 1e-2 (the bottom line), 1e0, 1e2 and 1e4.
        y_ticks=(
            (5702.7, 1e-2, "bottom"),
            (5530.5, 1e-1, "tick"),
            (5358.5, 1e0, "tick"),
            (5185.5, 1e1, "tick"),
            (5013.5, 1e2, "tick"),
            (4841.5, 1e3, "tick"),
            (4668.5, 1e4, "tick"),
        ),
        ends=(0.46, 0.86),
        grid=_grid(0.462, 0.858, 0.002),
        files={
            "hf": "Xe/lhuillier_1987__fig6_hf_Xe.csv",
            "rpae": "Xe/lhuillier_1987__fig6_rpae_Xe.csv",
            "ldrpa": "Xe/lhuillier_1987__fig6_ldrpa_Xe.csv",
        },
        # The solid line (1.4) and the dotted and dashed lines (0.26) meet the left
        # line.
        not_y_ticks=(5332.5, 5460.0),
        labels=(
            (4625, 872, 4679, 953),  # Xe
            (4912, 1535, 4966, 1599),  # 6s
            (4599, 1877, 4654, 2090),  # <-5d->
            (4518, 2217, 4573, 2280),  # 7s
            (4688, 2272, 4743, 2345),  # 6d
        ),
        marker_box=(5300, 780, 5560, 1420),
        markers=13,
        # Above 0.79 Ry (HF) and 0.765 Ry (RPAE) the resonances of the solid and
        # dashed lines (5d, 7s, 6d, 8s) overlap and cross: their pieces cannot be
        # told apart.
        ranges={"hf": (0.462, 0.790), "rpae": (0.462, 0.765), "ldrpa": (0.462, 0.858)},
        cross=(5410, 798, 5441, 829),
    ),
    Panel(
        figure="7",
        frame=(3011, 4603, 800, 1511),
        # The bottom ticks every 0.05 Ry, labelled 1.0 to 1.4 (the right line); the
        # left line is not labelled (the text: 0.98 Ry, checked).
        x_ticks=(
            (3086.5, 1.00, "tick"),
            (3275.0, 1.05, "tick"),
            (3464.0, 1.10, "tick"),
            (3653.0, 1.15, "tick"),
            (3843.5, 1.20, "tick"),
            (4032.5, 1.25, "tick"),
            (4221.5, 1.30, "tick"),
            (4411.5, 1.35, "tick"),
            (4602.5, 1.40, "right"),
        ),
        y_ticks=(
            (1510.9, 1e-2, "bottom"),
            (1333.5, 1e-1, "tick"),
            (1156.0, 1e0, "tick"),
            (980.0, 1e1, "tick"),
            (799.7, 1e2, "top"),
        ),
        ends=(0.98, 1.40),
        grid=_grid(0.985, 1.395, 0.005),
        files={
            "hf": "Xe/lhuillier_1987__fig7_hf_Xe.csv",
            "rpae": "Xe/lhuillier_1987__fig7_rpae_Xe.csv",
            "ldrpa": "Xe/lhuillier_1987__fig7_ldrpa_Xe.csv",
        },
        # The dashed line (1.1) and the solid line (0.03) meet the left line.
        not_y_ticks=(1147.5, 1424.5),
        labels=((935, 4306, 991, 4541),),  # Xe ATI
    ),
)


@dataclass(frozen=True)
class Line:
    """A fitted frame line: position = slope * along + intercept, in pixels."""

    slope: float
    intercept: float

    def at(self, along: float) -> float:
        """Return the line's position at a row (vertical) or column (horizontal)."""
        return self.slope * along + self.intercept

    def over(self, along: np.ndarray) -> np.ndarray:
        """Return the line's positions at an array of rows or columns."""
        return self.slope * along + self.intercept


def _spans(mask: np.ndarray) -> list[tuple[int, int]]:
    """Return the first and last index of each run of True."""
    index = np.flatnonzero(mask)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [(int(group[0]), int(group[-1])) for group in groups]


def _fit_line(dark: np.ndarray, vertical: bool, near: int, span: range) -> Line:
    """Fit a frame line near a column (vertical) or row, over a span of rows or columns.

    The centre of the dark run (at most 8 pixels wide) nearest the expected position
    within 12 pixels, in every other row or column; outliers beyond 1.5 pixels (a tick
    or a curve meeting the line) are rejected, three times.
    """
    points = []
    for along in span[::2]:
        cut = slice(near - 12, near + 13)
        line = dark[along, cut] if vertical else dark[cut, along]
        runs = [(a + b) / 2 + near - 12 for a, b in _spans(line) if b - a < 8]
        if runs:
            points.append((along, min(runs, key=lambda c: abs(c - near))))
    data = np.array(points)
    for _ in range(3):
        slope, intercept = np.polyfit(data[:, 0], data[:, 1], 1)
        keep = np.abs(data[:, 1] - (slope * data[:, 0] + intercept)) < 1.5
        data = data[keep]
    if len(data) < 0.6 * len(span[::2]):
        msg = f"a frame line near {near}: only {len(data)} points on it"
        raise ValueError(msg)
    slope, intercept = np.polyfit(data[:, 0], data[:, 1], 1)
    return Line(float(slope), float(intercept))


@dataclass(frozen=True)
class Frame:
    """A panel's frame lines and the shear correction they give."""

    left: Line
    right: Line
    top: Line
    bottom: Line

    @property
    def shear(self) -> float:
        """Columns per row of the vertical lines (the mean of the two)."""
        return (self.left.slope + self.right.slope) / 2

    @property
    def tilt(self) -> float:
        """Rows per column of the horizontal lines (the mean of the two)."""
        return (self.top.slope + self.bottom.slope) / 2

    @property
    def corner(self) -> tuple[float, float]:
        """The bottom-left corner (column, row), where left and bottom lines meet."""
        left, bottom = self.left, self.bottom
        column = (left.slope * bottom.intercept + left.intercept) / (
            1 - left.slope * bottom.slope
        )
        return column, bottom.at(column)

    def corrected(self, column: float, row: float) -> tuple[float, float]:
        """Return a pixel's corrected (column, row): sheared onto the axes."""
        c0, r0 = self.corner
        return column - self.shear * (row - r0), row - self.tilt * (column - c0)

    def column(self, corrected: float, row: float) -> float:
        """Return the image column of a corrected column at an image row."""
        return corrected + self.shear * (row - self.corner[1])

    def ends(self) -> dict[str, float]:
        """Return the corrected positions of the four lines."""
        c0, r0 = self.corner
        return {
            "left": c0,
            "right": self.corrected(self.right.at(r0), r0)[0],
            "bottom": r0,
            "top": self.corrected(c0, self.top.at(c0))[1],
        }

    def interior(self, rows: np.ndarray, columns: np.ndarray) -> np.ndarray:
        """Return True for pixels at least ``FRAME_CLEAR`` inside the frame lines."""
        r, c = np.meshgrid(rows, columns, indexing="ij")
        m = FRAME_CLEAR
        return (
            (c > self.left.over(r) + m)
            & (c < self.right.over(r) - m)
            & (r > self.top.over(c) + m)
            & (r < self.bottom.over(c) - m)
        )


def _frame(dark: np.ndarray, panel: Panel) -> Frame:
    """Fit the four lines of a panel's frame."""
    left, right, top, bottom = panel.frame
    rows = range(top + 15, bottom - 15)
    columns = range(left + 15, right - 15)
    return Frame(
        _fit_line(dark, True, left, rows),
        _fit_line(dark, True, right, rows),
        _fit_line(dark, False, top, columns),
        _fit_line(dark, False, bottom, columns),
    )


@dataclass(frozen=True)
class Ticks:
    """The tick marks found: image (column, row) of each, and corrected positions."""

    bottom: list[tuple[float, float]]
    left: list[tuple[float, float]]


def _find_ticks(dark: np.ndarray, frame: Frame, panel: Panel) -> Ticks:
    """Return the dark stretches just inside the bottom and left frame lines."""
    inner, outer = TICK_INSET
    left, right, top, bottom = panel.frame
    columns = range(left + 15, right - 15)
    found = []
    for column in columns:
        row = round(frame.bottom.at(column))
        found.append(bool(dark[row - outer : row - inner, column].all()))
    bottom_ticks = []
    for a, b in _spans(np.array(found)):
        column = (a + b) / 2 + columns.start
        bottom_ticks.append((column, frame.bottom.at(column)))
    rows = range(top + 15, bottom - 15)
    found = []
    for row in rows:
        column = round(frame.left.at(row))
        found.append(bool(dark[row, column + inner : column + outer].all()))
    left_ticks = []
    for a, b in _spans(np.array(found)):
        row = (a + b) / 2 + rows.start
        left_ticks.append((frame.left.at(row), row))
    return Ticks(bottom_ticks, left_ticks)


def _match(found: Sequence[float], expected: Sequence[float], what: str) -> None:
    """Raise unless the found positions are the expected ones, each within tolerance."""
    if len(found) != len(expected) or any(
        abs(f - e) > TICK_TOLERANCE
        for f, e in zip(sorted(found), sorted(expected), strict=True)
    ):
        found_ = [round(f, 1) for f in sorted(found)]
        msg = f"{what}: found {found_}, expected {sorted(expected)}"
        raise ValueError(msg)


def _axes(dark: np.ndarray, panel: Panel) -> tuple[Frame, Ticks, Axis, Axis]:
    """Fit a panel's frame, check its ticks and return its calibrated axes."""
    frame = _frame(dark, panel)
    ticks = _find_ticks(dark, frame, panel)
    name = f"Fig. {panel.figure}"
    ends = frame.ends()
    for position, _, where in (*panel.x_ticks, *panel.y_ticks):
        if where != "tick" and abs(ends[where] - position) > TICK_TOLERANCE:
            msg = f"{name}: the {where} line at {ends[where]:.1f}, not {position}"
            raise ValueError(msg)
    _match(
        [frame.corrected(c, r)[0] for c, r in ticks.bottom],
        [p for p, _, w in panel.x_ticks if w == "tick"] + list(panel.not_x_ticks),
        f"{name}: x ticks",
    )
    _match(
        [frame.corrected(c, r)[1] for c, r in ticks.left],
        [p for p, _, w in panel.y_ticks if w == "tick"] + list(panel.not_y_ticks),
        f"{name}: y ticks",
    )
    x_axis = Axis(
        tuple(p for p, _, _ in panel.x_ticks),
        tuple(v for _, v, _ in panel.x_ticks),
        tolerance=AXIS_TOLERANCE,
    )
    y_axis = Axis(
        tuple(p for p, _, _ in panel.y_ticks),
        tuple(v for _, v, _ in panel.y_ticks),
        log=True,
        tolerance=AXIS_TOLERANCE,
    )
    x_axis.check()
    y_axis.check()
    # The unlabelled frame lines are where the text puts the ends of the range.
    for where, end in zip(("left", "right"), panel.ends, strict=True):
        value = float(x_axis(ends[where]))
        if abs(value - end) > 0.001:
            msg = f"{name}: the {where} line at {value:.4f} Ry, not {end}"
            raise ValueError(msg)
    return frame, ticks, x_axis, y_axis


@dataclass(frozen=True)
class Ink:
    """A panel's ink inside the frame, frame, ticks, lettering and markers erased.

    ``curves`` maps ``hf``, ``rpae`` and ``ldrpa`` to boolean images of the box
    ``rows`` x ``columns``: the solid, dashed and dotted lines.
    """

    rows: range
    columns: range
    curves: dict[str, np.ndarray]
    interior: np.ndarray
    markers: list[tuple[int, int, int, int]]


def _erase_ticks(
    box: np.ndarray, origin: tuple[int, int], frame: Frame, ticks: Ticks, panel: Panel
) -> None:
    """Erase the tick marks (not the curves meeting the frame) from a box of ink."""
    top, left = origin
    tick_columns = [p for p, _, w in panel.x_ticks if w == "tick"]
    for column, row in ticks.bottom:
        if any(
            abs(frame.corrected(column, row)[0] - t) <= TICK_TOLERANCE
            for t in tick_columns
        ):
            c, r = round(column) - left, round(row) - top
            length = _tick_length(box[: r + 1, c][::-1])
            box[r - length : r + 1, c - TICK_CLEAR : c + TICK_CLEAR + 1] = False
    tick_rows = [p for p, _, w in panel.y_ticks if w == "tick"]
    for column, row in ticks.left:
        if any(
            abs(frame.corrected(column, row)[1] - t) <= TICK_TOLERANCE
            for t in tick_rows
        ):
            c, r = round(column) - left, round(row) - top
            length = _tick_length(box[r, c:])
            box[r - TICK_CLEAR : r + TICK_CLEAR + 1, c : c + length + 1] = False


def _erase_boxes(
    box: np.ndarray,
    origin: tuple[int, int],
    boxes: Sequence[tuple[int, int, int, int]],
    margin: int,
) -> None:
    """Erase boxes (top, left, bottom, right; image pixels) with a margin."""
    top, left = origin
    for t, west, b, east in boxes:
        box[
            t - margin - top : b + margin - top,
            west - margin - left : east + margin - left,
        ] = False


def _ink(dark: np.ndarray, frame: Frame, ticks: Ticks, panel: Panel) -> Ink:
    """Erase the frame, ticks, lettering and markers; sort the ink into the curves."""
    left, right, top, bottom = panel.frame
    rows = range(top - 20, bottom + 20)
    columns = range(left - 20, right + 20)
    origin = (rows.start, columns.start)
    box = dark[rows.start : rows.stop, columns.start : columns.stop].copy()
    interior = frame.interior(
        np.arange(rows.start, rows.stop), np.arange(columns.start, columns.stop)
    )
    box &= interior
    _erase_ticks(box, origin, frame, ticks, panel)
    lettering = [*panel.labels, *([panel.cross] if panel.cross else [])]
    for t, west, b, east in lettering:
        if (
            box[
                t - rows.start : b - rows.start,
                west - columns.start : east - columns.start,
            ].sum()
            < 200
        ):
            msg = f"Fig. {panel.figure}: no lettering in the box {(t, west, b, east)}"
            raise ValueError(msg)
    _erase_boxes(box, origin, lettering, 0)
    markers = _markers(dark, panel)
    _erase_boxes(box, origin, markers, MARKER_MARGIN)
    labels, _ = ndimage.label(box, structure=np.ones((3, 3), dtype=int))
    curves = {name: np.zeros_like(box) for name in CURVES}
    for index, where in enumerate(ndimage.find_objects(labels), start=1):
        if where is None:
            continue
        size = max(where[0].stop - where[0].start, where[1].stop - where[1].start)
        piece = labels[where] == index
        if size >= SOLID:
            curves["hf"][where] |= piece
        elif size <= DOT:
            if piece.sum() >= 4:
                curves["ldrpa"][where] |= piece
        elif piece.sum() >= 10:
            curves["rpae"][where] |= piece
    return Ink(rows, columns, curves, interior, markers)


def _markers(dark: np.ndarray, panel: Panel) -> list[tuple[int, int, int, int]]:
    """Return the boxes (top, left, bottom, right) of the filled markers of a panel."""
    if panel.marker_box is None:
        return []
    t, west, b, east = panel.marker_box
    opened = ndimage.binary_opening(
        dark[t:b, west:east], structure=np.ones((MARKER, MARKER), dtype=bool)
    )
    labels, count = ndimage.label(opened)
    if count != panel.markers:
        msg = f"Fig. {panel.figure}: {count} markers, not {panel.markers}"
        raise ValueError(msg)
    return [
        (s[0].start + t, s[1].start + west, s[0].stop + t, s[1].stop + west)
        for s in ndimage.find_objects(labels)
        if s is not None
    ]


def _tick_length(line: np.ndarray) -> int:
    """Return how far a tick runs from the frame: the dark pixels from the first one."""
    dark = np.flatnonzero(line)
    if dark.size == 0:
        return 0
    first = int(dark[0])
    white = np.flatnonzero(~line[first:])
    return first + (int(white[0]) if white.size else len(line) - first) + 2


def _column_of(axis: Axis, value: float) -> float:
    """Return the corrected column of a photon energy."""
    positions = np.arange(0.0, float(SHAPE[1]))
    return float(np.interp(value, np.asarray(axis(positions), dtype=float), positions))


def _run(mask: np.ndarray, column: int, ink: Ink) -> float | None:
    """Return the centre row of a curve's single short run in an image column."""
    c = column - ink.columns.start
    if not 0 <= c < mask.shape[1]:
        return None
    spans = _spans(mask[:, c])
    if len(spans) != 1 or spans[0][1] - spans[0][0] + 1 > RUN:
        return None
    # A run touching the erased frame is cut by it.
    first, last = spans[0]
    if not (ink.interior[first - 1, c] and ink.interior[last + 1, c]):
        return None
    return (spans[0][0] + spans[0][1]) / 2 + ink.rows.start


def _bridge(mask: np.ndarray, column: int, ink: Ink) -> tuple[float | None, str]:
    """Interpolate a dashed curve across a gap between its dashes.

    Only across a true gap: the nearest readable columns within ``GAP`` pixels either
    side, with no ink at all in the columns between them (an unreadable steep stretch
    or peak is not bridged).
    """

    def empty(c: int) -> bool:
        return not mask[:, c - ink.columns.start].any()

    def nearest(step: int) -> int | None:
        for c in range(column + step, column + step * (GAP + 1), step):
            if _run(mask, c, ink) is not None:
                return c
            if not empty(c):
                return None
        return None

    if not empty(column):
        return None, ""
    left, right = nearest(-1), nearest(1)
    if left is None or right is None:
        return None, ""
    r_left, r_right = _run(mask, left, ink), _run(mask, right, ink)
    if r_left is None or r_right is None:
        return None, ""
    row = r_left + (r_right - r_left) * (column - left) / (right - left)
    return row, f"between pixels ({left}, {r_left:.1f}) and ({right}, {r_right:.1f})"


@dataclass(frozen=True)
class Sample:
    """One value read from a curve: where in the image and what it is."""

    grid: float
    column: float
    row: float
    energy: float
    sigma: float
    how: str


def _sample(
    mask: np.ndarray,
    ink: Ink,
    calibration: tuple[Frame, Axis, Axis],
    energy: float,
    row_guess: float,
    bridge: bool,
) -> Sample | None:
    """Read a curve at a photon energy, or None where it cannot be read.

    The column is that of the energy at the curve's row (found by one iteration from a
    guess); a gap of the dashes or dots is bridged (``bridge``) by the nearest readable
    columns within ``GAP`` pixels either side.
    """
    frame, x_axis, y_axis = calibration
    corrected = _column_of(x_axis, energy)
    column = round(frame.column(corrected, row_guess))
    for _ in range(2):
        row = _run(mask, column, ink)
        how = ""
        if row is None and bridge:
            row, how = _bridge(mask, column, ink)
        if row is None:
            return None
        new = round(frame.column(corrected, row))
        if new == column:
            break
        column = new
    else:
        return None
    if not how:
        how = f"pixel ({column}, {row:.1f})"
    c, r = frame.corrected(column, row)
    return Sample(energy, column, row, float(x_axis(c)), float(y_axis(r)), how)


@dataclass(frozen=True)
class Reading:
    """A panel's samples, and what the checks need of it."""

    panel: Panel
    frame: Frame
    axes: tuple[Axis, Axis]
    ink: Ink
    samples: dict[str, list[Sample]]
    dark: np.ndarray


# A value is kept only if a neighbouring grid point was read too and the two agree:
# their rows differ by at most SLOPE pixels per pixel of column (the steepest stretch
# the RUN rule lets through) plus SLACK pixels. This drops a stray piece of another
# line read where the curve itself has a gap.
SLOPE = 1.5
SLACK = 4.0


def _consistent(samples: list[Sample], grid: Sequence[float]) -> list[Sample]:
    """Return the samples that agree with a sample at a neighbouring grid point."""
    index = {energy: k for k, energy in enumerate(grid)}
    by_index = {index[s.grid]: s for s in samples}

    def agrees(a: Sample, b: Sample) -> bool:
        return abs(a.row - b.row) <= SLOPE * abs(a.column - b.column) + SLACK

    kept = []
    for k, sample in sorted(by_index.items()):
        neighbours = [by_index[j] for j in (k - 1, k + 1) if j in by_index]
        if any(agrees(sample, n) for n in neighbours):
            kept.append(sample)
    return kept


def _read(dark: np.ndarray, panel: Panel) -> Reading:
    """Return the samples of a panel's curves."""
    frame, ticks, x_axis, y_axis = _axes(dark, panel)
    ink = _ink(dark, frame, ticks, panel)
    out: dict[str, list[Sample]] = {}
    for name in CURVES:
        mask = ink.curves[name]
        samples = []
        guess = (ink.rows.start + ink.rows.stop) / 2
        for energy in panel.grid:
            if panel.ranges is not None and not (
                panel.ranges[name][0] <= energy <= panel.ranges[name][1]
            ):
                continue
            sample = _sample(
                mask, ink, (frame, x_axis, y_axis), energy, guess, name != "hf"
            )
            if sample is not None:
                samples.append(sample)
                guess = sample.row
        out[name] = _consistent(samples, panel.grid)
    return Reading(panel, frame, (x_axis, y_axis), ink, out, dark)


# McCown et al. (1982): alpha = 4e-32 cm^4 W^-1 at 193 nm (abstract and Table I);
# sigma = alpha hbar omega, with hbar omega = h c / 193 nm, in 1e-49 cm^4 s; the
# photon energy of 193 nm in Ry (R_infinity hc = 109 737.316 cm^-1).
MCCOWN_ALPHA = 4e-32
HC_J_CM = 1.98644586e-23  # h c, J cm (CODATA 2018, exact)
MCCOWN_SIGMA = MCCOWN_ALPHA * HC_J_CM / 193e-7 / 1e-49
MCCOWN_RY = 1e7 / 193 / 109737.316
# How far the cross may lie from 193 nm (Ry), and its height over MCCOWN_SIGMA as read
# (the record's notes).
CROSS_ENERGY = 0.002
CROSS_RATIO = (0.93, 1.07)


def _check_mccown(reading: Reading) -> float:
    """Raise unless the cross is at 193 nm; return its height over McCown et al.'s.

    The cross's centre is the centroid of the dark pixels in its box (the cross
    alone there).
    """
    if reading.panel.cross is None:
        msg = f"Fig. {reading.panel.figure}: no cross"
        raise ValueError(msg)
    top, left, bottom, right = reading.panel.cross
    rows, columns = np.nonzero(reading.dark[top:bottom, left:right])
    if rows.size < 50:
        msg = "Fig. 6: no cross in its box"
        raise ValueError(msg)
    column, row = float(columns.mean()) + left, float(rows.mean()) + top
    frame, (x_axis, y_axis) = reading.frame, reading.axes
    c, rr = frame.corrected(column, row)
    energy, sigma = float(x_axis(c)), float(y_axis(rr))
    ratio = sigma / MCCOWN_SIGMA
    problems = []
    if abs(energy - MCCOWN_RY) > CROSS_ENERGY:
        problems.append(f"at {energy:.4f} Ry, not {MCCOWN_RY:.4f}")
    if not CROSS_RATIO[0] <= ratio <= CROSS_RATIO[1]:
        problems.append(f"{ratio:.3f} times McCown et al.'s {MCCOWN_SIGMA:.3f}")
    if problems:
        msg = "Fig. 6, McCown et al.'s cross: " + "; ".join(problems)
        raise ValueError(msg)
    return ratio


def readings() -> list[Reading]:
    """Return every panel's reading, the checks passed."""
    grey = image(DIRECTORY, PAGE)
    if grey.shape != SHAPE:
        msg = f"page {PAGE}: an image of {grey.shape}, not {SHAPE}"
        raise ValueError(msg)
    dark = grey < DARK
    out = []
    for panel in PANELS:
        reading = _read(dark, panel)
        if panel.figure == "6":
            _check_mccown(reading)
        out.append(reading)
    return out


def extract() -> dict[str, list[Point]]:
    """Return the samples of Figs. 6 and 7, against the photon energy (Ry)."""
    files: dict[str, list[Point]] = {}
    for reading in readings():
        for name, samples in reading.samples.items():
            files[reading.panel.files[name]] = [
                Point(
                    f"{s.grid:.3f}",
                    f"{s.energy:.4f}",
                    f"{s.sigma:.3e}",
                    "",
                    s.how,
                )
                for s in samples
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
