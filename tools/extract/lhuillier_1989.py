# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract L'Huillier, Tang and Lambropoulos (1989): MQDT two-photon cross sections.

The two-photon cross sections of Kr and Xe are in figures only (data rule 8, TC-699).
The PDF is a scan: every page is one 1-bit image, 5169 x 7129 pixels (600 dpi), read
here from its native pixels (``_figure.image``).

* **Fig. 6** (p. 8, library figure ``fig_31bfc6fa4487``): the two-photon
  autoionisation spectra between the 2P3/2 and 2P1/2 limits, sigma_2 in 1e-50 cm^4 s
  on linear axes against the **total two-photon energy** in 1e3 cm^-1; (a) Kr,
  115.0-116.5, (b) Xe, 103-106. The solid line is the total cross section, the
  dashed line the partial cross section to J = 0 final states.
* **Fig. 10** (p. 9, library figure ``fig_2942ab283ba9``): two-photon ionisation
  above the second limit, on logarithmic axes; (a) Xe, 110-130, (b) Kr, 120-140, in
  1e3 cm^-1 of total energy. Only the solid line is this paper's (the others are
  McGuire's, Gangopadhyay et al.'s and L'Huillier and Wendin's; not taken).

How each panel is read:

* **the frame**: its four lines are fitted to the image (the centre of the dark run
  near the expected position in every other row or column, outliers beyond 1.5
  pixels rejected). The scan is sheared (the vertical lines lean by up to 0.013 pixel
  per row in Fig. 10, the horizontal ones by 0.003): a pixel (column, row) is mapped to
  *corrected* coordinates, the column it would have on the bottom frame line and the
  row it would have on the left one, through the mean slopes of the frame's lines.
* **the ticks**: the inward tick marks are found in the image along the bottom and
  left frame lines (a row or column dark from 6 to 18 pixels inside the line) and
  must be where ``x_ticks`` and ``y_ticks`` write them, within ``TICK_TOLERANCE``
  pixels; dark stretches there that are not ticks (a curve meeting the frame) are
  written in ``not_ticks`` and must be found too. The frame lines are the end values
  of the x axis and of the linear y axes. The axes are calibrated from those
  positions (``_figure.Axis``) and every tick checked on the fitted scale.
* **the curves**: the frame lines and the tick marks are erased and the ink inside
  the frame split into connected pieces (8-connected). The **solid** curve is the
  pieces at least ``SOLID`` pixels long or high (the lettering, the dashes and the
  markers are smaller); the **dashed** J = 0 curve of Fig. 6 is the smaller pieces
  that lie wholly in the band near zero (``J0_BAND``). Each curve is sampled on a grid
  of total energies (``GRID``); at a sample's column the curve's ink must form a single
  run no taller than ``RUN`` pixels (a steep flank or a resonance's peak gives a taller
  run or two, and is not read), whose centre is the value. A dash's gap is bridged by
  the linear interpolation between the nearest readable columns within ``GAP`` pixels.
  A J = 0 run that touches the erased bottom frame (values below about 0.07e-50) is not
  read. Windows written in ``EXCLUDE`` (a marker on the solid line) are not read.

Every value's ``printed`` cell gives its pixel (column, row) in the page image and
the total energy at it; the data file's ``x`` is the **photon** wavenumber, half the
total two-photon energy (the record's decision). One pixel is 1.12 cm^-1 of total
energy and 0.011e-50 cm^4 s in Fig. 6(a), 2.25 cm^-1 and 0.017e-50 in Fig. 6(b), and
19 cm^-1 and 0.43 % in Fig. 10; the solid lines are 6-8 pixels wide, the dashes 5-6.

Checks (``_check_mccown``): the cross marking McCown et al.'s measurement in Fig. 6(b)
lies on the solid line, at a total energy within 40 cm^-1 of twice the ArF line centre
McCown et al. give (103 466 cm^-1); its height is compared (and reported) with
McCown et al.'s printed coefficient, 4e-32 cm^4 W^-1, as sigma = alpha hbar omega.
"""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from _figure import Axis, image
from _library import Point, write
from scipy import ndimage

DIRECTORY = "LHuillier_1989_MultiphotonIonizationRare"
SHAPE = (7129, 5169)
DARK = 128
# How far, pixels, inside a frame line a tick is looked for (a row or column dark over
# this whole stretch), how far from a tick it may be found, and how much of the frame
# and of each tick (across it) is erased before the curves are read.
TICK_INSET = (6, 18)
TICK_TOLERANCE = 2.0
FRAME_CLEAR = 5
TICK_CLEAR = 5
# A piece of ink this long or high (pixels) is part of a solid curve; lettering,
# dashes and markers are smaller.
SOLID = 150
# The tallest run (pixels) a value is read from: the solid lines are 6-8 pixels wide,
# the dashes 5-6; a taller run is a steep stretch, a peak or two lines together.
RUN = 12
# How far (pixels) either side a sample in a gap of the dashes looks for the line.
GAP = 25


@dataclass(frozen=True)
class Panel:
    """One panel of a figure: where it is, its ticks and how its curves are read.

    Tick positions are corrected pixel coordinates (the module's docstring): a column
    for the x axis, a row for the y axis. ``where`` is ``left``, ``right``, ``bottom``
    or ``top`` for a frame line that carries an end value, or ``tick``.
    """

    figure: str
    page: int
    gas: str
    frame: tuple[int, int, int, int]
    """The frame lines' approximate columns (left, right) and rows (top, bottom)."""
    x_ticks: tuple[tuple[float, float, str], ...]
    """(position, total energy in 1e3 cm^-1, where)."""
    y_ticks: tuple[tuple[float, float, str], ...]
    """(position, sigma_2 in the axis's unit, where)."""
    log: bool
    grid: tuple[float, ...]
    """The total energies (1e3 cm^-1) where the curves are sampled."""
    files: dict[str, str]
    """The data file of each curve read (``total``, and ``j0`` in Fig. 6)."""
    not_x_ticks: tuple[float, ...] = ()
    """Corrected columns of dark stretches at the bottom line that are not ticks."""
    not_y_ticks: tuple[float, ...] = ()
    """Corrected rows of dark stretches at the left line that are not ticks."""
    exclude: tuple[tuple[float, float], ...] = ()
    """Windows of total energy (1e3 cm^-1) where the solid curve is not read."""
    j0_band: float = 0.0
    """Fig. 6: the J = 0 curve's dashes lie wholly below this sigma_2."""


def _grid(first: float, last: float, step: float) -> tuple[float, ...]:
    """Return the grid from first to last (inclusive) in steps, rounded to 3 places."""
    count = round((last - first) / step)
    return tuple(round(first + step * k, 3) for k in range(count + 1))


# The x ticks of Fig. 6 lie 2-5 pixels (3-11 cm^-1) off the straight line through the
# frame lines and each other: the residuals of the fit, within the tolerance below.
FIG6_X_TOLERANCE = 4.0
# Fig. 6(b)'s y tick at 5 lies 3 pixels off the line through the frame and the 10 tick.
FIG6_Y_TOLERANCE = 3.5
# The minor ticks of Fig. 10's logarithmic axes are hand-placed: they scatter about
# the fitted scale by up to 6.2 pixels (2.7 %); its x ticks by up to 3 pixels.
FIG10_Y_TOLERANCE = 6.5
FIG10_X_TOLERANCE = 3.5

PANELS = (
    Panel(
        figure="6",
        page=8,
        gas="Kr",
        frame=(916, 2251, 775, 1676),
        x_ticks=(
            (917.1, 115.0, "left"),
            (1366.5, 115.5, "tick"),
            (1805.5, 116.0, "tick"),
            (2251.5, 116.5, "right"),
        ),
        y_ticks=((1674.2, 0.0, "bottom"), (1225.5, 5.0, "tick"), (774.5, 10.0, "top")),
        log=False,
        # Below 115.01 and above 116.49 the curves run into the frame lines.
        grid=_grid(115.01, 116.49, 0.01),
        files={
            "total": "Kr/lhuillier_1989__fig6_total_Kr.csv",
            "j0": "Kr/lhuillier_1989__fig6_J0_Kr.csv",
        },
        # The solid line (6.0e-50) and the J = 0 dashes (0.2e-50) meet the left line.
        not_y_ticks=(1136.5, 1655.5),
        j0_band=1.0,
    ),
    Panel(
        figure="6",
        page=8,
        gas="Xe",
        frame=(916, 2253, 1981, 2877),
        x_ticks=(
            (917.3, 103.0, "left"),
            (1369.0, 104.0, "tick"),
            (1810.5, 105.0, "tick"),
            (2254.5, 106.0, "right"),
        ),
        y_ticks=(
            (2876.6, 0.0, "bottom"),
            (2574.5, 5.0, "tick"),
            (2279.5, 10.0, "tick"),
            (1980.2, 15.0, "top"),
        ),
        log=False,
        grid=_grid(103.02, 105.98, 0.02),
        files={
            "total": "Xe/lhuillier_1989__fig6_total_Xe.csv",
            "j0": "Xe/lhuillier_1989__fig6_J0_Xe.csv",
        },
        # The steep flanks of the J = 0 peaks at 104.6 and 105.7 cross the bottom
        # line's tick zone.
        not_x_ticks=(1633.0, 2117.5),
        # The solid line (6.6e-50) and the J = 0 dashes meet the left line.
        not_y_ticks=(2481.5, 2860.0),
        # McCown et al.'s cross lies on the solid line here (the check).
        exclude=((103.55, 103.69),),
        j0_band=3.0,
    ),
    Panel(
        figure="10",
        page=9,
        gas="Xe",
        frame=(3242, 4296, 3387, 4277),
        x_ticks=(
            (3246.0, 110.0, "left"),
            (3769.5, 120.0, "tick"),
            (4302.1, 130.0, "right"),
        ),
        y_ticks=(
            (3409.0, 2e-49, "tick"),
            (3576.5, 1e-49, "tick"),
            (3601.0, 9e-50, "tick"),
            (3631.0, 8e-50, "tick"),
            (3661.5, 7e-50, "tick"),
            (3693.5, 6e-50, "tick"),
            (3742.5, 5e-50, "tick"),
            (3791.0, 4e-50, "tick"),
            (3956.0, 2e-50, "tick"),
            (4122.5, 1e-50, "tick"),
            (4143.5, 9e-51, "tick"),
            (4167.0, 8e-51, "tick"),
            (4195.0, 7e-51, "tick"),
            (4237.5, 6e-51, "tick"),
        ),
        log=True,
        grid=_grid(110.5, 129.5, 0.5),
        files={"total": "Xe/lhuillier_1989__fig10_total_Xe.csv"},
        # McGuire's dashed line (3544), this paper's solid line (3853.5, over the
        # 3e-50 tick) and the linear-RPAE short dashes (3900.5) meet the left line.
        not_y_ticks=(3544.0, 3853.5, 3900.5),
        # Gangopadhyay et al.'s cross at 129 lies on the solid line.
        exclude=((128.6, 129.4),),
    ),
    Panel(
        figure="10",
        page=9,
        gas="Kr",
        frame=(3256, 4312, 4553, 5736),
        x_ticks=(
            (3261.9, 120.0, "left"),
            (3787.5, 130.0, "tick"),
            (4320.3, 140.0, "right"),
        ),
        y_ticks=(
            # The top line carries the label 10^-49.
            (4552.2, 1e-49, "top"),
            (4578.5, 9e-50, "tick"),
            (4608.0, 8e-50, "tick"),
            (4637.5, 7e-50, "tick"),
            (4672.0, 6e-50, "tick"),
            (4717.0, 5e-50, "tick"),
            # Merged with the start of McGuire's dashed line.
            (4771.0, 4e-50, "tick"),
            (4831.5, 3e-50, "tick"),
            (4932.5, 2e-50, "tick"),
            (5091.5, 1e-50, "tick"),
            (5118.5, 9e-51, "tick"),
            (5150.0, 8e-51, "tick"),
            (5177.0, 7e-51, "tick"),
            (5211.5, 6e-51, "tick"),
            (5262.0, 5e-51, "tick"),
            (5313.0, 4e-51, "tick"),
            (5376.0, 3e-51, "tick"),
            (5476.5, 2e-51, "tick"),
            (5633.0, 1e-51, "tick"),
        ),
        log=True,
        grid=_grid(120.5, 139.5, 0.5),
        files={"total": "Kr/lhuillier_1989__fig10_total_Kr.csv"},
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
    name = f"Fig. {panel.figure} {panel.gas}"
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
    six = panel.figure == "6"
    x_axis = Axis(
        tuple(p for p, _, _ in panel.x_ticks),
        tuple(v for _, v, _ in panel.x_ticks),
        tolerance=FIG6_X_TOLERANCE if six else FIG10_X_TOLERANCE,
    )
    y_axis = Axis(
        tuple(p for p, _, _ in panel.y_ticks),
        tuple(v for _, v, _ in panel.y_ticks),
        log=panel.log,
        tolerance=FIG6_Y_TOLERANCE if six else FIG10_Y_TOLERANCE,
    )
    x_axis.check()
    y_axis.check()
    return frame, ticks, x_axis, y_axis


@dataclass(frozen=True)
class Ink:
    """A panel's ink inside the frame, frame and ticks erased, split into curves.

    ``solid`` and ``dashed`` are boolean images of the box ``rows`` x ``columns``.
    """

    rows: range
    columns: range
    solid: np.ndarray
    dashed: np.ndarray
    interior: np.ndarray


def _ink(
    dark: np.ndarray, frame: Frame, ticks: Ticks, panel: Panel, y_axis: Axis
) -> Ink:
    """Erase the frame and ticks, and sort the remaining ink into the curves."""
    left, right, top, bottom = panel.frame
    rows = range(top - 20, bottom + 20)
    columns = range(left - 20, right + 20)
    box = dark[rows.start : rows.stop, columns.start : columns.stop].copy()
    interior = frame.interior(
        np.arange(rows.start, rows.stop), np.arange(columns.start, columns.stop)
    )
    box &= interior
    tick_columns = [p for p, _, w in panel.x_ticks if w == "tick"]
    for column, row in ticks.bottom:
        if any(
            abs(frame.corrected(column, row)[0] - t) <= TICK_TOLERANCE
            for t in tick_columns
        ):
            c, r = round(column) - columns.start, round(row) - rows.start
            length = _tick_length(box[: r + 1, c][::-1])
            box[r - length : r + 1, c - TICK_CLEAR : c + TICK_CLEAR + 1] = False
    tick_rows = [p for p, _, w in panel.y_ticks if w == "tick"]
    for column, row in ticks.left:
        if any(
            abs(frame.corrected(column, row)[1] - t) <= TICK_TOLERANCE
            for t in tick_rows
        ):
            c, r = round(column) - columns.start, round(row) - rows.start
            length = _tick_length(box[r, c:])
            box[r - TICK_CLEAR : r + TICK_CLEAR + 1, c : c + length + 1] = False
    labels, _ = ndimage.label(box, structure=np.ones((3, 3), dtype=int))
    slices = ndimage.find_objects(labels)
    solid = np.zeros_like(box)
    dashed = np.zeros_like(box)
    band_row = None
    if panel.j0_band:
        band_row = _row_of(y_axis, panel.j0_band) - rows.start
    for index, where in enumerate(slices, start=1):
        if where is None:
            continue
        height = where[0].stop - where[0].start
        width = where[1].stop - where[1].start
        piece = labels[where] == index
        if max(height, width) >= SOLID:
            solid[where] |= piece
        elif band_row is not None and where[0].start > band_row and piece.sum() >= 6:
            dashed[where] |= piece
    return Ink(rows, columns, solid, dashed, interior)


def _tick_length(line: np.ndarray) -> int:
    """Return how far a tick runs from the frame: the dark pixels from the first one."""
    dark = np.flatnonzero(line)
    if dark.size == 0:
        return 0
    first = int(dark[0])
    white = np.flatnonzero(~line[first:])
    return first + (int(white[0]) if white.size else len(line) - first) + 2


def _row_of(axis: Axis, value: float) -> float:
    """Return the corrected position of a value on an axis (its inverse)."""
    positions = np.arange(0.0, float(SHAPE[0]))
    values = np.asarray(axis(positions), dtype=float)
    order = np.argsort(values)
    return float(np.interp(value, values[order], positions[order]))


def _column_of(axis: Axis, value: float) -> float:
    """Return the corrected column of a total energy."""
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
    """Read a curve at a total energy, or None where it cannot be read.

    The column is that of the energy at the curve's row (found by one iteration from a
    guess); a gap of the dashes is bridged (``bridge``) by the nearest readable columns
    within ``GAP`` pixels either side.
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


def _read(dark: np.ndarray, panel: Panel) -> Reading:
    """Return the samples of a panel's curves."""
    frame, ticks, x_axis, y_axis = _axes(dark, panel)
    ink = _ink(dark, frame, ticks, panel, y_axis)
    curves = {"total": (ink.solid, False)}
    if panel.j0_band:
        curves["j0"] = (ink.dashed, True)
    out: dict[str, list[Sample]] = {}
    for name, (mask, bridge) in curves.items():
        samples = []
        guess = (ink.rows.start + ink.rows.stop) / 2
        for energy in panel.grid:
            if name == "total" and any(a <= energy <= b for a, b in panel.exclude):
                continue
            sample = _sample(mask, ink, (frame, x_axis, y_axis), energy, guess, bridge)
            if sample is not None:
                samples.append(sample)
                guess = sample.row
        out[name] = samples
    return Reading(panel, frame, (x_axis, y_axis), ink, out)


# McCown et al. (1982): alpha = 4e-32 cm^4 W^-1 at 193 nm (abstract and Table I).
# sigma = alpha hbar omega, with hbar omega = h c / 193 nm. The cross is drawn at twice
# the photon wavenumber of 193 nm (103.627e3 cm^-1), not at twice McCown et al.'s ArF
# line centre (103 466 cm^-1, their p. 3320).
MCCOWN_ALPHA = 4e-32
MCCOWN_TOTAL = 2 * 1e7 / 193 / 1e3
HC_J_CM = 1.98644586e-23  # h c, J cm (CODATA 2018, exact)
MCCOWN_SIGMA = MCCOWN_ALPHA * HC_J_CM / 193e-7 / 1e-50  # 1e-50 cm^4 s
# How far the cross is drawn from McCown et al.'s total energy (1e3 cm^-1), and from
# the solid line (relative); and its height over MCCOWN_SIGMA as read (the record's
# notes).
CROSS_ENERGY = 0.040
CROSS_ON_LINE = 0.05
CROSS_RATIO = (1.30, 1.42)


def _cross(reading: Reading) -> tuple[float, float, float, float]:
    """Return the cross's centre in Fig. 6(b): (column, row, total energy, sigma_2).

    The solid line through the cross's window is interpolated between the samples on
    either side; the cross's arms are the solid ink in the window farther than 5 pixels
    from that line, within 40 pixels of it, and the centre is their centroid (the arms
    are symmetric about it; the part hidden in the line is symmetric too).
    """
    ((low, high),) = reading.panel.exclude
    total = reading.samples["total"]
    before = max((s for s in total if s.energy < low), key=lambda s: s.energy)
    after = min((s for s in total if s.energy > high), key=lambda s: s.energy)
    ink = reading.ink
    found = []
    for column in range(round(before.column) + 1, round(after.column)):
        line = before.row + (after.row - before.row) * (column - before.column) / (
            after.column - before.column
        )
        c = column - ink.columns.start
        for row in range(round(line) - 40, round(line) + 41):
            if abs(row - line) > 5 and ink.solid[row - ink.rows.start, c]:
                found.append((column, row))
    centre = np.mean(np.array(found, dtype=float), axis=0)
    c, r = reading.frame.corrected(float(centre[0]), float(centre[1]))
    x_axis, y_axis = reading.axes
    return float(centre[0]), float(centre[1]), float(x_axis(c)), float(y_axis(r))


def _check_mccown(reading: Reading) -> float:
    """Raise unless the cross is where McCown et al.'s point belongs; return its ratio.

    It must lie at McCown et al.'s total energy (within ``CROSS_ENERGY``) and on the
    solid line (within ``CROSS_ON_LINE``); its height over McCown et al.'s printed
    coefficient as sigma = alpha hbar omega must stay what the record states.
    """
    _, _, energy, sigma = _cross(reading)
    total = reading.samples["total"]
    before = max((s for s in total if s.energy < energy), key=lambda s: s.energy)
    after = min((s for s in total if s.energy > energy), key=lambda s: s.energy)
    line = before.sigma + (after.sigma - before.sigma) * (energy - before.energy) / (
        after.energy - before.energy
    )
    ratio = sigma / MCCOWN_SIGMA
    problems = []
    if abs(energy - MCCOWN_TOTAL) > CROSS_ENERGY:
        problems.append(f"at {energy:.4f}e3 cm-1, not {MCCOWN_TOTAL}")
    if abs(sigma / line - 1) > CROSS_ON_LINE:
        problems.append(f"{sigma:.3f} off the line's {line:.3f}")
    if not CROSS_RATIO[0] <= ratio <= CROSS_RATIO[1]:
        problems.append(f"{ratio:.3f} times McCown et al.'s {MCCOWN_SIGMA:.3f}")
    if problems:
        msg = "Fig. 6(b), McCown et al.'s cross: " + "; ".join(problems)
        raise ValueError(msg)
    return ratio


def readings() -> list[Reading]:
    """Return every panel's reading, the checks passed."""
    pages: dict[int, np.ndarray] = {}
    out = []
    for panel in PANELS:
        if panel.page not in pages:
            grey = image(DIRECTORY, panel.page)
            if grey.shape != SHAPE:
                msg = f"page {panel.page}: an image of {grey.shape}, not {SHAPE}"
                raise ValueError(msg)
            pages[panel.page] = grey < DARK
        reading = _read(pages[panel.page], panel)
        if panel.figure == "6" and panel.gas == "Xe":
            _check_mccown(reading)
        out.append(reading)
    return out


def extract() -> dict[str, list[Point]]:
    """Return the samples of Figs. 6 and 10, against the photon wavenumber."""
    files: dict[str, list[Point]] = {}
    for reading in readings():
        six = reading.panel.figure == "6"
        for name, samples in reading.samples.items():
            files[reading.panel.files[name]] = [
                Point(
                    f"{s.grid:.2f}" if six else f"{s.grid:.1f}",
                    # The photon wavenumber: half the total two-photon energy.
                    f"{s.energy * 1e3 / 2:.1f}" if six else f"{s.energy * 1e3 / 2:.0f}",
                    f"{s.sigma:.3f}" if six else f"{s.sigma:.3e}",
                    "",
                    f"{s.how}; total energy {s.energy:.4f}e3 cm-1",
                )
                for s in samples
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
