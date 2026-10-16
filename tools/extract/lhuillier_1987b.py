# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract L'Huillier and Wendin (1987, J. Phys. B): two-photon cross sections of Xe.

The two-photon one-electron ionisation cross sections of the Xe 5p shell are in
figures only (data rule 8, TC-699). The PDF is a scan: every page is one 1-bit image,
2076 x 2945 pixels (300 dpi), read here from its native pixels (``_figure.image``).
PDF page 2 is journal p. L37.

* **Fig. 2** (p. L40, PDF page 5, library figure ``fig_171e5219ba1c``): sigma_2 in
  1e-49 cm^4 s on a linear axis (0 to 4.5) against the photon energy, 0.31 to 1.40 Ry,
  LDA basis: the independent-electron approximation (solid), the LDRPA with every
  electron-photon interaction screened (dashed) and with only the first screened
  (dotted).
* **Fig. 3** (p. L41, PDF page 6, ``fig_fea006c5adce``): the same on an HF basis, 0.46
  to 0.86 Ry, logarithmic (1e-2 to 1e5): HF (solid), RPAE fully screened (dashed) and
  RPAE with only the first interaction screened (dotted). Only the dotted line is
  taken: the solid and dashed lines are the calculations lhuillier_1987 (TC-883) holds
  from its Fig. 6 (checked by overlay; the record's notes).
* **Fig. 4** (p. L42, PDF page 7, ``fig_1a8946c9a5c6``): RPAE (solid) and LDRPA
  (dashed), both fully screened, 0.31 to 1.40 Ry, logarithmic (1e-2 to 1e5), with
  Gangopadhyay et al.'s squares and McCown et al.'s cross (not taken).

How each panel is read (the method of lhuillier_1987, scaled to 300 dpi):

* **the frame**: its four lines are fitted to the image (the centre of the dark run
  near the expected position in every other row or column, outliers beyond 1.5 pixels
  rejected); a pixel (column, row) is mapped to *corrected* coordinates, the column it
  would have on the bottom frame line and the row it would have on the left one,
  through the mean slopes of the frame's lines.
* **the ticks**: inward tick marks lie on all four lines. Those of the calibrating
  sides (``x_side``, ``y_side``: Fig. 2's top and right lines, which no curve meets;
  the bottom and left lines elsewhere) are found as columns (rows) dark from
  ``TICK_INSET[0]`` to ``TICK_INSET[1]`` pixels inside the line and must be where
  ``x_ticks`` and ``y_ticks`` write them, within ``TICK_TOLERANCE``; dark stretches
  there that are not ticks (a curve meeting the frame) are written in ``not_x_ticks``
  and ``not_y_ticks`` and must be found too. Labelled frame lines are end values (the
  bottom line: 0 in Fig. 2, 1e-2 in Figs. 3 and 4). The axes are calibrated from
  those positions (``_figure.Axis``), every tick within ``AXIS_TOLERANCE`` of the
  fitted scale. Fig. 3's left and right lines must lie at the energies the text gives
  (0.46 and 0.86 Ry) within ``ENDS``.
* **the curves**: the frame lines and the ticks found on all four sides are erased,
  then the lettering and hatched threshold marks (boxes written in ``labels``),
  McCown et al.'s cross (``cross``) and Fig. 4's filled squares (the cores an erosion
  by a ``MARKER`` x ``MARKER`` square leaves inside ``marker_box``, grown back; their
  count is checked), with a margin of ``MARGIN`` pixels. The remaining ink is split into
  connected pieces (8-connected): pieces at least ``SOLID`` pixels long or high are the
  **solid** line, pieces no larger than ``DOT`` pixels the **dotted** line, the rest
  the **dashed** line (smaller specks are dust). Each curve is sampled on a grid of
  photon energies (``grid``) within its ``ranges``; at a sample's column the curve's
  ink must form a single run no taller than ``RUN`` pixels (a steep flank, a peak, or
  two pieces in one column are not read), whose centre is the value. A gap between
  dashes or dots is bridged by the linear interpolation between the nearest readable
  columns within ``GAP`` pixels, if no ink of that curve lies between. A value is kept
  only if a neighbouring grid point was read and agrees with it (``SLOPE``,
  ``SLACK``).

Every value's ``printed`` cell gives its pixel (column, row) in the page image. One
pixel is 0.00097 Ry and 0.0054 (1e-49 cm^4 s) in Fig. 2, 0.00035 Ry and 1.9 % in Fig.
3, 0.00093 Ry and 1.8 % in Fig. 4; the lines are 3-4 pixels wide.

**Checks** (:func:`checks`; each raises if missed): Fig. 4's six squares lie at the
photon energies (half the total) of the first six rows of gangopadhyay_1986's Table I
and at the sums of its two columns, within ``SQUARE_ENERGY`` and ``SQUARE_VALUE``;
Fig. 4's cross lies at the photon energy of 193 nm (the text: 'at 0.47 Ryd (193 nm)')
within ``CROSS_ENERGY``; Fig. 4's LDRPA (log) agrees with Fig. 2's LDRPA (linear), the
same calculation drawn twice, within ``REDRAWN`` where Fig. 2's value exceeds 0.3; Fig.
3's solid and dashed lines (read for this only) agree with lhuillier_1987's Fig. 6 HF
and RPAE within ``REPEATED``. Reported: the crosses' heights over McCown et al.'s alpha
hbar omega; Fig. 3's cross (drawn at 0.484 Ry, not 0.472); the hatched threshold
marks; Fig. 4's RPAE and LDRPA over lhuillier_1987's Figs. 6 and 7.
"""

import csv
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from _figure import Axis, image
from _library import DATA, Point, write
from scipy import ndimage

DIRECTORY = "LHuillier_1987_TwoPhotonOne"
SHAPE = (2945, 2076)
DARK = 128
TICK_INSET = (4, 8)
"""A tick is dark over this whole stretch inside its frame line, pixels (the minor
ticks are 9-10 pixels long, the major 15)."""
TICK_TOLERANCE = 2.0
"""A tick lies within this of where it is written, corrected pixels."""
AXIS_TOLERANCE = 1.5
"""The ticks' largest residual from the fitted scales, pixels."""
ENDS = 0.002
"""Fig. 3's frame lines lie at the text's 0.46 and 0.86 Ry within this, Ry (two
pixels)."""
FRAME_CLEAR = 3
TICK_CLEAR = 3
"""How much of the frame lines and of each tick (either side) is erased, pixels."""
SOLID = 75
DOT = 7
DUST = 8
"""A piece this long or high (pixels) is the solid line, one no larger than DOT a dot,
one with fewer than DUST pixels dust; the dashes are 19-28 pixels long, the dots 5-7
pixels across."""
MARKER = 12
MARGIN = 2
"""Fig. 4's squares (17-18 pixels) keep a core after an erosion by a MARKER-pixel
square, the lines (3-4 pixels wide away from the resonances) do not; a square, the
cross and the lettering are erased with MARGIN pixels around them."""
RUN = 7
"""The tallest run (pixels) a value is read from: the lines are 3-4 pixels wide."""
GAP = 14
"""How far (pixels) either side a sample in a gap of the dashes or dots looks."""
SLOPE = 1.5
SLACK = 3.0
"""Neighbouring values agree if their rows differ by at most SLOPE pixels per pixel of
column plus SLACK pixels."""
STYLES = ("solid", "dashed", "dotted")
SIDES = ("left", "right", "top", "bottom")


@dataclass(frozen=True)
class Panel:
    """One figure: where it is, its ticks and how its curves are read.

    Tick positions are corrected pixel coordinates: a column for the x axis, a row for
    the y axis. ``where`` is ``bottom`` for the labelled bottom frame line, or
    ``tick``.
    """

    figure: str
    page: int
    frame: tuple[int, int, int, int]
    """The frame lines' approximate columns (left, right) and rows (top, bottom)."""
    x_side: str
    y_side: str
    x_ticks: tuple[tuple[float, float, str], ...]
    """(position, photon energy in Ry, where)."""
    y_ticks: tuple[tuple[float, float, str], ...]
    """(position, sigma_2 in 1e-49 cm^4 s, where)."""
    log: bool
    grid: tuple[float, ...]
    files: dict[str, str]
    """The data file of each line style taken."""
    ranges: dict[str, tuple[float, ...]]
    """Where each style is read: pairs of photon energies (Ry), low and high."""
    ends: tuple[float, float] | None = None
    """The photon energies (Ry) the text gives for the left and right frame lines."""
    not_x_ticks: tuple[float, ...] = ()
    not_y_ticks: tuple[float, ...] = ()
    labels: tuple[tuple[int, int, int, int], ...] = ()
    """Boxes of lettering and hatched marks (top, left, bottom, right; image pixels)."""
    cross: tuple[int, int, int, int] | None = None
    """McCown et al.'s cross (top, left, bottom, right; image pixels)."""
    marker_box: tuple[int, int, int, int] | None = None
    markers: int = 0
    hatches: tuple[tuple[int, int, int, int], ...] = field(default=())
    """The hatched threshold marks among ``labels``, reported."""
    solid: int = SOLID
    """A piece this long or high (pixels) is the solid line."""


def _grid(first: float, last: float, step: float) -> tuple[float, ...]:
    """Return the grid from first to last (inclusive) in steps, rounded to 4 places."""
    count = round((last - first) / step)
    return tuple(round(first + step * k, 4) for k in range(count + 1))


def _file(key: str) -> str:
    return f"Xe/lhuillier_1987b__{key}_Xe.csv"


PANELS = (
    Panel(
        figure="2",
        page=5,
        frame=(594, 1715, 1612, 2444),
        x_side="top",
        y_side="right",
        # The top ticks every 0.05 Ry from 0.35 (labels 0.4 to 1.2 under the bottom
        # line; the top ticks lie 1-3 pixels left of the bottom ones: the scan's shear,
        # corrected).
        x_ticks=tuple(
            (p, round(0.35 + 0.05 * k, 2), "tick")
            for k, p in enumerate(
                (
                    *(635.5, 686.5, 738.5, 790.5, 841.5, 893.5, 944.5),
                    *(996.0, 1048.0, 1099.5, 1150.5, 1202.0, 1253.5, 1304.5),
                    *(1356.5, 1407.5, 1459.0, 1510.5, 1561.5, 1613.5, 1664.5),
                )
            )
        ),
        # The right ticks every 0.25 (labels 0 to 4 on the left); the bottom line is 0.
        y_ticks=(
            (2445.8, 0.0, "bottom"),
            *(
                (p, 0.25 * (k + 1), "tick")
                for k, p in enumerate(
                    (
                        *(2399.1, 2353.1, 2307.1, 2261.1, 2215.1, 2169.1, 2123.1),
                        *(2077.1, 2030.1, 1984.1, 1938.1, 1891.6, 1845.1, 1799.1),
                        *(1752.6, 1706.1, 1659.6),
                    )
                )
            ),
        ),
        log=False,
        grid=_grid(0.312, 1.398, 0.002),
        files={
            "solid": _file("fig2_lda"),
            "dashed": _file("fig2_ldrpa"),
            "dotted": _file("fig2_ldrpa_first"),
        },
        ranges={
            "solid": (0.312, 1.398),
            "dashed": (0.312, 1.398),
            "dotted": (0.312, 1.398),
        },
        # The hatched threshold mark meets the top line; the dashed (0.30) and dotted
        # (0.07) lines meet the right line.
        not_x_ticks=(913.5,),
        not_y_ticks=(2391.1, 2435.1),
        labels=((1619, 907, 1763, 934),),  # the hatched LDA threshold
        hatches=((1619, 907, 1763, 934),),
    ),
    Panel(
        figure="3",
        page=6,
        frame=(609, 1738, 1198, 2033),
        x_side="bottom",
        y_side="left",
        x_ticks=tuple(
            (p, round(0.475 + 0.025 * k, 3), "tick")
            for k, p in enumerate(
                (
                    *(652.0, 722.5, 793.5, 863.5, 934.5, 1004.5, 1075.5, 1145.5),
                    *(1216.0, 1286.5, 1357.0, 1427.0, 1497.5, 1567.5, 1638.5),
                    1709.0,
                )
            )
        ),
        y_ticks=(
            (2033.2, 1e-2, "bottom"),
            (1913.5, 1e-1, "tick"),
            (1793.5, 1e0, "tick"),
            (1674.5, 1e1, "tick"),
            (1555.0, 1e2, "tick"),
            (1436.0, 1e3, "tick"),
            (1316.5, 1e4, "tick"),
        ),
        log=True,
        ends=(0.46, 0.86),
        grid=_grid(0.462, 0.858, 0.002),
        files={"dotted": _file("fig3_rpae_first")},
        # Above 0.79 Ry the resonances of the three lines overlap. The solid and
        # dashed lines are read for the check against lhuillier_1987 only.
        ranges={
            "solid": (0.462, 0.790),
            "dashed": (0.462, 0.765),
            "dotted": (0.462, 0.790),
        },
        # The solid (1.4) and dashed (0.26) lines meet the left line.
        not_y_ticks=(1775.5, 1864.5),
        labels=(
            (1470, 1156, 1496, 1187),  # 6s
            (1303, 1364, 1336, 1529),  # <- 5d ->
        ),
        cross=(1824, 667, 1852, 691),
    ),
    Panel(
        figure="4",
        page=7,
        frame=(613, 1787, 284, 1157),
        x_side="bottom",
        y_side="left",
        x_ticks=tuple(
            (p, round(0.35 + 0.05 * k, 2), "tick")
            for k, p in enumerate(
                (
                    *(655.5, 709.5, 763.5, 817.5, 871.5, 925.5, 979.5, 1033.5),
                    *(1087.5, 1141.5, 1194.5, 1248.5, 1302.5, 1356.0, 1409.5),
                    *(1463.5, 1517.0, 1571.0, 1624.5, 1679.0, 1732.5),
                )
            )
        ),
        y_ticks=(
            (1157.2, 1e-2, "bottom"),
            (1032.5, 1e-1, "tick"),
            (907.5, 1e0, "tick"),
            (781.5, 1e1, "tick"),
            (656.5, 1e2, "tick"),
            (531.5, 1e3, "tick"),
            (406.5, 1e4, "tick"),
        ),
        log=True,
        grid=_grid(0.312, 1.398, 0.002),
        files={"solid": _file("fig4_rpae"), "dashed": _file("fig4_ldrpa")},
        # The solid line: the threshold region below the 6s resonance, and the
        # above-threshold part; between, the resonances overlap.
        ranges={"solid": (0.462, 0.790, 0.976, 1.398), "dashed": (0.312, 1.398)},
        # The dashed line (0.48) meets the left line.
        not_y_ticks=(947.5,),
        labels=(
            (331, 1069, 357, 1099),  # 5d
            (332, 1186, 357, 1216),  # 7s
            (496, 951, 521, 982),  # 6s
            (469, 1218, 495, 1249),  # 6d
            (578, 1221, 604, 1251),  # 8s
            (698, 767, 858, 793),  # the hatched two-photon threshold
            (698, 1251, 858, 1277),  # the hatched one-photon threshold
        ),
        hatches=((698, 767, 858, 793), (698, 1251, 858, 1277)),
        cross=(944, 775, 968, 801),
        marker_box=(880, 815, 1045, 960),
        markers=6,
        # The squares cut the solid line into pieces 40-60 pixels long; the dashes
        # are at most 28.
        solid=40,
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

    The centre of the dark run (at most 7 pixels wide) nearest the expected position
    within 10 pixels, in every other row or column; outliers beyond 1.5 pixels (a tick
    or a curve meeting the line) are rejected, three times.
    """
    points = []
    for along in span[::2]:
        cut = slice(near - 10, near + 11)
        line = dark[along, cut] if vertical else dark[cut, along]
        runs = [(a + b) / 2 + near - 10 for a, b in _spans(line) if b - a < 7]
        if runs:
            points.append((along, min(runs, key=lambda c: abs(c - near))))
    data = np.array(points)
    for _ in range(3):
        slope, intercept = np.polyfit(data[:, 0], data[:, 1], 1)
        data = data[np.abs(data[:, 1] - (slope * data[:, 0] + intercept)) < 1.5]
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


def _side_ticks(dark: np.ndarray, frame: Frame, panel: Panel, side: str) -> list:
    """Return the image (column, row) of the dark stretches inside one frame line."""
    inner, outer = TICK_INSET
    left, right, top, bottom = panel.frame
    horizontal = side in {"top", "bottom"}
    span = range(left + 10, right - 10) if horizontal else range(top + 10, bottom - 10)
    line = getattr(frame, side)
    inward = 1 if side in {"top", "left"} else -1
    found = []
    for along in span:
        at = round(line.at(along))
        a, b = sorted((at + inward * inner, at + inward * outer))
        strip = dark[a:b, along] if horizontal else dark[along, a:b]
        found.append(bool(strip.all()))
    out = []
    for a, b in _spans(np.array(found)):
        along = (a + b) / 2 + span.start
        out.append((along, line.at(along)) if horizontal else (line.at(along), along))
    return out


def _match(found: Sequence[float], expected: Sequence[float], what: str) -> None:
    """Raise unless the found positions are the expected ones, each within tolerance."""
    if len(found) != len(expected) or any(
        abs(f - e) > TICK_TOLERANCE
        for f, e in zip(sorted(found), sorted(expected), strict=True)
    ):
        found_ = [round(f, 1) for f in sorted(found)]
        msg = f"{what}: found {found_}, expected {sorted(expected)}"
        raise ValueError(msg)


@dataclass(frozen=True)
class Calibration:
    """A panel's frame, its calibrated axes and the ticks found on every side."""

    frame: Frame
    x: Axis
    y: Axis
    ticks: dict[str, list]


def _calibrate(dark: np.ndarray, panel: Panel) -> Calibration:
    """Fit a panel's frame, check its ticks and return its calibrated axes."""
    frame = _frame(dark, panel)
    ticks = {s: _side_ticks(dark, frame, panel, s) for s in SIDES}
    name = f"Fig. {panel.figure}"
    ends = frame.ends()
    for position, _, where in (*panel.x_ticks, *panel.y_ticks):
        if where != "tick" and abs(ends[where] - position) > TICK_TOLERANCE:
            msg = f"{name}: the {where} line at {ends[where]:.1f}, not {position}"
            raise ValueError(msg)
    _match(
        [frame.corrected(c, r)[0] for c, r in ticks[panel.x_side]],
        [p for p, _, w in panel.x_ticks if w == "tick"] + list(panel.not_x_ticks),
        f"{name}: x ticks",
    )
    _match(
        [frame.corrected(c, r)[1] for c, r in ticks[panel.y_side]],
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
        log=panel.log,
        tolerance=AXIS_TOLERANCE,
    )
    x_axis.check()
    y_axis.check()
    if panel.ends is not None:
        for where, end in zip(("left", "right"), panel.ends, strict=True):
            value = float(x_axis(ends[where]))
            if abs(value - end) > ENDS:
                msg = f"{name}: the {where} line at {value:.4f} Ry, not {end}"
                raise ValueError(msg)
    return Calibration(frame, x_axis, y_axis, ticks)


def _tick_length(line: np.ndarray) -> int:
    """Return how far a tick runs from the frame: the dark pixels from the first one."""
    dark = np.flatnonzero(line)
    if dark.size == 0:
        return 0
    first = int(dark[0])
    white = np.flatnonzero(~line[first:])
    return first + (int(white[0]) if white.size else len(line) - first) + 2


@dataclass(frozen=True)
class Ink:
    """A panel's ink inside the frame, sorted into the three line styles."""

    rows: range
    columns: range
    styles: dict[str, np.ndarray]
    interior: np.ndarray
    """True inside the frame, away from the frame lines, ticks, lettering and
    markers: a run that touches anything else was cut by an erasure."""


def _erase(box: np.ndarray, origin: tuple[int, int], where: tuple[int, ...]) -> None:
    """Erase a box (top, left, bottom, right; image pixels) with ``MARGIN``."""
    top, left = origin
    t, west, b, east = where
    box[
        max(t - MARGIN - top, 0) : b + MARGIN - top,
        max(west - MARGIN - left, 0) : east + MARGIN - left,
    ] = False


def _cores(dark: np.ndarray, panel: Panel) -> list[tuple[int, int, int, int]]:
    """Return the cores (top, left, bottom, right) of a panel's filled squares.

    An erosion by a ``MARKER`` x ``MARKER`` square inside ``marker_box`` leaves a core
    of each square (17-18 pixels) and nothing of the lines (3-4 pixels wide), and
    separates squares that touch; their count is checked.
    """
    if panel.marker_box is None:
        return []
    t, west, b, east = panel.marker_box
    eroded = ndimage.binary_erosion(
        dark[t:b, west:east], structure=np.ones((MARKER, MARKER), dtype=bool)
    )
    labels, count = ndimage.label(eroded, structure=np.ones((3, 3), dtype=int))
    if count != panel.markers:
        msg = f"Fig. {panel.figure}: {count} squares, not {panel.markers}"
        raise ValueError(msg)
    return [
        (s[0].start + t, s[1].start + west, s[0].stop + t, s[1].stop + west)
        for s in ndimage.find_objects(labels)
        if s is not None
    ]


def markers(dark: np.ndarray, panel: Panel) -> list[tuple[int, int, int, int]]:
    """Return the boxes of a panel's filled squares: their cores grown back."""
    grow = MARKER // 2
    return [
        (t - grow, west - grow, b + grow, east + grow)
        for t, west, b, east in _cores(dark, panel)
    ]


def _erase_ticks(box: np.ndarray, origin: tuple[int, int], cal: Calibration) -> None:
    """Erase the tick marks found on the four frame lines from a box of ink."""
    top, left = origin
    for side, found in cal.ticks.items():
        for column, row in found:
            c, r = round(column) - left, round(row) - top
            if side == "bottom":
                n = _tick_length(box[: r + 1, c][::-1])
                box[r - n : r + 1, c - TICK_CLEAR : c + TICK_CLEAR + 1] = False
            elif side == "top":
                n = _tick_length(box[r:, c])
                box[r : r + n + 1, c - TICK_CLEAR : c + TICK_CLEAR + 1] = False
            elif side == "left":
                n = _tick_length(box[r, c:])
                box[r - TICK_CLEAR : r + TICK_CLEAR + 1, c : c + n + 1] = False
            else:
                n = _tick_length(box[r, : c + 1][::-1])
                box[r - TICK_CLEAR : r + TICK_CLEAR + 1, c - n : c + 1] = False


def _ink(dark: np.ndarray, cal: Calibration, panel: Panel) -> Ink:
    """Erase the frame, ticks, lettering and markers; sort the ink into the styles."""
    left, right, top, bottom = panel.frame
    rows = range(top - 20, bottom + 20)
    columns = range(left - 20, right + 20)
    origin = (rows.start, columns.start)
    box = dark[rows.start : rows.stop, columns.start : columns.stop].copy()
    interior = cal.frame.interior(
        np.arange(rows.start, rows.stop), np.arange(columns.start, columns.stop)
    )
    box &= interior
    before = box.copy()
    _erase_ticks(box, origin, cal)
    lettering = [*panel.labels, *([panel.cross] if panel.cross else [])]
    for t, west, b, east in lettering:
        ink = box[
            t - rows.start : b - rows.start, west - columns.start : east - columns.start
        ]
        if ink.sum() < 40:
            msg = f"Fig. {panel.figure}: no lettering in the box {(t, west, b, east)}"
            raise ValueError(msg)
    for where in (*lettering, *markers(dark, panel)):
        _erase(box, origin, where)
    interior &= ~(before & ~box)
    labels, _ = ndimage.label(box, structure=np.ones((3, 3), dtype=int))
    styles = {name: np.zeros_like(box) for name in STYLES}
    for index, where in enumerate(ndimage.find_objects(labels), start=1):
        if where is None:
            continue
        size = max(where[0].stop - where[0].start, where[1].stop - where[1].start)
        piece = labels[where] == index
        if piece.sum() < DUST:
            continue
        if size >= panel.solid:
            styles["solid"][where] |= piece
        elif size <= DOT:
            styles["dotted"][where] |= piece
        else:
            styles["dashed"][where] |= piece
    return Ink(rows, columns, styles, interior)


def _column_of(axis: Axis, value: float) -> float:
    """Return the corrected column of a photon energy."""
    positions = np.arange(0.0, float(SHAPE[1]))
    return float(np.interp(value, np.asarray(axis(positions), dtype=float), positions))


def _run(mask: np.ndarray, column: int, ink: Ink) -> float | None:
    """Return the centre row of a style's single short run in an image column."""
    c = column - ink.columns.start
    if not 0 <= c < mask.shape[1]:
        return None
    spans = _spans(mask[:, c])
    if len(spans) != 1 or spans[0][1] - spans[0][0] + 1 > RUN:
        return None
    first, last = spans[0]
    if not (ink.interior[first - 1, c] and ink.interior[last + 1, c]):
        return None
    return (first + last) / 2 + ink.rows.start


def _bridge(mask: np.ndarray, column: int, ink: Ink) -> tuple[float | None, str]:
    """Interpolate a dashed or dotted line across a gap between its pieces.

    Only across a true gap: the nearest readable columns within ``GAP`` pixels either
    side, with no ink at all of the style in the columns between them.
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
    """One value read from a line: where in the image and what it is."""

    grid: float
    column: float
    row: float
    energy: float
    sigma: float
    how: str


def _sample(
    mask: np.ndarray,
    ink: Ink,
    cal: Calibration,
    energy: float,
    row_guess: float,
    bridge: bool,
) -> Sample | None:
    """Read a line at a photon energy, or None where it cannot be read."""
    corrected = _column_of(cal.x, energy)
    column = round(cal.frame.column(corrected, row_guess))
    how = ""
    row = None
    for _ in range(2):
        row = _run(mask, column, ink)
        how = ""
        if row is None and bridge:
            row, how = _bridge(mask, column, ink)
        if row is None:
            return None
        new = round(cal.frame.column(corrected, row))
        if new == column:
            break
        column = new
    else:
        return None
    if row is None:
        return None
    if not how:
        how = f"pixel ({column}, {row:.1f})"
    c, r = cal.frame.corrected(column, row)
    return Sample(energy, column, row, float(cal.x(c)), float(cal.y(r)), how)


def _consistent(samples: list[Sample], grid: Sequence[float]) -> list[Sample]:
    """Return the samples that agree with a sample at a neighbouring grid point."""
    index = {energy: k for k, energy in enumerate(grid)}
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
    """A panel's calibration, ink and samples."""

    panel: Panel
    calibration: Calibration
    ink: Ink
    samples: dict[str, list[Sample]]
    dark: np.ndarray


def _inside(energy: float, ranges: tuple[float, ...]) -> bool:
    pairs = zip(ranges[::2], ranges[1::2], strict=True)
    return any(low <= energy <= high for low, high in pairs)


def read(panel: Panel) -> Reading:
    """Return a panel's reading, its ticks checked."""
    grey = image(DIRECTORY, panel.page)
    if grey.shape != SHAPE:
        msg = f"page {panel.page}: an image of {grey.shape}, not {SHAPE}"
        raise ValueError(msg)
    dark = grey < DARK
    cal = _calibrate(dark, panel)
    ink = _ink(dark, cal, panel)
    out: dict[str, list[Sample]] = {}
    for style, ranges in panel.ranges.items():
        mask = ink.styles[style]
        samples = []
        guess = (ink.rows.start + ink.rows.stop) / 2
        for energy in panel.grid:
            if not _inside(energy, ranges):
                continue
            sample = _sample(mask, ink, cal, energy, guess, style != "solid")
            if sample is not None:
                samples.append(sample)
                guess = sample.row
        out[style] = _consistent(samples, panel.grid)
    return Reading(panel, cal, ink, out, dark)


def _centre(reading: Reading, where: tuple[int, int, int, int]) -> tuple[float, float]:
    """Return the (photon energy, sigma_2) of the centroid of the ink in a box."""
    top, left, bottom, right = where
    rows, columns = np.nonzero(reading.dark[top:bottom, left:right])
    cal = reading.calibration
    c, r = cal.frame.corrected(float(columns.mean()) + left, float(rows.mean()) + top)
    return float(cal.x(c)), float(cal.y(r))


def _core_centre(
    reading: Reading, core: tuple[int, int, int, int]
) -> tuple[float, float]:
    """Return the (photon energy, sigma_2) of a square: the centre of its core."""
    top, left, bottom, right = core
    cal = reading.calibration
    c, r = cal.frame.corrected((left + right - 1) / 2, (top + bottom - 1) / 2)
    return float(cal.x(c)), float(cal.y(r))


def _hatch(reading: Reading, where: tuple[int, int, int, int]) -> float:
    """Return the photon energy (Ry) of a hatched mark's vertical bar.

    The bar is the column of the box darkest over the box's height.
    """
    top, left, bottom, right = where
    counts = reading.dark[top:bottom, left:right].sum(axis=0)
    column = left + float(np.flatnonzero(counts == counts.max()).mean())
    c, _ = reading.calibration.frame.corrected(column, (top + bottom) / 2)
    return float(reading.calibration.x(c))


RYDBERG_EV = 13.605693122994
"""R_inf h c in eV (CODATA 2018)."""
MCCOWN_RY = 1e7 / 193 / 109737.31568160
"""The photon energy of 193 nm in Ry (R_inf = 109 737.316 cm^-1, CODATA 2018)."""
HC_J_CM = 1.98644586e-23
"""h c in J cm (CODATA 2018, exact)."""
MCCOWN_SIGMA = 4e-32 * HC_J_CM / 193e-7 / 1e-49
"""mccown_1982's alpha = 4e-32 cm^4 W^-1 as sigma = alpha hbar omega, 1e-49 cm^4 s."""
CROSS_ENERGY = 0.004
"""Fig. 4's cross lies within this of 193 nm, Ry (four pixels and the text's
'0.47')."""
SQUARE_ENERGY = 0.007
SQUARE_VALUE = 0.10
"""Fig. 4's squares lie within this of Table I's photon energies (Ry) and sums
(relative): hand-placed markers 17-18 pixels across on a scan (one pixel 0.00093 Ry
and 1.8 %); they are drawn 0.001-0.0055 Ry to the right of Table I's energies."""
REDRAWN = 0.12
"""Fig. 4's LDRPA agrees with Fig. 2's within this (relative), where Fig. 2's exceeds
0.3: a value of 0.3 is 55 pixels above Fig. 2's bottom line, so its pixel (0.0054)
and line width (3-4 pixels) are under 4 %, and Fig. 4's 1.8 % per pixel."""


REPEATED = 0.12
"""Fig. 3's HF and RPAE lines agree with lhuillier_1987's Fig. 6 readings of the same
calculations within this (relative): 1.9 % (here) and 1.3 % (there) per pixel, and
each reading's line width."""


def _gangopadhyay() -> dict[str, tuple[float, float]]:
    """Return gangopadhyay_1986's Table I: row -> (photon energy in Ry, sum)."""
    out = {}
    columns = []
    for key in ("tab1_jc32", "tab1_jc12"):
        path = DATA / "values" / "ionisation" / "Xe" / f"gangopadhyay_1986__{key}.csv"
        with path.open(newline="", encoding="utf-8") as handle:
            columns.append(list(csv.DictReader(handle)))
    for a, b in zip(*columns, strict=True):
        energy = float(a["x"]) / RYDBERG_EV
        out[a["row"]] = (energy, (float(a["value"]) + float(b["value"])) / 1e-49)
    return out


def _lhuillier_1987() -> dict[str, dict[float, float]]:
    """Return lhuillier_1987's three lines (Figs. 6, 7; 1e-49 cm^4 s) by grid point."""
    out: dict[str, dict[float, float]] = {"hf": {}, "rpae": {}, "ldrpa": {}}
    for name, values in out.items():
        for fig in ("6", "7"):
            path = (
                DATA
                / "values"
                / "ionisation"
                / "Xe"
                / f"lhuillier_1987__fig{fig}_{name}_Xe.csv"
            )
            with path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    values[round(float(row["row"]), 3)] = float(row["value"])
    return out


def _ratio_summary(ratios: list[float]) -> str:
    if not ratios:
        return "none"
    a = np.array(ratios)
    return (
        f"{len(a)} points, median {np.median(a):.3f}, range {a.min():.3f}-{a.max():.3f}"
    )


def checks(readings: dict[str, Reading]) -> dict[str, object]:
    """Check Fig. 4's squares, cross and LDRPA; return the checks and the reports."""
    out: dict[str, object] = {}
    failures = []
    four = readings["4"]
    table = _gangopadhyay()
    squares = sorted(_core_centre(four, box) for box in _cores(four.dark, four.panel))
    rows = list(table)[: len(squares)]
    for (energy, value), row in zip(squares, rows, strict=True):
        expected, total = table[row]
        ratio = value / total
        out[f"Fig. 4 square at {energy:.4f} Ry over Table I {row} eV"] = round(ratio, 3)
        if abs(energy - expected) > SQUARE_ENERGY or abs(ratio - 1) > SQUARE_VALUE:
            failures.append(f"Fig. 4 square at {energy:.4f} Ry: {ratio:.3f}")
    for fig in ("3", "4"):
        cross = readings[fig].panel.cross
        if cross is None:
            continue
        energy, value = _centre(readings[fig], cross)
        out[f"Fig. {fig} cross (Ry; over 193 nm's {MCCOWN_RY:.4f}; over McCown)"] = (
            round(energy, 4),
            round(energy - MCCOWN_RY, 4),
            round(value / MCCOWN_SIGMA, 3),
        )
        if fig == "4" and abs(energy - MCCOWN_RY) > CROSS_ENERGY:
            failures.append(f"Fig. 4's cross at {energy:.4f} Ry")
    for fig in ("2", "4"):
        for k, where in enumerate(readings[fig].panel.hatches, start=1):
            out[f"Fig. {fig} hatched mark {k} at (Ry)"] = round(
                _hatch(readings[fig], where), 4
            )
    linear = {s.grid: s.sigma for s in readings["2"].samples["dashed"]}
    redrawn = [
        s.sigma / linear[s.grid]
        for s in four.samples["dashed"]
        if s.grid in linear and linear[s.grid] > 0.3
    ]
    out["Fig. 4 LDRPA over Fig. 2 LDRPA (where > 0.3)"] = _ratio_summary(redrawn)
    if not redrawn or max(abs(r - 1) for r in redrawn) > REDRAWN:
        failures.append(f"Fig. 4's LDRPA departs from Fig. 2's: {redrawn}")
    pra = _lhuillier_1987()
    for style, name in (("solid", "hf"), ("dashed", "rpae")):
        ratios = [
            s.sigma / pra[name][round(s.grid, 3)]
            for s in readings["3"].samples[style]
            if round(s.grid, 3) in pra[name]
        ]
        out[f"Fig. 3 {name} over lhuillier_1987's Fig. 6"] = _ratio_summary(ratios)
        if not ratios or max(abs(r - 1) for r in ratios) > REPEATED:
            failures.append(f"Fig. 3's {name} departs from lhuillier_1987's")
    for style, name in (("solid", "rpae"), ("dashed", "ldrpa")):
        ratios = [
            s.sigma / pra[name][round(s.grid, 3)]
            for s in four.samples[style]
            if round(s.grid, 3) in pra[name]
        ]
        out[f"Fig. 4 {name} over lhuillier_1987"] = _ratio_summary(ratios)
    if failures:
        raise ValueError("; ".join(failures))
    return out


def readings() -> dict[str, Reading]:
    """Return every panel's reading."""
    return {panel.figure: read(panel) for panel in PANELS}


def _points(reading: Reading) -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for style, file in reading.panel.files.items():
        files[file] = [
            Point(f"{s.grid:.3f}", f"{s.energy:.4f}", f"{s.sigma:.3e}", "", s.how)
            for s in reading.samples[style]
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the samples of Figs. 2, 3 and 4 against the photon energy (Ry)."""
    found = readings()
    checks(found)
    files: dict[str, list[Point]] = {}
    for reading in found.values():
        files.update(_points(reading))
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks(readings()).items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
