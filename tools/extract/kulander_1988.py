# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kulander (1988): single-active-electron ionisation rates of xenon.

The library copy is a scan of Phys. Rev. A 38, 778 (1988): every page one 1-bit image
(5060 x 6898 pixels, about 610 dpi) with an OCR text layer; PDF page n is journal page
777 + n. The paper's cross sections are printed in the text (the record's inline
points); its ionisation rates against intensity are drawn only, on p. 781 (PDF page
4), and are read here from the page image at its native pixels (data rule 8, TC-699):

* **Fig. 1** (left, library figure ``fig_b14628317562``): the rates at 293 nm of the
  5p_1 electrons (squares, the four 5p_+-1), the 5p_0 electrons (triangles) and their
  sum, the p shell's total (circles on a heavy line), with an I^3 fit (dashed), from
  1e13 to 2e14 W/cm^2. Taken: the total.
* **Fig. 2** (right, ``fig_bfb3456aaf40``): the p shell's total rate at 293 nm (solid
  squares), 586 nm (solid circles) and 1064 nm (solid triangles), with experimental
  rates as open symbols. Taken: 586 and 1064 nm; the 293 nm squares check Fig. 1.

How the figures are read:

* **the frames and ticks**: both are logarithmic in both axes. Each frame's four
  lines are written below as two (column, row) points each, and every point is
  checked to be the centre of a dark run of the image (``FRAME_TOLERANCE``). The tick
  marks point inward; each is looked for along the line ``SCAN`` pixels inside its
  side and must lie where the script writes it, within ``TICK_TOLERANCE``. The scan is
  slightly rotated and sheared (Fig. 1's left line leans 9.5 pixels over its height,
  its top line 5 pixels over its width), so every decade is a straight line through
  its two ticks (:class:`Line`), and a value is interpolated in log10 between the two
  decade lines that bracket the pixel, at the pixel's own row or column
  (:class:`Panel`). The decades are not evenly printed (Fig. 2: 264 to 288 pixels per
  decade of rate down the axis): the piecewise scale is exact at every tick.
* **the decades**: the y labels print 10^9 to 10^15 (Fig. 2; the frame's top and
  bottom lines are 10^15 and 10^9) and 10^11 to 10^14 (Fig. 1, at ticks; above the
  10^14 tick the scale is extended with its decade below). In Fig. 1 the 10^14 tick
  of the right side and the 10^11 tick of the left lie under markers: those decade
  lines are drawn parallel to their neighbours (10^13 and 10^12). The x labels print
  10^13 and 10^14; the 10^14 tick is a line from the top to the bottom tick (they
  differ by 11 and 12.5 pixels); the 10^13 decade is the left frame line, where the
  first calculated and measured points of both figures are drawn (the labels' glyphs
  are set 30 pixels off it, and 6 to 17 off the 10^14 tick, so they do not place it).
  Checked: the right frame lines read 2.00 and 2.04e14 W/cm^2 (the text: 'between
  1e13 and 2e14'), every calculated marker lies within ``ROUND`` of a round intensity
  (1, 2, 3, 3.5, 4, 5, 6, 7, 7.5, 8, 10, 15 or 20 x 10^13), and the open symbols the
  text places at 2.8e13 and 4.2e13 W/cm^2 read there.
* **the markers**: the ink is opened with a disk (``Panel.opening``: radius 10 in Fig.
  1, whose heavy line is thick, and 6 in Fig. 2), which removes the lines and the frame
  and leaves the filled markers. Each marker read is seeded at a pixel found by
  inspecting the image (``SEEDS``), takes the opened region there, and must have the
  area of its kind (``AREAS``). Its centre is the region's centroid; a circle taller
  than ``CIRCLE_HEIGHT`` pixels (in Fig. 1 two total circles merge with the 5p_0
  triangle below them, at 3.5e13 and 5e13) is centred ``CIRCLE_HALF`` pixels below its
  top row, the half-height of the clean circles (checked), at the mean column of
  its top ``CIRCLE_HEIGHT`` rows. Markers drawn on top of one another (Fig. 2: 586 and
  1064 nm from 1e14 W/cm^2, 293 nm at 1.5e14) are not read. Each ``printed`` cell gives
  the marker's pixel and the round intensity it lies at.
* **the open symbols** (Fig. 2, not data): their centre is the centroid of the white
  region inside the outline.

Checked (:func:`checks`; the script raises if one is missed): Fig. 2's open square,
'Perry et al. ... measured an ionization rate of 4 x 10^11 s^-1 for a peak intensity
of 1 x 10^13 W/cm^2' (p. 781), within 10 %; Fig. 1's total against the sum of its 5p_0
and 5p_1 markers (within 10 %), and Fig. 2's 293 nm squares against Fig. 1's total at
the same intensities (within 10 %); the dashed I^3 line's slope (3 within 3 %) and its
generalised cross section, rate / (I / hbar omega)^3, against the printed '1 x
10^-82 cm^6 s^2' (one significant figure: 0.5 to 1.5e-82).

One pixel is 0.24 % (Fig. 1) and 0.37 % (Fig. 2) in rate along y, and 0.19 % in
intensity along x.
"""

import math
from dataclasses import dataclass
from functools import cache
from typing import Literal

import numpy as np
from _figure import image
from _library import Point, write
from scipy import ndimage

DIRECTORY = "Kulander_1988_TimeDependentTheory"
PAGE = 4
"""PDF page 4 = p. 781: Figs. 1 and 2."""
SHAPE = (6898, 5060)
DARK = 128
SCAN = (8, 25)
"""The ticks are looked for along the lines this many pixels inside each side."""
TICK_TOLERANCE = 2.0
FRAME_TOLERANCE = 1.0
ROUND = 0.03
"""The largest relative distance of a calculated marker from a round intensity."""
NOMINAL = (1.0, 2.0, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 7.5, 8.0, 10.0, 15.0, 20.0)
"""Round intensities, 1e13 W/cm^2."""
CIRCLE_HEIGHT = 44
CIRCLE_HALF = 20.5
"""Fig. 1's six clean total circles are 41-43 rows tall, their centroids 20.0-21.5
(mean 20.5) rows below their top (checked within ``CIRCLE_SLACK``)."""
CIRCLE_SLACK = 1.5
HC_EV_NM = 1239.841984
"""hc, eV nm (CODATA 2018)."""
EV_J = 1.602176634e-19
"""The electronvolt, J (exact, SI 2019)."""

Kind = Literal["circle", "square", "triangle"]


@dataclass(frozen=True)
class Line:
    """A straight line through two (column, row) points of the page image."""

    a: tuple[float, float]
    b: tuple[float, float]

    def row(self, column: float) -> float:
        """Return the line's row at a column."""
        (c0, r0), (c1, r1) = self.a, self.b
        return r0 + (r1 - r0) * (column - c0) / (c1 - c0)

    def column(self, row: float) -> float:
        """Return the line's column at a row."""
        (c0, r0), (c1, r1) = self.a, self.b
        return c0 + (c1 - c0) * (row - r0) / (r1 - r0)


def _between(position: float, lines: list[tuple[float, float]]) -> float:
    """Return the value at a position from (position, value) decades, piecewise.

    The decades are in order of increasing value; outside them the nearest interval
    is extended.
    """
    for k in range(len(lines) - 1):
        (p0, v0), (p1, v1) = lines[k], lines[k + 1]
        if k == len(lines) - 2 or (position - p1) * (p1 - p0) <= 0:
            return v0 + (v1 - v0) * (position - p0) / (p1 - p0)
    msg = "fewer than two decades"
    raise ValueError(msg)


@dataclass(frozen=True)
class Panel:
    """A figure's frame, its decade lines and how its markers are opened.

    Attributes
    ----------
    name : str
        'Fig. 1' or 'Fig. 2'.
    left, right, top, bottom : Line
        The frame's lines.
    x_decades : tuple of (float, Line)
        log10 of the intensity (W/cm^2) and its line, left to right.
    y_decades : tuple of (float, Line)
        log10 of the rate (s^-1) and its line, bottom to top.
    y_ticks : tuple of (str, tuple of (float, float))
        The ticks written per side ('left', 'right'): (row, log10 rate).
    x_ticks : tuple of (str, tuple of (float, float))
        The ticks written per side ('top', 'bottom'): (column, log10 intensity).
    opening : int
        The radius (pixels) of the disk the ink is opened with.
    """

    name: str
    left: Line
    right: Line
    top: Line
    bottom: Line
    x_decades: tuple[tuple[float, Line], ...]
    y_decades: tuple[tuple[float, Line], ...]
    y_ticks: tuple[tuple[str, tuple[tuple[float, float], ...]], ...]
    x_ticks: tuple[tuple[str, tuple[tuple[float, float], ...]], ...]
    opening: int

    def intensity(self, column: float, row: float) -> float:
        """Return the intensity (W/cm^2) at a pixel."""
        lines = [(line.column(row), value) for value, line in self.x_decades]
        return 10 ** _between(column, lines)

    def rate(self, column: float, row: float) -> float:
        """Return the rate (s^-1) at a pixel."""
        lines = [(line.row(column), value) for value, line in self.y_decades]
        return 10 ** _between(row, lines)


def _y_line(left: Line, right: Line, row_left: float, row_right: float) -> Line:
    """Return a decade line through its left and right ticks (on the frame lines)."""
    return Line((left.column(row_left), row_left), (right.column(row_right), row_right))


def _x_line(top: Line, bottom: Line, column_top: float, column_bottom: float) -> Line:
    """Return a decade line through its top and bottom ticks (on the frame lines)."""
    return Line(
        (column_top, top.row(column_top)), (column_bottom, bottom.row(column_bottom))
    )


# Fig. 1: the frame lines (dark-run centres, read from the image) and the ticks.
F1_LEFT = Line((750.0, 701.0), (759.5, 2240.0))
F1_RIGHT = Line((2352.0, 701.0), (2360.5, 2240.0))
F1_TOP = Line((814.0, 644.0), (2296.0, 639.0))
F1_BOTTOM = Line((814.0, 2301.5), (2296.0, 2298.0))
F1_Y_LEFT = ((944.0, 14.0), (1362.0, 13.0), (1774.5, 12.0))
F1_Y_RIGHT = ((1354.0, 13.0), (1767.5, 12.0), (2184.0, 11.0))
F1_X = (("top", ((1980.5, 14.0),)), ("bottom", ((1991.5, 14.0),)))
_F1_13 = _y_line(F1_LEFT, F1_RIGHT, 1362.0, 1354.0)
_F1_12 = _y_line(F1_LEFT, F1_RIGHT, 1774.5, 1767.5)
# The 10^14 tick of the right side and the 10^11 tick of the left lie under markers:
# those decade lines are parallel to their neighbours.
_F1_14 = _y_line(F1_LEFT, F1_RIGHT, 944.0, 944.0 + (1354.0 - 1362.0))
_F1_11 = _y_line(F1_LEFT, F1_RIGHT, 2184.0 - (1767.5 - 1774.5), 2184.0)
FIG1 = Panel(
    "Fig. 1",
    F1_LEFT,
    F1_RIGHT,
    F1_TOP,
    F1_BOTTOM,
    x_decades=((13.0, F1_LEFT), (14.0, _x_line(F1_TOP, F1_BOTTOM, 1980.5, 1991.5))),
    y_decades=((11.0, _F1_11), (12.0, _F1_12), (13.0, _F1_13), (14.0, _F1_14)),
    y_ticks=(("left", F1_Y_LEFT), ("right", F1_Y_RIGHT)),
    x_ticks=F1_X,
    opening=10,
)

# Fig. 2: the frame's top and bottom lines are the 10^15 and 10^9 decades.
F2_LEFT = Line((2965.5, 709.0), (2968.0, 2238.0))
F2_RIGHT = Line((4575.5, 1473.0), (4576.5, 2238.0))
F2_TOP = Line((3027.0, 647.0), (4515.0, 651.0))
F2_BOTTOM = Line((3027.0, 2297.0), (4515.0, 2300.0))
F2_Y_LEFT = (
    (913.0, 14.0),
    (1186.0, 13.0),
    (1453.0, 12.0),
    (1728.5, 11.0),
    (2011.0, 10.0),
)
F2_Y_RIGHT = (
    (917.0, 14.0),
    (1190.0, 13.0),
    (1456.0, 12.0),
    (1732.0, 11.0),
    (2013.5, 10.0),
)
F2_X = (("top", ((4201.0, 14.0),)), ("bottom", ((4188.5, 14.0),)))
FIG2 = Panel(
    "Fig. 2",
    F2_LEFT,
    F2_RIGHT,
    F2_TOP,
    F2_BOTTOM,
    x_decades=((13.0, F2_LEFT), (14.0, _x_line(F2_TOP, F2_BOTTOM, 4201.0, 4188.5))),
    y_decades=(
        (9.0, F2_BOTTOM),
        *(
            (left[1], _y_line(F2_LEFT, F2_RIGHT, left[0], right[0]))
            for left, right in reversed(tuple(zip(F2_Y_LEFT, F2_Y_RIGHT, strict=True)))
        ),
        (15.0, F2_TOP),
    ),
    y_ticks=(("left", F2_Y_LEFT), ("right", F2_Y_RIGHT)),
    x_ticks=F2_X,
    opening=6,
)

AREAS: dict[tuple[str, Kind], tuple[int, int]] = {
    ("Fig. 1", "circle"): (1150, 1700),
    ("Fig. 1", "triangle"): (300, 500),
    ("Fig. 1", "square"): (740, 920),
    ("Fig. 2", "circle"): (1150, 1280),
    ("Fig. 2", "triangle"): (450, 600),
    ("Fig. 2", "square"): (840, 980),
}
"""The opened area (pixels) of each kind of marker, per figure."""

SEEDS: dict[str, tuple[Panel, Kind, tuple[tuple[int, int], ...]]] = {
    "fig1_total": (
        FIG1,
        "circle",
        (
            (759, 1937),
            (1131, 1621),
            (1425, 1316),
            (1620, 1101),
            (1836, 952),
            (1990, 874),
            (2206, 798),
            (2351, 783),
        ),
    ),
    "fig1_5p0": (
        FIG1,
        "triangle",
        (
            (759, 1982),
            (1135, 1660),
            (1431, 1353),
            (1841, 994),
            (1994, 941),
            (2206, 865),
            (2351, 939),
        ),
    ),
    "fig1_5p1": (
        FIG1,
        "square",
        (
            (758, 2205),
            (1129, 1880),
            (1430, 1622),
            (1619, 1445),
            (1835, 1251),
            (1992, 1088),
            (2208, 1017),
            (2352, 889),
        ),
    ),
    "fig2_293": (
        FIG2,
        "square",
        (
            (3336, 1355),
            (3636, 1166),
            (3829, 1014),
            (4046, 929),
            (4197, 876),
            (4571, 821),
        ),
    ),
    "fig2_586": (
        FIG2,
        "circle",
        (
            (2969, 2068),
            (3332, 1691),
            (3551, 1465),
            (3632, 1337),
            (3707, 1275),
            (3826, 1181),
            (4042, 1015),
        ),
    ),
    "fig2_1064": (
        FIG2,
        "triangle",
        (
            (3329, 2062),
            (3551, 1682),
            (3704, 1471),
            (3826, 1346),
            (3925, 1238),
            (4013, 1160),
            (4085, 1092),
        ),
    ),
}
"""Each series read: its figure, marker kind and seeds (column, row) by inspection."""

OPEN_SYMBOLS = {
    "Perry et al., 293 nm (open square)": (2960, 1563),
    "Perry et al., 586 nm (open circle)": (3731, 1332),
    "Freeman et al., 1064 nm (open triangle)": (3511, 1878),
}
"""Fig. 2's experimental points read for the checks: a pixel inside each outline."""
PERRY_RATE = 4e11
"""'an ionization rate of 4 x 10^11 s^-1 for a peak intensity of 1 x 10^13' (p. 781)."""
PRINTED_INTENSITIES = {
    "Perry et al., 586 nm (open circle)": 4.2e13,
    "Freeman et al., 1064 nm (open triangle)": 2.8e13,
}
"""The intensities the text gives the open symbols (W/cm^2, pp. 781-782)."""
DASHES = (1460, 1830)
"""The columns over which Fig. 1's dashed I^3 line is the topmost ink."""
SIGMA3_PRINTED = 1e-82
"""'The generalized cross section corresponding to this line is 1 x 10^-82 cm6 s2'."""
WAVELENGTH_293 = 293.0
FILES = {
    "fig1_total": "Xe/kulander_1988__fig1_rate_293_Xe.csv",
    "fig2_586": "Xe/kulander_1988__fig2_rate_586_Xe.csv",
    "fig2_1064": "Xe/kulander_1988__fig2_rate_1064_Xe.csv",
}


def _runs(line: np.ndarray, offset: float) -> list[float]:
    """Return the centres of the runs of dark pixels along a line of the image."""
    index = np.flatnonzero(line)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [(g[0] + g[-1]) / 2 + offset for g in groups]


@cache
def _dark() -> np.ndarray:
    """Return the page's dark pixels."""
    grey = image(DIRECTORY, PAGE)
    if grey.shape != SHAPE:
        msg = f"p. 781: an image of {grey.shape}, not {SHAPE}"
        raise ValueError(msg)
    return grey < DARK


def _check_lines(panel: Panel) -> None:
    """Raise unless the frame lines are dark runs where the script writes them."""
    dark = _dark()
    for side in ("left", "right", "top", "bottom"):
        line: Line = getattr(panel, side)
        for column, row in (line.a, line.b):
            if side in {"left", "right"}:
                start = int(column) - 12
                found = _runs(dark[int(row), start : start + 25], start)
                at = column
            else:
                start = int(row) - 12
                found = _runs(dark[start : start + 25, int(column)], start)
                at = row
            if not any(abs(f - at) <= FRAME_TOLERANCE for f in found):
                msg = f"{panel.name}: no {side} line at ({column}, {row}): {found}"
                raise ValueError(msg)


def _check_ticks(panel: Panel) -> None:
    """Raise unless every tick written is a dark run ``SCAN`` pixels inside its side."""
    dark = _dark()
    for side, ticks in (*panel.y_ticks, *panel.x_ticks):
        line: Line = getattr(panel, side)
        sign = 1 if side in {"left", "top"} else -1
        for position, value in ticks:
            for inset in SCAN:
                start = int(position) - 15
                if side in {"left", "right"}:
                    column = round(line.column(position) + sign * inset)
                    found = _runs(dark[start : start + 31, column], start)
                else:
                    row = round(line.row(position) + sign * inset)
                    found = _runs(dark[row, start : start + 31], start)
                if not any(abs(f - position) <= TICK_TOLERANCE for f in found):
                    msg = f"{panel.name}: no {side} tick 10^{value:g} at {position}"
                    raise ValueError(msg)


def _check_frame(panel: Panel) -> None:
    """Raise unless the frame lines and the ticks are where the script writes them."""
    _check_lines(panel)
    _check_ticks(panel)


def _disk(radius: int) -> np.ndarray:
    """Return a disk structuring element."""
    y, x = np.mgrid[-radius : radius + 1, -radius : radius + 1]
    return x * x + y * y <= radius * radius


@cache
def _opened(panel: Panel) -> np.ndarray:
    """Return the labelled regions of the panel's ink opened with its disk."""
    dark = _dark()
    rows = slice(int(panel.top.row(panel.left.a[0])) - 60, int(panel.bottom.b[1]) + 60)
    columns = slice(int(panel.left.a[0]) - 60, int(panel.right.b[0]) + 60)
    opened = np.zeros_like(dark)
    opened[rows, columns] = ndimage.binary_opening(
        dark[rows, columns], structure=_disk(panel.opening)
    )
    labels, _ = ndimage.label(opened)
    return labels


def _region(
    panel: Panel, kind: Kind, seed: tuple[int, int]
) -> tuple[np.ndarray, np.ndarray]:
    """Return the rows and columns of the opened region at a marker's seed."""
    labels = _opened(panel)
    column, row = seed
    window = labels[row - 12 : row + 13, column - 12 : column + 13]
    found = set(window[window > 0].tolist())
    if len(found) != 1:
        msg = f"{panel.name}: {len(found)} regions at the {kind} seed {seed}"
        raise ValueError(msg)
    rows, columns = np.nonzero(labels == found.pop())
    low, high = AREAS[(panel.name, kind)]
    if not low <= rows.size <= high:
        msg = f"{panel.name}: the {kind} at {seed} has an area of {rows.size}"
        raise ValueError(msg)
    return rows, columns


def _marker(panel: Panel, kind: Kind, seed: tuple[int, int]) -> tuple[float, float]:
    """Return a marker's centre (column, row) from the opened region at its seed."""
    rows, columns = _region(panel, kind, seed)
    if kind == "circle" and rows.max() - rows.min() + 1 > CIRCLE_HEIGHT:
        top = rows.min()
        upper = rows < top + CIRCLE_HEIGHT
        return float(columns[upper].mean()), top + CIRCLE_HALF
    return float(columns.mean()), float(rows.mean())


def _check_circle_half(out: dict[str, str]) -> None:
    """Check that Fig. 1's clean circles are centred ``CIRCLE_HALF`` below their top."""
    panel, kind, seeds = SEEDS["fig1_total"]
    offsets = []
    for seed in seeds:
        rows, _ = _region(panel, kind, seed)
        if rows.max() - rows.min() + 1 <= CIRCLE_HEIGHT:
            offsets.append(float(rows.mean()) - rows.min())
    out["Fig. 1 clean circles: centre below the top (pixels)"] = ", ".join(
        f"{offset:.1f}" for offset in offsets
    )
    if max(abs(offset - CIRCLE_HALF) for offset in offsets) > CIRCLE_SLACK:
        msg = f"Fig. 1: the clean circles are centred {offsets} below their tops"
        raise ValueError(msg)


@dataclass(frozen=True)
class Read:
    """A marker read: its pixel, intensity (W/cm^2), rate (s^-1), and nominal."""

    column: float
    row: float
    intensity: float
    rate: float
    nominal: float


def _nominal(intensity: float) -> float:
    """Return the round intensity (W/cm^2) a marker lies at, or raise."""
    best = min(NOMINAL, key=lambda n: abs(math.log(intensity / (n * 1e13))))
    if abs(intensity / (best * 1e13) - 1) > ROUND:
        msg = f"a marker at {intensity:.3e} W/cm^2 is off every round intensity"
        raise ValueError(msg)
    return best * 1e13


@cache
def series(name: str) -> tuple[Read, ...]:
    """Return a series' markers, left to right."""
    panel, kind, seeds = SEEDS[name]
    _check_frame(panel)
    out = []
    for seed in seeds:
        column, row = _marker(panel, kind, seed)
        intensity = panel.intensity(column, row)
        out.append(
            Read(column, row, intensity, panel.rate(column, row), _nominal(intensity))
        )
    return tuple(out)


def _hole(seed: tuple[int, int]) -> tuple[float, float]:
    """Return the centroid of the white region enclosed around an open symbol."""
    dark = _dark()
    column, row = seed
    window = ~dark[row - 40 : row + 41, column - 40 : column + 41]
    labels, _ = ndimage.label(window)
    index = labels[40, 40]
    rows, columns = np.nonzero(labels == index)
    if index == 0 or rows.min() == 0 or columns.min() == 0 or rows.max() == 80:
        msg = f"no enclosed white region at {seed}"
        raise ValueError(msg)
    return float(columns.mean()) + column - 40, float(rows.mean()) + row - 40


def _dashed_line() -> tuple[float, float, int]:
    """Return Fig. 1's dashed line: its slope, log10 c of W = c I^3, and samples.

    In each column of ``DASHES`` the topmost dark run inside the frame is sampled;
    between the dashes it is a solid line below (a smaller rate), and those samples
    are dropped by refitting until every sample kept lies no more than 0.02 decade
    below the line. The slope is fitted; c is the mean of log10 W - 3 log10 I over
    the samples kept (W in s^-1, I in W/cm^2).
    """
    dark = _dark()
    samples = []
    for column in range(*DASHES):
        top = round(FIG1.top.row(column)) + 20
        found = _runs(dark[top : int(FIG1.bottom.row(column)) - 20, column], top)
        if found:
            row = found[0]
            samples.append(
                (
                    math.log10(FIG1.intensity(column, row)),
                    math.log10(FIG1.rate(column, row)),
                )
            )
    points = np.array(samples)
    keep = np.ones(len(points), dtype=bool)
    for _ in range(10):
        slope, intercept = np.polyfit(points[keep, 0], points[keep, 1], 1)
        residual = points[:, 1] - (slope * points[:, 0] + intercept)
        keep = residual > -0.02
    slope, _ = np.polyfit(points[keep, 0], points[keep, 1], 1)
    log_c = float(np.mean(points[keep, 1] - 3 * points[keep, 0]))
    return float(slope), log_c, int(keep.sum())


def _check_calibration(out: dict[str, str]) -> None:
    """Check the right frame lines and the open symbols against the text."""
    for panel in (FIG1, FIG2):
        _check_frame(panel)
        column = panel.right.b[0]
        middle = (panel.top.row(column) + panel.bottom.row(column)) / 2
        edge = panel.intensity(panel.right.column(middle), middle)
        out[f"{panel.name}: the right frame line (W/cm2)"] = f"{edge:.3e}"
        if abs(edge / 2e14 - 1) > 0.05:
            msg = f"{panel.name}: the right frame reads {edge:.3e}, not 2e14 W/cm2"
            raise ValueError(msg)
    for label, seed in OPEN_SYMBOLS.items():
        column, row = _hole(seed)
        intensity, rate = FIG2.intensity(column, row), FIG2.rate(column, row)
        out[f"Fig. 2, {label}"] = f"{intensity:.3e} W/cm2, {rate:.3e} s-1"
        if label in PRINTED_INTENSITIES:
            printed = PRINTED_INTENSITIES[label]
            if abs(intensity / printed - 1) > ROUND:
                msg = f"Fig. 2: {label} at {intensity:.3e}, the text {printed:.1e}"
                raise ValueError(msg)
        elif abs(rate / PERRY_RATE - 1) > 0.1:
            msg = f"Fig. 2: {label} reads {rate:.3e} s-1, the text 4e11"
            raise ValueError(msg)


def _check_series(out: dict[str, str]) -> None:
    """Check Fig. 1's total against its parts, and Fig. 2's 293 nm against Fig. 1."""
    total = {r.nominal: r for r in series("fig1_total")}
    p0 = {r.nominal: r for r in series("fig1_5p0")}
    p1 = {r.nominal: r for r in series("fig1_5p1")}
    for nominal in sorted(set(total) & set(p0) & set(p1)):
        ratio = total[nominal].rate / (p0[nominal].rate + p1[nominal].rate)
        out[f"Fig. 1 at {nominal:.2g}: total / (5p0 + 5p1)"] = f"{ratio:.3f}"
        if abs(ratio - 1) > 0.1:
            msg = f"Fig. 1 at {nominal:.2g}: the total is {ratio:.3f} x the sum"
            raise ValueError(msg)
    for square in series("fig2_293"):
        if square.nominal in total:
            ratio = square.rate / total[square.nominal].rate
            out[f"293 nm at {square.nominal:.2g}: Fig. 2 / Fig. 1"] = f"{ratio:.3f}"
            if abs(ratio - 1) > 0.1:
                msg = f"293 nm at {square.nominal:.2g}: Fig. 2 is {ratio:.3f} x Fig. 1"
                raise ValueError(msg)


def _check_dashed(out: dict[str, str]) -> None:
    """Check Fig. 1's dashed line: an I^3 law with the printed sigma_3."""
    slope, log_c, count = _dashed_line()
    photon = HC_EV_NM / WAVELENGTH_293 * EV_J
    # sigma_3 = W / Phi^3 with Phi = I / (hbar omega), I in W/cm^2: W = c I^3.
    sigma3 = 10**log_c * photon**3
    out["Fig. 1 dashed line: slope, samples"] = f"{slope:.3f}, {count}"
    out["Fig. 1 dashed line: sigma_3 (cm6 s2)"] = f"{sigma3:.3e}"
    if abs(slope / 3 - 1) > 0.03:
        msg = f"Fig. 1: the dashed line's slope is {slope:.3f}, not 3"
        raise ValueError(msg)
    if not 0.5 * SIGMA3_PRINTED <= sigma3 <= 1.5 * SIGMA3_PRINTED:
        msg = f"Fig. 1: the dashed line gives sigma_3 = {sigma3:.3e}, printed 1e-82"
        raise ValueError(msg)


def checks() -> dict[str, str]:
    """Return the checks of the reading against the text and between the figures."""
    out: dict[str, str] = {}
    _check_calibration(out)
    _check_circle_half(out)
    _check_series(out)
    _check_dashed(out)
    return out


def extract() -> dict[str, list[Point]]:
    """Return Fig. 1's 293 nm total and Fig. 2's 586 and 1064 nm rates."""
    files: dict[str, list[Point]] = {}
    for name, file in FILES.items():
        files[file] = [
            Point(
                f"{r.nominal:.3g}",
                f"{r.intensity:.4g}",
                f"{r.rate:.3e}",
                "",
                f"pixel ({r.column:.1f}, {r.row:.1f}); at {r.nominal:.3g} W/cm2",
            )
            for r in series(name)
        ]
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
