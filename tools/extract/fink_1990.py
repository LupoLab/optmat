# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Fink and Johnson (1990): relativistic TDDF two-photon cross sections.

The paper prints no table of cross sections: its two-photon ionisation cross sections
are given only in figures (data rule 8, TC-699), and the PDF is a scan (every page one
1-bit image, 5060 x 6898 pixels, about 600 dpi, with an OCR text layer that is not
used), so the figures are read here from the native pixels (``_figure.image``). Table I
(the Dirac-Fock binding energies, which are the calculation's thresholds, and the
experimental ones) is read from the library's table cells for the checks; it is not
written as data.

* **Fig. 3** (p. 13, library figure ``fig_0567865d2622``, its second panel crop): Ar,
  total sigma^(2) in 1e-50 cm^4 s on linear axes against the photon energy omega in
  atomic units, 0.289 to 0.56.
* **Fig. 4** (p. 13, ``fig_607b97ca1f36``): Kr, the same, 0.27 to 0.46.
* **Fig. 5** (p. 13, ``fig_3ff16879fe95``): Xe, the same on a logarithmic y axis (1e-1
  to 1e4), 0.25 to 0.35.
* **Fig. 6** (p. 14, ``fig_7bd44496e19f``): Xe, the partial cross sections to the
  2P3/2 and 2P1/2 ionic cores, same axes as Fig. 5.

In Figs. 3 to 6 each computed value is a marker (Figs. 3-5: + linear, diamond circular
polarisation; Fig. 6: + linear 2P3/2, x linear 2P1/2, diamond circular 2P3/2, square
circular 2P1/2), joined by straight lines (solid for +, dashed for diamonds). The
**markers** are read, not the lines.

How each panel is read:

* **the frame**: its four lines are fitted to the image (the centre of the dark run
  near the expected position in every other row or column, outliers beyond 1.5 pixels
  rejected); a pixel (column, row) is mapped to *corrected* coordinates through the
  mean slopes of the frame's lines (a shear of the scan of at most a few thousandths of
  a pixel per pixel here).
* **the ticks**: every tick has a stub inside the frame (the labelled ones also cross
  it), found as a run of rows or columns dark from ``TICK_INSET`` pixels inside the
  frame line's centre; each must be where ``x_ticks`` and ``y_ticks`` write it, within
  ``TICK_TOLERANCE`` pixels, and dark stretches there that are not ticks (a curve or a
  marker meeting the frame) are written in ``not_x_ticks`` and ``not_y_ticks``. The
  axes are calibrated from those positions (``_figure.Axis``) and every tick is checked
  on the fitted scale.
* **the markers**: each marker read is seeded at a pixel found by inspecting the image
  (``seeds``; nie_2022 reads its markers the same way), and placed by the script: with
  the frame lines erased, the shape's template (:func:`_regions`: the + and x arms, the
  diamond's and the square's sides, about 24 pixels across, strokes 4 to 6 pixels) is
  scored at every pixel within ``REACH`` of the seed (the smallest dark fraction of its
  regions; a + or an x must also end, its arms' ends light in one direction at least,
  so that a steep line is not taken for one), and the marker's centre is the mean of
  the pixels within 2 of the best whose score is within ``PLATEAU`` of it. The best
  score must reach the shape's ``THRESHOLD``, and no two markers of one curve may land
  within ``CLASH`` pixels. A marker hidden in a cluster of others (where the two
  polarisations' curves cross, and in Fig. 6's dense stretches), or on a steep stretch
  where it cannot be told from the line, is not seeded: the record lists the gaps.
  The first marker of a curve may sit on the left frame line, where the plot is
  clipped and only its right half is drawn (``edge``): its row is the middle of the
  dark rows just inside the line (:func:`_edge_row`; no tick of the axis within
  ``EDGE_CLEAR`` rows), its omega the frame line's (the marker's visible half is a
  full half, so its centre is on the line within about 2 pixels).

Every value's ``printed`` cell gives its marker's centre (column, row) in the page
image. One pixel is 2.1e-4 a.u. (0.006 eV) in omega in Fig. 3, 1.5e-4 a.u. (0.004 eV)
in Fig. 4 and 0.8e-4 a.u. (0.002 eV) in Figs. 5 and 6; in sigma^(2) it is 0.0035e-50
cm^4 s (Fig. 3), 0.0078e-50 (Fig. 4) or 0.9 % (Figs. 5 and 6). A marker's centre is
found within about 2 pixels (a + and a diamond drawn over each other, within 4).

Checks (:func:`checks`; the script fails if one is missed): the frame's and the ticks'
positions (above); Table I's fine-structure splittings are the differences of its
binding energies within the rounding (0.01 eV); each curve's first computed frequency
lies at or above half the *experimental* 2P1/2 binding energy of Table I (within
``SAME``) and its last below the 2P3/2 one-photon threshold; Fig. 6's partial cross
sections add up to Fig. 5's totals at the frequencies where all three are read
(within ``SUM``); Fig. 5's totals against the Fink-Johnson curves that Lagutin et al.
(2017) redraw in their Fig. 8 (lagutin_2017, TC-880; vector, read here with that
script's machinery), 6.80-8.10 and 8.70-9.28 eV, away from the 6s resonance (within
``REDRAWN``); Fig. 3 against Starace and Jiang's ESI totals (starace_1987, TC-878),
which the text says agree 'within drawing accuracy' (p. 3812), 8.1-9.4 eV (within
``ESI``). Reported only: the first frequencies against half the Dirac-Fock 2P1/2
energy; how many markers of two curves stand at one computed frequency (within
``SAME``); Fig. 3 against McKenna and van der Hart (2004, TC-856); Fig. 4 against
McGuire (1981, TC-879).
"""

import csv
import functools
import itertools
import math
from dataclasses import dataclass

import lagutin_2017
import numpy as np
from _figure import Axis, image
from _json_table import json_table
from _library import Point, path, write
from scipy import ndimage

from optmat.constants import AU_ENERGY, E_CHARGE, HBAR, C

DIRECTORY = "Fink_1990_ElectronCorrelationsSpin"
SOURCE = "fink_1990"
TABLE_I = "tbl_175a623d1be6"
SHAPE = (6898, 5060)
DARK = 128
TICK_INSET = (3, 8)
"""A tick darkens the rows (columns) this far inside the frame line's centre, pixels."""
TICK_TOLERANCE = 2.0
"""A tick is found within this of where the script writes it, pixels."""
FRAME_CLEAR = 4
"""Pixels inside a frame line's centre erased with it before the markers are read."""
REACH = 6
"""A marker's centre is looked for within this of its seed, pixels."""
PLATEAU = 0.03
"""Pixels near the best whose score is within this of it share the centre."""
CLASH = 7
"""Two markers of one curve closer than this (pixels) are one marker read twice."""
EDGE_CLEAR = 16
"""A marker on the left frame line is read only with no tick within this, rows."""
X_TOLERANCE = 1.0
"""The x ticks lie on their fitted scale within this, pixels (worst 0.8)."""
Y_TOLERANCE = 3.0
"""The y ticks lie on their fitted scale within this, pixels: the printed spacing of the
y ticks varies by a few pixels (Fig. 4's from 125.5 to 129.5 pixels per 1e-50 cm^4 s,
Fig. 5's decades from 252 to 259 pixels); the worst residual is 2.7."""
LIGHT = 0.25
"""A light region (a + or x arm's end) is at most this dark."""
THRESHOLD = {"plus": 0.75, "cross": 0.75, "diamond": 0.7, "square": 0.7}
"""The smallest best score of a seeded marker of each shape."""
SAME = 0.0006
"""Markers of two curves at one computed frequency lie within this, a.u. (4 pixels;
markers of neighbouring frequencies are at least 0.0009 apart where read)."""
SUM = 0.10
"""Fig. 6's partials add up to Fig. 5's total within this, relative."""
REDRAWN = 0.25
"""Fig. 5 against Lagutin et al.'s redrawing of it, relative (a hand redrawing)."""
ESI = 0.15
"""Fig. 3 against Starace and Jiang's ESI totals, relative ('drawing accuracy')."""
HARTREE_EV = AU_ENERGY / E_CHARGE
"""The hartree in eV (CODATA, scipy), for Table I's energies."""
HC_EV_ANGSTROM = 2 * math.pi * HBAR * C / E_CHARGE * 1e10
"""h c in eV angstrom, for McGuire's wavelengths."""


@dataclass(frozen=True)
class Panel:
    """One figure: where it is, its ticks, and its curves' markers.

    Tick positions are corrected pixel coordinates: a column for the x axis, a row for
    the y axis. ``where`` is ``left``, ``right``, ``bottom`` or ``top`` for a frame
    line that carries an end value, or ``tick``.
    """

    figure: str
    page: int
    gas: str
    frame: tuple[int, int, int, int]
    """The frame lines' approximate columns (left, right) and rows (top, bottom)."""
    x_ticks: tuple[tuple[float, float, str], ...]
    """(position, omega in a.u., where)."""
    y_ticks: tuple[tuple[float, float, str], ...]
    """(position, sigma^(2) in 1e-50 cm^4 s, where)."""
    log: bool
    curves: dict[str, str]
    """The data file of each marker shape's curve."""
    seeds: dict[str, tuple[tuple[int, int], ...]]
    """Each shape's seeds, (column, row) in the page image."""
    edge: dict[str, tuple[int, int]]
    """Each shape's marker on the left frame line: (seed row, half window), rows."""
    not_x_ticks: tuple[float, ...] = ()
    not_y_ticks: tuple[float, ...] = ()


def _file(gas: str, key: str) -> str:
    return f"{gas}/{SOURCE}__{key}.csv"


PANELS = {
    "fig3": Panel(
        figure="3",
        page=13,
        gas="Ar",
        frame=(887, 2165, 4649, 5926),
        # Labelled 0.32, 0.40, 0.48 and 0.56 (the right frame line); a long tick every
        # 0.04, a short one between.
        x_ticks=(
            (938.5, 0.30, "tick"),
            (1033.5, 0.32, "tick"),
            (1127.5, 0.34, "tick"),
            (1223.0, 0.36, "tick"),
            (1318.0, 0.38, "tick"),
            (1413.0, 0.40, "tick"),
            (1508.0, 0.42, "tick"),
            (1601.5, 0.44, "tick"),
            (1697.0, 0.46, "tick"),
            (1791.6, 0.48, "tick"),
            (1886.1, 0.50, "tick"),
            (1980.1, 0.52, "tick"),
            (2075.6, 0.54, "tick"),
            (2170.6, 0.56, "right"),
        ),
        # Labelled 0.0 (the bottom line) to 4.0; a tick every 0.5.
        y_ticks=(
            (5930.5, 0.0, "bottom"),
            (5789.5, 0.5, "tick"),
            (5650.5, 1.0, "tick"),
            (5509.5, 1.5, "tick"),
            (5365.5, 2.0, "tick"),
            (5223.5, 2.5, "tick"),
            (5080.4, 3.0, "tick"),
            (4938.4, 3.5, "tick"),
            (4796.9, 4.0, "tick"),
        ),
        log=False,
        curves={
            "plus": _file("Ar", "fig3_linear"),
            "diamond": _file("Ar", "fig3_circular"),
        },
        seeds={
            "plus": (
                (908, 5613),
                (935, 5616),
                (1031, 5629),
                (1126, 5646),
                (1173, 5656),
                (1221, 5668),
                (1314, 5700),
                (1411, 5713),
                (1458, 5719),
                (1486, 5715),
                (1503, 5700),
                (1530, 5635),
                (1537, 5590),
                (1549, 5376),
                (1692, 5286),
                (1713, 5487),
                (1741, 5608),
                (1769, 5674),
                (1789, 5709),
                (1816, 5747),
                (1837, 5772),
                (1857, 5795),
                (1866, 5806),
                (1876, 5813),
                (1886, 5821),
                (1893, 5829),
            ),
            "diamond": (
                (908, 5740),
                (936, 5730),
                (1031, 5705),
                (1126, 5692),
                (1173, 5691),
                (1220, 5693),
                (1314, 5700),
                (1408, 5710),
                (1458, 5702),
                (1484, 5684),
                (1505, 5648),
                (1527, 5530),
                (1537, 5458),
                (1549, 5156),
                (1678, 4790),
                (1687, 4989),
                (1711, 5266),
                (1740, 5458),
                (1769, 5565),
                (1789, 5620),
                (1816, 5677),
                (1836, 5711),
                (1856, 5743),
                (1875, 5769),
                (1893, 5791),
                (1910, 5808),
            ),
        },
        edge={},
    ),
    "fig4": Panel(
        figure="4",
        page=13,
        gas="Kr",
        frame=(3164, 4443, 631, 1913),
        # Labelled 0.28 to 0.44; a tick every 0.02.
        x_ticks=(
            (3232.0, 0.28, "tick"),
            (3366.0, 0.30, "tick"),
            (3501.0, 0.32, "tick"),
            (3635.5, 0.34, "tick"),
            (3771.5, 0.36, "tick"),
            (3905.0, 0.38, "tick"),
            (4040.5, 0.40, "tick"),
            (4175.0, 0.42, "tick"),
            (4309.5, 0.44, "tick"),
        ),
        # Labelled 0.0 (the bottom line) to 10.0 (the top line); a tick every 1.0. The
        # tick at 2.0 merges with the first + marker on the left line (not used).
        y_ticks=(
            (1911.7, 0.0, "bottom"),
            (1783.5, 1.0, "tick"),
            (1526.0, 3.0, "tick"),
            (1397.5, 4.0, "tick"),
            (1268.0, 5.0, "tick"),
            (1142.5, 6.0, "tick"),
            (1015.0, 7.0, "tick"),
            (888.0, 8.0, "tick"),
            (759.5, 9.0, "tick"),
            (630.3, 10.0, "top"),
        ),
        log=False,
        curves={
            "plus": _file("Kr", "fig4_linear"),
            "diamond": _file("Kr", "fig4_circular"),
        },
        seeds={
            "plus": (
                (3231, 1651),
                (3258, 1656),
                (3298, 1662),
                (3339, 1670),
                (3365, 1675),
                (3413, 1686),
                (3433, 1691),
                (3500, 1709),
                (3546, 1720),
                (3566, 1727),
                (3606, 1738),
                (3634, 1746),
                (3713, 1757),
                (3762, 1736),
                (3770, 1719),
                (3807, 1258),
                (3903, 804),
                (3925, 1239),
                (3962, 867),
                (4082, 1323),
                (4107, 1563),
                (4110, 1579),
                (4134, 1673),
                (4158, 1724),
                (4173, 1750),
                (4180, 1759),
                (4202, 1786),
                (4221, 1807),
                (4240, 1823),
                (4259, 1834),
                (4292, 1761),
            ),
            "diamond": (
                (3230, 1710),
                (3258, 1708),
                (3298, 1707),
                (3339, 1709),
                (3365, 1711),
                (3413, 1715),
                (3433, 1718),
                (3482, 1725),
                (3500, 1730),
                (3546, 1738),
                (3566, 1745),
                (3606, 1751),
                (3664, 1758),
                (3713, 1757),
                (3763, 1694),
                (3768, 1669),
                (3807, 1178),
                (3904, 915),
                (3926, 1318),
                (3961, 1337),
                (3971, 1172),
                (4082, 1093),
                (4107, 1393),
                (4109, 1415),
                (4134, 1556),
                (4157, 1636),
                (4173, 1678),
                (4180, 1690),
                (4201, 1732),
                (4222, 1764),
                (4240, 1789),
                (4256, 1806),
                (4276, 1808),
                (4292, 1708),
            ),
        },
        edge={"diamond": (1710, 10)},
        # The tick at 2.0 with the first + marker; the first diamond marker.
        not_y_ticks=(1647.0, 1707.0, 1724.5),
    ),
    "fig5": Panel(
        figure="5",
        page=13,
        gas="Xe",
        frame=(3113, 4394, 4462, 5740),
        # Labelled 0.26 to 0.34; a tick every 0.01.
        x_ticks=tuple(
            (position, round(0.26 + 0.01 * k, 2), "tick")
            for k, position in enumerate(
                (3240.0, 3369.0, 3496.5, 3625.0, 3753.5, 3881.0, 4009.5, 4138.0, 4265.5)
            )
        ),
        # Decades, labelled 10^-1 (the bottom line) to 10^4 (the top line).
        y_ticks=(
            (5736.3, 0.1, "bottom"),
            (5484.5, 1.0, "tick"),
            (5226.0, 10.0, "tick"),
            (4969.5, 100.0, "tick"),
            (4712.5, 1000.0, "tick"),
            (4457.9, 10000.0, "top"),
        ),
        log=True,
        curves={
            "plus": _file("Xe", "fig5_linear"),
            "diamond": _file("Xe", "fig5_circular"),
        },
        seeds={
            "plus": (
                (3189, 5338),
                (3341, 5352),
                (3418, 5364),
                (3495, 5374),
                (3572, 5392),
                (3624, 5395),
                (3787, 5242),
                (3803, 5174),
                (3818, 5098),
                (3833, 4993),
                (3846, 4853),
                (3857, 4650),
                (3933, 4918),
                (3944, 4988),
                (4227, 5144),
                (4240, 5111),
                (4288, 4904),
            ),
            "diamond": (
                (3189, 5352),
                (3265, 5357),
                (3341, 5366),
                (3418, 5378),
                (3495, 5390),
                (3572, 5395),
                (3624, 5395),
                (3689, 5382),
                (3739, 5340),
                (3800, 5200),
                (3818, 5150),
                (3830, 5085),
                (3842, 4985),
                (3857, 4811),
                (3933, 4950),
                (3946, 5005),
                (3953, 5040),
                (3971, 5084),
                (4006, 5153),
                (4080, 5218),
                (4242, 5168),
                (4288, 5022),
            ),
        },
        edge={"plus": (5333, 5), "diamond": (5350, 10)},
        # The first markers on the left line.
        not_y_ticks=(5332.0, 5346.0, 5351.0, 5353.0),
    ),
    "fig6": Panel(
        figure="6",
        page=14,
        gas="Xe",
        frame=(911, 2187, 678, 1958),
        x_ticks=tuple(
            (position, round(0.26 + 0.01 * k, 2), "tick")
            for k, position in enumerate(
                (1041.5, 1168.5, 1297.0, 1425.0, 1552.0, 1679.5, 1807.5, 1935.0, 2062.0)
            )
        ),
        y_ticks=(
            (1957.5, 0.1, "bottom"),
            (1700.5, 1.0, "tick"),
            (1445.0, 10.0, "tick"),
            (1193.0, 100.0, "tick"),
            (935.0, 1000.0, "tick"),
            (677.8, 10000.0, "top"),
        ),
        log=True,
        curves={
            "plus": _file("Xe", "fig6_linear_p32"),
            "cross": _file("Xe", "fig6_linear_p12"),
            "diamond": _file("Xe", "fig6_circular_p32"),
            "square": _file("Xe", "fig6_circular_p12"),
        },
        seeds={
            "plus": (
                (988, 1601),
                (1065, 1604),
                (1142, 1612),
                (1218, 1627),
                (1294, 1642),
                (1372, 1660),
                (1424, 1671),
                (1474, 1671),
                (1728, 1149),
                (1883, 1467),
                (1908, 1475),
                (1921, 1475),
                (1934, 1475),
                (2081, 1170),
            ),
            "cross": (
                (988, 1694),
                (1065, 1697),
                (1142, 1700),
                (1218, 1705),
                (1295, 1710),
                (1372, 1714),
                (1423, 1717),
                (1730, 1380),
                (1745, 1452),
                (1752, 1497),
                (1765, 1560),
                (1778, 1600),
                (1793, 1630),
                (1805, 1660),
                (1830, 1692),
                (1880, 1710),
                (1910, 1693),
                (1932, 1676),
                (1945, 1661),
                (1958, 1640),
                (1971, 1617),
                (1985, 1590),
                (1997, 1563),
                (2010, 1535),
                (2023, 1495),
                (2034, 1453),
            ),
            "diamond": (
                (988, 1601),
                (1065, 1604),
                (1142, 1612),
                (1218, 1627),
                (1294, 1642),
                (1372, 1651),
                (1423, 1652),
                (1474, 1637),
                (1728, 1174),
                (1882, 1434),
                (2081, 1246),
            ),
            "square": (
                (988, 1729),
                (1065, 1733),
                (1142, 1737),
                (1218, 1741),
                (1295, 1748),
                (1372, 1751),
                (1423, 1755),
                (1742, 1645),
                (1752, 1692),
                (1795, 1825),
                (1810, 1842),
                (1830, 1860),
                (1882, 1860),
                (1912, 1840),
                (2083, 1540),
            ),
        },
        edge={"square": (1736, 14)},
        # The first markers on the left line (an x arm touches the 10^0 tick's stub).
        not_y_ticks=(1594.0, 1689.5, 1696.0, 1725.0, 1736.5, 1746.5),
    ),
}


@dataclass(frozen=True)
class Line:
    """A fitted frame line: position = slope * along + intercept, in pixels."""

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

    def ends(self) -> dict[str, float]:
        """Return the corrected positions of the four lines."""
        c0, r0 = self.corner
        return {
            "left": c0,
            "right": self.corrected(self.right.at(r0), r0)[0],
            "bottom": r0,
            "top": self.corrected(c0, self.top.at(c0))[1],
        }


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


def _find_ticks(
    dark: np.ndarray, frame: Frame, panel: Panel
) -> tuple[list[float], list[float]]:
    """Return the corrected positions of the tick stubs of the bottom and left lines."""
    inner, outer = TICK_INSET
    left, right, top, bottom = panel.frame
    columns = range(left + 12, right - 11)
    found = []
    for column in columns:
        row = round(frame.bottom.at(column))
        found.append(bool(dark[row - outer : row - inner + 1, column].all()))
    x_found = []
    for a, b in _spans(np.array(found)):
        column = (a + b) / 2 + columns.start
        x_found.append(frame.corrected(column, frame.bottom.at(column))[0])
    rows = range(top + 12, bottom - 11)
    found = []
    for row in rows:
        column = round(frame.left.at(row))
        found.append(bool(dark[row, column + inner : column + outer + 1].all()))
    y_found = []
    for a, b in _spans(np.array(found)):
        row = (a + b) / 2 + rows.start
        y_found.append(frame.corrected(frame.left.at(row), row)[1])
    return x_found, y_found


def _match(found: list[float], expected: list[float], what: str) -> None:
    """Raise unless the found positions are the expected ones, each within tolerance."""
    if len(found) != len(expected) or any(
        abs(f - e) > TICK_TOLERANCE
        for f, e in zip(sorted(found), sorted(expected), strict=True)
    ):
        found_ = [round(f, 1) for f in sorted(found)]
        msg = f"{what}: found {found_}, expected {sorted(expected)}"
        raise ValueError(msg)


def _axes(dark: np.ndarray, panel: Panel) -> tuple[Frame, Axis, Axis]:
    """Fit a panel's frame, check its ticks and return its calibrated axes."""
    frame = _frame(dark, panel)
    x_found, y_found = _find_ticks(dark, frame, panel)
    name = f"Fig. {panel.figure}"
    ends = frame.ends()
    for position, _, where in (*panel.x_ticks, *panel.y_ticks):
        if where != "tick" and abs(ends[where] - position) > TICK_TOLERANCE:
            msg = f"{name}: the {where} line at {ends[where]:.1f}, not {position}"
            raise ValueError(msg)
    _match(
        x_found,
        [p for p, _, w in panel.x_ticks if w == "tick"] + list(panel.not_x_ticks),
        f"{name}: x ticks",
    )
    _match(
        y_found,
        [p for p, _, w in panel.y_ticks if w == "tick"] + list(panel.not_y_ticks),
        f"{name}: y ticks",
    )
    x_axis = Axis(
        tuple(p for p, _, _ in panel.x_ticks),
        tuple(v for _, v, _ in panel.x_ticks),
        tolerance=X_TOLERANCE,
    )
    y_axis = Axis(
        tuple(p for p, _, _ in panel.y_ticks),
        tuple(v for _, v, _ in panel.y_ticks),
        log=panel.log,
        tolerance=Y_TOLERANCE,
    )
    x_axis.check()
    y_axis.check()
    return frame, x_axis, y_axis


@functools.cache
def _regions() -> dict[str, tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...]]]:
    """Return each marker shape's regions: (dark ones, light ones).

    Offsets from the candidate centre in pixels (dx right, dy down), on a 37 x 37
    stencil. A + is a vertical and a horizontal bar about 22 pixels long and 4 wide;
    an x the same on the diagonals; a diamond is an outline about 24 pixels across, its
    strokes 4 to 6 pixels wide; a square's sides are 3 to 4 pixels wide, their centres
    10 to 11 pixels from the marker's. The light regions of a + or an x come in
    pairs, the ends of one bar (or diagonal) each: one pair at least must be light.
    """
    dy, dx = np.mgrid[-18:19, -18:19]
    ax, ay = np.abs(dx), np.abs(dy)
    bar = (ax <= 1, ay <= 1)
    plus_dark = (
        bar[0] & (dy <= -5) & (dy >= -9),
        bar[0] & (dy >= 5) & (dy <= 9),
        bar[1] & (dx <= -5) & (dx >= -9),
        bar[1] & (dx >= 5) & (dx <= 9),
    )
    plus_light = (
        bar[0] & (dy <= -14) & (dy >= -16),
        bar[0] & (dy >= 14) & (dy <= 16),
        bar[1] & (dx <= -14) & (dx >= -16),
        bar[1] & (dx >= 14) & (dx <= 16),
    )
    diagonal = (np.abs(dx - dy) <= 1, np.abs(dx + dy) <= 1)
    cross_dark = tuple(
        diagonal[sx * sy < 0] & (sx * dx >= 5) & (sx * dx <= 8) & (sy * dy >= 4)
        for sx in (1, -1)
        for sy in (1, -1)
    )
    cross_light = tuple(
        diagonal[sx * sy < 0] & (sx * dx >= 12) & (sx * dx <= 14) & (sy * dy >= 11)
        for sx, sy in ((1, 1), (-1, -1), (1, -1), (-1, 1))
    )
    diamond_dark = tuple(
        (ax + ay <= 12) & (ax + ay >= 7) & (sx * dx >= 3) & (sy * dy >= 3)
        for sx in (1, -1)
        for sy in (1, -1)
    )
    side = (ax <= 7, ay <= 7)
    square_dark = (
        side[0] & (dy <= -9) & (dy >= -12),
        side[0] & (dy >= 9) & (dy <= 12),
        side[1] & (dx <= -9) & (dx >= -12),
        side[1] & (dx >= 9) & (dx <= 12),
    )
    return {
        "plus": (plus_dark, plus_light),
        "cross": (cross_dark, cross_light),
        "diamond": (diamond_dark, ()),
        "square": (square_dark, ()),
    }


def _scores(box: np.ndarray, shape: str) -> np.ndarray:
    """Return a shape's score at every pixel of a box: its darkest-region minimum.

    For a + or an x, zero unless one of its pairs of light regions (its two bars'
    ends) is light (each at most :data:`LIGHT` dark).
    """
    dark, light = _regions()[shape]

    def fraction(region: np.ndarray) -> np.ndarray:
        kernel = region.astype(float) / region.sum()
        return ndimage.correlate(box, kernel, mode="constant")

    score = np.min([fraction(r) for r in dark], axis=0)
    if light:
        ends = [fraction(r) <= LIGHT for r in light]
        score = np.where((ends[0] & ends[1]) | (ends[2] & ends[3]), score, 0.0)
    return score


@dataclass(frozen=True)
class Marker:
    """One marker read: its shape, centre (image column, row), values and score."""

    shape: str
    column: float
    row: float
    omega: float
    sigma: float
    score: float


@dataclass(frozen=True)
class Reading:
    """A panel's markers, and its calibration."""

    panel: Panel
    frame: Frame
    axes: tuple[Axis, Axis]
    markers: dict[str, list[Marker]]


@functools.cache
def _page(page: int) -> np.ndarray:
    """Return a page's image as dark (True) and light pixels."""
    grey = image(DIRECTORY, page)
    if grey.shape != SHAPE:
        msg = f"page {page}: an image of {grey.shape}, not {SHAPE}"
        raise ValueError(msg)
    return grey < DARK


def _cleared(
    dark: np.ndarray, frame: Frame, panel: Panel
) -> tuple[np.ndarray, int, int]:
    """Return the panel's box as floats, all but the frame's inside erased.

    With the box's first row and column in the page image.
    """
    left, right, top, bottom = panel.frame
    r0, c0 = top - 30, left - 30
    box = dark[r0 : bottom + 31, c0 : right + 31]
    rows, columns = np.indices(box.shape).astype(float)
    rows, columns = rows + r0, columns + c0
    m = FRAME_CLEAR

    def line(fitted: Line, along: np.ndarray) -> np.ndarray:
        return fitted.slope * along + fitted.intercept

    inside = (
        (columns > line(frame.left, rows) + m)
        & (columns < line(frame.right, rows) - m)
        & (rows > line(frame.top, columns) + m)
        & (rows < line(frame.bottom, columns) - m)
    )
    return np.where(inside, box, False).astype(float), r0, c0


def _place(
    cleared: tuple[np.ndarray, int, int], shape: str, seed: tuple[int, int]
) -> tuple[float, float, float]:
    """Return a seeded marker's centre (column, row) and best score.

    Raises if the best score within :data:`REACH` of the seed is below the shape's
    :data:`THRESHOLD`.
    """
    column, row = seed
    page, r0, c0 = cleared
    half = REACH + 19
    r, c = row - r0, column - c0
    box = page[r - half : r + half + 1, c - half : c + half + 1]
    score = _scores(box, shape)[
        half - REACH : half + REACH + 1, half - REACH : half + REACH + 1
    ]
    best = float(score.max())
    if best < THRESHOLD[shape]:
        msg = f"the {shape} seeded at {seed}: best score {best:.2f}"
        raise ValueError(msg)
    rb, cb = np.unravel_index(int(np.argmax(score)), score.shape)
    near = np.zeros_like(score, dtype=bool)
    near[max(rb - 2, 0) : rb + 3, max(cb - 2, 0) : cb + 3] = True
    rr, cc = np.nonzero(near & (score >= best - PLATEAU))
    return (
        column - REACH + float(cc.mean()),
        row - REACH + float(rr.mean()),
        best,
    )


def _edge_row(
    dark: np.ndarray, frame: Frame, panel: Panel, edge: tuple[int, int]
) -> float:
    """Return the row of a marker on the left frame line, from its half inside.

    The plot is clipped at the frame: a marker centred on the left line shows only its
    right half. Its row is the middle of the dark rows within ``half`` of the seed row
    in the five columns from 5 to 9 pixels right of the line's centre (a diamond's
    right half, or a + bar's right arm with the line leaving it). Raises if a tick of
    the axis is within :data:`EDGE_CLEAR` rows (its stub would be read with it).
    """
    seed, half = edge
    near = [
        p for p, _, w in panel.y_ticks if w == "tick" and abs(p - seed) < EDGE_CLEAR
    ]
    if near:
        msg = f"Fig. {panel.figure}: a tick at {near} next to the edge marker at {seed}"
        raise ValueError(msg)
    line = round(frame.left.at(seed))
    patch = dark[seed - half : seed + half + 1, line + 5 : line + 10]
    rows = np.flatnonzero(patch.any(axis=1))
    if rows.size == 0 or rows[0] == 0 or rows[-1] == 2 * half:
        msg = f"Fig. {panel.figure}: the edge marker at row {seed} is not isolated"
        raise ValueError(msg)
    return seed - half + (rows[0] + rows[-1]) / 2


def _read(dark: np.ndarray, panel: Panel) -> Reading:
    """Return a panel's markers, sorted by omega, per shape."""
    frame, x_axis, y_axis = _axes(dark, panel)
    cleared = _cleared(dark, frame, panel)
    markers: dict[str, list[Marker]] = {shape: [] for shape in panel.curves}
    for shape, seeds in panel.seeds.items():
        for seed in seeds:
            column, row, score = _place(cleared, shape, seed)
            x, y = frame.corrected(column, row)
            markers[shape].append(
                Marker(shape, column, row, float(x_axis(x)), float(y_axis(y)), score)
            )
    for shape, seed in panel.edge.items():
        row = _edge_row(dark, frame, panel, seed)
        column = frame.left.at(row)
        x, y = frame.corrected(column, row)
        markers[shape].append(
            Marker(shape, column, row, float(x_axis(x)), float(y_axis(y)), 1.0)
        )
    for shape, found in markers.items():
        found.sort(key=lambda k: (k.omega, k.sigma))
        for a, b in itertools.pairwise(found):
            if math.hypot(a.column - b.column, a.row - b.row) < CLASH:
                msg = f"Fig. {panel.figure}: two {shape} markers at {a.column, a.row}"
                raise ValueError(msg)
    return Reading(panel, frame, (x_axis, y_axis), markers)


@functools.cache
def reading(key: str) -> Reading:
    """Return a panel's reading."""
    panel = PANELS[key]
    return _read(_page(panel.page), panel)


def table_i() -> dict[str, dict[str, float]]:
    """Return Table I's binding energies (eV) per gas: DF and experimental, by shell."""
    rows = json_table(DIRECTORY, TABLE_I)[1:]
    out: dict[str, dict[str, float]] = {}
    gas = ""
    for element, shell, df, df_split, expt, expt_split in rows:
        gas = element or gas
        entry = out.setdefault(gas, {})
        if shell:
            j = shell.split("_")[-1].strip("{}") if "_" in shell else ""
            entry[f"df{j}"] = float(df)
            entry[f"expt{j}"] = float(expt)
        else:
            entry["df_split"] = float(df_split)
            entry["expt_split"] = float(expt_split)
    return out


def extract() -> dict[str, list[Point]]:
    """Return the markers of Figs. 3 to 6, after the checks."""
    checks()
    files: dict[str, list[Point]] = {}
    for key, panel in PANELS.items():
        found = reading(key)
        for shape, file in panel.curves.items():
            files[file] = [
                Point(
                    f"{k:02d}",
                    f"{m.omega:.4f}",
                    f"{m.sigma:.3g}" if panel.log else f"{m.sigma:.3f}",
                    "",
                    f"marker centre at pixel ({m.column:.1f}, {m.row:.1f})",
                )
                for k, m in enumerate(found.markers[shape], 1)
            ]
    return files


def _near(markers: list[Marker], omega: float) -> Marker | None:
    """Return the marker at a computed frequency (within :data:`SAME`), if read."""
    found = [m for m in markers if abs(m.omega - omega) <= SAME]
    return min(found, key=lambda m: abs(m.omega - omega)) if found else None


def _redrawn(key: str) -> tuple[np.ndarray, np.ndarray]:
    """Return Lagutin et al.'s redrawing of Fig. 5 (their Fig. 8, blue dashes).

    Photon energies (eV) and values (1e-50 cm^4 s) of every vertex of the blue
    subpaths in the panel's frame, less the legend's sample, sorted by energy; read
    with lagutin_2017's machinery (TC-880).
    """
    panel = lagutin_2017.PANELS[key]
    box = lagutin_2017.frame(key)
    x_axis, y_axis = lagutin_2017.axes(key)
    strokes = [
        st
        for st in lagutin_2017._strokes(panel.page)
        if not st.filled
        and st.colour == lagutin_2017.BLUE
        and np.any(lagutin_2017._inside(st.vertices, box))
        and not any(np.all(lagutin_2017._inside(st.vertices, b)) for b in panel.boxes)
    ]
    vertices = np.concatenate(
        [st.vertices[lagutin_2017._inside(st.vertices, box, 0.0)] for st in strokes]
    )
    vertices = vertices[np.argsort(vertices[:, 0], kind="stable")]
    energy = np.asarray(x_axis(vertices[:, 0]), dtype=float)
    value = np.asarray(y_axis(vertices[:, 1]), dtype=float) * 10.0  # 1e-49 to 1e-50
    return energy, value


def checks() -> dict[str, object]:  # noqa: PLR0912, PLR0915 (one list of checks)
    """Return the checks of the reading; raise if one fails its margin.

    Table I's splittings; each curve's first frequency against Table I's thresholds;
    the two polarisations' common frequencies; Fig. 6's partials against Fig. 5's
    totals; Fig. 5 against Lagutin et al.'s redrawing of it; Fig. 3 against Starace and
    Jiang's ESI totals. Reported only: the first frequencies against half the
    Dirac-Fock 2P1/2 energy, Fig. 3 against McKenna and van der Hart (2004), Fig. 4
    against McGuire (1981).
    """
    out: dict[str, object] = {}
    failures: list[str] = []
    energies = table_i()
    for gas, entry in energies.items():
        if "df_split" not in entry:
            continue
        df = entry["df1/2"] - entry["df3/2"]
        expt = entry["expt1/2"] - entry["expt3/2"]
        out[f"Table I {gas} splittings (DF, expt)"] = (round(df, 2), round(expt, 2))
        if (
            abs(df - entry["df_split"]) > 0.0101
            or abs(expt - entry["expt_split"]) > 0.0101
        ):
            failures.append(f"Table I {gas}: splittings {df:.2f}, {expt:.2f}")
    for key, panel in PANELS.items():
        found = reading(key)
        entry = energies[panel.gas]
        half_expt = entry["expt1/2"] / 2 / HARTREE_EV
        half_df = entry["df1/2"] / 2 / HARTREE_EV
        one_photon = min(entry["df3/2"], entry["expt3/2"]) / HARTREE_EV
        for shape, markers in found.markers.items():
            first, last = markers[0].omega, markers[-1].omega
            out[f"Fig. {panel.figure} {shape} first omega, half 2P1/2 expt and DF"] = (
                round(first, 4),
                round(half_expt, 4),
                round(half_df, 4),
            )
            if first < half_expt - SAME or last >= one_photon:
                failures.append(
                    f"Fig. {panel.figure} {shape}: {first:.4f} to {last:.4f}"
                )
        shapes = list(found.markers)
        for a, b in itertools.pairwise(shapes):
            offsets = [
                min(abs(k.omega - m.omega) for k in found.markers[b])
                for m in found.markers[a]
            ]
            out[f"Fig. {panel.figure} {a}/{b}: markers, at one frequency"] = (
                len(offsets),
                sum(o <= SAME for o in offsets),
            )
    totals, partials = reading("fig5").markers, reading("fig6").markers
    for total, parts in (
        ("plus", ("plus", "cross")),
        ("diamond", ("diamond", "square")),
    ):
        ratios = []
        for m in totals[total]:
            found_parts = [_near(partials[p], m.omega) for p in parts]
            if all(found_parts):
                ratios.append(sum(f.sigma for f in found_parts if f) / m.sigma)
        out[f"Fig. 6 {'+'.join(parts)} / Fig. 5 {total}"] = (
            len(ratios),
            round(min(ratios), 3),
            round(max(ratios), 3),
        )
        if not ratios or max(abs(r - 1) for r in ratios) > SUM:
            failures.append(f"Fig. 6 partials against Fig. 5 {total}: {ratios}")
    for key, shape, windows in (
        ("fig8a", "plus", ((6.80, 8.10), (8.70, 9.28))),
        ("fig8b", "diamond", ((6.80, 8.10), (8.70, 9.28))),
    ):
        energy, value = _redrawn(key)
        ratios = []
        for m in totals[shape]:
            e = m.omega * HARTREE_EV
            if any(low <= e <= high for low, high in windows):
                redrawn = 10 ** np.interp(e, energy, np.log10(value))
                ratios.append(m.sigma / float(redrawn))
        out[f"Fig. 5 {shape} / Lagutin Fig. 8 redrawing"] = (
            len(ratios),
            round(min(ratios), 3),
            round(max(ratios), 3),
        )
        if not ratios or max(abs(r - 1) for r in ratios) > REDRAWN:
            failures.append(f"Fig. 5 {shape} against Lagutin's redrawing: {ratios}")
    argon = reading("fig3").markers
    for shape, file, window, margin in (
        ("plus", "Ar/starace_1987__fig7_esi_linear_Ar.csv", (8.1, 9.4), ESI),
        ("diamond", "Ar/starace_1987__fig6_esi_circular_Ar.csv", (8.1, 9.4), ESI),
        ("plus", "Ar/mckenna_2004__fig5_sigma2_Ar.csv", (8.0, 10.0), None),
    ):
        x, y = _curve(file)
        ratios = [
            m.sigma * 1e-50 / float(np.interp(m.omega * HARTREE_EV, x, y))
            for m in argon[shape]
            if window[0] <= m.omega * HARTREE_EV <= window[1]
        ]
        out[
            f"Fig. 3 {shape} / {file.split('__')[0][3:]} {file.split('__')[1][:-4]}"
        ] = (
            len(ratios),
            round(min(ratios), 3),
            round(max(ratios), 3),
        )
        if margin is not None and max(abs(r - 1) for r in ratios) > margin:
            failures.append(f"Fig. 3 {shape} against {file}: {ratios}")
    x, y = _curve("Kr/mcguire_1981__fig3_sigma2_Kr.csv")
    energy = HC_EV_ANGSTROM / x
    order = np.argsort(energy)
    ratios = [
        (
            round(m.omega * HARTREE_EV, 3),
            round(
                m.sigma
                * 1e-50
                / float(np.interp(m.omega * HARTREE_EV, energy[order], y[order])),
                2,
            ),
        )
        for m in reading("fig4").markers["diamond"]
        if energy.min() <= m.omega * HARTREE_EV <= 9.8
    ]
    out["Fig. 4 diamond / mcguire_1981 fig3_sigma2_Kr (eV, ratio)"] = ratios
    if failures:
        raise ValueError("; ".join(failures))
    return out


def _curve(file: str) -> tuple[np.ndarray, np.ndarray]:
    """Return another record's data file: its x column and values times its scale."""
    with path(file).open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    scale = 1e-50 if "lhuillier" in file else 1.0
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["value"]) for r in rows]) * scale
    return x, y


def main() -> None:
    """Write the data files, and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"{label}: {result}")


if __name__ == "__main__":
    main()
