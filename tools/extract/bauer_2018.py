# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bauer and Deng (2018): the TDSE ionisation of H(1s) in a flat-top pulse.

The paper prints no table of its calculations. The time-dependent Schrodinger equation
(TDSE) probabilities of atomic hydrogen after a '2-10-2' pulse (sin^2 turn-on and
turn-off of the vector potential over two cycles each, a ten-cycle plateau) are given
only in figures, all vector drawings placed as form XObjects (matplotlib's: the labels
are text, the curves polylines, the markers four-arc circles or squares):

* **Fig. 1** (journal p. 4, PDF page 5, library figure ``fig_917161e2757b``) and
  **Fig. 2** (``fig_c26d0459a537``): ground-state (black open circles), ionisation
  (blue filled circles) and excitation (red squares) probabilities against the peak
  intensity, linear axes, omega = 0.057 and 0.114 au: 'Each circle or square
  corresponds to a single TDSE simulation. B-splines are used to display the curves.'
* **Fig. 5** (journal p. 5, PDF page 6, ``fig_e778f524ee82``) and **Fig. 6**
  (``fig_79b9623bd88f``): the TDSE ionisation probability (grey filled circles on a grey
  B-spline curve) with five models' curves, (a) log-log, (b) linear, at 0.057 and
  0.114 au.
* **Fig. 7** (journal p. 5, PDF page 6, ``fig_ab49da534666``) and **Fig. 8** (journal
  p. 6, PDF page 7, ``fig_a558d83e4736``): the three probabilities against the
  frequency omega at F = 0.02 and 0.063 au, log-log, as B-spline curves only (blue solid
  ionisation, red dashed excitation, black dotted ground state); the simulations'
  frequencies are not shown.

What is read (:func:`_figure.painted`, :func:`tokens`):

* the frames: four stroked lines (bottom, left, top, right) of the axes' colour; the
  ticks are two-vertex lines of that colour perpendicular to the bottom and left edges
  with one end on them, inward; the top and right edges carry none (checked: any there
  must repeat one below or left). The dashed vertical line at I_BSI in Figs. 5 and 6
  also starts on the bottom edge: a tick is the only segment longer than :data:`DOT`
  at its position (a dotted curve's dots, shorter, may share it).
* the axes: the labels are text, read character by character from the page's text
  (:func:`tokens`); every label lies within :data:`LABEL_AT` of a long tick (x: the
  whole label's middle; y: the base's). Linear axes: the ticks are equally spaced in
  value, the labelled ones fix the step. Logarithmic axes: the long ticks are the
  decades, consecutive, and the short ticks 2 to 9 times the decade below (those before
  the first decade counted down from 9). One labelled decade suffices there: the
  counted values must lie on one log scale within :data:`TICK_TOLERANCE`
  (:meth:`_figure.Axis.check`), which a wrong count misses by points. Figs. 7 and 8's x
  axes are labelled at 10^-1 only; the excitation humps the text places at omega =
  0.375, 0.444, 0.114 and 0.228 au fall there (checked).
* the markers (:data:`SERIES`): each drawn twice, filled and stroked, at one centre
  (checked); the centre of the stroked path's box. No legend has a marker.
* the curves: Figs. 7 and 8's blue ionisation curve, polylines only (checked: no Bezier
  segment), in pieces that join end to start (checked); every vertex inside the frame
  is kept, in order, except a vertex repeating the one before it (a zero-length
  segment, or the joint of two pieces). Fig. 7's curve has 641 vertices, Fig. 8's
  about 4700 (0.03 to 0.16 pt apart: the B-spline drawn densely; kept, not resampled).

Every point's ``printed`` cell gives the panel and the marker's centre or the vertex in
page points.

**Checks** (:func:`checks`), against every number the paper prints for these quantities
and against the figures' other copies of the same data:

* Fig. 5(b) and Fig. 6(b) put their markers where Figs. 1 and 2 put theirs (the same
  intensities and probabilities, within :data:`LINEAR_INTENSITY` and
  :data:`LINEAR_AGREEMENT`), as the captions say.
* Fig. 6(a)'s markers are Fig. 2's above 7e12 W/cm^2 (the intensities within
  :data:`GRID`, the probabilities within :data:`LOG_LINEAR` where P >= 0.03).
* Fig. 5(a)'s markers are NOT Fig. 1's, though the caption says so: other intensities
  and, at the intensities both have, other probabilities (:func:`fig5a_against_fig1`,
  reported). The grey curves are smoothing B-splines: none passes through its panel's
  markers (:func:`on_curve`, reported: up to 6 % off where P >= 0.1, more below).
* in Figs. 1 and 2 the ground-state, excitation and ionisation probabilities sum to 1
  within :data:`SUM_TOLERANCE`.
* the dashed I_BSI lines of Figs. 5 and 6 lie at the printed 1.37 x 10^14 W/cm^2 (the
  intensity of F_BSI = 1/16 au, I = F^2 in atomic units), within :data:`BSI`; and the
  text's pairs of field and intensity (:data:`TEXT_FIELDS`) agree with I = F^2.
* the text's statements on the figures: excitation 'less than 10%' at 0.057 au (Fig. 1);
  'up to 30% near I = 1.5 x 10^14' at 0.114 au (Fig. 2); 'nearly fully ionized' at 4.0
  x 10^14 (Figs. 1 and 2, >= 0.97); Figs. 7 and 8's excitation humps near omega =
  0.375 and 0.444 (Fig. 7) and 0.114 and 0.228 au (Fig. 8) within :data:`HUMP`, and
  Fig. 8's excitation above the ionisation over '~0.37 a.u. <= omega <= ~0.47 a.u.'.

:func:`report` prints the numbers, and the comparisons of Figs. 7 and 8 at omega =
0.057 and 0.114 au with Figs. 5(a), 1 and 6(a) at the same field.
"""

import ctypes
import functools
import itertools
import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, Painted, painted, pdf
from _library import Point, write
from pypdfium2 import raw

from optmat.constants import AU_ENERGY, AU_FIELD, EPS0, HBAR, C

DIRECTORY = "Bauer_2018_NumericalTestsTheoretical"
SOURCE = "bauer_2018"
BLACK = (3, 3, 3)  # the axes, the ground state and Fig. 5's 'Popruzhenko et al.'
WHITE = (255, 255, 255)
BLUE = (34, 86, 166)
RED = (238, 41, 40)
GREY = (126, 128, 131)
I_AU = 0.5 * C * EPS0 * AU_FIELD**2 / 1e4  # W/cm^2 for I = F^2 = 1 in atomic units
AU_TIME_FS = HBAR / AU_ENERGY * 1e15
STRAIGHT = 0.01  # an edge or a tick is vertical or horizontal to 0.01 pt
EDGE = 40.0  # a frame's edge is longer than 40 pt; ticks are under 2.5 pt
ON_LINE = 0.03  # a tick ends on its edge within 0.03 pt
TICK_LONGEST = 4.0
LONG_TICK = 1.5  # a long tick is about twice the short ones (1.96 to 1.98 times)
MIRROR = 0.03
# A dash of the I_BSI line is 1.36 pt long; the dotted curves' dots 0.34 pt or less.
DOT = 0.5
# Every tick lies on its fitted scale within 0.05 pt (the drawing's coordinates are
# quantised in steps of about 0.03 pt; the largest residual is printed by report()).
TICK_TOLERANCE = 0.05
# A label's middle (x) or its base's (y) lies within 3 pt of its long tick (found: up to
# 2.5 pt, the x labels '10^13' of Fig. 5(a); the long ticks are 10 pt apart or more).
LABEL_AT = 3.0
LABEL_BESIDE = 16.0  # the labels lie within 16 pt below or left of the frame
ON_EDGE = 0.05  # a marker or vertex is in the frame if within 0.05 pt of it
TWIN = 0.005  # a marker's fill and stroke share a centre to 0.005 pt
# Two copies of a linear plot agree within 0.0005 in probability (the y scales are 101
# to 111 pt per unit, the quantum 0.03 pt, 0.0003) and 1e11 W/cm^2 in intensity (the x
# scales are 42 to 44 pt per 1e14 W/cm^2: 0.03 pt is 7e10).
LINEAR_AGREEMENT = 0.0005
LINEAR_INTENSITY = 1e11
# Fig. 6(a)'s intensities are Fig. 2's within 0.3 % or, at the lowest, Fig. 2's 1e11.
GRID = 0.003
# Fig. 6(a) (log) and Fig. 2 (linear) agree within 2 % where P >= 0.03 (Fig. 2's
# 0.0003 resolution is 1 % there).
LOG_LINEAR = 0.02
LOG_LINEAR_FROM = 0.03
# The grey curves are compared with the markers where P >= 0.02 in the linear panels
# (their resolution).
ON_CURVE_FROM = 0.02
SUM_TOLERANCE = 0.002  # the three probabilities sum to 1 within 0.002 (y offsets)
BSI = 0.005  # the I_BSI lines lie within 0.5 % of 1.37e14
EXCITED_RANGE = (0.37, 0.47)  # au, Fig. 8's excitation above the ionisation
HUMP = 0.05  # an excitation hump lies within 5 % of the text's omega ('near')
F_BSI = 1 / 16  # the barrier-suppression field of H(1s), E_B^2/4Z (Eq. 1), au
MARKER_LARGEST = 3.0  # pt; the markers are 2.4 pt across or less
NEARLY_FULL = 0.97  # 'nearly fully ionized' at 4e14 W/cm^2
EXCITATION_BELOW = 0.1  # 'less than 10%' (Fig. 1)
# 'up to 30% near I = 1.5e14' (Fig. 2): the largest is 0.34 at 1.7e14 (0.29 at 1.5e14)
EXCITATION_PEAK = (0.25, 0.35, 1.3e14, 1.8e14)
# The field of each intensity the text prints is the field it prints with it within half
# its last digit, 0.0005 au.
HALF_DIGIT = 0.0005


@dataclass(frozen=True)
class Panel:
    """A panel: name, page (1-based), frame's lower left corner (pt, to 1 pt), axes."""

    name: str
    page: int
    corner: tuple[float, float]
    x_log: bool
    y_log: bool


PANELS = {
    "fig1": Panel("Fig. 1", 5, (97.0, 633.2), x_log=False, y_log=False),
    "fig2": Panel("Fig. 2", 5, (97.9, 418.9), x_log=False, y_log=False),
    "fig5a": Panel("Fig. 5(a)", 6, (97.1, 544.8), x_log=True, y_log=True),
    "fig5b": Panel("Fig. 5(b)", 6, (96.9, 395.5), x_log=False, y_log=False),
    "fig6a": Panel("Fig. 6(a)", 6, (357.5, 487.9), x_log=True, y_log=True),
    "fig6b": Panel("Fig. 6(b)", 6, (357.4, 337.0), x_log=False, y_log=False),
    "fig7": Panel("Fig. 7", 6, (364.6, 139.7), x_log=True, y_log=True),
    "fig8": Panel("Fig. 8", 7, (103.7, 668.2), x_log=True, y_log=True),
}
# The linear panels' x axes are in units of 10^14 W/cm^2 (their titles).
X_SCALE = {"fig1": 1e14, "fig2": 1e14, "fig5b": 1e14, "fig6b": 1e14}


@dataclass(frozen=True)
class Series:
    """A marker series: stroke and fill colours, vertex count, size range (pt)."""

    stroke: tuple[int, int, int]
    fill: tuple[int, int, int]
    vertices: int
    size: tuple[float, float]


# The stroked paths' boxes (the line's width included): the circles 2.19 to 2.36 pt
# across, the squares and the ground state's open circles smaller.
IONISED = Series(BLUE, BLUE, 14, (2.1, 2.3))
GROUND = Series(BLACK, WHITE, 14, (1.0, 2.3))
EXCITED = Series(RED, RED, 5, (0.8, 2.3))
TDSE = Series(GREY, GREY, 14, (2.1, 2.4))
SERIES = {
    "fig1": IONISED,
    "fig2": IONISED,
    "fig5a": TDSE,
    "fig5b": TDSE,
    "fig6a": TDSE,
    "fig6b": TDSE,
}
DATASETS = {
    "fig5a": f"H/{SOURCE}__fig5a_tdse_H_800nm.csv",
    "fig1": f"H/{SOURCE}__fig1_tdse_H_800nm.csv",
    "fig6a": f"H/{SOURCE}__fig6a_tdse_H_400nm.csv",
    "fig7": f"H/{SOURCE}__fig7_tdse_H_spectrum_F002.csv",
    "fig8": f"H/{SOURCE}__fig8_tdse_H_spectrum_F0063.csv",
}
FIELD = {"fig7": 0.02, "fig8": 0.063}  # au, printed


@dataclass(frozen=True)
class Frame:
    """A panel's frame (l, b, r, t, pt), its axes, and the page's paths around it."""

    box: tuple[float, float, float, float]
    x: Axis
    y: Axis
    drawn: tuple[Painted, ...]


def _region(panel: Panel) -> tuple[float, float, float, float]:
    """Return the region (l, b, r, t, pt) that holds a panel and its labels."""
    x0, y0 = panel.corner
    return (x0 - 30, y0 - 20, x0 + 190, y0 + 130)


def _box(
    drawn: Sequence[Painted], corner: tuple[float, float], where: str
) -> tuple[float, float, float, float]:
    """Return the frame (l, b, r, t) whose lower left corner is within 1 pt of this."""
    edges = []
    for d in drawn:
        v = d.vertices
        if d.stroke != BLACK or d.fill is not None or len(v) != 2:
            continue
        extent = v.max(axis=0) - v.min(axis=0)
        if max(extent) > EDGE and min(extent) < STRAIGHT:
            edges.append(v)
    x0, y0 = corner
    near = [
        e for e in edges if abs(e[:, 0].min() - x0) < 1 and abs(e[:, 1].min() - y0) < 1
    ]
    bottoms = [e for e in near if np.ptp(e[:, 1]) < STRAIGHT]
    lefts = [e for e in near if np.ptp(e[:, 0]) < STRAIGHT]
    if len(bottoms) != 1 or len(lefts) != 1:
        msg = f"{where}: {len(bottoms)} bottom and {len(lefts)} left edges"
        raise ValueError(msg)
    left, right = float(bottoms[0][:, 0].min()), float(bottoms[0][:, 0].max())
    bottom, top = float(lefts[0][:, 1].min()), float(lefts[0][:, 1].max())
    for fixed, value in ((0, right), (1, top)):
        if not any(
            np.ptp(e[:, fixed]) < STRAIGHT and abs(e[0, fixed] - value) < STRAIGHT
            for e in edges
        ):
            msg = f"{where}: the {'right' if fixed == 0 else 'top'} edge is missing"
            raise ValueError(msg)
    return left, bottom, right, top


def _edge_ticks(
    drawn: Sequence[Painted],
    box: tuple[float, float, float, float],
    axis: int,
    *,
    far: bool,
) -> list[tuple[float, float]]:
    """Return an edge's ticks as (position, length), sorted.

    ``axis`` 0: the bottom (or, ``far``, top) edge's vertical ticks; 1: the left (or
    right) edge's horizontal ticks.
    """
    along, across = (0, 1) if axis == 0 else (1, 0)
    edge = box[across + 2] if far else box[across]
    low, high = box[along], box[along + 2]
    lines = [
        d.vertices
        for d in drawn
        if d.stroke == BLACK
        and len(d.vertices) == 2
        and abs(d.vertices[0, along] - d.vertices[1, along]) <= STRAIGHT
        and np.ptp(d.vertices[:, across]) > 0
    ]
    ticks = {}
    for v in lines:
        ends = sorted(v[:, across], key=lambda e: abs(e - edge))
        length = abs(ends[1] - ends[0])
        position = float(v[0, along])
        if not (
            abs(ends[0] - edge) < ON_LINE
            and 0 < length < TICK_LONGEST
            and low - ON_LINE <= position <= high + ON_LINE
        ):
            continue
        dashes = sum(
            1
            for w in lines
            if abs(float(w[0, along]) - position) < STRAIGHT
            and DOT < np.ptp(w[:, across]) < TICK_LONGEST
        )
        if dashes == 1:
            ticks[round(position, 3)] = round(float(length), 3)
    return sorted(ticks.items())


@functools.cache
def _characters(page: int) -> tuple[tuple[str, float, float, float, float], ...]:
    """Return a page's characters with their boxes (l, b, r, t, pt), in text order."""
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    text = document[page - 1].get_textpage()
    out = []
    for index in range(text.count_chars()):
        char = text.get_text_range(index, 1)
        if not char.strip():
            out.append((" ", 0.0, 0.0, 0.0, 0.0))
            continue
        left, bottom, right, top = text.get_charbox(index)
        out.append((char, left, bottom, right, top))
    return tuple(out)


@dataclass(frozen=True)
class Token:
    """A run of characters on one baseline, under 2 pt apart, and its box (pt)."""

    text: str
    left: float
    bottom: float
    right: float
    top: float

    @property
    def centre(self) -> tuple[float, float]:
        """The centre of the box."""
        return (self.left + self.right) / 2, (self.bottom + self.top) / 2


def tokens(page: int, region: tuple[float, float, float, float]) -> list[Token]:
    """Return the text runs of a page whose characters lie in a region (l, b, r, t).

    A superscript, or a minus sign (set higher), starts a run of its own.
    """
    left, bottom, right, top = region
    runs: list[list] = []
    previous = None
    for char, west, south, east, north in _characters(page):
        x, y = (west + east) / 2, (south + north) / 2
        if char == " " or not (left <= x <= right and bottom <= y <= top):
            previous = None
            continue
        if (
            previous is not None
            and abs(south - previous[0]) < 0.5
            and 0 <= west - previous[1] < 2.0
        ):
            run = runs[-1]
            run[0] += char
            run[2], run[3], run[4] = min(run[2], south), east, max(run[4], north)
        else:
            runs.append([char, west, south, east, north])
        previous = (south, east)
    return [Token(r[0].replace("\N{MINUS SIGN}", "-"), *r[1:]) for r in runs]


def _labels(
    page: int, box: tuple[float, float, float, float], axis: int, *, log: bool
) -> list[tuple[float, float]]:
    """Return an axis's labels as (position along the axis, value), sorted."""
    left, bottom, right, top = box
    if axis == 0:
        region = (left - 8, bottom - LABEL_BESIDE, right + 8, bottom - 0.5)
    else:
        region = (left - LABEL_BESIDE - 8, bottom - 4, left - 0.5, top + 4)
    found = sorted(tokens(page, region), key=lambda t: t.left)
    out = []
    if not log:
        for t in found:
            try:
                out.append((t.centre[axis], float(t.text)))
            except ValueError:
                continue
        return sorted(out)
    consumed: set[int] = set()
    for index, base in enumerate(found):
        if base.text != "10" or index in consumed:
            continue
        # The exponent: the raised runs that follow the base, each within 2 pt of the
        # one before (a minus sign is a run of its own).
        exponent, end = "", base.right
        for later, t in enumerate(found):
            if 0 <= t.left - end < 2.0 and 1 < t.bottom - base.bottom < 5:
                exponent, end = exponent + t.text, t.right
                consumed.add(later)
        if not exponent:
            msg = f"page {page}: the label 10 at {base.centre} has no exponent"
            raise ValueError(msg)
        where = (base.left + end) / 2 if axis == 0 else base.centre[1]
        out.append((where, 10.0 ** int(exponent)))
    return sorted(out)


def _calibrate(
    ticks: list[tuple[float, float]],
    labels: list[tuple[float, float]],
    *,
    log: bool,
    where: str,
) -> Axis:
    """Return the axis calibrated from its ticks and labels, checked on every tick."""
    positions = [p for p, _ in ticks]
    shortest = min(length for _, length in ticks)
    long = [p for p, length in ticks if length > LONG_TICK * shortest]
    labelled = {}
    for position, value in labels:
        nearest = min(long, key=lambda p: abs(p - position))
        if abs(nearest - position) > LABEL_AT:
            msg = f"{where}: the label {value} at {position:.2f} is not at a long tick"
            raise ValueError(msg)
        labelled[positions.index(nearest)] = value
    if len(labelled) < (1 if log else 2):
        msg = f"{where}: {len(labelled)} labelled ticks"
        raise ValueError(msg)
    if not log:
        (k0, v0), (k1, v1) = sorted(labelled.items())[:2]
        step = (v1 - v0) / (k1 - k0)
        values = [v0 + (k - k0) * step for k in range(len(positions))]
    else:
        decades = [positions.index(p) for p in long]
        k0, v0 = min(labelled.items())
        d0 = round(math.log10(v0)) - decades.index(k0)
        values = []
        for k in range(len(positions)):
            if k < decades[0]:
                values.append((10 - (decades[0] - k)) * 10.0 ** (d0 - 1))
                continue
            below = max(i for i, d in enumerate(decades) if d <= k)
            values.append((k - decades[below] + 1) * 10.0 ** (d0 + below))
    for k, value in labelled.items():
        if not math.isclose(values[k], value, rel_tol=1e-9, abs_tol=1e-12):
            msg = f"{where}: tick {k} counts to {values[k]}, labelled {value}"
            raise ValueError(msg)
    axis = Axis(tuple(positions), tuple(values), log=log, tolerance=TICK_TOLERANCE)
    axis.check()
    return axis


@functools.cache
def frame(key: str) -> Frame:
    """Return a panel's frame and calibrated axes, its labels checked."""
    panel = PANELS[key]
    drawn = painted(DIRECTORY, panel.page, _region(panel))
    box = _box(drawn, panel.corner, panel.name)
    axes = []
    for axis, log in ((0, panel.x_log), (1, panel.y_log)):
        where = f"{panel.name} {'xy'[axis]}"
        ticks = _edge_ticks(drawn, box, axis, far=False)
        positions = [p for p, _ in ticks]
        for p, _ in _edge_ticks(drawn, box, axis, far=True):
            if min(abs(p - q) for q in positions) > MIRROR:
                msg = f"{where}: a tick of the opposite edge has none on this one"
                raise ValueError(msg)
        labels = _labels(panel.page, box, axis, log=log)
        axes.append(_calibrate(ticks, labels, log=log, where=where))
    return Frame(box, axes[0], axes[1], tuple(drawn))


def residual(key: str) -> float:
    """Return the largest distance (pt) of a panel's ticks from their fitted scales."""
    found = frame(key)
    worst = 0.0
    for axis in (found.x, found.y):
        values = np.asarray(axis.values)
        y = np.log10(values) if axis.log else values
        slope, intercept = np.polyfit(axis.positions, y, 1)
        distance = np.abs((y - intercept) / slope - np.asarray(axis.positions))
        worst = max(worst, float(np.max(distance)))
    return worst


def _inside(x: float, y: float, box: tuple[float, float, float, float]) -> bool:
    """Return whether a point lies in a frame, within ON_EDGE."""
    left, bottom, right, top = box
    return (
        left - ON_EDGE <= x <= right + ON_EDGE
        and bottom - ON_EDGE <= y <= top + ON_EDGE
    )


@dataclass(frozen=True)
class Reading:
    """One marker or vertex as read: its position (pt) and its values."""

    x: float
    y: float
    abscissa: float
    value: float


def _read(key: str, points: Sequence[tuple[float, float]]) -> list[Reading]:
    """Return points read on a panel's axes (intensity in W/cm^2), in order given."""
    found = frame(key)
    scale = X_SCALE.get(key, 1.0)
    return [
        Reading(x, y, float(found.x(x)) * scale, float(found.y(y))) for x, y in points
    ]


def markers(key: str, series: Series) -> list[Reading]:
    """Return a panel's markers of a series, in order of intensity.

    Each marker is a stroked path with a filled twin at its centre (checked).
    """
    found = frame(key)
    stroked = [
        d
        for d in found.drawn
        if d.stroke == series.stroke
        and len(d.vertices) == series.vertices
        and _inside(*d.centre, found.box)
    ]
    filled = [
        d.centre
        for d in found.drawn
        if d.fill == series.fill
        and d.stroke is None
        and len(d.vertices) == series.vertices
        and _inside(*d.centre, found.box)
    ]
    name = PANELS[key].name
    if len(filled) != len(stroked):
        msg = f"{name}: {len(stroked)} stroked markers, {len(filled)} filled"
        raise ValueError(msg)
    for d in stroked:
        if not series.size[0] <= d.size <= series.size[1]:
            msg = f"{name}: a marker {d.size:.2f} pt across at {d.centre}"
            raise ValueError(msg)
        if min(math.dist(d.centre, c) for c in filled) > TWIN:
            msg = f"{name}: the marker at {d.centre} has no filled twin"
            raise ValueError(msg)
    centres = sorted({(round(d.centre[0], 3), round(d.centre[1], 3)) for d in stroked})
    if len(centres) != len(stroked):
        msg = f"{name}: two markers share a centre"
        raise ValueError(msg)
    return _read(key, centres)


def _bezier_segments(page: int, colour: tuple[int, int, int], shortest: int) -> int:
    """Return the Bezier segments of a page's stroked paths of a colour and length."""
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    count = 0
    for obj in document[page - 1].get_objects(max_depth=15):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        n = raw.FPDFPath_CountSegments(obj.raw)
        if n < shortest:
            continue
        channels = [ctypes.c_uint() for _ in range(4)]
        if not raw.FPDFPageObj_GetStrokeColor(
            obj.raw, *(ctypes.byref(c) for c in channels)
        ):
            continue
        if tuple(c.value for c in channels[:3]) != colour:
            continue
        count += sum(
            raw.FPDFPathSegment_GetType(raw.FPDFPath_GetPathSegment(obj.raw, i))
            == raw.FPDF_SEGMENT_BEZIERTO
            for i in range(n)
        )
    return count


def _pieces(key: str, colour: tuple[int, int, int]) -> list[np.ndarray]:
    """Return a panel's stroked polylines of a colour (20 vertices or more), sorted."""
    found = frame(key)
    pieces = [
        d.vertices
        for d in found.drawn
        if d.stroke == colour
        and d.fill is None
        and len(d.vertices) >= 20
        and _inside(*d.centre, found.box)
    ]
    pieces.sort(key=lambda v: float(v[0, 0]))
    for before, after in itertools.pairwise(pieces):
        if math.dist(before[-1], after[0]) > STRAIGHT:
            msg = f"{PANELS[key].name}: the curve's pieces do not join"
            raise ValueError(msg)
    return pieces


def curve(key: str, colour: tuple[int, int, int] = BLUE) -> list[Reading]:
    """Return a panel's curve: every vertex in the frame, in order, repeats dropped."""
    found = frame(key)
    kept: list[tuple[float, float]] = []
    for piece in _pieces(key, colour):
        for x, y in piece:
            point = (float(x), float(y))
            if _inside(*point, found.box) and (not kept or point != kept[-1]):
                kept.append(point)
    if any(b[0] < a[0] for a, b in itertools.pairwise(kept)):
        msg = f"{PANELS[key].name}: the curve turns back"
        raise ValueError(msg)
    return _read(key, kept)


def dots(key: str, colour: tuple[int, int, int]) -> list[Reading]:
    """Return a dashed or dotted curve's vertices inside the frame, in order of x.

    Its dashes are two-vertex strokes 2.5 pt or more from the frame's left, right and
    bottom edges (the ticks are nearer).
    """
    found = frame(key)
    left, bottom, right, top = found.box
    points = {
        (float(x), float(y))
        for d in found.drawn
        if d.stroke == colour and len(d.vertices) == 2 and np.ptp(d.vertices[:, 0]) < 3
        for x, y in d.vertices
        if left + 2.5 < x < right - 2.5 and bottom + 2.5 < y < top
    }
    return _read(key, sorted(points))


def bsi_line(key: str) -> float:
    """Return the intensity (W/cm^2) of a panel's dashed vertical I_BSI line.

    The line is the one position with 10 or more vertical dashes (a tick is one).
    """
    found = frame(key)
    positions = [
        round(float(d.vertices[0, 0]), 3)
        for d in found.drawn
        if d.stroke == BLACK
        and len(d.vertices) == 2
        and np.ptp(d.vertices[:, 0]) < STRAIGHT
        and DOT < np.ptp(d.vertices[:, 1]) < TICK_LONGEST
    ]
    lines = {x for x in positions if positions.count(x) >= 10}
    if len(lines) != 1:
        msg = f"{PANELS[key].name}: {len(lines)} dashed vertical lines"
        raise ValueError(msg)
    return float(found.x(lines.pop())) * X_SCALE.get(key, 1.0)


def interpolate(points: Sequence[Reading], at: float, *, log: bool) -> float:
    """Return a curve's value at an abscissa, linearly (in log-log if ``log``)."""
    x = np.array([p.abscissa for p in points])
    y = np.array([p.value for p in points])
    if log:
        return float(np.exp(np.interp(math.log(at), np.log(x), np.log(y))))
    return float(np.interp(at, x, y))


def grey_curve(key: str) -> list[Reading]:
    """Return Fig. 5's or 6's grey TDSE curve (its pieces' vertices in the frame).

    Its pieces are the grey stroked paths larger than a marker.
    """
    found = frame(key)
    kept = {
        (float(x), float(y))
        for d in found.drawn
        if d.stroke == GREY and d.fill is None and d.size > MARKER_LARGEST
        for x, y in d.vertices
        if _inside(float(x), float(y), found.box)
    }
    return _read(key, sorted(kept))


def on_curve(key: str, least: float = 0.0) -> float:
    """Return the largest relative distance of a panel's markers from its grey curve.

    Where P >= ``least``, and in (b) only where P >= ON_CURVE_FROM (the linear axis
    does not resolve less).
    """
    log = PANELS[key].y_log
    line = grey_curve(key)
    worst = 0.0
    for m in markers(key, TDSE):
        if not line[0].abscissa <= m.abscissa <= line[-1].abscissa:
            continue
        if m.value < least or (not log and m.value < ON_CURVE_FROM):
            continue
        worst = max(worst, abs(interpolate(line, m.abscissa, log=log) / m.value - 1))
    return worst


def same_markers(copy: str, original: str) -> tuple[int, float, float]:
    """Return a linear panel's markers against another's: count, worst dI, dP.

    Each of ``copy``'s markers is paired with ``original``'s nearest in intensity.
    """
    reference = markers(original, SERIES[original])
    worst_i = worst_p = 0.0
    found = markers(copy, SERIES[copy])
    for m in found:
        r = min(reference, key=lambda r: abs(r.abscissa - m.abscissa))
        worst_i = max(worst_i, abs(m.abscissa - r.abscissa))
        worst_p = max(worst_p, abs(m.value - r.value))
    return len(found), worst_i, worst_p


def fig6a_against_fig2() -> list[tuple[float, float, float]]:
    """Return Fig. 6(a)'s markers with Fig. 2's: (intensity, P 6(a), P Fig. 2)."""
    reference = markers("fig2", IONISED)
    out = []
    for m in markers("fig6a", TDSE):
        r = min(reference, key=lambda r: abs(math.log(r.abscissa / m.abscissa)))
        if abs(r.abscissa - m.abscissa) > max(GRID * m.abscissa, LINEAR_INTENSITY):
            msg = f"Fig. 6(a): the marker at {m.abscissa:.4e} has none in Fig. 2"
            raise ValueError(msg)
        out.append((m.abscissa, m.value, r.value))
    return out


def fig5a_against_fig1() -> list[tuple[float, float, float]]:
    """Return Fig. 5(a)'s markers at Fig. 1's intensities: (intensity, 5(a), Fig. 1).

    Only where both have a marker within GRID.
    """
    reference = markers("fig1", IONISED)
    out = []
    for m in markers("fig5a", TDSE):
        r = min(reference, key=lambda r: abs(math.log(r.abscissa / m.abscissa)))
        if abs(r.abscissa / m.abscissa - 1) <= GRID:
            out.append((m.abscissa, m.value, r.value))
    return out


def probabilities(key: str) -> list[tuple[float, float, float, float]]:
    """Return Fig. 1's or 2's (intensity, ground, excitation, ionisation), in order."""
    ground = markers(key, GROUND)
    excited = markers(key, EXCITED)
    ionised = markers(key, IONISED)
    out = []
    for g, e, i in zip(ground, excited, ionised, strict=True):
        if max(abs(g.abscissa - i.abscissa), abs(e.abscissa - i.abscissa)) > 2e11:
            msg = f"{PANELS[key].name}: the series' intensities differ at {i.abscissa}"
            raise ValueError(msg)
        out.append((i.abscissa, g.value, e.value, i.value))
    return out


def humps(key: str, ranges: Sequence[tuple[float, float]]) -> list[float]:
    """Return the omega (au) of the excitation curve's highest point in each range."""
    excited = dots(key, RED)
    out = []
    for low, high in ranges:
        within = [r for r in excited if low <= r.abscissa <= high]
        out.append(max(within, key=lambda r: r.value).abscissa)
    return out


def excitation_above_ionisation(key: str) -> list[tuple[float, float]]:
    """Return the omega ranges (au) over which the excitation exceeds the ionisation."""
    ionised = curve(key)
    excited = dots(key, RED)
    flags = [
        (r.abscissa, r.value > interpolate(ionised, r.abscissa, log=True))
        for r in excited
        if ionised[0].abscissa <= r.abscissa <= ionised[-1].abscissa
    ]
    ranges, start = [], None
    for (omega, above), (_, after) in zip(
        flags, [*flags[1:], (0.0, False)], strict=True
    ):
        if above and start is None:
            start = omega
        if start is not None and not after:
            ranges.append((start, omega))
            start = None
    return ranges


def intensity(field_au: float) -> float:
    """Return the intensity (W/cm^2) of a peak field (au): I = F^2 in atomic units."""
    return field_au**2 * I_AU


# The text's (F au, I W/cm^2) pairs, pp. 1-6 of the journal; '2.5 F_BSI = 0.156 a.u.'
# and 'F_BSI = 0.063 a.u.' as their exact values.
TEXT_FIELDS = (
    (2.5 * F_BSI, 8.6e14),
    (F_BSI, 1.37e14),
    (0.053, 1.0e14),
    (0.107, 4.0e14),
    (0.038, 5e13),
    (0.027, 2.5e13),
)
HUMPS = {
    "fig7": ((0.375, 0.444), ((0.33, 0.41), (0.41, 0.5))),
    "fig8": ((0.114, 0.228), ((0.09, 0.14), (0.18, 0.27))),
}


def checks() -> None:  # noqa: PLR0912 - one check after another
    """Raise unless the readings agree with every number the paper prints."""
    for copy, original in (("fig5b", "fig1"), ("fig6b", "fig2")):
        n, d_i, d_p = same_markers(copy, original)
        if d_i > LINEAR_INTENSITY or d_p > LINEAR_AGREEMENT:
            msg = f"{copy}: {n} markers, worst dI {d_i:.3g}, dP {d_p:.5f}"
            raise ValueError(msg)
    for i, p6, p2 in fig6a_against_fig2():
        if p2 >= LOG_LINEAR_FROM and abs(p6 / p2 - 1) > LOG_LINEAR:
            msg = f"Fig. 6(a) at {i:.4e}: {p6:.4g}, Fig. 2 {p2:.4g}"
            raise ValueError(msg)
    for key in ("fig5a", "fig5b", "fig6a", "fig6b"):
        if abs(bsi_line(key) / 1.37e14 - 1) > BSI:
            msg = f"{key}: the I_BSI line at {bsi_line(key):.4e}"
            raise ValueError(msg)
    for key in ("fig1", "fig2"):
        rows = probabilities(key)
        for i, g, e, p in rows:
            if abs(g + e + p - 1) > SUM_TOLERANCE:
                msg = f"{key} at {i:.3e}: the probabilities sum to {g + e + p:.4f}"
                raise ValueError(msg)
        if rows[-1][3] < NEARLY_FULL:
            msg = f"{key}: {rows[-1][3]:.3f} ionised at {rows[-1][0]:.3e}"
            raise ValueError(msg)
    if max(e for _, _, e, _ in probabilities("fig1")) >= EXCITATION_BELOW:
        msg = "Fig. 1: the excitation reaches 10 %"
        raise ValueError(msg)
    i, _, e, _ = max(probabilities("fig2"), key=lambda r: r[2])
    low, high, first, last = EXCITATION_PEAK
    if not (low <= e <= high and first <= i <= last):
        msg = f"Fig. 2: the excitation peaks at {e:.3f}, {i:.3e} W/cm^2"
        raise ValueError(msg)
    for key, (expected, ranges) in HUMPS.items():
        for omega, text in zip(humps(key, ranges), expected, strict=True):
            if abs(omega / text - 1) > HUMP:
                msg = f"{key}: an excitation hump at {omega:.4f}, the text's {text}"
                raise ValueError(msg)
    # 'the excitation is significant in a frequency range of ~0.37 a.u. <= omega <=
    # ~0.47 a.u. ... For these frequencies the ionization is much weaker'.
    wide = max(excitation_above_ionisation("fig8"), key=lambda r: r[1] - r[0])
    if not wide[0] <= EXCITED_RANGE[0] < EXCITED_RANGE[1] <= wide[1]:
        msg = f"Fig. 8: the excitation exceeds the ionisation over {wide}"
        raise ValueError(msg)
    for field, printed in TEXT_FIELDS:
        if abs(math.sqrt(printed / I_AU) - field) > HALF_DIGIT:
            msg = f"F = {field} au is {intensity(field):.3e} W/cm^2, not {printed}"
            raise ValueError(msg)
    for key in FIELD:
        if _bezier_segments(PANELS[key].page, BLUE, 20):
            msg = f"{key}: the ionisation curve has Bezier segments"
            raise ValueError(msg)


def _points(key: str, found: Sequence[Reading]) -> list[Point]:
    """Return a dataset's points, numbered in order of the abscissa."""
    name = PANELS[key].name
    if key in FIELD:
        what, x_format, y_format = "TDSE ionisation curve: vertex", "{:.5g}", "{:.4g}"
    elif PANELS[key].y_log:
        what, x_format, y_format = "TDSE: marker", "{:.4e}", "{:.4g}"
    else:  # a linear axis: the probability to 0.0001
        what, x_format, y_format = "ionization: marker", "{:.4e}", "{:.4f}"
    width = len(str(len(found)))
    return [
        Point(
            f"p{k:0{width}d}",
            x_format.format(r.abscissa),
            y_format.format(r.value),
            "",
            f"{name} {what} at ({r.x:.3f}, {r.y:.3f}) pt",
        )
        for k, r in enumerate(found, start=1)
    ]


def readings(key: str) -> list[Reading]:
    """Return a dataset's readings: the panel's markers, or Figs. 7 and 8's curve."""
    return curve(key) if key in FIELD else markers(key, SERIES[key])


def extract() -> dict[str, list[Point]]:
    """Return the datasets' points, one file per panel."""
    checks()
    return {relative: _points(key, readings(key)) for key, relative in DATASETS.items()}


def pulse_fs(omega_au: float, cycles: int = 14) -> float:
    """Return the duration (fs) of a number of cycles at a frequency (au)."""
    return cycles * 2 * math.pi / omega_au * AU_TIME_FS


def report() -> None:
    """Print the readings' and the checks' numbers."""
    for key, panel in PANELS.items():
        found = frame(key)
        print(
            f"{panel.name}: frame {np.round(found.box, 2)}, "
            f"ticks {len(found.x.positions)} x {len(found.y.positions)}, "
            f"worst tick residual {residual(key):.3f} pt"
        )
    for key in DATASETS:
        r = readings(key)
        values = [p.value for p in r]
        print(
            f"{key}: {len(r)} points, {r[0].abscissa:.4g}-{r[-1].abscissa:.4g}, "
            f"P {min(values):.3g}-{max(values):.3g}"
        )
    for copy, original in (("fig5b", "fig1"), ("fig6b", "fig2")):
        print(
            copy,
            "against",
            original,
            "(n, dI W/cm^2, dP)",
            same_markers(copy, original),
        )
    pairs = fig6a_against_fig2()
    print(
        "Fig. 6(a)/Fig. 2 - 1 where P >= 0.03:",
        [round(a / b - 1, 4) for _, a, b in pairs if b >= LOG_LINEAR_FROM],
    )
    print(
        "Fig. 5(a) against Fig. 1 at the same intensities (I, 5(a), Fig. 1):",
        [(f"{i:.3e}", round(a, 4), round(b, 4)) for i, a, b in fig5a_against_fig1()],
    )
    for key in ("fig5a", "fig5b", "fig6a", "fig6b"):
        print(
            f"{key}: the grey curve off the markers <= {on_curve(key):.4f}, "
            f"<= {on_curve(key, 0.1):.4f} where P >= 0.1; "
            f"I_BSI line {bsi_line(key):.4e}"
        )
    for key in ("fig1", "fig2"):
        rows = probabilities(key)
        sums = [g + e + p for _, g, e, p in rows]
        i, _, e, _ = max(rows, key=lambda r: r[2])
        print(
            f"{key}: sums {min(sums):.4f}-{max(sums):.4f}; excitation max {e:.4f} at "
            f"{i:.3e}; ionisation at 4e14 {rows[-1][3]:.4f}; at 3.5e14 excitation "
            f"{rows[-2][2]:.4f}"
        )
    for key, (_, ranges) in HUMPS.items():
        print(key, "excitation humps", [round(h, 4) for h in humps(key, ranges)])
    print("Fig. 8 excitation > ionisation", excitation_above_ionisation("fig8"))
    for field, printed in TEXT_FIELDS:
        print(f"F = {field} au: {intensity(field):.4e} W/cm^2 (printed {printed:.3g})")
    for key, field in FIELD.items():
        line = curve(key)
        i = intensity(field)
        for omega, other in ((0.057, "fig5a"), (0.057, "fig1"), (0.114, "fig6a")):
            mark = markers(other, SERIES[other])
            log = PANELS[other].y_log
            print(
                f"{key} (F = {field}, I = {i:.4e}) at omega {omega}: "
                f"{interpolate(line, omega, log=True):.4g}; {other} at that "
                f"intensity {interpolate(mark, i, log=log):.4g}"
            )
    for omega in (0.057, 0.114):
        print(f"14 cycles at omega {omega}: {pulse_fs(omega):.3f} fs")


def main() -> None:
    """Write the data files and print the checks' numbers."""
    write(extract())
    report()


if __name__ == "__main__":
    main()
