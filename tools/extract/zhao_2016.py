# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Zhao et al. (2016): one-electron TDSE ionisation of H and the rare gases.

The paper prints no table of its calculations: the time-dependent Schrodinger equation
(TDSE) probabilities are given only in **Fig. 2** (p. 3, library figure
``fig_962af80d9355``), atomic hydrogen at 600, 800, 1000 and 1200 nm, and **Fig. 3**
(p. 3, ``fig_2ad80701eb1d``), He and Ne at 400 nm and Ar, Kr and Xe at 800 nm, all for
a 15 fs (FWHM) Gaussian pulse (the captions). The correction factors TDSE/PPT of **Fig.
1** (p. 2, ``fig_1c7094c066f6``) are read too, but only to check the printed factors and
for the record's notes: no quantity holds them. The text prints the Xe model potential
(p. 2), an inline dataset of the record, not read here.

All three figures are vector drawings on the page (no image, no form) and are read
from the drawing (data rule 8, TC-699) with :func:`_figure.painted`. Their labels are
glyphs drawn as filled paths, not text.

* the frames: four stroked black lines per panel (bottom, left, top, right), and the
  ticks separate stroked black two-vertex lines perpendicular to an edge with one end
  on it; long and short ticks (the long ones about twice as long, or more). The top and
  right edges repeat the bottom and left ones' ticks (checked: every tick there has
  one at the same position below or left, within :data:`MIRROR`).
* the axes (:data:`PANELS`): the labels, read on the page, sit at the long ticks and
  only there (checked, :func:`_check_labels`).

  - linear (Fig. 1, Fig. 2's insets): one short tick halfway between long ones; every
    tick is assigned its value counting in half steps from the first long tick, and
    each long tick must get its label's value.
  - logarithmic in decades (the y axes of Figs. 2 and 3): the long ticks the decades,
    eight short ones between, 2 to 9 times the lower decade.
  - logarithmic in steps (the x axes of Figs. 2 and 3): the ticks are consecutive
    multiples of one step (10^13 or 10^14 W/cm^2) also beyond the decade, e.g. Fig.
    2's 2, 3, ..., 14 x 10^13 W/cm^2 with the long tick at 10^14, and Fig. 3(c)'s
    4 ... 24 x 10^13 with long ticks at its labels 10^14 and 2 x 10^14. That this is
    what the ticks are is checked: every tick lies on the fitted log scale within
    :data:`TICK_TOLERANCE` (Fig. 2's 13 ticks within 0.07 pt on 186 pt per decade),
    which the multiples of a decade would miss by several points; and Fig. 2's
    insets, a linear axis from 6 to 13 x 10^13 W/cm^2 labelled at four long ticks,
    put the same TDSE points at the same intensities (:func:`checks`).

  Every tick lies on its fitted scale within :data:`TICK_TOLERANCE`
  (:meth:`_figure.Axis.check`): the drawing's coordinates are quantised in steps of
  about 0.07 pt (Fig. 2's y ticks alternate 14.10 and 14.20 pt apart for equal steps).
* the markers: Fig. 2's TDSE points are blue open circles (stroked 13-vertex paths,
  7.36 pt across in the main panels, 3.96 pt in the insets, the legends' smaller or
  larger); the insets' PPT points red open stars (11-vertex, 5.2 pt; the legends'
  6.3 pt). Fig. 3's TDSE points blue open circles 4.44 pt across, the legend's one
  directly above the legend's PPT line sample (a black horizontal line about 10 pt
  long); Fig. 1's averages red open circles 4.87 pt (the legends' 6.09 pt). Their
  centres. A few of Fig. 3's markers are centred just outside the frame (up to 2.2 pt,
  the axis extended), the rest inside it or on it.
* the lines: Fig. 1(a)'s correction factor is one black polyline (1428 vertices), its
  'Fitting with a constant' red dashes; Fig. 1(b)'s and (c)'s 'Polynomial fit' blue
  polylines (20 vertices, extending beyond the frame); Fig. 2's 'PPTX<factor>' black
  polylines, one per main panel and inset.

Every point's ``printed`` cell gives the panel and the marker's centre in page points.
The rows are numbered in order of intensity.

**Checks** (:func:`checks`), against every number the paper prints for these
quantities:

* Fig. 1(a)'s dashed line lies at 1.49, the label beside it and the caption's 'The
  average of the correction factors is taken to be 1.49'.
* Fig. 1(a)'s curve averages 1.49 over its intensities (trapezoids; 1.490).
* Fig. 1(b)'s averages at 600, 800, 1000 and 1200 nm are the factors Fig. 2's legends
  print, 'PPTX1.42', 'PPTX1.49', 'PPTX1.36' and 'PPTX1.26' (within
  :data:`FIG1_TOLERANCE`).
* in each inset of Fig. 2, the 'PPTX' line over the PPT stars is the legend's factor
  applied to the rate: the line is 1 - (1 - P)^factor of the stars' P, not factor x P
  (within :data:`FACTOR_TOLERANCE`, where P >= :data:`FACTOR_FROM`).
* Fig. 2's insets and main panels give the same TDSE probabilities at the same
  intensities (within :data:`INSET_AGREEMENT`).
* Fig. 2's main panels put their TDSE points at 2.0, 2.5, ..., 14.0 x 10^13 W/cm^2
  (within :data:`GRID`), and the insets at 6.0 ... 13.0 x 10^13.
* the arrows marking the Keldysh parameter in Fig. 2 ('gamma=2', ... read on the
  page) lie at the intensity where hydrogen's gamma = omega sqrt(2 I_p)/F (I_p = 0.5
  hartree) takes that value (within :data:`ARROW`); in Fig. 3 the same with each
  gas's ionisation energy (optmat's default model's first channel).

The text's statements are reported (:func:`fig1`): Fig. 1(a)'s factor 'oscillates but
is within about a factor of 2'; Fig. 1(b)'s smooth curve 'is between 1.1 and 1.4'; Fig.
1(c)'s 'covers from 1.4 to 1.1'.
"""

import functools
import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from _figure import Axis, Painted, painted
from _library import Point, write

from optmat.constants import AU_ANGULAR_FREQUENCY, AU_FIELD, EPS0, C
from optmat.gas.gases import Gas
from optmat.gas.ionisation_rates import ionisation_rate_models

DIRECTORY = "Zhao_2016_AnalyticalModelCalibrating"
SOURCE = "zhao_2016"
BLACK, RED, BLUE = (0, 0, 0), (255, 0, 0), (0, 0, 255)
# The drawing's coordinates are quantised in steps of about 0.0708 pt: ticks of equal
# steps are drawn up to 0.1 pt apart. Every tick lies on its fitted scale within 0.07 pt
# (the largest residual, Fig. 2's x axis); 0.1 pt allows for it and no more.
TICK_TOLERANCE = 0.1
STRAIGHT = 0.01  # a tick or an edge is vertical or horizontal to 0.01 pt
# A tick ends on its edge within 0.03 pt (the edges' ends and the ticks are quantised
# separately: up to 0.011 pt apart).
ON_LINE = 0.03
TICK_LONGEST = 4.0  # ticks are under 3 pt long; the frames' edges 50 pt or more
LONG_TICK = 1.5  # a long tick is over 1.5 times the shortest (1.75 to 1.9 times)
MIRROR = 0.03  # the top and right edges' ticks are the bottom and left ones' to 0.02 pt
EDGE = 40.0  # a frame's edge is a stroked line longer than this (pt)
# A marker is in the frame if its centre is within 0.05 pt of it (Figs. 1 and 2 put
# the first and last markers on the edges, centred up to 0.03 pt outside); Fig. 3 draws
# some just beyond the axis's end, up to 2.2 pt (its markers are 4.44 pt across).
ON_EDGE = 0.05
BEYOND_EDGE = 2.5
# The labels' glyphs: centred within a panel's label_left of the frame (y; the axis
# title is further: 17 pt or more in Fig. 3, 22 in Fig. 1, 28 in Fig. 2) or 12 pt below
# it (x); a label's glyphs are less than 1.5 pt apart (its neighbours 2 pt or more), and
# its middle lies within 3 pt of its long tick (a superscript shifts it). In Figs. 1 to
# 3's main panels a label's glyphs are up to 3.5 pt apart (the 'x10' and its
# superscript), its neighbours 10 pt or more: label_gap 4 pt there, 1.5 in the insets.
LABEL_BELOW = 12.0
LABEL_GAP = 1.5
LABEL_AT = 3.0
# Fig. 1(b)'s averages at 600 ... 1200 nm are the legends' printed factors within
# 0.002: its y scale is 164.9 pt per unit, the drawing's quantum 0.07 pt (0.0004); they
# agree within 0.0005. The insets' factors (solved from the PPTX line and the stars)
# within 0.02 where P >= 0.02: there the 0.07 pt quantum is up to 1.2 % of P (Fig.
# 2(a)'s inset, 303 pt per unit), about 0.017 in the factor; below, it is not resolved.
FIG1_TOLERANCE = 0.002
FACTOR_TOLERANCE = 0.02
FACTOR_FROM = 0.02
# Inset and main panel agree within 3 %, or 0.0007 in probability where the inset's
# quantum (0.07 pt is 0.0003 to 0.0006 in probability) dominates.
INSET_AGREEMENT = 0.03
INSET_ABSOLUTE = 7e-4
GRID = 0.01  # the TDSE intensities are on the 0.5e13 grid within 1 %
INSET_RANGE = (5.9e13, 1.31e14)  # the insets' points, 6.0 to 13.0e13 W/cm^2
# The arrows lie within 5 % of the intensity of the printed gamma (the largest
# difference found 1.3 %, Fig. 3(e)'s gamma = 1.8; the arrows are placed by hand, and
# the model atoms' I_p are not printed).
ARROW = 0.05
I_AU = 0.5 * C * EPS0 * AU_FIELD**2 / 1e4  # W/cm^2 per atomic unit of intensity
H_IP_AU = 0.5  # hydrogen's ionisation energy in a one-electron TDSE (infinite mass)


@dataclass(frozen=True)
class Scale:
    """How an axis's ticks get their values.

    ``kind``: ``linear`` (``labels``: the long ticks' values in order; one short tick
    halfway between), ``decades`` (``labels``: the long ticks' decades, as powers of
    ten) or ``steps`` (``labels``: the long ticks' values; every tick a consecutive
    multiple of ``step``).
    """

    kind: str
    labels: tuple[float, ...]
    step: float = 0.0


@dataclass(frozen=True)
class Panel:
    """A panel: name, page, frame's lower left corner (pt, to 1 pt), axes.

    ``label_left``: how far left of the frame (pt) the y labels' glyphs are centred;
    ``label_gap``: the largest gap (pt) between the glyphs of one label.
    """

    name: str
    page: int
    corner: tuple[float, float]
    x: Scale
    y: Scale
    label_left: float
    label_gap: float = 1.5


FIG2_X = Scale("steps", (1e14,), 1e13)
FIG2_Y = Scale("decades", tuple(range(-9, 1)))
FIG3_Y = Scale("decades", tuple(range(-7, 1)))
INSET_X = Scale("linear", (6e13, 8e13, 1e14, 1.2e14))
PANELS = {
    # Fig. 1: x '6.0x10^13' ... '1.2x10^14' (W/cm^2), '600' ... '1600' (nm), '8' ...
    # '16' (eV); y '1.0' ... '2.4', '0.9' ... '1.5', '0.8' ... '1.8'.
    "fig1a": Panel(
        "Fig. 1(a)",
        2,
        (366.6, 635.6),
        Scale("linear", (6e13, 8e13, 1e14, 1.2e14)),
        Scale("linear", (1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4)),
        16.0,
        4.0,
    ),
    "fig1b": Panel(
        "Fig. 1(b)",
        2,
        (366.6, 509.0),
        Scale("linear", (600, 800, 1000, 1200, 1400, 1600)),
        Scale("linear", (0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)),
        16.0,
        4.0,
    ),
    "fig1c": Panel(
        "Fig. 1(c)",
        2,
        (366.6, 382.4),
        Scale("linear", tuple(range(8, 17))),
        Scale("linear", (0.8, 1.0, 1.2, 1.4, 1.6, 1.8)),
        16.0,
        4.0,
    ),
    # Fig. 2: x '10^14' (W/cm^2) under one long tick; y '10^-9' ... '10^0'.
    "fig2a": Panel("Fig. 2(a)", 3, (148.1, 611.2), FIG2_X, FIG2_Y, 20.0, 4.0),
    "fig2b": Panel("Fig. 2(b)", 3, (337.1, 611.2), FIG2_X, FIG2_Y, 20.0, 4.0),
    "fig2c": Panel("Fig. 2(c)", 3, (148.1, 472.6), FIG2_X, FIG2_Y, 20.0, 4.0),
    "fig2d": Panel("Fig. 2(d)", 3, (337.1, 472.6), FIG2_X, FIG2_Y, 20.0, 4.0),
    # The insets: x '6.0x10^13', '8.0x10^13', '1.0x10^14', '1.2x10^14'; y '0.00' ...
    "fig2a_inset": Panel(
        "Fig. 2(a) inset",
        3,
        (226.8, 621.0),
        INSET_X,
        Scale("linear", (0.0, 0.04, 0.08, 0.12, 0.16)),
        12.0,
    ),
    "fig2b_inset": Panel(
        "Fig. 2(b) inset",
        3,
        (419.0, 621.0),
        INSET_X,
        Scale("linear", (0.0, 0.02, 0.04, 0.06, 0.08, 0.10)),
        12.0,
    ),
    "fig2c_inset": Panel(
        "Fig. 2(c) inset",
        3,
        (231.6, 482.2),
        INSET_X,
        Scale("linear", (0.0, 0.02, 0.04, 0.06, 0.08)),
        12.0,
    ),
    "fig2d_inset": Panel(
        "Fig. 2(d) inset",
        3,
        (419.0, 482.2),
        INSET_X,
        Scale("linear", (0.0, 0.02, 0.04, 0.06)),
        12.0,
    ),
    # Fig. 3: x (a) '10^15', (b) '8x10^14', (c) '10^14' and '2x10^14', (d) '10^14',
    # (e) '8x10^13'; y '10^-7' ... '10^0'.
    "fig3a": Panel(
        "Fig. 3(a)", 3, (338.8, 313.3), Scale("steps", (1e15,), 1e14), FIG3_Y, 14.0, 4.0
    ),
    "fig3b": Panel(
        "Fig. 3(b)", 3, (456.5, 313.3), Scale("steps", (8e14,), 1e14), FIG3_Y, 14.0, 4.0
    ),
    "fig3c": Panel(
        "Fig. 3(c)",
        3,
        (338.8, 224.5),
        Scale("steps", (1e14, 2e14), 1e13),
        FIG3_Y,
        14.0,
        4.0,
    ),
    "fig3d": Panel(
        "Fig. 3(d)", 3, (456.5, 224.5), Scale("steps", (1e14,), 1e13), FIG3_Y, 14.0, 4.0
    ),
    "fig3e": Panel(
        "Fig. 3(e)", 3, (338.8, 138.6), Scale("steps", (8e13,), 1e13), FIG3_Y, 14.0, 4.0
    ),
}
# The insets' frames lie inside the main panels'; their markers are not the main ones.
INSET_OF = {f"fig2{p}": f"fig2{p}_inset" for p in "abcd"}


@dataclass(frozen=True)
class Marker:
    """A marker series: stroke colour, vertex count and the data markers' size (pt)."""

    stroke: tuple[int, int, int]
    vertices: int
    size: tuple[float, float]


TDSE_MAIN = Marker(BLUE, 13, (7.3, 7.45))
TDSE_INSET = Marker(BLUE, 13, (3.9, 4.05))
PPT_INSET = Marker(RED, 11, (5.15, 5.3))
TDSE_FIG3 = Marker(BLUE, 13, (4.4, 4.5))
AVERAGE_FIG1 = Marker(RED, 13, (4.8, 4.95))

# Fig. 2: each panel's wavelength (nm) and its legend's 'PPTX' factor, read on the page.
FIG2 = {"a": (600, 1.42), "b": (800, 1.49), "c": (1000, 1.36), "d": (1200, 1.26)}
# Fig. 3: gas and wavelength (nm) of each panel, the caption's.
FIG3 = {
    "a": ("He", 400),
    "b": ("Ne", 400),
    "c": ("Ar", 800),
    "d": ("Kr", 800),
    "e": ("Xe", 800),
}
# The Keldysh-parameter arrows' labels, read on the page, by panel, in order of x.
ARROWS = {
    "fig2a": (2.0, 1.3),
    "fig2b": (1.5, 1.0),
    "fig2c": (1.5, 1.0),
    "fig2d": (1.3, 1.0),
    "fig3a": (1.5, 1.0),
    "fig3b": (1.5, 1.0),
    "fig3c": (1.5, 1.0),
    "fig3d": (1.5, 1.0),
    "fig3e": (1.8, 1.2),
}
DATASETS = {
    **{
        f"fig2{p}": f"H/{SOURCE}__fig2{p}_tdse_H_{w}nm.csv"
        for p, (w, _) in FIG2.items()
    },
    **{
        f"fig3{p}": f"{g}/{SOURCE}__fig3{p}_tdse_{g}_{w}nm.csv"
        for p, (g, w) in FIG3.items()
    },
}


@dataclass(frozen=True)
class Frame:
    """A panel's frame (l, b, r, t, pt), its axes, and the page's paths around it."""

    box: tuple[float, float, float, float]
    x: Axis
    y: Axis
    x_long: tuple[float, ...]
    y_long: tuple[float, ...]
    drawn: tuple[Painted, ...]


def _edges(drawn: list[Painted]) -> list[np.ndarray]:
    """Return the frames' edges: stroked black straight lines longer than EDGE."""
    out = []
    for d in drawn:
        v = d.vertices
        if d.stroke != BLACK or d.fill is not None or not 2 <= len(v) <= 3:
            continue
        extent = v.max(axis=0) - v.min(axis=0)
        if max(extent) > EDGE and min(extent) < STRAIGHT:
            out.append(v)
    return out


def _box(
    drawn: list[Painted], corner: tuple[float, float], where: str
) -> tuple[float, float, float, float]:
    """Return the frame (l, b, r, t) whose lower left corner is within 1 pt of this."""
    edges = _edges(drawn)
    x0, y0 = corner
    bottoms = [
        e
        for e in edges
        if np.ptp(e[:, 1]) < STRAIGHT
        and abs(e[0, 1] - y0) < 1
        and abs(e[:, 0].min() - x0) < 1
    ]
    lefts = [
        e
        for e in edges
        if np.ptp(e[:, 0]) < STRAIGHT
        and abs(e[0, 0] - x0) < 1
        and abs(e[:, 1].min() - y0) < 1
    ]
    if len(bottoms) != 1 or len(lefts) != 1:
        msg = f"{where}: {len(bottoms)} bottom and {len(lefts)} left edges"
        raise ValueError(msg)
    left, right = float(bottoms[0][:, 0].min()), float(bottoms[0][:, 0].max())
    bottom, top = float(lefts[0][:, 1].min()), float(lefts[0][:, 1].max())
    for fixed, value, along in ((0, right, (bottom, top)), (1, top, (left, right))):
        found = [
            e
            for e in edges
            if np.ptp(e[:, fixed]) < STRAIGHT
            and abs(e[0, fixed] - value) < STRAIGHT * 5
            and abs(e[:, 1 - fixed].min() - along[0]) < 0.05
            and abs(e[:, 1 - fixed].max() - along[1]) < 0.05
        ]
        if len(found) != 1:
            msg = f"{where}: the {'right' if fixed == 0 else 'top'} edge is missing"
            raise ValueError(msg)
    return left, bottom, right, top


def _edge_ticks(
    drawn: list[Painted], box: tuple[float, float, float, float], axis: int, far: bool
) -> list[tuple[float, float]]:
    """Return an edge's ticks as (position, length), sorted.

    ``axis`` 0: the bottom (or, ``far``, top) edge's vertical ticks; 1: the left (or
    right) edge's horizontal ticks. Only ticks within the edge's span.
    """
    along, across = (0, 1) if axis == 0 else (1, 0)
    edge = box[across + 2] if far else box[across]
    low, high = box[along], box[along + 2]
    ticks = set()
    for d in drawn:
        v = d.vertices
        if (
            d.stroke != BLACK
            or len(v) != 2
            or abs(v[0, along] - v[1, along]) > STRAIGHT
        ):
            continue
        ends = sorted(v[:, across], key=lambda e: abs(e - edge))
        length = abs(ends[1] - ends[0])
        position = float(v[0, along])
        if (
            abs(ends[0] - edge) < ON_LINE
            and 0 < length < TICK_LONGEST
            and low - ON_LINE <= position <= high + ON_LINE
        ):
            ticks.add((round(position, 3), round(float(length), 3)))
    return sorted(ticks)


def _ticks(
    drawn: list[Painted], box: tuple[float, float, float, float], axis: int, where: str
) -> tuple[list[float], list[float]]:
    """Return an axis's tick positions and which of them are long.

    The far edge's ticks (top, right) must each have one on the near edge at the same
    position, within MIRROR.
    """
    near = _edge_ticks(drawn, box, axis, far=False)
    far = _edge_ticks(drawn, box, axis, far=True)
    positions = [p for p, _ in near]
    if any(min(abs(p - q) for q in positions) > MIRROR for p, _ in far):
        msg = f"{where}: a tick of the opposite edge has none on this one"
        raise ValueError(msg)
    shortest = min(length for _, length in near)
    long = [p for p, length in near if length > LONG_TICK * shortest]
    return positions, long


def _values(
    scale: Scale, positions: list[float], long: list[float], where: str
) -> list[float]:
    """Return each tick's value, from its scale and the long ticks' labels."""
    if len(long) != len(scale.labels):
        msg = f"{where}: {len(long)} long ticks, {len(scale.labels)} labels"
        raise ValueError(msg)
    first = positions.index(long[0])
    if scale.kind == "linear":
        half = (scale.labels[1] - scale.labels[0]) / 2
        values = [scale.labels[0] + (k - first) * half for k in range(len(positions))]
    elif scale.kind == "steps":
        n0 = round(scale.labels[0] / scale.step)
        values = [(n0 + k - first) * scale.step for k in range(len(positions))]
    else:  # decades: the long ticks, then 2 to 9 times the decade below each
        if positions[0] != long[0]:
            msg = f"{where}: a short tick below the lowest decade"
            raise ValueError(msg)
        values, decade, multiple = [], 0.0, 1
        for p in positions:
            if p in long:
                decade, multiple = scale.labels[long.index(p)], 1
            else:
                multiple += 1
            values.append(multiple * 10.0**decade)
    for p, label in zip(long, scale.labels, strict=True):
        got = values[positions.index(p)]
        expected = 10.0**label if scale.kind == "decades" else label
        if not math.isclose(got, expected, rel_tol=1e-9, abs_tol=1e-12):
            msg = (
                f"{where}: the long tick at {p} pt counts to {got}, labelled {expected}"
            )
            raise ValueError(msg)
    return values


def _axis(scale: Scale, positions: list[float], long: list[float], where: str) -> Axis:
    """Return the calibrated axis, checked against every tick."""
    values = _values(scale, positions, long, where)
    axis = Axis(
        tuple(positions),
        tuple(float(v) for v in values),
        log=scale.kind != "linear",
        tolerance=TICK_TOLERANCE,
    )
    axis.check()
    return axis


def _groups(spans: list[tuple[float, float]], gap: float) -> list[tuple[float, float]]:
    """Return the runs of glyph extents (low, high) less than ``gap`` apart."""
    groups: list[list[float]] = []
    for low, high in sorted(spans):
        if groups and low - groups[-1][1] < gap:
            groups[-1][1] = max(groups[-1][1], high)
        else:
            groups.append([low, high])
    return [(low, high) for low, high in groups]


def _check_labels(found: Frame, panel: Panel, inner: tuple | None) -> None:
    """Raise unless the label glyphs sit at the long ticks, one label each.

    The glyphs are filled black paths: the y labels' wholly left of the frame and
    centred within the panel's ``label_left`` of it, the x labels' wholly below it and
    centred within LABEL_BELOW, not a y label's; none over an inset's frame or its
    labels. A label is a
    run of glyphs less than the panel's ``label_gap`` apart, and its middle lies within
    LABEL_AT of its long tick. Where the x labels touch (Fig. 2(c) and (d)'s insets),
    they make fewer runs than there are long ticks: then the glyphs' extent must
    cover the first and last long ticks and its middle lie within LABEL_AT of theirs.
    """
    left, bottom, right, top = found.box
    glyphs = [d for d in found.drawn if d.fill == BLACK and d.stroke is None]

    def outside(x: float, y: float) -> bool:
        return inner is None or not (
            inner[0] - LABEL_BELOW < x < inner[2]
            and inner[1] - LABEL_BELOW < y < inner[3]
        )

    y_glyphs = [
        d
        for d in glyphs
        if left - panel.label_left < d.centre[0]
        and d.right < left
        and bottom - 3 < d.centre[1] < top + 3
        and outside(*d.centre)
    ]
    beside = [(d.bottom, d.top) for d in y_glyphs]
    under = [
        (d.left, d.right)
        for d in glyphs
        if bottom - LABEL_BELOW < d.centre[1]
        and d.top < bottom
        and left - 15 < d.centre[0] < right + 15
        and outside(*d.centre)
        and all(d is not g for g in y_glyphs)
    ]
    for name, spans, ticks in (("y", beside, found.y_long), ("x", under, found.x_long)):
        runs = _groups(spans, panel.label_gap)
        middles = [(low + high) / 2 for low, high in runs]
        if len(runs) == len(ticks) and all(
            abs(m - t) <= LABEL_AT for m, t in zip(middles, ticks, strict=True)
        ):
            continue
        low, high = min(s[0] for s in spans), max(s[1] for s in spans)
        touching = (
            name == "x"
            and 1 <= len(runs) < len(ticks)
            and low < ticks[0] < ticks[-1] < high
            and abs((low + high) / 2 - (ticks[0] + ticks[-1]) / 2) <= LABEL_AT
        )
        if not touching:
            msg = (
                f"{panel.name}: the {name} labels at {middles} are not at the long "
                f"ticks {ticks}"
            )
            raise ValueError(msg)


@functools.cache
def frame(key: str) -> Frame:
    """Return a panel's frame and calibrated axes, its labels checked."""
    panel = PANELS[key]
    x0, y0 = panel.corner
    region = (x0 - 40, y0 - 30, x0 + 200, y0 + 160)
    drawn = painted(DIRECTORY, panel.page, region)
    box = _box(drawn, panel.corner, panel.name)
    x_positions, x_long = _ticks(drawn, box, 0, f"{panel.name} x")
    y_positions, y_long = _ticks(drawn, box, 1, f"{panel.name} y")
    found = Frame(
        box,
        _axis(panel.x, x_positions, x_long, f"{panel.name} x"),
        _axis(panel.y, y_positions, y_long, f"{panel.name} y"),
        tuple(x_long),
        tuple(y_long),
        tuple(drawn),
    )
    inner = frame(INSET_OF[key]).box if key in INSET_OF else None
    _check_labels(found, panel, inner)
    return found


def _inside(x: float, y: float, box: tuple, margin: float) -> bool:
    """Return whether a point lies inside a box (l, b, r, t) widened by a margin."""
    return (
        box[0] - margin <= x <= box[2] + margin
        and box[1] - margin <= y <= box[3] + margin
    )


@dataclass(frozen=True)
class Reading:
    """One marker or vertex as read: its position (pt) and its values."""

    x: float
    y: float
    abscissa: float
    value: float


def _read(found: Frame, points: list[tuple[float, float]]) -> list[Reading]:
    """Return points read on a frame's axes, in order of x."""
    return [
        Reading(x, y, float(found.x(x)), float(found.y(y))) for x, y in sorted(points)
    ]


def _markers(
    key: str,
    series: Marker,
    margin: float = ON_EDGE,
    legend: Callable[[Painted], bool] | None = None,
) -> list[Reading]:
    """Return a panel's markers of a series, outside any inset, the legend's dropped.

    Markers of the series' stroke and vertex count but another size are the legend's
    (one); or, with ``legend``, the one marker it picks.
    """
    found = frame(key)
    inset = frame(INSET_OF[key]).box if key in INSET_OF else None
    candidates = [
        d
        for d in found.drawn
        if d.stroke == series.stroke
        and len(d.vertices) == series.vertices
        and _inside(*d.centre, found.box, margin)
        and not (inset is not None and _inside(*d.centre, inset, 0.5))
    ]
    if legend is None:
        data = [d for d in candidates if series.size[0] <= d.size <= series.size[1]]
        others = [d for d in candidates if d not in data]
    else:
        others = [d for d in candidates if legend(d)]
        data = [d for d in candidates if not legend(d)]
        if any(not series.size[0] <= d.size <= series.size[1] for d in data):
            msg = f"{PANELS[key].name}: a marker of another size"
            raise ValueError(msg)
    if len(others) != 1:
        msg = f"{PANELS[key].name}: {len(others)} legend markers, expected one"
        raise ValueError(msg)
    centres = {(round(d.centre[0], 3), round(d.centre[1], 3)) for d in data}
    if len(centres) != len(data):
        msg = f"{PANELS[key].name}: two markers share a centre"
        raise ValueError(msg)
    return _read(found, sorted(centres))


def _fig3_legend(key: str) -> Callable[[Painted], bool]:
    """Return the test for Fig. 3's legend marker: above the PPT line's sample."""
    found = frame(key)
    samples = [
        d
        for d in found.drawn
        if d.stroke == BLACK
        and len(d.vertices) == 2
        and np.ptp(d.vertices[:, 1]) < STRAIGHT
        and 9 < np.ptp(d.vertices[:, 0]) < 11
        and _inside(*d.centre, found.box, 0)
    ]
    if len(samples) != 1:
        msg = f"{PANELS[key].name}: {len(samples)} legend line samples"
        raise ValueError(msg)
    x, y = samples[0].centre

    def is_legend(d: Painted) -> bool:
        return abs(d.centre[0] - x) < 1 and 0 < d.centre[1] - y < 8

    return is_legend


def tdse(key: str) -> list[Reading]:
    """Return a main panel's TDSE points (Fig. 2 or 3), in order of intensity."""
    if key.startswith("fig2"):
        return _markers(key, TDSE_MAIN)
    return _markers(key, TDSE_FIG3, BEYOND_EDGE, _fig3_legend(key))


def _polyline(
    key: str, colour: tuple[int, int, int], vertices: tuple[int, int]
) -> tuple[Frame, list[tuple[float, float]]]:
    """Return the vertices of a panel's one polyline of a colour and vertex count."""
    found = frame(key)
    inset = frame(INSET_OF[key]).box if key in INSET_OF else None
    lines = [
        d
        for d in found.drawn
        if d.stroke == colour
        and d.fill is None
        and vertices[0] <= len(d.vertices) <= vertices[1]
        and _inside(*d.centre, found.box, 0)
        and not (inset is not None and _inside(*d.centre, inset, 0))
    ]
    if len(lines) != 1:
        msg = f"{PANELS[key].name}: {len(lines)} lines, expected one"
        raise ValueError(msg)
    v = lines[0].vertices
    left, _, right, _ = found.box
    kept = [
        (float(x), float(y)) for x, y in v if left - ON_EDGE <= x <= right + ON_EDGE
    ]
    return found, kept


def ppt_factor_line(key: str) -> list[Reading]:
    """Return a Fig. 2 panel's or inset's 'PPTX<factor>' line's vertices, read."""
    found, kept = _polyline(key, BLACK, (15, 30))
    return _read(found, kept)


Pairs = list[tuple[float, float]]


@dataclass(frozen=True)
class Fig1:
    """Fig. 1's readings.

    ``curve``: Fig. 1(a)'s vertices (intensity W/cm^2, factor); ``constant``: its
    dashed line's value; ``mean``: the curve's mean over intensity (trapezoids);
    ``wavelength`` and ``ionisation``: Fig. 1(b)'s and (c)'s averages (nm or eV,
    factor); ``wavelength_fit`` and ``ionisation_fit``: the blue polynomials' vertices
    within the frame.
    """

    curve: Pairs
    constant: float
    mean: float
    wavelength: Pairs
    wavelength_fit: Pairs
    ionisation: Pairs
    ionisation_fit: Pairs


def fig1() -> Fig1:
    """Return Fig. 1's readings: the curve, the averages and the fits."""
    found, kept = _polyline("fig1a", BLACK, (1000, 2000))
    curve = _read(found, kept)
    x = np.array([r.abscissa for r in curve])
    y = np.array([r.value for r in curve])
    mean = float(np.sum((y[1:] + y[:-1]) / 2 * np.diff(x)) / (x[-1] - x[0]))
    dashes = [
        d
        for d in found.drawn
        if d.stroke == RED
        and len(d.vertices) == 2
        and np.ptp(d.vertices[:, 1]) < STRAIGHT
        and _inside(*d.centre, found.box, 0)
    ]
    levels = sorted({round(float(d.vertices[0, 1]), 3) for d in dashes})
    # The legend's two dashes are higher, inside the legend at the top.
    constant_level = min(levels)
    if sum(round(float(d.vertices[0, 1]), 3) == constant_level for d in dashes) < 15:
        msg = "Fig. 1(a): the dashed constant is not found"
        raise ValueError(msg)
    series = {}
    for key, name in (("fig1b", "wavelength"), ("fig1c", "ionisation")):
        series[name] = [(r.abscissa, r.value) for r in _markers(key, AVERAGE_FIG1)]
        fit_frame, fit = _polyline(key, BLUE, (15, 30))
        series[f"{name}_fit"] = [(r.abscissa, r.value) for r in _read(fit_frame, fit)]
    return Fig1(
        [(r.abscissa, r.value) for r in curve],
        float(found.y(constant_level)),
        mean,
        series["wavelength"],
        series["wavelength_fit"],
        series["ionisation"],
        series["ionisation_fit"],
    )


def _arrows(key: str) -> list[float]:
    """Return the x (pt) of a panel's Keldysh-parameter arrows, in order."""
    found = frame(key)
    inset = frame(INSET_OF[key]).box if key in INSET_OF else None
    heads = {
        round(float(d.vertices[0, 0]), 3)
        for d in found.drawn
        if d.fill == BLACK
        and len(d.vertices) == 4
        and d.size < 5
        and _inside(*d.centre, found.box, 0)
        and not (inset is not None and _inside(*d.centre, inset, 0))
    }
    return sorted(heads)


def keldysh_intensity(gamma: float, wavelength_nm: float, ip_au: float) -> float:
    """Return the intensity (W/cm^2) at which gamma = omega sqrt(2 I_p)/F."""
    omega = 2 * math.pi * C / (wavelength_nm * 1e-9) / AU_ANGULAR_FREQUENCY
    field = omega * math.sqrt(2 * ip_au) / gamma
    return field**2 * I_AU


def arrows() -> dict[str, list[tuple[float, float, float]]]:
    """Return each arrow's (printed gamma, intensity read, intensity of that gamma).

    Fig. 2 with hydrogen's I_p; Fig. 3 with each gas's in optmat's default model.
    """
    models = ionisation_rate_models()["default"]
    out = {}
    for key, gammas in ARROWS.items():
        letter = key[-1]
        if key.startswith("fig2"):
            wavelength, ip = FIG2[letter][0], H_IP_AU
        else:
            gas, wavelength = FIG3[letter]
            ip = models[Gas(gas)].channels[0].ionization_energy_au
        found = frame(key)
        xs = _arrows(key)
        if len(xs) != len(gammas):
            msg = f"{PANELS[key].name}: {len(xs)} arrows, {len(gammas)} labels"
            raise ValueError(msg)
        out[key] = [
            (g, float(found.x(x)), keldysh_intensity(g, wavelength, ip))
            for g, x in zip(gammas, xs, strict=True)
        ]
    return out


def inset_comparison(letter: str) -> list[tuple[float, float, float, float]]:
    """Return a Fig. 2 inset's points beside the main panel's at the same intensity.

    As (inset intensity, main intensity, inset probability, main probability).
    """
    main = tdse(f"fig2{letter}")
    out = []
    for r in _markers(f"fig2{letter}_inset", TDSE_INSET):
        nearest = min(main, key=lambda m: abs(math.log(m.abscissa / r.abscissa)))
        out.append((r.abscissa, nearest.abscissa, r.value, nearest.value))
    return out


def inset_factors(letter: str) -> list[tuple[float, float, float]]:
    """Return a Fig. 2 inset's (intensity, PPT, factor) at each PPT star.

    The factor multiplies the rate, not the probability: the 'PPTX' line is
    1 - (1 - P)^factor of the stars' P (it is not factor x P, which it falls below by
    2 % at P = 0.07). The factor is solved from the line, interpolated linearly at the
    star's x: ln(1 - line)/ln(1 - P).
    """
    line = ppt_factor_line(f"fig2{letter}_inset")
    xs = np.array([r.x for r in line])
    ys = np.array([r.value for r in line])
    out = []
    for star in _markers(f"fig2{letter}_inset", PPT_INSET):
        if xs[0] <= star.x <= xs[-1]:
            above = float(np.interp(star.x, xs, ys))
            factor = math.log1p(-above) / math.log1p(-star.value)
            out.append((star.abscissa, star.value, factor))
    return out


def tdse_over_ppt(letter: str) -> list[tuple[float, float, float]]:
    """Return Fig. 2's TDSE over the PPT, as (intensity, ratio of P, ratio of rates).

    At each TDSE point of a main panel: the PPT is the 'PPTX' line (interpolated
    linearly in log-log at the point's intensity) with the factor taken out of the
    rate, P = 1 - (1 - P_line)^(1/factor); the ratio of rates is
    ln(1 - P_TDSE)/ln(1 - P). Reported, for the record's notes; the paper prints Fig.
    1(a)'s factor at 800 nm only.
    """
    factor = FIG2[letter][1]
    line = ppt_factor_line(f"fig2{letter}")
    xs = np.log([r.abscissa for r in line])
    ys = np.log([r.value for r in line])
    out = []
    for r in tdse(f"fig2{letter}"):
        if not xs[0] <= math.log(r.abscissa) <= xs[-1]:
            continue
        p_line = math.exp(float(np.interp(math.log(r.abscissa), xs, ys)))
        ppt = -math.expm1(math.log1p(-p_line) / factor)
        out.append((r.abscissa, r.value / ppt, math.log1p(-r.value) / math.log1p(-ppt)))
    return out


def fit_ends() -> dict[str, tuple[float, float]]:
    """Return Fig. 1(b)'s and (c)'s polynomial fits at the axes' ends.

    At 600 and 1600 nm, and 8 and 16 eV: interpolated linearly between the polyline's
    vertices, those beyond the frame included.
    """
    out = {}
    for key, name in (("fig1b", "wavelength"), ("fig1c", "ionisation")):
        found = frame(key)
        lines = [
            d
            for d in found.drawn
            if d.stroke == BLUE
            and 15 <= len(d.vertices) <= 30
            and _inside(*d.centre, found.box, 0)
        ]
        if len(lines) != 1:
            msg = f"{PANELS[key].name}: {len(lines)} polynomial fits"
            raise ValueError(msg)
        v = lines[0].vertices
        xs = np.asarray(found.x(v[:, 0]), dtype=float)
        ys = np.asarray(found.y(v[:, 1]), dtype=float)
        left, _, right, _ = found.box
        low, high = float(found.x(left)), float(found.x(right))
        out[name] = (float(np.interp(low, xs, ys)), float(np.interp(high, xs, ys)))
    return out


def checks() -> None:  # noqa: PLR0912 - one check after another
    """Raise unless the readings agree with every number the paper prints."""
    one = fig1()
    if abs(one.constant - 1.49) > 0.005:
        msg = f"Fig. 1(a): the dashed constant reads {one.constant}, not 1.49"
        raise ValueError(msg)
    averages = dict(one.wavelength)
    for wavelength, factor in FIG2.values():
        nearest = min(averages, key=lambda w: abs(w - wavelength))
        if (
            abs(nearest - wavelength) > 2
            or abs(averages[nearest] - factor) > FIG1_TOLERANCE
        ):
            msg = f"Fig. 1(b) at {wavelength} nm: {averages[nearest]}, legend {factor}"
            raise ValueError(msg)
    for letter, (wavelength, factor) in FIG2.items():
        for intensity, ppt, ratio in inset_factors(letter):
            if ppt >= FACTOR_FROM and abs(ratio - factor) > FACTOR_TOLERANCE:
                msg = (
                    f"Fig. 2({letter}) inset at {intensity:.3e}: factor {ratio:.4f}, "
                    f"legend {factor} ({wavelength} nm)"
                )
                raise ValueError(msg)
        for inset_i, main_i, inset_p, main_p in inset_comparison(letter):
            if abs(inset_i / main_i - 1) > GRID:
                msg = f"Fig. 2({letter}): inset point at {inset_i:.3e} has no main one"
                raise ValueError(msg)
            if abs(inset_p - main_p) > max(INSET_AGREEMENT * main_p, INSET_ABSOLUTE):
                msg = (
                    f"Fig. 2({letter}) at {main_i:.3e}: inset {inset_p:.4g}, "
                    f"main {main_p:.4g}"
                )
                raise ValueError(msg)
        main = tdse(f"fig2{letter}")
        grid = [2e13 + 5e12 * k for k in range(25)]
        if len(main) != len(grid) or any(
            abs(r.abscissa / g - 1) > GRID for r, g in zip(main, grid, strict=True)
        ):
            msg = f"Fig. 2({letter}): the TDSE points are not on the 0.5e13 grid"
            raise ValueError(msg)
    for key, found in arrows().items():
        for gamma, read, expected in found:
            if abs(read / expected - 1) > ARROW:
                msg = f"{key}: gamma = {gamma} at {read:.3e}, expected {expected:.3e}"
                raise ValueError(msg)
    if abs(one.mean - 1.49) > 0.005:
        msg = f"Fig. 1(a): the curve averages {one.mean}, not 1.49"
        raise ValueError(msg)


def _points(key: str, found: list[Reading]) -> list[Point]:
    """Return a dataset's points, numbered in order of intensity."""
    return [
        Point(
            f"p{k:02d}",
            f"{r.abscissa:.4e}",
            f"{r.value:.4g}",
            "",
            f"{PANELS[key].name} TDSE: marker at ({r.x:.3f}, {r.y:.3f}) pt",
        )
        for k, r in enumerate(found, start=1)
    ]


def extract() -> dict[str, list[Point]]:
    """Return the TDSE probabilities of Figs. 2 and 3, one file per panel."""
    checks()
    return {relative: _points(key, tdse(key)) for key, relative in DATASETS.items()}


def main() -> None:
    """Write the data files and print the checks' numbers."""
    write(extract())
    one = fig1()
    print(f"Fig. 1(a): constant {one.constant:.4f}, curve mean {one.mean:.4f}")
    curve = np.array(one.curve)
    print(f"  curve {curve[:, 1].min():.3f} to {curve[:, 1].max():.3f}")
    for name in ("wavelength", "ionisation"):
        print(name, [(round(a, 2), round(v, 3)) for a, v in getattr(one, name)])
        fit = getattr(one, f"{name}_fit")
        print(" fit", [(round(a, 2), round(v, 3)) for a, v in fit])
    for letter in FIG2:
        print(letter, "PPTX factor", [round(r, 4) for _, _, r in inset_factors(letter)])
        print(
            letter,
            "inset/main - 1",
            [round(i / m - 1, 4) for _, _, i, m in inset_comparison(letter)],
        )
    for key, found in arrows().items():
        print(
            key, [(g, f"{r:.3e}", f"{e:.3e}", round(r / e - 1, 3)) for g, r, e in found]
        )
    print("Fig. 1(b), (c) fits at the axes' ends", fit_ends())
    for letter, (wavelength, factor) in FIG2.items():
        ratios = tdse_over_ppt(letter)
        inset = [p for i, p, _ in ratios if INSET_RANGE[0] <= i <= INSET_RANGE[1]]
        print(
            f"{wavelength} nm TDSE/PPT (intensity, P ratio, rate ratio); mean over "
            f"6-13e13 {np.mean(inset):.3f}, legend {factor}:",
            [(f"{i:.2e}", round(p, 3), round(w, 3)) for i, p, w in ratios],
        )


if __name__ == "__main__":
    main()
