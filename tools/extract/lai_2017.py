# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Lai et al. (2017): measured ion yields of Ne, Ar, Kr and Xe, 0.4 to 3.9 um.

The paper prints no yield as a number: the measured singly charged ion yields against
the peak intensity are in Figs. 1 to 4 only (data rule 8, TC-699), all four on PDF
page 2, drawn as vectors (MATLAB's: the tick labels are text, the markers stroked
paths, the theory curves polylines). The measured points are read; the PPT and ADK
curves (the authors' calculations of others' formulas) are not. Na, K, Zn and Mg
(Figs. 5 to 8, 11, 15) are outside optmat's media and are not read.

* **Fig. 1** (library figure ``fig_31008d4c7a41``): Xe, (a) 0.8 um, linear (squares,
  'Expt.-LP') and circular (circles, 'Expt.-CP'); (b) 3.3 um, linear only.
* **Fig. 2** (``fig_0d4d59b034d0``): Ar, (a) 0.4, (b) 0.8, (c) 1.3 um, linear and
  circular (the legend is in panel (c) only).
* **Fig. 3** (``fig_4b3af4da8236``): Kr, (a) 0.8 um, linear and circular; (b) 3.3
  (squares), 3.6 (circles) and 3.9 um (triangles), 'displaced from each other
  arbitrarily along the vertical axis' (the caption); their polarisation is not stated.
* **Fig. 4** (``fig_73722d72eb42``): Ne, 0.8 um, linear and circular.

How the drawing is read (:func:`_figure.painted`, :func:`_figure.texts`):

* the panels: each plot area is a white-filled rectangle (:data:`PANELS` gives its
  lower-left corner to a point, and the script checks that exactly one area is there);
  the panel's letter is the text '(a)', '(b)' or '(c)' inside it (Fig. 4 has none).
* the axes: the ticks are black two-vertex lines across the area's bottom (x) and left
  (y) edges, one end on the edge; the top and right edges repeat them (not used). Both
  axes are logarithmic. The x axis's long ticks are the decades, labelled '10' with a
  superscript exponent (text objects, read here: :func:`_labels`); its short ticks are
  2 to 9 times the decade below (those before the first decade counted down from 9),
  a short tick also drawn at each decade being dropped. The y axis has long ticks only,
  every second decade, each labelled. Figs. 1(b) and 2(b) print no y labels: their
  ticks are at the left panel's positions (checked, :data:`SHARED`), whose labels
  they share. Every tick
  must lie on the fitted log scale within :data:`TICK_TOLERANCE`
  (:meth:`_figure.Axis.check`).
* the markers: black stroked, unfilled paths under :data:`MARKER_SIZE` across: squares
  (five vertices, closed), circles (13, four Bezier arcs) and triangles (four, closed,
  pointing up). The centre of a square or circle is the middle of its vertices' extent;
  of a triangle, the mean of its three corners (its centroid: MATLAB centres the '^'
  marker there; in Fig. 3(b)'s legend the triangle's centroid lies on its row,
  0.02 pt from the rows' even spacing, where the middle of its height lies 0.56 pt
  above). Markers outside the plot area are dropped. The legend's markers are the
  ones beside a legend text 'Expt...' (:data:`LEGEND_LEFT`, :data:`LEGEND_ROW`); each
  such text must have exactly one, of the shape :data:`PANELS` assigns that series
  (so the legend itself confirms which shape is which polarisation or wavelength).

Every point's ``printed`` cell gives the figure, the panel and the marker's centre in
page points. The rows are numbered in order of intensity.

**Checks**:

* (:func:`checks`) the extracted intensity range of each gas and wavelength against
  Table I's I_min and I_max (``tbl_73c01a54f6cc``, read from its cells,
  :func:`_json_table.json_table`), printed to two or three figures
  (:data:`RANGE_AGREEMENT`). Table I's ranges are the linear sets' (and Fig. 3(b)'s):
  the circular sets reach beyond them and are not compared.
* (:func:`fig9_check`) Fig. 9(a) (PDF p. 5) replots Xe's 0.8 um yields on a linear
  intensity axis: its 37 squares must be Fig. 1(a)'s 37 LP points, read from a second,
  independently calibrated drawing (:data:`FIG9_INTENSITY`, :data:`FIG9_YIELD`).
* (:func:`ratio_check`, reported by :func:`main`, not enforced) Fig. 13(a) (PDF p. 6):
  the ratio Y_CP/Y_LP at I_BSI (Table II, ``tbl_67a91b746408``) that the paper derives
  from these curves, against the same ratio from the extracted LP and CP points, each
  set's log yield fitted by a straight line in log intensity through the
  :data:`NEAREST` points nearest I_BSI. Fig. 13 is not taken as data.
"""

import math
from collections import defaultdict
from dataclasses import dataclass
from itertools import pairwise

import numpy as np
from _figure import Axis, Painted, painted, texts
from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Lai_2017_ExperimentalInvestigationStrong"
SOURCE = "lai_2017"
PAGE = 2
BLACK, WHITE = (0, 0, 0), (255, 255, 255)
MINUS = "\N{MINUS SIGN}"
TW = 1e12  # the figures' and Table I's TW/cm^2, in W/cm^2
# The worst tick of any panel lies 0.04 pt off its axis's fitted log scale (MATLAB
# rounds the tick positions; the y ticks are 22 to 30 pt per decade): 0.05 pt allowed.
TICK_TOLERANCE = 0.05
# A tick's end on the edge, and a tick's straightness, to the drawing's precision.
ON_EDGE = 0.01
# Long ticks are 3.2 to 5.3 pt, short ones 1.7 to 2.6 pt: 1.5 times the shortest.
LONG = 1.5
# Markers are 2.8 to 4.0 pt across their vertices; tick labels' and legend lines'
# paths are larger or of other vertex counts.
MARKER_SIZE = 6.0
# A legend marker lies 20 to 30 pt left of its text's centre, on its row.
LEGEND_LEFT = (10.0, 45.0)
LEGEND_ROW = 2.0
# A label '10' lies within 3 pt of its tick (along the axis) and within 22 pt of the
# edge (across it); its exponent within 4 pt of 5.5 pt right, 4 pt up of its centre.
LABEL_ALONG = 3.0
LABEL_ACROSS = 22.0
EXPONENT_OFFSET = (5.5, 4.0)
EXPONENT_NEAR = 4.0
# Figs. 1(b) and 2(b) print no y labels and share the left panel's y ticks: the same
# positions to 0.001 pt.
SHARES_Y = {"fig1b": "fig1a", "fig2b": "fig2a"}
SHARED = 0.001
# Table I prints I_min and I_max to two or three figures (e.g. 26, 240, 1030 TW/cm^2:
# up to 2 % rounding); the extracted ranges agree within 3 %, except Ar at 1.3 um,
# 59.6 to 239.9 against 58 to 230 (+2.8 %, +4.3 %): 5 % allowed.
RANGE_AGREEMENT = 0.05

# Fig. 3(b) draws every data marker twice, the copies within 0.002 pt of each other
# (the legend's once); the other panels once. Two different points of a set are never
# closer than 0.058 pt (repeated measurements at one intensity, apart in the yield).
DRAWN = {"fig3b": 2}
COPY = 0.01

SQUARE, CIRCLE, TRIANGLE = "square", "circle", "triangle"
VERTICES = {SQUARE: 5, CIRCLE: 13, TRIANGLE: 4}


@dataclass(frozen=True)
class Series:
    """One measured series of a panel: its marker, legend text and dataset."""

    shape: str
    legend: str
    wavelength_nm: int
    polarization: str


@dataclass(frozen=True)
class PanelSpec:
    """A panel as the paper draws it: figure, letter, gas, area corner, series."""

    figure: int
    letter: str
    gas: str
    corner: tuple[float, float]
    series: tuple[Series, ...]
    legend_panel: str


def _lp_cp(wavelength_nm: int) -> tuple[Series, ...]:
    """Return a panel's linear (squares) and circular (circles) series."""
    return (
        Series(SQUARE, f"Expt.{MINUS}LP", wavelength_nm, "linear"),
        Series(CIRCLE, f"Expt.{MINUS}CP", wavelength_nm, "circular"),
    )


PANELS = (
    PanelSpec(1, "(a)", "Xe", (80.7, 573.0), _lp_cp(800), "fig1a"),
    PanelSpec(
        1,
        "(b)",
        "Xe",
        (190.3, 573.0),
        (Series(SQUARE, f"Expt.{MINUS}LP", 3300, "linear"),),
        "fig1b",
    ),
    PanelSpec(2, "(a)", "Ar", (74.3, 229.1), _lp_cp(400), "fig2c"),
    PanelSpec(2, "(b)", "Ar", (184.6, 229.1), _lp_cp(800), "fig2c"),
    PanelSpec(2, "(c)", "Ar", (74.3, 107.8), _lp_cp(1300), "fig2c"),
    PanelSpec(3, "(a)", "Kr", (339.7, 580.0), _lp_cp(800), "fig3a"),
    PanelSpec(
        3,
        "(b)",
        "Kr",
        (457.4, 579.9),
        (
            Series(SQUARE, "Expt.(3.3)", 3300, "not stated"),
            Series(CIRCLE, "Expt.(3.6)", 3600, "not stated"),
            Series(TRIANGLE, "Expt.(3.9)", 3900, "not stated"),
        ),
        "fig3b",
    ),
    PanelSpec(4, "", "Ne", (363.0, 101.8), _lp_cp(800), "fig4"),
)
FIGURE_ID = {
    1: "fig_31008d4c7a41",
    2: "fig_0d4d59b034d0",
    3: "fig_4b3af4da8236",
    4: "fig_73722d72eb42",
}


def _name(spec: PanelSpec) -> str:
    """Return a panel's short name, ``fig1a`` (``fig4`` for the one-panel Fig. 4)."""
    return f"fig{spec.figure}{spec.letter.strip('()')}"


def dataset_id(spec: PanelSpec, series: Series) -> str:
    """Return a dataset's id: figure, panel, gas, wavelength and polarisation."""
    tag = {"linear": "LP", "circular": "CP", "not stated": ""}[series.polarization]
    wavelength = f"{series.wavelength_nm / 1000:g}um".replace(".", "p")
    return "_".join(p for p in (_name(spec), spec.gas, wavelength, tag) if p)


@dataclass(frozen=True)
class Panel:
    """A panel's plot area (l, b, r, t, pt) and its calibrated axes."""

    spec: PanelSpec
    area: tuple[float, float, float, float]
    x_axis: Axis
    y_axis: Axis


def _box(drawn: Painted) -> tuple[float, float, float, float]:
    """Return a path's vertices' extent (l, b, r, t), without the stroke."""
    xs, ys = drawn.vertices[:, 0], drawn.vertices[:, 1]
    return float(xs.min()), float(ys.min()), float(xs.max()), float(ys.max())


def _inside(x: float, y: float, box: tuple[float, float, float, float]) -> bool:
    """Return whether a point lies inside a box (l, b, r, t)."""
    return box[0] <= x <= box[2] and box[1] <= y <= box[3]


def _ticks(
    drawn: list[Painted], area: tuple[float, float, float, float], axis: int
) -> tuple[list[float], list[float]]:
    """Return the long and short tick positions on the area's bottom or left edge.

    ``axis`` 0 is x (vertical ticks on the bottom edge), 1 is y (horizontal ticks on
    the left edge).
    """
    edge = area[1] if axis == 0 else area[0]
    along, across = (0, 1) if axis == 0 else (1, 0)
    lo, hi = (area[0], area[2]) if axis == 0 else (area[1], area[3])
    found: list[tuple[float, float]] = []
    for d in drawn:
        if d.stroke != BLACK or d.fill is not None or len(d.vertices) != 2:
            continue
        line = d.vertices
        if abs(line[0, along] - line[1, along]) > ON_EDGE:
            continue
        ends = sorted(line[:, across])
        if abs(ends[0] - edge) > ON_EDGE or ends[1] - ends[0] > 10:
            continue
        position = float(line[0, along])
        if lo - ON_EDGE <= position <= hi + ON_EDGE:
            found.append((round(position, 3), float(ends[1] - ends[0])))
    shortest = min(length for _, length in found)
    if max(length for _, length in found) <= LONG * shortest:
        # One length only (the y axes): every tick is a labelled one.
        return sorted({p for p, _ in found}), []
    long = sorted({p for p, length in found if length > LONG * shortest})
    short = sorted({p for p, length in found if length <= LONG * shortest})
    short = [s for s in short if min(abs(s - p) for p in long) > ON_EDGE]
    return long, short


def _exponent(text: str) -> int:
    """Return a printed exponent ('2', '-10' with a Unicode minus) as an integer."""
    return int(number(text))


def _labels(
    words: list[tuple[str, float, float]],
    positions: list[float],
    area: tuple[float, float, float, float],
    axis: int,
) -> list[int]:
    """Return the decade printed at each long tick: '10' and its superscript.

    ``axis`` 0 is x (labels below the bottom edge), 1 is y (left of the left edge).
    """
    decades = []
    for p in positions:
        if axis == 0:
            bases = [
                (x, y)
                for w, x, y in words
                if w == "10"
                and abs(x - p) < LABEL_ALONG
                and 0 < area[1] - y < LABEL_ACROSS
            ]
        else:
            bases = [
                (x, y)
                for w, x, y in words
                if w == "10"
                and abs(y - p) < LABEL_ALONG
                and 0 < area[0] - x < LABEL_ACROSS
            ]
        if len(bases) != 1:
            msg = f"{len(bases)} labels at the tick {p:.3f} pt"
            raise ValueError(msg)
        x, y = bases[0]
        target = (x + EXPONENT_OFFSET[0], y + EXPONENT_OFFSET[1])
        exponents = [
            w
            for w, u, v in words
            if w.lstrip(MINUS).isdigit()
            and math.hypot(u - target[0], v - target[1]) < EXPONENT_NEAR
        ]
        if len(exponents) != 1:
            msg = f"{len(exponents)} exponents at the label at ({x:.1f}, {y:.1f}) pt"
            raise ValueError(msg)
        decades.append(_exponent(exponents[0]))
    return decades


def _x_axis(long: list[float], short: list[float], decades: list[int]) -> Axis:
    """Return a log x axis from its decades (long ticks) and the short ticks between.

    The short ticks between two decades are 2 to 9 times the lower; those before the
    first decade are counted down from 9 times the decade below it, those after the
    last up from 2 times it.
    """
    if any(b - a != 1 for a, b in pairwise(decades)):
        msg = f"the x decades {decades} are not consecutive"
        raise ValueError(msg)
    positions = list(long)
    values = [10.0**d for d in decades]
    below = sorted((s for s in short if s < long[0]), reverse=True)
    above = sorted(s for s in short if s > long[-1])
    if len(below) > 8 or len(above) > 8:
        msg = "more short ticks beyond the decades than a decade holds"
        raise ValueError(msg)
    for k, s in enumerate(below):
        positions.append(s)
        values.append((9 - k) * 10.0 ** (decades[0] - 1))
    for k, s in enumerate(above):
        positions.append(s)
        values.append((2 + k) * 10.0 ** decades[-1])
    for lower, upper, decade in zip(long, long[1:], decades, strict=False):
        between = sorted(s for s in short if lower < s < upper)
        if len(between) != 8:
            msg = f"{len(between)} short ticks in a decade"
            raise ValueError(msg)
        positions += between
        values += [k * 10.0**decade for k in range(2, 10)]
    order = np.argsort(positions)
    axis = Axis(
        tuple(float(positions[i]) for i in order),
        tuple(float(values[i]) for i in order),
        log=True,
        tolerance=TICK_TOLERANCE,
    )
    axis.check()
    return axis


def _y_axis(long: list[float], short: list[float], decades: list[int]) -> Axis:
    """Return a log y axis from its labelled long ticks (every second decade)."""
    if short:
        msg = f"unexpected short y ticks at {short}"
        raise ValueError(msg)
    axis = Axis(
        tuple(long),
        tuple(10.0**d for d in decades),
        log=True,
        tolerance=TICK_TOLERANCE,
    )
    axis.check()
    return axis


def _drawn() -> tuple[list[Painted], list[tuple[str, float, float]]]:
    """Return page 2's paths and text objects."""
    region = (0.0, 0.0, 610.0, 810.0)
    return painted(DIRECTORY, PAGE, region), texts(DIRECTORY, PAGE, region)


def panels() -> list[Panel]:
    """Return the eight panels, their areas found and their axes calibrated."""
    drawn, words = _drawn()
    areas = [
        _box(d) for d in drawn if d.fill == WHITE and d.stroke is None and d.size > 90
    ]
    out: list[Panel] = []
    y_ticks: dict[str, list[float]] = {}
    y_decades: dict[str, list[int]] = {}
    for spec in PANELS:
        mine = [
            a
            for a in set(areas)
            if abs(a[0] - spec.corner[0]) < 1 and abs(a[1] - spec.corner[1]) < 1
        ]
        if len(mine) != 1:
            msg = f"{_name(spec)}: {len(mine)} plot areas at {spec.corner}"
            raise ValueError(msg)
        area = mine[0]
        letters = [
            w for w, x, y in words if w in {"(a)", "(b)", "(c)"} and _inside(x, y, area)
        ]
        if letters != ([spec.letter] if spec.letter else []):
            msg = f"{_name(spec)}: the panel's letter is {letters}"
            raise ValueError(msg)
        x_long, x_short = _ticks(drawn, area, 0)
        y_long, y_short = _ticks(drawn, area, 1)
        x_axis = _x_axis(x_long, x_short, _labels(words, x_long, area, 0))
        if _name(spec) in SHARES_Y:
            # Figs. 1(b) and 2(b) print no y labels: they share the left panel's.
            left = SHARES_Y[_name(spec)]
            shared = y_ticks[left]
            if len(shared) != len(y_long) or any(
                abs(a - b) > SHARED for a, b in zip(shared, y_long, strict=True)
            ):
                msg = f"{_name(spec)}: the y ticks are not {left}'s"
                raise ValueError(msg)
            decades = y_decades[left]
        else:
            decades = _labels(words, y_long, area, 1)
        y_ticks[_name(spec)], y_decades[_name(spec)] = y_long, decades
        out.append(Panel(spec, area, x_axis, _y_axis(y_long, y_short, decades)))
    return out


@dataclass(frozen=True)
class Marker:
    """A marker: its shape and its centre in page points."""

    shape: str
    x: float
    y: float


def _shape(d: Painted) -> str | None:
    """Return a stroked path's marker shape, or None if it is not a marker."""
    if d.stroke != BLACK or d.fill is not None or d.size > MARKER_SIZE:
        return None
    for shape, count in VERTICES.items():
        if len(d.vertices) == count and np.allclose(d.vertices[0], d.vertices[-1]):
            return shape
    return None


def _centre(d: Painted, shape: str) -> tuple[float, float]:
    """Return a marker's centre: a triangle's centroid, else its extent's middle."""
    if shape == TRIANGLE:
        corners = d.vertices[:3]
        if not corners[2, 1] > corners[0, 1]:
            msg = f"a triangle at {d.centre} does not point up"
            raise ValueError(msg)
        return float(corners[:, 0].mean()), float(corners[:, 1].mean())
    left, bottom, right, top = _box(d)
    return (left + right) / 2, (bottom + top) / 2


def _markers(
    panel: Panel,
    drawn: list[Painted],
    words: list[tuple[str, float, float]],
) -> dict[str, list[Marker]]:
    """Return a panel's data markers by series legend, the legend's own removed."""
    found = []
    for d in drawn:
        shape = _shape(d)
        if shape is None:
            continue
        x, y = _centre(d, shape)
        if _inside(x, y, panel.area):
            found.append(Marker(shape, x, y))
    legend = [
        (w, x, y)
        for w, x, y in words
        if w.startswith("Expt") and _inside(x, y, panel.area)
    ]
    in_legend: set[Marker] = set()
    shapes: dict[str, str] = {}
    for text, x, y in legend:
        beside = [
            m
            for m in found
            if LEGEND_LEFT[0] < x - m.x < LEGEND_LEFT[1] and abs(m.y - y) < LEGEND_ROW
        ]
        if len(beside) != 1:
            msg = f"{_name(panel.spec)}: {len(beside)} markers beside '{text}'"
            raise ValueError(msg)
        in_legend.add(beside[0])
        shapes[text] = beside[0].shape
    by_series: dict[str, list[Marker]] = defaultdict(list)
    for m in found:
        if m not in in_legend:
            by_series[m.shape].append(m)
    times = DRAWN.get(_name(panel.spec), 1)
    for shape, markers in by_series.items():
        by_series[shape] = _once(markers, times, _name(panel.spec))
    out = {}
    for series in panel.spec.series:
        if legend and shapes.get(series.legend) != series.shape:
            msg = f"{_name(panel.spec)}: '{series.legend}' is not a {series.shape}"
            raise ValueError(msg)
        out[series.legend] = sorted(by_series.pop(series.shape, []), key=lambda m: m.x)
    if by_series:
        msg = f"{_name(panel.spec)}: markers of no series: {dict(by_series)}"
        raise ValueError(msg)
    return out


def _once(markers: list[Marker], times: int, where: str) -> list[Marker]:
    """Return markers drawn ``times`` times each at one place, once each.

    Raises unless every marker has exactly ``times`` copies within :data:`COPY`.
    """
    groups: list[list[Marker]] = []
    for m in sorted(markers, key=lambda m: (m.x, m.y)):
        for group in groups:
            if math.hypot(m.x - group[0].x, m.y - group[0].y) < COPY:
                group.append(m)
                break
        else:
            groups.append([m])
    if any(len(group) != times for group in groups):
        counts = sorted({len(group) for group in groups})
        msg = f"{where}: markers drawn {counts} times, expected {times}"
        raise ValueError(msg)
    return [group[0] for group in groups]


def legends() -> dict[str, dict[str, str]]:
    """Return each legend's text-to-shape assignment, as read from the drawing."""
    drawn, words = _drawn()
    out = {}
    for panel in panels():
        name = _name(panel.spec)
        if panel.spec.legend_panel != name:
            continue
        _markers(panel, drawn, words)  # raises unless the legend matches PANELS
        out[name] = {s.legend: s.shape for s in panel.spec.series}
    return out


@dataclass(frozen=True)
class Reading:
    """One measured point: the marker's centre (pt) and its intensity and yield."""

    x: float
    y: float
    intensity: float
    yield_: float


def readings() -> dict[str, list[Reading]]:
    """Return every dataset's points, by dataset id, in order of intensity."""
    drawn, words = _drawn()
    out = {}
    for panel in panels():
        markers = _markers(panel, drawn, words)
        for series in panel.spec.series:
            out[dataset_id(panel.spec, series)] = [
                Reading(
                    m.x,
                    m.y,
                    float(panel.x_axis(m.x)) * TW,
                    float(panel.y_axis(m.y)),
                )
                for m in markers[series.legend]
            ]
    return out


TABLE1_HEADER = [
    "",
    "\\lambda(μm)",
    "I_{\\text{min}}(TW/cm2)",
    "I_{\\text{max}}(TW/cm2)",
    "\\gamma_{\\text{BSI}}",
    "\\gamma_{\\text{max}}",
]


def table1() -> dict[tuple[str, int], tuple[float, float]]:
    """Return Table I's I_min and I_max (W/cm^2) by gas and wavelength (nm)."""
    grid = json_table(DIRECTORY, "tbl_73c01a54f6cc")
    if grid[0] != TABLE1_HEADER:
        msg = f"unexpected header of Table I: {grid[0]!r}"
        raise ValueError(msg)
    out = {}
    gas = ""
    for cells in grid[1:]:
        gas = cells[0] or gas
        wavelength = round(float(number(cells[1])) * 1000)
        low, high = float(number(cells[2])), float(number(cells[3]))
        out[gas, wavelength] = (low * TW, high * TW)
    return out


def ranges() -> list[tuple[str, int, float, float, float, float]]:
    """Return, per gas and wavelength, the extracted and Table I's intensity ranges.

    Table I's ranges are the linearly polarised sets' (and Fig. 3(b)'s): the circular
    sets reach beyond them (Xe 0.8 um to 188, Ar 1.3 um to 282 TW/cm^2) and are not
    compared.
    """
    found: dict[tuple[str, int], list[float]] = defaultdict(list)
    for spec in PANELS:
        for series in spec.series:
            if series.polarization == "circular":
                continue
            for r in readings_cache()[dataset_id(spec, series)]:
                found[spec.gas, series.wavelength_nm].append(r.intensity)
    printed = table1()
    out = []
    for key, values in found.items():
        low, high = printed[key]
        out.append((key[0], key[1], min(values), max(values), low, high))
    return out


_CACHE: dict[str, list[Reading]] = {}


def readings_cache() -> dict[str, list[Reading]]:
    """Return :func:`readings`, read once."""
    if not _CACHE:
        _CACHE.update(readings())
    return _CACHE


def checks() -> None:
    """Raise unless each set's intensity range agrees with Table I's."""
    for gas, wavelength, low, high, t_low, t_high in ranges():
        for read, printed in ((low, t_low), (high, t_high)):
            if abs(read / printed - 1) > RANGE_AGREEMENT:
                msg = (
                    f"{gas} {wavelength} nm: the extracted range {low:.3g} to "
                    f"{high:.3g} W/cm^2 is not Table I's {t_low:.3g} to {t_high:.3g}"
                )
                raise ValueError(msg)


# Fig. 13(a) (PDF p. 6, library figure fig_55247f967b87): Y_CP/Y_LP at I_BSI, the
# noble gases; filled markers, coloured by wavelength (the caption: blue 0.4 um, red
# 0.8, green 1.3), shaped by gas (the legend's grey markers: circle Ar, left triangle
# Xe, up triangle Kr, diamond Ne). The plot area's left edge and its y ticks, 0 to 0.6.
FIG13_PAGE = 6
FIG13_REGION = (60.0, 255.0, 300.0, 345.0)
FIG13_LEFT = 98.394
FIG13_TICK = 1.8  # the y ticks are 1.70 pt long
FIG13_LABEL = 12.0  # the labels' centres lie 4 to 8 pt left of the edge
BLUE, RED, GREEN, GREY = (0, 0, 255), (255, 0, 0), (0, 204, 0), (128, 128, 128)
FIG13: dict[tuple[tuple[int, int, int] | None, int, str], str] = {
    # (fill, vertices, pointing): the dataset ids whose LP and CP sets give the ratio
    (BLUE, 13, ""): "fig2a_Ar_0p4um",
    (RED, 13, ""): "fig2b_Ar_0p8um",
    (GREEN, 13, ""): "fig2c_Ar_1p3um",
    (RED, 4, "left"): "fig1a_Xe_0p8um",
    (RED, 4, "up"): "fig3a_Kr_0p8um",
    (RED, 5, ""): "fig4_Ne_0p8um",
}
# Table II's I_BSI (TW/cm^2) by gas (tbl_67a91b746408).
STAR = "\N{ASTERISK OPERATOR}"
TABLE2_HEADER = ["", "Ip (eV)", f"n{STAR}", f"cn2{STAR}l{STAR}", "IBSI (TW/cm2)"]
# The ratio from the curves: each set's log yield against log intensity fitted with a
# straight line through the NEAREST points to I_BSI (by log intensity).
NEAREST = 5


def _fig13_axis() -> Axis:
    """Return Fig. 13(a)'s linear y axis, from its tick marks and their labels."""
    drawn = painted(DIRECTORY, FIG13_PAGE, FIG13_REGION)
    words = texts(DIRECTORY, FIG13_PAGE, FIG13_REGION)
    ticks = sorted(
        round(float(d.vertices[0, 1]), 3)
        for d in drawn
        if d.stroke == BLACK
        and len(d.vertices) == 2
        and abs(d.vertices[0, 1] - d.vertices[1, 1]) < ON_EDGE
        and abs(min(d.vertices[:, 0]) - FIG13_LEFT) < ON_EDGE
        and max(d.vertices[:, 0]) - min(d.vertices[:, 0]) < FIG13_TICK
    )
    values = []
    for t in ticks:
        labels = [
            w
            for w, x, y in words
            if w and abs(y - t) < LABEL_ALONG and 0 < FIG13_LEFT - x < FIG13_LABEL
        ]
        if len(labels) != 1:
            msg = f"Fig. 13(a): {len(labels)} labels at the tick {t:.3f} pt"
            raise ValueError(msg)
        values.append(float(number(labels[0])))
    axis = Axis(tuple(ticks), tuple(values), tolerance=TICK_TOLERANCE)
    axis.check()
    return axis


def _pointing(d: Painted) -> str:
    """Return where a filled triangle points, 'up', 'down', 'left' or 'right'.

    The base is the side whose two corners share a coordinate; '' if not a triangle.
    """
    if len(d.vertices) != 4:
        return ""
    corners = d.vertices[:3]
    for axis, (low, high) in ((0, ("left", "right")), (1, ("down", "up"))):
        for k in range(3):
            base = [corners[j, axis] for j in range(3) if j != k]
            if abs(base[0] - base[1]) < ON_EDGE:
                return low if corners[k, axis] < base[0] else high
    return ""


def fig13() -> dict[str, float]:
    """Return Fig. 13(a)'s ratios, by the stem of the LP and CP sets they come from."""
    axis = _fig13_axis()
    drawn = painted(DIRECTORY, FIG13_PAGE, FIG13_REGION)
    out = {}
    for d in drawn:
        key = (d.fill, len(d.vertices), _pointing(d))
        if d.stroke is not None or key not in FIG13:
            continue
        if len(d.vertices) == 4:
            y = float(d.vertices[:3, 1].mean())
        else:
            y = float((d.vertices[:, 1].min() + d.vertices[:, 1].max()) / 2)
        if FIG13[key] in out:
            msg = f"Fig. 13(a): two markers for {FIG13[key]}"
            raise ValueError(msg)
        out[FIG13[key]] = float(axis(y))
    if len(out) != len(FIG13):
        msg = f"Fig. 13(a): {len(out)} markers, expected {len(FIG13)}"
        raise ValueError(msg)
    return out


def table2() -> dict[str, float]:
    """Return Table II's I_BSI (W/cm^2) by gas."""
    grid = table(DIRECTORY, "tbl_67a91b746408")
    if grid[0] != TABLE2_HEADER:
        msg = f"unexpected header of Table II: {grid[0]!r}"
        raise ValueError(msg)
    return {cells[0]: float(number(cells[4])) * TW for cells in grid[1:]}


def _at(found: list[Reading], intensity: float) -> float:
    """Return a set's yield at an intensity, from a line through the nearest points."""
    x = np.log([r.intensity for r in found])
    y = np.log([r.yield_ for r in found])
    nearest = np.argsort(np.abs(x - math.log(intensity)))[:NEAREST]
    slope, intercept = np.polyfit(x[nearest], y[nearest], 1)
    return float(np.exp(slope * math.log(intensity) + intercept))


def ratio_check() -> list[tuple[str, float, float]]:
    """Return, per set, Y_CP/Y_LP at I_BSI from the extracted curves and Fig. 13(a)'s.

    Reported, not enforced: the paper's ratios are not taken as data (the record's
    notes), and they differ from the curves' by more than the reading.
    """
    found = readings_cache()
    bsi = table2()
    out = []
    for stem, printed in fig13().items():
        gas = stem.split("_")[1]
        lp, cp = found[f"{stem}_LP"], found[f"{stem}_CP"]
        ratio = _at(cp, bsi[gas]) / _at(lp, bsi[gas])
        out.append((stem, ratio, printed))
    return out


# Fig. 9(a) (PDF p. 5, library figure fig_826039bb622f): Xe's 0.8 um yield again,
# 'Expt.' squares, against a linear intensity axis (its labels under panel (b), whose
# ticks it shares). Its markers must be Fig. 1(a)'s LP points, one for one.
FIG9_PAGE = 5
FIG9_REGION = (40.0, 520.0, 300.0, 760.0)
FIG9_SET = "fig1a_Xe_0p8um_LP"
# The two drawings agree within 0.29 % in intensity and 0.76 % in yield (the worst
# point: 0.1 pt of Fig. 9's 30 pt per decade); 0.5 % and 1.5 % allowed.
FIG9_INTENSITY = 0.005
FIG9_YIELD = 0.015


def fig9() -> list[tuple[float, float]]:
    """Return Fig. 9(a)'s markers, (intensity W/cm^2, yield), in order of intensity."""
    drawn = painted(DIRECTORY, FIG9_PAGE, FIG9_REGION)
    words = texts(DIRECTORY, FIG9_PAGE, FIG9_REGION)
    areas = sorted(
        {
            _box(d)
            for d in drawn
            if d.fill == WHITE and d.stroke is None and d.size > 90
        },
        key=lambda a: -a[1],
    )
    if len(areas) != 2:
        msg = f"Fig. 9: {len(areas)} plot areas"
        raise ValueError(msg)
    top, bottom = areas
    x_long, _ = _ticks(drawn, bottom, 0)
    top_long = _ticks(drawn, top, 0)[0]
    if len(top_long) != len(x_long) or any(
        abs(a - b) > 2 * SHARED for a, b in zip(top_long, x_long, strict=True)
    ):
        msg = "Fig. 9: panel (a)'s x ticks are not panel (b)'s"
        raise ValueError(msg)
    x_values = []
    for p in x_long:
        labels = [
            w
            for w, x, y in words
            if w.isdigit()
            and abs(x - p) < LABEL_ALONG
            and 0 < bottom[1] - y < LABEL_ACROSS
        ]
        if len(labels) != 1:
            msg = f"Fig. 9: {len(labels)} labels at the x tick {p:.3f} pt"
            raise ValueError(msg)
        x_values.append(float(labels[0]))
    x_axis = Axis(tuple(x_long), tuple(x_values), tolerance=TICK_TOLERANCE)
    x_axis.check()
    y_long, y_short = _ticks(drawn, top, 1)
    y_axis = _y_axis(y_long, y_short, _labels(words, y_long, top, 1))
    found = [
        _centre(d, SQUARE)
        for d in drawn
        if _shape(d) == SQUARE and _inside(*_centre(d, SQUARE), top)
    ]
    legend = [(x, y) for w, x, y in words if w == "Expt." and _inside(x, y, top)]
    data = [
        (x, y)
        for x, y in found
        if not any(
            LEGEND_LEFT[0] < u - x < LEGEND_LEFT[1] and abs(v - y) < LEGEND_ROW
            for u, v in legend
        )
    ]
    if len(found) - len(data) != len(legend):
        msg = "Fig. 9: the legend's marker is not found"
        raise ValueError(msg)
    return sorted((float(x_axis(x)) * TW, float(y_axis(y))) for x, y in data)


def fig9_check() -> tuple[float, float]:
    """Return the worst relative differences, intensity and yield, Fig. 9(a) to 1(a).

    Raises unless Fig. 9(a) has as many markers as Fig. 1(a)'s LP set and each pair
    (in order of intensity) agrees within :data:`FIG9_INTENSITY` and
    :data:`FIG9_YIELD`.
    """
    mine = readings_cache()[FIG9_SET]
    other = fig9()
    if len(other) != len(mine):
        msg = f"Fig. 9(a) has {len(other)} markers, Fig. 1(a) LP {len(mine)}"
        raise ValueError(msg)
    worst_i = max(
        abs(i / r.intensity - 1) for (i, _), r in zip(other, mine, strict=True)
    )
    worst_y = max(abs(y / r.yield_ - 1) for (_, y), r in zip(other, mine, strict=True))
    if worst_i > FIG9_INTENSITY or worst_y > FIG9_YIELD:
        msg = f"Fig. 9(a) differs from 1(a) by {worst_i:.3%} in I, {worst_y:.3%} in Y"
        raise ValueError(msg)
    return worst_i, worst_y


def _points(spec: PanelSpec, found: list[Reading]) -> list[Point]:
    """Return a dataset's points, numbered in order of intensity."""
    where = f"Fig. {spec.figure}{spec.letter}"
    return [
        Point(
            f"p{k:02d}",
            f"{r.intensity:.4e}",
            f"{r.yield_:.4e}",
            "",
            f"{where}: marker at ({r.x:.3f}, {r.y:.3f}) pt",
        )
        for k, r in enumerate(found, start=1)
    ]


def extract() -> dict[str, list[Point]]:
    """Return Figs. 1 to 4's measured ion yields, one file per series."""
    legends()
    checks()
    fig9_check()
    found = readings_cache()
    out = {}
    for spec in PANELS:
        for series in spec.series:
            key = dataset_id(spec, series)
            out[f"{spec.gas}/{SOURCE}__{key}.csv"] = _points(spec, found[key])
    return out


def main() -> None:
    """Write the data files, and report Fig. 13(a) against the curves."""
    write(extract())
    worst_i, worst_y = fig9_check()
    print(f"Fig. 9(a) against Fig. 1(a) LP: {worst_i:.2%} in I, {worst_y:.2%} in Y")
    for stem, ratio, printed in ratio_check():
        print(
            f"{stem}: Y_CP/Y_LP at I_BSI {ratio:.3f} (curves), {printed:.3f} (Fig. 13)"
        )


if __name__ == "__main__":
    main()
