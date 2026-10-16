# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Petrov et al. (2016): two-photon 3p ionisation cross sections of argon.

The paper prints no table of cross sections and no supplement exists: its generalised
two-photon ionisation cross sections (G2PICS) are given only in figures, all vector
drawings (no image on any page), read from the drawing (data rule 8, TC-699). The PDF
pages are the journal's pages 033408-1 to -15.

* **Fig. 10** (p. 8, library figure ``fig_8cec5ab9f473``): the total 3p cross section,
  log scale, 8 to 15 eV: four curves, many-electron correlations (corr) and in addition
  core polarisation (corr + CP), each in the length (L) and velocity (V) form.
* **Fig. 6** (p. 7, ``fig_ef2a6676c813``): the partial cross section to the
  3p5 ep(1S) channel, linear scale: the same four curves.
* **Fig. 7** (p. 7, ``fig_6bb337bc395a``): the ep(1D) channel, length form: corr and
  corr + CP (the LOPT curve is not read).
* **Fig. 8** (p. 8, ``fig_b0d6534b02ab``): the ef(1D) channel, length form: corr and
  corr + CP (the LOPT curve is not read).
* **Fig. 9** (p. 8, ``fig_88024d0eaf1f``): Fig. 10's four curves on a linear scale, read
  for the checks only (:func:`checks`), not written.

How the drawing is read:

* the frames: four stroked edges (0.70 to 0.72 pt lines, the corners joined into
  three-vertex paths with a tick), found near the position :data:`FIGURES` gives; the
  ticks are stroked two-vertex segments of the same width, perpendicular to an edge
  with one end on it: outside the frame on the bottom and left edges, inside on the
  top and right, long (3.7 to 3.8 pt) at the labelled values, short (1.9 pt) between.
  The top and right edges repeat the bottom and left ones' ticks (checked).
* the axes: x linear in photon energy (eV), long ticks at the printed 8, 9, ..., 15,
  one short tick between; y linear (Figs. 6 to 9) with a short tick between labels, or
  logarithmic (Fig. 10) with long ticks at the decades 10^1 to 10^5 and short ones at
  2, ..., 9 times each, starting from 2 at the frame's bottom edge; the printed unit is
  10^-51 cm^4 s. The label values are written in :data:`FIGURES` as read on the
  rendered page; every tick, short and long, lies on its fitted scale within
  :data:`TICK_TOLERANCE`.
* the curves: the coloured strokes (red and blue, 0.69 to 0.72 pt) with at least half
  their vertices in the frame (a curve of the plot below may cross it), less the
  legend's samples (inside :attr:`Figure.legend`). A solid curve is one polyline (it
  runs far beyond the frame at the resonances, clipped when drawn); a dashed one has
  been flattened into its dashes, each a short stroke of two or three vertices lying on
  the plotted line (a dash that spans a vertex of the line is split there or bent
  there). The curves are drawn one after another, each starting at the frame's left
  edge (8 eV), in the legend's order; :func:`curves` splits the strokes there and
  checks each curve's colour and style against :attr:`Figure.curves` (solid; dashed;
  dash-dot: about one dot per 3.7 pt dash; dash-dot-dot: about two per 3.0 pt dash).
* the order: every stroke runs forward in photon energy (checked); the dashes of a
  curve come in no fixed order and are put in order by :func:`_ordered` (by the x
  where each starts, and within one 0.05 pt column of a resonance's flank in the
  direction they run). Between dashes the line is not drawn (gaps of 1.9 to 4.7 pt
  along it); a gap wider than :data:`GAP` times the curve's median gap is allowed only
  where the curve leaves the frame (checked). Strokes of no length are dropped.
* the values: every vertex and dash end inside the frame (vertices above the top or
  below the bottom edge are where the curve leaves the plot and are dropped), in
  order along the curve; the files hold them in cm^4 s. Rows are numbered in that
  order; each ``printed`` cell gives the point's position in page points.

**Checks** (:func:`checks`): the numbers the text prints for these curves (pp. 6-8;
the conclusions' 17 % at the threshold, p. 10, matches the ep(1S) partial's 1.167,
not the total's 1.13, and is reported, not checked), the partial cross sections
summed against the total (Eq. 19), Fig. 9 against Fig. 10 (the same calculation drawn
twice), and Fig. 10's vertical lines against the NIST levels of Ar I they mark; the
script fails if one is missed by more than its stated margin. :func:`mckenna_ratios`
reports the comparison with mckenna_2004's R-matrix Floquet total (not a check).
"""

import csv
import ctypes
import functools
import itertools
import math
from dataclasses import dataclass

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, pdf
from _library import DATA, Point, write
from pypdfium2 import raw

from optmat.constants import E_CHARGE, HBAR, C

DIRECTORY = "Petrov_2016_CorrelationPolarizationEffects"
SOURCE = "petrov_2016"
BLACK = (0, 0, 0)
RED = (255, 0, 0)
BLUE = (0, 0, 255)
WIDTH = (0.68, 0.73)
"""Frames, ticks and curves are stroked between these widths, pt."""
EDGE = 0.03
"""A tick's end lies on its frame edge within this, pt."""
TICK_MAX = 4.5
"""Ticks are shorter than this, pt (the long ones 3.7 to 3.8)."""
LONG = 3.0
"""Ticks longer than this are the long, labelled ones, pt (the short ones 1.9)."""
TICK_TOLERANCE = 0.06
"""Every tick lies on its fitted scale within this, pt."""
MIRROR = 0.03
"""A top or right tick repeats a bottom or left one within this, pt."""
FRAME_FIND = 3.0
"""The frame's edges lie within this of :data:`FIGURES`' position, pt."""
OUTSIDE = 0.05
"""A vertex further than this outside the frame is dropped, pt."""
GAP = 1.5
"""The widest undrawn gap between consecutive dashes inside the frame, as a multiple
of the curve's median gap (its dash pattern's: 1.9 to 2.0 pt for the dash-dot lines,
3.5 to 4.2 pt for the dashed ones)."""


@dataclass(frozen=True)
class Figure:
    """One plot: where it is, what its labels print and which curves it draws."""

    name: str
    page: int
    frame: tuple[float, float, float, float]
    """The frame's left, bottom, right, top, pt, to within :data:`FRAME_FIND`."""
    y_labels: tuple[float, ...]
    """The values the y labels print, 10^-51 cm^4 s (the decades for a log axis)."""
    legend: tuple[float, float, float, float]
    """A box (l, b, r, t, pt) holding the legend's line samples and nothing else."""
    curves: tuple[tuple[str, tuple[int, int, int], str], ...]
    """Each curve in drawing (and legend) order: its key, colour and style."""
    log: bool = False
    x_labels: tuple[float, ...] = (8, 9, 10, 11, 12, 13, 14, 15)


FIGURES = {
    "fig6": Figure(
        "Fig. 6",
        7,
        (339.13, 584.44, 552.95, 734.57),
        (0, 2, 4, 6, 8, 10, 12),
        (345.0, 688.0, 372.0, 727.0),
        (
            ("corr_L", RED, "dashdot"),
            ("corr_V", RED, "dashdotdot"),
            ("corrcp_L", BLUE, "solid"),
            ("corrcp_V", BLUE, "dash"),
        ),
    ),
    "fig7": Figure(
        "Fig. 7",
        7,
        (341.0, 145.7, 553.0, 295.0),
        (0, 2, 4, 6, 8),
        (347.0, 254.0, 374.0, 284.0),
        (("corr_L", RED, "dashdot"), ("corrcp_L", BLUE, "solid")),
    ),
    "fig8": Figure(
        "Fig. 8",
        8,
        (77.0, 586.0, 289.6, 734.7),
        (0, 2, 4, 6, 8, 10, 12),
        (80.0, 698.0, 108.0, 728.0),
        (("corr_L", RED, "dashdot"), ("corrcp_L", BLUE, "solid")),
    ),
    "fig9": Figure(
        "Fig. 9",
        8,
        (78.2, 157.0, 289.6, 304.4),
        (0, 5, 10, 15, 20),
        (82.0, 262.0, 110.0, 298.0),
        (
            ("corr_L", RED, "dashdot"),
            ("corr_V", RED, "dashdotdot"),
            ("corrcp_L", BLUE, "solid"),
            ("corrcp_V", BLUE, "dash"),
        ),
    ),
    "fig10": Figure(
        "Fig. 10",
        8,
        (345.69, 587.63, 552.83, 733.06),
        (10.0, 100.0, 1000.0, 10000.0, 100000.0),
        (350.0, 686.0, 377.0, 725.0),
        (
            ("corr_L", RED, "dashdot"),
            ("corr_V", RED, "dashdotdot"),
            ("corrcp_L", BLUE, "solid"),
            ("corrcp_V", BLUE, "dash"),
        ),
        log=True,
    ),
}
LOG_BOTTOM = 2.0
"""Fig. 10's lowest tick, on the frame's bottom edge: 2 x 10^-51 cm^4 s."""

FILES = {
    ("fig10", "corrcp_L"): "Ar/petrov_2016__fig10_corrcp_L_sigma2_Ar.csv",
    ("fig10", "corrcp_V"): "Ar/petrov_2016__fig10_corrcp_V_sigma2_Ar.csv",
    ("fig10", "corr_L"): "Ar/petrov_2016__fig10_corr_L_sigma2_Ar.csv",
    ("fig10", "corr_V"): "Ar/petrov_2016__fig10_corr_V_sigma2_Ar.csv",
    ("fig6", "corrcp_L"): "Ar/petrov_2016__fig6_corrcp_L_sigma2_ep1S_Ar.csv",
    ("fig6", "corrcp_V"): "Ar/petrov_2016__fig6_corrcp_V_sigma2_ep1S_Ar.csv",
    ("fig6", "corr_L"): "Ar/petrov_2016__fig6_corr_L_sigma2_ep1S_Ar.csv",
    ("fig6", "corr_V"): "Ar/petrov_2016__fig6_corr_V_sigma2_ep1S_Ar.csv",
    ("fig7", "corrcp_L"): "Ar/petrov_2016__fig7_corrcp_L_sigma2_ep1D_Ar.csv",
    ("fig7", "corr_L"): "Ar/petrov_2016__fig7_corr_L_sigma2_ep1D_Ar.csv",
    ("fig8", "corrcp_L"): "Ar/petrov_2016__fig8_corrcp_L_sigma2_ef1D_Ar.csv",
    ("fig8", "corr_L"): "Ar/petrov_2016__fig8_corr_L_sigma2_ef1D_Ar.csv",
}


@dataclass(frozen=True)
class Stroke:
    """A stroked path object: its drawing index, colour, width and vertices (pt)."""

    index: int
    colour: tuple[int, int, int]
    width: float
    vertices: np.ndarray

    def inside(self, box: tuple[float, float, float, float], margin: float) -> bool:
        """Whether every vertex lies within a box (l, b, r, t), within a margin."""
        left, bottom, right, top = box
        x, y = self.vertices[:, 0], self.vertices[:, 1]
        return bool(
            np.all(x >= left - margin)
            and np.all(x <= right + margin)
            and np.all(y >= bottom - margin)
            and np.all(y <= top + margin)
        )

    def touches(self, box: tuple[float, float, float, float]) -> bool:
        """Whether any vertex lies within a box (l, b, r, t)."""
        left, bottom, right, top = box
        x, y = self.vertices[:, 0], self.vertices[:, 1]
        return bool(np.any((x >= left) & (x <= right) & (y >= bottom) & (y <= top)))


@functools.cache
def _strokes(page: int) -> tuple[Stroke, ...]:
    """Return a page's stroked paths in drawing order, each as one polyline.

    Raises if a path lies in a form, has a Bezier segment or more than one subpath
    (none does; checked).
    """
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    found = []
    for index, obj in enumerate(document[page - 1].get_objects(max_depth=15)):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        fill, stroked = ctypes.c_int(), ctypes.c_int()
        raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill), ctypes.byref(stroked))
        if not stroked.value:
            continue
        if obj.container is not None:
            msg = f"p. {page}: a stroked path inside a form, which is not mapped here"
            raise ValueError(msg)
        width = ctypes.c_float()
        raw.FPDFPageObj_GetStrokeWidth(obj.raw, ctypes.byref(width))
        channels = [ctypes.c_uint() for _ in range(4)]
        raw.FPDFPageObj_GetStrokeColor(obj.raw, *(ctypes.byref(c) for c in channels))
        colour = (channels[0].value, channels[1].value, channels[2].value)
        matrix = obj.get_matrix()
        points = []
        for k in range(raw.FPDFPath_CountSegments(obj.raw)):
            segment = raw.FPDFPath_GetPathSegment(obj.raw, k)
            kind = raw.FPDFPathSegment_GetType(segment)
            if kind == raw.FPDF_SEGMENT_BEZIERTO:
                msg = f"p. {page}: a Bezier segment, which this script does not read"
                raise ValueError(msg)
            if kind == raw.FPDF_SEGMENT_MOVETO and k > 0:
                msg = f"p. {page}: a path of several subpaths, not expected here"
                raise ValueError(msg)
            x, y = ctypes.c_float(), ctypes.c_float()
            raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
            point = matrix.on_point(x.value, y.value)
            points.append((float(point[0]), float(point[1])))
        found.append(
            Stroke(
                index,
                colour,
                round(float(width.value), 3),
                np.array(points, dtype=float).reshape(-1, 2),
            )
        )
    return tuple(found)


def _framed(stroke: Stroke) -> bool:
    """Whether a stroke is drawn in the frames' and ticks' colour and width."""
    return stroke.colour == BLACK and WIDTH[0] <= stroke.width <= WIDTH[1]


def _segments(page: int) -> list[tuple[np.ndarray, np.ndarray]]:
    """Return every straight segment of the page's black frame-width strokes."""
    return [
        (a, b)
        for stroke in _strokes(page)
        if _framed(stroke)
        for a, b in zip(stroke.vertices[:-1], stroke.vertices[1:], strict=True)
    ]


@functools.cache
def _frame(key: str) -> tuple[float, float, float, float]:
    """Return a plot's frame (l, b, r, t, pt): its four long edges."""
    figure = FIGURES[key]
    left, bottom, right, top = figure.frame
    height, width = top - bottom, right - left
    vertical, horizontal = [], []
    for a, b in _segments(figure.page):
        # An edge spans the plot (another plot on the page may share its line).
        if (
            abs(a[0] - b[0]) < EDGE
            and abs(a[1] - b[1]) > 0.9 * height
            and min(a[1], b[1]) < bottom + FRAME_FIND
            and max(a[1], b[1]) > top - FRAME_FIND
        ):
            vertical.append(float(a[0]))
        if (
            abs(a[1] - b[1]) < EDGE
            and abs(a[0] - b[0]) > 0.9 * width
            and min(a[0], b[0]) < left + FRAME_FIND
            and max(a[0], b[0]) > right - FRAME_FIND
        ):
            horizontal.append(float(a[1]))
    out = []
    for edges, approx in (
        (vertical, left),
        (horizontal, bottom),
        (vertical, right),
        (horizontal, top),
    ):
        near = [e for e in edges if abs(e - approx) < FRAME_FIND]
        if not near or max(near) - min(near) > EDGE:
            msg = f"{figure.name}: {len(near)} frame edges near {approx}: {near}"
            raise ValueError(msg)
        out.append(float(np.mean(near)))
    return out[0], out[1], out[2], out[3]


@dataclass(frozen=True)
class Ticks:
    """A plot's tick positions along one edge, pt: all and long ones."""

    every: tuple[float, ...]
    long: tuple[float, ...]


def _ticks(key: str) -> tuple[Ticks, Ticks]:
    """Return the x ticks (bottom edge) and y ticks (left edge) of a plot.

    Checks that the top and right edges repeat them within :data:`MIRROR`.
    """
    figure = FIGURES[key]
    left, bottom, right, top = _frame(key)
    edges: dict[str, list[tuple[float, float]]] = {
        "bottom": [],
        "top": [],
        "left": [],
        "right": [],
    }
    for (x0, y0), (x1, y1) in _segments(figure.page):
        length = math.hypot(x1 - x0, y1 - y0)
        if length > TICK_MAX:
            continue
        if abs(x0 - x1) < EDGE and left - EDGE <= x0 <= right + EDGE:
            for name, edge in (("bottom", bottom), ("top", top)):
                if min(abs(y0 - edge), abs(y1 - edge)) < EDGE:
                    edges[name].append((float(x0), length))
        if abs(y0 - y1) < EDGE and bottom - EDGE <= y0 <= top + EDGE:
            for name, edge in (("left", left), ("right", right)):
                if min(abs(x0 - edge), abs(x1 - edge)) < EDGE:
                    edges[name].append((float(y0), length))
    out = []
    for near, far in (("bottom", "top"), ("left", "right")):
        every = sorted({round(p, 3) for p, _ in edges[near]})
        long = sorted({round(p, 3) for p, length in edges[near] if length > LONG})
        mirror = sorted({round(p, 3) for p, _ in edges[far]})
        if len(mirror) != len(every) or any(
            min(abs(m - p) for p in every) > MIRROR for m in mirror
        ):
            msg = f"{figure.name}: the {far} ticks do not repeat the {near} ones"
            raise ValueError(msg)
        out.append(Ticks(tuple(every), tuple(long)))
    return out[0], out[1]


def _linear_values(ticks: Ticks, labels: tuple[float, ...], name: str) -> list[float]:
    """Return the values of every tick of a linear axis, from its labelled long ones."""
    if len(ticks.long) != len(labels):
        msg = f"{name}: {len(ticks.long)} long ticks for {len(labels)} labels"
        raise ValueError(msg)
    steps = (len(ticks.every) - 1) / (len(ticks.long) - 1)
    if steps != int(steps) or ticks.every[0] != ticks.long[0]:
        msg = f"{name}: the short ticks do not divide the labels' intervals evenly"
        raise ValueError(msg)
    step = (labels[1] - labels[0]) / steps
    return [labels[0] + k * step for k in range(len(ticks.every))]


def _log_values(ticks: Ticks, decades: tuple[float, ...], name: str) -> list[float]:
    """Return the values of every tick of Fig. 10's log axis, from :data:`LOG_BOTTOM`.

    Raises unless the long ticks are the decades and the ticks number as the values.
    """
    if len(ticks.long) != len(decades):
        msg = f"{name}: {len(ticks.long)} long ticks for {len(decades)} decades"
        raise ValueError(msg)
    below = [decades[0] / 10 * k for k in range(int(LOG_BOTTOM), 10)]
    values = below + [d * k for d in decades[:-1] for k in range(1, 10)]
    values.append(decades[-1])
    if len(ticks.every) != len(values):
        msg = f"{name}: {len(ticks.every)} ticks for {len(values)} values"
        raise ValueError(msg)
    return values


@functools.cache
def axes(key: str) -> tuple[Axis, Axis]:
    """Return a plot's photon-energy (eV) and cross-section (10^-51 cm^4 s) axes.

    Every tick is checked on the fitted scale, and each long tick at its label.
    """
    figure = FIGURES[key]
    x_ticks, y_ticks = _ticks(key)
    x_axis = Axis(
        x_ticks.every,
        tuple(_linear_values(x_ticks, figure.x_labels, figure.name)),
        tolerance=TICK_TOLERANCE,
    )
    y_values = (
        _log_values(y_ticks, figure.y_labels, figure.name)
        if figure.log
        else _linear_values(y_ticks, figure.y_labels, figure.name)
    )
    y_axis = Axis(
        y_ticks.every, tuple(y_values), log=figure.log, tolerance=TICK_TOLERANCE
    )
    for axis in (x_axis, y_axis):
        axis.check()
    for ticks, labels, axis in (
        (x_ticks, figure.x_labels, x_axis),
        (y_ticks, figure.y_labels, y_axis),
    ):
        for position, label in zip(ticks.long, labels, strict=True):
            read = float(axis(position))
            # Within the tick tolerance, in the axis's (relative, if log) units.
            margin = (
                TICK_TOLERANCE * axis.resolution() * (abs(label) if axis.log else 1)
            )
            if abs(read - label) > margin + 1e-9:
                msg = f"{figure.name}: the long tick at {position} reads {read}"
                raise ValueError(msg)
    return x_axis, y_axis


def residuals(key: str) -> tuple[float, float]:
    """Return the worst tick residual of a plot's x and y axes, pt."""
    out = []
    for axis in axes(key):
        y = np.log10(axis.values) if axis.log else np.asarray(axis.values)
        slope, intercept = np.polyfit(axis.positions, y, 1)
        out.append(float(np.max(np.abs((y - intercept) / slope - axis.positions))))
    return out[0], out[1]


SOLID = 10
"""A stroke of more than this many vertices is a solid curve (the dashes have two or
three)."""


def _style(strokes: list[Stroke]) -> str:
    """Return a curve's line style from its strokes: solid, or the dash pattern.

    Solid: one polyline of more than :data:`SOLID` vertices. Dashed: the strokes'
    lengths along the line: the dash-dot pattern has long dashes of 3.7 pt (3.6 to
    4.1) and one dot (0.6 to 0.9 pt) each; the dash-dot-dot pattern dashes of 3.0 pt
    (2.5 to 3.3) and two dots each; the dashed pattern dashes of 3.0 pt and no dots. A
    dash split at a vertex of the drawn polyline counts in neither, so the ratios are
    approximate.
    """
    if len(strokes) == 1 and len(strokes[0].vertices) > SOLID:
        return "solid"
    if any(len(s.vertices) > SOLID for s in strokes):
        return "mixed"
    lengths = np.array(
        [float(np.sum(np.hypot(*np.diff(s.vertices, axis=0).T))) for s in strokes]
    )
    dots = int(np.sum((lengths > 0.6) & (lengths < 0.9)))
    long = int(np.sum((lengths > 3.6) & (lengths < 4.1)))
    short = int(np.sum((lengths > 2.5) & (lengths < 3.3)))
    if long > 2 * short and 0.8 * long < dots < 2.0 * long:
        return "dashdot"
    if short > 2 * long and 1.6 * short < dots < 3.5 * short:
        return "dashdotdot"
    if short > 2 * long and dots < 0.5 * short:
        return "dash"
    return f"unknown ({dots} dots, {long} long, {short} short dashes)"


def _in_frame(stroke: Stroke, frame: tuple[float, float, float, float]) -> bool:
    """Whether at least half a stroke's vertices lie in a frame.

    A curve of another plot on the page may cross this one's frame (the solid lines
    run far beyond their frames at the resonances, clipped when drawn); a dash crossing
    the frame's edge has one of its two ends inside.
    """
    left, bottom, right, top = frame
    x, y = stroke.vertices[:, 0], stroke.vertices[:, 1]
    inside = (
        (x >= left - OUTSIDE)
        & (x <= right + OUTSIDE)
        & (y >= bottom - OUTSIDE)
        & (y <= top + OUTSIDE)
    )
    return bool(np.mean(inside) >= 0.5)


COLUMN = 0.02
"""Strokes that start within this in x (pt) are in one column of the drawing, whose
quantum is about 0.05 pt (the coordinates differ in the third decimal within one)."""
EXIT = 4.0
"""A dash end this near the frame's top or bottom (pt) may be where the curve leaves
the plot: the next dash, beyond the edge, is clipped away (a dash is 3.0 or 3.7 pt)."""
BACK = 0.06
"""The ordered points may step back in x by up to this (pt): one quantum of the
drawing, where a dash ends one column on from where the next one, in its column,
starts."""


@functools.cache
def curves(key: str) -> dict[str, list[Stroke]]:
    """Return a plot's curves: each one's strokes, keyed as :attr:`Figure.curves`.

    The coloured strokes in the frame (:func:`_in_frame`), less the legend's, in
    drawing order. A curve starts with a stroke on the frame's left edge after one that
    is not, and a solid curve (a polyline of more than :data:`SOLID` vertices) is a
    curve of its own. Raises unless the curves are those :attr:`Figure.curves` lists,
    in its order, colours and styles, and the legend holds a sample of each colour.
    """
    figure = FIGURES[key]
    frame = _frame(key)
    left = frame[0]
    coloured = [
        s
        for s in _strokes(figure.page)
        if s.colour != BLACK and WIDTH[0] <= s.width <= WIDTH[1]
    ]
    legend = [s for s in coloured if s.inside(figure.legend, 0.0)]
    if {s.colour for s in legend} != {colour for _, colour, _ in figure.curves}:
        msg = f"{figure.name}: the legend's samples are not the curves' colours"
        raise ValueError(msg)
    inside = [
        s for s in coloured if _in_frame(s, frame) and not s.touches(figure.legend)
    ]
    groups: list[list[Stroke]] = []
    previous_on_edge, previous_solid = False, False
    for stroke in inside:
        on_edge = bool(stroke.vertices[:, 0].min() <= left + 0.1)
        solid = len(stroke.vertices) > SOLID
        if (on_edge and not previous_on_edge) or solid or previous_solid:
            groups.append([])
        if not groups:
            msg = f"{figure.name}: a curve does not start at the left edge"
            raise ValueError(msg)
        groups[-1].append(stroke)
        previous_on_edge, previous_solid = on_edge, solid
    found = [({s.colour for s in group}, _style(group)) for group in groups]
    expected = [({colour}, style) for _, colour, style in figure.curves]
    if found != expected:
        msg = f"{figure.name}: curves {found}, expected {expected}"
        raise ValueError(msg)
    return {
        name: group for (name, _, _), group in zip(figure.curves, groups, strict=True)
    }


def _ordered(segments: list[np.ndarray]) -> list[np.ndarray]:
    """Put a curve's strokes in order along it, each running forward (checked).

    By the photon energy where each starts; strokes starting within :data:`COLUMN` of
    the first in a column are in that column. The strokes of one column, on the
    near-vertical flanks of a resonance, are put in the direction they run: the rising
    ones upwards, the falling ones downwards; where one column holds both (the top of
    a peak, or the bottom of a dip, narrower than 0.05 pt), the group whose first
    stroke starts nearer the last point placed comes first (the rising flank of a
    peak, the falling flank of a dip).
    """
    columns: list[list[np.ndarray]] = []
    for v in sorted(segments, key=lambda v: float(v[0, 0])):
        if columns and float(v[0, 0]) - float(columns[-1][0][0, 0]) < COLUMN:
            columns[-1].append(v)
        else:
            columns.append([v])
    out: list[np.ndarray] = []
    for column in columns:
        rising = sorted(
            (v for v in column if v[-1, 1] >= v[0, 1]), key=lambda v: float(v[0, 1])
        )
        falling = sorted(
            (v for v in column if v[-1, 1] < v[0, 1]), key=lambda v: -float(v[0, 1])
        )
        if rising and falling and out:
            end = out[-1][-1]
            up = float(np.hypot(*(rising[0][0] - end)))
            down = float(np.hypot(*(falling[0][0] - end)))
            out.extend(rising + falling if up <= down else falling + rising)
        else:
            out.extend(rising + falling)
    return out


@functools.cache
def curve(key: str, name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return one curve: its points (pt), photon energies (eV), values (1e-51 cm4 s).

    Every vertex and dash end inside the frame, in order along the curve
    (:func:`_ordered`); see the module docstring.
    """
    figure = FIGURES[key]
    left, bottom, right, top = _frame(key)
    strokes = curves(key)[name]
    for stroke in strokes:
        # Each stroke runs forward along the curve, never back in photon energy: the
        # solid polyline as drawn, and each dash as the line it was cut from.
        if np.any(np.diff(stroke.vertices[:, 0]) < -1e-9):
            msg = f"{figure.name} {name}: a stroke steps back in photon energy"
            raise ValueError(msg)

    def within(point: np.ndarray) -> bool:
        return bool(
            left - OUTSIDE <= point[0] <= right + OUTSIDE
            and bottom - OUTSIDE <= point[1] <= top + OUTSIDE
        )

    def edge(point: np.ndarray) -> bool:
        # Outside the frame, or so near its top or bottom that the dash beyond it is
        # clipped away: where the curve leaves the plot.
        return not within(point) or min(abs(point[1] - top), abs(point[1] - bottom)) < (
            EXIT
        )

    # A stroke of no length (the dasher leaves a few) is a point of the line that
    # another stroke already holds; it is dropped, as its column order is undefined.
    segments = _ordered(
        [s.vertices for s in strokes if np.ptp(s.vertices, axis=0).max() > 0.01]
    )
    gaps = [
        float(np.hypot(*(after[0] - before[-1])))
        for before, after in itertools.pairwise(segments)
    ]
    drawn = [g for g in gaps if g > 0.01]
    widest = GAP * float(np.median(drawn)) if drawn else 0.0
    for gap, before, after in zip(gaps, segments[:-1], segments[1:], strict=True):
        if gap > widest and not (edge(before[-1]) and edge(after[0])):
            msg = f"{figure.name} {name}: a gap of {gap:.2f} pt at {after[0]}"
            raise ValueError(msg)
    line = np.vstack(segments)
    keep = np.ones(len(line), dtype=bool)
    keep[1:] = np.any(np.abs(np.diff(line, axis=0)) > 1e-6, axis=1)
    line = line[keep]
    line = line[np.array([within(point) for point in line], dtype=bool)]
    if np.any(np.diff(line[:, 0]) < -BACK):
        msg = f"{figure.name} {name}: the curve steps back in photon energy"
        raise ValueError(msg)
    x_axis, y_axis = axes(key)
    energy = np.asarray(x_axis(line[:, 0]), dtype=float)
    value = np.asarray(y_axis(line[:, 1]), dtype=float)
    return line, energy, value


EDGE_EV = 0.002
"""At the frame's left edge (8 eV) the first point is taken within this (eV): one
quantum of the drawing (0.05 pt at 30 pt/eV)."""


def at(key: str, name: str, energy: float) -> float:
    """Return a curve's value (10^-51 cm^4 s) at a photon energy (eV), interpolated.

    Linearly between the two points that bracket it (log-linearly on Fig. 10's axis).
    """
    _, energies, values = curve(key, name)
    if abs(energy - energies[0]) <= EDGE_EV:
        return float(values[0])
    above = int(np.argmax(energies >= energy))
    if energies[above] < energy or above == 0:
        msg = f"{FIGURES[key].name} {name}: {energy} eV is outside the curve"
        raise ValueError(msg)
    e0, e1 = energies[above - 1], energies[above]
    v0, v1 = values[above - 1], values[above]
    t = 0.0 if e1 == e0 else (energy - e0) / (e1 - e0)
    if FIGURES[key].log:
        return float(10 ** ((1 - t) * math.log10(v0) + t * math.log10(v1)))
    return float((1 - t) * v0 + t * v1)


def peak(key: str, name: str, low: float, high: float) -> tuple[float, float]:
    """Return a curve's highest point in a window: photon energy (eV) and value."""
    _, energies, values = curve(key, name)
    mask = (energies >= low) & (energies <= high)
    k = int(np.argmax(values[mask]))
    return float(energies[mask][k]), float(values[mask][k])


def resonance_lines() -> list[float]:
    """Return the photon energies (eV) of Fig. 10's vertical lines.

    The black lines, as long as the frame is high and crossing it, which mark the
    experimental 3p5_3/2 nl'[3/2]_1 energies (NIST ASD, the paper's Ref. 30).
    """
    left, bottom, right, top = _frame("fig10")
    x_axis, _ = axes("fig10")
    found = set()
    for (x0, y0), (x1, y1) in _segments(8):
        if (
            abs(x0 - x1) < EDGE
            and left + 1 < x0 < right - 1
            and min(y0, y1) <= bottom + EDGE
            and max(y0, y1) > top
        ):
            found.add(round(float(x_axis(x0)), 4))
    return sorted(found)


# The reference marks: the repository's copy of the NIST ASD levels of Ar I
# (nist_asd_2026), cm^-1; 1 cm^-1 = 100 h c / e eV.
NIST_LEVELS = DATA / "external" / "nist_asd" / "Ar_I_levels.tsv"
EV_PER_WAVENUMBER = 100 * 2 * math.pi * HBAR * C / E_CHARGE
LEVEL = 0.002
"""A vertical line is at its level within this (eV): 0.05 pt at 29.6 pt/eV."""


def nist_levels() -> dict[str, float]:
    """Return the Ar I 3p5(2P3/2) nl [3/2] J = 1 odd levels below 15 eV, eV, by nl."""
    out = {}
    for line in NIST_LEVELS.read_text(encoding="utf-8").splitlines()[1:]:
        cells = [c.strip('"') for c in line.split("\t")]
        configuration, term, j, level = cells[0], cells[1], cells[2], cells[4]
        if "(2P*<3/2>)" not in configuration or term != "2[3/2]*" or j != "1":
            continue
        energy = float(level) * EV_PER_WAVENUMBER
        if energy < FIGURES["fig10"].x_labels[-1]:
            out[configuration.split(".")[-1]] = energy
    return out


# How far the readings may stray from the numbers the text prints: the text's
# percentages are given to two figures and at an energy or over a range; the drawing's
# quantum (0.05 pt) is 0.4 % of a value on Fig. 10's axis (31 pt per decade) and 0.03
# to 0.1 units of 10^-51 on the linear ones (7 to 12 pt per unit).
PERCENT = 0.03
"""A percentage the text prints is matched within this (absolute, in the ratio)."""
TWICE = 0.03
"""Fig. 9 and Fig. 10, the same calculation, agree within this (relative)."""
SUM = 0.04
"""The partial cross sections sum to the total within this (relative): three readings
of linear plots against one of the log plot, each within about 1 %."""
OFF_RESONANCE = (8.0, 8.5, 9.0, 9.5, 10.0, 10.5)
"""Photon energies (eV) where every curve is flat: below the 4s resonance's wing."""


def checks() -> dict[str, object]:  # noqa: PLR0912 (one list of checks)
    """Return the checks of the reading against the text, between figures and sums.

    Raises if one fails its margin (:data:`PERCENT`, :data:`TWICE`, :data:`SUM`).
    """
    out: dict[str, object] = {}
    failures: list[str] = []

    def ratio(label: str, read: float, printed: float, margin: float) -> None:
        out[label] = (round(read, 4), printed)
        if abs(read - printed) > margin:
            failures.append(f"{label}: read {read:.3f}, printed {printed:.3f}")

    # p. 8: Fig. 9, core polarisation increases the total 'in length ... and velocity
    # ... forms by 18% and 15% respectively at omega = 9.5 eV'.
    for form, printed in (("L", 1.18), ("V", 1.15)):
        for key in ("fig9", "fig10"):
            read = at(key, f"corrcp_{form}", 9.5) / at(key, f"corr_{form}", 9.5)
            ratio(f"{key} corr+CP/corr {form} at 9.5 eV", read, printed, PERCENT)
    # Abstract: core polarisation 'increases the photoionization cross section at the
    # photon energy from 8 eV to 10 eV by approximately 15%' (reported over the range,
    # the length and velocity forms); conclusions: 'by approximately 17% at the
    # threshold' (the length form at 8 eV, the lowest energy drawn; reported).
    for form in ("L", "V"):
        out[f"fig9 corr+CP/corr {form}, 8-10 eV"] = [
            round(at("fig9", f"corrcp_{form}", e) / at("fig9", f"corr_{form}", e), 3)
            for e in (8.0, 9.0, 10.0)
        ]
    # p. 6, Fig. 4: with the correlations the length and velocity totals differ by
    # '8% at the 8-11-eV photon energy region': the corr curves of Figs. 9 and 10.
    for key in ("fig9", "fig10"):
        spread = [
            at(key, "corr_V", e) / at(key, "corr_L", e) for e in (8.0, 9.0, 10.0, 11.0)
        ]
        out[f"{key} corr V/L at 8, 9, 10, 11 eV"] = [round(s, 3) for s in spread]
        if max(abs(s - 1) for s in spread) > 0.08 + PERCENT:
            failures.append(f"{key}: corr V/L strays beyond 8 % between 8 and 11 eV")
    # p. 7, Fig. 6: core polarisation increases the ep(1S) partial 'by 17% for the
    # length form ... and by 11% for the velocity form' (no energy given: reported at
    # 8, 9 and 10 eV; the text's figure is matched at the threshold, 8 eV).
    for form, printed in (("L", 1.17), ("V", 1.11)):
        out[f"fig6 corr+CP/corr {form} at 8, 9, 10 eV"] = [
            round(at("fig6", f"corrcp_{form}", e) / at("fig6", f"corr_{form}", e), 3)
            for e in (8.0, 9.0, 10.0)
        ]
        read = at("fig6", f"corrcp_{form}", 8.0) / at("fig6", f"corr_{form}", 8.0)
        ratio(f"fig6 corr+CP/corr {form} at 8 eV", read, printed, PERCENT)
    # p. 8, Fig. 8: core polarisation increases the ef(1D) partial 'by 16% in the 8-11
    # eV photon energy range'.
    increase = [
        at("fig8", "corrcp_L", e) / at("fig8", "corr_L", e) for e in (8.0, 9.0, 10.0)
    ]
    out["fig8 corr+CP/corr L at 8, 9, 10 eV"] = [round(r, 3) for r in increase]
    if max(abs(r - 1.16) for r in increase) > PERCENT + 0.02:
        failures.append("Fig. 8: corr+CP/corr strays from 16 % between 8 and 10 eV")
    # p. 7, Fig. 7: core polarisation 'induces an additional shift of the theoretical
    # minimum by 0.2 eV ..., but the threshold G2PICSs remained almost unchanged'.
    minima = []
    for name in ("corr_L", "corrcp_L"):
        _, energies, values = curve("fig7", name)
        mask = (energies > 9.0) & (energies < 11.5)
        minima.append(float(energies[mask][int(np.argmin(values[mask]))]))
    out["fig7 ep(1D) minimum, corr and corr+CP (eV)"] = [round(m, 3) for m in minima]
    if abs((minima[0] - minima[1]) - 0.2) > 0.1:
        failures.append(f"Fig. 7: the minimum moves {minima[0] - minima[1]:.3f} eV")
    threshold = at("fig7", "corrcp_L", 8.0) / at("fig7", "corr_L", 8.0)
    out["fig7 corr+CP/corr L at 8 eV"] = round(threshold, 3)
    # Eq. (19): the total is the sum of the partials, ep(1S) + ep(1D) + ef(1D), in the
    # length form (the 1D partials are drawn only in it), for corr and corr + CP.
    for name in ("corr_L", "corrcp_L"):
        for energy in OFF_RESONANCE:
            parts = sum(at(k, name, energy) for k in ("fig6", "fig7", "fig8"))
            total = at("fig9", name, energy)
            out[f"sum of partials / Fig. 9 total, {name} at {energy} eV"] = round(
                parts / total, 4
            )
            if abs(parts / total - 1) > SUM:
                failures.append(
                    f"{name} at {energy} eV: partials {parts:.3f}, total {total:.3f}"
                )
    # The same calculation drawn twice: Fig. 9 (linear) against Fig. 10 (log).
    for name, _, _ in FIGURES["fig10"].curves:
        for energy in (*OFF_RESONANCE, 13.0):
            read = at("fig10", name, energy) / at("fig9", name, energy)
            out[f"Fig. 10 / Fig. 9, {name} at {energy} eV"] = round(read, 4)
            if abs(read - 1) > TWICE:
                failures.append(f"Fig. 10 over Fig. 9, {name} at {energy}: {read:.4f}")
    # Fig. 10's vertical lines, 'the experimental energies of the 3p5_3/2 nl'[3/2]_1
    # resonances [30]' (p. 8): each at a NIST level of Ar I (nist_asd_2026).
    lines = resonance_lines()
    levels = nist_levels()
    out["Fig. 10 vertical lines (eV)"] = lines
    out["NIST 3p5_3/2 nl[3/2]_1 levels (eV)"] = {
        k: round(v, 4) for k, v in levels.items()
    }
    if len(lines) != len(levels) or any(
        abs(line - level) > LEVEL
        for line, level in zip(lines, sorted(levels.values()), strict=True)
    ):
        failures.append(f"Fig. 10's lines {lines} are not the NIST levels {levels}")
    # The 4s resonance: drawn peaks against the experimental line (11.62 eV, p. 5).
    for name in ("corr_L", "corr_V", "corrcp_L", "corrcp_V"):
        out[f"Fig. 10 {name} 4s, highest drawn point (eV, 1e-51 cm4 s)"] = tuple(
            round(v, 4) for v in peak("fig10", name, 11.0, 12.5)
        )
    if failures:
        raise ValueError("; ".join(failures))
    return out


UNIT = 1e-51
"""The axes' printed unit, 10^-51 cm^4 s; the data files are in cm^4 s."""


def _points(key: str, name: str) -> list[Point]:
    """Return a curve as data points (cm^4 s), one per vertex or dash end."""
    line, energy, value = curve(key, name)
    return [
        Point(
            f"p{k:04d}",
            f"{e:.5f}",
            f"{v * UNIT:.5e}",
            "",
            f"({x:.3f}, {y:.3f}) pt",
        )
        for k, ((x, y), e, v) in enumerate(zip(line, energy, value, strict=True), 1)
    ]


def extract() -> dict[str, list[Point]]:
    """Return every taken curve, after the checks."""
    checks()
    return {path: _points(key, name) for (key, name), path in FILES.items()}


# A comparison the register quotes (not a check): McKenna and van der Hart's R-matrix
# Floquet total at 1e12 W/cm2 (mckenna_2004, Fig. 5), read from its committed file.
MCKENNA = "ionisation/Ar/mckenna_2004__fig5_sigma2_Ar.csv"


def mckenna_ratios() -> dict[float, tuple[float, float, float]]:
    """Return mckenna_2004's, corr+CP L and V totals at 8 to 11 eV, 1e-51 cm4 s."""
    with (DATA / "values" / MCKENNA).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    energies = np.array([float(r["x"]) for r in rows])
    values = np.array([float(r["value"]) for r in rows]) / 1e-51
    out = {}
    for energy in (8.0, 9.0, 10.0, 11.0):
        # Log-interpolated, as mckenna_2004's script does; the curve is smooth here.
        reference = float(10 ** np.interp(energy, energies, np.log10(values)))
        out[energy] = (
            reference,
            at("fig9", "corrcp_L", energy),
            at("fig9", "corrcp_V", energy),
        )
    return out


def main() -> None:
    """Write the data files, and print the readings' extent and the checks."""
    write(extract())
    for key, figure in FIGURES.items():
        x_res, y_res = residuals(key)
        for name, _, _ in figure.curves:
            line, energies, values = curve(key, name)
            print(
                f"{figure.name} {name}: {len(line)} points, {energies[0]:.4f} to "
                f"{energies[-1]:.4f} eV, {values.min():.4g} to {values.max():.4g}; "
                f"ticks within {x_res:.3f} (x) and {y_res:.3f} (y) pt"
            )
    for label, result in checks().items():
        print(f"check {label}: {result}")
    for energy, (reference, length, velocity) in mckenna_ratios().items():
        print(
            f"{energy} eV: mckenna_2004 {reference:.3f}, corr+CP L {length:.3f} "
            f"({length / reference:.3f}), V {velocity:.3f} ({velocity / reference:.3f})"
        )


if __name__ == "__main__":
    main()
