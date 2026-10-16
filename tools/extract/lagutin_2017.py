# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Lagutin et al. (2017): CIPFCP two-photon ionisation cross sections of Xe 5p.

The paper prints no table of its cross sections (Table I holds geometric factors only)
and no supplement exists: the generalised two-photon ionisation cross sections (its
G2PICS) are given only in figures, vector drawings on the page (Origin; the one image
of the PDF, on p. 5, is Fig. 1), read from the drawing (data rule 8, TC-699). The
journal pages 063414-1 to -10 are the PDF pages 1 to 10.

* **Fig. 3** (p. 6, library figure ``fig_2941080626c2``): the total cross section,
  length form, linear polarisation, in three approximations (HF, PF, PFCP); top panel
  7.0 to 11.0 eV, bottom panel its enlargement from 10.3 to 10.8 eV.
* **Fig. 4** (p. 7, ``fig_8ed2e2069b66``): the total cross section, linear
  polarisation, PF and CIPFCP, each in the length and velocity forms, 6.3 to 11.5 eV.
* **Fig. 6** (p. 7, ``fig_b52931f6435b``): the present total cross section with (W)
  and without (WO) "nonlinear screening", correlation (IV) excluded; L'Huillier and
  Wendin's curves in the same plot are not read.
* **Fig. 7** (p. 8, ``fig_2ce602ac73f9``): the partial cross section of the
  5p3/2^5 ep3/2 (J = 0) channel, CIPF with and without the shake-up correlation (IV),
  length and velocity forms, on a linear axis.
* **Fig. 8** (p. 8, ``fig_edb8d8774747``): the CIPFCP total cross section, length form,
  for linear (top panel) and circular (bottom panel) polarisation; the other
  calculations drawn (Fink and Johnson, Gangopadhyay et al., McGuire) and McCown et
  al.'s point are not read as data.

How the drawing is read:

* the frames: four stroked lines of the tick width per plot, found near the position
  :data:`PANELS` gives; the ticks are stroked two-vertex segments of the same width,
  perpendicular to the bottom (x) or left (y) edge with one end on it, short (about
  1.6 to 1.9 pt) and long (about 3.2 to 3.8 pt) at the labelled values.
* the axes: the labels are text objects; their values are written in :data:`PANELS`
  as printed and checked against the text objects of the page under the bottom and
  left of each frame (:func:`label_check`). The long ticks are fitted to the labels;
  each short tick is then given the nearest value of the axis's grid (the label
  interval divided by the number of short ticks between the first two labels, or 1 to
  9 times a decade on a log axis) and every tick, short and long, must lie on the
  scale fitted to all of them within :data:`TICK_TOLERANCE`. The top panel of Fig. 8
  prints no x labels: its ticks must coincide with the bottom panel's, whose labels
  it shares.
* the curves: a curve is the set of stroked subpaths of one width and colour that
  have a vertex inside the frame, and more there than inside any other plot of the
  page (the resonances of the bottom panels of Figs. 3 and 8 run up into the top
  panels' frames), less those lying wholly inside a box of
  :data:`PANELS` (the legend's samples) or of the curve (Fig. 7's hatched arrows,
  drawn in the stroke of its solid curve). Solid curves are
  long polylines; dashed ones are drawn dash by dash, each dash (and each dot of a
  dash-dot line) a two-vertex subpath lying on the plotted line. Where two curves share
  a width and a colour (Fig. 4's PF L and V and CIPFCP L and V, Fig. 6's W and WO,
  Fig. 7's CIPF shake-up excluded L and V) each is drawn in one block, in the order
  :data:`CURVES` lists them, told apart by their pattern: a solid line, or the length
  of its full dashes. The split between the blocks is where the most subpaths match
  their own block's pattern (a dot, or a dash shortened at a bend, may match neither or
  the other's), at the largest jump on the page among equally good places
  (:func:`_blocks`); the misfits are counted and reported.
* the values: every vertex of the curve's subpaths that lies inside the frame (within
  :data:`EDGE`), less exact repeats, sorted by photon energy (in drawing order where
  two share one). Vertices beyond the frame (the resonances run far above it, and
  Fig. 8's circular curve below it near 10.2 eV) are in the PDF but clipped from the
  printed figure, and are not taken; the checks use them only to locate resonances.
  Each ``printed`` cell gives the vertex's position in page points.

**Checks** (:func:`checks`), against every number the paper prints for these curves
and against the experimental reference marks it draws; the script fails if one is
missed by more than its stated margin.
"""

import ctypes
import functools
import math
import re
from dataclasses import dataclass

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, pdf, texts
from _library import DATA, Point, write
from pypdfium2 import raw

from optmat.constants import AU_ENERGY, E_CHARGE, HBAR, C

DIRECTORY = "Lagutin_2017_RelativisticCorrelationPolarization"
SOURCE = "lagutin_2017"
EDGE = 0.05
"""A frame line, a tick's end on it, or a vertex inside it, within this, pt."""
FRAME_FIND = 1.0
"""The frame's lines lie within this of :data:`PANELS`' position, pt."""
TICK_MAX = 5.0
"""Ticks are shorter than this, pt (the long ones 3.2 to 3.8)."""
LONG = 2.5
"""Ticks longer than this are the long, labelled ones, pt (the short 1.6 to 1.9)."""
TICK_TOLERANCE = 0.1
"""Every tick lies on its fitted scale within this, pt: the ticks' spacing alternates
by about 0.05 pt (Fig. 8's x ticks are 4.714 and 4.769 pt apart)."""
SOLID = 6.0
"""A two-vertex subpath longer than this is a piece of a solid line, pt (the dashes
are 2.4 to 4.4 pt)."""
DASH = 0.15
"""A full dash matches its curve's dash length within this, pt."""
MATCH = 0.1
"""A label's text object lies within this of its long tick, pt, after the offset
common to the axis."""


@dataclass(frozen=True)
class Panel:
    """One plot: where it is, what its axes print, and what in it is not a curve."""

    name: str
    page: int
    frame: tuple[float, float, float, float]
    """The frame's left, bottom, right, top, pt, to within :data:`FRAME_FIND`."""
    tick_width: float
    x_labels: tuple[str, ...]
    """The photon energies the x labels print, eV (those of ``x_from`` if empty)."""
    y_labels: tuple[str, ...]
    """The values the y labels print, in the axis's unit (10^-49 cm^4 s)."""
    log: bool = True
    boxes: tuple[tuple[float, float, float, float], ...] = ()
    """Boxes (l, b, r, t, pt) whose wholly enclosed subpaths are not curves."""
    x_from: str | None = None
    """The panel whose x labels (and tick positions) this one shares."""


PANELS = {
    "fig3a": Panel(
        "Fig. 3, top",
        6,
        (348.74, 559.03, 549.88, 664.11),
        0.669,
        ("7.0", "7.5", "8.0", "8.5", "9.0", "9.5", "10.0", "10.5", "11.0"),
        ("0.01", "0.1", "1", "10", "100", "1000", "10000"),
        boxes=((352.0, 628.0, 381.0, 657.0),),
    ),
    "fig3b": Panel(
        "Fig. 3, bottom",
        6,
        (348.74, 417.19, 549.88, 522.27),
        0.669,
        ("10.3", "10.4", "10.5", "10.6", "10.7", "10.8"),
        ("0.1", "1", "10", "100", "1000", "10000"),
        boxes=((352.0, 485.0, 382.0, 511.0),),
    ),
    "fig4": Panel(
        "Fig. 4",
        7,
        (88.89, 553.55, 267.92, 662.64),
        0.595,
        (
            "6.5",
            "7.0",
            "7.5",
            "8.0",
            "8.5",
            "9.0",
            "9.5",
            "10.0",
            "10.5",
            "11.0",
            "11.5",
        ),
        ("0.01", "0.1", "1", "10", "100"),
        boxes=((98.0, 614.0, 116.0, 655.0),),
    ),
    "fig6": Panel(
        "Fig. 6",
        7,
        (342.62, 585.13, 556.5, 734.1),
        0.711,
        ("6.5", "7.0", "7.5", "8.0", "8.5", "9.0", "9.5", "10.0", "10.5"),
        ("0.01", "0.1", "1", "10", "100", "1000"),
        boxes=((347.0, 690.0, 368.0, 726.0),),
    ),
    "fig7": Panel(
        "Fig. 7",
        8,
        (81.87, 587.8, 293.24, 735.01),
        0.703,
        ("7.0", "7.5", "8.0", "8.5", "9.0"),
        ("0.0", "0.1", "0.2", "0.3"),
        log=False,
        boxes=((89.0, 665.0, 110.0, 700.0),),
    ),
    "fig8a": Panel(
        "Fig. 8, top (linear polarisation)",
        8,
        (342.44, 585.19, 556.28, 734.12),
        0.711,
        (),
        ("0.01", "0.1", "1", "10", "100", "1000"),
        boxes=((348.0, 688.0, 369.0, 726.0),),
        x_from="fig8b",
    ),
    "fig8b": Panel(
        "Fig. 8, bottom (circular polarisation)",
        8,
        (342.44, 425.13, 556.33, 574.06),
        0.711,
        ("6.5", "7.0", "7.5", "8.0", "8.5", "9.0", "9.5", "10.0", "10.5"),
        ("0.01", "0.1", "1", "10", "100", "1000"),
        boxes=((348.0, 538.0, 374.0, 563.0),),
    ),
}
UNIT = 1e-49
"""Every y axis prints its values in units of 10^-49 cm^4 s."""


@dataclass(frozen=True)
class Curve:
    """One drawn curve: its plot, stroke and pattern, and its data file."""

    panel: str
    width: float
    colour: tuple[int, int, int]
    pattern: float | None
    """The length of its full dashes, pt, or None for a solid line."""
    file: str
    boxes: tuple[tuple[float, float, float, float], ...] = ()
    """Boxes (l, b, r, t, pt): the subpaths of its stroke wholly inside are not it."""


def _file(key: str) -> str:
    return f"Xe/{SOURCE}__{key}.csv"


RED, BLUE, BLACK = (255, 0, 0), (0, 0, 255), (0, 0, 0)
CURVES = {
    "fig3a_hf": Curve("fig3a", 0.926, RED, None, _file("fig3a_hf")),
    "fig3a_pf": Curve("fig3a", 0.926, BLUE, None, _file("fig3a_pf")),
    "fig3a_pfcp": Curve("fig3a", 0.926, BLACK, None, _file("fig3a_pfcp")),
    "fig3b_hf": Curve("fig3b", 0.926, RED, None, _file("fig3b_hf")),
    "fig3b_pf": Curve("fig3b", 0.926, BLUE, None, _file("fig3b_pf")),
    "fig3b_pfcp": Curve("fig3b", 0.926, BLACK, None, _file("fig3b_pfcp")),
    # Fig. 4: the two red curves are drawn in one block each, PF L (dash-dot, dashes
    # 3.17 pt) first; the black, CIPFCP L (solid) first, then V (dashes 3.67 pt).
    "fig4_pf_length": Curve("fig4", 0.595, RED, 3.17, _file("fig4_pf_length")),
    "fig4_pf_velocity": Curve("fig4", 0.595, RED, 2.39, _file("fig4_pf_velocity")),
    "fig4_cipfcp_length": Curve(
        "fig4", 0.824, BLACK, None, _file("fig4_cipfcp_length")
    ),
    "fig4_cipfcp_velocity": Curve(
        "fig4", 0.824, BLACK, 3.67, _file("fig4_cipfcp_velocity")
    ),
    # Fig. 6: WO (dashes 4.39 pt) first, then W (solid).
    "fig6_present_wo": Curve("fig6", 0.985, BLACK, 4.39, _file("fig6_present_wo")),
    "fig6_present_w": Curve("fig6", 0.985, BLACK, None, _file("fig6_present_w")),
    # Fig. 7: shake-up excluded L (dash-dot, dashes 3.74 pt) first, then V (dash-dot-
    # dot, 2.82 pt); CIPF V dashed (0.973 pt), CIPF L solid (0.487 pt, the one curve
    # of that width once the hatched arrows are boxed out).
    "fig7_noshakeup_length": Curve(
        "fig7", 0.703, RED, 3.74, _file("fig7_noshakeup_length")
    ),
    "fig7_noshakeup_velocity": Curve(
        "fig7", 0.703, RED, 2.82, _file("fig7_noshakeup_velocity")
    ),
    "fig7_cipf_length": Curve(
        "fig7",
        0.487,
        BLACK,
        None,
        _file("fig7_cipf_length"),
        boxes=((200.0, 600.0, 230.0, 664.0),),  # the hatched arrows' outlines
    ),
    "fig7_cipf_velocity": Curve(
        "fig7", 0.973, BLACK, None, _file("fig7_cipf_velocity")
    ),
    "fig8_cipfcp_linear": Curve(
        "fig8a", 0.985, BLACK, None, _file("fig8_cipfcp_linear")
    ),
    "fig8_cipfcp_circular": Curve(
        "fig8b", 0.985, BLACK, None, _file("fig8_cipfcp_circular")
    ),
}


@dataclass(frozen=True)
class Stroke:
    """A stroked subpath (page points), its width, colour and drawing index."""

    index: int
    width: float
    colour: tuple[int, int, int]
    vertices: np.ndarray
    filled: bool = False
    curved: bool = False
    """Whether the path has a Bezier segment (glyphs, the legend's marks)."""

    @property
    def length(self) -> float:
        """The length of the subpath, pt."""
        return float(np.sum(np.hypot(*np.diff(self.vertices, axis=0).T)))


@functools.cache
def _strokes(page: int) -> tuple[Stroke, ...]:
    """Return a page's stroked or filled paths, split into subpaths at each move.

    Raises if a path lies in a form (not mapped here). A path with a Bezier segment is
    marked ``curved`` (its control points kept as vertices): no curve, frame or tick
    of the figures is one (checked where they are selected).
    """
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    found = []
    for index, obj in enumerate(document[page - 1].get_objects(max_depth=15)):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        fill, stroked = ctypes.c_int(), ctypes.c_int()
        raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill), ctypes.byref(stroked))
        if obj.container is not None:
            msg = f"p. {page}: a path inside a form, which is not mapped here"
            raise ValueError(msg)
        width = ctypes.c_float()
        raw.FPDFPageObj_GetStrokeWidth(obj.raw, ctypes.byref(width))
        channels = [ctypes.c_uint() for _ in range(4)]
        get = (
            raw.FPDFPageObj_GetStrokeColor
            if stroked.value
            else raw.FPDFPageObj_GetFillColor
        )
        get(obj.raw, *(ctypes.byref(c) for c in channels))
        colour = (channels[0].value, channels[1].value, channels[2].value)
        matrix = obj.get_matrix()
        subpaths: list[list[tuple[float, float]]] = []
        curved = False
        for k in range(raw.FPDFPath_CountSegments(obj.raw)):
            segment = raw.FPDFPath_GetPathSegment(obj.raw, k)
            kind = raw.FPDFPathSegment_GetType(segment)
            curved = curved or kind == raw.FPDF_SEGMENT_BEZIERTO
            x, y = ctypes.c_float(), ctypes.c_float()
            raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
            point = matrix.on_point(x.value, y.value)
            if kind == raw.FPDF_SEGMENT_MOVETO or not subpaths:
                subpaths.append([])
            subpaths[-1].append((float(point[0]), float(point[1])))
        found.extend(
            Stroke(
                index,
                round(float(width.value), 3),
                colour,
                np.array(s, dtype=float).reshape(-1, 2),
                filled=not stroked.value,
                curved=curved,
            )
            for s in subpaths
        )
    return tuple(found)


def _lines(page: int, width: float) -> list[Stroke]:
    """Return a page's stroked two-vertex subpaths of a width (frames and ticks)."""
    return [
        s
        for s in _strokes(page)
        if not s.filled
        and not s.curved
        and len(s.vertices) == 2
        and abs(s.width - width) < 1e-3
    ]


@functools.cache
def frame(key: str) -> tuple[float, float, float, float]:
    """Return a plot's frame (l, b, r, t, pt): its four lines of the tick width."""
    panel = PANELS[key]
    left, bottom, right, top = panel.frame
    found = {"bottom": [], "top": [], "left": [], "right": []}
    for line in _lines(panel.page, panel.tick_width):
        (x0, y0), (x1, y1) = line.vertices
        if (
            abs(y0 - y1) < EDGE
            and abs(min(x0, x1) - left) < FRAME_FIND
            and abs(max(x0, x1) - right) < FRAME_FIND
        ):
            for name, edge in (("bottom", bottom), ("top", top)):
                if abs(y0 - edge) < FRAME_FIND:
                    found[name].append(y0)
        if (
            abs(x0 - x1) < EDGE
            and abs(min(y0, y1) - bottom) < FRAME_FIND
            and abs(max(y0, y1) - top) < FRAME_FIND
        ):
            for name, edge in (("left", left), ("right", right)):
                if abs(x0 - edge) < FRAME_FIND:
                    found[name].append(x0)
    if any(len({round(v, 2) for v in values}) != 1 for values in found.values()):
        msg = f"{panel.name}: no single frame near {panel.frame}: {found}"
        raise ValueError(msg)
    return (found["left"][0], found["bottom"][0], found["right"][0], found["top"][0])


@dataclass(frozen=True)
class Ticks:
    """A plot's tick positions along one edge, pt: all, and the long ones."""

    every: tuple[float, ...]
    long: tuple[float, ...]


def _ticks(key: str) -> tuple[Ticks, Ticks]:
    """Return the x ticks (bottom edge) and y ticks (left edge) of a plot."""
    panel = PANELS[key]
    left, bottom, right, top = frame(key)
    x_ticks: list[tuple[float, float]] = []
    y_ticks: list[tuple[float, float]] = []
    for line in _lines(panel.page, panel.tick_width):
        (x0, y0), (x1, y1) = line.vertices
        length = line.length
        if length > TICK_MAX:
            continue
        if (
            abs(x0 - x1) < EDGE
            and left - EDGE <= x0 <= right + EDGE
            and min(abs(y0 - bottom), abs(y1 - bottom)) < EDGE
        ):
            x_ticks.append((x0, length))
        if (
            abs(y0 - y1) < EDGE
            and bottom - EDGE <= y0 <= top + EDGE
            and min(abs(x0 - left), abs(x1 - left)) < EDGE
        ):
            y_ticks.append((y0, length))
    out = []
    for ticks in (x_ticks, y_ticks):
        every = sorted({round(p, 3) for p, _ in ticks})
        long = sorted({round(p, 3) for p, length in ticks if length > LONG})
        out.append(Ticks(tuple(every), tuple(long)))
    return out[0], out[1]


def _grid(ticks: Ticks, labels: tuple[str, ...], log: bool, name: str) -> Axis:
    """Return an axis calibrated on every tick: the long ones at their labels.

    Each short tick takes the nearest value of the axis's grid on the scale fitted to
    the long ticks; every tick must then lie on the scale fitted to all of them within
    :data:`TICK_TOLERANCE`.
    """
    values = [float(v) for v in labels]
    if len(ticks.long) != len(values):
        msg = f"{name}: {len(ticks.long)} long ticks for {len(values)} labels"
        raise ValueError(msg)
    coarse = Axis(ticks.long, tuple(values), log=log, tolerance=TICK_TOLERANCE)
    coarse.check()
    if log:
        exponents = range(math.floor(math.log10(values[0])) - 1, 6)
        grid = np.array([k * 10.0**e for e in exponents for k in range(1, 10)])
        grid_y = np.log10(grid)
    else:
        between = [p for p in ticks.every if ticks.long[0] < p < ticks.long[1]]
        step = (values[1] - values[0]) / (len(between) + 1)
        grid = values[0] + step * np.arange(-50, 51)
        grid_y = grid
    slope, intercept = coarse._fit()
    assigned = []
    for position in ticks.every:
        nearest = int(np.argmin(np.abs(grid_y - (slope * position + intercept))))
        assigned.append(float(grid[nearest]))
    if len(set(assigned)) != len(assigned):
        msg = f"{name}: two ticks take one grid value"
        raise ValueError(msg)
    axis = Axis(ticks.every, tuple(assigned), log=log, tolerance=TICK_TOLERANCE)
    axis.check()
    return axis


@functools.cache
def axes(key: str) -> tuple[Axis, Axis]:
    """Return a plot's photon-energy (eV) and cross-section axes, checked on every tick.

    The cross section is in the axes' printed unit, 10^-49 cm^4 s.
    """
    panel = PANELS[key]
    x_ticks, y_ticks = _ticks(key)
    if panel.x_from is not None:
        other, _ = _ticks(panel.x_from)
        if len(other.every) != len(x_ticks.every) or any(
            abs(a - b) > TICK_TOLERANCE
            for a, b in zip(other.every, x_ticks.every, strict=True)
        ):
            msg = f"{panel.name}: the x ticks are not those of {panel.x_from}"
            raise ValueError(msg)
        x_labels = PANELS[panel.x_from].x_labels
    else:
        x_labels = panel.x_labels
    x_axis = _grid(x_ticks, x_labels, log=False, name=f"{panel.name} x")
    y_axis = _grid(y_ticks, panel.y_labels, log=panel.log, name=f"{panel.name} y")
    return x_axis, y_axis


def residuals(key: str) -> tuple[float, float]:
    """Return the worst tick residual of a plot's x and y axes, pt."""
    out = []
    for axis in axes(key):
        y = np.log10(axis.values) if axis.log else np.asarray(axis.values)
        slope, intercept = np.polyfit(axis.positions, y, 1)
        out.append(float(np.max(np.abs((y - intercept) / slope - axis.positions))))
    return out[0], out[1]


_NUMBER = re.compile(r"^-?\d+(\.\d+)?$")


def label_check(key: str) -> float:
    """Check a plot's printed labels against the page's text objects.

    The x labels are the numbers of the text below the frame (one object or several,
    read left to right); the y labels the numbers of the text objects just left of
    the frame, bottom to top, each within :data:`MATCH` of its long tick after the
    offset common to the axis (the glyphs' centres sit a fixed distance from the
    ticks). Returns that offset, pt.
    """
    panel = PANELS[key]
    left, bottom, right, top = frame(key)
    _, y_ticks = _ticks(key)
    found = texts(
        DIRECTORY, panel.page, (left - 20.0, bottom - 15.0, right + 20.0, top + 5.0)
    )
    if panel.x_labels:
        below = sorted((x, t) for t, x, y in found if bottom - 15.0 < y < bottom - 3.0)
        words = [w for _, t in below for w in t.split() if _NUMBER.match(w)]
        if tuple(words) != panel.x_labels:
            msg = f"{panel.name}: the x labels print {words}, not {panel.x_labels}"
            raise ValueError(msg)
    beside = sorted(
        (y, t)
        for t, x, y in found
        if left - 17.5 < x < left - 2.0 and y > bottom - 3.0 and _NUMBER.match(t)
    )
    if tuple(t for _, t in beside) != panel.y_labels:
        msg = f"{panel.name}: the y labels print {beside}, not {panel.y_labels}"
        raise ValueError(msg)
    offsets = np.array([y for y, _ in beside]) - np.array(y_ticks.long)
    y_offset = float(np.median(offsets))
    if np.max(np.abs(offsets - y_offset)) > MATCH:
        msg = f"{panel.name}: the y labels are not at their long ticks: {offsets}"
        raise ValueError(msg)
    return y_offset


def _inside(
    vertices: np.ndarray, box: tuple[float, float, float, float], margin: float = EDGE
) -> np.ndarray:
    """Return which vertices lie inside a box (l, b, r, t), within a margin, pt."""
    left, bottom, right, top = box
    return (
        (vertices[:, 0] >= left - margin)
        & (vertices[:, 0] <= right + margin)
        & (vertices[:, 1] >= bottom - margin)
        & (vertices[:, 1] <= top + margin)
    )


def _style(key: str) -> list[Stroke]:
    """Return a curve's stroke group: the subpaths of its width and colour in its plot.

    Those with a vertex inside the frame, less those wholly inside a box.
    """
    curve = CURVES[key]
    panel = PANELS[curve.panel]
    # A subpath is drawn in the plot holding most of its vertices: the resonances of
    # Fig. 8's bottom panel run up into the top panel's frame (and Fig. 3's likewise).
    others = [frame(k) for k, p in PANELS.items() if p.page == panel.page]
    box = frame(curve.panel)
    found = [
        s
        for s in _strokes(panel.page)
        if not s.filled
        and abs(s.width - curve.width) < 1e-3
        and s.colour == curve.colour
        and np.any(_inside(s.vertices, box))
        and np.sum(_inside(s.vertices, box))
        == max(np.sum(_inside(s.vertices, other)) for other in others)
        and not any(np.all(_inside(s.vertices, b)) for b in panel.boxes + curve.boxes)
    ]
    if any(s.curved for s in found):
        msg = f"{key}: a Bezier segment in the curve, which this script does not read"
        raise ValueError(msg)
    return found


def _matches(stroke: Stroke, pattern: float | None) -> bool:
    """Whether a subpath is unmistakably of a pattern: solid, or a full dash."""
    if pattern is None:
        return len(stroke.vertices) > 2 or stroke.length > SOLID
    return len(stroke.vertices) == 2 and abs(stroke.length - pattern) <= DASH


MISFIT = 0.1
"""At most this fraction of a block's subpaths may match the other curve's pattern."""


@functools.cache
def _blocks(key: str) -> tuple[tuple[str, ...], int, tuple[int, int]]:
    """Return the curves sharing a curve's stroke, where they split, and the misfits.

    At most two curves share a stroke here, drawn one block after the other in the
    order of :data:`CURVES`. The split is the index that gives the most subpaths
    matching their block's pattern (:func:`_matches`; a dash cut short at a bend or a
    piece's end may match neither, or the other's); among equally good ones, the one
    across which the drawing jumps farthest on the page (the end of one curve to the
    start of the next). Raises if more than :data:`MISFIT` of a block matches only the
    other pattern. Returns the sharing curves, the split, and each block's misfits.
    """
    curve = CURVES[key]
    sharing = tuple(
        k
        for k, c in CURVES.items()
        if (c.panel, c.width, c.colour) == (curve.panel, curve.width, curve.colour)
    )
    strokes = _style(key)
    if len(sharing) == 1:
        return sharing, len(strokes), (0, 0)
    if len(sharing) != 2:
        msg = f"{key}: {len(sharing)} curves share its stroke"
        raise ValueError(msg)
    first = np.array([_matches(s, CURVES[sharing[0]].pattern) for s in strokes])
    second = np.array([_matches(s, CURVES[sharing[1]].pattern) for s in strokes])
    only_first, only_second = first & ~second, second & ~first
    scores = [
        int(only_first[:k].sum() + only_second[k:].sum())
        for k in range(len(strokes) + 1)
    ]
    best = [k for k, score in enumerate(scores) if score == max(scores)]

    def jump(k: int) -> float:
        if k in (0, len(strokes)):
            return -1.0
        return float(np.hypot(*(strokes[k].vertices[0] - strokes[k - 1].vertices[-1])))

    split = max(best, key=jump)
    misfits = (int(only_second[:split].sum()), int(only_first[split:].sum()))
    if misfits[0] > MISFIT * split or misfits[1] > MISFIT * (len(strokes) - split):
        msg = f"{key}: the curves sharing its stroke do not split ({misfits} misfits)"
        raise ValueError(msg)
    return sharing, split, misfits


def _split(key: str) -> tuple[Stroke, ...]:
    """Return the subpaths of one curve, split from those sharing its stroke."""
    sharing, split, _ = _blocks(key)
    strokes = _style(key)
    return tuple(strokes[:split] if sharing.index(key) == 0 else strokes[split:])


@functools.cache
def curve(key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a curve: its page positions (pt), photon energies (eV), values (cm^4 s).

    Every vertex inside the frame (on its edge at most), less exact repeats, sorted by
    photon energy.
    """
    box = frame(CURVES[key].panel)
    seen: set[tuple[float, float]] = set()
    kept = []
    for stroke in _split(key):
        for vertex in stroke.vertices[_inside(stroke.vertices, box, 0.0)]:
            mark = (round(float(vertex[0]), 3), round(float(vertex[1]), 3))
            if mark not in seen:
                seen.add(mark)
                kept.append(vertex)
    line = np.array(kept, dtype=float).reshape(-1, 2)
    line = line[np.argsort(line[:, 0], kind="stable")]
    x_axis, y_axis = axes(CURVES[key].panel)
    energy = np.asarray(x_axis(line[:, 0]), dtype=float)
    value = np.asarray(y_axis(line[:, 1]), dtype=float) * UNIT
    return line, energy, value


def _points(key: str) -> list[Point]:
    """Return a curve as data points, one per vertex."""
    line, energy, value = curve(key)
    return [
        Point(
            f"v{k:04d}",
            f"{e:.5f}",
            f"{v:.4e}",
            "",
            f"vertex at ({x:.3f}, {y:.3f}) pt",
        )
        for k, ((x, y), e, v) in enumerate(zip(line, energy, value, strict=True), 1)
    ]


def at(key: str, energy: float) -> float:
    """Return a curve's cross section (cm^4 s) at a photon energy (eV).

    Interpolated in the logarithm of the value (linearly on Fig. 7's linear axis)
    between the vertices that bracket the energy.
    """
    _, energies, values = curve(key)
    if not energies[0] <= energy <= energies[-1]:
        msg = f"{key}: {energy} eV is outside the curve"
        raise ValueError(msg)
    if PANELS[CURVES[key].panel].log:
        return float(10 ** np.interp(energy, energies, np.log10(values)))
    return float(np.interp(energy, energies, values))


def resonance(key: str, low: float, high: float) -> float:
    """Return the photon energy (eV) of a curve's highest drawn point in a window.

    From every vertex of the curve's subpaths, those beyond the frame (clipped from
    the figure) included, mapped through the plot's calibration: the mean energy of
    the vertices within :data:`TOP` of the highest, which is the tip of a solid line's
    peak, or the middle of a dashed line's two flanks where its dashes stop above the
    frame.
    """
    x_axis, _ = axes(CURVES[key].panel)
    vertices = np.vstack([s.vertices for s in _split(key)])
    energies = np.asarray(x_axis(vertices[:, 0]), dtype=float)
    mask = (energies >= low) & (energies <= high)
    if not np.any(mask):
        msg = f"{key}: nothing drawn between {low} and {high} eV"
        raise ValueError(msg)
    heights = vertices[mask, 1]
    return float(np.mean(energies[mask][heights >= heights.max() - TOP]))


TOP = 1.0
"""The vertices within this of a peak's highest are its tip, pt (:func:`resonance`)."""


# The reference marks. NIST levels: the repository's copy of the NIST ASD levels of Xe I
# (nist_asd_2026), cm^-1. 1 cm^-1 = 100 h c / e eV.
EV_PER_WAVENUMBER = 100 * 2 * math.pi * HBAR * C / E_CHARGE
RYDBERG_EV = AU_ENERGY / 2 / E_CHARGE
NIST_LEVELS = DATA / "external" / "nist_asd" / "Xe_I_levels.tsv"


@functools.cache
def odd_j1_levels() -> dict[str, float]:
    """Return the odd-parity J = 1 levels of Xe I up to 11.5 eV, eV, by name."""
    out = {}
    for line in NIST_LEVELS.read_text(encoding="utf-8").splitlines()[1:]:
        cells = [c.strip('"') for c in line.split("\t")]
        configuration, term, j, level = cells[0], cells[1], cells[2], cells[4]
        if j != "1" or not term.endswith("*"):
            continue
        energy = float(level) * EV_PER_WAVENUMBER
        if energy <= 11.5:
            name = configuration.split(".")[-1] + " " + term + "1"
            core = "'" if "1/2" in configuration else ""
            out[name.replace(" ", core + " ", 1)] = energy
    return out


def _bars(key: str, bottom: float, top: float) -> list[float]:
    """Return the x (eV) of the coloured vertical bars of a plot's 'Exper.' box.

    The stroked subpaths, not black, all of whose vertices share one x, starting in
    the box (between ``bottom`` and ``top``, pt), through the plot's x calibration;
    one value per bar (a dashed bar's dashes merged).
    """
    panel = PANELS[key]
    left, _, right, _ = frame(key)
    xs = sorted(
        float(s.vertices[0, 0])
        for s in _strokes(panel.page)
        if not s.filled
        and s.colour != BLACK
        and np.ptp(s.vertices[:, 0]) < EDGE
        and left < s.vertices[0, 0] < right
        and bottom - EDGE <= float(np.min(s.vertices[:, 1])) <= top
    )
    merged: list[list[float]] = []
    for x in xs:
        if merged and x - merged[-1][-1] < 0.5:
            merged[-1].append(x)
        else:
            merged.append([x])
    x_axis, _ = axes(key)
    return [float(x_axis(float(np.mean(m)))) for m in merged]


def _dotted_lines(key: str, width: float) -> list[float]:
    """Return the x (eV) of a plot's dotted vertical lines: dots of one width merged."""
    panel = PANELS[key]
    box = frame(key)
    xs = sorted(
        round(float(s.vertices[0, 0]), 2)
        for s in _strokes(panel.page)
        if not s.filled
        and abs(s.width - width) < 1e-3
        and np.ptp(s.vertices[:, 0]) < EDGE
        and np.all(_inside(s.vertices, box))
    )
    merged: list[list[float]] = []
    for x in xs:
        if merged and x - merged[-1][-1] < 0.5:
            merged[-1].append(x)
        else:
            merged.append([x])
    x_axis, _ = axes(key)
    return [float(x_axis(float(np.mean(m)))) for m in merged if len(m) > 10]


def star(key: str) -> tuple[float, float]:
    """Return the photon energy (eV) and value (cm^4 s) of a plot's star marker.

    The filled ten-pointed (closed, eleven-vertex) path inside the frame and outside
    its boxes (the legend's star is boxed out): its centre, the mean of its ten
    corners.
    """
    panel = PANELS[key]
    box = frame(key)
    found = [
        s
        for s in _strokes(panel.page)
        if s.filled
        and len(s.vertices) == 11
        and np.all(_inside(s.vertices, box))
        and not any(np.all(_inside(s.vertices, b)) for b in panel.boxes)
    ]
    if len(found) != 1:
        msg = f"{panel.name}: {len(found)} stars"
        raise ValueError(msg)
    centre = found[0].vertices[:10].mean(axis=0)
    x_axis, y_axis = axes(key)
    return float(x_axis(centre[0])), float(y_axis(centre[1])) * UNIT


def _nearest_level(energy: float) -> tuple[str, float]:
    """Return the NIST level nearest a photon energy: its name and energy (eV)."""
    levels = odd_j1_levels()
    name = min(levels, key=lambda n: abs(levels[n] - energy))
    return name, levels[name]


# Margins. TEXT: numbers the text prints to two figures (7.7 %, 'a factor of 2').
# SHIFT: resonance shifts printed to 0.01 eV, located on a drawing whose energy grid
# near a resonance is not given (0.02 eV). MARK: a reference mark drawn at a NIST
# energy, 0.01 eV (half a point at 48 pt/eV). TWICE: one calculation drawn in two
# figures, off resonance, where 0.1 pt in energy does not matter (3 %).
TEXT = 0.15
SHIFT = 0.02
MARK = 0.01
TWICE = 0.03
MCCOWN_COEFFICIENT_CM4_PER_W = 4e-32
"""McCown et al.'s two-photon coefficient at 193 nm, cm^4 W^-1 (mccown_1982), times
the photon energy in J gives a cross section in cm^4 s."""
MCCOWN_NM = 193.0
LV_THRESHOLD = {"CIPFCP": 0.077, "PF": 0.80}
"""p. 6: |sigma_L - sigma_V| / mean at the threshold, 80 % (PF) and 7.7 % (CIPFCP)."""
SHIFTS = {
    "6s rel": ("fig3a_hf", (8.6, 9.0), "fig3a_pf", (8.6, 9.0), -0.11),
    "6s CP": ("fig3a_pf", (8.6, 9.0), "fig3a_pfcp", (8.2, 8.6), -0.36),
    "5d rel": ("fig3b_hf", (10.55, 10.645), "fig3b_pf", (10.55, 10.63), -0.02),
    "5d CP": ("fig3b_pf", (10.55, 10.63), "fig3b_pfcp", (10.48, 10.55), -0.08),
}
"""pp. 5-6: the shifts of the 5p^5 6s and 5d resonances by relativistic effects (HF to
PF) and by core polarisation (PF to PFCP): -0.11, -0.36, -0.02 and -0.08 eV. The
windows (eV) bracket each curve's resonance in Fig. 3."""


def checks() -> dict[str, object]:  # noqa: PLR0912 (one list of checks)
    """Return the checks of the reading against the text and the reference marks.

    Raises if one fails its margin (:data:`TEXT`, :data:`SHIFT`, :data:`MARK`,
    :data:`TWICE`).
    """
    out: dict[str, object] = {}
    failures = []
    for key, panel in PANELS.items():
        out[f"{panel.name} y label offset (pt)"] = label_check(key)
        out[f"{panel.name} worst tick residuals (pt)"] = residuals(key)
    # p. 6: the threshold L/V difference, PF 80 %, CIPFCP 7.7 %, at the curves' common
    # first energy (Fig. 4).
    for name, printed in LV_THRESHOLD.items():
        low = name.lower()
        start = max(
            curve(f"fig4_{low}_length")[1][0], curve(f"fig4_{low}_velocity")[1][0]
        )
        length, velocity = (
            at(f"fig4_{low}_length", start),
            at(f"fig4_{low}_velocity", start),
        )
        read = abs(length - velocity) / ((length + velocity) / 2)
        out[f"{name} L/V difference at {start:.3f} eV"] = (read, printed)
        if abs(read - printed) > TEXT * printed:
            failures.append(f"{name} L/V difference {read:.3f}, printed {printed}")
    # The CIPFCP curves start at the 5p1/2 two-photon threshold, 0.9875 Ry / 2 (p. 3,
    # the experimental ionisation energy used in the energy denominators).
    threshold = 0.9875 * RYDBERG_EV / 2
    for key in ("fig4_cipfcp_length", "fig4_cipfcp_velocity", "fig8_cipfcp_linear"):
        first = float(curve(key)[1][0])
        out[f"{key} first energy (eV)"] = (first, threshold)
        if abs(first - threshold) > MARK:
            failures.append(f"{key} starts at {first:.4f} eV, not {threshold:.4f}")
    # pp. 5-6: the resonance shifts in Fig. 3.
    for name, (before, first, after, second, printed) in SHIFTS.items():
        shift = resonance(after, *second) - resonance(before, *first)
        out[f"Fig. 3 {name} shift (eV)"] = (shift, printed)
        if abs(shift - printed) > SHIFT:
            failures.append(f"Fig. 3 {name} shift {shift:.3f} eV, printed {printed}")
    # p. 7: 'at a photon energy of 8 eV ... the respective factor equals only to 2'
    # (nonlinear screening, Fig. 6).
    factor = at("fig6_present_wo", 8.0) / at("fig6_present_w", 8.0)
    out["Fig. 6 WO / W at 8 eV"] = (factor, 2.0)
    if abs(factor / 2.0 - 1) > TEXT:
        failures.append(f"Fig. 6 WO / W at 8 eV is {factor:.3f}, printed 2")
    # McCown et al.'s point, drawn in Figs. 4 and 8: 4e-32 cm^4 W^-1 times the photon
    # energy at 193 nm.
    photon_ev = 2 * math.pi * HBAR * C / E_CHARGE / (MCCOWN_NM * 1e-9)
    expected = MCCOWN_COEFFICIENT_CM4_PER_W * photon_ev * E_CHARGE
    for key in ("fig4", "fig8a"):
        energy, value = star(key)
        out[f"{PANELS[key].name} McCown star (eV, cm^4 s)"] = (
            (energy, photon_ev),
            (value, expected),
        )
        if abs(energy - photon_ev) > 2 * MARK or abs(value / expected - 1) > TWICE:
            failures.append(f"{PANELS[key].name}: the star at {energy}, {value}")
    # The experimental resonance energies drawn: Figs. 3 and 4's 'Exper.' bars and
    # Fig. 8's dotted lines, each at a NIST J = 1 odd level.
    marks = {
        "Fig. 3 bars": _bars("fig3a", 693.85, 711.0),
        "Fig. 4 bars": _bars("fig4", 678.24, 693.8),
        "Fig. 8 dotted lines": _dotted_lines("fig8b", 0.492),
    }
    for name, energies in marks.items():
        named = []
        for energy in energies:
            level, nist = _nearest_level(energy)
            named.append((level, round(energy, 4), round(nist, 4)))
            if abs(energy - nist) > MARK:
                failures.append(f"{name}: a mark at {energy:.4f} eV, no NIST level")
        out[name] = named
    # The resonances of the CIPFCP calculation (Fig. 4, length form) against the
    # experimental 6s and 6s' energies.
    for label, window in (("6s", (8.2, 8.7)), ("6s'", (9.4, 9.9))):
        out[f"CIPFCP L {label} resonance (eV)"] = resonance(
            "fig4_cipfcp_length", *window
        )
    # The same calculation drawn twice, off resonance: Fig. 4's CIPFCP L and Fig. 8's
    # linear curve; Fig. 3's PF and Fig. 4's PF L (both relativistic LOPT, length
    # form); Fig. 3's top panel and its enlargement. At each vertex of the first in
    # the windows (eV), against the second interpolated there (the top panel of Fig.
    # 3 draws its curves with few vertices: 0.03 to 0.07 eV apart).
    pairs = (
        ("fig4_cipfcp_length", "fig8_cipfcp_linear", ((6.75, 8.2), (8.8, 9.4))),
        ("fig3a_pf", "fig4_pf_length", ((7.05, 8.5), (9.0, 9.4))),
        ("fig3a_pfcp", "fig3b_pfcp", ((10.31, 10.43), (10.72, 10.79))),
        ("fig3a_pf", "fig3b_pf", ((10.31, 10.43), (10.72, 10.79))),
        ("fig3a_hf", "fig3b_hf", ((10.31, 10.43), (10.75, 10.79))),
    )
    for first, second, windows in pairs:
        _, energies, values = curve(first)
        ratios = [
            v / at(second, e)
            for e, v in zip(energies, values, strict=True)
            if any(low <= e <= high for low, high in windows)
        ]
        worst = max(ratios, key=lambda r: abs(r - 1))
        out[f"{first} / {second}, {len(ratios)} vertices, worst"] = worst
        if abs(worst - 1) > TWICE:
            failures.append(f"{first} / {second}: {worst:.3f} off resonance")
    # Reported: core polarisation's effect at the threshold ('22% higher', p. 9; Fig.
    # 3, PFCP over PF), and the circular-to-linear ratio of Fig. 8.
    out["Fig. 3 PFCP / PF at 7.1 eV"] = at("fig3a_pfcp", 7.1) / at("fig3a_pf", 7.1)
    for energy in (6.8, 7.0, 7.5, 8.0, 9.0, 10.0):
        out[f"Fig. 8 circular / linear at {energy} eV"] = at(
            "fig8_cipfcp_circular", energy
        ) / at("fig8_cipfcp_linear", energy)
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return every curve, after the checks."""
    checks()
    return {c.file: _points(key) for key, c in CURVES.items()}


def main() -> None:
    """Write the data files, and print the checks."""
    write(extract())
    for key in CURVES:
        line, energies, values = curve(key)
        print(
            f"{key}: {len(line)} vertices, {energies[0]:.4f} to {energies[-1]:.4f} eV, "
            f"{values.min():.3e} to {values.max():.3e} cm4 s"
        )
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
