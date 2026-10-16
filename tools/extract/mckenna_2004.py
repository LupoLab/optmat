# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract McKenna and van der Hart (2004): R-matrix Floquet cross sections, Ne, Ar.

The paper prints no table: its generalised two- and three-photon ionisation cross
sections are given only in figures, all vector drawings on the page (no image, no
form; the PDF's one image is the download banner of page 1), read from the drawing
(data rule 8, TC-699). The PDF's page 1 is that banner: journal p. 457 is PDF page 2,
so the journal page is the PDF page plus 455.

* **Fig. 1** (PDF p. 6, journal p. 461, library figure ``fig_f6941eec6a0c``): Ne,
  two-photon, 10.8 to 21.6 eV.
* **Fig. 2** (p. 7/462, ``fig_4da3adbbfca0``): its enlargement from 19.5 to 21.1 eV.
* **Fig. 3** (p. 8/463, ``fig_0af720bedf15``): Ne, three-photon, 7.2 to 10.8 eV.
* **Fig. 5** (p. 10/465, ``fig_bdcfada5769c``): Ar, two-photon, 7.9 to 14.1 eV.
* **Fig. 7** (p. 12/467, ``fig_c9dbe2d04620``): Ar, two-photon, 14.25 to 14.65 eV.
* **Fig. 8** (p. 13/468, ``fig_714c65223ddb``): Ar, three-photon, 5.3 to 7.8 eV, and
  its inset, 7.12 to 7.14 eV on a linear axis in units of 10^-84 cm^6 s^2.

Figs. 4 and 6 are ionisation rates against the photon energy at one intensity (the
resonances' anatomy), which no quantity of the vocabulary holds, and are not read.

How the drawing is read:

* the frames: a stroked rectangle (one closed five-vertex subpath, 0.35 pt lines) per
  plot, found near the position :data:`FIGURES` gives; the ticks are stroked two-vertex
  subpaths of the same width, perpendicular to an edge with one end on it, short (about
  1.3 pt) at every tick and long (about 2.7 pt) at the labelled ones. The top and right
  edges repeat the bottom and left ones' ticks (checked, within :data:`MIRROR`).
* the axes: the labels are glyphs drawn as filled paths, not text; their values are
  written in :data:`FIGURES` as read on the rendered page. x linear: the long ticks are
  the labels, the short ones divide each interval between labels into equal steps
  (2 in every figure; checked to be a whole number). y logarithmic in decades: the long
  ticks the decades, the short ones 1, 2, ..., 9 times each decade; the inset's y is
  linear, every tick long and labelled. Every tick, short and long, lies on its fitted
  scale within :data:`TICK_TOLERANCE` (:meth:`_figure.Axis.check`); the worst residual
  of each axis is reported by :func:`checks`.
* the curves: the stroked paths drawn with the thick line (0.66 to 0.69 pt; ticks and
  frames are 0.35 pt, the dashed threshold lines and the label arrows 0.40 pt), inside
  the plot's frame (for Fig. 8, the inset's curve is the one inside the inset's frame,
  the main curve the others). Each curve is drawn as many polyline subpaths (no
  Bezier segment), split over several path objects, in no particular order or
  direction; :func:`_chain` joins them end to end into one polyline, each next piece
  the one whose end is nearest the polyline's end (within :data:`JOIN`, and unique),
  turned to continue it. Where two pieces overlap (one piece of Fig. 8 restarts 1.1 pt
  before the previous one's end), the later piece's vertices before the end are dropped.
  The polyline must then never step back in photon energy: the curve is a function of
  it.
* the values: every vertex of the polyline, in drawing order, less exact repeats of
  the previous one. The drawing's coordinates are quantised in steps of about 0.13 pt:
  on the flat parts the line is a staircase of horizontal and vertical steps, and in
  the narrowest resonances several vertices share one photon energy (the line runs up
  and down inside one 0.13 pt column). Rows are numbered in drawing order; each
  ``printed`` cell gives the vertex's position in page points.

**Checks** (:func:`checks`), against every number the paper prints for these curves
(pp. 460, 462, 465, 466, 468, 469) and between the figures that draw the same
calculation twice; the script fails if one is missed by more than its stated margin.
"""

import ctypes
import functools
import itertools
import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, pdf
from _library import Point, write
from pypdfium2 import raw

DIRECTORY = "McKenna_2003_MultiphotonIonizationCross"
SOURCE = "mckenna_2004"
TICK_WIDTH = 0.35
"""The stroke width of the frames and ticks, pt."""
CURVE_WIDTH = 0.55
"""Strokes wider than this are the curves (0.66 to 0.69 pt), pt."""
EDGE = 0.02
"""A tick's end lies on its frame edge within this, pt."""
TICK_MAX = 4.0
"""Ticks are shorter than this, pt (the long ones 2.7)."""
LONG = 2.0
"""Ticks longer than this are the long, labelled ones, pt."""
TICK_TOLERANCE = 0.13
"""Every tick lies on its fitted scale within this, pt: the drawing's quantum."""
MIRROR = 0.02
"""A top or right tick repeats a bottom or left one within this, pt."""
FRAME_FIND = 3.0
"""The frame's corners lie within this of :data:`FIGURES`' position, pt."""
JOIN = 0.2
"""A curve's next piece starts within this of what is drawn, pt (:func:`_gap`)."""
TAIL = 2.0
"""How far back in x a next piece may start on the drawn curve, pt."""
STEP = 1.0
"""The largest undrawn vertical step between two pieces, pt."""


@dataclass(frozen=True)
class Figure:
    """One plot: where it is and what its tick labels print."""

    name: str
    page: int
    frame: tuple[float, float, float, float]
    """The frame's left, bottom, right, top, pt, to within :data:`FRAME_FIND`."""
    x_labels: tuple[float, ...]
    """The photon energies the x labels print, eV."""
    y_labels: tuple[float, ...]
    """The values the y labels print (decades, or the inset's linear values)."""
    log: bool = True
    inside: str | None = None
    """For an inset, the key of the plot it is drawn inside."""


FIGURES = {
    "fig1": Figure(
        "Fig. 1",
        6,
        (188.04, 551.75, 415.52, 741.2),
        (10, 12, 14, 16, 18, 20, 22),
        (1e-52, 1e-51, 1e-50, 1e-49, 1e-48),
    ),
    "fig2": Figure(
        "Fig. 2",
        7,
        (188.57, 548.34, 419.64, 740.9),
        (19.50, 19.75, 20.00, 20.25, 20.50, 20.75, 21.00),
        (1e-52, 1e-51, 1e-50, 1e-49, 1e-48),
    ),
    "fig3": Figure(
        "Fig. 3",
        8,
        (187.87, 553.87, 412.11, 740.73),
        (7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0),
        (1e-87, 1e-86, 1e-85, 1e-84, 1e-83, 1e-82, 1e-81, 1e-80),
    ),
    "fig5": Figure(
        "Fig. 5",
        10,
        (188.82, 548.31, 419.64, 740.66),
        (7.7, 8.7, 9.7, 10.7, 11.7, 12.7, 13.7),
        (1e-51, 1e-50, 1e-49, 1e-48, 1e-47),
    ),
    "fig7": Figure(
        "Fig. 7",
        12,
        (187.44, 555.75, 410.03, 741.24),
        (14.25, 14.35, 14.45, 14.55, 14.65),
        (1e-51, 1e-50, 1e-49, 1e-48, 1e-47),
    ),
    "fig8": Figure(
        "Fig. 8",
        13,
        (187.95, 543.03, 414.41, 741.17),
        (5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0),
        (1e-84, 1e-83, 1e-82, 1e-81, 1e-80, 1e-79),
    ),
    "fig8inset": Figure(
        "Fig. 8, inset",
        13,
        (239.83, 658.71, 298.33, 714.04),
        (7.12, 7.13, 7.14),
        (5.0, 5.5, 6.0, 6.5, 7.0),
        log=False,
        inside="fig8",
    ),
}
INSET_SCALE = 1e-84
"""The inset's y axis prints its values in units of 10^-84 cm^6 s^2 ('x 10^-84')."""

FILES = {
    "fig1": "Ne/mckenna_2004__fig1_sigma2_Ne.csv",
    "fig2": "Ne/mckenna_2004__fig2_sigma2_Ne.csv",
    "fig3": "Ne/mckenna_2004__fig3_sigma3_Ne.csv",
    "fig5": "Ar/mckenna_2004__fig5_sigma2_Ar.csv",
    "fig7": "Ar/mckenna_2004__fig7_sigma2_Ar.csv",
    "fig8": "Ar/mckenna_2004__fig8_sigma3_Ar.csv",
    "fig8inset": "Ar/mckenna_2004__fig8_inset_sigma3_Ar.csv",
}


@dataclass(frozen=True)
class Stroke:
    """A stroked path object's subpaths (page points) and its stroke width."""

    width: float
    subpaths: tuple[np.ndarray, ...]

    def inside(self, box: tuple[float, float, float, float], margin: float) -> bool:
        """Whether every vertex lies within a box (l, b, r, t), within a margin."""
        left, bottom, right, top = box
        points = np.vstack(self.subpaths)
        return bool(
            np.all(points[:, 0] >= left - margin)
            and np.all(points[:, 0] <= right + margin)
            and np.all(points[:, 1] >= bottom - margin)
            and np.all(points[:, 1] <= top + margin)
        )


@functools.cache
def _strokes(page: int) -> tuple[Stroke, ...]:
    """Return a page's stroked paths, split into subpaths at each move.

    Raises if a path lies in a form (not mapped here) or has a Bezier segment (the
    figures are polylines; checked).
    """
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    found = []
    for obj in document[page - 1].get_objects(max_depth=15):
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
        matrix = obj.get_matrix()
        subpaths: list[list[tuple[float, float]]] = []
        for index in range(raw.FPDFPath_CountSegments(obj.raw)):
            segment = raw.FPDFPath_GetPathSegment(obj.raw, index)
            kind = raw.FPDFPathSegment_GetType(segment)
            if kind == raw.FPDF_SEGMENT_BEZIERTO:
                msg = f"p. {page}: a Bezier segment, which this script does not read"
                raise ValueError(msg)
            x, y = ctypes.c_float(), ctypes.c_float()
            raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
            point = matrix.on_point(x.value, y.value)
            if kind == raw.FPDF_SEGMENT_MOVETO or not subpaths:
                subpaths.append([])
            subpaths[-1].append((float(point[0]), float(point[1])))
        found.append(
            Stroke(
                round(float(width.value), 3),
                tuple(np.array(s, dtype=float).reshape(-1, 2) for s in subpaths),
            )
        )
    return tuple(found)


def _frame(figure: Figure) -> tuple[float, float, float, float]:
    """Return a plot's frame (l, b, r, t, pt): its stroked closed rectangle."""
    found = []
    for stroke in _strokes(figure.page):
        if abs(stroke.width - TICK_WIDTH) > 1e-3:
            continue
        for sub in stroke.subpaths:
            if len(sub) != 5 or not np.allclose(sub[0], sub[-1]):
                continue
            box = (
                float(sub[:, 0].min()),
                float(sub[:, 1].min()),
                float(sub[:, 0].max()),
                float(sub[:, 1].max()),
            )
            if all(
                abs(a - b) < FRAME_FIND for a, b in zip(box, figure.frame, strict=True)
            ):
                found.append(box)
    if len({tuple(round(v, 3) for v in box) for box in found}) != 1:
        msg = f"{figure.name}: {len(found)} frames near {figure.frame}"
        raise ValueError(msg)
    return found[0]


@dataclass(frozen=True)
class Ticks:
    """A plot's tick positions along one edge and its mirror, pt: all and long ones."""

    every: tuple[float, ...]
    long: tuple[float, ...]


def _ticks(figure: Figure) -> tuple[Ticks, Ticks]:
    """Return the x ticks (bottom edge) and y ticks (left edge) of a plot.

    Checks that the top and right edges repeat them within :data:`MIRROR`.
    """
    left, bottom, right, top = _frame(figure)
    edges: dict[str, list[tuple[float, float]]] = {
        "bottom": [],
        "top": [],
        "left": [],
        "right": [],
    }
    segments = [
        sub
        for stroke in _strokes(figure.page)
        if abs(stroke.width - TICK_WIDTH) <= 1e-3
        for sub in stroke.subpaths
        if len(sub) == 2
    ]
    for (x0, y0), (x1, y1) in segments:
        length = math.hypot(x1 - x0, y1 - y0)
        if length > TICK_MAX:
            continue
        if abs(x0 - x1) < EDGE and left - EDGE <= x0 <= right + EDGE:
            for name, edge in (("bottom", bottom), ("top", top)):
                if min(abs(y0 - edge), abs(y1 - edge)) < EDGE:
                    edges[name].append((x0, length))
        if abs(y0 - y1) < EDGE and bottom - EDGE <= y0 <= top + EDGE:
            for name, edge in (("left", left), ("right", right)):
                if min(abs(x0 - edge), abs(x1 - edge)) < EDGE:
                    edges[name].append((y0, length))
    out = []
    for near, far in (("bottom", "top"), ("left", "right")):
        every = sorted({round(p, 3) for p, length in edges[near] if length <= LONG})
        long = sorted({round(p, 3) for p, length in edges[near] if length > LONG})
        mirror = sorted({round(p, 3) for p, _ in edges[far]})
        if len(mirror) != len(set(every) | set(long)) or any(
            min(abs(m - p) for p in every + long) > MIRROR for m in mirror
        ):
            msg = f"{figure.name}: the {far} ticks do not repeat the {near} ones"
            raise ValueError(msg)
        out.append(Ticks(tuple(every), tuple(long)))
    return out[0], out[1]


def _linear_values(ticks: Ticks, labels: Sequence[float], name: str) -> list[float]:
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


def _log_values(ticks: Ticks, decades: Sequence[float], name: str) -> list[float]:
    """Return the values of every tick of a log axis: 1, 2, ..., 9 times each decade."""
    if len(ticks.long) != len(decades):
        msg = f"{name}: {len(ticks.long)} long ticks for {len(decades)} decades"
        raise ValueError(msg)
    values = [d * k for d in decades[:-1] for k in range(1, 10)] + [decades[-1]]
    if len(ticks.every) != len(values):
        msg = f"{name}: {len(ticks.every)} short ticks for {len(values)} values"
        raise ValueError(msg)
    return values


@functools.cache
def axes(key: str) -> tuple[Axis, Axis]:
    """Return a plot's photon-energy (eV) and cross-section axes, checked on every tick.

    The inset's y axis is in its printed unit, 10^-84 cm^6 s^2.
    """
    figure = FIGURES[key]
    x_ticks, y_ticks = _ticks(figure)
    x_values = _linear_values(x_ticks, figure.x_labels, figure.name)
    x_axis = Axis(x_ticks.every, tuple(x_values), tolerance=TICK_TOLERANCE)
    if figure.log:
        y_values = _log_values(y_ticks, figure.y_labels, figure.name)
    else:
        y_values = _linear_values(y_ticks, figure.y_labels, figure.name)
    y_axis = Axis(
        y_ticks.every, tuple(y_values), log=figure.log, tolerance=TICK_TOLERANCE
    )
    for axis in (x_axis, y_axis):
        axis.check()
        # The long ticks are among the ticks checked: each must be at its label.
    for position, label in zip(x_ticks.long, figure.x_labels, strict=True):
        if abs(float(x_axis(position)) - label) > 1e-6 * max(1.0, abs(label)) + (
            TICK_TOLERANCE * x_axis.resolution()
        ):
            msg = f"{figure.name}: the long x tick at {position} is not at {label}"
            raise ValueError(msg)
    return x_axis, y_axis


def residuals(key: str) -> tuple[float, float]:
    """Return the worst tick residual of a plot's x and y axes, pt."""
    out = []
    for axis in axes(key):
        slope, intercept = np.polyfit(
            axis.positions,
            np.log10(axis.values) if axis.log else np.asarray(axis.values),
            1,
        )
        y = np.log10(axis.values) if axis.log else np.asarray(axis.values)
        out.append(float(np.max(np.abs((y - intercept) / slope - axis.positions))))
    return out[0], out[1]


def _pieces(key: str) -> list[np.ndarray]:
    """Return the subpaths of a plot's curve, as drawn."""
    figure = FIGURES[key]
    frame = _frame(figure)
    insets = [
        _frame(f) for f in FIGURES.values() if f.inside == key and f.page == figure.page
    ]
    pieces = []
    for stroke in _strokes(figure.page):
        if stroke.width <= CURVE_WIDTH or not stroke.inside(frame, 0.01):
            continue
        if any(stroke.inside(box, 0.01) for box in insets):
            continue
        pieces.extend(stroke.subpaths)
    if not pieces:
        msg = f"{figure.name}: no curve"
        raise ValueError(msg)
    return pieces


def _to_segment(point: np.ndarray, start: np.ndarray, end: np.ndarray) -> float:
    """Return the distance from a point to a segment, pt."""
    along = end - start
    length = float(along @ along)
    t = 0.0 if length == 0 else float(np.clip((point - start) @ along / length, 0, 1))
    return float(np.hypot(*(point - start - t * along)))


def _gap(point: np.ndarray, line: np.ndarray) -> float:
    """Return how far a piece's end lies from continuing a polyline, pt.

    Its distance to the polyline's last :data:`TAIL` in x (a piece that starts on
    what is drawn, or retraces it), or, directly above or below the polyline's end,
    the vertical step between them (a step of the staircase left undrawn), if that is
    within :data:`STEP`.
    """
    end = line[-1]
    if abs(point[0] - end[0]) < EDGE and abs(point[1] - end[1]) <= STEP:
        return 0.0
    tail = line[line[:, 0] >= end[0] - TAIL]
    tail = line[-2:] if len(tail) < 2 else tail
    return min(_to_segment(point, a, b) for a, b in itertools.pairwise(tail))


def _chain(pieces: list[np.ndarray], name: str) -> np.ndarray:
    """Join a curve's pieces end to end into one polyline.

    Starts from the leftmost end of any piece; the next piece is the one with an end
    that continues the polyline (:func:`_gap` within :data:`JOIN`), which must be the
    only one, turned to start there; its vertices before the polyline's end in x are
    dropped (an overlap). Most pieces start exactly where the last ended; in Fig. 8 one
    retraces the last 1.1 pt of the previous one; at Fig. 5's threshold and at the
    tops of a few narrow peaks the next starts directly below the last one's end.
    """
    left = min(range(len(pieces)), key=lambda k: pieces[k][:, 0].min())
    first = pieces[left]
    line = first if first[0, 0] <= first[-1, 0] else first[::-1]
    rest = [p for k, p in enumerate(pieces) if k != left]
    while rest:
        gaps = sorted(
            (_gap(p[-1 if backwards else 0], line), k, backwards)
            for k, p in enumerate(rest)
            for backwards in (False, True)
        )
        gap, k, backwards = gaps[0]
        if gap > JOIN or (len(gaps) > 1 and gaps[1][0] <= JOIN and gaps[1][1] != k):
            msg = f"{name}: no single piece continues the curve at {line[-1]}"
            raise ValueError(msg)
        piece = rest.pop(k)
        piece = piece[::-1] if backwards else piece
        line = np.vstack([line, piece[piece[:, 0] >= line[-1, 0] - 1e-9]])
    keep = np.ones(len(line), dtype=bool)
    keep[1:] = np.any(np.diff(line, axis=0) != 0, axis=1)
    line = line[keep]
    if np.any(np.diff(line[:, 0]) < -1e-9):
        msg = f"{name}: the joined curve steps back in photon energy"
        raise ValueError(msg)
    return line


@functools.cache
def curve(key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a plot's curve: the page positions (pt), photon energies (eV), values.

    The values are in cm^2N s^(N-1), the inset's too (its printed unit times
    :data:`INSET_SCALE`).
    """
    figure = FIGURES[key]
    line = _chain(_pieces(key), figure.name)
    # A vertex on the frame's bottom or top edge is where the axis cuts the line (the
    # threshold's rise in Figs. 1 and 5): its value is only bounded, and is dropped.
    _, bottom, _, top = _frame(figure)
    line = line[
        (np.abs(line[:, 1] - bottom) > EDGE) & (np.abs(line[:, 1] - top) > EDGE)
    ]
    x_axis, y_axis = axes(key)
    energy = np.asarray(x_axis(line[:, 0]), dtype=float)
    value = np.asarray(y_axis(line[:, 1]), dtype=float)
    if not figure.log:
        value = value * INSET_SCALE
    return line, energy, value


def _points(key: str) -> list[Point]:
    """Return a plot's curve as data points, one per vertex."""
    line, energy, value = curve(key)
    return [
        Point(
            f"v{k:04d}",
            f"{e:.6f}",
            f"{v:.4e}",
            "",
            f"vertex {k} at ({x:.3f}, {y:.3f}) pt",
        )
        for k, ((x, y), e, v) in enumerate(zip(line, energy, value, strict=True), 1)
    ]


def at(key: str, energy: float) -> float:
    """Return a plot's cross section at a photon energy (eV), log-interpolated.

    Between the two vertices that bracket it, first in drawing order.
    """
    _, energies, values = curve(key)
    above = int(np.argmax(energies >= energy))
    if energies[above] < energy or above == 0:
        msg = f"{FIGURES[key].name}: {energy} eV is outside the curve"
        raise ValueError(msg)
    e0, e1 = energies[above - 1], energies[above]
    v0, v1 = values[above - 1], values[above]
    t = 0.0 if e1 == e0 else (energy - e0) / (e1 - e0)
    return float(10 ** ((1 - t) * math.log10(v0) + t * math.log10(v1)))


def _within(key: str, low: float, high: float) -> np.ndarray:
    """Return a plot's values between two photon energies (eV)."""
    _, energies, values = curve(key)
    return values[(energies >= low) & (energies <= high)]


def _local_minima(key: str, low: float, high: float) -> list[float]:
    """Return the photon energies (eV) of a plot's deepest point in a window."""
    _, energies, values = curve(key)
    mask = (energies >= low) & (energies <= high)
    return [float(energies[mask][int(np.argmin(values[mask]))])]


# How far the drawing may stray from the numbers the text prints (relative), and the
# figures drawing the same calculation from each other: the drawing's quantum (0.13
# pt, 0.3 to 0.8 % of a value on these log axes) and the text's own rounding (two or
# three figures), with room for the text's 'within 8 %' and 'within 10 %' being
# statements about ranges.
TEXT = 0.05
TWICE = 0.03
ENERGY = 0.01
"""Positions the text prints (eV) are matched within this: the text's two decimals and
the flat bottom of Fig. 2's 20.08 eV minimum (0.007 eV wide in the drawing)."""


def checks() -> dict[str, object]:
    """Return the checks of the reading against the text and between figures.

    Raises if one fails its margin (:data:`TEXT`, :data:`TWICE`, :data:`ENERGY`).
    """
    out: dict[str, object] = {}
    failures = []

    def compare(label: str, read: float, printed: float, margin: float) -> None:
        out[label] = (read, printed, read / printed - 1)
        if abs(read / printed - 1) > margin:
            failures.append(f"{label}: read {read:.4g}, printed {printed:.4g}")

    # p. 460: Ne, two-photon, 'Between 11 and 16 eV the cross sections remain fairly
    # constant, within 10% of 1.65 x 10^-51 cm4 s'.
    ne2 = _within("fig1", 11.0, 16.0)
    out["Ne 2-photon 11-16 eV range"] = (float(ne2.min()), float(ne2.max()))
    if ne2.min() < 1.65e-51 * (1 - 0.10 - TEXT) or ne2.max() > 1.65e-51 * (
        1 + 0.10 + TEXT
    ):
        failures.append("Fig. 1 strays beyond 10 % of 1.65e-51 between 11 and 16 eV")
    # p. 465: Ar, two-photon, 'For photon energies less than 11.5 eV ... within 8% of
    # 1.2 x 10^-50 cm4 s': from the top of the threshold's rise (8.0 eV) to 11.4 eV;
    # the 4s 1Po resonance's wing takes the drawn curve to 1.08 times 1.2e-50 at 11.43
    # eV and 1.17 times at 11.5 eV (reported).
    ar2 = _within("fig5", 8.0, 11.4)
    out["Ar 2-photon 8.0-11.4 eV range"] = (float(ar2.min()), float(ar2.max()))
    out["Ar 2-photon at 11.5 eV"] = at("fig5", 11.5)
    if ar2.min() < 1.2e-50 * (1 - 0.08 - TEXT) or ar2.max() > 1.2e-50 * (
        1 + 0.08 + TEXT
    ):
        failures.append("Fig. 5 strays beyond 8 % of 1.2e-50 below 11.5 eV")
    # p. 465: 'The next resonant state, at approximately 12.77 eV' (4s 1Se): the drawn
    # peak.
    _, e_ar2, v_ar2 = curve("fig5")
    window = (e_ar2 > 12.4) & (e_ar2 < 13.2)
    read = float(e_ar2[window][int(np.argmax(v_ar2[window]))])
    out["Ar 2-photon 4s 1Se peak (eV)"] = (read, 12.77)
    if abs(read - 12.77) > ENERGY:
        failures.append(f"Fig. 5's 4s 1Se peak at {read:.4f} eV, printed 12.77")
    # p. 466: 'The cross section of the resonance at 13.9 eV rises to 3.6 x 10^-48'.
    peak = _within("fig5", 13.7, 14.1)
    compare("Ar 2-photon 13.9 eV peak", float(peak.max()), 3.6e-48, TEXT)
    # p. 468: Ar, three-photon, 'Between 5.5 and 6.5 eV, the cross sections rise up to
    # 1.0 x 10^-82 cm6 s2 from a minimum value of 6.75 x 10^-84 cm6 s2'.
    ar3 = _within("fig8", 5.5, 6.5)
    compare("Ar 3-photon minimum 5.5-6.5 eV", float(ar3.min()), 6.75e-84, TEXT)
    # 'rise up to 1.0 x 10^-82' by 6.5 eV: where the drawn curve reaches 1.0e-82 (it
    # rises 3.3 decades per eV there, so a value at 6.5 eV would test the text's
    # rounding of the energy, not the reading).
    _, e_ar3, v_ar3 = curve("fig8")
    reach = float(e_ar3[int(np.argmax((e_ar3 > 6.0) & (v_ar3 >= 1.0e-82)))])
    out["Ar 3-photon reaches 1.0e-82 at (eV)"] = (reach, 6.5)
    if abs(reach - 6.5) > ENERGY:
        failures.append(f"Fig. 8 reaches 1.0e-82 at {reach:.4f} eV, not 6.5")
    # p. 462: Fig. 2's interference minima 'at 20.74 eV', 'at 20.08 eV', 'at 21.07 eV'.
    for printed, window in ((20.08, (20.03, 20.20)), (20.74, (20.68, 20.80))):
        (read,) = _local_minima("fig2", *window)
        out[f"Ne Fig. 2 minimum near {printed} eV"] = (read, printed)
        if abs(read - printed) > ENERGY:
            failures.append(f"Fig. 2's minimum at {read:.4f} eV, printed {printed}")
    # 21.07 eV: a dip in the 5d shoulder above the 6s peak; the curve's lowest point
    # between the 6s peak and the frame's end.
    (read,) = _local_minima("fig2", 21.03, 21.10)
    out["Ne Fig. 2 minimum near 21.07 eV"] = (read, 21.07)
    if abs(read - 21.07) > ENERGY:
        failures.append(f"Fig. 2's minimum at {read:.4f} eV, printed 21.07")
    # p. 469: the inset's 4d 1De and 6s 1Se features enhance the cross section 'by
    # nearly 20%' and 'by 2%': the peaks over the smooth background drawn under them
    # (reported, not a pass/fail check: the background is not printed).
    _, e_in, v_in = curve("fig8inset")
    peak_4d = float(v_in[(e_in > 7.125) & (e_in < 7.130)].max())
    dip = float(v_in[(e_in > 7.129) & (e_in < 7.133)].min())
    out["Ar inset 4d peak / dip after it"] = peak_4d / dip
    # The same calculation drawn twice: Fig. 1 against Fig. 2 (19.5-21.1 eV) and Fig. 8
    # against its inset (7.12-7.14 eV), off the resonances, where the drawn curve is
    # flat enough for 0.13 pt in energy not to matter.
    for big, small, energies in (
        ("fig1", "fig2", (19.55, 19.60, 20.20, 20.30, 20.40, 20.85, 20.90)),
        ("fig8", "fig8inset", (7.122, 7.124, 7.135, 7.138)),
    ):
        for energy in energies:
            compare(
                f"{FIGURES[big].name} over {FIGURES[small].name} at {energy} eV",
                at(big, energy),
                at(small, energy),
                TWICE,
            )
    if failures:
        raise ValueError("; ".join(failures))
    return out


# The photon energies the record's notes quote, eV: 193 nm (ArF) and 248.6 nm (KrF),
# 1239.842 eV nm over the wavelength.
QUOTED = {"Ar 3-photon at 193 nm": ("fig8", 1239.842 / 193.0)}


def quoted() -> dict[str, tuple[float, float]]:
    """Return the cross sections the record's notes quote: (photon energy, value)."""
    return {name: (energy, at(key, energy)) for name, (key, energy) in QUOTED.items()}


def extract() -> dict[str, list[Point]]:
    """Return every plot's curve, after the checks."""
    checks()
    return {FILES[key]: _points(key) for key in FILES}


def main() -> None:
    """Write the data files, and print the checks and the quoted values."""
    write(extract())
    for key in FILES:
        line, energies, values = curve(key)
        x_res, y_res = residuals(key)
        print(
            f"{FIGURES[key].name}: {len(line)} vertices, {energies[0]:.4f} to "
            f"{energies[-1]:.4f} eV, {values.min():.3e} to {values.max():.3e}; "
            f"ticks within {x_res:.3f} (x) and {y_res:.3f} (y) pt"
        )
    for label, result in checks().items():
        print(f"check {label}: {result}")
    for label, (energy, value) in quoted().items():
        print(f"{label} ({energy:.4f} eV): {value:.3e}")


if __name__ == "__main__":
    main()
