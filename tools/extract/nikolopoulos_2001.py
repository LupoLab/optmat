# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Nikolopoulos and Lambropoulos (2001): two-photon ionisation of He, 12-47 eV.

The paper prints no table: its two-photon single- and double-ionisation cross sections
are given in figures, vector drawings on the page (the PDF's one image is the download
banner of page 1), read from the drawing (data rule 8, TC-699). The journal page is the
PDF page plus 543 (p. 545 is PDF page 2).

* **Fig. 4** (PDF p. 14, journal p. 557, library figure ``fig_521c07c70805``): the
  partial two-photon cross sections to the S and D waves (He+, single ionisation),
  12.35 to 38.0 eV, two polylines; the markers of Bachau et al (1990) and Proulx and
  Shakeshaft (1993) are not taken.
* **Fig. 5** (p. 15/558, ``fig_47d9f579f9f5``): the total two-photon single (He+) and
  double (He++) ionisation cross sections, 42.7 to 46.5 eV, as markers joined by lines;
  the markers are the calculation's points. The (+) points, pseudoresonances the
  authors replaced by interpolation (p. 558), are not taken.
* **Fig. 3** (p. 14/557, ``fig_9a66eb16b035``): the D wave, 12.75 to 23.70 eV, 'this
  work, R = 40 a.u.', against Saenz and Lambropoulos (1999), whose dotted curve is not
  taken. Its curve is not Fig. 4's D wave drawn again: its energy grid differs, and
  Fig. 4's is lower by up to a quarter (reported by :func:`checks`).

How the drawing is read:

* the paths: stroked or filled, straight segments only for the frames, ticks and
  curves (the circles of Fig. 5 are Bezier outlines, whose bounding box is used); no
  form. Colours: black (35, 31, 32) for the frames, ticks, Fig. 4's S wave, Fig. 3's
  own curve and Fig. 5's He+; red (236, 42, 45) for Fig. 4's D wave and Fig. 5's He++.
* the frames: a stroked closed rectangle (one five-vertex subpath) near
  :data:`FIGURES`' position. The ticks: stroked two-vertex subpaths perpendicular to an
  edge with one end on it, short (about 1.53 pt) and long (3.06 pt, at the labels); the
  top and right edges repeat the bottom and left ones (checked). The labels are glyphs
  drawn as filled paths, not text: their values are written in :data:`FIGURES` as read
  on the rendered page. x linear: long ticks at the labels, short ones every
  :attr:`Figure.x_minor` eV; y logarithmic: long ticks at the decades, short ones at 2,
  ..., 9 times them. Each short tick is given its value from the long ticks' scale
  (within a fifth of a step), then every tick lies on the fitted scale within
  :data:`TICK_TOLERANCE` (:meth:`_figure.Axis.check`). The drawing's coordinates are
  quantised in steps of about 0.15 pt.
* the curves: the stroked polylines of a colour and width (:data:`CURVES`) with more
  than :data:`MIN_VERTICES` vertices inside the frame (the legend's samples have two).
  Every vertex, in drawing order, less those on the frame's bottom or top edge (where
  the frame clips the line at a deep minimum; the value there is only bounded). A
  polyline split into subpaths (Fig. 4's D wave, clipped at 21.2 eV) is joined in
  drawing order; it must never step back in photon energy.
* the markers (Fig. 5): the filled paths of the marker's colour and size inside the
  frame, outside the legend's box; each is drawn twice (filled, then outlined; the
  count is checked). The value is the marker's centre. The joining line runs through
  every marker's centre (checked, within :data:`ON_LINE`), and from the first marker to
  the frame's left edge (that vertex is not a value).

**Checks** (:func:`checks`), against every number the paper prints for these curves
(p. 558). The comparisons of Fig. 3 with Fig. 4, and with the values Sato et al. (2011)
attribute to this calculation, are reported, not checked (the record's notes).
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

DIRECTORY = "Nikolopoulos_2001_MultichannelTheoryTwo"
SOURCE = "nikolopoulos_2001"
BLACK = (35, 31, 32)
RED = (236, 42, 45)
FRAME_FIND = 0.5
"""The frame's corners lie within this of :data:`FIGURES`' position, pt."""
EDGE = 0.02
"""A tick's end lies on its frame edge within this, pt."""
SHORT = (1.4, 1.7)
"""The length of a short tick, pt (1.53 to 1.54 drawn)."""
LONG = (2.9, 3.2)
"""The length of a long tick, pt (3.06 to 3.09 drawn)."""
MIRROR = 0.02
"""A top or right tick repeats a bottom or left one within this, pt."""
TICK_TOLERANCE = 0.12
"""Every tick lies on its fitted scale within this, pt: the drawing's quantum."""
ASSIGN = 0.2
"""A short tick lies within this fraction of a step of the value it is given."""
MIN_VERTICES = 10
"""A curve has more vertices than this (a legend sample has two)."""
ON_LINE = 0.1
"""A marker's centre lies on its joining line, and each of the line's vertices on a
marker's centre, within this, pt (a circle's box, from its Bezier outline, is centred
0.08 pt left of the line's vertex)."""


@dataclass(frozen=True)
class Figure:
    """One plot: where it is and what its tick labels print."""

    name: str
    page: int
    frame: tuple[float, float, float, float]
    """The frame's left, bottom, right, top, pt, to within :data:`FRAME_FIND`."""
    x_labels: tuple[float, ...]
    """The photon energies the x labels print, eV (at the long ticks)."""
    x_minor: float
    """The step of the x ticks, eV."""
    y_decades: tuple[float, ...]
    """The decades the y labels print (at the long ticks), cm^4 s."""
    legend: tuple[float, float, float, float] | None = None
    """The legend's box inside the frame (l, b, r, t, pt): its markers are not data."""


FIGURES = {
    "fig3": Figure(
        "Fig. 3",
        14,
        (200.24, 523.46, 459.70, 740.75),
        (12.75, 17.75, 22.75),
        1.0,
        (1e-55, 1e-54, 1e-53, 1e-52, 1e-51, 1e-50, 1e-49, 1e-48),
    ),
    "fig4": Figure(
        "Fig. 4",
        14,
        (200.30, 225.89, 457.14, 439.92),
        (10, 20, 30),
        5.0,
        (1e-55, 1e-54, 1e-53, 1e-52, 1e-51, 1e-50, 1e-49, 1e-48, 1e-47),
        (349.52, 398.95, 461.27, 431.98),
    ),
    "fig5": Figure(
        "Fig. 5",
        15,
        (197.73, 522.41, 459.81, 740.81),
        (42.5, 44.5, 46.5),
        1.0,
        (1e-55, 1e-54, 1e-53, 1e-52),
        (383.52, 705.40, 444.36, 734.10),
    ),
}


@dataclass(frozen=True)
class Curve:
    """A drawn polyline: its plot, stroke colour and width (pt)."""

    figure: str
    colour: tuple[int, int, int]
    width: float


CURVES = {
    "fig4_S": Curve("fig4", BLACK, 0.459),
    "fig4_D": Curve("fig4", RED, 0.764),
    "fig3_D": Curve("fig3", BLACK, 0.463),
    "fig5_He+_line": Curve("fig5", BLACK, 0.780),
    "fig5_He++_line": Curve("fig5", RED, 0.780),
}


@dataclass(frozen=True)
class Marker:
    """A kind of marker: its plot, fill colour, size range (pt) and joining line."""

    figure: str
    colour: tuple[int, int, int]
    size: tuple[float, float]
    line: str
    bezier: bool
    """Whether its outline is drawn with Bezier segments (a circle) or lines."""


MARKERS = {
    "fig5_He+": Marker("fig5", BLACK, (4.5, 4.8), "fig5_He+_line", bezier=True),
    "fig5_He++": Marker("fig5", RED, (3.9, 4.2), "fig5_He++_line", bezier=False),
}

FILES = {
    "fig3_D": "He/nikolopoulos_2001__fig3_sigma2_D_He.csv",
    "fig4_S": "He/nikolopoulos_2001__fig4_sigma2_S_He.csv",
    "fig4_D": "He/nikolopoulos_2001__fig4_sigma2_D_He.csv",
    "fig5_He+": "He/nikolopoulos_2001__fig5_sigma2_single_He.csv",
    "fig5_He++": "He/nikolopoulos_2001__fig5_sigma2_double_He.csv",
}


@dataclass(frozen=True)
class Drawn:
    """A path object: paint, width (pt), subpaths (page points), whether curved."""

    stroke: tuple[int, int, int] | None
    fill: tuple[int, int, int] | None
    width: float
    subpaths: tuple[np.ndarray, ...]
    bezier: bool

    @property
    def box(self) -> tuple[float, float, float, float]:
        """The bounding box of the vertices (l, b, r, t, pt)."""
        points = np.vstack(self.subpaths)
        return (
            float(points[:, 0].min()),
            float(points[:, 1].min()),
            float(points[:, 0].max()),
            float(points[:, 1].max()),
        )


def _colour(obj: pdfium.PdfObject, stroke: bool) -> tuple[int, int, int] | None:
    """Return a path's stroke or fill colour, or None if it is not so painted."""
    fill_mode, stroked = ctypes.c_int(), ctypes.c_int()
    raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill_mode), ctypes.byref(stroked))
    if not (stroked.value if stroke else fill_mode.value):
        return None
    channels = [ctypes.c_uint() for _ in range(4)]
    get = raw.FPDFPageObj_GetStrokeColor if stroke else raw.FPDFPageObj_GetFillColor
    get(obj.raw, *(ctypes.byref(c) for c in channels))
    return channels[0].value, channels[1].value, channels[2].value


@functools.cache
def _drawn(page: int) -> tuple[Drawn, ...]:
    """Return a page's path objects, in drawing order. Raises on a path in a form."""
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    found = []
    for obj in document[page - 1].get_objects(max_depth=15):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        if obj.container is not None:
            msg = f"p. {page}: a path inside a form, which is not mapped here"
            raise ValueError(msg)
        matrix = obj.get_matrix()
        subpaths: list[list[tuple[float, float]]] = []
        bezier = False
        for index in range(raw.FPDFPath_CountSegments(obj.raw)):
            segment = raw.FPDFPath_GetPathSegment(obj.raw, index)
            kind = raw.FPDFPathSegment_GetType(segment)
            bezier = bezier or kind == raw.FPDF_SEGMENT_BEZIERTO
            x, y = ctypes.c_float(), ctypes.c_float()
            raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
            point = matrix.on_point(x.value, y.value)
            if kind == raw.FPDF_SEGMENT_MOVETO or not subpaths:
                subpaths.append([])
            subpaths[-1].append((float(point[0]), float(point[1])))
        width = ctypes.c_float()
        raw.FPDFPageObj_GetStrokeWidth(obj.raw, ctypes.byref(width))
        found.append(
            Drawn(
                _colour(obj, stroke=True),
                _colour(obj, stroke=False),
                round(float(width.value), 3),
                tuple(np.array(s, dtype=float).reshape(-1, 2) for s in subpaths),
                bezier,
            )
        )
    return tuple(found)


def _to_segment(point: np.ndarray, start: np.ndarray, end: np.ndarray) -> float:
    """Return the distance from a point to a segment, pt."""
    along = end - start
    length = float(along @ along)
    t = 0.0 if length == 0 else float(np.clip((point - start) @ along / length, 0, 1))
    return float(np.hypot(*(point - start - t * along)))


def _inside(
    point: np.ndarray, box: tuple[float, float, float, float], margin: float
) -> bool:
    """Whether a point lies within a box (l, b, r, t), within a margin."""
    left, bottom, right, top = box
    return bool(
        left - margin <= point[0] <= right + margin
        and bottom - margin <= point[1] <= top + margin
    )


@functools.cache
def frame(key: str) -> tuple[float, float, float, float]:
    """Return a plot's frame (l, b, r, t, pt): its one stroked closed rectangle."""
    figure = FIGURES[key]
    found = set()
    for drawn in _drawn(figure.page):
        if drawn.stroke is None or drawn.bezier:
            continue
        for sub in drawn.subpaths:
            if len(sub) != 5 or not np.allclose(sub[0], sub[-1]):
                continue
            box = (
                round(float(sub[:, 0].min()), 3),
                round(float(sub[:, 1].min()), 3),
                round(float(sub[:, 0].max()), 3),
                round(float(sub[:, 1].max()), 3),
            )
            if all(
                abs(a - b) < FRAME_FIND for a, b in zip(box, figure.frame, strict=True)
            ):
                found.add(box)
    if len(found) != 1:
        msg = f"{figure.name}: {len(found)} frames near {figure.frame}"
        raise ValueError(msg)
    return found.pop()


@dataclass(frozen=True)
class Ticks:
    """A plot's tick positions along one edge, pt: every one and the long ones."""

    every: tuple[float, ...]
    long: tuple[float, ...]


def _tick(
    sub: np.ndarray, box: tuple[float, float, float, float]
) -> tuple[str, float, bool] | None:
    """Return (edge, position, long) if a subpath is a tick of a frame, else None."""
    if len(sub) != 2:
        return None
    left, bottom, right, top = box
    (x0, y0), (x1, y1) = sub
    length = math.hypot(x1 - x0, y1 - y0)
    if not (SHORT[0] <= length <= SHORT[1] or LONG[0] <= length <= LONG[1]):
        return None
    is_long = length >= LONG[0]
    if abs(x0 - x1) < EDGE and left - EDGE <= x0 <= right + EDGE:
        for name, edge in (("bottom", bottom), ("top", top)):
            if min(abs(y0 - edge), abs(y1 - edge)) < EDGE:
                return name, round(float(x0), 3), is_long
    if abs(y0 - y1) < EDGE and bottom - EDGE <= y0 <= top + EDGE:
        for name, edge in (("left", left), ("right", right)):
            if min(abs(x0 - edge), abs(x1 - edge)) < EDGE:
                return name, round(float(y0), 3), is_long
    return None


def _ticks(key: str) -> tuple[Ticks, Ticks]:
    """Return the x ticks (bottom edge) and y ticks (left edge) of a plot.

    Checks that the top and right edges repeat them within :data:`MIRROR`.
    """
    figure = FIGURES[key]
    left, bottom, right, top = frame(key)
    edges: dict[str, list[tuple[float, bool]]] = {
        "bottom": [],
        "top": [],
        "left": [],
        "right": [],
    }
    for drawn in _drawn(figure.page):
        if drawn.stroke != BLACK or drawn.bezier:
            continue
        for sub in drawn.subpaths:
            found = _tick(sub, (left, bottom, right, top))
            if found is not None:
                edges[found[0]].append((found[1], found[2]))
    out = []
    for near, far in (("bottom", "top"), ("left", "right")):
        every = sorted({p for p, _ in edges[near]})
        long = sorted({p for p, is_long in edges[near] if is_long})
        mirror = sorted({p for p, _ in edges[far]})
        if len(mirror) != len(every) or any(
            abs(m - p) > MIRROR for m, p in zip(mirror, every, strict=True)
        ):
            msg = f"{figure.name}: the {far} ticks do not repeat the {near} ones"
            raise ValueError(msg)
        out.append(Ticks(tuple(every), tuple(long)))
    return out[0], out[1]


def _linear(ticks: Ticks, labels: Sequence[float], step: float, name: str) -> Axis:
    """Return a linear axis: the labelled long ticks, then every tick at its value.

    The ticks lie every ``step`` from the first label (eV).
    """
    if len(ticks.long) != len(labels):
        msg = f"{name}: {len(ticks.long)} long x ticks for {len(labels)} labels"
        raise ValueError(msg)
    rough = Axis(ticks.long, tuple(float(v) for v in labels))
    values = []
    for position in ticks.every:
        value = float(rough(position))
        nearest = labels[0] + round((value - labels[0]) / step) * step
        if abs(value - nearest) > ASSIGN * step:
            msg = f"{name}: the x tick at {position} is not on a {step} eV step"
            raise ValueError(msg)
        values.append(round(nearest, 9))
    if len(set(values)) != len(values):
        msg = f"{name}: two x ticks at one value"
        raise ValueError(msg)
    return Axis(ticks.every, tuple(values), tolerance=TICK_TOLERANCE)


def _log(ticks: Ticks, decades: Sequence[float], name: str) -> Axis:
    """Return a log axis: the decades' long ticks, then k x 10^n at every tick."""
    if len(ticks.long) != len(decades):
        msg = f"{name}: {len(ticks.long)} long y ticks for {len(decades)} decades"
        raise ValueError(msg)
    rough = Axis(ticks.long, tuple(float(v) for v in decades), log=True)
    values = []
    for position in ticks.every:
        value = math.log10(float(rough(position)))
        decade = math.floor(value + 1e-6)
        k = round(10 ** (value - decade))
        k, decade = (1, decade + 1) if k == 10 else (k, decade)
        if k < 1 or abs(value - decade - math.log10(k)) > ASSIGN * math.log10(
            (k + 1) / k
        ):
            msg = f"{name}: the y tick at {position} is not at k x 10^n"
            raise ValueError(msg)
        values.append(float(f"{k}e{decade}"))
    if len(set(values)) != len(values):
        msg = f"{name}: two y ticks at one value"
        raise ValueError(msg)
    return Axis(ticks.every, tuple(values), log=True, tolerance=TICK_TOLERANCE)


@functools.cache
def axes(key: str) -> tuple[Axis, Axis]:
    """Return a plot's photon-energy (eV) and cross-section (cm^4 s) axes, checked."""
    figure = FIGURES[key]
    x_ticks, y_ticks = _ticks(key)
    x_axis = _linear(x_ticks, figure.x_labels, figure.x_minor, figure.name)
    y_axis = _log(y_ticks, figure.y_decades, figure.name)
    x_axis.check()
    y_axis.check()
    return x_axis, y_axis


def residuals(key: str) -> tuple[float, float]:
    """Return the worst tick residual of a plot's x and y axes, pt."""
    out = []
    for axis in axes(key):
        y = np.log10(axis.values) if axis.log else np.asarray(axis.values)
        slope, intercept = np.polyfit(axis.positions, y, 1)
        out.append(float(np.max(np.abs((y - intercept) / slope - axis.positions))))
    return out[0], out[1]


@functools.cache
def polyline(name: str) -> np.ndarray:
    """Return a curve's vertices (pt) in drawing order, its subpaths joined."""
    curve = CURVES[name]
    figure = FIGURES[curve.figure]
    box = frame(curve.figure)
    found = [
        d
        for d in _drawn(figure.page)
        if d.stroke == curve.colour
        and abs(d.width - curve.width) < 1e-3
        and not d.bezier
        and d.fill is None
        and sum(len(s) for s in d.subpaths) > MIN_VERTICES
        and all(_inside(p, box, EDGE) for s in d.subpaths for p in s)
    ]
    if len(found) != 1:
        msg = f"{figure.name}: {len(found)} polylines for {name}"
        raise ValueError(msg)
    line = np.vstack(found[0].subpaths)
    keep = np.ones(len(line), dtype=bool)
    keep[1:] = np.any(np.diff(line, axis=0) != 0, axis=1)
    line = line[keep]
    if np.any(np.diff(line[:, 0]) < -EDGE):
        msg = f"{figure.name}: {name} steps back in photon energy"
        raise ValueError(msg)
    return line


@functools.cache
def curve(name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a curve's page positions (pt), photon energies (eV) and values (cm^4 s).

    Vertices on the frame's bottom or top edge are dropped.
    """
    key = CURVES[name].figure
    _, bottom, _, top = frame(key)
    line = polyline(name)
    line = line[
        (np.abs(line[:, 1] - bottom) > EDGE) & (np.abs(line[:, 1] - top) > EDGE)
    ]
    x_axis, y_axis = axes(key)
    return (
        line,
        np.asarray(x_axis(line[:, 0]), dtype=float),
        np.asarray(y_axis(line[:, 1]), dtype=float),
    )


@functools.cache
def markers(name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a marker set's centres (pt), photon energies (eV) and values (cm^4 s).

    Checks that each is drawn twice (filled, then outlined) and that each centre is a
    vertex of the joining line, and the line has no other vertex inside the frame.
    """
    marker = MARKERS[name]
    figure = FIGURES[marker.figure]
    box = frame(marker.figure)
    filled, outlined = [], []
    for d in _drawn(figure.page):
        if d.bezier != marker.bezier:
            continue
        left, bottom, right, top = d.box
        size = max(right - left, top - bottom)
        centre = np.array([(left + right) / 2, (bottom + top) / 2])
        if not marker.size[0] <= size <= marker.size[1]:
            continue
        if not _inside(centre, box, 0) or (
            figure.legend is not None and _inside(centre, figure.legend, 0)
        ):
            continue
        if d.fill == marker.colour and d.stroke is None:
            filled.append(centre)
        elif d.stroke == marker.colour and d.fill is None:
            outlined.append(centre)
    centres = np.array(sorted(filled, key=lambda c: c[0]))
    if len(outlined) != len(filled) or not np.allclose(
        centres, sorted(outlined, key=lambda c: c[0]), atol=1e-3
    ):
        msg = f"{figure.name}: {name}'s markers are not each filled and outlined"
        raise ValueError(msg)
    line = polyline(marker.line)
    inner = line[line[:, 0] > box[0] + EDGE]
    far = max(
        min(_to_segment(c, a, b) for a, b in itertools.pairwise(line)) for c in centres
    )
    unmarked = max(float(np.min(np.hypot(*(centres - v).T))) for v in inner)
    if far > ON_LINE or unmarked > ON_LINE:
        msg = f"{figure.name}: {name}'s markers are not on the vertices of its line"
        raise ValueError(msg)
    x_axis, y_axis = axes(marker.figure)
    return (
        centres,
        np.asarray(x_axis(centres[:, 0]), dtype=float),
        np.asarray(y_axis(centres[:, 1]), dtype=float),
    )


def _points(name: str) -> list[Point]:
    """Return a curve's vertices or a marker set's centres as data points."""
    if name in MARKERS:
        positions, energies, values = markers(name)
        what = "marker"
    else:
        positions, energies, values = curve(name)
        what = "vertex"
    return [
        Point(
            f"{what[0]}{k:03d}",
            f"{e:.5f}",
            f"{v:.4e}",
            "",
            f"{what} {k} at ({x:.3f}, {y:.3f}) pt",
        )
        for k, ((x, y), e, v) in enumerate(
            zip(positions, energies, values, strict=True), 1
        )
    ]


def at(name: str, energy: float) -> float:
    """Return a curve's value (cm^4 s) at a photon energy (eV), log-interpolated."""
    if name in MARKERS:
        _, energies, values = markers(name)
    else:
        _, energies, values = curve(name)
    above = int(np.argmax(energies >= energy))
    if energies[above] < energy or above == 0:
        msg = f"{name}: {energy} eV is outside the curve"
        raise ValueError(msg)
    e0, e1 = energies[above - 1], energies[above]
    v0, v1 = values[above - 1], values[above]
    t = 0.0 if e1 == e0 else (energy - e0) / (e1 - e0)
    return float(10 ** ((1 - t) * math.log10(v0) + t * math.log10(v1)))


# The margin. TEXT: the text's '~ 1.2 x 10^-53' is two figures, rounded or cut: the
# drawn value must lie from 1.15 to 1.30e-53, widened by the digitising uncertainty
# (1 %, the record's notes).
TEXT = (1.15e-53 * 0.99, 1.30e-53 * 1.01)

# The values Sato et al., J. Phys. B 44, 161001 (2011), Table 2, attribute to this
# calculation ('Theory (LOPT)', their ref. [12]): their wavelengths (nm) and cross
# sections (1e-52 cm4 s). A secondary quotation, compared here, never written as data.
SATO = ((61.4, 8.0), (58.4, 360.0), (56.0, 2.7))
HC = 1239.84198
"""h c in eV nm (CODATA 2018), for the photon energy of a wavelength."""


def checks() -> dict[str, object]:
    """Return the checks against the text, and the comparisons reported.

    Raises if one fails its margin.
    """
    out: dict[str, object] = {}
    failures = []

    def compare(label: str, read: float, printed: float, margin: float) -> None:
        out[label] = (read, printed, read / printed - 1)
        if abs(read / printed - 1) > margin:
            failures.append(f"{label}: read {read:.4g}, printed {printed:.4g}")

    # p. 558: 'for the two-photon He+ ionization cross section we obtain sigma_2(He+)
    # ~ 1.2 x 10^-53 cm4 s', at omega = 45 eV.
    he = at("fig5_He+", 45.0)
    out["He+ at 45 eV, Fig. 5"] = (he, 1.2e-53, he / 1.2e-53 - 1)
    if not TEXT[0] <= he <= TEXT[1]:
        failures.append(f"He+ at 45 eV: read {he:.4g}, printed ~1.2e-53")
    # p. 558: 'the He2+ cross section is about 10% of the He+ cross section': the
    # ratio of the markers, reported, and required between 4 and 20 %.
    ratio = markers("fig5_He++")[2] / markers("fig5_He+")[2]
    out["He++ / He+ (Fig. 5), min and max"] = (float(ratio.min()), float(ratio.max()))
    if ratio.min() < 0.04 or ratio.max() > 0.2:
        failures.append(f"He++ / He+ from {ratio.min():.3f} to {ratio.max():.3f}")
    # Fig. 3's own curve and Fig. 4's D wave, reported: not one calculation drawn
    # twice (different energy grids; Fig. 4 is lower by up to a quarter).
    _, energies, values = curve("fig3_D")
    ratios = [
        at("fig4_D", float(energy)) / value
        for energy, value in zip(energies, values, strict=True)
        if 13.0 <= energy <= 20.5 or 21.8 <= energy <= 22.8
    ]
    out["Fig. 4 D over Fig. 3 D, off resonance, min and max"] = (
        min(ratios),
        max(ratios),
    )
    # Reported: Sato et al. (2011)'s values for this calculation, against Fig. 3's D
    # wave and against Fig. 4's S + D.
    for wavelength, quoted in SATO:
        energy = HC / wavelength
        d_wave = at("fig3_D", energy)
        total = at("fig4_S", energy) + at("fig4_D", energy)
        out[f"Sato's at {wavelength} nm ({energy:.3f} eV): Fig. 3 D, Fig. 4 S + D"] = (
            quoted * 1e-52,
            d_wave / (quoted * 1e-52),
            total / (quoted * 1e-52),
        )
    # Reported: the 1s2p resonance peak in each figure (eV).
    for name in ("fig3_D", "fig4_D", "fig4_S"):
        _, energies, values = curve(name)
        window = (energies > 20.5) & (energies < 21.5)
        out[f"{name} peak near 21 eV"] = float(
            energies[window][int(np.argmax(values[window]))]
        )
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return every curve and marker set, after the checks."""
    checks()
    return {path: _points(name) for name, path in FILES.items()}


def main() -> None:
    """Write the data files, and print the readings and the checks."""
    write(extract())
    for name in FILES:
        if name in MARKERS:
            _, energies, values = markers(name)
            key = MARKERS[name].figure
        else:
            _, energies, values = curve(name)
            key = CURVES[name].figure
        x_res, y_res = residuals(key)
        print(
            f"{name}: {len(energies)} points, {energies[0]:.4f} to {energies[-1]:.4f} "
            f"eV, {values.min():.3e} to {values.max():.3e}; ticks within "
            f"{x_res:.3f} (x) and {y_res:.3f} (y) pt"
        )
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
