# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Mihelic and Horvat (2021): multiphoton ionisation cross sections of He.

The paper prints no table: its generalised two-, three- and four-photon partial
ionisation cross sections are given only in figures, drawn by matplotlib as vectors
(data rule 8, TC-699). The pages are numbered 043108-1 to -11, PDF pages 1 to 11.

* **Fig. 3** (PDF p. 6, library figure ``fig_5f3321672a12``): He, partial by the total
  orbital angular momentum L of the final state, linear polarisation. (a) two photons,
  L = 0, 2, 12.30 to 24.60 eV, below the one-photon threshold (BTI); (b) the same,
  24.60 to 43.54 eV, above it (ATI); (c), (d) three photons, L = 1, 3, 8.20 to 12.30
  and 12.30 to 43.54 eV; (e), (f) four photons, L = 0, 2, 4, 6.15 to 8.20 and 8.20 to
  43.54 eV (the frames). The dotted black curves of (a), (c), (e), Saenz and
  Lambropoulos's (saenz_1999), are not taken.
* **Fig. 4** (p. 7, ``fig_9c2c2686d9cd``): He, the total three-photon cross section
  from 25.25 to 27.35 eV, two calculations (channels with n_a <= 10, solid; n_a <= 5,
  dashed).
* **Fig. 5** (p. 7, ``fig_1b0120cd87bd``): He, two photons, the partial cross sections
  of the n_a l_a epsilon l L channels from 59.05 to 64.30 eV, L = 0 (total, 1s es,
  2s es, 2p ep) and L = 2 (total, 1s ed, 2s ed, 2p ep, 2p ef).

Fig. 2 (hydrogen) is drawn differently, its lines converted to filled outlines; it is
not read (the record's notes). Figs. 1 and 6 hold no cross section.

How the drawing is read:

* the frames: each plot's axes rectangle, drawn as four stroked black segments (the
  spines) and used as the clip path of its curves; both are checked against
  :data:`FIGURES`.
* the ticks: stroked black segments perpendicular to an edge with one end on it, major
  (about 2.1 to 2.3 pt) or minor (1.3 pt); top and bottom, left and right are merged
  (where both edges have ticks, they agree within :data:`MIRROR`). The labels are
  glyphs drawn as filled paths, not text: their values are written in :data:`FIGURES`
  as read on the rendered page (every label of each axis). x linear: the major ticks
  at :attr:`Figure.x_major`, the minor every :attr:`Figure.x_minor` eV. y logarithmic:
  the major ticks at :attr:`Figure.y_major` (every decade or every second one), the
  minor ones (Fig. 4 only) at k times a decade. Every tick lies on the fitted scale
  within :data:`TICK_TOLERANCE` (:meth:`_figure.Axis.check`).
* the curves: stroked polylines of the colour and width of :data:`CURVES`, clipped to
  the plot's frame, with more than :data:`MIN_VERTICES` vertices (a legend's sample has
  two). A curve drawn as several path objects or subpaths is joined in drawing order.
  matplotlib simplified the paths before writing them (vertices on a straight stretch
  are dropped), so the vertices are not the calculation's energy grid. Vertices
  outside the frame (beyond :data:`EDGE`) are not drawn (the clip hides them) and are
  not taken; the curve must never step back in photon energy.
* the values: every vertex inside the frame, in drawing order; each ``printed`` cell
  gives the vertex's position in page points.

**Checks** (:func:`checks`), against every number the paper prints for these curves:
the positions of the core-excited peak (40.8 eV) and of the thresholds where the
cross sections jump (32.7, 21.8 and 16.35 eV; p. 043108-7), and the dotted line of Fig.
4 at 1 a.u. (27.2 eV, p. 043108-7); and between Fig. 4's n_a <= 5 calculation and Fig.
3(d), which it enlarges. The script fails if one is missed.
"""

import ctypes
import functools
import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, pdf
from _library import Point, write
from pypdfium2 import raw

DIRECTORY = "Mihelic_2021_CalculationMultiphotonIonization"
SOURCE = "mihelic_2021"
BLACK = (0, 0, 0)
BLUE = (31, 119, 180)
ORANGE = (255, 127, 14)
GREEN = (44, 160, 44)
PURPLE = (116, 71, 157)
RED = (223, 39, 40)
GREY = (187, 187, 187)
FRAME_FIND = 0.01
"""The frame (spines and clip path) lies within this of :data:`FIGURES`' box, pt."""
EDGE = 0.02
"""A tick's end lies on its frame edge within this, pt."""
MINOR = (1.2, 1.4)
"""The length of a minor tick, pt (1.27 to 1.31 drawn)."""
MAJOR = (2.0, 2.4)
"""The length of a major tick, pt (2.12 to 2.29 drawn)."""
MIRROR = 0.01
"""A tick on the far edge repeats one on the near edge within this, pt."""
TICK_TOLERANCE = 0.01
"""Every tick lies on its fitted scale within this, pt."""
ASSIGN = 0.2
"""A minor tick lies within this fraction of a step of the value it is given."""
MIN_VERTICES = 3
"""A curve has more vertices than this (a legend's sample has two)."""
BACK = 0.005
"""A curve may step back in x by this much (the drawing's rounding), pt."""


@dataclass(frozen=True)
class Figure:
    """One plot: where it is, what its tick labels print, and its unit."""

    name: str
    page: int
    frame: tuple[float, float, float, float]
    """The frame's left, bottom, right, top, pt (the spines and the clip path)."""
    x_major: tuple[float, ...]
    """The photon energies at the major x ticks, eV (each labelled)."""
    x_minor: float | None
    """The step of the minor x ticks, eV, or None if there are none."""
    y_major: tuple[float, ...]
    """The values at the major y ticks (labelled each, or every other), axis unit."""
    unit: float
    """The axis's unit: 10^-52 cm^4 s, 10^-84 cm^6 s^2, 10^-118 cm^8 s^3."""


def _decades(first: int, last: int, step: int = 1) -> tuple[float, ...]:
    """Return 10^first, ..., 10^last in steps of ``step`` decades."""
    return tuple(float(f"1e{k}") for k in range(first, last + 1, step))


FIGURES = {
    "3a": Figure(
        "Fig. 3(a)",
        6,
        (93.073, 661.608, 237.198, 734.101),
        (12.5, 15.0, 17.5, 20.0, 22.5),
        0.5,
        _decades(-3, 3),
        1e-52,
    ),
    "3b": Figure(
        "Fig. 3(b)",
        6,
        (242.603, 661.608, 530.854, 734.101),
        (25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5),
        0.5,
        _decades(-1, 6),
        1e-52,
    ),
    "3c": Figure(
        "Fig. 3(c)",
        6,
        (93.073, 570.992, 237.198, 643.485),
        (9.0, 10.0, 11.0, 12.0),
        0.2,
        _decades(-3, 3),
        1e-84,
    ),
    "3d": Figure(
        "Fig. 3(d)",
        6,
        (242.603, 570.992, 530.854, 643.485),
        (15.0, 20.0, 25.0, 30.0, 35.0, 40.0),
        1.0,
        _decades(-4, 4, 2),
        1e-84,
    ),
    "3e": Figure(
        "Fig. 3(e)",
        6,
        (93.073, 480.38, 237.198, 552.87),
        (6.5, 7.0, 7.5, 8.0),
        0.1,
        _decades(-3, 3),
        1e-118,
    ),
    "3f": Figure(
        "Fig. 3(f)",
        6,
        (242.603, 480.38, 530.854, 552.87),
        (10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0),
        1.0,
        _decades(-4, 5),
        1e-118,
    ),
    "4": Figure(
        "Fig. 4",
        7,
        (343.56, 632.36, 556.56, 738.17),
        (25.5, 26.0, 26.5, 27.0),
        0.1,
        _decades(0, 2),
        1e-84,
    ),
    "5a": Figure(
        "Fig. 5, L = 0",
        7,
        (351.61, 199.32, 554.94, 266.65),
        (60.0, 61.0, 62.0, 63.0, 64.0),
        None,
        _decades(-3, 2),
        1e-52,
    ),
    "5b": Figure(
        "Fig. 5, L = 2",
        7,
        (351.61, 131.98, 554.94, 199.32),
        (60.0, 61.0, 62.0, 63.0, 64.0),
        None,
        _decades(-7, 1, 2),
        1e-52,
    ),
}


@dataclass(frozen=True)
class Curve:
    """A drawn curve: its plot, stroke colour and width (pt), and its dataset id."""

    figure: str
    colour: tuple[int, int, int]
    width: float
    dataset: str

    @property
    def file(self) -> str:
        """The data file, as the source record names it."""
        return f"He/{SOURCE}__{self.dataset}_He.csv"


CURVES = {
    "3a_L0": Curve("3a", BLUE, 0.981, "fig3a_sigma2_L0"),
    "3a_L2": Curve("3a", ORANGE, 0.654, "fig3a_sigma2_L2"),
    "3b_L0": Curve("3b", BLUE, 0.981, "fig3b_sigma2_L0"),
    "3b_L2": Curve("3b", ORANGE, 0.654, "fig3b_sigma2_L2"),
    "3c_L1": Curve("3c", GREEN, 0.981, "fig3c_sigma3_L1"),
    "3c_L3": Curve("3c", PURPLE, 0.654, "fig3c_sigma3_L3"),
    "3d_L1": Curve("3d", GREEN, 0.981, "fig3d_sigma3_L1"),
    "3d_L3": Curve("3d", PURPLE, 0.654, "fig3d_sigma3_L3"),
    "3e_L0": Curve("3e", BLUE, 0.981, "fig3e_sigma4_L0"),
    "3e_L2": Curve("3e", ORANGE, 0.654, "fig3e_sigma4_L2"),
    "3e_L4": Curve("3e", RED, 0.981, "fig3e_sigma4_L4"),
    "3f_L0": Curve("3f", BLUE, 0.981, "fig3f_sigma4_L0"),
    "3f_L2": Curve("3f", ORANGE, 0.654, "fig3f_sigma4_L2"),
    "3f_L4": Curve("3f", RED, 0.981, "fig3f_sigma4_L4"),
    "4_na10": Curve("4", BLACK, 0.763, "fig4_sigma3_na10"),
    "4_na5": Curve("4", BLUE, 0.763, "fig4_sigma3_na5"),
    "5a_total": Curve("5a", BLACK, 1.063, "fig5_sigma2_L0"),
    "5a_1ses": Curve("5a", BLUE, 0.729, "fig5_sigma2_L0_1ses"),
    "5a_2ses": Curve("5a", ORANGE, 0.729, "fig5_sigma2_L0_2ses"),
    "5a_2pep": Curve("5a", RED, 0.729, "fig5_sigma2_L0_2pep"),
    "5b_total": Curve("5b", BLACK, 1.063, "fig5_sigma2_L2"),
    "5b_1sed": Curve("5b", BLUE, 0.729, "fig5_sigma2_L2_1sed"),
    "5b_2sed": Curve("5b", ORANGE, 0.729, "fig5_sigma2_L2_2sed"),
    "5b_2pep": Curve("5b", RED, 0.729, "fig5_sigma2_L2_2pep"),
    "5b_2pef": Curve("5b", GREEN, 0.729, "fig5_sigma2_L2_2pef"),
}


@dataclass(frozen=True)
class Drawn:
    """A stroked path object: colour, width (pt), subpaths and clip box (page pt)."""

    stroke: tuple[int, int, int]
    width: float
    subpaths: tuple[np.ndarray, ...]
    clip: tuple[float, float, float, float] | None


def _clip(obj: pdfium.PdfObject) -> tuple[float, float, float, float] | None:
    """Return the box of an object's clip path, if it is one rectangle, else None."""
    clip = raw.FPDFPageObj_GetClipPath(obj.raw)
    if not clip or raw.FPDFClipPath_CountPaths(clip) != 1:
        return None
    points = []
    for index in range(raw.FPDFClipPath_CountPathSegments(clip, 0)):
        segment = raw.FPDFClipPath_GetPathSegment(clip, 0, index)
        x, y = ctypes.c_float(), ctypes.c_float()
        raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
        points.append((x.value, y.value))
    array = np.array(points, dtype=float)
    return (
        round(float(array[:, 0].min()), 3),
        round(float(array[:, 1].min()), 3),
        round(float(array[:, 0].max()), 3),
        round(float(array[:, 1].max()), 3),
    )


@functools.cache
def _drawn(page: int) -> tuple[Drawn, ...]:
    """Return a page's stroked straight-segment paths, in drawing order."""
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    found = []
    for obj in document[page - 1].get_objects(max_depth=15):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        fill, stroked = ctypes.c_int(), ctypes.c_int()
        raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill), ctypes.byref(stroked))
        if not stroked.value:
            continue  # glyphs, and the outlined legend samples of Fig. 5
        if obj.container is not None:
            msg = f"p. {page}: a stroked path inside a form, which is not mapped here"
            raise ValueError(msg)
        matrix = obj.get_matrix()
        subpaths: list[list[tuple[float, float]]] = []
        curved = False
        for index in range(raw.FPDFPath_CountSegments(obj.raw)):
            segment = raw.FPDFPath_GetPathSegment(obj.raw, index)
            kind = raw.FPDFPathSegment_GetType(segment)
            curved = curved or kind == raw.FPDF_SEGMENT_BEZIERTO
            x, y = ctypes.c_float(), ctypes.c_float()
            raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
            point = matrix.on_point(x.value, y.value)
            if kind == raw.FPDF_SEGMENT_MOVETO or not subpaths:
                subpaths.append([])
            subpaths[-1].append((float(point[0]), float(point[1])))
        if curved:
            continue  # no stroked curve is drawn with Bezier segments
        channels = [ctypes.c_uint() for _ in range(4)]
        raw.FPDFPageObj_GetStrokeColor(obj.raw, *(ctypes.byref(c) for c in channels))
        width = ctypes.c_float()
        raw.FPDFPageObj_GetStrokeWidth(obj.raw, ctypes.byref(width))
        found.append(
            Drawn(
                (channels[0].value, channels[1].value, channels[2].value),
                round(float(width.value), 3),
                tuple(np.array(s, dtype=float).reshape(-1, 2) for s in subpaths),
                _clip(obj),
            )
        )
    return tuple(found)


def _near(a: Sequence[float], b: Sequence[float], tolerance: float) -> bool:
    """Whether two boxes agree within a tolerance."""
    return all(abs(x - y) <= tolerance for x, y in zip(a, b, strict=True))


@functools.cache
def frame(key: str) -> tuple[float, float, float, float]:
    """Return a plot's frame (l, b, r, t, pt), checked against its four spines."""
    figure = FIGURES[key]
    left, bottom, right, top = figure.frame
    wanted = [
        ((left, bottom), (left, top)),
        ((right, bottom), (right, top)),
        ((left, bottom), (right, bottom)),
        ((left, top), (right, top)),
    ]
    found = 0
    for (x0, y0), (x1, y1) in wanted:
        for drawn in _drawn(figure.page):
            if drawn.stroke != BLACK or len(drawn.subpaths) != 1:
                continue
            sub = drawn.subpaths[0]
            if len(sub) == 2 and _near(
                (*sorted(sub[:, 0]), *sorted(sub[:, 1])),
                (min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)),
                FRAME_FIND,
            ):
                found += 1
                break
    if found != len(wanted):
        msg = f"{figure.name}: {found} spines at {figure.frame}"
        raise ValueError(msg)
    return figure.frame


@dataclass(frozen=True)
class Ticks:
    """A plot's tick positions along one axis, pt: every one and the major ones."""

    every: tuple[float, ...]
    major: tuple[float, ...]


def _tick(
    sub: np.ndarray, box: tuple[float, float, float, float]
) -> tuple[str, float, bool] | None:
    """Return (axis, position, major) if a subpath is a tick of a frame, else None."""
    if len(sub) != 2:
        return None
    left, bottom, right, top = box
    (x0, y0), (x1, y1) = sub
    length = math.hypot(x1 - x0, y1 - y0)
    if not (MINOR[0] <= length <= MINOR[1] or MAJOR[0] <= length <= MAJOR[1]):
        return None
    major = length >= MAJOR[0]
    if abs(x0 - x1) < EDGE and left - EDGE <= x0 <= right + EDGE:
        for edge in (bottom, top):
            if min(abs(y0 - edge), abs(y1 - edge)) < EDGE:
                return "x", round(float(x0), 3), major
    if abs(y0 - y1) < EDGE and bottom - EDGE <= y0 <= top + EDGE:
        for edge in (left, right):
            if min(abs(x0 - edge), abs(x1 - edge)) < EDGE:
                return "y", round(float(y0), 3), major
    return None


def _merge(found: list[tuple[float, bool]], name: str) -> Ticks:
    """Merge one axis's ticks from its two edges, checking that they agree."""
    found = sorted(found)
    every: list[float] = []
    major: list[float] = []
    for position, is_major in found:
        if every and abs(position - every[-1]) <= MIRROR:
            continue  # the same tick on the far edge, or drawn twice
        every.append(position)
        if is_major:
            major.append(position)
    for position, is_major in found:
        nearest = min(every, key=lambda p: abs(p - position))
        if abs(nearest - position) > MIRROR or (nearest in major) != is_major:
            msg = f"{name}: the ticks of the two edges disagree at {position}"
            raise ValueError(msg)
    return Ticks(tuple(every), tuple(major))


def _ticks(key: str) -> tuple[Ticks, Ticks]:
    """Return a plot's x and y ticks, from both edges of each axis."""
    figure = FIGURES[key]
    box = frame(key)
    found: dict[str, list[tuple[float, bool]]] = {"x": [], "y": []}
    for drawn in _drawn(figure.page):
        if drawn.stroke != BLACK:
            continue
        for sub in drawn.subpaths:
            tick = _tick(sub, box)
            if tick is not None:
                found[tick[0]].append((tick[1], tick[2]))
    return _merge(found["x"], figure.name), _merge(found["y"], figure.name)


def _linear(ticks: Ticks, labels: Sequence[float], step: float | None) -> Axis:
    """Return a linear axis: the major ticks at the labels, the minor every step."""
    if len(ticks.major) != len(labels):
        msg = f"{len(ticks.major)} major x ticks for {len(labels)} labels"
        raise ValueError(msg)
    if step is None:
        if len(ticks.every) != len(ticks.major):
            msg = "minor x ticks where none are expected"
            raise ValueError(msg)
        return Axis(ticks.major, tuple(labels), tolerance=TICK_TOLERANCE)
    rough = Axis(ticks.major, tuple(labels))
    values = []
    for position in ticks.every:
        value = float(rough(position))
        nearest = labels[0] + round((value - labels[0]) / step) * step
        if abs(value - nearest) > ASSIGN * step:
            msg = f"the x tick at {position} is not on a {step} eV step"
            raise ValueError(msg)
        values.append(round(nearest, 9))
    return Axis(ticks.every, tuple(values), tolerance=TICK_TOLERANCE)


def _log(ticks: Ticks, majors: Sequence[float]) -> Axis:
    """Return a log axis: the major ticks at their decades, minor at k x 10^n."""
    if len(ticks.major) != len(majors):
        msg = f"{len(ticks.major)} major y ticks for {len(majors)} values"
        raise ValueError(msg)
    rough = Axis(ticks.major, tuple(majors), log=True)
    values = []
    for position in ticks.every:
        value = math.log10(float(rough(position)))
        decade = math.floor(value + 1e-6)
        k = round(10 ** (value - decade))
        k, decade = (1, decade + 1) if k == 10 else (k, decade)
        if k < 1 or abs(value - decade - math.log10(k)) > ASSIGN * math.log10(
            (k + 1) / k
        ):
            msg = f"the y tick at {position} is not at k x 10^n"
            raise ValueError(msg)
        values.append(float(f"{k}e{decade}"))
    return Axis(ticks.every, tuple(values), log=True, tolerance=TICK_TOLERANCE)


@functools.cache
def axes(key: str) -> tuple[Axis, Axis]:
    """Return a plot's photon-energy (eV) and value (axis unit) axes, checked."""
    figure = FIGURES[key]
    x_ticks, y_ticks = _ticks(key)
    try:
        x_axis = _linear(x_ticks, figure.x_major, figure.x_minor)
        y_axis = _log(y_ticks, figure.y_major)
        x_axis.check()
        y_axis.check()
    except ValueError as error:
        msg = f"{figure.name}: {error}"
        raise ValueError(msg) from error
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
    """Return a curve's vertices inside its frame (pt), in drawing order."""
    curve = CURVES[name]
    figure = FIGURES[curve.figure]
    box = frame(curve.figure)
    found = [
        d
        for d in _drawn(figure.page)
        if d.stroke == curve.colour
        and abs(d.width - curve.width) < 1e-3
        and sum(len(s) for s in d.subpaths) > MIN_VERTICES
        and d.clip is not None
        and _near(d.clip, box, FRAME_FIND)
    ]
    if not found:
        msg = f"{figure.name}: no curve for {name}"
        raise ValueError(msg)
    line = np.vstack([s for d in found for s in d.subpaths])
    left, bottom, right, top = box
    inside = (
        (line[:, 0] >= left - EDGE)
        & (line[:, 0] <= right + EDGE)
        & (line[:, 1] >= bottom - EDGE)
        & (line[:, 1] <= top + EDGE)
    )
    line = line[inside]
    keep = np.ones(len(line), dtype=bool)
    keep[1:] = np.any(np.diff(line, axis=0) != 0, axis=1)
    line = line[keep]
    if np.any(np.diff(line[:, 0]) < -BACK):
        msg = f"{figure.name}: {name} steps back in photon energy"
        raise ValueError(msg)
    return line


@functools.cache
def curve(name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a curve's page positions (pt), photon energies (eV) and values.

    The values are in cm^2N s^(N-1): the axis's printed unit times
    :attr:`Figure.unit`.
    """
    key = CURVES[name].figure
    line = polyline(name)
    x_axis, y_axis = axes(key)
    return (
        line,
        np.asarray(x_axis(line[:, 0]), dtype=float),
        np.asarray(y_axis(line[:, 1]), dtype=float) * FIGURES[key].unit,
    )


def _points(name: str) -> list[Point]:
    """Return a curve's vertices as data points."""
    line, energies, values = curve(name)
    return [
        Point(
            f"v{k:04d}",
            f"{e:.5f}",
            f"{v:.4e}",
            "",
            f"vertex {k} at ({x:.3f}, {y:.3f}) pt",
        )
        for k, ((x, y), e, v) in enumerate(zip(line, energies, values, strict=True), 1)
    ]


def at(name: str, energy: float) -> float:
    """Return a curve's value at a photon energy (eV), log-interpolated."""
    _, energies, values = curve(name)
    above = int(np.argmax(energies >= energy))
    if energies[above] < energy or above == 0:
        msg = f"{name}: {energy} eV is outside the curve"
        raise ValueError(msg)
    e0, e1 = energies[above - 1], energies[above]
    v0, v1 = values[above - 1], values[above]
    t = 0.0 if e1 == e0 else (energy - e0) / (e1 - e0)
    return float(10 ** ((1 - t) * math.log10(v0) + t * math.log10(v1)))


def _peak(name: str, low: float, high: float) -> float:
    """Return the photon energy (eV) of a curve's highest vertex in a window."""
    _, energies, values = curve(name)
    window = (energies >= low) & (energies <= high)
    return float(energies[window][int(np.argmax(values[window]))])


def _jump(name: str, low: float, high: float) -> tuple[float, float]:
    """Return where (eV) and by how much (decades) a curve changes most in a window.

    Between consecutive vertices, at the middle of the pair.
    """
    _, energies, values = curve(name)
    window = np.flatnonzero((energies >= low) & (energies <= high))
    steps = np.abs(np.diff(np.log10(values[window])))
    k = int(np.argmax(steps))
    middle = (energies[window[k]] + energies[window[k + 1]]) / 2
    return float(middle), float(steps[k])


PEAK = 0.1
"""The core-excited peak is printed as 'approximately 40.8 eV': within this, eV."""
THRESHOLD = 0.1
"""A threshold printed to 0.1 eV or 0.05 eV is matched within this, eV."""
TWICE = 0.03
"""One calculation drawn in two figures agrees within this: two readings between
vertices of simplified paths (0.1 pt at 8 to 33 pt per decade, up to 3 %)."""
AU_EV = 27.211386
"""1 a.u. of energy in eV (CODATA 2018), where Fig. 4's dotted line is drawn."""


def _check_peaks(out: dict[str, object], failures: list[str]) -> None:
    """Check the core-excited peak, 'At hbar omega ~ 40.8 eV' (p. 043108-7)."""
    for name in ("3b_L0", "3b_L2", "3d_L1", "3d_L3", "3f_L0", "3f_L2", "3f_L4"):
        read = _peak(name, 39.5, 42.0)
        out[f"{name} core-excited peak (eV)"] = (read, 40.8)
        if abs(read - 40.8) > PEAK:
            failures.append(f"{name}: the peak at {read:.3f} eV, printed 40.8")


def _check_jumps(out: dict[str, object], failures: list[str]) -> None:
    """Check the thresholds where the cross sections jump (pp. 043108-6, -7).

    'sharp jumps in sigma_L^(K) mark higher ionization thresholds ... the 2 l_a
    epsilon_2 l channels at 32.7, 21.8, and 16.35 eV in Figs. 3(b), 3(d), and 3(f)': the
    largest step between vertices within 0.3 eV of each, which must be there in (b) and
    in (d)'s L = 1. In (d)'s L = 3 and in (f) the step is smaller than the resonances
    beside it (reported).
    """
    for name, printed, required in (
        ("3b_L0", 32.7, True),
        ("3b_L2", 32.7, True),
        ("3d_L1", 21.8, True),
        ("3d_L3", 21.8, False),
        ("3f_L0", 16.35, False),
        ("3f_L2", 16.35, False),
        ("3f_L4", 16.35, False),
    ):
        where, size = _jump(name, printed - 0.3, printed + 0.3)
        out[f"{name} largest step near {printed} eV (eV, decades)"] = (where, size)
        if required and abs(where - printed) > THRESHOLD:
            failures.append(f"{name}: the jump at {where:.3f} eV, printed {printed}")


def _check_fig4(out: dict[str, object], failures: list[str]) -> None:
    """Check Fig. 4's dotted line at 1 a.u. and its n_a <= 5 curve against Fig. 3(d).

    'the peak positions converge to the limit of 1 a.u. (27.2 eV), which is marked with
    a dotted vertical line in Fig. 4'; Fig. 4 is 'an enlarged view' of Fig. 3(d) from
    25.2 to 27.2 eV, whose n_a <= 5 calculation is Fig. 3(d)'s L = 1 and 3 summed (p.
    043108-7).
    """
    x_axis, _ = axes("4")
    lines = [
        d
        for d in _drawn(FIGURES["4"].page)
        if d.stroke == GREY and len(d.subpaths) == 1 and len(d.subpaths[0]) == 2
    ]
    if len(lines) != 1:
        failures.append(f"Fig. 4: {len(lines)} grey dotted lines")
    else:
        read = float(x_axis(float(lines[0].subpaths[0][0, 0])))
        out["Fig. 4 dotted line (eV)"] = (read, AU_EV)
        if abs(read - AU_EV) > THRESHOLD / 10:
            failures.append(f"Fig. 4's dotted line at {read:.4f} eV, not 1 a.u.")
    for energy in (25.3, 25.7, 25.85, 26.3, 26.6, 27.0, 27.3):
        total = at("3d_L1", energy) + at("3d_L3", energy)
        read = at("4_na5", energy)
        out[f"Fig. 4 n_a <= 5 over Fig. 3(d) L = 1 + 3 at {energy} eV"] = read / total
        if abs(read / total - 1) > TWICE:
            failures.append(f"Fig. 4 and Fig. 3(d) differ at {energy} eV")


def _report_fig5(out: dict[str, object]) -> None:
    """Report the largest excess of a channel of Fig. 5 over its total."""
    for total, parts in (
        ("5a_total", ("5a_1ses", "5a_2ses", "5a_2pep")),
        ("5b_total", ("5b_1sed", "5b_2sed", "5b_2pep", "5b_2pef")),
    ):
        _, energies, _ = curve(total)
        worst = 0.0
        for part in parts:
            _, part_energies, part_values = curve(part)
            for energy, value in zip(part_energies, part_values, strict=True):
                if energies[0] < energy < energies[-1]:
                    worst = max(worst, value / at(total, float(energy)) - 1)
        out[f"{total}: the largest excess of a channel over it"] = worst


def checks() -> dict[str, object]:
    """Return the checks against the numbers the text prints. Raises on a miss."""
    out: dict[str, object] = {}
    failures: list[str] = []
    _check_peaks(out, failures)
    _check_jumps(out, failures)
    _check_fig4(out, failures)
    _report_fig5(out)
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return every curve, after the checks."""
    checks()
    return {curve.file: _points(name) for name, curve in CURVES.items()}


def main() -> None:
    """Write the data files, and print the readings and the checks."""
    write(extract())
    for name, spec in CURVES.items():
        line, energies, values = curve(name)
        x_res, y_res = residuals(spec.figure)
        print(
            f"{name}: {len(line)} vertices, {energies[0]:.4f} to {energies[-1]:.4f} "
            f"eV, {values.min():.3e} to {values.max():.3e}; ticks within "
            f"{x_res:.3f} (x) and {y_res:.3f} (y) pt"
        )
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
