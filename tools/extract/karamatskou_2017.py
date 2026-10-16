# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Karamatskou and Santra (2017): TDCIS two-photon ionisation of Xe 5p.

The paper prints no table and no supplement exists: its generalised two-photon
ionisation cross sections are given only in figures, vector drawings on the page
(gnuplot; the PDF's images are the arXiv stamp and the page furniture), read from the
drawing (data rule 8, TC-699). The library copy is the arXiv preprint
(arXiv:1701.04789v1); its page numbers are the PDF's.

* **Fig. 1** (p. 5, library figure ``fig_f3d1dc935443``): the interchannel and
  intrachannel cross sections, linear axis, 9.0 to 12.55 eV.
* **Fig. 2** (p. 6, ``fig_551a2ef821b9``): the interchannel cross section again (left
  axis, linear), with the one-photon excitation cross section (right axis, 10^3 Mb),
  which is not read.
* **Fig. 3** (p. 8, ``fig_92f6093163b2``): the interchannel cross section on a log axis,
  shifted in energy by the authors ('When shifting the 5p binding energy about 0.3 eV
  down to the energy found within RPAE', p. 8), with Wendin, Jonsson and L'Huillier's
  RPAE curve (J. Opt. Soc. Am. B 4, 833 (1987)), which is not read.

How the drawing is read:

* the frames: a stroked closed rectangle (five vertices) per plot, found near the
  position :data:`PANELS` gives; the ticks are stroked two-vertex segments of the same
  stroke, perpendicular to the bottom (x) or left (y) edge with one end on it, long
  (4.08 pt) at the labelled values and short (2.01 pt) between.
* the axes: the labels are text objects; their values are written in :data:`PANELS`
  as printed and checked against the text objects under the bottom and left of each
  frame (:func:`label_check`). The long ticks are fitted to the labels; each short tick
  takes the nearest value of the axis's grid (the label interval divided by the number
  of short ticks between two labels, or 1 to 9 times a decade on a log axis), and every
  tick must lie on the scale fitted to all of them within :data:`TICK_TOLERANCE`.
* the curves: each is one stroked path of its colour (:data:`CURVES`), a polyline
  through the computed points (gnuplot joins the data points with straight lines),
  plus the legend's two-vertex sample, checked to be horizontal and outside the plotted
  data's reach (above :data:`LEGEND` of the frame's height) and dropped. A vertex on an
  edge of the frame is where gnuplot clipped the line, not a computed point, and is
  dropped (one at each end of Fig. 1's interchannel curve and of Fig. 2's, one at the
  start of Fig. 3's).
* the values: every remaining vertex, in drawing order (increasing photon energy,
  checked). The drawing's coordinates are quantised in steps of 0.065 pt. Each
  ``printed`` cell gives the vertex's position in page points.

**Checks** (:func:`checks`). The paper prints no value of the cross section. The same
calculation is drawn in Figs. 1, 2 and 3: each point of Figs. 2 and 3 must pair with a
point of Fig. 1 at one energy offset (within :data:`SPREAD`) and with its value
(within :data:`VALUE` or :data:`RELATIVE`); the offset must be small for Fig. 2
(:data:`OFFSET`; it is 0.0137 eV, a finding) and 'about 0.3 eV' down for Fig. 3.
"""

import ctypes
import functools
import math
import re
from dataclasses import dataclass

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, pdf, texts
from _library import Point, write
from pypdfium2 import raw

DIRECTORY = "Karamatskou_2017_TimeDependentConfiguration"
SOURCE = "karamatskou_2017"
EDGE = 0.02
"""A tick's end on its frame edge, or a clip vertex on it, within this, pt."""
FRAME_FIND = 1.0
"""The frame's corners lie within this of :data:`PANELS`' position, pt."""
TICK_MAX = 5.0
"""Ticks are shorter than this, pt (the long ones 4.08)."""
LONG = 3.0
"""Ticks longer than this are the long, labelled ones, pt (the short ones 2.01)."""
TICK_TOLERANCE = 0.065
"""Every tick lies on its fitted scale within this, pt: the drawing's quantum."""
MATCH = 0.2
"""A y label's text object lies within this of its long tick, after the offset common
to the axis, pt."""
LEGEND = 0.85
"""The legend's sample lies higher than this fraction of the frame's height."""
UNIT = 1e-48
"""Every y axis prints its values in units of 10^-48 cm^4 s."""


@dataclass(frozen=True)
class Panel:
    """One plot: where it is and what its axes print."""

    name: str
    page: int
    frame: tuple[float, float, float, float]
    """The frame's left, bottom, right, top, pt, to within :data:`FRAME_FIND`."""
    x_labels: tuple[str, ...]
    y_labels: tuple[str, ...]
    log: bool = False


PANELS = {
    "fig1": Panel(
        "Fig. 1",
        5,
        (127.56, 447.34, 529.65, 711.0),
        ("9", "9.5", "10", "10.5", "11", "11.5", "12", "12.5"),
        ("0", "50", "100", "150", "200"),
    ),
    "fig2": Panel(
        "Fig. 2",
        6,
        (127.56, 444.1, 494.3, 707.8),
        ("9.5", "10", "10.5", "11", "11.5", "12"),
        ("0", "50", "100", "150"),
    ),
    "fig3": Panel(
        "Fig. 3",
        8,
        (141.1, 481.0, 520.0, 705.4),
        ("9", "9.5", "10", "10.5", "11", "11.5", "12"),
        ("0.1", "1", "10", "100", "1000"),
        log=True,
    ),
}


@dataclass(frozen=True)
class Curve:
    """One drawn curve: its plot, its stroke colour and its data file."""

    panel: str
    colour: tuple[int, int, int]
    width: float
    """The stroke width (in the path's own units; Fig. 2's left axis line has the
    two-photon curve's colour but not its width)."""
    file: str


CURVES = {
    "fig1_interchannel": Curve(
        "fig1", (0, 0, 140), 12.95, "Xe/karamatskou_2017__fig1_interchannel.csv"
    ),
    "fig1_intrachannel": Curve(
        "fig1", (0, 191, 255), 12.95, "Xe/karamatskou_2017__fig1_intrachannel.csv"
    ),
    "fig2_interchannel": Curve(
        "fig2", (26, 26, 143), 12.95, "Xe/karamatskou_2017__fig2_interchannel.csv"
    ),
    "fig3_interchannel_shifted": Curve(
        "fig3",
        (0, 0, 128),
        11.021,
        "Xe/karamatskou_2017__fig3_interchannel_shifted.csv",
    ),
}


@dataclass(frozen=True)
class Stroke:
    """A stroked subpath (page points), its stroke width and colour."""

    width: float
    colour: tuple[int, int, int]
    vertices: np.ndarray

    @property
    def length(self) -> float:
        """The length of the subpath, pt."""
        return float(np.sum(np.hypot(*np.diff(self.vertices, axis=0).T)))


@functools.cache
def _strokes(page: int) -> tuple[Stroke, ...]:
    """Return a page's stroked paths, split into subpaths at each move.

    Raises if a stroked path lies in a form (not mapped here) or has a Bezier segment
    (the figures are polylines; checked).
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
        channels = [ctypes.c_uint() for _ in range(4)]
        raw.FPDFPageObj_GetStrokeColor(obj.raw, *(ctypes.byref(c) for c in channels))
        colour = (channels[0].value, channels[1].value, channels[2].value)
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
        found.extend(
            Stroke(
                round(float(width.value), 3),
                colour,
                np.array(s, dtype=float).reshape(-1, 2),
            )
            for s in subpaths
        )
    return tuple(found)


@functools.cache
def _frame_stroke(key: str) -> tuple[tuple[float, float, float, float], float]:
    """Return a plot's frame (l, b, r, t, pt) and its stroke width.

    The frame is the stroked closed five-vertex rectangle near :data:`PANELS`'
    position; it is drawn twice (once before the curves, once after), and both must
    agree.
    """
    panel = PANELS[key]
    found = set()
    for stroke in _strokes(panel.page):
        sub = stroke.vertices
        if len(sub) != 5 or not np.allclose(sub[0], sub[-1]):
            continue
        box = (
            float(sub[:, 0].min()),
            float(sub[:, 1].min()),
            float(sub[:, 0].max()),
            float(sub[:, 1].max()),
        )
        if all(abs(a - b) < FRAME_FIND for a, b in zip(box, panel.frame, strict=True)):
            found.add((tuple(round(v, 3) for v in box), stroke.width))
    if len(found) != 1:
        msg = f"{panel.name}: {len(found)} distinct frames near {panel.frame}"
        raise ValueError(msg)
    ((box, width),) = found
    return (box[0], box[1], box[2], box[3]), width


def frame(key: str) -> tuple[float, float, float, float]:
    """Return a plot's frame: left, bottom, right, top, pt."""
    return _frame_stroke(key)[0]


@dataclass(frozen=True)
class Ticks:
    """A plot's tick positions along one edge, pt: all, and the long ones."""

    every: tuple[float, ...]
    long: tuple[float, ...]


def _ticks(key: str) -> tuple[Ticks, Ticks]:
    """Return the x ticks (bottom edge) and y ticks (left edge) of a plot."""
    panel = PANELS[key]
    (left, bottom, right, top), width = _frame_stroke(key)
    x_ticks: list[tuple[float, float]] = []
    y_ticks: list[tuple[float, float]] = []
    for stroke in _strokes(panel.page):
        if len(stroke.vertices) != 2 or abs(stroke.width - width) > 1e-3:
            continue
        (x0, y0), (x1, y1) = stroke.vertices
        length = stroke.length
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

    The cross section is in the axes' printed unit, 10^-48 cm^4 s.
    """
    panel = PANELS[key]
    x_ticks, y_ticks = _ticks(key)
    x_axis = _grid(x_ticks, panel.x_labels, log=False, name=f"{panel.name} x")
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

    The x labels are the numbers of the text objects just below the frame, left to
    right; the y labels those just left of it, bottom to top, each within
    :data:`MATCH` of its long tick after the offset common to the axis. Returns that
    offset, pt.
    """
    panel = PANELS[key]
    left, bottom, right, top = frame(key)
    x_ticks, y_ticks = _ticks(key)
    found = texts(DIRECTORY, panel.page, (left - 30.0, bottom - 20.0, right, top + 5))
    below = sorted(
        (x, t)
        for t, x, y in found
        if bottom - 20.0 < y < bottom - 5.0 and _NUMBER.match(t)
    )
    if tuple(t for _, t in below) != panel.x_labels:
        msg = f"{panel.name}: the x labels print {below}, not {panel.x_labels}"
        raise ValueError(msg)
    offsets = np.array([x for x, _ in below]) - np.array(x_ticks.long)
    if np.max(np.abs(offsets - np.median(offsets))) > 2.0:
        msg = f"{panel.name}: the x labels are not under their long ticks: {offsets}"
        raise ValueError(msg)
    beside = sorted(
        (y, t)
        for t, x, y in found
        if left - 30.0 < x < left - 2.0 and y > bottom - 5.0 and _NUMBER.match(t)
    )
    if tuple(t for _, t in beside) != panel.y_labels:
        msg = f"{panel.name}: the y labels print {beside}, not {panel.y_labels}"
        raise ValueError(msg)
    offsets = np.array([y for y, _ in beside]) - np.array(y_ticks.long)
    offset = float(np.median(offsets))
    if np.max(np.abs(offsets - offset)) > MATCH:
        msg = f"{panel.name}: the y labels are not at their long ticks: {offsets}"
        raise ValueError(msg)
    return offset


@functools.cache
def curve(key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a curve: its page positions (pt), photon energies (eV), values (cm^4 s).

    The longest subpath of the curve's stroke inside the frame; the stroke's other
    subpaths must be the legend's horizontal sample, high in the frame. Vertices on
    the frame's edges (gnuplot's clip points) are dropped; the rest must increase in
    photon energy.
    """
    spec = CURVES[key]
    panel = PANELS[spec.panel]
    left, bottom, right, top = frame(spec.panel)
    mine = [
        s
        for s in _strokes(panel.page)
        if s.colour == spec.colour
        and abs(s.width - spec.width) < 1e-3
        and np.all(s.vertices[:, 0] >= left - EDGE)
        and np.all(s.vertices[:, 0] <= right + EDGE)
        and np.all(s.vertices[:, 1] >= bottom - EDGE)
        and np.all(s.vertices[:, 1] <= top + EDGE)
    ]
    mine.sort(key=lambda s: -len(s.vertices))
    if not mine or len(mine[0].vertices) < 10:
        msg = f"{key}: no curve of colour {spec.colour}"
        raise ValueError(msg)
    for other in mine[1:]:
        flat = len(other.vertices) == 2 and np.ptp(other.vertices[:, 1]) < EDGE
        if not flat or other.vertices[0, 1] < bottom + LEGEND * (top - bottom):
            msg = f"{key}: a second subpath of its colour, not the legend's sample"
            raise ValueError(msg)
    line = mine[0].vertices
    if line[0, 0] > line[-1, 0]:
        line = line[::-1]
    on_edge = (
        (np.abs(line[:, 0] - left) < EDGE)
        | (np.abs(line[:, 0] - right) < EDGE)
        | (np.abs(line[:, 1] - bottom) < EDGE)
        | (np.abs(line[:, 1] - top) < EDGE)
    )
    line = line[~on_edge]
    if np.any(np.diff(line[:, 0]) <= 0):
        msg = f"{key}: the curve does not increase in photon energy"
        raise ValueError(msg)
    x_axis, y_axis = axes(spec.panel)
    energy = np.asarray(x_axis(line[:, 0]), dtype=float)
    value = np.asarray(y_axis(line[:, 1]), dtype=float) * UNIT
    return line, energy, value


def clipped(key: str) -> int:
    """Return how many of a curve's vertices lie on the frame's edges (dropped)."""
    spec = CURVES[key]
    panel = PANELS[spec.panel]
    longest = max(
        (
            s
            for s in _strokes(panel.page)
            if s.colour == spec.colour and abs(s.width - spec.width) < 1e-3
        ),
        key=lambda s: len(s.vertices),
    )
    return len(longest.vertices) - len(curve(key)[0])


def _points(key: str) -> list[Point]:
    """Return a curve as data points, one per computed point."""
    line, energy, value = curve(key)
    return [
        Point(
            f"p{k:03d}",
            f"{e:.4f}",
            f"{v:.4e}",
            "",
            f"vertex at ({x:.3f}, {y:.3f}) pt",
        )
        for k, ((x, y), e, v) in enumerate(zip(line, energy, value, strict=True), 1)
    ]


def peaks(key: str, floor: float) -> list[tuple[float, float]]:
    """Return a curve's local maxima above a value (cm^4 s): (eV, cm^4 s)."""
    _, energies, values = curve(key)
    return [
        (float(energies[k]), float(values[k]))
        for k in range(1, len(values) - 1)
        if values[k] >= values[k - 1]
        and values[k] > values[k + 1]
        and values[k] > floor
    ]


@dataclass(frozen=True)
class Pairing:
    """How one drawn curve's points map onto another's."""

    shift: float
    """The median energy difference of the paired points (second minus first), eV."""
    spread: float
    """The largest departure of a pair's energy difference from the median, eV."""
    offset: int
    """The second curve's point i is paired with the first's point i + offset."""
    pairs: int
    worst: float
    """The largest value difference of a pair, cm^4 s, beyond the margin allowed it."""


def pairing(first: str, second: str, absolute: float, relative: float) -> Pairing:
    """Return how a curve's points map onto another's, drawn from one calculation.

    Each point of the second is paired with a point of the first at a fixed offset in
    the list, the offset (within 10) that makes the energy differences most nearly
    constant. ``worst`` is the largest amount by which a pair's values differ beyond
    ``absolute`` (cm^4 s) or ``relative`` of the larger, whichever is greater (0 if
    every pair is within its margin).
    """
    _, e1, v1 = curve(first)
    _, e2, v2 = curve(second)
    best: Pairing | None = None
    for offset in range(-10, 11):
        pairs = [(i + offset, i) for i in range(len(e2)) if 0 <= i + offset < len(e1)]
        de = np.array([e2[i] - e1[j] for j, i in pairs])
        spread = float(np.max(np.abs(de - np.median(de))))
        if best is None or spread < best.spread:
            excess = [
                abs(v2[i] - v1[j]) - max(absolute, relative * max(v1[j], v2[i]))
                for j, i in pairs
            ]
            best = Pairing(
                float(np.median(de)), spread, offset, len(pairs), max(0.0, *excess)
            )
    if best is None:  # pragma: no cover (the loop always runs)
        msg = f"{first}, {second}: nothing to pair"
        raise ValueError(msg)
    return best


# Margins. SPREAD: two drawings of one calculation keep the spacing of its points
# within their quanta (0.065 pt: 0.0006 eV). VALUE: their values agree within the
# linear axes' quantum (0.065 pt is 0.06 x 10^-48 cm^4 s in Fig. 1, 0.04 in Fig. 2),
# here 0.1 x 10^-48 cm^4 s, or within 3 % (Fig. 3's log axis: 0.065 pt is 0.3 %, with
# the tick residuals 0.6 %; a computed point near a peak is steep in energy).
# OFFSET: Figs. 1 and 2 may differ by a small offset in energy (a finding; 0.0137 eV
# read), not by more than 0.02 eV. ABOUT: the text's 'about 0.3 eV' (p. 8).
SPREAD = 0.002
VALUE = 0.1 * UNIT
RELATIVE = 0.03
OFFSET = 0.02
ABOUT = 0.05


def checks() -> dict[str, object]:
    """Return the checks of the reading. Raises if one fails its margin."""
    out: dict[str, object] = {}
    failures = []
    for key, panel in PANELS.items():
        out[f"{panel.name} y label offset (pt)"] = label_check(key)
        out[f"{panel.name} worst tick residuals (pt)"] = residuals(key)
    for key in CURVES:
        out[f"{key} clip vertices dropped"] = clipped(key)
    # Figs. 1 and 2 draw one interchannel calculation (p. 7: 'the interchannel
    # two-photon ionization cross section shown in Fig. 1 is compared to the
    # one-photon excitation cross section'); Fig. 3 draws it moved in energy ('When
    # shifting the 5p binding energy about 0.3 eV down', p. 8).
    for second, low, high in (
        ("fig2_interchannel", -OFFSET, OFFSET),
        ("fig3_interchannel_shifted", -0.3 - ABOUT, -0.3 + ABOUT),
    ):
        found = pairing("fig1_interchannel", second, VALUE, RELATIVE)
        out[f"{second} against fig1_interchannel"] = found
        if found.spread > SPREAD or found.worst > 0 or not low <= found.shift <= high:
            failures.append(f"{second} does not map onto Fig. 1: {found}")
    for key in CURVES:
        out[f"{key} peaks above 2e-48 (eV, cm4 s)"] = peaks(key, 2 * UNIT)
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
            f"{key}: {len(line)} points, {energies[0]:.4f} to {energies[-1]:.4f} eV, "
            f"{values.min():.3e} to {values.max():.3e} cm4 s"
        )
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
