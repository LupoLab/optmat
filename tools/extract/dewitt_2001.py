# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract DeWitt, Wells & Jones (2001): O2:Xe and N2:Ar ion-yield ratios at 790 nm.

The paper prints the ratios against intensity only in two figures; the text gives
round numbers (the record's inline datasets). Both figures are vector drawings on
p. 3 (PDF page 3), read from the drawing (data rule 8, TC-699):

**Fig. 1** (library figure ``fig_4f3a1feb7459``): the ratio of the molecular to the
Xe+ yield against the intensity (W/cm^2, log axis). Only the O2 points of this work
('O2: 790 nm, 100 fs', open squares) are taken: not the S2 points (filled squares,
filled circles) nor Guo et al.'s O2 points (open circles; shown 'without their
detection efficiency correction', the caption).

**Fig. 2** (library figure ``fig_9c453cfd77c2``): the ratio to the Ar+ yield; only the
N2 points of this work ('N2: 790 nm, 100 fs', open squares) are taken: not F2
(filled circles) nor Guo et al.'s N2 (open triangles).

How the drawing is read:

* the paths: every path object of the page is split into its subpaths (at each
  move-to), in page points; one object can hold many markers or many bars.
* the axes: the frame is a stroked rectangle; the tick marks are its stroked
  two-point subpaths perpendicular to the bottom (x) or left (y) edge, starting on
  it, under 4 pt long; long ticks (over 2.5 pt) are labelled. The labels are text
  objects: the x labels '10' with an exponent (the text to its upper right), matched
  to the long tick within 2.5 pt of the label's centre; the y labels '0.5' to '3.5',
  matched to the long tick within 2 pt of the label's height. On the log x axis a
  short tick between the decades 10^n and 10^(n+1) is k x 10^n, k counted from 2
  (below the first decade counted down from 9, above the last up from 2); on the
  linear y axis the short ticks divide the step between labels evenly (and continue
  it beyond them). Every tick must lie on the fitted scale within ``TICK_TOLERANCE``
  (the drawing's coordinates are quantised: the y ticks' spacing alternates 5.30 and
  5.41 pt).
* the markers: closed four-segment straight subpaths, axis-aligned squares 3.5-3.9 pt
  a side, stroked and not filled (the S2 squares of Fig. 1 are filled), inside the
  frame and outside the legend box; their centres.
* the error bars: stroked vertical two-point subpaths at a marker's x (within
  ``ON_MARKER`` pt) with one end inside the marker's square; the other end is the
  bar's end, and a stroked horizontal cap must be drawn there across the bar. In
  Fig. 2 every N2 marker has a bar up and a bar down; in Fig. 1 only one O2 marker
  has bars ('Errors smaller than symbol size not shown', the legend), and the others'
  ``sigma`` is left empty.

``sigma`` is half the bar's full length in the ratio (the axis is linear); the bars'
definition is not stated. Every ``printed`` cell gives the marker's centre and the
bars' ends in page points. The rows are numbered in order of intensity (then ratio).

:func:`checks` compares the points with the text: the O2:Xe ratio 'between 0.15 and
0.2' at low intensity (p. 2) and the N2:Ar ratio 'approximately 1.7' (p. 3).
"""

import ctypes
import math
from dataclasses import dataclass

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, pdf, texts
from _library import Point, write
from pypdfium2 import raw

DIRECTORY = "DeWitt_2001_RatiometricComparisonIntense"
SOURCE = "dewitt_2001"
PAGE = 3
# The largest residual of a tick from the fitted scale, pt: the drawing's coordinates
# are quantised to about 0.1 pt (the y ticks' spacing alternates 5.30 and 5.41 pt).
TICK_TOLERANCE = 0.1
# A tick is under 4 pt long; a long (labelled) tick is over 2.5 pt (they are 3.1 pt,
# the short ones 1.6 pt).
TICK_LONGEST = 4.0
TICK_LONG = 2.5
# A tick starts on the frame's edge, a bar at the marker's x: within 0.05 pt.
ON_EDGE = 0.05
ON_MARKER = 0.05
# The data markers are squares 3.65-3.68 pt a side (the legend's 5.1 pt).
SQUARE = (3.5, 3.9)


@dataclass(frozen=True)
class Subpath:
    """A subpath of a path object: its vertices (pt) and how the object is painted."""

    vertices: np.ndarray
    closed: bool
    curved: bool
    stroked: bool
    filled: bool


def _matrices(obj: pdfium.PdfObject, points: np.ndarray) -> np.ndarray:
    """Map an object's points to the page through it and its enclosing forms."""
    out = np.array([obj.get_matrix().on_point(x, y) for x, y in points], dtype=float)
    container = obj.container
    while container is not None:
        matrix = container.get_matrix()
        out = np.array([matrix.on_point(x, y) for x, y in out], dtype=float)
        container = container.container
    return out.reshape(-1, 2)


def subpaths(
    directory: str, page: int, region: tuple[float, float, float, float]
) -> list[Subpath]:
    """Return the subpaths of the path objects touching a region (l, b, r, t, pt)."""
    left, bottom, right, top = region
    document = pdfium.PdfDocument(str(pdf(directory)))
    found: list[Subpath] = []
    for obj in document[page - 1].get_objects(max_depth=15):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        west, south, east, north = obj.get_bounds()
        if not (east >= left and west <= right and north >= bottom and south <= top):
            continue
        fill, stroke = ctypes.c_int(), ctypes.c_int()
        raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill), ctypes.byref(stroke))
        groups: list[list[tuple[float, float, int, bool]]] = []
        for index in range(raw.FPDFPath_CountSegments(obj.raw)):
            segment = raw.FPDFPath_GetPathSegment(obj.raw, index)
            x, y = ctypes.c_float(), ctypes.c_float()
            raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
            kind = raw.FPDFPathSegment_GetType(segment)
            closes = bool(raw.FPDFPathSegment_GetClose(segment))
            if kind == raw.FPDF_SEGMENT_MOVETO or not groups:
                groups.append([])
            groups[-1].append((x.value, y.value, kind, closes))
        for group in groups:
            points = _matrices(obj, np.array([(g[0], g[1]) for g in group]))
            found.append(
                Subpath(
                    points,
                    closed=group[-1][3],
                    curved=any(g[2] == raw.FPDF_SEGMENT_BEZIERTO for g in group),
                    stroked=bool(stroke.value),
                    filled=bool(fill.value),
                )
            )
    return found


def _segments(drawn: list[Subpath]) -> list[np.ndarray]:
    """Return the stroked straight two-point subpaths, as 2x2 arrays."""
    return [
        s.vertices
        for s in drawn
        if s.stroked and not s.curved and not s.closed and len(s.vertices) == 2
    ]


@dataclass(frozen=True)
class Figure:
    """A figure: its frame and legend box (l, b, r, t, pt) and its axes' labels.

    Attributes
    ----------
    name : str
        The figure, for messages.
    frame : tuple of float
        The plot frame, from the drawing.
    legend : tuple of float
        The legend's box, from the drawing (its markers are not points).
    y_labels : tuple of str
        The y axis's labels, as printed.
    """

    name: str
    frame: tuple[float, float, float, float]
    legend: tuple[float, float, float, float]
    y_labels: tuple[str, ...]


FIG1 = Figure(
    "Fig. 1",
    (88.39, 553.29, 295.75, 739.91),
    (120.41, 663.42, 287.34, 734.61),
    ("0.5", "1.0", "1.5", "2.0", "2.5", "3.0"),
)
FIG2 = Figure(
    "Fig. 2",
    (347.97, 576.57, 553.17, 735.26),
    (364.27, 687.95, 523.18, 729.79),
    ("0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5"),
)


def _inside(x: float, y: float, box: tuple[float, float, float, float]) -> bool:
    return box[0] <= x <= box[2] and box[1] <= y <= box[3]


def _ticks(
    segments: list[np.ndarray], figure: Figure, axis: int
) -> list[tuple[float, bool]]:
    """Return a figure's ticks on the bottom (0) or left (1) edge: (position, long)."""
    left, bottom, right, top = figure.frame
    edge = bottom if axis == 0 else left
    along, across = (0, 1) if axis == 0 else (1, 0)
    low, high = (left, right) if axis == 0 else (bottom, top)
    ticks: dict[float, bool] = {}
    for s in segments:
        if abs(s[0, along] - s[1, along]) > ON_EDGE:
            continue
        ends = sorted(s[:, across])
        length = ends[1] - ends[0]
        if abs(ends[0] - edge) > ON_EDGE or not 0 < length < TICK_LONGEST:
            continue
        position = round(float(s[:, along].mean()), 3)
        if low - ON_EDGE <= position <= high + ON_EDGE:
            ticks[position] = ticks.get(position, False) or length > TICK_LONG
    return sorted(ticks.items())


def _x_axis(segments: list[np.ndarray], figure: Figure) -> Axis:
    """Return the log intensity axis, W/cm^2: decades from the labels, k x 10^n."""
    ticks = _ticks(segments, figure, 0)
    words = texts(
        DIRECTORY,
        PAGE,
        (
            figure.frame[0] - 20,
            figure.frame[1] - 20,
            figure.frame[2] + 20,
            figure.frame[1],
        ),
    )
    labels = []
    for word, x, y in words:
        if word != "10":
            continue
        exponents = [
            (w, ex, ey)
            for w, ex, ey in words
            if w.isdigit() and 0 < ex - x < 10 and 0 < ey - y < 6
        ]
        if len(exponents) != 1:
            msg = f"{figure.name}: the label '10' at {x:.1f} has no exponent"
            raise ValueError(msg)
        labels.append(((x + exponents[0][1]) / 2, int(exponents[0][0])))
    decades: dict[float, int] = {}
    for position, long in ticks:
        if not long:
            continue
        near = [n for c, n in labels if abs(c - position) < 2.5]
        if len(near) != 1:
            msg = f"{figure.name}: {len(near)} labels at the long x tick {position}"
            raise ValueError(msg)
        decades[position] = near[0]
    if len(decades) != len(labels):
        msg = f"{figure.name}: {len(labels)} x labels, {len(decades)} long ticks"
        raise ValueError(msg)
    first = min(decades)
    values = []
    for position, _ in ticks:
        below = [p for p in decades if p <= position]
        if below:
            base = max(below)
            rank = sum(1 for p, _ in ticks if base < p <= position)
            values.append((rank + 1) * 10.0 ** decades[base])
        else:
            rank = sum(1 for p, _ in ticks if position <= p < first)
            values.append((10 - rank) * 10.0 ** (decades[first] - 1))
    axis = Axis(
        tuple(p for p, _ in ticks), tuple(values), log=True, tolerance=TICK_TOLERANCE
    )
    axis.check()
    return axis


def _y_axis(segments: list[np.ndarray], figure: Figure) -> Axis:
    """Return the linear ratio axis: the labels on the long ticks, the rest between."""
    ticks = _ticks(segments, figure, 1)
    left = figure.frame[0]
    words = texts(
        DIRECTORY, PAGE, (left - 20, figure.frame[1] - 5, left, figure.frame[3] + 5)
    )
    majors: dict[float, float] = {}
    for position, long in ticks:
        if not long:
            continue
        near = [
            w for w, _, y in words if w in figure.y_labels and abs(y - position) < 2
        ]
        if len(near) != 1:
            msg = f"{figure.name}: {len(near)} labels at the long y tick {position}"
            raise ValueError(msg)
        majors[position] = float(near[0])
    if sorted(majors.values()) != [float(v) for v in figure.y_labels]:
        msg = f"{figure.name}: y labels {sorted(majors.values())}"
        raise ValueError(msg)
    keys = sorted(majors)
    between = [p for p, _ in ticks if keys[0] < p < keys[1]]
    step = (majors[keys[1]] - majors[keys[0]]) / (len(between) + 1)
    values = []
    for position, _ in ticks:
        nearest = min(keys, key=lambda k: abs(k - position))
        rank = sum(
            1 for p, _ in ticks if min(nearest, position) < p <= max(nearest, position)
        )
        sign = 1 if position >= nearest else -1
        values.append(majors[nearest] + sign * rank * step)
    axis = Axis(tuple(p for p, _ in ticks), tuple(values), tolerance=TICK_TOLERANCE)
    axis.check()
    return axis


@dataclass(frozen=True)
class Marker:
    """A marker's centre and the ends of its bars (pt; None where not drawn)."""

    x: float
    y: float
    half: float
    top: float | None
    bottom: float | None


def _markers(drawn: list[Subpath], figure: Figure) -> list[Marker]:
    """Return a figure's open squares inside the frame, with their error bars."""
    squares = []
    for s in drawn:
        if not (s.stroked and not s.filled and s.closed and not s.curved):
            continue
        if len(s.vertices) != 5:
            continue
        xs, ys = s.vertices[:, 0], s.vertices[:, 1]
        width, height = float(np.ptp(xs)), float(np.ptp(ys))
        if not (SQUARE[0] <= width <= SQUARE[1] and abs(width - height) < 0.02):
            continue
        centre = (float(xs.min() + xs.max()) / 2, float(ys.min() + ys.max()) / 2)
        if _inside(*centre, figure.frame) and not _inside(*centre, figure.legend):
            squares.append((centre, height / 2))
    segments = _segments(drawn)
    markers = []
    for (x, y), half in sorted(set(squares)):
        ends: dict[int, float] = {}
        for s in segments:
            if (
                abs(s[0, 0] - s[1, 0]) > ON_MARKER
                or abs(s[:, 0].mean() - x) > ON_MARKER
            ):
                continue
            near, far = sorted(s[:, 1], key=lambda v: abs(v - y))
            if abs(near - y) > half + ON_MARKER or abs(far - y) <= half:
                continue
            caps = [
                c
                for c in segments
                if abs(c[0, 1] - far) < ON_MARKER
                and abs(c[1, 1] - far) < ON_MARKER
                and min(c[:, 0]) < x < max(c[:, 0])
            ]
            if not caps:
                msg = f"{figure.name}: the bar of ({x:.2f}, {y:.2f}) has no cap"
                raise ValueError(msg)
            direction = 1 if far > y else -1
            if direction in ends:
                msg = f"{figure.name}: two bars at ({x:.2f}, {y:.2f})"
                raise ValueError(msg)
            ends[direction] = float(far)
        if len(ends) == 2 and abs(ends[1] + ends[-1] - 2 * y) > BAR_SYMMETRY:
            msg = f"{figure.name}: the bars of ({x:.2f}, {y:.2f}) are not symmetric"
            raise ValueError(msg)
        markers.append(Marker(x, y, half, ends.get(1), ends.get(-1)))
    return markers


def _points(markers: list[Marker], x_axis: Axis, y_axis: Axis) -> list[Point]:
    """Return the markers' intensities, ratios and sigmas, in order of intensity."""
    read = []
    for m in markers:
        intensity = float(x_axis(m.x))
        ratio = float(y_axis(m.y))
        where = f"marker at ({m.x:.2f}, {m.y:.2f}) pt"
        if m.top is not None and m.bottom is not None:
            sigma = f"{(float(y_axis(m.top)) - float(y_axis(m.bottom))) / 2:.3f}"
            printed = f"{where}; bar ends {m.top:.2f} and {m.bottom:.2f}"
        elif m.top is None and m.bottom is None:
            sigma, printed = "", f"{where}; no bar"
        else:
            msg = f"a bar on one side only at ({m.x:.2f}, {m.y:.2f})"
            raise ValueError(msg)
        read.append((intensity, ratio, sigma, printed))
    read.sort()
    return [
        Point(str(k), f"{i:.4e}", f"{r:.3f}", s, p)
        for k, (i, r, s, p) in enumerate(read, start=1)
    ]


def readings() -> dict[str, tuple[list[Marker], Axis, Axis]]:
    """Return each figure's markers and calibrated axes."""
    out = {}
    for name, figure in (("fig1", FIG1), ("fig2", FIG2)):
        drawn = subpaths(DIRECTORY, PAGE, figure.frame)
        segments = _segments(drawn)
        out[name] = (
            _markers(drawn, figure),
            _x_axis(segments, figure),
            _y_axis(segments, figure),
        )
    return out


def extract() -> dict[str, list[Point]]:
    """Return the data files: O2:Xe (Fig. 1) and N2:Ar (Fig. 2) against intensity."""
    read = readings()
    o2 = _points(*read["fig1"])
    n2 = _points(*read["fig2"])
    if len(o2) != 18 or sum(1 for p in o2 if p.sigma) != 1:
        msg = f"Fig. 1: {len(o2)} O2 markers, expected 18 (one with bars)"
        raise ValueError(msg)
    if len(n2) != 34 or not all(p.sigma for p in n2):
        msg = f"Fig. 2: {len(n2)} N2 markers, expected 34, all with bars"
        raise ValueError(msg)
    checks(o2, n2)
    return {
        f"O2/{SOURCE}__fig1_O2_Xe.csv": o2,
        f"N2/{SOURCE}__fig2_N2_Ar.csv": n2,
    }


# The text's numbers, and how far the points may stray from them (the record's notes):
# * the O2:Xe ratio 'At low intensities ... remains nearly constant at a value between
#   0.15 and 0.2' (p. 2), below the 1e14 W/cm^2 above which 'saturation of the Xe yield
#   results in an increase of the ionization ratio': the median of the points below
#   1e14 W/cm^2 must lie in that range, and every point within O2_SCATTER of it;
# * the N2:Ar ratio 'nearly constant at approximately 1.7, before falling towards
#   unity near the saturation intensity' (p. 3): the median within N2_SPREAD of 1.7;
# * 'shifting the intensity scale of our N2+ data by 10% changes the N2:Ar ionization
#   ratio at 1 x 10^14 W cm^-2 from 1.90 to 0.93' (p. 1): the ratio interpolated (in
#   log intensity) at 1e14 W/cm^2 within N2_AT_1E14_TOLERANCE of 1.90, a quarter of
#   the bars' half-length there (0.18).
O2_RANGE = (0.15, 0.2)
O2_SCATTER = 0.05
O2_BELOW = 1e14
N2_PLATEAU = 1.7
N2_SPREAD = 0.2
N2_AT_1E14 = 1.90
N2_AT_1E14_TOLERANCE = 0.05
# The bars are symmetric: their halves agree within the drawing's quantisation, pt.
BAR_SYMMETRY = 0.15


def n2_at_1e14(n2: list[Point]) -> float:
    """Return the N2:Ar ratio interpolated in log intensity at 1e14 W/cm^2."""
    xs = np.log10([float(p.x) for p in n2])
    return float(np.interp(14.0, xs, [float(p.value) for p in n2]))


def checks(o2: list[Point], n2: list[Point]) -> None:
    """Raise unless the points agree with the numbers the text prints."""
    low = [float(p.value) for p in o2 if float(p.x) < O2_BELOW]
    median = float(np.median(low))
    if not O2_RANGE[0] <= median <= O2_RANGE[1]:
        msg = f"Fig. 1: the median O2:Xe below 1e14 W/cm2 is {median:.3f}"
        raise ValueError(msg)
    if not all(O2_RANGE[0] - O2_SCATTER <= v <= O2_RANGE[1] + O2_SCATTER for v in low):
        msg = f"Fig. 1: O2:Xe below 1e14 W/cm2 spans {min(low)}-{max(low)}"
        raise ValueError(msg)
    median = float(np.median([float(p.value) for p in n2]))
    if not math.isclose(median, N2_PLATEAU, abs_tol=N2_SPREAD):
        msg = f"Fig. 2: the N2:Ar median {median:.2f} is not about 1.7"
        raise ValueError(msg)
    at = n2_at_1e14(n2)
    if not math.isclose(at, N2_AT_1E14, abs_tol=N2_AT_1E14_TOLERANCE):
        msg = f"Fig. 2: N2:Ar at 1e14 W/cm2 is {at:.3f}, the text 1.90"
        raise ValueError(msg)


def main() -> None:
    """Write the data files and print what the checks compare."""
    files = extract()
    for name, points in files.items():
        values = [float(p.value) for p in points]
        xs = [float(p.x) for p in points]
        print(
            f"{name}: {len(points)} points, {min(xs):.3e}-{max(xs):.3e} W/cm2, "
            f"ratio {min(values):.3f}-{max(values):.3f}, "
            f"median {float(np.median(values)):.3f}"
        )
    write(files)


if __name__ == "__main__":
    main()
