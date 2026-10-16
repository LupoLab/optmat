# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Loriot, Hertz, Lavorel & Faucher (2008): N2 ionisation, N2+:Ar+, 800 nm.

The paper prints the N2:Ar ratio as one number (0.8 +- 0.1) and the O2:Xe ratio as a
bound, in the text (the record's inline datasets); the points are in two vector
figures, read from the drawing (data rule 8, TC-699). The PDF's page 1 is the
publisher's cover sheet: the article's page n is PDF page n + 1.

**Fig. 4** (article p. 4, PDF page 5, library figure ``fig_bea52e231748``): the
ionisation probability of N2 against the peak intensity (TW/cm^2), both axes
logarithmic: the full circles ('parallel (from our previous work [24])', Loriot et al.
2006) and the blue triangles ('perpendicular pump and probe polarization'). The
semi-empirical curve and the vertical lines bounding regions (A) and (B) are not
data; the lines' intensities are read for the record's notes.

**Fig. 7** (article p. 5, PDF page 6, library figure ``fig_0f8df2057eca``): the
N2+:Ar+ ratio against the peak intensity (TW/cm^2), linear axes: the black circles
('the present technique (800 nm, 100 fs)'). Not taken: DeWitt et al.'s points (blue
squares, curated from the original as dewitt_2001; read here only to compare,
:func:`compare_dewitt`) and Guo et al.'s (red triangles).

How the drawing is read:

* the paths: every path object of the page split into its subpaths (at each
  move-to), with the object's fill colour, in page points.
* the axes: the frames' tick marks are stroked two-point subpaths across the bottom
  (x) or left (y) edge, one end on it, under 4 pt long; long ticks (over 1.8 pt: they
  are 2.4-2.5 pt, the short ones 1.2 pt) carry the labels. The labels are drawn as
  glyph outlines, not text, so their values are written here (``AxisSpec.long``) as
  printed, in order, and the count of long ticks must match. The short ticks divide
  the step between two long ones evenly in value (Fig. 4's intensity axis has eight
  between 20 and 40 TW/cm^2: steps of 20/9; its probability axis eight per decade,
  2 to 9 x 10^n); beyond the last long tick they continue with ``AxisSpec.above``.
  Every tick must lie on the fitted scale within ``TICK_TOLERANCE``.
* the markers: Fig. 4's and Fig. 7's circles are closed curved subpaths filled black,
  3.0-3.2 pt across (their centre: the middle of their vertices' extent); Fig. 4's
  triangles closed three-sided subpaths filled blue (0, 0, 255), pointing up; their
  centre is taken at the middle of their vertical extent (the centroid lies a sixth
  of the height, 0.51 pt, lower: the digitising uncertainty covers half of that).
  Markers outside the frame are dropped.

Every ``printed`` cell gives the marker's centre in page points. The rows are numbered
in order of intensity.
"""

import csv
import ctypes
import math
from dataclasses import dataclass

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, pdf
from _library import DATA, Point, write
from pypdfium2 import raw

DIRECTORY = "Loriot_2007_FieldFreeMolecular"
SOURCE = "loriot_2008"
FIG4_PAGE, FIG7_PAGE = 5, 6
TICK_TOLERANCE = 0.05  # pt: the worst tick is 0.042 pt off (Fig. 4, x)
TICK_LONGEST = 4.0
TICK_LONG = 1.8
ON_EDGE = 0.02
BLACK, BLUE = (0, 0, 0), (0, 0, 255)
CIRCLE = (2.9, 3.3)  # the circles' vertices span 3.06 pt (Fig. 4) and 3.0 (Fig. 7)
TRIANGLE = (3.3, 3.8)  # the triangles' width, 3.56 pt
DEWITT_SQUARE = (2.4, 2.7)  # Fig. 7's filled blue squares, 2.58 pt a side


@dataclass(frozen=True)
class Subpath:
    """A subpath: its vertices (pt) and how its object is painted."""

    vertices: np.ndarray
    closed: bool
    curved: bool
    stroked: bool
    fill: tuple[int, int, int] | None


def _fill(obj: pdfium.PdfObject) -> tuple[int, int, int] | None:
    """Return a path object's fill colour, or None if it is not filled."""
    fill, stroke = ctypes.c_int(), ctypes.c_int()
    raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill), ctypes.byref(stroke))
    if not fill.value:
        return None
    channels = [ctypes.c_uint() for _ in range(4)]
    raw.FPDFPageObj_GetFillColor(obj.raw, *(ctypes.byref(c) for c in channels))
    return channels[0].value, channels[1].value, channels[2].value


def subpaths(page: int, region: tuple[float, float, float, float]) -> list[Subpath]:
    """Return the subpaths of the path objects touching a region (l, b, r, t, pt)."""
    left, bottom, right, top = region
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    found: list[Subpath] = []
    for obj in document[page - 1].get_objects(max_depth=15):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        if obj.container is not None:  # this paper's figures are drawn on the page
            msg = f"p. {page}: a path inside a form, which this script does not map"
            raise ValueError(msg)
        west, south, east, north = obj.get_bounds()
        if not (east >= left and west <= right and north >= bottom and south <= top):
            continue
        fill, stroke = ctypes.c_int(), ctypes.c_int()
        raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill), ctypes.byref(stroke))
        matrix = obj.get_matrix()
        groups: list[list[tuple[float, float, int, bool]]] = []
        for index in range(raw.FPDFPath_CountSegments(obj.raw)):
            segment = raw.FPDFPath_GetPathSegment(obj.raw, index)
            x, y = ctypes.c_float(), ctypes.c_float()
            raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
            kind = raw.FPDFPathSegment_GetType(segment)
            if kind == raw.FPDF_SEGMENT_MOVETO or not groups:
                groups.append([])
            x_page, y_page = matrix.on_point(x.value, y.value)
            closes = bool(raw.FPDFPathSegment_GetClose(segment))
            groups[-1].append((x_page, y_page, kind, closes))
        for group in groups:
            found.append(
                Subpath(
                    np.array([(g[0], g[1]) for g in group], dtype=float),
                    closed=group[-1][3],
                    curved=any(g[2] == raw.FPDF_SEGMENT_BEZIERTO for g in group),
                    stroked=bool(stroke.value),
                    fill=_fill(obj),
                )
            )
    return found


@dataclass(frozen=True)
class AxisSpec:
    """An axis's printed labels and how its short ticks continue beyond them.

    Attributes
    ----------
    long : tuple of float
        The values printed at the long ticks, in order (the glyphs are outlines).
    above : float
        The step of the short ticks beyond the last long tick.
    below : float
        The step of the short ticks before the first long tick.
    log : bool
        True for a logarithmic axis.
    """

    long: tuple[float, ...]
    above: float
    below: float
    log: bool


@dataclass(frozen=True)
class Figure:
    """A figure: its page, frame (l, b, r, t, pt) and axes."""

    name: str
    page: int
    frame: tuple[float, float, float, float]
    x: AxisSpec
    y: AxisSpec


FIG4 = Figure(
    "Fig. 4",
    FIG4_PAGE,
    (78.03, 591.29, 292.62, 733.25),
    AxisSpec((20.0, 40.0, 60.0, 80.0, 100.0, 120.0), 20 / 9, 20 / 9, log=True),
    AxisSpec(tuple(10.0**n for n in range(-7, -1)), 1e-2, 1e-8, log=True),
)
FIG7 = Figure(
    "Fig. 7",
    FIG7_PAGE,
    (333.45, 660.97, 549.14, 733.41),
    AxisSpec((40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0), 5.0, 5.0, log=False),
    AxisSpec((0.0, 0.5, 1.0, 1.5, 2.0), 0.25, 0.25, log=False),
)


def _ticks(drawn: list[Subpath], figure: Figure, axis: int) -> list[tuple[float, bool]]:
    """Return the ticks across the bottom (0) or left (1) edge: (position, long)."""
    left, bottom, right, top = figure.frame
    edge = bottom if axis == 0 else left
    along, across = (0, 1) if axis == 0 else (1, 0)
    low, high = (left, right) if axis == 0 else (bottom, top)
    ticks: dict[float, bool] = {}
    for s in drawn:
        if not s.stroked or s.curved or s.closed or len(s.vertices) != 2:
            continue
        v = s.vertices
        if abs(v[0, along] - v[1, along]) > ON_EDGE:
            continue
        length = abs(v[0, across] - v[1, across])
        on_edge = min(abs(v[0, across] - edge), abs(v[1, across] - edge)) < ON_EDGE
        if not on_edge or not 0 < length < TICK_LONGEST:
            continue
        position = round(float(v[:, along].mean()), 3)
        if low - ON_EDGE <= position <= high + ON_EDGE:
            ticks[position] = ticks.get(position, False) or length > TICK_LONG
    return sorted(ticks.items())


def _axis(ticks: list[tuple[float, bool]], spec: AxisSpec, where: str) -> Axis:
    """Return an axis: the long ticks' printed values, the short ones evenly between."""
    longs = [p for p, is_long in ticks if is_long]
    if len(longs) != len(spec.long):
        msg = f"{where}: {len(longs)} long ticks, {len(spec.long)} labels"
        raise ValueError(msg)
    value_at = dict(zip(longs, spec.long, strict=True))
    positions = [p for p, _ in ticks]
    values = []
    for position in positions:
        if position in value_at:
            values.append(value_at[position])
            continue
        lower = [p for p in longs if p < position]
        upper = [p for p in longs if p > position]
        if lower and upper:
            a, b = max(lower), min(upper)
            between = [p for p in positions if a < p < b]
            rank = between.index(position) + 1
            values.append(
                value_at[a] + (value_at[b] - value_at[a]) * rank / (len(between) + 1)
            )
        elif lower:
            rank = sum(1 for p in positions if max(lower) < p <= position)
            values.append(value_at[max(lower)] + rank * spec.above)
        else:
            rank = sum(1 for p in positions if position <= p < min(upper))
            values.append(value_at[min(upper)] - rank * spec.below)
    axis = Axis(tuple(positions), tuple(values), log=spec.log, tolerance=TICK_TOLERANCE)
    axis.check()
    return axis


def axes(figure: Figure) -> tuple[Axis, Axis, list[Subpath]]:
    """Return a figure's calibrated x and y axes and its subpaths."""
    drawn = subpaths(figure.page, figure.frame)
    x_axis = _axis(_ticks(drawn, figure, 0), figure.x, f"{figure.name} x")
    y_axis = _axis(_ticks(drawn, figure, 1), figure.y, f"{figure.name} y")
    return x_axis, y_axis, drawn


def _inside(x: float, y: float, frame: tuple[float, float, float, float]) -> bool:
    return frame[0] <= x <= frame[2] and frame[1] <= y <= frame[3]


def _centres(
    drawn: list[Subpath],
    figure: Figure,
    fill: tuple[int, int, int],
    curved: bool,
    span: tuple[float, float],
    corners: int,
) -> list[tuple[float, float]]:
    """Return the centres of a kind of marker inside the frame (pt), each once."""
    found = set()
    for s in drawn:
        if s.fill != fill or s.curved != curved or not s.closed:
            continue
        if not curved and len({tuple(np.round(v, 2)) for v in s.vertices}) != corners:
            continue
        xs, ys = s.vertices[:, 0], s.vertices[:, 1]
        if not span[0] <= float(np.ptp(xs)) <= span[1]:
            continue
        centre = (
            round(float(xs.min() + xs.max()) / 2, 3),
            round(float(ys.min() + ys.max()) / 2, 3),
        )
        if _inside(*centre, figure.frame):
            found.add(centre)
    return sorted(found)


def _points(
    centres: list[tuple[float, float]], x_axis: Axis, y_axis: Axis, value: str
) -> list[Point]:
    """Return markers' intensities (W/cm^2, from TW/cm^2) and values, in order."""
    read = sorted(
        (float(x_axis(x)) * 1e12, float(y_axis(y)), f"marker at ({x:.2f}, {y:.2f}) pt")
        for x, y in centres
    )
    return [
        Point(str(k), f"{i:.4e}", format(v, value), "", p)
        for k, (i, v, p) in enumerate(read, start=1)
    ]


def region_lines(drawn: list[Subpath], x_axis: Axis) -> list[float]:
    """Return the intensities (TW/cm^2) of Fig. 4's vertical lines across the plot."""
    left, bottom, right, top = FIG4.frame
    found = []
    for s in drawn:
        if not s.stroked or s.fill is not None or len(s.vertices) != 2:
            continue
        (x0, y0), (x1, y1) = s.vertices
        across = abs(x0 - x1) < ON_EDGE and abs(y1 - y0) > 0.9 * (top - bottom)
        if across and left + 1 < x0 < right - 1:
            found.append(round(float(x_axis(x0)), 2))
    return sorted(found)


def extract() -> dict[str, list[Point]]:
    """Return the data files: N2's probability (Fig. 4) and N2+:Ar+ (Fig. 7)."""
    x4, y4, drawn4 = axes(FIG4)
    circles = _centres(drawn4, FIG4, BLACK, curved=True, span=CIRCLE, corners=0)
    triangles = _centres(drawn4, FIG4, BLUE, curved=False, span=TRIANGLE, corners=3)
    x7, y7, drawn7 = axes(FIG7)
    ratios = _centres(drawn7, FIG7, BLACK, curved=True, span=CIRCLE, corners=0)
    if (len(circles), len(triangles), len(ratios)) != (12, 6, 7):
        msg = f"markers: {len(circles)}, {len(triangles)}, {len(ratios)}"
        raise ValueError(msg)
    return {
        f"N2/{SOURCE}__fig4_parallel.csv": _points(circles, x4, y4, ".3e"),
        f"N2/{SOURCE}__fig4_perpendicular.csv": _points(triangles, x4, y4, ".3e"),
        f"N2/{SOURCE}__fig7_N2_Ar.csv": _points(ratios, x7, y7, ".3f"),
    }


# Loriot et al. redraw DeWitt et al.'s N2:Ar points in Fig. 7; each must match one of
# the points read from DeWitt's own Fig. 2 (dewitt_2001) within these: 1 % in
# intensity and 0.02 in the ratio, about the two drawings' precision (0.06 pt of 22.6
# pt per 0.5 in DeWitt's Fig. 2, 0.01 pt of 30.2 pt per 0.5 here; the redrawing's own).
MATCH_INTENSITY = 0.01
MATCH_RATIO = 0.02


def compare_dewitt() -> list[tuple[float, float, float, float]]:
    """Return Fig. 7's DeWitt points matched to dewitt_2001's: (I, R here, I, R there).

    Raises unless every redrawn point matches one of DeWitt's own.
    """
    x7, y7, drawn7 = axes(FIG7)
    squares = _centres(drawn7, FIG7, BLUE, curved=False, span=DEWITT_SQUARE, corners=4)
    path = DATA / "values" / "ionisation" / "N2" / "dewitt_2001__fig2_N2_Ar.csv"
    with path.open(encoding="utf-8") as handle:
        theirs = [(float(r["x"]), float(r["value"])) for r in csv.DictReader(handle)]
    pairs = []
    for x, y in squares:
        intensity, ratio = float(x7(x)) * 1e12, float(y7(y))
        nearest = min(theirs, key=lambda t: abs(math.log(t[0] / intensity)))
        if (
            abs(nearest[0] / intensity - 1) > MATCH_INTENSITY
            or abs(nearest[1] - ratio) > MATCH_RATIO
        ):
            msg = f"Fig. 7: DeWitt's point at {intensity:.3e}, {ratio:.3f} unmatched"
            raise ValueError(msg)
        pairs.append((intensity, ratio, nearest[0], nearest[1]))
    return pairs


def main() -> None:
    """Write the data files and print the checks the record quotes."""
    files = extract()
    for name, points in files.items():
        xs = [float(p.x) for p in points]
        values = [float(p.value) for p in points]
        print(
            f"{name}: {len(points)} points, {min(xs):.3e}-{max(xs):.3e} W/cm2, "
            f"{min(values):.3g}-{max(values):.3g}"
        )
    x4, _, drawn4 = axes(FIG4)
    lines = region_lines(drawn4, x4)
    print("Fig. 4 region lines (TW/cm2):", lines)
    # The text's nonlinearity N_L = 7.5 +- 0.5 (from the fitted F, ref. [24]): the
    # slope of log P against log I of the parallel points inside region (A).
    parallel = files[f"N2/{SOURCE}__fig4_parallel.csv"]
    inside = [p for p in parallel if lines[1] <= float(p.x) / 1e12 <= lines[2]]
    slope = np.polyfit(
        np.log([float(p.x) for p in inside]),
        np.log([float(p.value) for p in inside]),
        1,
    )[0]
    print(f"Fig. 4 parallel, region (A): {len(inside)} points, slope {slope:.2f}")
    ratios = [float(p.value) for p in files[f"N2/{SOURCE}__fig7_N2_Ar.csv"]]
    print(f"Fig. 7 N2:Ar mean {np.mean(ratios):.3f}, sd {np.std(ratios, ddof=1):.3f}")
    pairs = compare_dewitt()
    worst_i = max(abs(a / c - 1) for a, _, c, _ in pairs)
    worst_r = max(abs(b - d) for _, b, _, d in pairs)
    print(
        f"DeWitt's points redrawn in Fig. 7: {len(pairs)}, worst {worst_i:.2%} in "
        f"intensity, {worst_r:.3f} in the ratio"
    )
    write(files)


if __name__ == "__main__":
    main()
