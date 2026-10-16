# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Le et al. (2016): one-electron TDSE ionisation of atomic hydrogen, Fig. 12.

The tutorial prints no table of its calculations: the time-dependent Schrodinger
equation (TDSE) probabilities of hydrogen are given only in **Fig. 12** (PDF page 19,
printed p. 18, library figure ``fig_d7530dc1bd66``), at 400 nm (a), 800 nm (b) and
1600 nm (c), against the laser intensity, beside the ADK, the PPT (Eq. 91) and the SFA
and SFA-CC, the last three multiplied by the factors the legends print (section 2.5.3).
Only the TDSE curves are read: the others are models.

The figure is a vector drawing on the page (no image, no form) and is read from the
drawing (data rule 8, TC-699). Its labels are glyphs drawn as filled paths, not text.

* the frames: a stroked black rectangle per panel (:data:`PANELS` gives its lower left
  corner), its edges drawn again as lines; the ticks separate stroked black two-vertex
  lines perpendicular to an edge with one end on it. The top and right edges repeat the
  bottom and left ones' ticks (checked: every tick there has one at the same position
  below or left, within :data:`MIRROR`).
* the axes: x linear, long ticks (4.2 pt) at the labels '0.2', '0.4', ..., '1.4' (in
  10^14 W/cm^2, the axis title) and short ones (2.1 pt) between and at the ends, so
  every tick 0.1 x 10^14 from the next: (a) and (b) 0.1 ... 1.5, (c) 0.2 ... 1.5. y
  logarithmic, every tick long and labelled, two decades apart: (a) '10^-14' ...
  '10^0', (b) '10^-14' ... '10^-2', (c) '10^-10' ... '10^0'. Every tick lies on its
  fitted scale within :data:`TICK_TOLERANCE` (:meth:`_figure.Axis.check`); the y labels'
  glyphs sit one label at each tick, the x labels' glyphs span the long ticks
  (:func:`_check_labels`).
* the curves: each panel's TDSE is the one stroked black path of over ten segments
  (the legend's first entry, 'TDSE', a black line sample, read on the page; the
  samples' colours in the legends' order are checked, :func:`_check_legend`). The
  other curves are polylines; the TDSE paths alone mix straight segments with cubic
  Bezier segments: (a) 2 cubics (from 1.2 to 7.9 and from 8.1 to 14.8 x 10^13 W/cm^2)
  and 5 straight segments; (b) 8 cubics, each over 13 or 14 steps of 10^12 W/cm^2,
  and 30 straight segments; (c) 3 cubics of 9 steps and 103 straight segments. Only the
  vertices that the drawn curve passes through are read (the ends of the straight
  segments and of the cubics; :func:`curve`), not the cubics' control points: these
  are the data, the cubics between them a drawing's interpolation. Each read vertex
  lies on a multiple of 10^12 W/cm^2, within :data:`GRID` (the SFA curve is drawn at
  every multiple of 10^12, 141 vertices in (a) and (b), 132 in (c): the calculations'
  grid). Where the path's end is drawn twice at the frame's edge (a vertex within
  :data:`DUPLICATE` of the next), the one nearer the grid is kept.

Every point's ``printed`` cell gives the panel, whether the vertex ends a straight
segment or a cubic, and its position in page points. The rows are numbered in order of
intensity.

**Checks** (:func:`checks`), against every number the paper prints for this figure:

* the vertical dotted lines, 'the position with gamma = 1' (caption), lie at the
  intensity where hydrogen's gamma = omega sqrt(2 I_p)/F (I_p = 0.5 hartree) is 1, with
  the intensity c eps0 F^2/2, within :data:`GAMMA`; (a) has none (gamma = 1 at 4.6 x
  10^14 W/cm^2 is beyond its axis).
* the TDSE curves of (a) and (b) span the text's '1 x 10^13 to 1.5 x 10^14 W cm^-2'
  and (c) its axis's 2 x 10^13 to 1.5 x 10^14, within :data:`GRID`.
The text's numbers that are not checks are reported (:func:`report`), for the record's
notes: 'at 4 x 10^13 W cm^-2 the ADK underestimates by almost two orders of magnitude'
(800 nm; the drawn ADK over the drawn TDSE, :func:`adk_ratio`), and the PPT curves over
the TDSE with each legend's factor taken out.
"""

import ctypes
import functools
import math
from dataclasses import dataclass

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, Painted, painted, pdf
from _library import Point, write
from pypdfium2 import raw

from optmat.constants import AU_ANGULAR_FREQUENCY, AU_FIELD, EPS0, C

DIRECTORY = "Le_2016_StrongFieldApproximation"
SOURCE = "le_2016"
PAGE = 19  # the PDF's page; printed p. 18 (the PDF has a cover page)
BLACK, RED, GREEN, BLUE, MAGENTA = (
    (0, 0, 0),
    (255, 0, 0),
    (0, 255, 0),
    (0, 0, 255),
    (255, 0, 255),
)
# The drawing's coordinates are quantised in steps of about 0.02 pt (the x ticks, 15.09
# to 15.11 pt apart for equal steps); every tick lies on its fitted scale within 0.009
# pt (the largest residual, the x axes of (a) and (b)). 0.03 pt allows for it.
TICK_TOLERANCE = 0.03
STRAIGHT = 0.01  # a tick is vertical or horizontal to 0.01 pt
ON_LINE = 0.01  # a tick ends on its edge within 0.01 pt (0.002 found)
TICK_LONGEST = 5.0  # the ticks are 2.1 and 4.2 pt long; the frames' edges 140 or more
LONG_TICK = 1.5  # a long tick is over 1.5 times the shortest (2.0 times)
MIRROR = 0.01  # the top and right edges' ticks are the bottom and left ones' to 0.002
# The labels' glyphs: the y labels' wholly left of the frame, beginning within 18 pt of
# it (the axis title begins 28 pt left); a label's glyphs ('10' and its exponent)
# overlap in height, its neighbours are 15 pt or more apart; its middle lies within
# 0.3 pt of its tick (1.0 allowed). The x labels' glyphs are centred within 10 pt below
# the frame (the axis title 12 pt or more).
LABEL_LEFT = 18.0
LABEL_GAP = 1.0
LABEL_AT = 1.0
LABEL_BELOW = 10.0
X_LABEL_AT = 3.0  # the x labels' extent is centred on the long ticks' within 3 pt
# A read vertex lies within 0.03 x 10^12 W/cm^2 of a multiple of 10^12 (0.045 pt; the
# largest offset found 0.013, the drawing's quantum); the curves' ends are drawn twice,
# 0.021 pt apart, at the frame's edge.
STEP = 1e12
GRID = 0.03
DUPLICATE = 0.05
# The gamma = 1 lines lie within 1 % of hydrogen's gamma = 1 intensity (0.3 % found;
# the axis's quantum 0.02 pt is 0.013e12 W/cm^2, 0.05 % of 2.85e13).
GAMMA = 0.01
# 'at 4 x 10^13 W cm^-2 the ADK underestimates by almost two orders of magnitude'
# (reported, not checked: the drawing gives 1/17 there, 1/78 at 3 x 10^13).
ADK_AT = 4e13
I_AU = 0.5 * C * EPS0 * AU_FIELD**2 / 1e4  # W/cm^2 per atomic unit of field squared
H_IP_AU = 0.5  # hydrogen's ionisation energy in a one-electron TDSE (infinite mass)
# A curve belongs to the panel whose frame holds its vertices within 1 pt (the ADK of
# (b) ends 0.35 pt below the frame; the panels are 30 pt or more apart).
CURVE_MARGIN = 1.0
SAMPLES = 64  # points per cubic when the drawn curve is evaluated (checks, report)


@dataclass(frozen=True)
class Panel:
    """A panel: name, frame's lower left corner (pt, to 1 pt), axes, pulse, factors.

    ``x_labels``: the long x ticks' values (W/cm^2); every tick is 0.1 x 10^14 from
    the next. ``y_labels``: the y ticks' decades, every tick labelled. ``wavelength``
    (nm) and ``cycles``: the text's; ``ppt``: the legend's 'PPT (x...)' factor, read on
    the page.
    """

    name: str
    corner: tuple[float, float]
    x_labels: tuple[float, ...]
    y_labels: tuple[int, ...]
    wavelength: int
    cycles: int
    ppt: float


X_LABELS = tuple(k * 2e13 for k in range(1, 8))  # '0.2' ... '1.4' x 10^14 W/cm^2
X_STEP = 1e13
PANELS = {
    "a": Panel(
        "Fig. 12(a)", (78.9, 620.4), X_LABELS, tuple(range(-14, 1, 2)), 400, 12, 0.33
    ),
    "b": Panel(
        "Fig. 12(b)", (78.9, 439.0), X_LABELS, tuple(range(-14, -1, 2)), 800, 6, 0.33
    ),
    "c": Panel(
        "Fig. 12(c)", (78.9, 256.1), X_LABELS, tuple(range(-10, 1, 2)), 1600, 6, 0.5
    ),
}
# The legends' entries from the top, read on the page: 'TDSE', 'ADK', 'PPT (x0.33)',
# 'SFA (x12)', 'SFA-CC (x0.051)' in (a); x0.33, x15, x0.064 in (b); x0.5, x20, x0.073
# in (c).
LEGEND = (BLACK, RED, GREEN, BLUE, MAGENTA)
# The intensities the TDSE curves span: the text's 1e13 to 1.5e14 W/cm^2 in (a) and
# (b); (c)'s axis begins at 2e13.
SPAN = {"a": (1e13, 1.5e14), "b": (1e13, 1.5e14), "c": (2e13, 1.5e14)}
DATASETS = {
    p: f"H/{SOURCE}__fig12{p}_tdse_H_{panel.wavelength}nm.csv"
    for p, panel in PANELS.items()
}


@dataclass(frozen=True)
class Frame:
    """A panel's frame (l, b, r, t, pt), its axes, and the page's paths around it."""

    box: tuple[float, float, float, float]
    x: Axis
    y: Axis
    x_long: tuple[float, ...]
    y_ticks: tuple[float, ...]
    drawn: tuple[Painted, ...]


def _box(drawn: list[Painted], panel: Panel) -> tuple[float, float, float, float]:
    """Return the frame: the stroked black rectangle whose lower left is the corner.

    From its vertices (the bounding box includes half the stroke's width).
    """
    x0, y0 = panel.corner
    found = [
        (
            float(d.vertices[:, 0].min()),
            float(d.vertices[:, 1].min()),
            float(d.vertices[:, 0].max()),
            float(d.vertices[:, 1].max()),
        )
        for d in drawn
        if d.stroke == BLACK
        and d.fill is None
        and len(d.vertices) == 5
        and abs(float(d.vertices[:, 0].min()) - x0) < 1
        and abs(float(d.vertices[:, 1].min()) - y0) < 1
    ]
    if len(found) != 1:
        msg = f"{panel.name}: {len(found)} frames"
        raise ValueError(msg)
    return found[0]


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
) -> list[tuple[float, float]]:
    """Return an axis's ticks (position, length), each mirrored on the far edge."""
    near = _edge_ticks(drawn, box, axis, far=False)
    far = _edge_ticks(drawn, box, axis, far=True)
    positions = [p for p, _ in near]
    if len(far) != len(near) or any(
        min(abs(p - q) for q in positions) > MIRROR for p, _ in far
    ):
        msg = f"{where}: the opposite edge's ticks are not this one's"
        raise ValueError(msg)
    return near


def _x_axis(ticks: list[tuple[float, float]], panel: Panel) -> tuple[Axis, list[float]]:
    """Return the x axis, every tick 0.1 x 10^14 from the next, and its long ticks."""
    shortest = min(length for _, length in ticks)
    positions = [p for p, _ in ticks]
    long = [p for p, length in ticks if length > LONG_TICK * shortest]
    if len(long) != len(panel.x_labels):
        msg = f"{panel.name} x: {len(long)} long ticks, {len(panel.x_labels)} labels"
        raise ValueError(msg)
    first = positions.index(long[0])
    values = [panel.x_labels[0] + (k - first) * X_STEP for k in range(len(positions))]
    for p, label in zip(long, panel.x_labels, strict=True):
        if not math.isclose(values[positions.index(p)], label, rel_tol=1e-9):
            msg = f"{panel.name} x: the long tick at {p} pt is not labelled {label}"
            raise ValueError(msg)
    axis = Axis(tuple(positions), tuple(values), tolerance=TICK_TOLERANCE)
    axis.check()
    return axis, long


def _y_axis(ticks: list[tuple[float, float]], panel: Panel) -> Axis:
    """Return the y axis, a labelled decade at every tick."""
    if len(ticks) != len(panel.y_labels):
        msg = f"{panel.name} y: {len(ticks)} ticks, {len(panel.y_labels)} labels"
        raise ValueError(msg)
    lengths = [length for _, length in ticks]
    if max(lengths) > LONG_TICK * min(lengths):
        msg = f"{panel.name} y: long and short ticks"
        raise ValueError(msg)
    axis = Axis(
        tuple(p for p, _ in ticks),
        tuple(10.0**e for e in panel.y_labels),
        log=True,
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


def _check_labels(found: Frame, panel: Panel) -> None:
    """Raise unless the y labels sit one at each tick and the x labels span the long.

    The glyphs are filled black paths. A y label is a run of glyphs wholly left of the
    frame, beginning within LABEL_LEFT of it, overlapping in height; its middle lies
    within LABEL_AT of its tick. The x labels' glyphs, wholly below the frame and
    centred within LABEL_BELOW of it, are drawn several to a path: their extent must
    cover the first and last long ticks and be centred on them within X_LABEL_AT.
    """
    left, bottom, right, top = found.box
    glyphs = [d for d in found.drawn if d.fill == BLACK and d.stroke is None]
    beside = [
        (d.bottom, d.top)
        for d in glyphs
        if left - LABEL_LEFT < d.left
        and d.right < left
        and bottom - 5 < d.centre[1] < top + 5
    ]
    middles = [(low + high) / 2 for low, high in _groups(beside, LABEL_GAP)]
    if len(middles) != len(found.y_ticks) or any(
        abs(m - t) > LABEL_AT for m, t in zip(middles, found.y_ticks, strict=True)
    ):
        msg = f"{panel.name}: the y labels at {middles} are not at {found.y_ticks}"
        raise ValueError(msg)
    under = [
        (d.left, d.right)
        for d in glyphs
        if bottom - LABEL_BELOW < d.centre[1]
        and d.top < bottom
        and left - 15 < d.centre[0] < right + 15
    ]
    low, high = min(s[0] for s in under), max(s[1] for s in under)
    ticks = found.x_long
    if not (
        low < ticks[0] < ticks[-1] < high
        and abs((low + high) / 2 - (ticks[0] + ticks[-1]) / 2) <= X_LABEL_AT
    ):
        msg = f"{panel.name}: the x labels ({low}, {high}) do not span {ticks}"
        raise ValueError(msg)


@functools.cache
def frame(key: str) -> Frame:
    """Return a panel's frame and calibrated axes, its labels checked."""
    panel = PANELS[key]
    x0, y0 = panel.corner
    drawn = painted(DIRECTORY, PAGE, (x0 - 40, y0 - 30, x0 + 230, y0 + 160))
    box = _box(drawn, panel)
    x_axis, x_long = _x_axis(_ticks(drawn, box, 0, f"{panel.name} x"), panel)
    y_ticks = _ticks(drawn, box, 1, f"{panel.name} y")
    found = Frame(
        box,
        x_axis,
        _y_axis(y_ticks, panel),
        tuple(x_long),
        tuple(p for p, _ in y_ticks),
        tuple(drawn),
    )
    _check_labels(found, panel)
    return found


def _inside(x: float, y: float, box: tuple, margin: float) -> bool:
    """Return whether a point lies inside a box (l, b, r, t) widened by a margin."""
    return (
        box[0] - margin <= x <= box[2] + margin
        and box[1] - margin <= y <= box[3] + margin
    )


@dataclass(frozen=True)
class Segment:
    """A path's segment: its kind (``move``, ``line``, ``cubic``) and points (pt).

    A cubic's points are its two control points and its end.
    """

    kind: str
    points: tuple[tuple[float, float], ...]


KINDS = {
    raw.FPDF_SEGMENT_MOVETO: "move",
    raw.FPDF_SEGMENT_LINETO: "line",
    raw.FPDF_SEGMENT_BEZIERTO: "cubic",
}


def _stroke(obj: pdfium.PdfObject) -> tuple[int, int, int] | None:
    """Return a path's stroke colour, or None if it is not stroked."""
    fill, stroked = ctypes.c_int(), ctypes.c_int()
    raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill), ctypes.byref(stroked))
    channels = [ctypes.c_uint() for _ in range(4)]
    if not stroked.value or not raw.FPDFPageObj_GetStrokeColor(
        obj.raw, *(ctypes.byref(c) for c in channels)
    ):
        return None
    return channels[0].value, channels[1].value, channels[2].value


@functools.cache
def _paths(colour: tuple[int, int, int]) -> tuple[tuple[Segment, ...], ...]:
    """Return the page's stroked paths of a colour with over ten segments, by segment.

    The figure is drawn on the page (no form; checked); each path is one subpath.
    """
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    found = []
    for obj in document[PAGE - 1].get_objects(max_depth=15):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        count = raw.FPDFPath_CountSegments(obj.raw)
        if count <= 10 or _stroke(obj) != colour:
            continue
        if obj.container is not None:
            msg = f"p. {PAGE}: a curve inside a form, which this script does not map"
            raise ValueError(msg)
        matrix = obj.get_matrix()
        kinds, points = [], []
        for index in range(count):
            segment = raw.FPDFPath_GetPathSegment(obj.raw, index)
            x, y = ctypes.c_float(), ctypes.c_float()
            raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
            kinds.append(KINDS[raw.FPDFPathSegment_GetType(segment)])
            points.append(matrix.on_point(x.value, y.value))
        if kinds[0] != "move" or "move" in kinds[1:]:
            msg = f"p. {PAGE}: a curve of several subpaths"
            raise ValueError(msg)
        segments, k = [], 0
        while k < count:
            # pdfium gives a cubic as three BEZIERTO points: two controls, the end.
            width = 3 if kinds[k] == "cubic" else 1
            if kinds[k : k + width] != [kinds[k]] * width:
                msg = f"p. {PAGE}: a cubic of fewer than three points"
                raise ValueError(msg)
            segments.append(Segment(kinds[k], tuple(points[k : k + width])))
            k += width
        found.append(tuple(segments))
    return tuple(found)


def _path(key: str, colour: tuple[int, int, int]) -> tuple[Segment, ...]:
    """Return a panel's one path of a colour (over ten segments), within its frame."""
    box = frame(key).box
    inside = [
        p
        for p in _paths(colour)
        if all(_inside(*s.points[-1], box, CURVE_MARGIN) for s in p)
    ]
    if len(inside) != 1:
        msg = f"{PANELS[key].name}: {len(inside)} curves of colour {colour}"
        raise ValueError(msg)
    return inside[0]


@dataclass(frozen=True)
class Reading:
    """A vertex as read: its kind, position (pt) and values."""

    kind: str
    x: float
    y: float
    abscissa: float
    value: float


def _off_grid(intensity: float) -> float:
    """Return how far an intensity is from a multiple of STEP, in STEP."""
    return abs(intensity / STEP - round(intensity / STEP))


@functools.cache
def curve(key: str) -> tuple[Reading, ...]:
    """Return a panel's TDSE: the vertices its drawn curve passes through, read.

    The ends of its straight segments and cubics (the first vertex a move), not the
    cubics' control points. Of two vertices within DUPLICATE in x, the one nearer a
    multiple of STEP is kept.
    """
    found = frame(key)
    ends = [
        ("curve end" if s.kind == "cubic" else "line", *s.points[-1])
        for s in _path(key, BLACK)
    ]
    kept: list[tuple[str, float, float]] = []
    for end in ends:
        if kept and abs(end[1] - kept[-1][1]) < DUPLICATE:
            if _off_grid(float(found.x(end[1]))) < _off_grid(
                float(found.x(kept[-1][1]))
            ):
                kept[-1] = end
            continue
        kept.append(end)
    return tuple(
        Reading(kind, x, y, float(found.x(x)), float(found.y(y))) for kind, x, y in kept
    )


def drawn_curve(key: str, colour: tuple[int, int, int]) -> np.ndarray:
    """Return a panel's path of a colour as drawn, (intensity, value) rows.

    The cubics evaluated at SAMPLES points each: what the reader of the figure sees,
    for the checks and the report, not the data.
    """
    found = frame(key)
    out: list[tuple[float, float]] = []
    for s in _path(key, colour):
        if s.kind != "cubic":
            out.append(s.points[-1])
            continue
        p0 = np.array(out[-1])
        c1, c2, p3 = (np.array(p) for p in s.points)
        for t in np.linspace(0, 1, SAMPLES + 1)[1:]:
            point = (
                (1 - t) ** 3 * p0
                + 3 * (1 - t) ** 2 * t * c1
                + 3 * (1 - t) * t**2 * c2
                + t**3 * p3
            )
            out.append((float(point[0]), float(point[1])))
    xy = np.array(sorted(out))
    return np.column_stack(
        [np.asarray(found.x(xy[:, 0])), np.asarray(found.y(xy[:, 1]))]
    )


def at(rows: np.ndarray, intensity: float) -> float:
    """Return a drawn curve's value at an intensity (linear in log of the value)."""
    return float(np.exp(np.interp(intensity, rows[:, 0], np.log(rows[:, 1]))))


def keldysh_intensity(gamma: float, wavelength_nm: float, ip_au: float) -> float:
    """Return the intensity (W/cm^2) at which gamma = omega sqrt(2 I_p)/F."""
    omega = 2 * math.pi * C / (wavelength_nm * 1e-9) / AU_ANGULAR_FREQUENCY
    field = omega * math.sqrt(2 * ip_au) / gamma
    return field**2 * I_AU


def gamma_lines(key: str) -> list[float]:
    """Return the intensities of a panel's dotted vertical lines (W/cm^2).

    Stroked black two-vertex vertical lines inside the frame, nearly its height, drawn
    wider than the frame (0.78 pt against 0.32).
    """
    found = frame(key)
    box = found.box
    edge_width = min(
        d.width for d in found.drawn if d.stroke == BLACK and len(d.vertices) == 5
    )
    xs = sorted(
        float(d.vertices[0, 0])
        for d in found.drawn
        if d.stroke == BLACK
        and len(d.vertices) == 2
        and np.ptp(d.vertices[:, 0]) < STRAIGHT
        and np.ptp(d.vertices[:, 1]) > 0.9 * (box[3] - box[1])
        and box[0] + 1 < d.vertices[0, 0] < box[2] - 1
        and d.width > edge_width * 1.5
    )
    return [float(found.x(x)) for x in xs]


def _check_legend(key: str) -> None:
    """Raise unless the legend's line samples are, from the top, LEGEND's colours.

    The samples: stroked two-vertex horizontal lines 15 to 20 pt long inside the
    frame. The black sample ('TDSE') is drawn as wide as the black curve.
    """
    found = frame(key)
    samples = sorted(
        (
            d
            for d in found.drawn
            if d.stroke is not None
            and len(d.vertices) == 2
            and np.ptp(d.vertices[:, 1]) < STRAIGHT
            and 15 < np.ptp(d.vertices[:, 0]) < 20
            and _inside(*d.centre, found.box, 0)
        ),
        key=lambda d: -d.centre[1],
    )
    if tuple(d.stroke for d in samples) != LEGEND:
        msg = f"{PANELS[key].name}: the legend's samples {[d.stroke for d in samples]}"
        raise ValueError(msg)
    curves = [
        d
        for d in found.drawn
        if d.stroke == BLACK and len(d.vertices) > 10 and d.fill is None
    ]
    if len(curves) != 1 or not math.isclose(curves[0].width, samples[0].width):
        msg = f"{PANELS[key].name}: the TDSE curve is not drawn as its legend's sample"
        raise ValueError(msg)


def checks() -> None:
    """Raise unless the readings agree with every number the paper prints."""
    for key, panel in PANELS.items():
        _check_legend(key)
        read = curve(key)
        for r in read:
            if _off_grid(r.abscissa) > GRID:
                msg = f"{panel.name}: a vertex at {r.abscissa:.4e}, off the 1e12 grid"
                raise ValueError(msg)
        low, high = SPAN[key]
        if (
            abs(read[0].abscissa - low) > GRID * STEP
            or abs(read[-1].abscissa - high) > GRID * STEP
        ):
            msg = f"{panel.name}: the TDSE spans {read[0].abscissa}-{read[-1].abscissa}"
            raise ValueError(msg)
        lines = gamma_lines(key)
        expected = keldysh_intensity(1.0, panel.wavelength, H_IP_AU)
        beyond = expected > float(frame(key).x(frame(key).box[2]))
        if (beyond and lines) or (
            not beyond and (len(lines) != 1 or abs(lines[0] / expected - 1) > GAMMA)
        ):
            msg = f"{panel.name}: gamma = 1 lines at {lines}, expected {expected:.4e}"
            raise ValueError(msg)


def adk_ratio() -> float:
    """Return Fig. 12(b)'s drawn ADK over its drawn TDSE at ADK_AT."""
    return at(drawn_curve("b", RED), ADK_AT) / at(drawn_curve("b", BLACK), ADK_AT)


def _points(key: str) -> list[Point]:
    """Return a dataset's points, numbered in order of intensity."""
    return [
        Point(
            f"p{k:03d}",
            f"{r.abscissa:.4e}",
            f"{r.value:.4g}",
            "",
            f"{PANELS[key].name} TDSE: {r.kind} vertex at ({r.x:.3f}, {r.y:.3f}) pt",
        )
        for k, r in enumerate(curve(key), start=1)
    ]


def extract() -> dict[str, list[Point]]:
    """Return the TDSE probabilities of Fig. 12, one file per panel."""
    checks()
    return {relative: _points(key) for key, relative in DATASETS.items()}


def report() -> None:
    """Print the checks' numbers and the PPT curves over the TDSE."""
    for key, panel in PANELS.items():
        found = frame(key)
        read = curve(key)
        kinds = [r.kind for r in read]
        print(
            f"{panel.name}: {len(read)} vertices ({kinds.count('curve end')} cubic "
            f"ends), {read[0].abscissa:.4e}-{read[-1].abscissa:.4e} W/cm^2, "
            f"{read[0].value:.4g}-{read[-1].value:.4g}; largest offset from the grid "
            f"{max(_off_grid(r.abscissa) for r in read):.4f}e12"
        )
        x_ticks = np.array(found.x.positions)
        print(
            f"  x {1e14 / found.x.resolution():.2f} pt per 1e14 W/cm^2, ticks "
            f"{np.round(np.diff(x_ticks), 3).tolist()}"
        )
        print(f"  y {math.log(10) / found.y.resolution():.3f} pt per decade")
        expected = keldysh_intensity(1.0, panel.wavelength, H_IP_AU)
        print(f"  gamma = 1 lines {gamma_lines(key)}, expected {expected:.4e}")
        ppt = drawn_curve(key, GREEN)
        ratios = [
            (r.abscissa, at(ppt, r.abscissa) / panel.ppt / r.value)
            for r in read
            if ppt[0, 0] <= r.abscissa <= ppt[-1, 0]
        ]
        print(
            f"  PPT (legend's x{panel.ppt} taken out) over the TDSE at its vertices:",
            [(f"{i:.2e}", round(q, 3)) for i, q in ratios],
        )
    print(f"Fig. 12(b): ADK/TDSE at {ADK_AT:.0e} W/cm^2 {adk_ratio():.4f}")


def main() -> None:
    """Write the data files and print the checks' numbers."""
    write(extract())
    report()


if __name__ == "__main__":
    main()
