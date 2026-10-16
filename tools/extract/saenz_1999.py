# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Saenz and Lambropoulos (1999): two- to four-photon ionisation of He.

The paper prints no table: its generalised N-photon ionisation cross sections, partial
by the final state's symmetry, are given only in figures, all vector drawings on the
page (the PDF's one image is the download banner of page 1), read from the drawing
(data rule 8, TC-699). The PDF's page 1 is that banner: journal p. 5629 is PDF page 2,
so the journal page is the PDF page plus 5627.

* **Fig. 1** (PDF p. 6, journal p. 5633, library figure ``fig_c3ccd1aa8c7b``): two
  photons, (a) the 1S and (b) the 1D final states, 12.65 to 23.85 eV.
* **Fig. 2** (p. 7/5634, ``fig_4538c5a2639e``): three photons, (a) 1P and (b) 1F,
  8.5 to 11.95 eV.
* **Fig. 3** (p. 8/5635, ``fig_f558c85055c3``, and its continuation p. 9/5636,
  ``fig_38efdb3d55c0``): four photons, (a) 1S, (b) 1D and (c) 1G, 6.35 to 8.15 eV.
* **Fig. 4** (p. 9/5636, ``fig_f06112f44e83``): two photons, 1D, the length (dashed)
  against the velocity (solid) gauge, 12.65 to 23.85 eV.

The photon-energy axes are as drawn: the authors moved each order's energies so that
the calculated ionisation threshold falls at the true one (their ground state is about
0.69 eV off; p. 5632). Nothing here shifts them back (the record's notes).

How the drawing is read:

* every path is a stroked straight segment of two vertices (no Bezier, no form), all of
  one width per plot (0.81 to 0.83 pt): frames, ticks, the curves and the axis labels'
  glyphs alike. The frames: four segments forming the rectangle :data:`FIGURES` gives
  (within :data:`FRAME_FIND`).
* the ticks: segments perpendicular to an edge with one end on it (within
  :data:`EDGE`), pointing inwards, :data:`SHORT` (about 2.09 pt) or :data:`LONG` (4.17
  pt) long; the top and right edges repeat the bottom and left ones (checked, within
  :data:`MIRROR`). The labels are glyphs drawn as strokes, not text: their values are
  written in :data:`FIGURES` as read on the rendered page. x linear: the long ticks are
  the labels, the short ones every :attr:`Figure.x_minor` eV (also outside the first
  and last labels); y logarithmic: the long ticks are the decades, the short ones 2,
  ..., 9 times a decade (also below the first and above the last). Each short tick is
  assigned its value from the long ticks' scale (it must lie within a fifth of a step,
  or of the log step, of it), then every tick lies on the fitted scale within
  :data:`TICK_TOLERANCE` (:meth:`_figure.Axis.check`).
* the curves: the other segments inside the frame. The program drew each curve as a
  chain of segments, each starting where the last ended, in drawing order; a chain
  breaks where the curve leaves the frame (and is clipped there) and at the dashes of
  Fig. 4's dashed line. In Figs. 1-3 every chain belongs to the one curve, and the
  chains are taken in order of their first photon energy; in Fig. 4 the solid line is
  the one chain that runs from the left edge of the frame to the right, the dashed line
  every other chain (its dashes, in drawing order). Each curve must never step back in
  photon energy.
* the values: every vertex, less exact repeats, and less the vertices on the frame's
  bottom or top edge, where the frame clips the line (the value there is only bounded).
  Rows are numbered in drawing order; each ``printed`` cell gives the vertex's position
  in page points (and the dash for Fig. 4's dashed line).

**Checks** (:func:`checks`), against every number the paper prints for these curves
(pp. 5632, 5634-5635) and between Figs. 1(b) and 4, which draw the same velocity-gauge
calculation (established here); the script fails if one is missed by more than its
stated margin.
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

DIRECTORY = "Saenz_1999_TheoreticalTwoThree"
SOURCE = "saenz_1999"
FRAME_FIND = 0.5
"""The frame's edges lie within this of :data:`FIGURES`' position, pt."""
EDGE = 0.02
"""A tick's end lies on its frame edge within this, pt."""
SHORT = (1.9, 2.3)
"""The length of a short tick, pt (2.086 drawn)."""
LONG = (3.9, 4.4)
"""The length of a long tick, pt (4.17 drawn)."""
MIRROR = 0.02
"""A top or right tick repeats a bottom or left one within this, pt."""
TICK_TOLERANCE = 0.05
"""Every tick lies on its fitted scale within this, pt."""
ASSIGN = 0.2
"""A short tick lies within this fraction of a step of the value it is given."""
JOIN = 0.005
"""A chain's next segment starts within this of the last one's end, pt (the drawing
rounds its coordinates to about 0.001 pt), and a curve may step back by as much."""


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
    """The decades the y labels print (at the long ticks), in the axis unit."""
    unit: float
    """The y axis's unit: '10^-52 cm^4 sec^1' and so on."""


FIGURES = {
    "fig1a": Figure(
        "Fig. 1(a)",
        6,
        (186.025, 558.038, 461.415, 744.007),
        (14, 16, 18, 20, 22),
        1.0,
        (1e-4, 1e-3, 0.01, 0.1, 1, 10, 100, 1000),
        1e-52,
    ),
    "fig1b": Figure(
        "Fig. 1(b)",
        6,
        (186.025, 329.661, 461.415, 509.191),
        (14, 16, 18, 20, 22),
        1.0,
        (1e-3, 0.01, 0.1, 1, 10, 100, 1000, 1e4),
        1e-52,
    ),
    "fig2a": Figure(
        "Fig. 2(a)",
        7,
        (186.025, 557.829, 461.415, 743.797),
        (9, 10, 11),
        0.2,
        (1e-3, 0.01, 0.1, 1, 10, 100, 1000, 1e4, 1e5),
        1e-84,
    ),
    "fig2b": Figure(
        "Fig. 2(b)",
        7,
        (186.025, 331.22, 461.415, 517.188),
        (9, 10, 11),
        0.2,
        (0.01, 0.1, 1, 10),
        1e-84,
    ),
    "fig3a": Figure(
        "Fig. 3(a)",
        8,
        (186.236, 558.128, 461.627, 744.097),
        (6.5, 7.0, 7.5, 8.0),
        0.1,
        (1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 1, 10, 100),
        1e-116,
    ),
    "fig3b": Figure(
        "Fig. 3(b)",
        8,
        (186.236, 326.517, 461.626, 512.486),
        (6.5, 7.0, 7.5, 8.0),
        0.1,
        (1e-5, 1e-4, 1e-3, 0.01, 0.1, 1, 10, 100, 1000),
        1e-116,
    ),
    "fig3c": Figure(
        "Fig. 3(c)",
        9,
        (186.025, 557.826, 461.415, 743.794),
        (6.5, 7.0, 7.5, 8.0),
        0.1,
        (1e-3, 0.01, 0.1),
        1e-116,
    ),
    "fig4": Figure(
        "Fig. 4",
        9,
        (186.034, 301.845, 461.424, 483.596),
        (14, 16, 18, 20, 22),
        1.0,
        (1e-3, 0.01, 0.1, 1, 10, 100, 1000, 1e4),
        1e-52,
    ),
}

# (figure, curve) -> file. Fig. 4's curves are "solid" (velocity) and "dashed"
# (length); the other figures have one curve, "line".
FILES = {
    ("fig1a", "line"): "He/saenz_1999__fig1a_sigma2_1S_He.csv",
    ("fig1b", "line"): "He/saenz_1999__fig1b_sigma2_1D_He.csv",
    ("fig2a", "line"): "He/saenz_1999__fig2a_sigma3_1P_He.csv",
    ("fig2b", "line"): "He/saenz_1999__fig2b_sigma3_1F_He.csv",
    ("fig3a", "line"): "He/saenz_1999__fig3a_sigma4_1S_He.csv",
    ("fig3b", "line"): "He/saenz_1999__fig3b_sigma4_1D_He.csv",
    ("fig3c", "line"): "He/saenz_1999__fig3c_sigma4_1G_He.csv",
    ("fig4", "solid"): "He/saenz_1999__fig4_sigma2_1D_velocity_He.csv",
    ("fig4", "dashed"): "He/saenz_1999__fig4_sigma2_1D_length_He.csv",
}


@functools.cache
def _segments(page: int) -> tuple[tuple[np.ndarray, ...], tuple[float, ...]]:
    """Return a page's stroked segments (2 x 2 arrays, pt) in drawing order, and widths.

    Raises unless every stroked path is one straight segment outside any form.
    """
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    found, widths = [], []
    for obj in document[page - 1].get_objects(max_depth=15):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        fill, stroked = ctypes.c_int(), ctypes.c_int()
        raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill), ctypes.byref(stroked))
        if not stroked.value:
            continue  # the filled glyphs of the panel letters, (a), (b), (c)
        count = raw.FPDFPath_CountSegments(obj.raw)
        kinds = [
            raw.FPDFPathSegment_GetType(raw.FPDFPath_GetPathSegment(obj.raw, k))
            for k in range(count)
        ]
        if (
            obj.container is not None
            or count != 2
            or kinds != [raw.FPDF_SEGMENT_MOVETO, raw.FPDF_SEGMENT_LINETO]
        ):
            msg = f"p. {page}: a stroked path that is not one straight segment"
            raise ValueError(msg)
        matrix = obj.get_matrix()
        points = []
        for k in range(count):
            segment = raw.FPDFPath_GetPathSegment(obj.raw, k)
            x, y = ctypes.c_float(), ctypes.c_float()
            raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
            points.append(matrix.on_point(x.value, y.value))
        found.append(np.array(points, dtype=float))
        width = ctypes.c_float()
        raw.FPDFPageObj_GetStrokeWidth(obj.raw, ctypes.byref(width))
        widths.append(float(width.value))
    return tuple(found), tuple(widths)


def _is_frame(seg: np.ndarray, frame: tuple[float, float, float, float]) -> bool:
    """Whether a segment is one of the frame's edges: from one corner to another."""
    left, bottom, right, top = frame
    corners = np.array([[left, bottom], [right, bottom], [right, top], [left, top]])
    ends = [int(np.argmin(np.hypot(*(corners - end).T))) for end in seg]
    return ends[0] != ends[1] and all(
        float(np.hypot(*(corners[k] - end))) < FRAME_FIND
        for k, end in zip(ends, seg, strict=True)
    )


def _frame(key: str) -> tuple[float, float, float, float]:
    """Return a plot's frame (l, b, r, t, pt), checked against its four drawn edges."""
    figure = FIGURES[key]
    edges = sum(_is_frame(seg, figure.frame) for seg in _segments(figure.page)[0])
    if edges != 4:
        msg = f"{figure.name}: {edges} frame edges near {figure.frame}"
        raise ValueError(msg)
    return figure.frame


def _tick(
    seg: np.ndarray, frame: tuple[float, float, float, float]
) -> tuple[str, float, bool] | None:
    """Return (edge, position, long) if a segment is a tick of the frame, else None."""
    left, bottom, right, top = frame
    (x0, y0), (x1, y1) = seg
    length = math.hypot(x1 - x0, y1 - y0)
    if SHORT[0] <= length <= SHORT[1]:
        is_long = False
    elif LONG[0] <= length <= LONG[1]:
        is_long = True
    else:
        return None
    if abs(x0 - x1) < EDGE and left - EDGE <= x0 <= right + EDGE:
        for name, edge, inward in (("bottom", bottom, 1), ("top", top, -1)):
            for end, other in ((y0, y1), (y1, y0)):
                if abs(end - edge) < EDGE and (other - end) * inward > 0:
                    return name, round(float(x0), 3), is_long
    if abs(y0 - y1) < EDGE and bottom - EDGE <= y0 <= top + EDGE:
        for name, edge, inward in (("left", left, 1), ("right", right, -1)):
            for end, other in ((x0, x1), (x1, x0)):
                if abs(end - edge) < EDGE and (other - end) * inward > 0:
                    return name, round(float(y0), 3), is_long
    return None


@dataclass(frozen=True)
class Ticks:
    """A plot's tick positions along one edge, pt: every one and the long ones."""

    every: tuple[float, ...]
    long: tuple[float, ...]


def _ticks(key: str) -> tuple[Ticks, Ticks]:
    """Return the x ticks (bottom edge) and y ticks (left edge) of a plot.

    Checks that the top and right edges repeat them within :data:`MIRROR`.
    """
    figure = FIGURES[key]
    frame = _frame(key)
    edges: dict[str, list[tuple[float, bool]]] = {
        "bottom": [],
        "top": [],
        "left": [],
        "right": [],
    }
    for seg in _segments(figure.page)[0]:
        found = _tick(seg, frame)
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
        values.append(k * 10.0**decade)
    if len(set(values)) != len(values):
        msg = f"{name}: two y ticks at one value"
        raise ValueError(msg)
    return Axis(ticks.every, tuple(values), log=True, tolerance=TICK_TOLERANCE)


@functools.cache
def axes(key: str) -> tuple[Axis, Axis]:
    """Return a plot's photon-energy (eV) and value (axis unit) axes, checked."""
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


def _chains(key: str) -> list[np.ndarray]:
    """Return a plot's curve chains, in drawing order: polylines of joined segments."""
    figure = FIGURES[key]
    frame = _frame(key)
    left, bottom, right, top = frame
    chains: list[list[np.ndarray]] = []
    last = -2
    segments, _ = _segments(figure.page)
    for index, seg in enumerate(segments):
        inside = (
            np.all(seg[:, 0] >= left - EDGE)
            and np.all(seg[:, 0] <= right + EDGE)
            and np.all(seg[:, 1] >= bottom - EDGE)
            and np.all(seg[:, 1] <= top + EDGE)
        )
        if not inside or _is_frame(seg, frame) or _tick(seg, frame) is not None:
            continue
        if (
            chains
            and index == last + 1
            and np.allclose(chains[-1][-1][1], seg[0], 0, JOIN)
        ):
            chains[-1].append(seg)
        else:
            chains.append([seg])
        last = index
    return [np.vstack([c[0][:1]] + [s[1:] for s in c]) for c in chains]


def _monotonic(line: np.ndarray, name: str) -> np.ndarray:
    """Drop exact repeats; raise if the polyline steps back in photon energy."""
    keep = np.ones(len(line), dtype=bool)
    keep[1:] = np.any(np.diff(line, axis=0) != 0, axis=1)
    line = line[keep]
    if np.any(np.diff(line[:, 0]) < -JOIN):
        msg = f"{name}: the curve steps back in photon energy"
        raise ValueError(msg)
    return line


@functools.cache
def _vertices(key: str, curve: str) -> tuple[tuple[np.ndarray, int], ...]:
    """Return a curve's vertices (pt) with the number of the chain each is in."""
    figure = FIGURES[key]
    chains = _chains(key)
    left, _, right, _ = _frame(key)
    if key == "fig4":
        spans = [
            k
            for k, c in enumerate(chains)
            if c[0, 0] - left < EDGE and right - c[-1, 0] < EDGE
        ]
        if len(spans) != 1:
            msg = f"{figure.name}: {len(spans)} chains span the frame"
            raise ValueError(msg)
        chosen = (
            spans
            if curve == "solid"
            else [k for k in range(len(chains)) if k != spans[0]]
        )
    else:
        chosen = sorted(range(len(chains)), key=lambda k: chains[k][0, 0])
    out = [(vertex, k) for k in chosen for vertex in chains[k]]
    line = _monotonic(np.array([v for v, _ in out]), f"{figure.name} {curve}")
    if len(line) != len(out):  # repeats dropped: rebuild with the chain numbers
        kept, previous = [], None
        for vertex, k in out:
            if previous is None or np.any(vertex != previous):
                kept.append((vertex, k))
            previous = vertex
        out = kept
    return tuple(out)


@functools.cache
def curve(key: str, name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return a curve: page positions (pt), photon energies (eV), values, chain numbers.

    The values are in cm^2N s^(N-1): the axis's printed unit times
    :attr:`Figure.unit`. Vertices on the frame's bottom or top edge are dropped.
    """
    figure = FIGURES[key]
    _, bottom, _, top = _frame(key)
    found = [
        (v, k)
        for v, k in _vertices(key, name)
        if abs(v[1] - bottom) > EDGE and abs(v[1] - top) > EDGE
    ]
    line = np.array([v for v, _ in found])
    chain = np.array([k for _, k in found])
    x_axis, y_axis = axes(key)
    energy = np.asarray(x_axis(line[:, 0]), dtype=float)
    value = np.asarray(y_axis(line[:, 1]), dtype=float) * figure.unit
    return line, energy, value, chain


def _points(key: str, name: str) -> list[Point]:
    """Return a curve as data points, one per vertex."""
    line, energy, value, chain = curve(key, name)
    dashed = name == "dashed"
    return [
        Point(
            f"v{k:04d}",
            f"{e:.6f}",
            f"{v:.4e}",
            "",
            (f"dash {c + 1}, " if dashed else "")
            + f"vertex {k} at ({x:.3f}, {y:.3f}) pt",
        )
        for k, ((x, y), e, v, c) in enumerate(
            zip(line, energy, value, chain, strict=True), 1
        )
    ]


def at(key: str, name: str, energy: float) -> float:
    """Return a curve's value at a photon energy (eV), log-interpolated.

    Between the two vertices that bracket it, first in drawing order; for the dashed
    line, only within a dash (raises in a gap).
    """
    _, energies, values, chain = curve(key, name)
    above = int(np.argmax(energies >= energy))
    if energies[above] < energy or above == 0 or chain[above] != chain[above - 1]:
        msg = f"{FIGURES[key].name} {name}: {energy} eV is not on a drawn stretch"
        raise ValueError(msg)
    e0, e1 = energies[above - 1], energies[above]
    v0, v1 = values[above - 1], values[above]
    t = 0.0 if e1 == e0 else (energy - e0) / (e1 - e0)
    return float(10 ** ((1 - t) * math.log10(v0) + t * math.log10(v1)))


def _dash_near(key: str, energy: float) -> float:
    """Return the dashed line's value at the drawn vertex nearest a photon energy."""
    _, energies, values, _ = curve(key, "dashed")
    return float(values[int(np.argmin(np.abs(energies - energy)))])


# The margins. DRAWN: the digitising uncertainty of a value (0.5 %, the record's notes)
# with log-interpolation between vertices, and the dashed line read at its nearest
# vertex (within 0.03 eV of 13.5 eV, where the curve changes 1 % per 0.1 eV). The
# printed numbers' own rounding is added to it: 1.3 and 1.35 to half a unit of their
# last digit (4 % and 0.4 %), 0.2 and 0.15 to one significant figure (25 % and 3.3 %).
DRAWN = 0.03
SAME = 0.03
"""Figs. 1(b) and 4 (solid) drawn from one calculation agree within this."""


def checks() -> dict[str, object]:
    """Return the checks against the text and between Figs. 1(b) and 4.

    Raises if one fails its margin.
    """
    out: dict[str, object] = {}
    failures = []

    def compare(label: str, read: float, printed: float, margin: float) -> None:
        out[label] = (read, printed, read / printed - 1)
        if abs(read / printed - 1) > margin:
            failures.append(f"{label}: read {read:.4g}, printed {printed:.4g}")

    # p. 5634: at 13.5 eV, 1D: 'sigma_D = 1.3 x 10^-51 cm4 s (velocity gauge) and 1.35 x
    # 10^-51 cm4 s (length gauge)'; 1S: 'sigma_S = 0.2 x 10^-51 cm4 s (velocity) and
    # 0.15 x 10^-51 cm4 s (length)'.
    compare(
        "1D velocity 13.5 eV, Fig. 4 solid",
        at("fig4", "solid", 13.5),
        1.3e-51,
        DRAWN + 0.04,
    )
    compare(
        "1D length 13.5 eV, Fig. 4 dashed",
        _dash_near("fig4", 13.5),
        1.35e-51,
        DRAWN + 0.004 + 0.02,
    )
    compare(
        "1D 13.5 eV, Fig. 1(b) against velocity",
        at("fig1b", "line", 13.5),
        1.3e-51,
        DRAWN + 0.04,
    )
    compare(
        "1S 13.5 eV, Fig. 1(a) against velocity",
        at("fig1a", "line", 13.5),
        0.2e-51,
        DRAWN + 0.25,
    )
    # Reported, not a pass/fail check: Fig. 1(a) against the length-gauge 0.15.
    s = at("fig1a", "line", 13.5)
    out["1S 13.5 eV, Fig. 1(a) against length 0.15e-51"] = (
        s,
        0.15e-51,
        s / 0.15e-51 - 1,
    )
    # Fig. 1(b) and Fig. 4's solid line: the same velocity-gauge calculation, off the
    # resonances.
    for energy in (13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0):
        compare(
            f"Fig. 1(b) over Fig. 4 solid at {energy} eV",
            at("fig1b", "line", energy),
            at("fig4", "solid", energy),
            SAME,
        )
    # p. 5632: 'From figure 1 ... yielding the value sigma(2) ~ 10^-51 cm4 s' at the
    # ninth harmonic (13.5 eV, p. 5630): the sum of Figs. 1(a) and 1(b), reported.
    total = at("fig1a", "line", 13.5) + at("fig1b", "line", 13.5)
    out["sigma2 total 13.5 eV, Fig. 1 (a) + (b), against 1e-51"] = (total, 1e-51)
    # p. 5635: 'According to figure 3 the generalized cross section is sigma(4) ~
    # 10^-118 cm8 s3, or perhaps larger by a factor of about two' at 7.5 eV: the sum of
    # the three partial cross sections must lie between 1e-118 and four times it (the
    # text's factor of two, with the reading of an order-of-magnitude statement).
    parts = [at(k, "line", 7.5) for k in ("fig3a", "fig3b", "fig3c")]
    total4 = sum(parts)
    out["sigma4 7.5 eV, 1S, 1D, 1G"] = tuple(parts)
    out["sigma4 7.5 eV, Fig. 3 total, against 1e-118"] = (total4, 1e-118)
    if not 1e-118 <= total4 <= 4e-118:
        failures.append(f"sigma4 at 7.5 eV: {total4:.3g}, not 1-4 x 1e-118")
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return every curve, after the checks."""
    checks()
    return {path: _points(key, name) for (key, name), path in FILES.items()}


def main() -> None:
    """Write the data files, and print the readings and the checks."""
    write(extract())
    for (key, name), _ in FILES.items():
        line, energies, values, chain = curve(key, name)
        x_res, y_res = residuals(key)
        print(
            f"{FIGURES[key].name} {name}: {len(line)} vertices in "
            f"{len(set(chain.tolist()))} chains, {energies[0]:.4f} to "
            f"{energies[-1]:.4f} eV, {values.min():.3e} to {values.max():.3e}; "
            f"ticks within {x_res:.3f} (x) and {y_res:.3f} (y) pt"
        )
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
