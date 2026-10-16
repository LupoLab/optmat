# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Sato et al. (2018): single and double ionisation of Ar by TD-OCC and others.

The paper prints no number of its ionisation probabilities: they are given only in
**Fig. 3** (PDF p. 5, journal p. 051101-4, library figure ``fig_5d38508bb252``), the
probabilities of finding one (a) and two (b) electrons outside a sphere of radius
R0 = 20 a.u. against time in optical cycles, 1.0 to 3.0, for TDHF, TD-OCCD, TD-OCCDT
and TD-CASSCF, Ar at 800 nm, 6 x 10^14 W/cm^2, three cycles foot to foot. The panels
are labelled '(a) Ar+' and '(b) Ar2+'. The figure is a vector drawing on the page (no
image; the PDF's one image is the download banner of page 1), read from the drawing
(data rule 8, TC-699) with :func:`_figure.painted`. Its labels are glyphs drawn as
filled paths, not text: their values are written in :data:`PANELS` and
:data:`LEGEND` as read on the rendered page, and the script checks that a label's
glyphs sit at each long tick and beside each legend sample.

How the drawing is read:

* the frames: a black stroked closed rectangle (five vertices, drawn twice) per panel,
  near the position :data:`PANELS` gives; the ticks black stroked two-vertex lines of
  3.63 pt, perpendicular to an edge with one end on it (within :data:`EDGE`), all
  labelled (on the left, and at the bottom of panel (b)). The top and right edges
  repeat the bottom and left ones' ticks (checked, within :data:`MIRROR`); panel (a)'s
  x ticks are panel (b)'s (checked).
* the axes: linear, the y labels 0.0 to 0.5 in (a) and 0.00 to 0.04 in (b), the x
  labels 1.0 to 3.0 cycles. Each tick has a filled glyph path centred on it within
  :data:`LABEL`, outside the frame (checked; the count of ticks is the count of
  labels). Every tick lies on its fitted scale within :data:`TICK_TOLERANCE`
  (:meth:`_figure.Axis.check`): the drawing's coordinates step by about 0.058 pt.
* the curves (:data:`LEGEND`): TDHF red (1.152 pt), TD-OCCD blue (1.152 pt), TD-OCCDT
  green (0.576 pt) and TD-CASSCF black (1.727 pt), stroked paths of more than two
  vertices inside the frame, two or three pieces each, joined end to end (each next
  piece starts within :data:`JOIN` of the last one's end). The legend's samples (two
  vertices) lie top to bottom in the order of :data:`LEGEND`, as its labels on the
  rendered page, each with a filled glyph path to its left at its height (checked). A
  polyline must not step back in time by more than :data:`BACKSTEP`, and must end at
  its rightmost vertex.
* the value: each curve's last vertex, on the 3.0-cycle tick (checked, within
  :data:`EDGE`): the end of the pulse.

The digitising uncertainty is :data:`DIGITISING` pt on the y axis (one quantum, which
bounds half a quantum plus the worst tick residual, checked): 2.9e-4 in (a) and 2.4e-5
in (b), absolute.

**Checks** (:func:`checks`). The text prints no number for these curves. The script
checks the calibration's zero: every curve starts at t = 1.0 cycle (before the field
ionises; the pulse starts at t = 0) at 0 within :data:`START` pt. It reports the range
of each double-ionisation curve over the last :data:`TAIL` cycle, as the curves
oscillate there.
"""

import functools
import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from _figure import Axis, Painted, painted
from _library import Point, write

DIRECTORY = "Sato_2018_CommunicationTimeDependent"
SOURCE = "sato_2018"
PAGE = 5
BLACK = (0, 0, 0)
RED = (255, 0, 0)
BLUE = (0, 0, 255)
GREEN = (0, 255, 0)
TICK_MAX = 5.0
"""Ticks are shorter than this, pt (3.63)."""
EDGE = 0.02
"""A tick's end lies on its frame edge within this, pt."""
MIRROR = 0.02
"""A top or right tick (or the other panel's x tick) repeats one within this, pt."""
FRAME_FIND = 2.0
"""A frame's corners lie within this of :data:`PANELS`' position, pt."""
LABEL = 1.0
"""A label's glyphs are centred within this of their tick or sample, across, pt."""
LABEL_REACH = 25.0
"""A label's glyphs lie within this of the frame or the sample, outside it, pt."""
TICK_TOLERANCE = 0.06
"""Every tick lies on its fitted scale within this, pt: the drawing's quantum."""
QUANTUM = 0.058
"""The step of the drawing's coordinates, pt (curve vertices 739.905 and 739.848)."""
DIGITISING = 0.06
"""The digitising uncertainty of a vertex on the y axis, pt: one quantum, which bounds
half the quantum plus the worst tick residual (checked)."""
JOIN = 0.005
"""A curve's next piece starts within this of the last one's end, pt (0.001 seen)."""
START = 0.12
"""Every curve starts at t = 1 cycle at 0 within this, pt: two quanta."""
TAIL = 0.2
"""The double-ionisation curves' range is reported over this last part, cycles."""
BACKSTEP = 1.0
"""A curve may step back in time by at most this, pt: the drawn lines step back at a
few kinks, by 0.12 to 0.88 pt (at most 0.022 cycle, TD-OCCD in (b) at 2.45 cycles;
reported by :func:`checks`), none at the end."""


@dataclass(frozen=True)
class Panel:
    """One plot: where it is and what its tick labels print."""

    name: str
    frame: tuple[float, float, float, float]
    """The frame's left, bottom, right, top, pt, to within :data:`FRAME_FIND`."""
    y_labels: tuple[float, ...]
    """The values the y labels print, bottom to top."""
    ion: str
    """The panel's label: the ion whose probability it draws."""


PANELS = {
    "fig3a": Panel(
        "Fig. 3(a)",
        (363.01, 632.73, 524.55, 742.38),
        (0.0, 0.1, 0.2, 0.3, 0.4, 0.5),
        "Ar+",
    ),
    "fig3b": Panel(
        "Fig. 3(b)",
        (363.01, 517.51, 524.55, 627.16),
        (0.00, 0.01, 0.02, 0.03, 0.04),
        "Ar2+",
    ),
}
X_LABELS = (1.0, 1.5, 2.0, 2.5, 3.0)
"""The x labels below panel (b), optical cycles."""
X_FROM = "fig3b"
"""The panel whose x axis is labelled."""


@dataclass(frozen=True)
class Method:
    """One curve: its legend label, stroke colour and width (pt), and file key."""

    label: str
    colour: tuple[int, int, int]
    width: float
    key: str


LEGEND = (
    Method("TDHF", RED, 1.152, "tdhf"),
    Method("TD-OCCD", BLUE, 1.152, "occd"),
    Method("TD-OCCDT", GREEN, 0.576, "occdt"),
    Method("TD-CASSCF", BLACK, 1.727, "casscf"),
)
"""The legend, top to bottom, as on the rendered page (both panels)."""
INTENSITY = "6e14"
"""The peak intensity, W/cm^2, as the text prints it."""


def file(panel: str, method: Method) -> str:
    """Return a curve's data file, as the record names it."""
    return f"Ar/{SOURCE}__{panel}_{method.key}_Ar_800nm.csv"


@functools.cache
def _paths() -> tuple[Painted, ...]:
    """Return the page's painted paths (the whole page)."""
    return tuple(painted(DIRECTORY, PAGE, (0.0, 0.0, 1000.0, 1000.0)))


@functools.cache
def frame(key: str) -> tuple[float, float, float, float]:
    """Return a panel's frame (l, b, r, t, pt): its black stroked closed rectangle."""
    panel = PANELS[key]
    found: list[tuple[float, float, float, float]] = []
    for path in _paths():
        v = path.vertices
        if path.stroke != BLACK or len(v) != 5 or not np.allclose(v[0], v[-1]):
            continue
        box = (
            float(v[:, 0].min()),
            float(v[:, 1].min()),
            float(v[:, 0].max()),
            float(v[:, 1].max()),
        )
        if all(abs(a - b) < FRAME_FIND for a, b in zip(box, panel.frame, strict=True)):
            found.append(box)
    # The frame is drawn twice, the copies 0.001 pt apart: one frame within EDGE.
    if not found or np.ptp(np.array(found), axis=0).max() > EDGE:
        msg = f"{panel.name}: {len(found)} frames near {panel.frame}"
        raise ValueError(msg)
    left, bottom, right, top = found[0]
    return left, bottom, right, top


@functools.cache
def ticks(key: str) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Return the x ticks (bottom edge) and y ticks (left edge) of a panel, pt.

    Checks that the top and right edges repeat them within :data:`MIRROR`.
    """
    panel = PANELS[key]
    left, bottom, right, top = frame(key)
    edges: dict[str, list[float]] = {"bottom": [], "top": [], "left": [], "right": []}
    for path in _paths():
        if path.stroke != BLACK or len(path.vertices) != 2:
            continue
        (x0, y0), (x1, y1) = path.vertices
        if math.hypot(x1 - x0, y1 - y0) > TICK_MAX:
            continue
        if abs(x0 - x1) < EDGE and left - EDGE <= x0 <= right + EDGE:
            for name, edge in (("bottom", bottom), ("top", top)):
                if min(abs(y0 - edge), abs(y1 - edge)) < EDGE:
                    edges[name].append(float(x0))
        if abs(y0 - y1) < EDGE and bottom - EDGE <= y0 <= top + EDGE:
            for name, edge in (("left", left), ("right", right)):
                if min(abs(x0 - edge), abs(x1 - edge)) < EDGE:
                    edges[name].append(float(y0))
    out = []
    for near, far in (("bottom", "top"), ("left", "right")):
        every = sorted({round(p, 3) for p in edges[near]})
        mirror = sorted({round(p, 3) for p in edges[far]})
        if (
            not every
            or len(mirror) != len(every)
            or any(abs(m - p) > MIRROR for m, p in zip(mirror, every, strict=True))
        ):
            msg = f"{panel.name}: the {far} ticks do not repeat the {near} ones"
            raise ValueError(msg)
        out.append(tuple(every))
    return out[0], out[1]


def _glyphs() -> list[Painted]:
    """Return the page's filled black paths: the glyphs of the figure's labels."""
    return [p for p in _paths() if p.fill == BLACK and p.stroke is None]


def _check_labels(key: str) -> None:
    """Check that a label's glyphs sit at every tick of a panel, outside its frame."""
    panel = PANELS[key]
    left, bottom, _, _ = frame(key)
    x_ticks, y_ticks = ticks(key)
    glyphs = _glyphs()
    for y in y_ticks:
        if not any(
            abs(g.centre[1] - y) < LABEL and 0 < left - g.right < LABEL_REACH
            for g in glyphs
        ):
            msg = f"{panel.name}: no label at the y tick {y}"
            raise ValueError(msg)
    if key == X_FROM:
        for x in x_ticks:
            if not any(
                abs(g.centre[0] - x) < LABEL and 0 < bottom - g.top < LABEL_REACH
                for g in glyphs
            ):
                msg = f"{panel.name}: no label at the x tick {x}"
                raise ValueError(msg)


@functools.cache
def axes(key: str) -> tuple[Axis, Axis]:
    """Return a panel's time (optical cycles) and probability axes, checked."""
    panel = PANELS[key]
    x_ticks, y_ticks = ticks(key)
    shared, _ = ticks(X_FROM)
    if len(shared) != len(x_ticks) or any(
        abs(a - b) > MIRROR for a, b in zip(shared, x_ticks, strict=True)
    ):
        msg = f"{panel.name}: its x ticks are not {PANELS[X_FROM].name}'s"
        raise ValueError(msg)
    if len(shared) != len(X_LABELS) or len(y_ticks) != len(panel.y_labels):
        msg = f"{panel.name}: the ticks do not match the labels"
        raise ValueError(msg)
    _check_labels(key)
    x_axis = Axis(shared, X_LABELS, tolerance=TICK_TOLERANCE)
    y_axis = Axis(y_ticks, panel.y_labels, tolerance=TICK_TOLERANCE)
    for axis in (x_axis, y_axis):
        axis.check()
    return x_axis, y_axis


def residuals(key: str) -> tuple[float, float]:
    """Return the worst tick residual of a panel's x and y axes, pt."""
    out = []
    for axis in axes(key):
        slope, intercept = np.polyfit(axis.positions, axis.values, 1)
        positions = (np.asarray(axis.values) - intercept) / slope
        out.append(float(np.max(np.abs(positions - np.asarray(axis.positions)))))
    return out[0], out[1]


def _inside(path: Painted, box: tuple[float, float, float, float]) -> bool:
    """Whether a path's vertices lie within a frame (l, b, r, t), within EDGE."""
    left, bottom, right, top = box
    v = path.vertices
    return bool(
        np.all(v[:, 0] >= left - EDGE)
        and np.all(v[:, 0] <= right + EDGE)
        and np.all(v[:, 1] >= bottom - EDGE)
        and np.all(v[:, 1] <= top + EDGE)
    )


def _strokes(key: str, method: Method) -> list[Painted]:
    """Return a panel's stroked paths of a curve's colour and width."""
    return [
        p
        for p in _paths()
        if p.stroke == method.colour
        and abs(p.width - method.width) < 1e-3
        and _inside(p, frame(key))
    ]


def _check_legend(key: str) -> None:
    """Check the legend: its samples top to bottom in :data:`LEGEND`'s order."""
    panel = PANELS[key]
    heights = []
    for method in LEGEND:
        samples = [p for p in _strokes(key, method) if len(p.vertices) == 2]
        if len(samples) != 1:
            msg = f"{panel.name}: {len(samples)} legend samples for {method.label}"
            raise ValueError(msg)
        (x0, y0), _ = samples[0].vertices
        if not any(
            abs(g.centre[1] - y0) < LABEL and 0 < x0 - g.right < LABEL_REACH
            for g in _glyphs()
        ):
            msg = f"{panel.name}: no label beside the {method.label} sample"
            raise ValueError(msg)
        heights.append(float(y0))
    if heights != sorted(heights, reverse=True):
        msg = f"{panel.name}: the legend's samples are not in the order {LEGEND}"
        raise ValueError(msg)


def _chain(pieces: Sequence[np.ndarray], name: str) -> np.ndarray:
    """Join a curve's pieces end to end: each next one starts where the last ended."""
    rest = sorted(pieces, key=lambda p: float(p[0, 0]))
    line = rest.pop(0)
    while rest:
        nexts = [k for k, p in enumerate(rest) if np.hypot(*(p[0] - line[-1])) < JOIN]
        if len(nexts) != 1:
            msg = f"{name}: {len(nexts)} pieces continue the curve at {line[-1]}"
            raise ValueError(msg)
        line = np.vstack([line, rest.pop(nexts[0])[1:]])
    back = np.diff(line[:, 0])
    if back.min() < -BACKSTEP or line[-1, 0] < line[:, 0].max():
        msg = f"{name}: the curve steps back in time by {-back.min():.3f} pt"
        raise ValueError(msg)
    return line


def backsteps(key: str, method: Method) -> list[tuple[float, float]]:
    """Return where a curve steps back in time: (time in cycles, step in pt)."""
    line = curve(key, method)
    x_axis, _ = axes(key)
    steps = np.diff(line[:, 0])
    return [
        (float(x_axis(line[k, 0])), float(steps[k]))
        for k in np.flatnonzero(steps < -1e-9)
    ]


@functools.cache
def curve(key: str, method: Method) -> np.ndarray:
    """Return a panel's curve of one method, page points, in drawing order."""
    _check_legend(key)
    name = f"{PANELS[key].name} {method.label}"
    pieces = [p.vertices for p in _strokes(key, method) if len(p.vertices) > 2]
    if not pieces:
        msg = f"{name}: no curve"
        raise ValueError(msg)
    return _chain(pieces, name)


def end(key: str, method: Method) -> tuple[float, float, float, float]:
    """Return a curve at its end: x, y (pt), time (cycles), probability."""
    x, y = curve(key, method)[-1]
    x_axis, y_axis = axes(key)
    return float(x), float(y), float(x_axis(x)), float(y_axis(y))


def sigma(key: str) -> float:
    """Return the digitising uncertainty of a panel's probabilities."""
    _, y_axis = axes(key)
    return DIGITISING * y_axis.resolution()


def tail(key: str, method: Method) -> tuple[float, float]:
    """Return a curve's lowest and highest value over its last :data:`TAIL` cycle."""
    line = curve(key, method)
    x_axis, y_axis = axes(key)
    values = np.asarray(y_axis(line[:, 1]), dtype=float)
    times = np.asarray(x_axis(line[:, 0]), dtype=float)
    last = values[times >= times[-1] - TAIL]
    return float(last.min()), float(last.max())


def checks() -> dict[str, object]:
    """Return the checks of the reading; raise if one fails."""
    out: dict[str, object] = {}
    failures = []
    right = ticks(X_FROM)[0][-1]
    for key, panel in PANELS.items():
        x_res, y_res = residuals(key)
        out[f"{panel.name} worst tick residuals (pt)"] = (x_res, y_res)
        if QUANTUM / 2 + y_res > DIGITISING:
            failures.append(f"{panel.name}: y ticks beyond {DIGITISING} pt")
        zero = ticks(key)[1][0]
        for method in LEGEND:
            line = curve(key, method)
            x0, y0 = line[0]
            x_end = line[-1, 0]
            out[f"{panel.name} {method.label} start (pt from 1 cycle, 0)"] = (
                float(x0 - ticks(X_FROM)[0][0]),
                float(y0 - zero),
            )
            if abs(x0 - ticks(X_FROM)[0][0]) > EDGE or abs(y0 - zero) > START:
                failures.append(f"{panel.name} {method.label}: starts at {line[0]}")
            if abs(x_end - right) > EDGE:
                failures.append(f"{panel.name} {method.label}: ends at {x_end} pt")
            if key == "fig3b":
                out[f"{panel.name} {method.label} last {TAIL} cycle"] = tail(
                    key, method
                )
            if backsteps(key, method):
                out[f"{panel.name} {method.label} steps back (cycles, pt)"] = backsteps(
                    key, method
                )
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return each curve's value at the end of the pulse, after the checks."""
    checks()
    files = {}
    for key, panel in PANELS.items():
        for method in LEGEND:
            x, y, cycles, value = end(key, method)
            files[file(key, method)] = [
                Point(
                    "I6e14",
                    INTENSITY,
                    f"{value:.4g}",
                    "",
                    f"{panel.name} ({panel.ion}) {method.label}: the curve's last "
                    f"vertex at ({x:.3f}, {y:.3f}) pt, {cycles:.4f} optical cycles",
                )
            ]
    return files


def main() -> None:
    """Write the data files, and print the readings and the checks."""
    write(extract())
    for key, panel in PANELS.items():
        for method in LEGEND:
            _, _, cycles, value = end(key, method)
            print(
                f"{panel.name} {method.label}: {value:.5f} at {cycles:.4f} cycles "
                f"+- {sigma(key):.2g}"
            )
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
