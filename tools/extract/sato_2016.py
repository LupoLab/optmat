# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Sato et al. (2016): the TD-CASSCF survival probability of He at 800 nm.

The paper prints no number of its 800 nm survival probabilities: they are given only
in **Fig. 6** (PDF p. 7, library figure ``fig_13d7c9bcd599``), the probability of
finding both electrons inside a sphere of radius R0 = 20 a.u. against time in optical
cycles, for 4, 8 and 14 x 10^14 W/cm^2, 0 to 6 cycles, with the TD-CASSCF (0,0,14)
curve (MCTDHF for He: all electrons active) and the authors' ADK and PPT curves.
**Fig. 4(c)** (PDF p. 6, ``fig_bd7f3ee6e2fd``) draws the same survival at 8 x 10^14
W/cm^2 against time in atomic units for (0,0,1) to (0,0,14) with L_max = 59; it is
read only as a check of Fig. 6(b) and of the pulse's length. Both are vector drawings
on the page (no image), read from the drawing (data rule 8, TC-699) with
:func:`_figure.painted`; their tick labels are text, read with :func:`_figure.texts`.

How the drawing is read:

* the frames: a black stroked closed rectangle (five vertices) per panel, near the
  position :data:`PANELS` gives; the ticks black stroked two-vertex lines shorter than
  :data:`TICK_MAX`, perpendicular to an edge with one end on it (within :data:`EDGE`).
  The top and right edges repeat the bottom and left ones' ticks (checked, within
  :data:`MIRROR`). Fig. 6's ticks are all long (3.46 pt) and labelled on the left
  and, below panel (c) only, at the bottom; panels (a) and (b) share (c)'s x ticks
  (checked, within :data:`MIRROR`). Fig. 4(c)'s long ticks (3.49 pt) are labelled and
  four short ones (1.72 pt) divide each interval (checked to be a whole, constant
  number).
* the axes: each long tick's label is the text object centred within :data:`LABEL`
  of it along the axis, outside the frame, and must be the only one; the short ticks
  take the values of the equal steps between. Every tick lies on its fitted linear
  scale within :data:`TICK_TOLERANCE` (:meth:`_figure.Axis.check`), the drawing's
  quantum (its coordinates step by about 0.058 pt).
* the curves: the (0,0,14) curve is the black stroked line of width
  :data:`FIG6_WIDTH` (Fig. 6) with more than two vertices, inside the frame; the
  legend's sample of that width (two vertices) is beside the label '(0,0,14)' (checked,
  within :data:`LABEL`). In Fig. 4(c) it is the only thick black line
  (:data:`FIG4C_WIDTH`, drawn in pieces joined end to end; the legend's sample beside
  '(0,0,14)' is 1.108 pt wide in the drawing's units, the other black line, (0,0,9),
  0.346 pt and dashed). The polyline must never step back in time.
* the value: the curve's last vertex. In Fig. 6 it lies on the 6-cycle tick, the end
  of the pulse (tau = 6T, Eq. 27); in Fig. 4(c) at 661.6 au, 6T = 661.9 au (checked,
  within :data:`END_AU`). The data file holds 1 - S, the probability that at least one
  electron is outside R0 at t = 6T, with S (the survival, as drawn) in ``printed``.

The digitising uncertainty of S is :data:`DIGITISING` pt on the y axis (one quantum,
which bounds half the quantum plus the worst tick residual, checked): 2.7e-7 in (a),
1.3e-5 in (b) and 1.3e-4 in (c), as the axes span 4e-4, 0.02 and 0.2 on the same
90.4 pt. It is the same absolute uncertainty on 1 - S, the ``sigma`` of each row, and
0.08 to 0.09 % of it at every intensity.

**Checks** (:func:`checks`; the script fails if one is missed by more than its
margin), against every number the paper prints for these curves and between the two
figures:

* the caption's multipliers: 'the ADK rates must be scaled by 2.55, 1.35, and 0.96 ...
  The analogous factors for the PPT rates are 1.24, 1.27, and 1.24' for 'the survival
  probabilities to match those of TD-CASSCF'. With S = exp(-integral of the rate),
  the multiplier is ln S(TD-CASSCF) / ln S(model) at t = 6T: read from the red (ADK)
  and green dashed (PPT) curves, whose legend samples sit beside 'ADK' and 'PPT'
  (within :data:`MULTIPLIER`, relative). The ADK and PPT curves are read for this
  check only.
* Fig. 4(c)'s (0,0,14) at its end against Fig. 6(b) at 6 cycles (within
  :data:`TWICE` in S; L_max = 59 there, 71 in Fig. 6). The time axis of Fig. 4(c)
  puts the curve's end at 6T, the pulse's length the text gives for 800 nm and 8 x
  10^14 W/cm^2 (p. 5), which Fig. 6's caption does not restate.
"""

import functools
import itertools
import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from _figure import Axis, Painted, painted, texts
from _library import Point, write

DIRECTORY = "Sato_2016_TimeDependentComplete"
SOURCE = "sato_2016"
BLACK = (0, 0, 0)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
TICK_MAX = 5.0
"""Ticks are shorter than this, pt (Fig. 6's 3.46, Fig. 4(c)'s 3.49 and 1.72)."""
LONG = 2.5
"""Ticks longer than this are the long, labelled ones, pt."""
EDGE = 0.02
"""A tick's end lies on its frame edge within this, pt."""
MIRROR = 0.02
"""A top or right tick (or another panel's x tick) repeats one within this, pt."""
FRAME_FIND = 2.0
"""A frame's corners lie within this of :data:`PANELS`' position, pt."""
LABEL = 1.0
"""A label's centre lies within this of its tick along the axis, pt."""
LABEL_REACH = 30.0
"""A label lies within this of the frame, outside it, pt."""
TICK_TOLERANCE = 0.06
"""Every tick lies on its fitted scale within this, pt: the drawing's quantum."""
QUANTUM = 0.058
"""The step of the drawing's coordinates, pt: 0.055 in Fig. 6, 0.058 in Fig. 4(c)
(whose ticks are 3.657 and 3.599 pt apart for equal steps)."""
DIGITISING = 0.06
"""The digitising uncertainty of a vertex on the y axis, pt: one quantum, which bounds
half the quantum plus the worst tick residual (checked; 0.029 + 0.027 in Fig. 4(c),
0.028 + 0.001 in Fig. 6), and covers the curves' last vertices stepping one quantum
up and down on the flat end."""
JOIN = 1e-3
"""A curve's next piece starts where the last ended within this, pt."""
FIG6_WIDTH = 1.1
"""The stroke width of Fig. 6's (0,0,14) curve, pt (ADK and PPT 0.825, ticks 0.55)."""
FIG4C_WIDTH = 1.385
"""The stroke width of Fig. 4(c)'s (0,0,14) curve, in the drawing's units."""
OMEGA_AU = 45.563353 / 800.0
"""The photon energy at 800 nm, hartree: 1e7 / (800 nm) cm^-1 over 219474.63 cm^-1."""
SIX_T_AU = 6 * 2 * math.pi / OMEGA_AU
"""The pulse's length, tau = 6T, au (661.9)."""
END_AU = 1.0
"""Fig. 4(c)'s curve ends at 6T within this, au (its x quantum is 0.2 au)."""
MULTIPLIER = 0.02
"""The caption's multipliers are reproduced within this, relative: their two
decimals (0.2 to 0.5 %) and the reading of S near 1 in (a) (2.2e-7 in 1.2e-4 for ADK,
0.2 %), with room for the ADK/PPT curves' time shift ('shifted in time', the
caption), which is not described."""
TWICE = 3e-4
"""Fig. 4(c) and Fig. 6(b) agree within this in S (2 % of 1 - S): L_max 59 against
71 and the two drawings' quanta (5e-5 at most)."""


@dataclass(frozen=True)
class Panel:
    """One plot: where it is, what it holds."""

    name: str
    page: int
    frame: tuple[float, float, float, float]
    """The frame's left, bottom, right, top, pt, to within :data:`FRAME_FIND`."""
    width: float
    """The stroke width of the (0,0,14) curve."""
    x_from: str | None = None
    """The key of the panel whose labelled x axis this one shares, or None."""
    intensity: str = ""
    """The peak intensity, W/cm^2, as the caption prints it."""


PANELS = {
    "fig6a": Panel(
        "Fig. 6(a)", 7, (356.08, 644.80, 552.53, 735.63), FIG6_WIDTH, "fig6c", "4e14"
    ),
    "fig6b": Panel(
        "Fig. 6(b)", 7, (356.08, 545.52, 552.53, 636.34), FIG6_WIDTH, "fig6c", "8e14"
    ),
    "fig6c": Panel(
        "Fig. 6(c)", 7, (356.08, 446.23, 552.53, 537.06), FIG6_WIDTH, None, "14e14"
    ),
    "fig4c": Panel("Fig. 4(c)", 6, (357.20, 330.79, 555.16, 422.32), FIG4C_WIDTH),
}
FIG6 = ("fig6a", "fig6b", "fig6c")
ADK = {"fig6a": 2.55, "fig6b": 1.35, "fig6c": 0.96}
"""The caption's ADK multipliers, in ascending intensity."""
PPT = {"fig6a": 1.24, "fig6b": 1.27, "fig6c": 1.24}
"""The caption's PPT multipliers."""
FILE = "He/sato_2016__fig6_casscf_He_800nm.csv"
ROWS = {"fig6a": "I4e14", "fig6b": "I8e14", "fig6c": "I14e14"}


@functools.cache
def _paths(page: int) -> tuple[Painted, ...]:
    """Return a page's painted paths (the whole page)."""
    return tuple(painted(DIRECTORY, page, (0.0, 0.0, 1000.0, 1000.0)))


@functools.cache
def _texts(page: int) -> tuple[tuple[str, float, float], ...]:
    """Return a page's text objects: text and centre, pt."""
    return tuple(texts(DIRECTORY, page, (0.0, 0.0, 1000.0, 1000.0)))


@functools.cache
def frame(key: str) -> tuple[float, float, float, float]:
    """Return a panel's frame (l, b, r, t, pt): its black stroked closed rectangle."""
    panel = PANELS[key]
    found = set()
    for path in _paths(panel.page):
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
            found.add(tuple(round(c, 3) for c in box))
    if len(found) != 1:
        msg = f"{panel.name}: {len(found)} frames near {panel.frame}"
        raise ValueError(msg)
    left, bottom, right, top = found.pop()
    return left, bottom, right, top


@dataclass(frozen=True)
class Ticks:
    """A panel's tick positions along one edge, pt: every one and the long ones."""

    every: tuple[float, ...]
    long: tuple[float, ...]


@functools.cache
def ticks(key: str) -> tuple[Ticks, Ticks]:
    """Return the x ticks (bottom edge) and y ticks (left edge) of a panel.

    Checks that the top and right edges repeat them within :data:`MIRROR`.
    """
    panel = PANELS[key]
    left, bottom, right, top = frame(key)
    edges: dict[str, list[tuple[float, float]]] = {
        "bottom": [],
        "top": [],
        "left": [],
        "right": [],
    }
    for path in _paths(panel.page):
        if path.stroke != BLACK or len(path.vertices) != 2:
            continue
        (x0, y0), (x1, y1) = path.vertices
        length = math.hypot(x1 - x0, y1 - y0)
        if length > TICK_MAX:
            continue
        if abs(x0 - x1) < EDGE and left - EDGE <= x0 <= right + EDGE:
            for name, edge in (("bottom", bottom), ("top", top)):
                if min(abs(y0 - edge), abs(y1 - edge)) < EDGE:
                    edges[name].append((float(x0), length))
        if abs(y0 - y1) < EDGE and bottom - EDGE <= y0 <= top + EDGE:
            for name, edge in (("left", left), ("right", right)):
                if min(abs(x0 - edge), abs(x1 - edge)) < EDGE:
                    edges[name].append((float(y0), length))
    out = []
    for near, far in (("bottom", "top"), ("left", "right")):
        every = sorted({round(p, 3) for p, _ in edges[near]})
        long = sorted({round(p, 3) for p, length in edges[near] if length > LONG})
        mirror = sorted({round(p, 3) for p, _ in edges[far]})
        if (
            not long
            or len(mirror) != len(every)
            or any(min(abs(m - p) for p in every) > MIRROR for m in mirror)
        ):
            msg = f"{panel.name}: the {far} ticks do not repeat the {near} ones"
            raise ValueError(msg)
        out.append(Ticks(tuple(every), tuple(long)))
    return out[0], out[1]


def _label(key: str, position: float, along: str) -> float:
    """Return the value printed at a long tick: the one text centred on it, outside."""
    panel = PANELS[key]
    left, bottom, _, _ = frame(key)
    found = []
    for text, x, y in _texts(panel.page):
        try:
            value = float(text)
        except ValueError:
            continue
        if along == "x" and abs(x - position) < LABEL and 0 < bottom - y < LABEL_REACH:
            found.append(value)
        if along == "y" and abs(y - position) < LABEL and 0 < left - x < LABEL_REACH:
            found.append(value)
    if len(found) != 1:
        msg = f"{panel.name}: {len(found)} labels at the {along} tick {position}"
        raise ValueError(msg)
    return found[0]


def _values(key: str, ticks_: Ticks, along: str) -> list[float]:
    """Return every tick's value: the long ones' labels, the short ones' equal steps."""
    name = PANELS[key].name
    labels = [_label(key, p, along) for p in ticks_.long]
    between = [
        sum(1 for p in ticks_.every if a < p < b)
        for a, b in itertools.pairwise(ticks_.long)
    ]
    if len(set(between)) > 1:
        msg = f"{name}: the short {along} ticks do not divide the intervals evenly"
        raise ValueError(msg)
    steps = (between[0] if between else 0) + 1
    step = (labels[1] - labels[0]) / steps
    if any(
        abs((b - a) - steps * step) > 1e-9 * max(1.0, abs(step))
        for a, b in itertools.pairwise(labels)
    ):
        msg = f"{name}: the {along} labels are not equally spaced"
        raise ValueError(msg)
    # Each tick's count of steps from the first long tick, from its position: the
    # long ticks' spacing is the scale, and every tick is then checked on it.
    pitch = (ticks_.long[1] - ticks_.long[0]) / steps
    return [
        labels[0] + round((p - ticks_.long[0]) / pitch) * step for p in ticks_.every
    ]


@functools.cache
def axes(key: str) -> tuple[Axis, Axis]:
    """Return a panel's time and survival axes, checked on every tick.

    The time axis is in optical cycles (Fig. 6) or atomic units (Fig. 4(c)).
    """
    panel = PANELS[key]
    x_ticks, y_ticks = ticks(key)
    if panel.x_from is not None:
        shared, _ = ticks(panel.x_from)
        if len(shared.every) != len(x_ticks.every) or any(
            abs(a - b) > MIRROR
            for a, b in zip(shared.every, x_ticks.every, strict=True)
        ):
            msg = f"{panel.name}: its x ticks are not {PANELS[panel.x_from].name}'s"
            raise ValueError(msg)
        x_axis, _ = axes(panel.x_from)
    else:
        x_axis = Axis(
            x_ticks.every,
            tuple(_values(key, x_ticks, "x")),
            tolerance=TICK_TOLERANCE,
        )
    y_axis = Axis(
        y_ticks.every, tuple(_values(key, y_ticks, "y")), tolerance=TICK_TOLERANCE
    )
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
    if np.any(np.diff(line[:, 0]) < -1e-9):
        msg = f"{name}: the curve steps back in time"
        raise ValueError(msg)
    return line


def _inside(path: Painted, box: tuple[float, float, float, float]) -> bool:
    """Whether a path lies within a frame (l, b, r, t), within :data:`EDGE`."""
    left, bottom, right, top = box
    return (
        path.left >= left - EDGE - path.width
        and path.right <= right + EDGE + path.width
        and path.bottom >= bottom - EDGE - path.width
        and path.top <= top + EDGE + path.width
    )


def _strokes(key: str, colour: tuple[int, int, int], width: float) -> list[Painted]:
    """Return a panel's stroked paths of one colour and width inside its frame."""
    panel = PANELS[key]
    return [
        p
        for p in _paths(panel.page)
        if p.stroke == colour and abs(p.width - width) < 1e-3 and _inside(p, frame(key))
    ]


def _legend(key: str, colour: tuple[int, int, int], width: float, label: str) -> None:
    """Check that a curve's legend sample (two vertices) sits beside its label."""
    panel = PANELS[key]
    samples = [
        float(p.vertices[0, 1])
        for p in _strokes(key, colour, width)
        if len(p.vertices) == 2
    ]
    if colour == GREEN:  # the dashed PPT line: the sample is its path's first dash
        samples = [
            float(p.vertices[0, 1])
            for p in _strokes(key, colour, width)
            if abs(p.vertices[0, 1] - p.vertices[1, 1]) < EDGE
            and p.vertices[1, 0] - p.vertices[0, 0] > 15.0
        ]
    labels = [y for text, _, y in _texts(panel.page) if text == label]
    if len(samples) != 1 or not any(abs(y - samples[0]) < LABEL for y in labels):
        msg = f"{panel.name}: no legend sample beside '{label}'"
        raise ValueError(msg)


@functools.cache
def curve(key: str) -> np.ndarray:
    """Return a panel's (0,0,14) curve, page points, in drawing order."""
    panel = PANELS[key]
    pieces = [
        p.vertices for p in _strokes(key, BLACK, panel.width) if len(p.vertices) > 2
    ]
    if not pieces:
        msg = f"{panel.name}: no (0,0,14) curve"
        raise ValueError(msg)
    if key in FIG6:
        _legend(key, BLACK, panel.width, "(0,0,14)")
    return _chain(pieces, panel.name)


def end(key: str) -> tuple[float, float, float, float]:
    """Return a panel's curve at its end: x, y (pt), time, survival S."""
    x, y = curve(key)[-1]
    x_axis, y_axis = axes(key)
    return float(x), float(y), float(x_axis(x)), float(y_axis(y))


def _model_end(key: str, colour: tuple[int, int, int], label: str) -> float:
    """Return the survival S of the ADK (red) or PPT (green) curve at 6 cycles."""
    _legend(key, colour, 0.825, label)
    vertices = np.vstack([p.vertices for p in _strokes(key, colour, 0.825)])
    last = vertices[vertices[:, 0] >= vertices[:, 0].max() - EDGE]
    if np.ptp(last[:, 1]) > EDGE:
        msg = f"{PANELS[key].name}: the {label} curve ends at more than one value"
        raise ValueError(msg)
    _, y_axis = axes(key)
    return float(y_axis(float(last[0, 1])))


def sigma(key: str) -> float:
    """Return the digitising uncertainty of a panel's survival (and of 1 - S)."""
    _, y_axis = axes(key)
    return DIGITISING * y_axis.resolution()


def checks() -> dict[str, object]:
    """Return the checks of the reading; raise if one fails its margin."""
    out: dict[str, object] = {}
    failures = []
    for key, panel in PANELS.items():
        x_res, y_res = residuals(key)
        out[f"{panel.name} worst tick residuals (pt)"] = (x_res, y_res)
        if QUANTUM / 2 + y_res > DIGITISING:
            failures.append(f"{panel.name}: y ticks beyond {DIGITISING} pt")
    for key in FIG6:
        x, _, cycles, survival = end(key)
        x_ticks, _ = ticks("fig6c")
        if abs(x - x_ticks.every[-1]) > EDGE or abs(cycles - 6.0) > 1e-3:
            failures.append(f"{PANELS[key].name}: the curve ends at {cycles} cycles")
        for name, colour, printed in (("ADK", RED, ADK), ("PPT", GREEN, PPT)):
            model = _model_end(key, colour, name)
            read = math.log(survival) / math.log(model)
            out[f"{PANELS[key].name} {name} multiplier"] = (read, printed[key], model)
            if abs(read / printed[key] - 1) > MULTIPLIER:
                failures.append(
                    f"{PANELS[key].name}: {name} multiplier {read:.3f}, "
                    f"printed {printed[key]}"
                )
    _, _, t_end, s4 = end("fig4c")
    out["Fig. 4(c) curve ends at (au)"] = (t_end, SIX_T_AU)
    if abs(t_end - SIX_T_AU) > END_AU:
        failures.append(f"Fig. 4(c) ends at {t_end:.1f} au, not 6T")
    _, _, _, s6 = end("fig6b")
    out["Fig. 4(c) S over Fig. 6(b) S"] = (s4, s6, s4 - s6)
    if abs(s4 - s6) > TWICE:
        failures.append(f"Fig. 4(c) S {s4:.6f} against Fig. 6(b) {s6:.6f}")
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return the TD-CASSCF ionisation probability 1 - S at t = 6T, after the checks."""
    checks()
    points = []
    for key in FIG6:
        x, y, _, survival = end(key)
        points.append(
            Point(
                ROWS[key],
                PANELS[key].intensity,
                f"{1 - survival:.4g}",
                f"{sigma(key):.2g}",
                f"{PANELS[key].name} (0,0,14): survival {survival:.7f} at t = 6T, "
                f"the curve's last vertex at ({x:.3f}, {y:.3f}) pt",
            )
        )
    return {FILE: points}


def main() -> None:
    """Write the data file, and print the readings and the checks."""
    write(extract())
    for key in FIG6:
        _, _, _, survival = end(key)
        print(
            f"{PANELS[key].name}: S = {survival:.7f}, 1 - S = {1 - survival:.4e} "
            f"+- {sigma(key):.2g} ({sigma(key) / (1 - survival):.2%})"
        )
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
