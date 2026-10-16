# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Pathak, Sato and Ishikawa (2020): single-ionisation probabilities of Ar.

The paper prints no number of its time-dependent calculations: the probabilities are
only in **Fig. 2** (PDF p. 7, journal p. 034110-6, library figure ``fig_37c2c1bdbe02``,
whose crops lost panel (a); the PDF has all three), "Time evolution of single ionization
probability of Ar irradiated by a laser pulse with a wavelength of 1200 nm and a peak
intensity of 1 x 10^14 W/cm^2 (a), 2 x 10^14 W/cm^2 (b), and 4 x 10^14 W/cm^2 (c)
calculated with TDHF, TD-OMP2, TD-OCCD, TD-CASSCF, and SAE methods", a vector drawing
(the page has no image), read from the drawing (data rule 8, TC-699).

How the drawing is read:

* the frames: per panel four thin filled rectangles (0.8 pt) found near
  :data:`PANELS`' boxes; the ticks are two-vertex paths (each with a 0.8 pt filled
  rectangle around it) with one end on a frame line, inward, on all four sides; the top
  and right ticks repeat the bottom and left ones (checked, within :data:`MIRROR`).
* the axes: the tick labels are glyphs drawn as filled paths, not text; their values
  are written in :data:`X_LABELS` and :data:`PANELS` as read on the rendered page. Time
  (a.u.) 0 to 500 in steps of 100; probability (a) 0.00000 to 0.00075 in steps of
  0.00025, (b) 0.00 to 0.03 in 0.01, (c) 0.0 to 0.4 in 0.1, all linear. Every tick lies
  on the fitted scale within :data:`TOLERANCE`.
* the curves: each method's curve is one filled outline (the stroke converted to an
  outline polygon of lines and cubic Bezier segments), coloured as the legend of panel
  (c) shows: TDHF blue, TD-OMP2 green, TD-OCCD red, TD-CASSCF black (thick, one closed
  outline), SAE black dotted (one path of 58 to 71 closed squares, one per dot). The
  legend's sample lines are checked to stand in that order, top to bottom
  (:data:`LEGEND`). The Beziers are flattened (:data:`BEZIER_STEPS` per segment);
  the line's centre at a time is the mean of the outline's two crossings of that
  vertical (the midpoint of a stroke's vertical chord is on its centre line); the
  dotted line's centre is interpolated between the centres of the dots on either side.
* the times: the end of the pulse, t = 3T = 496.44 au (Eq. 29; T = 2 pi / omega_0 at
  1200 nm, 165.48 au), and the end of each drawn curve, which the script finds from the
  drawing: the solid outlines end half a line width (a projecting cap) beyond the line,
  the dots end with the line; every curve ends at 546.08 au = 3.3T within
  :data:`END_TOLERANCE` (the paper: "The simulations are continued after the end of
  pulse for further 6000 time steps", at 20 000 steps per cycle: 0.3T).

The value at 3T is the dataset's value; the value at the curve's end, 3.3T, is in each
point's ``printed`` cell. Each curve also starts at zero (checked at t = 20 au, within
:data:`ZERO`).

Checks against the text (p. 7), the only statements about Fig. 2: "a systematic
overestimation by TD-OMP2 and underestimation by TDHF in comparison to the TD-CASSCF
result" (at every intensity, both times), and at the highest intensity "the SAE result
is as accurate as the TD-OCCD one" (there SAE and TD-OCCD are both within 15 % of
TD-CASSCF, while TDHF is 40 % below it).
"""

import ctypes
import itertools
from dataclasses import dataclass
from typing import Final

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, Drawn, painted, pdf
from _library import Point, write
from pypdfium2 import raw

DIRECTORY: Final = "Pathak_2020_TimeDependentOptimized"
SOURCE: Final = "pathak_2020"
PAGE: Final = 7
FIGURE: Final = "fig_37c2c1bdbe02"

# Time: the period of 1200 nm light in atomic units of time (CODATA 2018).
SPEED_OF_LIGHT: Final = 299_792_458.0
AU_TIME_S: Final = 2.4188843265857e-17
PERIOD_AU: Final = 1200e-9 / SPEED_OF_LIGHT / AU_TIME_S
PULSE_END_AU: Final = 3 * PERIOD_AU
CURVE_END_AU: Final = 3.3 * PERIOD_AU


@dataclass(frozen=True)
class Panel:
    """One panel of Fig. 2: its intensity, frame and y tick labels."""

    intensity: str
    """The peak intensity, W/cm^2, as the caption prints it."""
    box: tuple[float, float, float, float]
    """The frame's outer edges (left, bottom, right, top), page points, approximate."""
    y_labels: tuple[float, ...]
    """The y tick labels, bottom to top, as read on the rendered page."""


PANELS: Final = {
    "a": Panel(
        "1e14", (356.36, 620.26, 536.35, 704.49), (0.0, 0.00025, 0.0005, 0.00075)
    ),
    "b": Panel("2e14", (356.36, 520.14, 536.35, 604.37), (0.0, 0.01, 0.02, 0.03)),
    "c": Panel("4e14", (356.36, 420.01, 536.35, 504.25), (0.0, 0.1, 0.2, 0.3, 0.4)),
}
X_LABELS: Final = (0.0, 100.0, 200.0, 300.0, 400.0, 500.0)
# How far a frame may be from where PANELS says, and a tick's end from its frame line.
NEAR: Final = 1.0
# The largest residual of a tick on its axis's fitted scale, and between a tick and its
# mirror on the opposite edge, page points.
TOLERANCE: Final = 0.01
MIRROR: Final = 0.01

# The methods, in the order of the legend (panel (c)), top to bottom, with their paint.
LEGEND: Final = ("TDHF", "TD-OMP2", "TD-OCCD", "TD-CASSCF", "SAE")
BLUE: Final = (0, 0, 255)
GREEN: Final = (0, 128, 0)
RED: Final = (255, 0, 0)
BLACK: Final = (0, 0, 0)
COLOURS: Final = {"TDHF": BLUE, "TD-OMP2": GREEN, "TD-OCCD": RED}
# A curve spans the time axis; the legend's sample lines are 17 pt long.
CURVE_WIDTH: Final = 100.0
# The legend's sample lines: left and right ends, page points (panel (c)).
LEGEND_SPAN: Final = (362.0, 382.0)
# The dotted line has dozens of dots; the solid outline is one closed subpath.
MIN_DOTS: Final = 20
# The dots a dotted line's period and dot length are measured over: its last ten.
DOTS_PERIOD: Final = 10
BEZIER_STEPS: Final = 16
END_TOLERANCE: Final = 0.5  # au
ZERO: Final = 0.05  # pt, the curves' distance from the zero tick at t = 20 au
ZERO_TIME: Final = 20.0
# The text's statements, made numbers: at 4e14, SAE and TD-OCCD within 15 % of
# TD-CASSCF, TDHF more than 30 % below it.
NEAR_CASSCF: Final = 0.15
TDHF_BELOW: Final = 0.30

DATASETS: Final = {
    "TDHF": f"Ar/{SOURCE}__fig2_tdhf_Ar.csv",
    "TD-OMP2": f"Ar/{SOURCE}__fig2_tdomp2_Ar.csv",
    "TD-OCCD": f"Ar/{SOURCE}__fig2_tdoccd_Ar.csv",
    "TD-CASSCF": f"Ar/{SOURCE}__fig2_tdcasscf_Ar.csv",
    "SAE": f"Ar/{SOURCE}__fig2_sae_Ar.csv",
}


# --- the drawing ---------------------------------------------------------------


@dataclass(frozen=True)
class Outline:
    """A filled path: its fill colour and its subpaths, flattened, in page points."""

    fill: tuple[int, int, int] | None
    subpaths: tuple[np.ndarray, ...]

    @property
    def left(self) -> float:
        """The leftmost x of the outline."""
        return min(float(s[:, 0].min()) for s in self.subpaths)

    @property
    def right(self) -> float:
        """The rightmost x of the outline."""
        return max(float(s[:, 0].max()) for s in self.subpaths)

    @property
    def bottom(self) -> float:
        """The lowest y of the outline."""
        return min(float(s[:, 1].min()) for s in self.subpaths)

    @property
    def top(self) -> float:
        """The highest y of the outline."""
        return max(float(s[:, 1].max()) for s in self.subpaths)


def _container_to_page(obj: pdfium.PdfObject, points: np.ndarray) -> np.ndarray:
    """Map points from an object's container (form XObjects) to the page."""
    out = np.array(points, dtype=float).reshape(-1, 2)
    container = obj.container
    while container is not None:
        matrix = container.get_matrix()
        out = np.array([matrix.on_point(x, y) for x, y in out], dtype=float)
        out = out.reshape(-1, 2)
        container = container.container
    return out


def _fill(obj: pdfium.PdfObject) -> tuple[int, int, int] | None:
    """Return a path's fill colour, or None if it is not filled."""
    fill_mode, stroked = ctypes.c_int(), ctypes.c_int()
    raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill_mode), ctypes.byref(stroked))
    if not fill_mode.value:
        return None
    channels = [ctypes.c_uint() for _ in range(4)]
    if not raw.FPDFPageObj_GetFillColor(obj.raw, *(ctypes.byref(c) for c in channels)):
        return None
    return channels[0].value, channels[1].value, channels[2].value


def _flatten(obj: pdfium.PdfObject) -> tuple[np.ndarray, ...]:
    """Return a path's subpaths as polylines, the Bezier segments flattened."""
    matrix = obj.get_matrix()
    subpaths: list[list[tuple[float, float]]] = []
    pending: list[tuple[float, float]] = []  # control points of a Bezier segment
    for index in range(raw.FPDFPath_CountSegments(obj.raw)):
        segment = raw.FPDFPath_GetPathSegment(obj.raw, index)
        x, y = ctypes.c_float(), ctypes.c_float()
        raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
        point = matrix.on_point(x.value, y.value)
        kind = raw.FPDFPathSegment_GetType(segment)
        if kind == raw.FPDF_SEGMENT_MOVETO:
            subpaths.append([point])
        elif kind == raw.FPDF_SEGMENT_LINETO:
            subpaths[-1].append(point)
        else:  # a cubic Bezier: two control points, then the end point
            pending.append(point)
            if len(pending) == 3:
                start = np.array(subpaths[-1][-1])
                c1, c2, end = (np.array(p) for p in pending)
                for t in np.linspace(0, 1, BEZIER_STEPS + 1)[1:]:
                    b = (
                        (1 - t) ** 3 * start
                        + 3 * (1 - t) ** 2 * t * c1
                        + 3 * (1 - t) * t**2 * c2
                        + t**3 * end
                    )
                    subpaths[-1].append((float(b[0]), float(b[1])))
                pending = []
    return tuple(
        _container_to_page(obj, np.array(s, dtype=float)) for s in subpaths if s
    )


def outlines(region: tuple[float, float, float, float]) -> list[Outline]:
    """Return the filled paths of the page that touch a region (l, b, r, t, pt)."""
    left, bottom, right, top = region
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    found = []
    for obj in document[PAGE - 1].get_objects(max_depth=15):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        fill = _fill(obj)
        if fill is None:
            continue
        outline = Outline(fill, _flatten(obj))
        if not outline.subpaths:
            continue
        if (
            outline.right >= left
            and outline.left <= right
            and outline.top >= bottom
            and outline.bottom <= top
        ):
            found.append(outline)
    return found


# --- the axes ------------------------------------------------------------------


@dataclass(frozen=True)
class Frame:
    """A panel's calibrated axes, and its frame lines (page points)."""

    x: Axis
    y: Axis
    zero: float
    """The zero tick's y, page points."""


def _ticks(drawn: list[Drawn], box: tuple[float, float, float, float]) -> Frame | str:
    """Return a panel's axes from its tick marks, or what is wrong with them."""
    left, bottom, right, top = box
    # The frame lines' centres: 0.4 pt inside the outer edges of the 0.8 pt rectangles.
    edges = {"left": left + 0.4, "bottom": bottom + 0.4, "right": right - 0.4}
    edges["top"] = top - 0.4
    marks: dict[str, list[float]] = {name: [] for name in edges}
    for d in drawn:
        if len(d.vertices) != 2:
            continue
        (x0, y0), (x1, y1) = d.vertices
        if abs(x0 - x1) < 1e-3:  # vertical: on the bottom or top line
            for name in ("bottom", "top"):
                if min(abs(y0 - edges[name]), abs(y1 - edges[name])) < NEAR / 4:
                    marks[name].append(float(x0))
        elif abs(y0 - y1) < 1e-3:  # horizontal: on the left or right line
            for name in ("left", "right"):
                if min(abs(x0 - edges[name]), abs(x1 - edges[name])) < NEAR / 4:
                    marks[name].append(float(y0))
    for found in marks.values():
        found.sort()
    for a, b in (("bottom", "top"), ("left", "right")):
        if len(marks[a]) != len(marks[b]) or np.any(
            np.abs(np.subtract(marks[a], marks[b])) > MIRROR
        ):
            return f"the {a} and {b} ticks differ: {marks[a]} and {marks[b]}"
    return Frame(
        Axis(tuple(marks["bottom"]), X_LABELS, tolerance=TOLERANCE),
        Axis(tuple(marks["left"]), (), tolerance=TOLERANCE),
        marks["left"][0] if marks["left"] else float("nan"),
    )


def frame(letter: str) -> Frame:
    """Return a panel's axes, calibrated from its ticks, checked against every one."""
    panel = PANELS[letter]
    left, bottom, right, top = panel.box
    drawn: list[Drawn] = list(
        painted(DIRECTORY, PAGE, (left - 2, bottom - 2, right + 2, top + 2))
    )
    lines = [
        d
        for d in drawn
        if len(d.vertices) == 5
        and min(d.right - d.left, d.top - d.bottom) < 1
        and max(d.right - d.left, d.top - d.bottom) > 50
    ]
    found = (
        min(d.left for d in lines),
        min(d.bottom for d in lines),
        max(d.right for d in lines),
        max(d.top for d in lines),
    )
    if len(lines) != 4 or np.any(np.abs(np.subtract(found, panel.box)) > NEAR):
        msg = f"Fig. 2({letter}): frame {found} ({len(lines)} lines), not {panel.box}"
        raise ValueError(msg)
    axes = _ticks(drawn, found)
    if isinstance(axes, str):
        msg = f"Fig. 2({letter}): {axes}"
        raise ValueError(msg)
    if len(axes.x.positions) != len(X_LABELS) or len(axes.y.positions) != len(
        panel.y_labels
    ):
        msg = (
            f"Fig. 2({letter}): {len(axes.x.positions)} x and "
            f"{len(axes.y.positions)} y ticks, not {len(X_LABELS)} and "
            f"{len(panel.y_labels)}"
        )
        raise ValueError(msg)
    y_axis = Axis(axes.y.positions, panel.y_labels, tolerance=TOLERANCE)
    axes.x.check()
    y_axis.check()
    return Frame(axes.x, y_axis, axes.y.positions[0])


def residuals(axis: Axis) -> float:
    """Return the largest distance of a tick from its axis's fitted scale, pt."""
    slope, intercept = np.polyfit(axis.positions, axis.values, 1)
    predicted = (np.asarray(axis.values) - intercept) / slope
    return float(np.max(np.abs(predicted - np.asarray(axis.positions))))


def x_of(axis: Axis, time: float) -> float:
    """Return the page x of a time, au."""
    slope, intercept = np.polyfit(axis.positions, axis.values, 1)
    return float((time - intercept) / slope)


# --- the curves ----------------------------------------------------------------


def _crossings(polygon: np.ndarray, x0: float) -> list[float]:
    """Return the y of every crossing of a closed polygon's edges with x = x0."""
    closed = np.vstack([polygon, polygon[:1]])
    found = []
    for (xa, ya), (xb, yb) in itertools.pairwise(closed):
        if xa == xb or not min(xa, xb) <= x0 < max(xa, xb):
            continue
        found.append(float(ya + (yb - ya) * (x0 - xa) / (xb - xa)))
    return sorted(found)


@dataclass(frozen=True)
class Reading:
    """One curve's centre at one time: page coordinates and the value."""

    time: float
    x: float
    y: float
    value: float


@dataclass(frozen=True)
class Curve:
    """A method's curve in one panel: its drawing, kind and end (au)."""

    method: str
    outline: Outline
    dotted: bool

    def dots(self) -> list[tuple[float, float, float]]:
        """Return the dotted line's dots, (left, right, centre y), along the line."""
        return sorted(
            (
                float(s[:, 0].min()),
                float(s[:, 0].max()),
                float(s[:, 1].min() + s[:, 1].max()) / 2,
            )
            for s in self.outline.subpaths
        )

    def centre(self, x0: float) -> float:
        """Return the line's centre y at page x0."""
        if self.dotted:
            dots = self.dots()
            xs = [(a + b) / 2 for a, b, _ in dots]
            ys = [y for _, _, y in dots]
            # Beyond the last dot's centre (the line ends in a dot or a gap): that
            # dot's centre, the line being flat there.
            return float(np.interp(x0, xs, ys))
        (polygon,) = self.outline.subpaths
        found = _crossings(polygon, x0)
        if len(found) != 2:
            msg = f"{self.method}: {len(found)} crossings at x = {x0:.3f}, not 2"
            raise ValueError(msg)
        return (found[0] + found[1]) / 2

    def width(self, x0: float) -> float:
        """Return the solid line's vertical width at page x0 (its width where flat)."""
        (polygon,) = self.outline.subpaths
        found = _crossings(polygon, x0)
        return found[-1] - found[0]

    def end(self) -> tuple[float, float]:
        """Return the page x where the drawn line ends, and the latest it can end.

        A solid outline ends half the line's width (a projecting cap) beyond the line.
        A dotted line ends at its last dot or in the gap after it.
        """
        if self.dotted:
            dots = self.dots()[-DOTS_PERIOD:]
            period = float(np.median(np.diff([a for a, _, _ in dots])))
            length = float(np.median([b - a for a, b, _ in dots]))
            return dots[-1][1], dots[-1][1] + period - length
        right = self.outline.right
        end = right - self.width(right - 1.0) / 2
        return end, end


def curves(letter: str) -> dict[str, Curve]:
    """Return the five curves of a panel, identified by paint and kind."""
    left, bottom, right, top = PANELS[letter].box
    found: dict[str, Curve] = {}
    for outline in outlines((left, bottom, right, top)):
        if outline.right - outline.left < CURVE_WIDTH:
            continue
        if not (left < outline.left and outline.right < right):
            continue  # the frame and the legend's box
        if not (bottom < outline.bottom and outline.top < top):
            continue
        method = next((m for m, c in COLOURS.items() if c == outline.fill), None)
        dotted = False
        if outline.fill == BLACK:
            dotted = len(outline.subpaths) >= MIN_DOTS
            if not dotted and len(outline.subpaths) != 1:
                msg = (
                    f"Fig. 2({letter}): a black curve of {len(outline.subpaths)} parts"
                )
                raise ValueError(msg)
            method = "SAE" if dotted else "TD-CASSCF"
        if method is None:
            continue
        if method in found:
            msg = f"Fig. 2({letter}): two curves for {method}"
            raise ValueError(msg)
        found[method] = Curve(method, outline, dotted)
    if set(found) != set(LEGEND):
        msg = f"Fig. 2({letter}): curves {sorted(found)}, not {sorted(LEGEND)}"
        raise ValueError(msg)
    return found


def check_legend() -> None:
    """Raise unless panel (c)'s legend samples stand in LEGEND's order, top down."""
    left, bottom, right, top = PANELS["c"].box
    samples = [
        o
        for o in outlines((left, bottom, right, top))
        if LEGEND_SPAN[0] < o.left
        and o.right < LEGEND_SPAN[1]
        and o.right - o.left > 10
        and o.top - o.bottom < 2
    ]
    samples.sort(key=lambda o: -o.top)
    kinds = []
    for o in samples:
        method = next((m for m, c in COLOURS.items() if c == o.fill), None)
        if o.fill == BLACK:
            method = "SAE" if len(o.subpaths) > 1 else "TD-CASSCF"
        kinds.append(method)
    if tuple(kinds) != LEGEND:
        msg = f"Fig. 2(c): legend samples {kinds}, not {LEGEND}"
        raise ValueError(msg)


@dataclass(frozen=True)
class Result:
    """One method's readings in one panel."""

    pulse_end: Reading
    curve_end: Reading


def read(letter: str) -> dict[str, Result]:
    """Return each method's value at 3T and at the end of its curve, in a panel."""
    axes = frame(letter)
    x_pulse_end = x_of(axes.x, PULSE_END_AU)
    x_zero_time = x_of(axes.x, ZERO_TIME)
    results = {}
    for method, curve in curves(letter).items():
        zero = curve.centre(x_zero_time)
        if abs(zero - axes.zero) > ZERO:
            msg = f"Fig. 2({letter}) {method}: {zero - axes.zero:.3f} pt off zero"
            raise ValueError(msg)
        x_end, x_latest = curve.end()
        end_time = float(axes.x(x_end))
        latest = float(axes.x(x_latest))
        if not end_time - END_TOLERANCE <= CURVE_END_AU <= latest + END_TOLERANCE:
            msg = f"Fig. 2({letter}) {method}: ends at {end_time:.2f} au, not 3.3T"
            raise ValueError(msg)
        readings = []
        for time, x0 in ((PULSE_END_AU, x_pulse_end), (end_time, x_end)):
            y = curve.centre(x0)
            readings.append(Reading(time, x0, y, float(axes.y(y))))
        results[method] = Result(*readings)
    return results


def checks(found: dict[str, dict[str, Result]]) -> None:
    """Raise unless the readings agree with what the text says of Fig. 2."""
    for letter, results in found.items():
        for which in ("pulse_end", "curve_end"):
            value = {m: getattr(r, which).value for m, r in results.items()}
            if not value["TD-OMP2"] > value["TD-CASSCF"] > value["TDHF"]:
                msg = f"Fig. 2({letter}) {which}: not OMP2 > CASSCF > TDHF: {value}"
                raise ValueError(msg)
            if letter == "c":
                casscf = value["TD-CASSCF"]
                for method in ("SAE", "TD-OCCD"):
                    if abs(value[method] / casscf - 1) > NEAR_CASSCF:
                        msg = f"Fig. 2(c): {method} not within 15 % of TD-CASSCF"
                        raise ValueError(msg)
                if value["TDHF"] / casscf > 1 - TDHF_BELOW:
                    msg = "Fig. 2(c): TDHF not 30 % below TD-CASSCF"
                    raise ValueError(msg)


def readings() -> dict[str, dict[str, Result]]:
    """Return every panel's readings, checked."""
    check_legend()
    found = {letter: read(letter) for letter in PANELS}
    checks(found)
    return found


def extract() -> dict[str, list[Point]]:
    """Return each method's probabilities at 3T against the intensity, one file each."""
    found = readings()
    files: dict[str, list[Point]] = {}
    for method, relative in DATASETS.items():
        points = []
        for letter, panel in PANELS.items():
            result = found[letter][method]
            at, end = result.pulse_end, result.curve_end
            printed = (
                f"Fig. 2({letter}) {method}: line centre at ({at.x:.3f}, "
                f"{at.y:.3f}) pt, "
                f"t = 3T = {at.time:.2f} au; at the curve's last drawn point, t = "
                f"{end.time:.2f} au (the line ends at 3.3T = {CURVE_END_AU:.2f} au), "
                f"{end.value:.4g} at ({end.x:.3f}, {end.y:.3f}) pt"
            )
            points.append(
                Point(f"fig2{letter}", panel.intensity, f"{at.value:.4g}", "", printed)
            )
        files[relative] = points
    return files


def main() -> None:
    """Write the data files and print the calibration and the readings."""
    write(extract())
    for letter in PANELS:
        axes = frame(letter)
        print(
            f"Fig. 2({letter}): tick residuals x {residuals(axes.x):.4f} pt, "
            f"y {residuals(axes.y):.4f} pt; {1 / axes.x.resolution():.4f} pt/au, "
            f"{axes.y.resolution():.4g} per pt"
        )
        for method, curve in curves(letter).items():
            if not curve.dotted:
                x0 = x_of(axes.x, 520.0)
                print(f"  {method}: line width {curve.width(x0):.3f} pt at 520 au")
    for letter, results in readings().items():
        print(f"Fig. 2({letter}) {PANELS[letter].intensity} W/cm^2")
        for method, result in results.items():
            at, end = result.pulse_end, result.curve_end
            print(
                f"  {method:10s} 3T: {at.value:.4g}  end ({end.time:.2f} au): "
                f"{end.value:.4g}  end/3T - 1: {end.value / at.value - 1:+.4f}"
            )


if __name__ == "__main__":
    main()
