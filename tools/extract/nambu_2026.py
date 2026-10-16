# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Nambu et al. (2026): measured ionisation phase shifts of Ar, Kr and Xe.

**Figs. 6 and 7** (pdf pp. 7 and 8, library figures ``fig_119981293e1f`` and
``fig_63d440c06a97``): the phase shift of the probe, normalised to 1 atm (rad/atm),
against the pump's peak intensity (W/cm^2), in Ar, Kr and Xe ionised by 800 nm (Fig. 6)
and 400 nm (Fig. 7) pulses, linearly (LP, blue squares) and circularly (CP, red
circles) polarised, with error bars. Only the measured points are taken: not the PPT
and TDSE curves, nor Fig. 8 (267 nm).

Both figures are vector drawings, each placed on its page as one form XObject, and are
read from the drawing (data rule 8, TC-699) with :func:`_figure.painted`, which maps
the form's paths to the page:

* the panels: each plot area is a white-filled rectangle about 87.5 by 105.5 pt; the
  panel's letter, gas ('Ar', 'Kr', 'Xe') and wavelength ('800 nm', '400 nm') are the
  text objects inside it. The legend is a smaller white-filled rectangle inside it.
* the axes: the tick marks are black lines on the plot area's bottom (x) and left (y)
  edges, major ticks 1.05 pt long and minor ones 0.53 pt. The x majors are the
  labelled decades (10^14 W/cm^2 in Fig. 6; 10^13 and 10^14 in Fig. 7), the y majors
  10^-3 to 10^3 rad/atm, and the minor ticks the multiples 2 to 9 of the decades,
  counted outward from the majors. Fig. 6(c) draws no x major (the gamma = 1 line is
  drawn over the 10^14 tick); its minor ticks are read as the only run of consecutive
  multiples that lies on one log scale. Every tick is checked against the fitted log
  scale (:meth:`_figure.Axis.check`, within 0.01 pt).
* the markers: LP, the blue-filled squares (five-vertex paths, 1.17 pt wide; each
  also drawn as a stroked outline, not used); CP, the red-filled circles (13-vertex
  paths, 1.56 pt wide; each also drawn as a 26-vertex ring, not used). Those inside
  the legend are dropped. Their centres.
* the error bars: each is drawn as two vertical lines from the marker's centre, one
  down and one up, with a horizontal cap at each end, in the marker's colour (where
  one marker's bar ends at the next marker's centre, the two lines are one polyline,
  split here; where a cap crosses Fig. 6's gamma = 1 line, it is drawn as a filled
  outline). Each line is assigned to the marker whose centre its first vertex is
  nearest, unambiguously (:data:`BAR_TO_MARKER`, :data:`BAR_MARGIN`); every marker has
  one line each way and a cap at each end. The bars are symmetric about the value in
  rad, not on the log axis (checked, :data:`SYMMETRY`): they are +-sigma, and sigma is
  half the bar's length in rad. Bars shorter than the marker (hidden in print) are
  in the drawing and read like the others.

Every point's ``printed`` cell gives the panel, the marker's centre and the bar's two
ends in page points. The row is the point's number in its panel, in order of
intensity then phase (several markers can share an intensity).

**Check** (:func:`checks`): Fig. 6's vertical dash-dot lines mark where the Keldysh
parameter is 1 (its caption); read on the x axis, they must agree with
I = I_p c eps0 m_e omega^2 / e^2 at 800 nm, from the paper's ionisation energies
(p. 8: Ar 15.76, Kr 14.00, Xe 12.13 eV).

**Parameters**: the record's pump and probe pulse durations and spot sizes are checked
against Table 2 (``tbl_454c4abab809``), the probe that goes with each pump being Nie
et al. 2022's (the record's notes).
"""

import math
import tomllib
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
from _figure import Axis, Painted, painted, texts
from _library import DATA, Point, table, write

from optmat.constants import E_CHARGE, EPS0, M_E, C

DIRECTORY = "Nambu_2026_StrongFieldPhysics"
SOURCE = "nambu_2026"
BLUE, RED, BLACK, WHITE = (0, 0, 255), (255, 0, 0), (0, 0, 0), (255, 255, 255)
KELDYSH_GREY = (23, 23, 23)
Y_DECADES = tuple(range(-3, 4))
# A bar's first vertex lies within 0.007 pt of its marker's centre (the marker's
# outline and the bar are rounded separately), and two markers of a panel can be
# 0.0105 pt apart (Fig. 6(b), LP): a bar's nearest marker must be 0.003 pt nearer than
# the next.
BAR_TO_MARKER = 0.008
BAR_MARGIN = 0.003
# A tick or bar is vertical or horizontal to the drawing's precision (0.001 pt) and a
# little more: bars start up to 0.003 pt either side of the marker's x.
STRAIGHT = 0.005
# A stroke drawn as its filled outline is 0.526 pt wide; the outline's mid-line lies
# up to 0.006 pt from the bar's end.
CAP_OUTLINE = 0.6
CAP_TO_END = 0.01
# Major ticks are 1.05 pt long, minor ticks 0.53 pt.
MAJOR_TICK = 0.8
# The bars are +-sigma in rad, not symmetric on the log axis: a bar of half-length h
# decades would be asymmetric in rad by about h ln 10 of sigma (20 % for the longest,
# 0.09 decade). The halves' difference is allowed 1 % of sigma, or 0.005 pt of the
# axis (the marker's centre and the bar's ends are each rounded to 0.001-0.003 pt),
# which is what bars shorter than about 0.3 pt reach (up to 2.7 % of sigma).
SYMMETRY = 0.01
BAR_ROUNDING = 0.005
# How far Fig. 6's gamma = 1 lines may stray from I_p's intensity (the lines are drawn
# at the paper's own value of the constants and I_p; they agree to 0.08 %).
KELDYSH_AGREEMENT = 0.002
IONISATION_EV = {"Ar": 15.76, "Kr": 14.00, "Xe": 12.13}


@dataclass(frozen=True)
class Figure:
    """One figure: its page, its region on the page and its axes' labelled decades."""

    name: str
    page: int
    region: tuple[float, float, float, float]
    label: str
    x_decades: tuple[int, ...]


FIGURES = (
    Figure("fig6", 7, (320.0, 145.0, 545.0, 400.0), "800 nm", (14,)),
    Figure("fig7", 8, (58.0, 478.0, 280.0, 735.0), "400 nm", (13, 14)),
)
MARKERS = {"LP": (BLUE, 5), "CP": (RED, 13)}


@dataclass(frozen=True)
class Panel:
    """A panel of a figure: its gas, plot area, legend box and axes."""

    gas: str
    letter: str
    area: tuple[float, float, float, float]
    legend: tuple[float, float, float, float]
    x_axis: Axis
    y_axis: Axis


def _box(drawn: Painted) -> tuple[float, float, float, float]:
    """Return a rectangle's vertices' extent (l, b, r, t), without the stroke."""
    xs, ys = drawn.vertices[:, 0], drawn.vertices[:, 1]
    return float(xs.min()), float(ys.min()), float(xs.max()), float(ys.max())


def _inside(x: float, y: float, box: tuple[float, float, float, float]) -> bool:
    """Return whether a point lies inside a box (l, b, r, t)."""
    return box[0] <= x <= box[2] and box[1] <= y <= box[3]


def _segments(drawn: list[Painted], colour: tuple[int, int, int]) -> list[np.ndarray]:
    """Return the stroked lines of a colour, as 2x2 arrays, polylines split.

    A polyline is split into its segments: where one marker's bar ends at the centre
    of the next marker, the PDF joins the two bars into one three-vertex line.
    """
    return [
        d.vertices[k : k + 2]
        for d in drawn
        if d.stroke == colour and d.fill is None and len(d.vertices) <= 4
        for k in range(len(d.vertices) - 1)
    ]


def _ticks(
    lines: list[np.ndarray], area: tuple[float, float, float, float], axis: int
) -> tuple[list[float], list[float]]:
    """Return the major and minor tick positions on the area's bottom or left edge.

    ``axis`` 0 is x (vertical ticks on the bottom edge), 1 is y (horizontal ticks on
    the left edge).
    """
    edge = area[1] if axis == 0 else area[0]
    along, across = (0, 1) if axis == 0 else (1, 0)
    major, minor = [], []
    for line in lines:
        if abs(line[0, along] - line[1, along]) > STRAIGHT:
            continue
        ends = sorted(line[:, across])
        if abs(ends[0] - edge) > STRAIGHT:
            continue
        length = ends[1] - ends[0]
        if length > 2 * MAJOR_TICK:
            continue
        position = float(line[:, along].mean())
        lo, hi = (area[0], area[2]) if axis == 0 else (area[1], area[3])
        if not lo - STRAIGHT <= position <= hi + STRAIGHT:
            continue
        (major if length > MAJOR_TICK else minor).append(position)
    return sorted(set(np.round(major, 3))), sorted(set(np.round(minor, 3)))


def _log_axis(
    major: list[float], minor: list[float], decades: tuple[int, ...], where: str
) -> Axis:
    """Return a log axis from its major ticks (the decades) and its minor ticks.

    The minor ticks between two majors are 2 to 9 times the lower decade; those beyond
    the first major are counted down from it (9, 8, ... times the decade below), those
    beyond the last up from it (2, 3, ...). Minor ticks on a major are dropped.
    """
    if len(major) != len(decades):
        msg = f"{where}: {len(major)} major ticks, expected {len(decades)}"
        raise ValueError(msg)
    minor = [m for m in minor if min(abs(m - p) for p in major) > STRAIGHT]
    positions = list(major)
    values = [10.0**d for d in decades]
    below = sorted((m for m in minor if m < major[0]), reverse=True)
    above = sorted(m for m in minor if m > major[-1])
    for k, m in enumerate(below):
        positions.append(m)
        values.append((9 - k) * 10.0 ** (decades[0] - 1))
    for k, m in enumerate(above):
        positions.append(m)
        values.append((2 + k) * 10.0 ** decades[-1])
    for lower, upper, decade in zip(major, major[1:], decades, strict=False):
        between = sorted(m for m in minor if lower < m < upper)
        if len(between) != 8:
            msg = f"{where}: {len(between)} minor ticks in a decade"
            raise ValueError(msg)
        positions += between
        values += [k * 10.0**decade for k in range(2, 10)]
    if len(below) > 8 or len(above) > 8:
        msg = f"{where}: more minor ticks beyond the majors than a decade holds"
        raise ValueError(msg)
    order = np.argsort(positions)
    axis = Axis(
        tuple(float(positions[i]) for i in order),
        tuple(float(values[i]) for i in order),
        log=True,
        tolerance=0.01,
    )
    axis.check()
    return axis


def _without_major(minor: list[float], decade: int, where: str) -> Axis:
    """Return a log axis from minor ticks alone, the decade's major tick not drawn.

    The ticks are consecutive multiples 2, ..., 9 of 10^(decade - 1), then of
    10^decade (the decade itself not drawn). Each possible first multiple is tried;
    exactly one must put every tick on one log scale.
    """
    fits = []
    for first in range(2, 10):
        values = []
        k, power = first, decade - 1
        for _ in minor:
            values.append(k * 10.0**power)
            k += 1
            if k == 10:
                k, power = 2, power + 1
        axis = Axis(tuple(minor), tuple(values), log=True, tolerance=0.01)
        try:
            axis.check()
        except ValueError:
            continue
        fits.append(axis)
    if len(fits) != 1:
        msg = f"{where}: {len(fits)} readings of the minor ticks fit a log scale"
        raise ValueError(msg)
    return fits[0]


def _panels(figure: Figure, drawn: list[Painted]) -> list[Panel]:
    """Return a figure's three panels, with their axes calibrated from the ticks.

    A panel whose x major tick is not drawn (Fig. 6(c): the gamma = 1 line is drawn
    over the 10^14 tick) is calibrated from its minor ticks alone
    (:func:`_without_major`).
    """
    areas = [_box(d) for d in drawn if d.fill == WHITE and d.size > 80]
    legends = [
        _box(d)
        for d in drawn
        if d.fill == WHITE and d.stroke is None and 30 < d.size < 40
    ]
    words = texts(DIRECTORY, figure.page, figure.region)
    black = _segments(drawn, BLACK)
    panels: list[Panel] = []
    for area in sorted(set(areas), key=lambda a: (-a[3], a[0])):
        gases = [w for w, x, y in words if w in IONISATION_EV and _inside(x, y, area)]
        letters = [w for w, x, y in words if w.startswith("(") and _inside(x, y, area)]
        labels = [w for w, x, y in words if w == figure.label and _inside(x, y, area)]
        boxes = [b for b in legends if _inside(b[0], b[1], area)]
        if len(gases) != 1 or len(letters) != 1 or len(labels) != 1 or len(boxes) != 1:
            msg = f"{figure.name}: a panel's gas, letter, label or legend is not found"
            raise ValueError(msg)
        where = f"{figure.name}{letters[0]}"
        x_major, x_minor = _ticks(black, area, 0)
        y_major, y_minor = _ticks(black, area, 1)
        if not x_major and len(figure.x_decades) == 1:
            x_axis = _without_major(x_minor, figure.x_decades[0], f"{where} x")
        else:
            x_axis = _log_axis(x_major, x_minor, figure.x_decades, f"{where} x")
        panels.append(
            Panel(
                gases[0],
                letters[0],
                area,
                boxes[0],
                x_axis,
                _log_axis(y_major, y_minor, Y_DECADES, f"{where} y"),
            )
        )
    if [p.letter for p in panels] != ["(a)", "(b)", "(c)"]:
        msg = f"{figure.name}: panels {[p.letter for p in panels]}"
        raise ValueError(msg)
    return panels


@dataclass(frozen=True)
class Reading:
    """One measured point as read: marker centre and bar ends (pt), and values."""

    x: float
    y: float
    low: float
    high: float
    intensity: float
    phase: float
    sigma: float
    asymmetry: float
    hidden: bool


def _readings(
    panel: Panel, drawn: list[Painted], polarization: str, where: str
) -> list[Reading]:
    """Return a panel's markers of one polarisation, with their error bars."""
    colour, vertices = MARKERS[polarization]
    markers = [
        d
        for d in drawn
        if len(d.vertices) == vertices
        and d.fill == colour
        and d.stroke is None
        and d.size < 2
        and _inside(*d.centre, panel.area)
        and not _inside(*d.centre, panel.legend)
    ]
    centres = np.array([m.centre for m in markers])
    lines = [
        s
        for s in _segments(drawn, colour)
        if _inside(s[0, 0], s[0, 1], panel.area)
        and not _inside(s[0, 0], s[0, 1], panel.legend)
    ]
    bars = [s for s in lines if abs(s[0, 0] - s[1, 0]) < STRAIGHT]
    caps = [s for s in lines if abs(s[0, 1] - s[1, 1]) < STRAIGHT]
    # Where a bar's cap crosses Fig. 6's gamma = 1 line, the PDF draws it as a filled
    # rectangle, the stroke's outline (0.53 pt high): its mid-line is the cap.
    for d in drawn:
        if (
            len(d.vertices) == 5
            and d.fill == colour
            and d.stroke is None
            and d.top - d.bottom < CAP_OUTLINE
            and _inside(*d.centre, panel.area)
            and not _inside(*d.centre, panel.legend)
        ):
            left, bottom, right, top = _box(d)
            middle = (top + bottom) / 2
            caps.append(np.array([[left, middle], [right, middle]]))
    owned: dict[int, list[np.ndarray]] = defaultdict(list)
    for bar in bars:
        distance = np.sort(np.hypot(*(centres - bar[0]).T))
        nearest = int(np.argmin(np.hypot(*(centres - bar[0]).T)))
        if distance[0] > BAR_TO_MARKER or distance[1] - distance[0] < BAR_MARGIN:
            msg = f"{where}: the bar from {bar[0]} pt has no unambiguous marker"
            raise ValueError(msg)
        owned[nearest].append(bar)
    readings = []
    for index, marker in enumerate(markers):
        mine = owned[index]
        ups = [bar for bar in mine if bar[1, 1] >= bar[0, 1]]
        if len(mine) != 2 or len(ups) != 1:
            msg = f"{where}: the marker at {marker.centre} has {len(mine)} bar lines"
            raise ValueError(msg)
        x, y = marker.centre
        low, high = sorted(float(bar[1, 1]) for bar in mine)
        for end in (low, high):
            if not any(
                abs(cap[0, 1] - end) < CAP_TO_END
                and min(cap[:, 0]) < x < max(cap[:, 0])
                for cap in caps
            ):
                msg = f"{where}: no cap at the bar end ({x:.3f}, {end:.3f})"
                raise ValueError(msg)
        phase = float(panel.y_axis(y))
        lower, upper = float(panel.y_axis(low)), float(panel.y_axis(high))
        sigma = (upper - lower) / 2
        asymmetry = ((upper - phase) - (phase - lower)) / 2
        rounding = phase * BAR_ROUNDING * panel.y_axis.resolution()
        allowed = max(SYMMETRY * sigma, rounding)
        if abs(asymmetry) > allowed:
            msg = f"{where}: the bar at ({x:.3f}, {y:.3f}) is not +-sigma in rad"
            raise ValueError(msg)
        half_size = marker.size / 2
        readings.append(
            Reading(
                x,
                y,
                low,
                high,
                float(panel.x_axis(x)),
                phase,
                sigma,
                asymmetry / sigma if sigma else 0.0,
                high - y < half_size and y - low < half_size,
            )
        )
    if len(caps) != 2 * len(markers):
        msg = f"{where}: {len(caps)} caps for {len(markers)} markers"
        raise ValueError(msg)
    return sorted(readings, key=lambda r: (r.x, r.y))


def readings() -> dict[str, list[Reading]]:
    """Return every panel's readings, by dataset id (``fig6_Ar_LP``, ...)."""
    out = {}
    for figure in FIGURES:
        drawn = painted(DIRECTORY, figure.page, figure.region)
        for panel in _panels(figure, drawn):
            for polarization in MARKERS:
                where = f"{figure.name}{panel.letter} {polarization}"
                key = f"{figure.name}_{panel.gas}_{polarization}"
                out[key] = _readings(panel, drawn, polarization, where)
    return out


def _points(name: str, found: list[Reading], letter: str) -> list[Point]:
    """Return a dataset's points, numbered in order of intensity."""
    points = []
    for k, r in enumerate(found, start=1):
        printed = (
            f"{name} {letter}: marker at ({r.x:.3f}, {r.y:.3f}) pt; "
            f"bar from ({r.x:.3f}, {r.low:.3f}) to ({r.x:.3f}, {r.high:.3f}) pt"
        )
        points.append(
            Point(
                f"p{k:02d}",
                f"{r.intensity:.4e}",
                f"{r.phase:.4e}",
                f"{r.sigma:.3e}",
                printed,
            )
        )
    return points


def keldysh_lines() -> list[tuple[str, float, float]]:
    """Return Fig. 6's gamma = 1 lines: the gas, the line's intensity and I_p's."""
    figure = FIGURES[0]
    drawn = painted(DIRECTORY, figure.page, figure.region)
    omega = 2 * math.pi * C / 800e-9
    out = []
    for panel in _panels(figure, drawn):
        # The line's dashes, each a filled rectangle 0.526 pt wide (the stroke's
        # outline); pieces merged with a crossing marker or bar are left out.
        dashes = [
            _box(d)
            for d in drawn
            if d.fill == KELDYSH_GREY
            and len(d.vertices) == 5
            and _inside(*d.centre, panel.area)
        ]
        xs = [(box[0] + box[2]) / 2 for box in dashes]
        if not dashes or max(xs) - min(xs) > STRAIGHT:
            msg = f"fig6{panel.letter}: the gamma = 1 line is not where expected"
            raise ValueError(msg)
        line = float(np.mean(xs))
        # gamma^2 = I_p / (2 U_p), U_p = e^2 E^2 / (4 m omega^2), I = c eps0 E^2 / 2.
        expected = IONISATION_EV[panel.gas] * E_CHARGE * C * EPS0 * M_E * omega**2
        expected /= E_CHARGE**2 * 1e4  # W/m^2 to W/cm^2
        out.append((panel.gas, float(panel.x_axis(line)), expected))
    return out


def checks() -> None:
    """Raise unless Fig. 6's gamma = 1 lines sit at I_p's intensity."""
    for gas, read, expected in keldysh_lines():
        if abs(read / expected - 1) > KELDYSH_AGREEMENT:
            msg = f"Fig. 6 {gas}: gamma = 1 at {read:.3e}, I_p gives {expected:.3e}"
            raise ValueError(msg)


TABLE2_HEADER = ["", "PulseDuration(fs)", "SpotSize(µm)"]
BEAMS = {
    # wavelength_nm: (pump row, probe row of the beam that probes this pump)
    800: ("800nmpump", "400nmprobe"),
    400: ("400nmpump", "800nmprobe"),
}


def check_parameters() -> None:
    """Raise unless the record's pulse durations and spot sizes are Table 2's."""
    grid = table(DIRECTORY, "tbl_454c4abab809")
    if grid[0] != TABLE2_HEADER:
        msg = f"unexpected header of Table 2: {grid[0]!r}"
        raise ValueError(msg)
    rows = {cells[0]: (float(cells[1]), float(cells[2])) for cells in grid[1:]}
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        record = tomllib.load(handle)
    for entry in record["datasets"]:
        p = entry["parameters"]
        pump, probe = BEAMS[p["wavelength_nm"]]
        stated = (p["pulse_fs"], p["pump_waist_um"])
        stated_probe = (p["probe_pulse_fs"], p["probe_waist_um"])
        if stated != rows[pump] or stated_probe != rows[probe]:
            msg = f"{entry['id']}: pulse or spot size differs from Table 2"
            raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return Figs. 6 and 7's measured phase shifts, one file per gas and pump."""
    checks()
    check_parameters()
    found = readings()
    out = {}
    for figure in FIGURES:
        drawn = painted(DIRECTORY, figure.page, figure.region)
        letters = {p.gas: p.letter for p in _panels(figure, drawn)}
        for gas, letter in letters.items():
            for polarization in MARKERS:
                key = f"{figure.name}_{gas}_{polarization}"
                name = "Fig. 6" if figure.name == "fig6" else "Fig. 7"
                out[f"{gas}/{SOURCE}__{key}.csv"] = _points(name, found[key], letter)
    return out


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
