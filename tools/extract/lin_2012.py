# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Lin et al. (2012): O2+/Xe+ and N2+/Ar+ yield ratios at 800 and 2000 nm.

The paper prints the ratios only in **Fig. 2** (p. 2, library figure
``fig_31c497b1a2f7``): '(a)' O2+/Xe+ and '(b)' N2+/Ar+ against the 'Laser intensity
(TW / cm^2)', each at 800 nm (black filled squares) and 2000 nm (red filled circles),
without error bars. The text gives no number for them. Fig. 1 (library figure
``fig_47d4eab47394``) plots the four yields themselves, 'Ion signal [Counts / (shot
Pa)]', at 2000 nm (main panels) and 800 nm (insets); it is read here only to check
Fig. 2 (:func:`fig1_quotients`), not curated (the record's notes).

Both figures are vector drawings on the page (no form, no image) and are read from the
drawing (data rule 8, TC-699) with :func:`_figure.painted`. The labels are glyphs
drawn as filled paths, not text.

* the frames: each panel's frame and ticks are stroked black paths of move-to/line-to
  pairs (checked: every pair is a horizontal or vertical line), two per panel (one
  for the inset of Fig. 1(a)), the same box. The ticks are the pairs perpendicular to
  an edge with one end on it, under 3 pt long: on the bottom (x) and left (y) edges,
  mirrored on the top and right ones in Fig. 2 and on the top in Fig. 1's main panels
  (checked: the same positions). Long ticks (about
  twice the short ones) are the labelled decades.
* the axes are logarithmic: the decades are those the labels print on the page
  (:data:`PANELS`), the short ticks 2 to 9 times a decade, counted outward from the
  long ticks (as in ``nambu_2026``). Every tick lies on the fitted log scale within
  :data:`TICK_TOLERANCE` (:meth:`_figure.Axis.check`). In Fig. 2 each long tick of the
  y axis, and only it, has a label to its left, and the x axis's only label lies
  under its one long tick (10^2 TW/cm^2) (:func:`_check_labels`). Fig. 2's x axis runs
  from 30 to 400 TW/cm^2 (80.9 pt per decade); the y axis of (a) from 0.025 to 1 (73.5
  pt per decade), of (b) from 0.1 to 10 (58.9 pt per decade).
* the markers: filled, unstroked paths of a series' colour, vertex count (5: a
  square; 13: a circle) and size (:data:`PANELS`), inside the frame, outside the
  legend box (:data:`PANELS`, from the drawing: the legend's markers and its labels'
  glyphs) and, in Fig. 1's main panels, outside the inset's frame. Each is also drawn
  as a stroked outline (not used); Fig. 1's open circles are white-filled with a red
  outline. Their centres. Exactly one marker of each series lies in the legend box,
  and no data marker touches it.

Every point's ``printed`` cell gives the panel, the series and the marker's centre in
page points. The rows are numbered in order of intensity. The intensity is converted
from TW/cm^2 to W/cm^2.

**Checks** (:func:`checks`):

* Fig. 2(a) is the quotient of Fig. 1(a)'s O2+ and Xe+ yields: every Fig. 2 point has a
  Fig. 1 pair (an O2+ and a Xe+ marker at the same x) at its intensity, within
  :data:`SAME_INTENSITY`, and their quotient agrees with it within :data:`QUOTIENT`,
  save the one point listed in :data:`QUOTIENT_EXCEPTIONS`. For N2+/Ar+ only some of
  Fig. 2(b)'s intensities are Fig. 1(b)'s, and at most of those the quotients differ;
  that is reported (:func:`n2_ar_comparison`), not checked.
* the text (p. 2): 'the ratio of N2+/Ar+ keeps almost constant, to be around 1'; 'The
  two ratios [O2+/Xe+] almost coincide at about the lowest intensity used in the
  experiment' and 'for the same intensities above 5-6 x 10^13 W/cm2, the longer the
  wavelength, the less pronounced ionization suppression of O2'; and (p. 3) 'the ratio
  of 800 nm increases by about 4 times when the laser intensity increases by about 4
  times' (:func:`text_numbers`).
"""

import math
from dataclasses import dataclass

import numpy as np
from _figure import Axis, Painted, painted
from _library import Point, write

DIRECTORY = "Lin_2012_IonizationSuppressionDiatomic"
SOURCE = "lin_2012"
PAGE = 2
FIG1_REGION = (320.0, 420.0, 570.0, 775.0)
FIG2_REGION = (320.0, 100.0, 560.0, 270.0)
BLACK, RED, WHITE = (0, 0, 0), (255, 0, 0), (255, 255, 255)
TW_CM2 = 1e12  # W/cm^2 per TW/cm^2
# The drawing's coordinates are given to 0.001 pt: a tick lies on its edge and is
# straight to that. The ticks lie on the fitted log scales within 0.028 pt (Fig. 2),
# 0.037 pt (Fig. 1's main panels) and 0.009 pt (its insets): the plotting program
# places them no better.
STRAIGHT = 0.002
TICK_TOLERANCE = 0.04
# Ticks are under 3 pt long (the longest, Fig. 2's y majors, 2.27 pt); a long tick is
# over 1.5 times the shortest of its axis (they are twice as long).
TICK_LONGEST = 3.0
LONG_TICK = 1.5
# A label's glyphs: within 14 pt left of the frame (the axis title is 17 pt or more
# away), within 1.5 pt of its tick's height; the x label within 12 pt below the frame
# and within 8 pt of its tick (the superscript sits to the right).
LABEL_LEFT = 14.0
LABEL_LEVEL = 1.5
LABEL_BELOW = 12.0
LABEL_ACROSS = 8.0
# A molecule's and its atom's markers in Fig. 1 are at the same x within 0.003 pt.
SAME_X = 0.01
# Fig. 2's intensities are Fig. 1's: within 0.5 % (Fig. 1's x scale is 38-41 pt per
# decade; 0.5 % is 0.09 pt).
SAME_INTENSITY = 0.005
# Fig. 2(a)'s ratio is Fig. 1(a)'s quotient within 1 %: Fig. 1's y scale is 19.6 pt
# per decade, and its markers placed as its ticks, to about 0.03 pt or 0.4 % of a
# yield; the quotients agree within 0.7 % (0.4 % in the 800 nm inset).
QUOTIENT = 0.01
# The one exception: at 800 nm, 110.5 TW/cm^2, Fig. 2(a) has 0.156 and the inset's
# yields give 0.133; the neighbouring points agree within 0.5 %.
QUOTIENT_EXCEPTIONS = {("800", 110.5)}
EXCEPTION_WIDTH = 0.01  # relative, to recognise the exception's intensity


@dataclass(frozen=True)
class Series:
    """A series of markers: its name, fill colour, vertex count and size range (pt)."""

    name: str
    fill: tuple[int, int, int]
    vertices: int
    size: tuple[float, float]


@dataclass(frozen=True)
class PanelSpec:
    """A panel: its name, frame's left edge, decades, series and legend box."""

    name: str
    region: tuple[float, float, float, float]
    left: float
    x_decades: tuple[int, ...]
    y_decades: tuple[int, ...]
    series: tuple[Series, ...]
    legend: tuple[float, float, float, float] | None = None
    inset_left: float | None = None


SQUARE_800 = Series("800", BLACK, 5, (2.0, 2.3))
CIRCLE_2000 = Series("2000", RED, 13, (2.3, 2.5))
PANELS = {
    # Fig. 2: decades read on the page: x '10^2' (TW/cm^2); y (a) '0.1', '1'; (b) '0.1',
    # '1', '10'.
    "fig2a": PanelSpec(
        "Fig. 2(a)",
        FIG2_REGION,
        349.352,
        (2,),
        (-1, 0),
        (SQUARE_800, CIRCLE_2000),
        legend=(358.5, 226.3, 402.6, 241.0),
    ),
    "fig2b": PanelSpec(
        "Fig. 2(b)",
        FIG2_REGION,
        459.184,
        (2,),
        (-1, 0, 1),
        (SQUARE_800, CIRCLE_2000),
        legend=(470.9, 226.4, 515.0, 241.1),
    ),
    # Fig. 1: x '10^1', '10^2' (and '10^3' in the inset of (a)); y '10^-1' to '10^5'
    # ('10^0' to '10^5' in the inset of (b)).
    "fig1a": PanelSpec(
        "Fig. 1(a)",
        FIG1_REGION,
        366.086,
        (1, 2),
        (-1, 0, 1, 2, 3, 4, 5),
        (
            Series("Xe", BLACK, 5, (2.0, 2.3)),
            Series("O2", WHITE, 13, (2.3, 2.5)),
        ),
        legend=(375.0, 714.5, 394.5, 732.5),
        inset_left=413.148,
    ),
    "fig1b": PanelSpec(
        "Fig. 1(b)",
        FIG1_REGION,
        461.445,
        (1, 2),
        (-1, 0, 1, 2, 3, 4, 5),
        (
            Series("Ar", BLACK, 5, (2.3, 2.6)),
            Series("N2", WHITE, 13, (2.7, 2.9)),
        ),
        legend=(471.0, 714.5, 489.5, 732.5),
        inset_left=508.741,
    ),
    "fig1a_inset": PanelSpec(
        "Fig. 1(a) inset",
        FIG1_REGION,
        413.148,
        (1, 2, 3),
        (-1, 0, 1, 2, 3, 4, 5),
        (
            Series("Xe", BLACK, 5, (0.9, 1.1)),
            Series("O2", WHITE, 13, (1.0, 1.3)),
        ),
    ),
    "fig1b_inset": PanelSpec(
        "Fig. 1(b) inset",
        FIG1_REGION,
        508.741,
        (1, 2),
        (0, 1, 2, 3, 4, 5),
        (
            Series("Ar", BLACK, 5, (0.9, 1.1)),
            Series("N2", WHITE, 13, (1.0, 1.2)),
        ),
    ),
}
# The curated datasets: (panel, series) -> data file.
DATASETS = {
    ("fig2a", "800"): f"O2/{SOURCE}__fig2a_O2_Xe_800nm.csv",
    ("fig2a", "2000"): f"O2/{SOURCE}__fig2a_O2_Xe_2000nm.csv",
    ("fig2b", "800"): f"N2/{SOURCE}__fig2b_N2_Ar_800nm.csv",
    ("fig2b", "2000"): f"N2/{SOURCE}__fig2b_N2_Ar_2000nm.csv",
}


@dataclass(frozen=True)
class Frame:
    """A panel's frame: its box (l, b, r, t, pt) and its stroked line segments."""

    box: tuple[float, float, float, float]
    segments: tuple[np.ndarray, ...]


def _frame(drawn: list[Painted], left: float, where: str) -> Frame:
    """Return the frame whose left edge is at ``left`` (pt), with its tick segments.

    The frame is drawn as one or two stroked black paths of move-to/line-to pairs;
    each pair must be a horizontal or vertical line.
    """
    paths = [
        d
        for d in drawn
        if d.stroke == BLACK
        and d.fill is None
        and len(d.vertices) >= 20
        and abs(float(d.vertices[:, 0].min()) - left) < STRAIGHT
    ]
    if not paths:
        msg = f"{where}: no frame"
        raise ValueError(msg)
    segments = []
    for path in paths:
        if len(path.vertices) % 2:
            msg = f"{where}: a frame path of {len(path.vertices)} vertices"
            raise ValueError(msg)
        for pair in path.vertices.reshape(-1, 2, 2):
            dx, dy = np.abs(pair[1] - pair[0])
            if min(dx, dy) > STRAIGHT:
                msg = f"{where}: a frame segment is not horizontal or vertical"
                raise ValueError(msg)
            segments.append(pair)
    every = np.concatenate([p.vertices for p in paths])
    box = (
        float(every[:, 0].min()),
        float(every[:, 1].min()),
        float(every[:, 0].max()),
        float(every[:, 1].max()),
    )
    return Frame(box, tuple(segments))


def _edge_ticks(frame: Frame, axis: int, far: bool) -> list[tuple[float, float]]:
    """Return the ticks on an edge as (position, length), sorted by position.

    ``axis`` 0: the bottom (or, ``far``, top) edge's vertical ticks; 1: the left (or
    right) edge's horizontal ticks.
    """
    along, across = (0, 1) if axis == 0 else (1, 0)
    edge = frame.box[across + 2] if far else frame.box[across]
    ticks = []
    for segment in frame.segments:
        if abs(segment[0, along] - segment[1, along]) > STRAIGHT:
            continue
        ends = sorted(segment[:, across], key=lambda e: abs(e - edge))
        length = abs(ends[1] - ends[0])
        if abs(ends[0] - edge) > STRAIGHT or not 0 < length < TICK_LONGEST:
            continue
        ticks.append((round(float(segment[0, along]), 3), round(float(length), 3)))
    return sorted(set(ticks))


def _ticks(frame: Frame, axis: int, where: str) -> tuple[list[float], list[float]]:
    """Return an axis's long and short tick positions.

    Where the opposite edge carries ticks too (Fig. 2's four edges, Fig. 1's main
    panels' top edges), they must be at the same positions.
    """
    near, far = _edge_ticks(frame, axis, False), _edge_ticks(frame, axis, True)
    if far and [p for p, _ in near] != [p for p, _ in far]:
        msg = f"{where}: the ticks of opposite edges differ"
        raise ValueError(msg)
    shortest = min(length for _, length in near)
    long = sorted({p for p, length in near if length > LONG_TICK * shortest})
    short = sorted(
        {p for p, length in near if length <= LONG_TICK * shortest} - set(long)
    )
    return long, short


def _log_axis(
    major: list[float], minor: list[float], decades: tuple[int, ...], where: str
) -> Axis:
    """Return a log axis from its long ticks (the decades) and its short ticks.

    The short ticks between two long ones are 2 to 9 times the lower decade; those
    below the first long tick are counted down from it (9, 8, ... times the decade
    below), those above the last up from it (2, 3, ...).
    """
    if len(major) != len(decades):
        msg = f"{where}: {len(major)} long ticks, expected {len(decades)}"
        raise ValueError(msg)
    below = sorted((m for m in minor if m < major[0]), reverse=True)
    above = sorted(m for m in minor if m > major[-1])
    if len(below) > 8 or len(above) > 8:
        msg = f"{where}: more short ticks beyond the long ones than a decade holds"
        raise ValueError(msg)
    positions = list(major)
    values = [10.0**d for d in decades]
    for k, m in enumerate(below):
        positions.append(m)
        values.append((9 - k) * 10.0 ** (decades[0] - 1))
    for k, m in enumerate(above):
        positions.append(m)
        values.append((2 + k) * 10.0 ** decades[-1])
    for lower, upper, decade in zip(major, major[1:], decades, strict=False):
        between = sorted(m for m in minor if lower < m < upper)
        if len(between) != 8:
            msg = f"{where}: {len(between)} short ticks in a decade"
            raise ValueError(msg)
        positions += between
        values += [k * 10.0**decade for k in range(2, 10)]
    order = np.argsort(positions)
    axis = Axis(
        tuple(float(positions[i]) for i in order),
        tuple(float(values[i]) for i in order),
        log=True,
        tolerance=TICK_TOLERANCE,
    )
    axis.check()
    return axis


def _inside(x: float, y: float, box: tuple[float, float, float, float]) -> bool:
    """Return whether a point lies inside a box (l, b, r, t)."""
    return box[0] <= x <= box[2] and box[1] <= y <= box[3]


def _touches(d: Painted, box: tuple[float, float, float, float]) -> bool:
    """Return whether a path's bounding box overlaps a box (l, b, r, t)."""
    return (
        d.right >= box[0]
        and d.left <= box[2]
        and d.top >= box[1]
        and d.bottom <= box[3]
    )


def _is_marker(d: Painted, series: Series) -> bool:
    """Return whether a path is a marker of a series (its fill, shape and size)."""
    return (
        d.fill == series.fill
        and d.stroke is None
        and len(d.vertices) == series.vertices
        and series.size[0] <= d.size <= series.size[1]
    )


def _check_labels(
    drawn: list[Painted],
    frame: Frame,
    x_long: list[float],
    y_long: list[float],
    where: str,
) -> None:
    """Raise unless the axes' label glyphs sit at the long ticks, and only there.

    The y labels: glyphs (black-filled paths) whose centre lies within
    :data:`LABEL_LEFT` pt left of the frame, grouped by height; each group within
    :data:`LABEL_LEVEL` pt of a long tick, one each. The x labels: glyphs within
    :data:`LABEL_BELOW` pt below the frame, all within :data:`LABEL_ACROSS` pt of the
    one long tick.
    """
    left, bottom, _, top = frame.box
    glyphs = [d for d in drawn if d.fill == BLACK and d.stroke is None]
    beside = sorted(
        d.centre[1]
        for d in glyphs
        if left - LABEL_LEFT < d.centre[0] < left and bottom - 3 < d.centre[1] < top + 3
    )
    groups: list[list[float]] = []
    for y in beside:
        if groups and y - groups[-1][-1] < 2 * LABEL_LEVEL:
            groups[-1].append(y)
        else:
            groups.append([y])
    levels = [float(np.mean(g)) for g in groups]
    if len(levels) != len(y_long) or any(
        abs(level - tick) > LABEL_LEVEL
        for level, tick in zip(levels, y_long, strict=True)
    ):
        msg = f"{where}: the y labels at {levels} are not at the long ticks {y_long}"
        raise ValueError(msg)
    under = [
        d.centre[0]
        for d in glyphs
        if bottom - LABEL_BELOW < d.centre[1] < bottom
        and left < d.centre[0] < frame.box[2]
    ]
    if (
        len(x_long) != 1
        or not under
        or any(abs(x - x_long[0]) > LABEL_ACROSS for x in under)
    ):
        msg = f"{where}: the x labels at {under} are not at the long tick {x_long}"
        raise ValueError(msg)


@dataclass(frozen=True)
class Reading:
    """One marker as read: its centre (pt) and its values."""

    x: float
    y: float
    intensity: float
    value: float


def _panel(key: str) -> dict[str, list[Reading]]:
    """Return a panel's markers by series, read on its calibrated axes."""
    spec = PANELS[key]
    drawn = painted(DIRECTORY, PAGE, spec.region)
    frame = _frame(drawn, spec.left, spec.name)
    x_long, x_short = _ticks(frame, 0, f"{spec.name} x")
    y_long, y_short = _ticks(frame, 1, f"{spec.name} y")
    x_axis = _log_axis(x_long, x_short, spec.x_decades, f"{spec.name} x")
    y_axis = _log_axis(y_long, y_short, spec.y_decades, f"{spec.name} y")
    if key.startswith("fig2"):
        _check_labels(drawn, frame, x_long, y_long, spec.name)
    inset = None
    if spec.inset_left is not None:
        inset = _frame(drawn, spec.inset_left, f"{spec.name}'s inset").box
    out = {}
    for series in spec.series:
        candidates = [
            d
            for d in drawn
            if _is_marker(d, series)
            and _inside(*d.centre, frame.box)
            and not (inset is not None and _inside(*d.centre, inset))
        ]
        if spec.legend is not None:
            in_legend = [d for d in candidates if _inside(*d.centre, spec.legend)]
            if len(in_legend) != 1:
                msg = f"{spec.name} {series.name}: {len(in_legend)} legend markers"
                raise ValueError(msg)
            candidates = [d for d in candidates if d is not in_legend[0]]
            if any(_touches(d, spec.legend) for d in candidates):
                msg = f"{spec.name} {series.name}: a data marker touches the legend"
                raise ValueError(msg)
        centres = sorted(
            {(round(d.centre[0], 3), round(d.centre[1], 3)) for d in candidates}
        )
        if len(centres) != len(candidates):
            msg = f"{spec.name} {series.name}: two markers share a centre"
            raise ValueError(msg)
        out[series.name] = [
            Reading(x, y, float(x_axis(x)), float(y_axis(y))) for x, y in centres
        ]
    return out


def readings() -> dict[tuple[str, str], list[Reading]]:
    """Return Fig. 2's markers, by (panel, series), in order of intensity."""
    out = {}
    for key in ("fig2a", "fig2b"):
        for name, found in _panel(key).items():
            out[key, name] = found
    return out


def _pairs(key: str, molecule: str, atom: str) -> list[tuple[float, float]]:
    """Return a Fig. 1 panel's (intensity, molecule/atom) at the x both are drawn."""
    found = _panel(key)
    out = []
    for m in found[molecule]:
        partners = [a for a in found[atom] if abs(a.x - m.x) < SAME_X]
        if len(partners) == 1:
            out.append((m.intensity, m.value / partners[0].value))
    return out


def fig1_quotients() -> list[tuple[str, float, float, float]]:
    """Return each Fig. 2(a) point with Fig. 1(a)'s quotient at its intensity.

    As (wavelength, intensity in TW/cm^2, Fig. 2's ratio, Fig. 1's O2+/Xe+); raises
    if a Fig. 2 point has no Fig. 1 pair within :data:`SAME_INTENSITY`.
    """
    found = readings()
    out = []
    for wavelength, key in (("2000", "fig1a"), ("800", "fig1a_inset")):
        pairs = _pairs(key, "O2", "Xe")
        intensities = np.array([i for i, _ in pairs])
        for r in found["fig2a", wavelength]:
            nearest = int(np.argmin(np.abs(np.log(intensities / r.intensity))))
            if abs(intensities[nearest] / r.intensity - 1) > SAME_INTENSITY:
                msg = f"Fig. 2(a) {wavelength} nm at {r.intensity:.1f}: no Fig. 1 pair"
                raise ValueError(msg)
            out.append((wavelength, r.intensity, r.value, pairs[nearest][1]))
    return out


def n2_ar_comparison() -> dict[str, dict[str, float]]:
    """Return, by wavelength, Fig. 1(b)'s and Fig. 2(b)'s N2+/Ar+ statistics.

    Fig. 1(b)'s quotients (where both markers are drawn at one x) and Fig. 2(b)'s
    points: their numbers, minima, medians and maxima; how many of Fig. 2(b)'s
    intensities are Fig. 1(b)'s (within :data:`SAME_INTENSITY`), how many of those
    agree with Fig. 1(b)'s quotient within :data:`QUOTIENT`, the highest intensity
    (TW/cm^2) at which one does not, and the largest relative difference.
    """
    found = readings()
    out = {}
    for wavelength, key in (("2000", "fig1b"), ("800", "fig1b_inset")):
        pairs = _pairs(key, "N2", "Ar")
        one = np.array([q for _, q in pairs])
        intensities = np.array([i for i, _ in pairs])
        two = found["fig2b", wavelength]
        values = np.array([r.value for r in two])
        shared = []
        for r in two:
            nearest = int(np.argmin(np.abs(np.log(intensities / r.intensity))))
            if abs(intensities[nearest] / r.intensity - 1) <= SAME_INTENSITY:
                shared.append((r.intensity, pairs[nearest][1] / r.value - 1))
        differing = [i for i, d in shared if abs(d) > QUOTIENT]
        out[wavelength] = {
            "fig1_n": len(one),
            "fig1_min": float(one.min()),
            "fig1_median": float(np.median(one)),
            "fig1_max": float(one.max()),
            "fig2_n": len(values),
            "fig2_min": float(values.min()),
            "fig2_median": float(np.median(values)),
            "fig2_max": float(values.max()),
            "fig2_at_fig1_intensity": len(shared),
            "agreeing": sum(abs(d) <= QUOTIENT for _, d in shared),
            "last_differing_TW_cm2": max(differing, default=float("nan")),
            "largest_difference": max(abs(d) for _, d in shared),
        }
    return out


def _power_law(found: list[Reading]) -> tuple[float, float]:
    """Return the exponent and prefactor of a power law fitted in log-log."""
    slope, intercept = np.polyfit(
        np.log10([r.intensity for r in found]), np.log10([r.value for r in found]), 1
    )
    return float(slope), float(10**intercept)


def text_numbers() -> dict[str, float]:
    """Return the numbers the text's statements are checked against.

    ``n2_ar_*``: the N2+/Ar+ points' range and median (both wavelengths). ``coincide``:
    the 2000 nm O2+/Xe+ over the 800 nm one at the lowest 800 nm intensity (log-log
    interpolated); ``above_60``: the smallest of that quotient at the 800 nm
    intensities above 60 TW/cm^2 that the 2000 nm points reach. ``slope_800``: the
    800 nm O2+/Xe+ power-law exponent; ``four_times``: the rise of that power law for
    four times the intensity; ``slope_2000``: the 2000 nm exponent.
    """
    found = readings()
    n2 = [r.value for w in ("800", "2000") for r in found["fig2b", w]]
    low, high = found["fig2a", "800"], found["fig2a", "2000"]
    xs = np.log([r.intensity for r in high])
    ys = np.log([r.value for r in high])

    def interpolated(intensity: float) -> float:
        return float(np.exp(np.interp(math.log(intensity), xs, ys)))

    above = [
        interpolated(r.intensity) / r.value
        for r in low
        if 60 <= r.intensity <= high[-1].intensity
    ]
    slope, _ = _power_law(low)
    slope_2000, _ = _power_law(high)
    return {
        "n2_ar_min": min(n2),
        "n2_ar_median": float(np.median(n2)),
        "n2_ar_max": max(n2),
        "coincide": interpolated(low[0].intensity) / low[0].value,
        "above_60": min(above),
        "slope_800": slope,
        "four_times": 4.0**slope,
        "slope_2000": slope_2000,
    }


def checks() -> None:
    """Raise unless Fig. 2(a) is Fig. 1(a)'s quotient and the text's statements hold.

    The text's bounds are loose (its words are 'around 1', 'almost coincide', 'about 4
    times'): N2+/Ar+ within 0.8 to 1.3 with a median within 10 % of 1; O2+/Xe+ at the
    lowest intensity within 20 % between the wavelengths, and 2000 nm above 800 nm at
    every intensity above 60 TW/cm^2; the 800 nm ratio rising by 3 to 5 times for 4
    times the intensity.
    """
    for wavelength, intensity, ratio, quotient in fig1_quotients():
        excepted = any(
            wavelength == w and abs(intensity / i - 1) < EXCEPTION_WIDTH
            for w, i in QUOTIENT_EXCEPTIONS
        )
        if not excepted and abs(quotient / ratio - 1) > QUOTIENT:
            msg = (
                f"Fig. 2(a) {wavelength} nm at {intensity:.1f} TW/cm^2: {ratio:.4f}, "
                f"Fig. 1(a)'s quotient {quotient:.4f}"
            )
            raise ValueError(msg)
    n = text_numbers()
    if (
        not (n["n2_ar_min"] > 0.8 and n["n2_ar_max"] < 1.3)
        or abs(n["n2_ar_median"] - 1) > 0.1
    ):
        msg = f"N2+/Ar+ is not 'around 1': {n}"
        raise ValueError(msg)
    if abs(n["coincide"] - 1) > 0.2 or n["above_60"] <= 1:
        msg = f"the O2+/Xe+ ratios of the two wavelengths do not behave as stated: {n}"
        raise ValueError(msg)
    if not 3 < n["four_times"] < 5:
        msg = f"the 800 nm O2+/Xe+ ratio does not rise 'about 4 times': {n}"
        raise ValueError(msg)


def _points(name: str, series: str, found: list[Reading]) -> list[Point]:
    """Return a dataset's points, numbered in order of intensity."""
    return [
        Point(
            f"p{k:02d}",
            f"{r.intensity * TW_CM2:.4e}",
            f"{r.value:.4g}",
            "",
            f"{name} {series} nm: marker at ({r.x:.3f}, {r.y:.3f}) pt",
        )
        for k, r in enumerate(found, start=1)
    ]


def extract() -> dict[str, list[Point]]:
    """Return Fig. 2's ratios, one file per panel and wavelength."""
    checks()
    found = readings()
    return {
        relative: _points(PANELS[key].name, series, found[key, series])
        for (key, series), relative in DATASETS.items()
    }


def main() -> None:
    """Write the data files and print the checks' numbers."""
    write(extract())
    for wavelength, intensity, ratio, quotient in fig1_quotients():
        print(f"{wavelength} nm {intensity:6.1f} TW/cm2: {ratio:.4f} {quotient:.4f}")
    print(n2_ar_comparison())
    print(text_numbers())


if __name__ == "__main__":
    main()
