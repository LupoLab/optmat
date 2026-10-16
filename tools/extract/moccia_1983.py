# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Moccia, Rahman and Rizzo (1983): two-photon ionisation of Ne and Ar.

The paper prints its total two-photon ionisation cross sections (approximate RPA sum
over intermediate states, static-exchange continuum, length and velocity forms) and
their circular/linear ratios only in figures; no table of them, no supplement (the
differential coefficients are 'available upon request', p. 11/2746). The PDF is a scan:
every page is one embedded 1-bit image, 2076 x 2945 pixels (about 250 dpi), behind an
IOPscience cover page (PDF p. 1). The figures are read from those native pixels (data
rule 8, TC-699). Journal page = PDF page + 2735.

* **Fig. 1** (p. 7/2742, ``fig_b5118ea4e934``): the total cross section of Ne for
  linearly polarised light, y in units of 10^-52 cm^4 s (the axis is labelled 10^52
  sigma), 0 to 15: velocity (full curve) and length (chain curve).
* **Fig. 2** (p. 8/2743, ``fig_8d179928b5fc``): the same for Ar, 10^-51 cm^4 s, 0 to 36.
* **Fig. 4** (p. 10/2745, ``fig_31f952447a04``): sigma^C / sigma^L of Ne, 0 to 2,
  velocity (full) and length (chain); a full horizontal line marks 1.5.
* **Fig. 5** (p. 10/2745, ``fig_1d8254490151``): the same for Ar.

The abscissa is the photon energy in atomic units (hartree). Its scale is enlarged
after the second resonance energy (p. 8/2743): two linear scales, joined at a thick bar
on the axis. Vertical broken lines mark the ionisation thresholds and the resonance
energies used (Table 1's experimental values).

How the images are read:

* the axes: the tick marks below the x axis and left of the y axis (:data:`FIGURES`),
  each found in the image within ``TICK_FOUND`` pixels of the position written there
  (:func:`_check_ticks`); each x segment's ticks, and the y ticks with the axis line as
  0, lie on one linear scale within ``TICK_TOLERANCE`` pixels (:class:`Calibration`).
  The axes are drawn square to within 3.5 pixels; no rotation is applied.
* the curves: each stretch (:data:`SEGMENTS`) is traced column by column from a seed
  pixel on it (:func:`trace`), and its line style checked (:func:`style`): the full
  curve has no gap, the chain curve has gaps. That tells the two forms apart.
* the values: at each sample (every ``STEP`` hartree within a stretch), the centre of
  the curve's run at that photon energy's column where the run is clear: no wider than
  a line at its slope, no other dark pixel within ``CLEAR`` pixels, and so for the
  line's runs ``WINDOW`` columns either side, no steeper than ``MAX_SLOPE`` rows per
  column; across a gap of the chain, the linear interpolation between the nearest clear
  columns (:func:`_read`). The samples left out (another curve, a resonance line or a
  tick nearby, or the line too steep) are listed in :data:`LEFT_OUT`; where the two
  curves of a figure are read within ``LINE`` rows of each other the reading cannot tell
  whose line it is, and the sample is left out of both (:data:`TOGETHER`). The script
  fails if either list changes. The ``printed`` cell of each value gives its pixel
  (column, row).

**Checks** (:func:`checks`, run by :func:`main` and :func:`extract`; it fails on a
miss): the ticks as above; the drawn R = 1.5 line of Figs. 4 and 5 against the y scale
(``RATIO_LINE``). Reported: the pixel sizes; the vertical broken lines found, against
the thresholds (half the ionisation potentials of p. 8/2743) and Table 1's experimental
energies on each figure's scale; the length/velocity ratios; the comparisons with
McKenna and van der Hart 2004 (TC-856), Pan, Gao and Starace 1990 (pan_1990b) and
Starace and Jiang 1987 (TC-878) (:func:`comparisons`).
"""

import itertools
import json
from dataclasses import dataclass

import numpy as np
from _figure import Axis, image
from _library import LIBRARY, Point, write

from optmat.datasets import DATA_DIRECTORY, Dataset, load_source

DIRECTORY = "Moccia_1983_TwoPhotonIonisation"
SOURCE = "moccia_1983"
DARK = 128
HARTREE_EV = 27.211386245988


@dataclass(frozen=True)
class Figure:
    """One plot on a scanned page: its axes, ticks and break, in image pixels.

    ``axis_row`` is the x axis line's row and ``axis_column`` the y axis line's column.
    ``x_low`` and ``x_high`` are the tick marks (column, hartree) of the two x scales,
    left and right of the thick bar at ``bar``; ``y_ticks`` (row, printed value) are
    those of the y axis, the axis line counting as 0. ``unit`` multiplies the printed
    y value. ``material`` and ``threshold`` (the ionisation potential, hartree) name the
    gas. ``ratio_line`` is two pixels (column, row) on the full horizontal line drawn at
    R = 1.5 in Figs. 4 and 5.
    """

    page: int
    table: str
    unit: float
    material: str
    threshold: float
    axis_row: int
    axis_column: int
    top: int
    bar: float
    x_low: tuple[tuple[float, float], ...]
    x_high: tuple[tuple[float, float], ...]
    y_ticks: tuple[tuple[float, float], ...]
    ratio_line: tuple[tuple[int, int], tuple[int, int]] | None = None

    def line_row(self, column: float) -> float:
        """Return the row of the drawn R = 1.5 line at a column (Figs. 4, 5)."""
        if self.ratio_line is None:
            return -1e9
        (c0, r0), (c1, r1) = self.ratio_line
        return r0 + (column - c0) * (r1 - r0) / (c1 - c0)


# The x ticks are probed ``X_DEPTH`` rows below the axis line, the y ticks ``Y_DEPTH``
# columns left of it.
X_DEPTH = 10
Y_DEPTH = 6
TICK_FOUND = 1.5
TICK_TOLERANCE = 6.0
# The ionisation potentials 'employed in the calculations' (p. 8/2743), hartree.
IP_NE = 0.7925
IP_AR = 0.5792

FIGURES = {
    "fig1": Figure(
        7,
        "fig_b5118ea4e934",
        1e-52,
        "Ne",
        IP_NE,
        2378,
        501,
        1600,
        1166.5,
        (
            (521.0, 0.40),
            (621.5, 0.45),
            (719.5, 0.50),
            (820.5, 0.55),
            (917.0, 0.60),
            (1017.5, 0.65),
            (1106.0, 0.70),
        ),
        ((1280.0, 0.74), (1447.0, 0.76), (1614.5, 0.78), (1785.5, 0.80)),
        (
            (1698.5, 15.0),
            (1833.5, 12.0),
            (1971.0, 9.0),
            (2107.5, 6.0),
            (2242.5, 3.0),
        ),
    ),
    "fig2": Figure(
        8,
        "fig_8d179928b5fc",
        1e-51,
        "Ar",
        IP_AR,
        1049,
        483,
        250,
        1189.5,
        ((556.0, 0.30), (700.5, 0.35), (842.0, 0.40), (986.5, 0.45), (1129.0, 0.50)),
        ((1227.5, 0.52), (1409.0, 0.54), (1588.5, 0.56), (1770.5, 0.58)),
        ((323.0, 36.0), (505.5, 27.0), (688.0, 18.0), (868.0, 9.0)),
    ),
    "fig4": Figure(
        10,
        "fig_31f952447a04",
        1.0,
        "Ne",
        IP_NE,
        1414,
        502,
        650,
        1167.5,
        (
            (522.5, 0.40),
            (620.0, 0.45),
            (719.0, 0.50),
            (819.5, 0.55),
            (916.5, 0.60),
            (1017.5, 0.65),
            (1107.5, 0.70),
        ),
        ((1281.5, 0.74), (1449.5, 0.76), (1617.5, 0.78), (1787.5, 0.80)),
        ((732.5, 2.0), (902.5, 1.5), (1074.5, 1.0), (1242.0, 0.5)),
        ((625, 903), (1682, 908)),
    ),
    "fig5": Figure(
        10,
        "fig_1d8254490151",
        1.0,
        "Ar",
        IP_AR,
        2471,
        491,
        1700,
        1194.0,
        ((562.0, 0.30), (706.5, 0.35), (847.0, 0.40), (992.0, 0.45), (1135.0, 0.50)),
        ((1233.0, 0.52), (1414.5, 0.54), (1593.0, 0.56), (1774.0, 0.58)),
        ((1748.0, 2.0), (1931.0, 1.5), (2113.0, 1.0), (2292.0, 0.5)),
        ((634, 1929), (1749, 1929)),
    ),
}

Run = tuple[int, int]


def _spans(mask: np.ndarray) -> list[Run]:
    """Return the first and last index of each run of True."""
    index = np.flatnonzero(mask)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [(int(group[0]), int(group[-1])) for group in groups]


def _centre(run: Run) -> float:
    """Return a run's centre."""
    return (run[0] + run[1]) / 2


def _axis_row(dark: np.ndarray, figure: Figure, column: int) -> float:
    """Return the centre of the x axis line at a column."""
    low = figure.axis_row - 4
    runs = _spans(dark[low : figure.axis_row + 5, column])
    if not runs:
        msg = f"p. {figure.page}: no x axis at column {column}"
        raise ValueError(msg)
    return low + _centre(
        min(runs, key=lambda r: abs(low + _centre(r) - figure.axis_row))
    )


class Calibration:
    """Map a pixel (column, row) to (hartree, value) through a figure's ticks."""

    def __init__(self, figure: Figure, dark: np.ndarray) -> None:
        self.figure = figure
        zero = _axis_row(dark, figure, figure.axis_column + 30)
        y = (*figure.y_ticks, (zero, 0.0))

        def scale(ticks: tuple[tuple[float, float], ...]) -> Axis:
            positions = tuple(p for p, _ in ticks)
            values = tuple(v for _, v in ticks)
            return Axis(positions, values, tolerance=TICK_TOLERANCE)

        self.axes = {
            "low": scale(figure.x_low),
            "high": scale(figure.x_high),
            "y": scale(y),
        }
        for axis in self.axes.values():
            axis.check()

    def energy(self, column: float) -> float:
        """Return the photon energy (hartree) at a column, on its side of the bar."""
        side = "low" if column < self.figure.bar else "high"
        return float(self.axes[side](column))

    def column(self, energy: float) -> float:
        """Return the (fractional) column of a photon energy (hartree)."""
        for side in ("low", "high"):
            axis = self.axes[side]
            a = float(axis(0.0))
            column = (energy - a) / (float(axis(1.0)) - a)
            if (column < self.figure.bar) == (side == "low"):
                return column
        msg = f"{energy} hartree falls in the bar"
        raise ValueError(msg)

    def value(self, row: float) -> float:
        """Return the y value, in the axis's printed unit, at a row."""
        return float(self.axes["y"](row))

    def residual(self) -> float:
        """Return the worst tick residual of the three scales, pixels."""
        worst = 0.0
        for axis in self.axes.values():
            a = float(axis(0.0))
            slope = float(axis(1.0)) - a
            for position, value in zip(axis.positions, axis.values, strict=True):
                worst = max(worst, abs((value - a) / slope - position))
        return worst


def _check_ticks(dark: np.ndarray, figure: Figure) -> None:
    """Raise unless every listed tick is found in the image."""
    below = dark[
        figure.axis_row + X_DEPTH, figure.axis_column : figure.axis_column + 1400
    ]
    found = [figure.axis_column + _centre(run) for run in _spans(below)]
    for position, value in (*figure.x_low, *figure.x_high, (figure.bar, 0.0)):
        if not any(abs(f - position) <= TICK_FOUND for f in found):
            msg = f"p. {figure.page}: no x tick {value} at {position}"
            raise ValueError(msg)
    left = dark[figure.top - 100 : figure.axis_row, figure.axis_column - Y_DEPTH]
    found = [figure.top - 100 + _centre(run) for run in _spans(left)]
    for position, value in figure.y_ticks:
        if not any(abs(f - position) <= TICK_FOUND for f in found):
            msg = f"p. {figure.page}: no y tick {value} at {position}"
            raise ValueError(msg)


@dataclass(frozen=True)
class Segment:
    """A stretch of one curve: data file, line style, a seed pixel, span (hartree)."""

    file: str
    figure: str
    style: str
    seed: tuple[int, int]
    span: tuple[float, float]


def _file(name: str, material: str) -> str:
    """Return a dataset's data file."""
    return f"{material}/{SOURCE}__{name}_{material}.csv"


def _segments(
    figure: str, name: str, material: str, rows: tuple[tuple, ...]
) -> tuple[Segment, ...]:
    """Return a curve's stretches: (style, seed, span) each."""
    return tuple(
        Segment(_file(name, material), figure, style, seed, span)
        for style, seed, span in rows
    )


SEGMENTS: tuple[Segment, ...] = (
    *_segments(
        "fig1",
        "fig1_sigma2_V",
        "Ne",
        (
            ("full", (625, 1940), (0.400, 0.615)),
            ("full", (1060, 2147), (0.625, 0.720)),
            ("full", (1307, 2185), (0.739, 0.758)),
        ),
    ),
    *_segments(
        "fig1",
        "fig1_sigma2_L",
        "Ne",
        (
            ("chain", (625, 2258), (0.400, 0.580)),
            ("chain", (1060, 2289), (0.625, 0.720)),
            ("chain", (1307, 2283), (0.739, 0.758)),
        ),
    ),
    *_segments(
        "fig2",
        "fig2_sigma2_V",
        "Ar",
        (
            ("full", (628, 878), (0.290, 0.430)),
            ("full", (1009, 867), (0.440, 0.510)),
            ("full", (1363, 910), (0.525, 0.545)),
            ("full", (1576, 936), (0.553, 0.563)),
        ),
    ),
    *_segments(
        "fig2",
        "fig2_sigma2_L",
        "Ar",
        (
            ("chain", (628, 822), (0.290, 0.355)),
            ("chain", (828, 875), (0.385, 0.430)),
            ("chain", (1009, 913), (0.440, 0.510)),
            ("chain", (1363, 942), (0.525, 0.545)),
            ("chain", (1576, 910), (0.553, 0.563)),
        ),
    ),
    *_segments(
        "fig4",
        "fig4_ratio_V",
        "Ne",
        (
            ("full", (625, 1277), (0.400, 0.615)),
            ("full", (1060, 1218), (0.625, 0.720)),
            ("full", (1353, 1201), (0.739, 0.758)),
            ("full", (1493, 1290), (0.763, 0.772)),
        ),
    ),
    *_segments(
        "fig4",
        "fig4_ratio_L",
        "Ne",
        (
            ("chain", (625, 985), (0.420, 0.455)),
            ("chain", (719, 1056), (0.465, 0.590)),
            ("chain", (1060, 1093), (0.625, 0.720)),
            ("chain", (1353, 1055), (0.739, 0.758)),
            ("chain", (1493, 1032), (0.763, 0.772)),
        ),
    ),
    *_segments(
        "fig5",
        "fig5_ratio_V",
        "Ar",
        (
            ("full", (634, 1995), (0.290, 0.430)),
            ("full", (1020, 1943), (0.445, 0.470)),
        ),
    ),
    *_segments(
        "fig5",
        "fig5_ratio_L",
        "Ar",
        (
            ("chain", (634, 2022), (0.290, 0.430)),
            ("chain", (1020, 1970), (0.445, 0.500)),
        ),
    ),
)

# The samples (hartree) that cannot be read cleanly, by file.
LEFT_OUT: dict[str, tuple[float, ...]] = {
    _file("fig1_sigma2_V", "Ne"): (
        0.4,
        0.405,
        0.415,
        0.42,
        0.425,
        0.6,
        0.605,
        0.61,
        0.615,
        0.625,
        0.63,
        0.635,
        0.64,
        0.645,
        0.65,
        0.655,
        0.66,
        0.665,
        0.67,
        0.72,
        0.739,
        0.74,
        0.741,
        0.742,
        0.756,
        0.757,
        0.758,
    ),
    _file("fig1_sigma2_L", "Ne"): (
        0.4,
        0.405,
        0.41,
        0.415,
        0.42,
        0.425,
        0.43,
        0.435,
        0.625,
        0.63,
        0.635,
        0.64,
        0.645,
        0.65,
        0.655,
        0.739,
        0.74,
        0.741,
        0.742,
        0.757,
        0.758,
    ),
    _file("fig2_sigma2_V", "Ar"): (
        0.38,
        0.385,
        0.43,
        0.44,
        0.445,
        0.51,
        0.525,
        0.526,
        0.527,
        0.528,
        0.542,
        0.543,
        0.544,
        0.545,
        0.553,
        0.554,
        0.563,
    ),
    _file("fig2_sigma2_L", "Ar"): (
        0.29,
        0.385,
        0.43,
        0.44,
        0.525,
        0.526,
        0.527,
        0.528,
        0.545,
        0.553,
        0.554,
        0.555,
        0.556,
        0.562,
        0.563,
    ),
    _file("fig4_ratio_V", "Ne"): (
        0.585,
        0.59,
        0.595,
        0.6,
        0.605,
        0.61,
        0.615,
        0.625,
        0.63,
        0.715,
        0.72,
        0.739,
        0.74,
        0.752,
        0.753,
        0.754,
        0.755,
        0.756,
        0.757,
        0.758,
        0.763,
        0.771,
        0.772,
    ),
    _file("fig4_ratio_L", "Ne"): (
        0.59,
        0.625,
        0.71,
        0.715,
        0.72,
        0.749,
        0.75,
        0.751,
        0.752,
        0.753,
        0.754,
        0.755,
        0.756,
        0.757,
        0.758,
        0.771,
        0.772,
    ),
    _file("fig5_ratio_V", "Ar"): (
        0.29,
        0.295,
        0.3,
        0.305,
        0.39,
        0.395,
        0.4,
        0.405,
        0.41,
        0.415,
        0.42,
        0.425,
        0.43,
        0.445,
        0.47,
    ),
    _file("fig5_ratio_L", "Ar"): (
        0.29,
        0.295,
        0.3,
        0.305,
        0.42,
        0.425,
        0.43,
        0.445,
        0.5,
    ),
}
# The samples read on both curves of a figure within ``LINE`` rows: left out of both.
TOGETHER: dict[str, tuple[float, ...]] = {}

# The tracing (:func:`trace`), pixels: the centres a prediction is fitted through, the
# fewest for a fit, how far a run may lie from the prediction, the widest run that
# moves the trace.
FIT = 10
MIN_FIT = 4
STEP_PIXELS = 4
MAX_WIDTH = 7
# The line style (:func:`style`), columns: the shortest unbroken run of a full curve,
# the longest dash of a chain curve, the fewest gaps of a chain curve per 100 columns.
FULL_RUN = 60
DASH = 50
CHAIN_GAPS = 4.0
# The reading (:func:`_read`), pixels: the clearance from any other dark pixel, a line's
# width and slack, the columns either side that must be clear too, the widest gap
# interpolated across, the steepest slope (rows per column); the samples' step, hartree,
# on the left and on the right of the bar.
CLEAR = 4
LINE = 4
SLACK = 2
WINDOW = 4
GAP = 24
MAX_SLOPE = 2.0
# Figs. 4 and 5: a value within this many rows of the drawn R = 1.5 line cannot be told
# from it.
NEAR = 6
STEP = 0.005
STEP_HIGH = 0.001


def trace(
    dark: np.ndarray, seed: tuple[int, int], columns: tuple[int, int]
) -> dict[int, Run | None]:
    """Return the run a curve passes through at each column, traced from a seed.

    From the seed forward, then back, column by column: the position predicted by a
    straight line through the last ``FIT`` traced centres; the run whose nearest pixel
    is within ``STEP_PIXELS`` of it is the curve's (the nearest, if several), and one no
    wider than ``MAX_WIDTH`` moves the trace. Columns with no such run are None.
    """
    out: dict[int, Run | None] = {}
    first, last = columns
    start = [(float(seed[0]), float(seed[1]))]
    for direction in (1, -1):
        history = list(start)
        column = seed[0]
        while first <= column <= last:
            recent = history[-FIT:]
            if len(recent) >= MIN_FIT:
                xs = np.array([x for x, _ in recent])
                ys = np.array([y for _, y in recent])
                predicted = float(np.polyval(np.polyfit(xs, ys, 1), column))
            else:
                predicted = recent[-1][1]
            low = int(predicted) - STEP_PIXELS - 3 * MAX_WIDTH
            pixels = dark[
                low : int(predicted) + STEP_PIXELS + 3 * MAX_WIDTH + 1, column
            ]
            near = [
                (abs((a + b) / 2 - predicted), (a, b))
                for a, b in ((low + a, low + b) for a, b in _spans(pixels))
                if abs(min(max(predicted, a), b) - predicted) <= STEP_PIXELS
            ]
            out[column] = min(near)[1] if near else None
            run = out[column]
            if run is not None and run[1] - run[0] + 1 <= MAX_WIDTH:
                history.append((float(column), _centre(run)))
            column += direction
        start = history[:FIT][::-1]
    return out


def style(found: list[bool]) -> str:
    """Return a traced stretch's style from which of its columns have a run.

    ``full``: an unbroken run of at least ``FULL_RUN`` columns (the curve's ends, where
    it turns steep into a resonance, may be lost). ``chain``: no run longer than
    ``DASH`` and, between the first and the last run, at least ``CHAIN_GAPS`` gaps per
    100 columns. Anything else is ``unknown``.
    """
    pieces = [(flag, len(list(group))) for flag, group in itertools.groupby(found)]
    while pieces and not pieces[0][0]:
        pieces.pop(0)
    while pieces and not pieces[-1][0]:
        pieces.pop()
    runs = [n for flag, n in pieces if flag]
    if not runs:
        return "unknown"
    if max(runs) >= FULL_RUN:
        return "full"
    gaps = sum(1 for flag, _ in pieces if not flag)
    width = sum(n for _, n in pieces)
    if max(runs) <= DASH and gaps >= CHAIN_GAPS * width / 100:
        return "chain"
    return "unknown"


def _clear(dark: np.ndarray, column: int, run: Run, slope: float) -> bool:
    """Return True if a run is thin and no other dark pixel is within ``CLEAR``."""
    a, b = run
    if b - a + 1 > LINE * np.hypot(1.0, slope) + SLACK:
        return False
    above = dark[a - CLEAR : a, column]
    below = dark[b + 1 : b + 1 + CLEAR, column]
    return not (above.any() or below.any())


@dataclass(frozen=True)
class Reading:
    """A value read from a figure: the sample, the pixel, the energy and the value."""

    wanted: float
    column: int
    row: float
    energy: float
    value: float
    printed: str


def _usable(dark: np.ndarray, line: dict[int, Run | None]) -> dict[int, float]:
    """Return the centres of a traced line's runs that are clear, with their window."""

    def slope(c: int) -> float:
        before, after = line.get(c - 3), line.get(c + 3)
        if before is None or after is None:
            return 0.0
        return (_centre(after) - _centre(before)) / 6

    clear = {
        c: _clear(dark, c, run, slope(c)) for c, run in line.items() if run is not None
    }
    open_ = {c: run is None or clear[c] for c, run in line.items()}
    return {
        c: _centre(run)
        for c, run in line.items()
        if run is not None
        and clear[c]
        and all(open_.get(d, True) for d in range(c - WINDOW, c + WINDOW + 1))
    }


def _read(
    calibration: Calibration, usable: dict[int, float], wanted: float
) -> Reading | None:
    """Return the curve's value at a photon energy, or None if it cannot be read."""
    column = round(calibration.column(wanted))
    window = [c for c in range(column - 6, column + 7) if c in usable]
    if len(window) >= 4:
        slope = abs(float(np.polyfit(window, [usable[c] for c in window], 1)[0]))
        if slope > MAX_SLOPE:
            return None
    if column in usable:
        row = usable[column]
        printed = f"pixel ({column}, {row:.1f})"
    else:
        left = next((c for c in range(column - 1, column - GAP, -1) if c in usable), 0)
        right = next((c for c in range(column + 1, column + GAP) if c in usable), 0)
        if not (left and right) or right - left > GAP:
            return None
        if (left < calibration.figure.bar) != (right < calibration.figure.bar):
            return None
        weight = (column - left) / (right - left)
        row = (1 - weight) * usable[left] + weight * usable[right]
        printed = (
            f"between pixels ({left}, {usable[left]:.1f}) and "
            f"({right}, {usable[right]:.1f})"
        )
    return Reading(
        wanted,
        column,
        row,
        calibration.energy(column),
        calibration.value(row),
        printed,
    )


def _grid(calibration: Calibration, span: tuple[float, float]) -> list[float]:
    """Return a span's samples: every ``STEP`` left of the bar, else ``STEP_HIGH``."""
    left = calibration.column(span[0]) < calibration.figure.bar
    step = STEP if left else STEP_HIGH
    first, last = (round(x / step) for x in span)
    return [round(k * step, 4) for k in range(first, last + 1)]


def calibrations() -> dict[str, tuple[Calibration, np.ndarray]]:
    """Return every figure's calibration and dark pixels, after checking its ticks."""
    pages: dict[int, np.ndarray] = {}
    out = {}
    for name, figure in FIGURES.items():
        if figure.page not in pages:
            pages[figure.page] = image(DIRECTORY, figure.page) < DARK
        _check_ticks(pages[figure.page], figure)
        out[name] = (Calibration(figure, pages[figure.page]), pages[figure.page])
    return out


def _follow(
    dark: np.ndarray, calibration: Calibration, segment: Segment, strict: bool
) -> dict[int, Run | None]:
    """Return a segment's traced line; raise if its style is not the one written."""
    first, last = (round(calibration.column(e)) for e in segment.span)
    line = trace(dark, segment.seed, (first - WINDOW, last + WINDOW))
    found = [line[c] is not None for c in range(first, last + 1)]
    drawn = style(found)
    if drawn != segment.style:
        msg = f"{segment.file} {segment.span}: a {drawn} line, not {segment.style}"
        if strict:
            raise ValueError(msg)
        print(msg)
    return line


def _together(out: dict[str, list[Reading]]) -> dict[str, list[float]]:
    """Return the samples where two curves of a figure are read within ``LINE`` rows."""
    figure_of = {segment.file: segment.figure for segment in SEGMENTS}
    found: dict[str, set[float]] = {}
    for a, b in itertools.combinations(sorted(out), 2):
        if figure_of[a] != figure_of[b]:
            continue
        rows_b = {r.wanted: r.row for r in out[b]}
        for r in out[a]:
            if r.wanted in rows_b and abs(r.row - rows_b[r.wanted]) < LINE:
                found.setdefault(a, set()).add(r.wanted)
                found.setdefault(b, set()).add(r.wanted)
    return {name: sorted(samples) for name, samples in found.items()}


def read_all(strict: bool = True) -> dict[str, list[Reading]]:
    """Return every segment's readings, by data file, in photon-energy order.

    Raises if a segment's line style is not the one written, or (if ``strict``) the
    samples left out are not exactly those of :data:`LEFT_OUT` and :data:`TOGETHER`.
    """
    scales = calibrations()
    out: dict[str, list[Reading]] = {}
    missing: dict[str, list[float]] = {}
    for segment in SEGMENTS:
        calibration, dark = scales[segment.figure]
        usable = _usable(dark, _follow(dark, calibration, segment, strict))
        for wanted in _grid(calibration, segment.span):
            reading = _read(calibration, usable, wanted)
            figure = calibration.figure
            if reading and abs(reading.row - figure.line_row(reading.column)) < NEAR:
                reading = None
            if reading is None:
                missing.setdefault(segment.file, []).append(wanted)
            else:
                out.setdefault(segment.file, []).append(reading)
    together = _together(out)
    for name, samples in together.items():
        out[name] = [r for r in out[name] if r.wanted not in samples]
    for name in set(missing) | set(LEFT_OUT) | set(together) | set(TOGETHER):
        if strict and tuple(missing.get(name, ())) != LEFT_OUT.get(name, ()):
            msg = f"{name}: left out {missing.get(name)}, not {LEFT_OUT.get(name)}"
            raise ValueError(msg)
        if strict and tuple(together.get(name, ())) != TOGETHER.get(name, ()):
            msg = f"{name}: together {together.get(name)}, not {TOGETHER.get(name)}"
            raise ValueError(msg)
    if not strict:
        print("LEFT_OUT", {k: tuple(v) for k, v in missing.items()})
        print("TOGETHER", {k: tuple(v) for k, v in together.items()})
    return {name: sorted(rows, key=lambda r: r.wanted) for name, rows in out.items()}


# Table 1 (p. 8/2743): the experimental excitation energies (hartree) the resonances
# were put at, by gas; with the two-photon thresholds (half the ionisation potentials)
# and the one-photon thresholds, the energies the vertical broken lines mark.
TABLE_1 = "tbl_c884f4d96ac8"
# The drawn R = 1.5 line's value on the y scale may stray this far from 1.5.
RATIO_LINE = 0.03
# A column whose dark pixels fill at least this share of the plot's height, in runs no
# longer than ``BROKEN_RUN`` rows, is a vertical broken line.
BROKEN_SHARE = 0.4
BROKEN_RUN = 40
COMPARE_EV = {"Ne": (11.0, 12.0, 13.0, 14.0, 15.0, 16.0), "Ar": (8.0, 9.0, 10.0, 11.0)}
OTHERS = {
    "Ne": {"mckenna_2004": "mckenna_2004/fig1_sigma2_Ne"},
    "Ar": {
        "mckenna_2004": "mckenna_2004/fig5_sigma2_Ar",
        "pan_1990b (all processes, L)": "pan_1990b/fig4c_all_L_Ar",
        "starace_1987 (ESI, L)": "starace_1987/fig7_esi_linear_Ar",
    },
}


def _table_1() -> dict[str, dict[str, float]]:
    """Return Table 1's experimental energies (hartree) by gas and excitation."""
    path = LIBRARY / DIRECTORY / "tables" / f"{TABLE_1}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))["cells"]
    grid = {(c["row"], c["column"]): c["text"].strip() for c in cells}
    out: dict[str, dict[str, float]] = {"Ne": {}, "Ar": {}}
    for row in range(2, 9):
        for gas, label, value in (("Ne", 0, 2), ("Ar", 3, 5)):
            if grid.get((row, value)):
                out[gas][grid[(row, label)]] = float(grid[(row, value)])
    return out


def broken_lines(calibration: Calibration, dark: np.ndarray) -> list[float]:
    """Return the photon energies (hartree) of the vertical broken lines found."""
    figure = calibration.figure
    top = round(min(r for r, _ in figure.y_ticks))
    columns = []
    for column in range(figure.axis_column + 3, 1800):
        runs = _spans(dark[top : figure.axis_row - 3, column])
        filled = sum(b - a + 1 for a, b in runs)
        if (
            filled >= BROKEN_SHARE * (figure.axis_row - 3 - top)
            and runs
            and max(b - a + 1 for a, b in runs) <= BROKEN_RUN
        ):
            columns.append(column)
    groups = np.split(np.array(columns), np.flatnonzero(np.diff(columns) > 3) + 1)
    return [round(calibration.energy(float(g.mean())), 4) for g in groups if g.size]


def _values(readings: dict[str, list[Reading]]) -> dict[str, dict[float, float]]:
    """Return each dataset's values (cm^4 s, or the ratio) by sample, by short name."""
    units = {segment.file: FIGURES[segment.figure].unit for segment in SEGMENTS}
    return {
        name.split("__")[1].removesuffix(".csv"): {
            r.wanted: r.value * units[name] for r in rows
        }
        for name, rows in readings.items()
    }


def _dataset(key: str) -> Dataset:
    """Return another record's dataset by its key, ``<source>/<id>``."""
    source, _ = key.split("/")
    record = load_source(DATA_DIRECTORY / "sources" / f"{source}.toml")
    return next(d for d in record.datasets if d.key == key)


def comparisons(values: dict[str, dict[float, float]]) -> dict[str, object]:
    """Return Figs. 1 and 2 against the other calculations of the library.

    At the samples nearest :data:`COMPARE_EV` (this paper over the other; the photon
    energy in eV is the hartree times ``HARTREE_EV``).
    """
    out: dict[str, object] = {}
    for gas, others in OTHERS.items():
        figure = "fig1" if gas == "Ne" else "fig2"
        for gauge in ("V", "L"):
            ours = values[f"{figure}_sigma2_{gauge}_{gas}"]
            for label, key in others.items():
                other = _dataset(key)
                energy, sigma = other.abscissa_values(), other.scaled_values()
                at = []
                nearest = {
                    min(ours, key=lambda x: abs(x * HARTREE_EV - wanted))
                    for wanted in COMPARE_EV[gas]
                }
                for e in sorted(nearest):
                    ev = e * HARTREE_EV
                    if energy.min() <= ev <= energy.max():
                        theirs = float(np.interp(ev, energy, sigma))
                        at.append((round(ev, 2), round(ours[e] / theirs, 2)))
                out[f"{figure} {gauge} / {label}"] = at
    return out


def checks(readings: dict[str, list[Reading]]) -> dict[str, object]:
    """Return the checks of the reading; raise if one fails its margin.

    The drawn R = 1.5 line of Figs. 4 and 5 on the y scale (:data:`RATIO_LINE`).
    Reported only: the tick residuals and pixel sizes; the vertical broken lines found
    against Table 1 and the thresholds; length over velocity; :func:`comparisons`.
    """
    out: dict[str, object] = {}
    failures: list[str] = []
    table_1 = _table_1()
    for name, (calibration, dark) in calibrations().items():
        figure = calibration.figure
        x = 1 / (
            float(calibration.axes["low"](1.0)) - float(calibration.axes["low"](0.0))
        )
        out[f"{name} ticks: worst residual, pixel x (hartree, low side), y"] = (
            round(calibration.residual(), 1),
            float(f"{1 / x:.2e}"),
            float(f"{calibration.axes['y'].resolution() * figure.unit:.2e}"),
        )
        expected = {
            "threshold": figure.threshold / 2,
            "ionisation potential": figure.threshold,
            **table_1[figure.material],
        }
        found = []
        for energy in broken_lines(calibration, dark):
            label = min(expected, key=lambda k: abs(expected[k] - energy))
            found.append((energy, label, round(energy - expected[label], 4)))
        out[f"{name} broken lines (found, nearest, found - it)"] = found
        if figure.ratio_line is not None:
            for column, row in figure.ratio_line:
                if not dark[row - 2 : row + 3, column].any():
                    failures.append(f"{name}: no R = 1.5 line at ({column}, {row})")
                value = calibration.value(row)
                out[f"{name} R = 1.5 line at column {column}"] = round(value, 3)
                if abs(value - 1.5) > RATIO_LINE:
                    failures.append(f"{name}: the R = 1.5 line reads {value:.3f}")
    values = _values(readings)
    for figure, gas in (("fig1", "Ne"), ("fig2", "Ar")):
        length, velocity = (values[f"{figure}_sigma2_{g}_{gas}"] for g in ("L", "V"))
        shared = sorted(set(length) & set(velocity))
        ratios = [length[e] / velocity[e] for e in shared]
        out[f"{figure} L/V"] = (
            shared[0],
            shared[-1],
            round(min(ratios), 2),
            round(max(ratios), 2),
            len(shared),
        )
    out.update(comparisons(values))
    if failures:
        raise ValueError("; ".join(failures))
    return out


def _points(readings: dict[str, list[Reading]]) -> dict[str, list[Point]]:
    """Return the data files' points: the value against photon energy (hartree)."""
    units = {segment.file: FIGURES[segment.figure].unit for segment in SEGMENTS}
    return {
        name: [
            Point(
                f"{r.wanted:.4f}",
                f"{r.energy:.4f}",
                f"{r.value * units[name]:.3e}"
                if "sigma2" in name
                else f"{r.value:.3f}",
                "",
                r.printed,
            )
            for r in rows
        ]
        for name, rows in readings.items()
    }


def extract() -> dict[str, list[Point]]:
    """Return the cross sections and ratios read from Figs. 1, 2, 4 and 5."""
    readings = read_all()
    checks(readings)
    return _points(readings)


def main() -> None:
    """Write the data files and print the checks."""
    readings = read_all()
    results = checks(readings)
    write(_points(readings))
    for label, result in results.items():
        print(f"{label}: {result}")


if __name__ == "__main__":
    main()
