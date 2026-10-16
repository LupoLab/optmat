# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract van der Hart (2006): R-matrix Floquet rates of He, Ne and Ar at 390 nm.

The rates are in figures only. **Fig. 6** (p. 8, library figure ``fig_1efe76bb16f6``)
draws those of He (solid line), Ne (dashed) and Ar (dot-dashed) against the intensity
on one set of axes; it is an image embedded in the PDF (900 x 635 pixels, 300 dpi),
read here from its native pixels (data rule 8, TC-699, TC-701):

* the axes: the frame's lines and the inward tick marks, found in the image and
  checked against the positions written below (``X_TICKS``, ``Y_TICKS``): x from 0 to
  2.5e14 W/cm^2 in 0.25 steps (the frame's left and right lines are 0 and 2.5); y, a
  logarithmic axis, 1e14 (the frame's top line), 1e12, 1e10 and 1e8 s^-1.
* each curve is traced column by column from a seed pixel on it (``SEEDS``), following
  the dark run nearest to the previous one within ``STEP`` pixels and bridging the gaps
  of the dashes; the trace is sampled at ``SAMPLES``, intensities on the curve's smooth
  part, away from the resonance spikes, each value from the centre of the dark run at
  that column. The ``printed`` cell of each value gives its pixel (column, row).

One pixel is 0.0033e14 W/cm^2 in x and 3.4 % in the rate; the lines are 2-4 pixels
wide. The reading is checked against van der Hart 2007's vector copy of the Ar curve
(``vanderhart_2007/fig4_floquet_Ar``, TC-700).
"""

import numpy as np
from _figure import Axis, image
from _library import Point, write

from optmat.datasets import datasets

DIRECTORY = "VanderHart_2006_IonizationRatesHe"
PAGE = 8
# Tick positions in the native image, (column or row, value): the frame's centre lines
# and the inward tick marks.
X_TICKS = (
    (107.0, 0.0),
    (185.0, 0.25),
    (262.0, 0.5),
    (339.0, 0.75),
    (416.0, 1.0),
    (492.5, 1.25),
    (569.5, 1.5),
    (646.5, 1.75),
    (723.5, 2.0),
    (800.5, 2.25),
    (878.5, 2.5),
)
Y_TICKS = ((16.5, 1e14), (153.0, 1e12), (289.5, 1e10), (425.5, 1e8))
DARK = 128
# The interior of the frame, clear of the frame lines and the tick marks.
INSIDE_ROWS = (22, 548)
INSIDE_COLUMNS = (116, 870)
STEP = 12
# The widest run, pixels, a value is read from: wider is a steep stretch of the line.
MAX_WIDTH = 6
# The seeds: at 1.30e14 W/cm^2 the three curves cross the column once each, Ar highest
# and He lowest; each is seeded at its run there, by that order.
SEED_INTENSITY = 1.30
SEEDS = {"Ar": 0, "Ne": 1, "He": 2}
# How far, pixels, a sample in a gap of the dashes (or on a steep stretch) looks for the
# traced line on either side, to interpolate between.
GAP = 10
# Where each curve is sampled, 1e14 W/cm^2: its smooth part, clear of the spikes (He
# below 1.05, Ne below 0.90, Ar below 0.65) and of the frame (He at 2.5).
SAMPLES = {
    "He": tuple(round(0.05 * k, 2) for k in range(21, 50)),
    "Ne": tuple(round(0.05 * k, 2) for k in range(18, 43) if k not in {20, 29}),
    "Ar": tuple(round(0.05 * k, 2) for k in range(13, 39) if k != 18),
}
# Left out: the peaks of Ne's resonances at 1.01 and 1.45e14 (2p5 4s) and of Ar's
# 3s2 3p5 3d (0.90e14), where the curve changes by some 40 % within a pixel's width
# (TC-701).
FILES = {
    "He": "He/vanderhart_2006__fig6_floquet_He.csv",
    "Ne": "Ne/vanderhart_2006__fig6_floquet_Ne.csv",
    "Ar": "Ar/vanderhart_2006__fig6_floquet_Ar.csv",
}


def _spans(mask: np.ndarray) -> list[tuple[int, int]]:
    """Return the first and last index of each run of True."""
    index = np.flatnonzero(mask)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [(int(group[0]), int(group[-1])) for group in groups]


def _runs(mask: np.ndarray) -> list[float]:
    """Return the centres of the runs of True in a column."""
    return [(first + last) / 2 for first, last in _spans(mask)]


def _check_ticks(dark: np.ndarray) -> None:
    """Raise unless the frame and ticks are where ``X_TICKS`` and ``Y_TICKS`` say."""
    near_bottom = dark[551]
    found_x = {round(c * 2) / 2 for c in _runs(near_bottom[116:870])}
    expected_x = {c - 116 for c, _ in X_TICKS[1:-1]}
    if found_x != expected_x:
        msg = f"Fig. 6: x ticks at {sorted(found_x)}, not {sorted(expected_x)}"
        raise ValueError(msg)
    near_left = dark[:, 114]
    found_y = {round(r * 2) / 2 for r in _runs(near_left[22:548])}
    expected_y = {r - 22 for r, _ in Y_TICKS[1:]}
    if found_y != expected_y:
        msg = f"Fig. 6: y ticks at {sorted(found_y)}, not {sorted(expected_y)}"
        raise ValueError(msg)
    for row in (15, 16, 17, 18, 555, 556, 557, 558):
        if dark[row, 107:879].mean() < 0.95:
            msg = f"Fig. 6: no frame line at row {row}"
            raise ValueError(msg)


def _trace(dark: np.ndarray, column: int, row: float) -> dict[int, tuple[int, int]]:
    """Return the run a curve passes through at each column, traced from a seed.

    From the seed pixel both ways, the run whose nearest pixel is closest to the
    previous position within ``STEP`` pixels; the position moves to that nearest pixel,
    so that a steep stretch (a long run) is followed up or down. Columns with no such
    run (a gap of the dashes) are skipped.
    """
    top, bottom = INSIDE_ROWS
    trace: dict[int, tuple[int, int]] = {}
    for direction in (1, -1):
        current = float(row)
        c = column
        while INSIDE_COLUMNS[0] <= c <= INSIDE_COLUMNS[1]:
            spans = [(top + a, top + b) for a, b in _spans(dark[top:bottom, c])]
            nearest = [(min(max(current, a), b), (a, b)) for a, b in spans]
            near = [(p, span) for p, span in nearest if abs(p - current) <= STEP]
            if near:
                current, span = min(near, key=lambda item: abs(item[0] - current))
                trace[c] = span
            c += direction
    return trace


def extract() -> dict[str, list[Point]]:
    """Return the Floquet rates of He, Ne and Ar sampled on Fig. 6's curves."""
    grey = image(DIRECTORY, PAGE)
    dark = grey < DARK
    if dark.shape != (635, 900):
        msg = f"Fig. 6: an image of {dark.shape}, not 635 x 900"
        raise ValueError(msg)
    _check_ticks(dark)
    x_axis = Axis(
        tuple(c for c, _ in X_TICKS), tuple(v for _, v in X_TICKS), tolerance=1.0
    )
    y_axis = Axis(
        tuple(r for r, _ in Y_TICKS),
        tuple(v for _, v in Y_TICKS),
        log=True,
        tolerance=1.0,
    )
    x_axis.check()
    y_axis.check()
    slope = (X_TICKS[-1][0] - X_TICKS[0][0]) / (X_TICKS[-1][1] - X_TICKS[0][1])
    seed = round(X_TICKS[0][0] + SEED_INTENSITY * slope)
    top, bottom = INSIDE_ROWS
    crossing = [top + (a + b) / 2 for a, b in _spans(dark[top:bottom, seed])]
    if len(crossing) != len(SEEDS):
        msg = f"Fig. 6: {len(crossing)} lines at the seed column, not {len(SEEDS)}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for gas, rank in SEEDS.items():
        trace = _trace(dark, seed, crossing[rank])
        points = []
        for wanted in SAMPLES[gas]:
            column = round(X_TICKS[0][0] + wanted * slope)
            centre, printed = _read(trace, column)
            points.append(
                Point(
                    f"{wanted:.2f}",
                    f"{float(x_axis(column)):.4f}e14",
                    f"{float(y_axis(centre)):.3e}",
                    "",
                    printed,
                )
            )
        files[FILES[gas]] = points
    _check_against_vector(files[FILES["Ar"]])
    return files


# How far the Ar curve read from this image may stray from van der Hart 2007's vector
# copy of it (it is 0.97 to 1.04 away from the resonance peaks).
AGREEMENT = 0.05


def _check_against_vector(argon: list[Point]) -> None:
    """Raise unless the Ar samples match the vector copy of the curve (TC-700)."""
    copy = next(d for d in datasets() if d.key == "vanderhart_2007/fig4_floquet_Ar")
    x = copy.abscissa_values()
    log_rate = np.log(copy.scaled_values())
    for point in argon:
        intensity = float(point.x) * 1e4  # W/cm2 to W/m2
        near = np.abs(x - intensity) < 0.004e18
        if not near.any():
            continue
        reference = float(np.exp(np.median(log_rate[near])))
        if abs(float(point.value) / reference - 1) > AGREEMENT:
            msg = (
                f"Fig. 6 reads Ar at {point.row}e14 as {point.value}, the vector copy "
                f"as {reference:.3e}"
            )
            raise ValueError(msg)


def _read(trace: dict[int, tuple[int, int]], column: int) -> tuple[float, str]:
    """Return a traced line's row at a column, and how it was read.

    The centre of the line's run there if it is no wider than ``MAX_WIDTH``; else (a gap
    or a steep stretch) the linear interpolation between the nearest such runs within
    ``GAP`` pixels on either side.
    """

    def usable(c: int) -> bool:
        return c in trace and trace[c][1] - trace[c][0] + 1 <= MAX_WIDTH

    def centre(c: int) -> float:
        return (trace[c][0] + trace[c][1]) / 2

    if usable(column):
        return centre(column), f"pixel ({column}, {centre(column):.1f})"
    left = next((c for c in range(column - 1, column - GAP - 1, -1) if usable(c)), None)
    right = next((c for c in range(column + 1, column + GAP + 1) if usable(c)), None)
    if left is None or right is None:
        msg = f"Fig. 6: no line within {GAP} pixels either side of column {column}"
        raise ValueError(msg)
    weight = (column - left) / (right - left)
    row = (1 - weight) * centre(left) + weight * centre(right)
    between = f"({left}, {centre(left):.1f}) and ({right}, {centre(right):.1f})"
    return row, f"between pixels {between}"


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
