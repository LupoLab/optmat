# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Dura et al. (2012): O2/Xe and N2/Ar ionisation ratios against wavelength.

The measured ratios are printed only in Fig. 4 (the text gives round numbers, quoted in
the record); they are read from the figure (data rule 8, TC-699). Fig. 4 (p. 2665, PDF
page 4, library figure ``fig_13703c126cd6``) is an image embedded in the PDF (1659 x
1004 pixels, 8-bit grey, lossless; 600 pixels per inch on the page), read at its native
pixels. Taken: the experimental points, 'gray diamonds, O2/Xe ratio; black circles,
N2/Ar ratio; error bars, standard deviation' (the caption), at 3 x 10^13 W/cm^2. Not
taken: the S-matrix curves and the open symbols ('experimental results from literature':
Guo et al., ref. 12, at 0.8 um; Walsh et al., ref. 11, at 10.6 um).

How the image is read:

* ink: a pixel darker than ``INK`` (grey level, 0 black to 255 white); the diamonds
  and their bars are mid-grey (about 110), the circles and theirs black.
* the axes: the frame is black, 7 pixels wide. The ratio axis (logarithmic) is
  calibrated from the tick marks along column ``Y_TICKS_COLUMN``, just inside the left
  frame: 29 runs of ink between the frame's top and bottom lines, 3 to 0.002 (the
  values ``Y_VALUES``; the long ones, black at column ``Y_LONG_COLUMN``, must be the
  labelled 1, 0.1 and 0.01). The wavelength axis (logarithmic, broken after the left
  panel) from the ticks along row ``X_TICKS_ROW`` of the left panel: four, read
  as 800, 900, 1000 and 2000 nm (only 1000 is labelled, and is the only long tick,
  black at row ``X_LONG_ROW``; 800 and 900 are the minor ticks below it, 2000 the
  panel's end; that the four lie on one log scale within ``TICK_TOLERANCE`` checks
  the reading). Each tick's position is the ink-weighted centre of its run.
* the markers: each is seeded at a pixel found by inspecting the image
  (``MARKERS``); the bar's column is the seed's. The marker's rows are the contiguous
  rows about the seed whose run of ink through the bar's column (within
  ``WINDOW`` pixels) is at least ``MARKER_WIDTH`` wide (the bar alone is 4); the
  centre row is the middle of the first and last, and the centre column the mean of
  the runs' ink-weighted centres in those two rows. A circle or diamond partly covered
  by a line or by the caps keeps its outer rows: at 1.4 um the N2/Ar circle lies on the
  N2/Ar curve and its extent is that of the circle's top and bottom (the curve is
  narrower than the circle there).
* the error bars: from the centre outward along the bar's three middle columns while
  they hold ink, the rows whose run is ``CAP_WIDTH`` wide (a cap is 21-23 pixels, a
  curve crossing it longer); the end is the mean row of the outermost such group.
  Where a cap cannot be told from a curve it lies on (the record's notes say which),
  that end is not read (``MARKERS``), and ``sigma`` is the half-length of the other
  half-bar in the ratio; where neither, empty.

``sigma`` is half the bar's full length in the ratio (the bars are 'standard
deviation'); on the logarithmic axis a bar symmetric in the ratio is longer below the
marker than above (checked: :func:`bar_symmetry`). Each ``printed`` cell gives the
marker's pixel (column, row) and the rows of its bar's ends.
"""

import math
from dataclasses import dataclass

import numpy as np
from _figure import Axis, image
from _library import Point, write

DIRECTORY = "Dura_2011_WavelengthDependenceSuppressed"
SOURCE = "dura_2012"
PAGE, INDEX, SHAPE = 4, 0, (1004, 1659)
INK = 200
# The ratio's ticks: column 250, rows 10-856 (inside the frame's top and bottom lines,
# rows 0-6 and 858-864); the long ones reach column 262.
Y_TICKS_COLUMN, Y_LONG_COLUMN, Y_ROWS = 250, 262, (10, 856)
Y_VALUES = (
    3.0,
    2.0,
    1.0,
    *(k / 10 for k in range(9, 0, -1)),
    *(k / 100 for k in range(9, 0, -1)),
    *(k / 1000 for k in range(9, 1, -1)),
)
Y_LONG = (1.0, 0.1, 0.01)
# The wavelength's ticks: row 850, columns 242-1340 (the left panel); the long one
# reaches row 836 (the short ones 846). The dotted grid lines continue the ticks up
# the plot in grey: a long tick is black (``BLACK``) at row 840.
X_TICKS_ROW, X_LONG_ROW, X_COLUMNS = 850, 840, (242, 1340)
BLACK = 60
X_VALUES = (800.0, 900.0, 1000.0, 2000.0)
X_LONG = (1000.0,)
TICK_TOLERANCE = 0.6  # pixels
WINDOW = 16  # pixels either side of the bar's column searched for a run
MARKER_WIDTH = 7  # pixels: a row of the marker (the bar alone is 4 wide)
CAP_WIDTH = (18, 25)  # pixels: a cap's run (21-23), shorter than a crossing curve


@dataclass(frozen=True)
class Seed:
    """A marker seeded by inspection, and which of its bar's ends are read.

    Attributes
    ----------
    column, row : int
        A pixel on the marker, on its bar's column.
    top, bottom : bool
        Whether the upper and the lower end of the bar are read.
    """

    column: int
    row: int
    top: bool
    bottom: bool


MARKERS: dict[str, tuple[Seed, ...]] = {
    # N2/Ar, black circles. At 1.4 um the caps lie on the circle and on the N2/Ar curve
    # (not read); at 1.5 um the upper cap lies on the curve (not read).
    "N2": (
        Seed(325, 171, top=True, bottom=True),
        Seed(861, 170, top=True, bottom=True),
        Seed(942, 148, top=False, bottom=False),
        Seed(1018, 208, top=False, bottom=True),
    ),
    # O2/Xe, grey diamonds.
    "O2": (
        Seed(325, 362, top=True, bottom=True),
        Seed(860, 286, top=True, bottom=True),
        Seed(942, 272, top=True, bottom=True),
    ),
}


def _runs(line: np.ndarray, start: int) -> list[float]:
    """Return the ink-weighted centres of the runs of ink along a line."""
    index = np.flatnonzero(line < INK)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    weight = 255 - line
    return [float(start + (g * weight[g]).sum() / weight[g].sum()) for g in groups]


def axes(grey: np.ndarray) -> tuple[Axis, Axis]:
    """Return the wavelength (nm) and ratio axes, calibrated from the tick marks."""
    first, last = Y_ROWS
    rows = _runs(grey[first:last, Y_TICKS_COLUMN], first)
    if len(rows) != len(Y_VALUES):
        msg = f"Fig. 4: {len(rows)} ratio ticks, not {len(Y_VALUES)}"
        raise ValueError(msg)
    long = [
        v
        for r, v in zip(rows, Y_VALUES, strict=True)
        if grey[round(r), Y_LONG_COLUMN] < BLACK
    ]
    if tuple(long) != Y_LONG:
        msg = f"Fig. 4: long ratio ticks at {long}"
        raise ValueError(msg)
    first, last = X_COLUMNS
    columns = _runs(grey[X_TICKS_ROW, first:last], first)
    if len(columns) != len(X_VALUES):
        msg = f"Fig. 4: {len(columns)} wavelength ticks, not {len(X_VALUES)}"
        raise ValueError(msg)
    long = [
        v
        for c, v in zip(columns, X_VALUES, strict=True)
        if grey[X_LONG_ROW, round(c)] < BLACK
    ]
    if tuple(long) != X_LONG:
        msg = f"Fig. 4: long wavelength ticks at {long}"
        raise ValueError(msg)
    x_axis = Axis(tuple(columns), X_VALUES, log=True, tolerance=TICK_TOLERANCE)
    y_axis = Axis(tuple(rows), Y_VALUES, log=True, tolerance=TICK_TOLERANCE)
    x_axis.check()
    y_axis.check()
    return x_axis, y_axis


def _run(grey: np.ndarray, row: int, column: int) -> tuple[int, float]:
    """Return the width and ink-weighted centre of the run of ink through a pixel."""
    line = grey[row, column - WINDOW : column + WINDOW + 1]
    middle = WINDOW
    if line[middle] >= INK:
        return 0, float(column)
    left = middle
    while left > 0 and line[left - 1] < INK:
        left -= 1
    right = middle
    while right < len(line) - 1 and line[right + 1] < INK:
        right += 1
    weight = 255 - line[left : right + 1]
    centre = (
        column
        - WINDOW
        + float((np.arange(left, right + 1) * weight).sum() / weight.sum())
    )
    return right - left + 1, centre


def _bar_end(grey: np.ndarray, column: int, row: float, direction: int) -> float:
    """Return the mean row of the outermost cap from a marker's centre (up -1)."""
    r = round(row)
    caps: list[int] = []
    while 0 <= r < grey.shape[0] and min(grey[r, column - 1 : column + 2]) < INK:
        width, _ = _run(grey, r, column)
        if CAP_WIDTH[0] <= width <= CAP_WIDTH[1]:
            caps.append(r)
        r += direction
    if not caps:
        msg = f"Fig. 4: no cap from ({column}, {row:.1f}) in direction {direction}"
        raise ValueError(msg)
    group = [caps[-1]]
    for c in reversed(caps[:-1]):
        if abs(c - group[-1]) != 1:
            break
        group.append(c)
    return float(np.mean(group))


@dataclass(frozen=True)
class Marker:
    """A marker read from the image: centre and bar ends (pixels; None: not read)."""

    x: float
    y: float
    top: float | None
    bottom: float | None


def _marker(grey: np.ndarray, seed: Seed) -> Marker:
    """Return a marker's centre and the ends of its bar that are read."""
    column = seed.column

    def wide(r: int) -> bool:
        return _run(grey, r, column)[0] >= MARKER_WIDTH

    if not wide(seed.row):
        msg = f"Fig. 4: the seed ({column}, {seed.row}) is not on a marker"
        raise ValueError(msg)
    first = last = seed.row
    while wide(first - 1):
        first -= 1
    while wide(last + 1):
        last += 1
    y = (first + last) / 2
    x = (_run(grey, first, column)[1] + _run(grey, last, column)[1]) / 2
    top = _bar_end(grey, column, y, -1) if seed.top else None
    bottom = _bar_end(grey, column, y, +1) if seed.bottom else None
    return Marker(x, y, top, bottom)


def readings() -> tuple[dict[str, list[Marker]], Axis, Axis]:
    """Return the markers of each pair, and the calibrated axes."""
    grey = image(DIRECTORY, PAGE, INDEX)
    if grey.shape != SHAPE:
        msg = f"Fig. 4: an image of {grey.shape}, not {SHAPE}"
        raise ValueError(msg)
    x_axis, y_axis = axes(grey)
    markers = {gas: [_marker(grey, s) for s in seeds] for gas, seeds in MARKERS.items()}
    return markers, x_axis, y_axis


def _point(k: int, m: Marker, x_axis: Axis, y_axis: Axis) -> Point:
    """Return a marker's wavelength (nm), ratio and sigma, with its pixels."""
    wavelength = float(x_axis(m.x))
    ratio = float(y_axis(m.y))
    where = f"pixel ({m.x:.1f}, {m.y:.1f})"
    if m.top is not None and m.bottom is not None:
        sigma = (float(y_axis(m.top)) - float(y_axis(m.bottom))) / 2
        printed = f"{where}; bar ends rows {m.top:.1f} and {m.bottom:.1f}"
    elif m.top is not None or m.bottom is not None:
        end = m.top if m.top is not None else m.bottom
        assert end is not None  # one of the two, by the branch
        sigma = abs(float(y_axis(end)) - ratio)
        hidden = "upper" if m.top is None else "lower"
        printed = f"{where}; one bar end, row {end:.1f} (the {hidden} not read)"
    else:
        return Point(
            str(k), f"{wavelength:.0f}", f"{ratio:.3f}", "", f"{where}; no bar end read"
        )
    return Point(str(k), f"{wavelength:.0f}", f"{ratio:.3f}", f"{sigma:.3f}", printed)


def bar_symmetry(markers: list[Marker], y_axis: Axis) -> list[float]:
    """Return, for bars with both ends read, the upper half over the lower in ratio."""
    out = []
    for m in markers:
        if m.top is None or m.bottom is None:
            continue
        value = float(y_axis(m.y))
        out.append((float(y_axis(m.top)) - value) / (value - float(y_axis(m.bottom))))
    return out


def extract() -> dict[str, list[Point]]:
    """Return the data files: O2/Xe and N2/Ar against wavelength."""
    markers, x_axis, y_axis = readings()
    return {
        f"{gas}/{SOURCE}__fig4_{gas}_{companion}.csv": [
            _point(k, m, x_axis, y_axis)
            for k, m in enumerate(sorted(markers[gas], key=lambda m: m.x), start=1)
        ]
        for gas, companion in (("O2", "Xe"), ("N2", "Ar"))
    }


def hidden_cap_row(m: Marker, y_axis: Axis) -> float:
    """Return the row where the unread end of a one-ended bar would be, if symmetric."""
    end = m.bottom if m.top is None else m.top
    assert end is not None  # a one-ended bar
    value = float(y_axis(m.y))
    other = 2 * value - float(y_axis(end))
    slope, intercept = np.polyfit(y_axis.positions, np.log10(y_axis.values), 1)
    return float((math.log10(other) - intercept) / slope)


def main() -> None:
    """Write the data files and print the checks the record quotes."""
    files = extract()
    markers, x_axis, y_axis = readings()
    for name, points in files.items():
        for p in points:
            print(name, p.row, p.x, p.value, p.sigma, p.printed)
    for gas in markers:
        print(
            gas,
            "bar symmetry (upper/lower, in ratio):",
            [f"{s:.2f}" for s in bar_symmetry(markers[gas], y_axis)],
        )
    for m in markers["N2"]:
        if (m.top is None) != (m.bottom is None):
            print(
                f"N2 one-ended bar at column {m.x:.0f}: the other end, if symmetric, "
                f"at row {hidden_cap_row(m, y_axis):.1f}"
            )
    print(
        f"pixels: {x_axis.resolution():.2%} of the wavelength, "
        f"{y_axis.resolution():.2%} of the ratio"
    )
    write(files)


if __name__ == "__main__":
    main()
