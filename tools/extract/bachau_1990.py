# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bachau, Lambropoulos and Tang (1990): two-photon cross sections of He.

The PDF is a scan (one 5060 x 6898 image per page, about 600 dpi) with a text layer.

* **Table I** (p. 5803, ``tbl_c78f5248bff6``, CSV; every cell checked against the
  table image, all agree): the generalised two-photon ionisation cross sections
  sigma_2 of He, cm^4 s, 'multiplied by 10^-53', at hbar omega = 0.5, 1.0, 1.1, 1.2
  and 1.3 au, in two columns: the authors' L^2-basis calculation and a single-channel
  quantum-defect Green's-function calculation (H. Rudolph's). 0.5 au is below the
  one-photon threshold; 1.0 to 1.3 au are above it, 'two-photon ATI' (the caption).
  One file per column; the row ids are the photon energies; ``printed`` keeps the row.
* **Fig. 1(a)** (p. 5802, ``fig_7d7bc8abfffe``): W / I^2 of two-photon ionisation of
  He, 10^-16 s^-1 W^-2 cm^4, against the photon energy in Ry, below the one-photon
  threshold, the authors' results the long-dashed line among four other authors'
  curves and markers. Read from the page image at its native pixels (data rule 8,
  TC-699): only where the long dashes stand clear of the other curves, 0.90 to 1.10
  Ry (the three dashes from the 2ph threshold's dotted line; beyond 1.10 Ry they run
  into the solid line and Victor's circles).

How Fig. 1(a) is read:

* the axes: the inward tick marks of the frame's bottom and left edges, written in
  :data:`X_TICKS` and :data:`Y_TICKS` and checked against the image (each found within
  :data:`FOUND` pixels; the long ticks at the labelled values 1, 1.5, 2 Ry and 2, 4, 6).
  x: every 0.1 Ry from 0.7 to 2.4 except 0.9 and 1.8, hidden by the dotted threshold
  lines; y: every 0.5 from 0 (the frame's bottom line) to 6.5. The scan is slightly
  sheared: the top edge's ticks (:data:`TOP_TICKS`) lie 6.5 to 10 pixels left of the
  bottom ones, so a pixel's column is moved to the bottom edge along the tick lines
  before the x calibration. The ticks are irregularly spaced (77 to 87 pixels per 0.1
  Ry, the same on both edges: the drawing's): the linear x scale fits them within
  :data:`X_TOLERANCE` pixels (0.011 Ry), the y scale within :data:`Y_TOLERANCE`.
* the long dashes: the connected dark regions (8-connected, grey below :data:`DARK`)
  inside :data:`WINDOW` that are :data:`DASH_HEIGHT` pixels high and
  :data:`DASH_WIDTH` wide (the short dashes are under 30 pixels high, the solid line
  with its circles is one region of some 180; checked: exactly three). The window
  starts right of the 2ph dotted line. Every :data:`EVERY` th column of a dash is
  sampled where the dash crosses it in one run of :data:`RUN` pixels (the dash's stroke,
  about 13 pixels on its slope): the rounded ends, and the columns where Kulander's
  '+' marker crosses the second dash at 1.0 Ry, are left out. The value is the run's
  centre; ``printed`` gives the pixel (column, row) and the run's length.

**Checks** (:func:`checks`): the text's 'largest discrepancy of a factor of 3.6
between the two results' (p. 5803); Fig. 1(a) at 1 Ry against Table I's 0.5 au entry
converted to W / I^2 = sigma_2 / (hbar omega)^2 (W = sigma_2 Phi^2 with Phi =
I / hbar omega; the paper prints no definition, this is the one the units imply).
Reported (:func:`nikolopoulos_2001`): Nikolopoulos and Lambropoulos (2001) draw Table
I's L^2 column as the squares of their Fig. 4 (nikolopoulos_2001, TC-875), read from
that vector drawing with nikolopoulos_2001's calibrated axes.
"""

import functools

import nikolopoulos_2001 as nikolopoulos
import numpy as np
from _figure import Axis, image
from _library import Point, table, write
from scipy import ndimage

DIRECTORY = "Bachau_1990_MethodCalculationMultiphoton"
SOURCE = "bachau_1990"
TABLE_I = "tbl_c78f5248bff6"
HEADER = ["\\hbar\\omega (a.u.)", "L^{2} basis", "Green's function"]
FILES = {
    1: "He/bachau_1990__tbl1_sigma2_L2_He.csv",
    2: "He/bachau_1990__tbl1_sigma2_qd_He.csv",
}
FIGURE_FILE = "He/bachau_1990__fig1a_W2_L2_He.csv"
SCALE = 1e-53
"""'The entries for sigma_2 must be multiplied by 10^-53' (Table I's caption)."""


@functools.cache
def _rows() -> tuple[tuple[str, ...], ...]:
    """Return Table I's rows, after checking its header."""
    grid = table(DIRECTORY, TABLE_I)
    if grid[0] != HEADER:
        msg = f"unexpected Table I header: {grid[0]!r}"
        raise ValueError(msg)
    for cells in grid[1:]:
        for cell in cells:
            float(cell)
    return tuple(tuple(cells) for cells in grid[1:])


def _points(column: int) -> list[Point]:
    """Return one column of Table I as data points (row id: the photon energy)."""
    return [
        Point(cells[0], cells[0], cells[column], "", " | ".join(cells))
        for cells in _rows()
    ]


# Fig. 1(a), on PDF page 2's image (native pixels: column, row).
PAGE = 2
DARK = 128
BOTTOM = 4113.0
"""The row of the frame's bottom line (y = 0) at the ticks' columns."""
TOP = 3077.0
"""The row of the frame's top line."""
X_TICKS = (
    (2974.0, 0.7),
    (3050.0, 0.8),
    (3222.0, 1.0),
    (3306.5, 1.1),
    (3393.0, 1.2),
    (3472.0, 1.3),
    (3558.5, 1.4),
    (3643.5, 1.5),
    (3720.0, 1.6),
    (3804.5, 1.7),
    (3966.0, 1.9),
    (4052.0, 2.0),
    (4130.0, 2.1),
    (4208.5, 2.2),
    (4294.5, 2.3),
    (4373.0, 2.4),
)
"""The bottom edge's ticks, column and Ry (0.9 and 1.8 lie under the dotted lines of
the experimental thresholds, 'expt.' 0.9 and 1.81)."""
TOP_TICKS = (
    2967.5,
    3045.0,
    3215.5,
    3299.5,
    3386.5,
    3465.0,
    3551.0,
    3635.5,
    3713.0,
    3797.0,
    3958.0,
    4043.5,
    4121.0,
    4199.0,
    4284.5,
    4363.5,
)
"""The top edge's ticks, column, at the same values as :data:`X_TICKS`."""
Y_TICKS = (
    (BOTTOM, 0.0),
    (4041.0, 0.5),
    (3965.5, 1.0),
    (3893.0, 1.5),
    (3820.5, 2.0),
    (3747.5, 2.5),
    (3671.0, 3.0),
    (3603.5, 3.5),
    (3529.5, 4.0),
    (3459.5, 4.5),
    (3384.5, 5.0),
    (3313.0, 5.5),
    (3237.0, 6.0),
    (3166.0, 6.5),
)
"""The left edge's ticks, row and W / I^2 (10^-16 s^-1 W^-2 cm^4)."""
LONG_X = (1.0, 1.5, 2.0)
LONG_Y = (2.0, 4.0, 6.0)
FOUND = 1.0
"""Each written tick is found in the image within this, pixels."""
X_TOLERANCE = 10.0
"""The x ticks lie on the linear scale within this, pixels (the drawing's spacing)."""
Y_TOLERANCE = 5.0
WINDOW = (3540, 3830, 3140, 3320)
"""Rows and columns (first, last) searched for the long dashes."""
DASH_HEIGHT = (50, 90)
DASH_WIDTH = (30, 70)
RUN = (9, 15)
"""A sampled column crosses the dash in one run of this many pixels."""
EVERY = 4


def _bands(dark: np.ndarray) -> dict[str, list[tuple[float, int]]]:
    """Return the tick marks found along the bottom, top and left edges.

    (position, length in pixels), the position the dark pixels' weighted centre.
    """
    found: dict[str, list[tuple[float, int]]] = {}
    for name, band, axis, start in (
        ("bottom", dark[4088:4103, 2905:4450], 0, 2905),
        ("top", dark[3088:3100, 2905:4450], 0, 2905),
        ("left", dark[3085:4108, 2903:2925], 1, 3085),
    ):
        counts = band.sum(axis=axis)
        need = 0.8 * band.shape[axis]
        labels, count = ndimage.label(counts >= need)
        marks = []
        for k in range(1, count + 1):
            index = np.flatnonzero(labels == k)
            centre = float(np.average(index, weights=counts[index])) + start
            if name == "left":
                row = dark[round(centre), 2895:3000]
                length = int(np.flatnonzero(~row[8:])[0]) + 8
            elif name == "bottom":
                column = dark[3900 : int(BOTTOM), round(centre)]
                length = int(np.flatnonzero(~column[::-1])[0])
            else:
                length = 0
            marks.append((centre, length))
        found[name] = marks
    return found


def _check_ticks(dark: np.ndarray) -> None:
    """Raise unless every written tick is in the image, the long ones at the labels."""
    found = _bands(dark)
    written = {
        "bottom": [(c, v) for c, v in X_TICKS],
        "top": [(c, v) for c, (_, v) in zip(TOP_TICKS, X_TICKS, strict=True)],
        "left": [(r, v) for r, v in Y_TICKS[1:]],
    }
    for name, ticks in written.items():
        for position, value in ticks:
            near = [m for m in found[name] if abs(m[0] - position) <= FOUND]
            if len(near) != 1:
                msg = f"Fig. 1(a): no {name} tick at {position} ({value})"
                raise ValueError(msg)
            length = near[0][1]
            long = value in (LONG_Y if name == "left" else LONG_X)
            if name != "top" and (length > 40) != long:
                msg = f"Fig. 1(a): the {name} tick at {value} is {length} pixels long"
                raise ValueError(msg)


@functools.cache
def figure_axes() -> tuple[Axis, Axis]:
    """Return Fig. 1(a)'s x axis (bottom-edge column, Ry) and y axis, checked."""
    dark = image(DIRECTORY, PAGE) < DARK
    _check_ticks(dark)
    x_axis = Axis(
        tuple(c for c, _ in X_TICKS),
        tuple(v for _, v in X_TICKS),
        tolerance=X_TOLERANCE,
    )
    y_axis = Axis(
        tuple(r for r, _ in Y_TICKS),
        tuple(v for _, v in Y_TICKS),
        tolerance=Y_TOLERANCE,
    )
    x_axis.check()
    y_axis.check()
    return x_axis, y_axis


def to_bottom(column: float, row: float) -> float:
    """Return the bottom-edge column of the tick line through a pixel (the shear)."""
    bottoms = np.array([c for c, _ in X_TICKS])
    shift = bottoms - np.array(TOP_TICKS)
    along = float(np.interp(column, bottoms, shift))
    return column + along * (BOTTOM - row) / (BOTTOM - TOP)


def residuals() -> tuple[float, float]:
    """Return the worst tick residual of Fig. 1(a)'s x and y axes, pixels."""
    out = []
    for axis in figure_axes():
        slope, intercept = np.polyfit(axis.positions, axis.values, 1)
        predicted = (np.asarray(axis.values) - intercept) / slope
        out.append(float(np.max(np.abs(predicted - np.asarray(axis.positions)))))
    return out[0], out[1]


@functools.cache
def dashes() -> tuple[tuple[tuple[int, float, int], ...], ...]:
    """Return the long dashes' samples: per dash, (column, run centre, run length)."""
    dark = image(DIRECTORY, PAGE) < DARK
    top, bottom, left, right = WINDOW
    labels, _ = ndimage.label(
        dark[top : bottom + 1, left : right + 1], structure=np.ones((3, 3))
    )
    found = []
    for k, box in enumerate(ndimage.find_objects(labels), 1):
        if box is None:
            continue
        height = box[0].stop - box[0].start
        width = box[1].stop - box[1].start
        if not (
            DASH_HEIGHT[0] <= height <= DASH_HEIGHT[1]
            and DASH_WIDTH[0] <= width <= DASH_WIDTH[1]
        ):
            continue
        mask = labels == k
        samples = []
        for c in range(box[1].start, box[1].stop, EVERY):
            rows = np.flatnonzero(mask[:, c])
            if rows.size == 0 or rows[-1] - rows[0] + 1 != rows.size:
                continue
            if RUN[0] <= rows.size <= RUN[1]:
                centre = (rows[0] + rows[-1]) / 2
                samples.append((c + left, float(centre + top), int(rows.size)))
        found.append(tuple(samples))
    found.sort(key=lambda s: s[0][0])
    if len(found) != 3:
        msg = f"Fig. 1(a): {len(found)} long dashes in the window, not 3"
        raise ValueError(msg)
    return tuple(found)


def figure_points() -> list[Point]:
    """Return Fig. 1(a)'s long-dashed line as data points, one per sampled column."""
    x_axis, y_axis = figure_axes()
    points = []
    for number, dash in enumerate(dashes(), 1):
        for column, row, run in dash:
            x = float(x_axis(to_bottom(column, row)))
            y = float(y_axis(row))
            points.append(
                Point(
                    f"d{number}c{column}",
                    f"{x:.4f}",
                    f"{y:.3f}",
                    "",
                    f"dash {number}, pixel ({column}, {row:.1f}), run {run}",
                )
            )
    return points


FACTOR = 3.6
"""p. 5803: 'The largest discrepancy of a factor of 3.6 between the two results'."""
HARTREE_J = 4.3597447222071e-18
"""The hartree in J (CODATA 2018)."""
AT_ONE_RY = 0.05
"""Fig. 1(a) at 1 Ry against Table I's 0.5 au entry: the reading's uncertainty there
(the x scale's 2-pixel residual at the 1 Ry tick, 0.6 %; the y scale's, 1 %) and the
interpolation across the '+' (the quadratic's residual, 0.4 %), with the drawing's."""


def at_one_rydberg() -> tuple[float, float]:
    """Return Fig. 1(a)'s long-dashed line at 1 Ry and Table I's 0.5 au entry, W / I^2.

    The line at 1 Ry, under Kulander's '+', interpolated by a quadratic in x fitted to
    the three dashes' samples.
    """
    points = figure_points()
    x = np.array([float(p.x) for p in points])
    y = np.array([float(p.value) for p in points])
    drawn = float(np.polyval(np.polyfit(x, y, 2), 1.0))
    sigma = float(_rows()[0][1]) * SCALE
    printed = sigma / (0.5 * HARTREE_J) ** 2 / 1e-16
    return drawn, printed


def checks() -> dict[str, object]:
    """Return the checks of Table I and Fig. 1(a); raise if one is missed."""
    ratios = {cells[0]: float(cells[1]) / float(cells[2]) for cells in _rows()}
    largest = max(ratios.values())
    # The text gives two figures: the largest ratio rounds to 3.6.
    if round(largest, 1) != FACTOR:
        msg = f"the largest L2 / Green's-function ratio is {largest:.3f}, not 3.6"
        raise ValueError(msg)
    drawn, printed = at_one_rydberg()
    if abs(drawn / printed - 1) > AT_ONE_RY:
        msg = f"Fig. 1(a) at 1 Ry is {drawn:.3f}, Table I's 0.5 au entry {printed:.3f}"
        raise ValueError(msg)
    return {
        "L2 over Green's function": ratios,
        "largest": largest,
        "Fig. 1(a) at 1 Ry, Table I at 0.5 au (W / I^2), ratio": (
            drawn,
            printed,
            drawn / printed,
        ),
    }


HARTREE_EV = 27.211386245988
"""The hartree in eV (CODATA 2018): nikolopoulos_2001's Fig. 4 is in eV."""
SQUARE = (51, 59, 151)
"""The fill of nikolopoulos_2001's Fig. 4 squares ('He+, Bachau et al (1990)')."""


def nikolopoulos_2001() -> dict[str, tuple[float, float, float]]:
    """Return nikolopoulos_2001's Fig. 4 squares against Table I's L^2 column.

    Per square, by Table I's photon energy (au): (its photon energy, eV; its value,
    cm^4 s; Table I's value at that energy). The squares are each drawn filled and
    outlined; the legend's square is left out.
    """
    figure = nikolopoulos.FIGURES["fig4"]
    left, bottom, right, top = nikolopoulos.frame("fig4")
    legend = figure.legend
    if legend is None:
        msg = "nikolopoulos_2001's Fig. 4 has no legend box"
        raise ValueError(msg)
    centres = set()
    for drawn in nikolopoulos._drawn(figure.page):
        if drawn.fill != SQUARE:
            continue
        west, south, east, north = drawn.box
        x, y = (west + east) / 2, (south + north) / 2
        inside = left <= x <= right and bottom <= y <= top
        if inside and not (legend[0] <= x <= legend[2] and legend[1] <= y <= legend[3]):
            centres.add((round(x, 3), round(y, 3)))
    x_axis, y_axis = nikolopoulos.axes("fig4")
    table_i = {float(cells[0]): float(cells[1]) * SCALE for cells in _rows()}
    out = {}
    for x, y in sorted(centres):
        energy = float(x_axis(x))
        omega = min(table_i, key=lambda w: abs(w * HARTREE_EV - energy))
        out[f"{omega} au"] = (energy, float(y_axis(y)), table_i[omega])
    if len(out) != len(centres):
        msg = "two of nikolopoulos_2001's squares match one row of Table I"
        raise ValueError(msg)
    return out


def extract() -> dict[str, list[Point]]:
    """Return Table I's two columns and Fig. 1(a)'s long dashes, after the checks."""
    checks()
    files = {relative: _points(column) for column, relative in FILES.items()}
    files[FIGURE_FILE] = figure_points()
    return files


def main() -> None:
    """Write the data files, and print the checks and the comparison."""
    write(extract())
    x_res, y_res = residuals()
    print(f"Fig. 1(a) ticks within {x_res:.1f} (x) and {y_res:.1f} (y) pixels")
    values = [(float(p.x), float(p.value)) for p in figure_points()]
    print(
        f"Fig. 1(a): {len(values)} samples, {values[0][0]:.3f} to {values[-1][0]:.3f} "
        f"Ry, {values[0][1]:.3f} to {values[-1][1]:.3f}"
    )
    for label, result in checks().items():
        print(f"check {label}: {result}")
    for label, (energy, drawn, printed) in nikolopoulos_2001().items():
        print(
            f"nikolopoulos_2001 Fig. 4 square at {energy:.3f} eV ({label}): "
            f"{drawn:.4e}, Table I {printed:.4e}, ratio {drawn / printed:.4f}"
        )


if __name__ == "__main__":
    main()
