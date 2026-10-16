# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Pathak, Sato and Ishikawa (2022): single-ionisation probabilities of Kr.

The paper prints no number of its time-dependent calculations: the probabilities are
only in **Fig. 2** (PDF p. 4, library figure ``fig_d7a5c6e5e701``), "Time evolution of
single ionization probability of Kr irradiated by a laser pulse with a wavelength of 800
nm and a peak intensity of 2 x 10^14 W/cm^2 calculated with TDHF, TD-OCCD, TD-OCCD(T),
and TD-CASSCF methods". The figure is an image embedded in the PDF (927 x 642 pixels,
JPEG, 300 dpi; its tick labels are part of the image), read here from its native pixels
(data rule 8, TC-699). The data availability statement: "The raw data supporting the
conclusion of this article will be made available by the authors, without undue
reservation"; they have not been asked.

How the image is read:

* the axes: the frame's lines and the inward tick marks, at the pixel positions written
  in :data:`X_TICKS` and :data:`Y_TICKS` (the darkness-weighted centres of the marks;
  the script finds them again in the image and checks each within :data:`TICK_CHECK`
  pixels). Time in optical cycles, 1.0 (the frame's left line) to 3.0 (its right line)
  in steps of 0.5; probability 0.0000 to 0.1200 in steps of 0.0200, linear. The zero
  tick is taken from the right edge (on the left the TDHF line covers it) and the 0.12
  tick from the left (on the right the TD-CASSCF line covers it); the other ticks are
  the same on both edges (checked).
* the curves: by colour, as the legend shows (TD-CASSCF red, TD-OCCD(T) black,
  TD-OCCD green, TDHF blue; the legend's sample lines are checked to stand in that
  order, top to bottom, :data:`LEGEND`). A pixel's weight for a colour is how far that
  channel exceeds the others (red, green, blue), or for black its darkness where the
  channels are equal within :data:`GREY`; weights below :data:`FLOOR` (JPEG noise) are
  dropped. Each curve is traced column by column from a seed row at column
  :data:`SEED_COLUMN` (:data:`SEEDS`; checked to be on the curve) up to column
  :data:`LAST_COLUMN`, the last one clear of the frame's right line (3.0 cycles is
  under that line, at column 907.4), its row in each column the weighted centre of the
  colour within :data:`WINDOW` rows of the previous one.
* the value at the end of the pulse, 3.0 cycles = 3T (the axis ends there; the
  simulation ran on to 3.3T, Table 2's note, which the figure does not show): the traced
  row at :data:`LAST_COLUMN`, 2.9939 cycles, 2.4 pixels (0.7 au) before the 3.0 tick,
  where the curves rise by less than 1e-5. As a check, a straight line through the
  traced rows of the last :data:`FIT_COLUMNS` columns (0.05 cycles), extrapolated to
  3.0 cycles, gives within 1 pixel the same (the TD-OCCD(T) line wiggles there by about
  1 pixel). The ``printed`` cell gives both rows.

One pixel is 2.16e-4 in probability and 0.0025 cycles in time; the lines are 3 to 4
pixels wide. The flat TDHF line before the pulse ionises (1.04 to 1.24 cycles, where
every curve is at zero) is drawn 1.5 pixels (3.2e-4) above the zero tick: the drawing's
own offset, not corrected, within the stated uncertainty (checked, within
:data:`ZERO_OFFSET`).
"""

from dataclasses import dataclass
from typing import Final

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, pdf
from _library import Point, write

DIRECTORY: Final = "Pathak_2022_TimeDependentOptimized"
SOURCE: Final = "pathak_2022"
PAGE: Final = 4
FIGURE: Final = "fig_d7a5c6e5e701"
SHAPE: Final = (642, 927)

# Tick positions in the native image, (column or row, value).
X_TICKS: Final = ((116.0, 1.0), (313.8, 1.5), (511.7, 2.0), (709.5, 2.5), (907.4, 3.0))
Y_TICKS: Final = (
    (14.0, 0.12),
    (106.85, 0.10),
    (199.55, 0.08),
    (292.25, 0.06),
    (385.1, 0.04),
    (477.75, 0.02),
    (570.55, 0.0),
)
TOLERANCE: Final = 0.3  # pixels: a tick's largest residual on the fitted scale
TICK_CHECK: Final = 0.25  # pixels: a written tick against the image's
# Where the ticks are measured: the rows of the bottom and top marks, the columns of
# the left and right marks, and the frame lines (the 1.0 and 3.0 ticks).
BOTTOM_MARKS: Final = (573, 578)
TOP_MARKS: Final = (7, 12)
LEFT_MARKS: Final = (118, 124)
RIGHT_MARKS: Final = (899, 905)
FRAME_ROWS: Final = (100, 500)
# Which edge each y tick is measured on (the other is covered by a curve there).
LEFT_ONLY: Final = frozenset({0.12})
RIGHT_ONLY: Final = frozenset({0.0})
MIRROR: Final = 0.25  # pixels: a tick on the left edge against its mirror

GREY: Final = 40  # the largest spread of the channels of a black (grey) pixel
FLOOR: Final = 40  # the smallest colour weight counted
LEGEND: Final = ("TD-CASSCF", "TD-OCCD(T)", "TD-OCCD", "TDHF")
COLOURS: Final = {
    "TD-CASSCF": "red",
    "TD-OCCD(T)": "black",
    "TD-OCCD": "green",
    "TDHF": "blue",
}
LEGEND_COLUMN: Final = 780
LEGEND_ROWS: Final = (250, 400)
SEED_COLUMN: Final = 880
SEEDS: Final = {"TD-CASSCF": 12.5, "TD-OCCD(T)": 25.5, "TD-OCCD": 37.0, "TDHF": 172.5}
SEED_CHECK: Final = 1.5  # pixels: a seed against the curve's centre there
LAST_COLUMN: Final = 905
FIT_COLUMNS: Final = 20
WINDOW: Final = 5
# The zero check: the TDHF line's rows before the pulse ionises.
ZERO_COLUMNS: Final = (130, 210)
ZERO_OFFSET: Final = 2.0  # pixels

DATASETS: Final = {
    "TDHF": f"Kr/{SOURCE}__fig2_tdhf_Kr.csv",
    "TD-OCCD": f"Kr/{SOURCE}__fig2_tdoccd_Kr.csv",
    "TD-OCCD(T)": f"Kr/{SOURCE}__fig2_tdoccdpt_Kr.csv",
    "TD-CASSCF": f"Kr/{SOURCE}__fig2_tdcasscf_Kr.csv",
}
INTENSITY: Final = "2e14"


def colour_image() -> np.ndarray:
    """Return the page's embedded image as RGB, native size, floats 0-255."""
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    images = [
        obj
        for obj in document[PAGE - 1].get_objects(max_depth=15)
        if isinstance(obj, pdfium.PdfImage)
    ]
    if len(images) != 1:
        msg = f"p. {PAGE}: {len(images)} images, not one"
        raise ValueError(msg)
    bitmap = images[0].get_bitmap(render=False)
    if bitmap.mode != "BGR":
        msg = f"Fig. 2: a {bitmap.mode} image, not BGR"
        raise ValueError(msg)
    rgb = np.asarray(bitmap.to_numpy(), dtype=float)[..., ::-1]
    if rgb.shape[:2] != SHAPE:
        msg = f"Fig. 2: an image of {rgb.shape[:2]}, not {SHAPE}"
        raise ValueError(msg)
    return rgb


def weights(rgb: np.ndarray) -> dict[str, np.ndarray]:
    """Return each colour's weight at every pixel (zero below :data:`FLOOR`)."""
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    high = np.maximum(np.maximum(r, g), b)
    low = np.minimum(np.minimum(r, g), b)
    found = {
        "red": r - np.maximum(g, b),
        "green": g - np.maximum(r, b),
        "blue": b - np.maximum(r, g),
        "black": np.where(high - low < GREY, 255 - high, 0.0),
    }
    return {k: np.where(v >= FLOOR, v, 0.0) for k, v in found.items()}


def _centre(profile: np.ndarray, first: int, last: int) -> float:
    """Return the weighted centre of a profile's darkness between two indices."""
    index = np.arange(first, last + 1)
    w = profile[first : last + 1] - np.median(profile)
    w = np.clip(w, 0, None)
    return float((index * w).sum() / w.sum())


def check_ticks(grey: np.ndarray) -> None:
    """Raise unless the frame and ticks are where the tick tables say."""
    dark = 255 - grey
    bottom = dark[BOTTOM_MARKS[0] : BOTTOM_MARKS[1] + 1].mean(axis=0)
    top = dark[TOP_MARKS[0] : TOP_MARKS[1] + 1].mean(axis=0)
    frame = dark[FRAME_ROWS[0] : FRAME_ROWS[1]].mean(axis=0)
    for column, value in X_TICKS:
        c = round(column)
        if value in {X_TICKS[0][1], X_TICKS[-1][1]}:  # the frame's lines
            found = [_centre(frame, c - 3, c + 3)]
        else:
            found = [_centre(bottom, c - 4, c + 4), _centre(top, c - 4, c + 4)]
        if any(abs(f - column) > TICK_CHECK for f in found):
            msg = f"Fig. 2: the {value} x tick at {found}, not {column}"
            raise ValueError(msg)
    left = dark[:, LEFT_MARKS[0] : LEFT_MARKS[1] + 1].mean(axis=1)
    right = dark[:, RIGHT_MARKS[0] : RIGHT_MARKS[1] + 1].mean(axis=1)
    for row, value in Y_TICKS:
        r = round(row)
        on_left = _centre(left, r - 4, r + 4)
        on_right = _centre(right, r - 4, r + 4)
        found = on_right if value in RIGHT_ONLY else on_left
        if abs(found - row) > TICK_CHECK:
            msg = f"Fig. 2: the {value} y tick at {found:.2f}, not {row}"
            raise ValueError(msg)
        if value not in LEFT_ONLY | RIGHT_ONLY and abs(on_left - on_right) > MIRROR:
            msg = f"Fig. 2: the {value} y tick at {on_left:.2f} and {on_right:.2f}"
            raise ValueError(msg)


def axes() -> tuple[Axis, Axis]:
    """Return the time (cycles) and probability axes, checked."""
    x_axis = Axis(
        tuple(c for c, _ in X_TICKS), tuple(v for _, v in X_TICKS), tolerance=TOLERANCE
    )
    y_axis = Axis(
        tuple(r for r, _ in Y_TICKS), tuple(v for _, v in Y_TICKS), tolerance=TOLERANCE
    )
    x_axis.check()
    y_axis.check()
    return x_axis, y_axis


def residuals(axis: Axis) -> float:
    """Return the largest distance of a tick from its axis's fitted scale, pixels."""
    slope, intercept = np.polyfit(axis.positions, axis.values, 1)
    predicted = (np.asarray(axis.values) - intercept) / slope
    return float(np.max(np.abs(predicted - np.asarray(axis.positions))))


def _row(weight: np.ndarray, column: int, near: float) -> float:
    """Return the weighted centre row of a colour in a column, near a row."""
    first, last = round(near) - WINDOW, round(near) + WINDOW
    w = weight[first : last + 1, column]
    if w.sum() == 0:
        msg = f"Fig. 2: no line at column {column} near row {near:.1f}"
        raise ValueError(msg)
    return float((np.arange(first, last + 1) * w).sum() / w.sum())


def check_legend(colour: dict[str, np.ndarray]) -> None:
    """Raise unless the legend's sample lines stand in LEGEND's order, top to bottom."""
    first, last = LEGEND_ROWS
    rows = {}
    for method, name in COLOURS.items():
        hits = np.flatnonzero(colour[name][first:last, LEGEND_COLUMN] > 100) + first
        if hits.size == 0 or np.ptp(hits) > WINDOW:
            msg = f"Fig. 2: the legend's {method} sample at rows {hits.tolist()}"
            raise ValueError(msg)
        rows[method] = float(hits.mean())
    order = tuple(sorted(rows, key=rows.__getitem__))
    if order != LEGEND:
        msg = f"Fig. 2: legend samples {order}, not {LEGEND}"
        raise ValueError(msg)


def check_zero(colour: dict[str, np.ndarray]) -> float:
    """Return the TDHF line's offset from the zero tick before the pulse, pixels."""
    zero = Y_TICKS[-1][0]
    rows = [
        _row(colour["blue"], c, zero) for c in range(ZERO_COLUMNS[0], ZERO_COLUMNS[1])
    ]
    offset = float(np.mean(rows)) - zero
    if abs(offset) > ZERO_OFFSET:
        msg = f"Fig. 2: the TDHF line before the pulse {offset:.2f} pixels off zero"
        raise ValueError(msg)
    return offset


@dataclass(frozen=True)
class Reading:
    """One curve's value at 3.0 cycles, and the trace it comes from."""

    row: float
    """The traced row at :data:`LAST_COLUMN`."""
    value: float
    fit_row: float
    """The row at 3.0 cycles of a line through the last :data:`FIT_COLUMNS` rows."""
    fit_value: float


def trace(weight: np.ndarray, seed: float) -> dict[int, float]:
    """Return a curve's centre row at each column from the seed's to the last."""
    rows = {}
    current = seed
    for column in range(SEED_COLUMN, LAST_COLUMN + 1):
        current = _row(weight, column, current)
        rows[column] = current
    return rows


def readings() -> dict[str, Reading]:
    """Return each method's probability at 3.0 cycles."""
    rgb = colour_image()
    check_ticks(rgb.mean(axis=2))
    x_axis, y_axis = axes()
    colour = weights(rgb)
    check_legend(colour)
    check_zero(colour)
    end_column = X_TICKS[-1][0]
    if abs(float(x_axis(end_column)) - 3.0) > 1e-3:
        msg = "Fig. 2: the 3.0 tick is not at 3.0 cycles"
        raise ValueError(msg)
    found = {}
    for method, name in COLOURS.items():
        seed = SEEDS[method]
        if abs(_row(colour[name], SEED_COLUMN, seed) - seed) > SEED_CHECK:
            msg = f"Fig. 2: the {method} seed is not on the {name} line"
            raise ValueError(msg)
        rows = trace(colour[name], seed)
        columns = np.arange(LAST_COLUMN - FIT_COLUMNS + 1, LAST_COLUMN + 1)
        fit = np.polyfit(columns, [rows[int(c)] for c in columns], 1)
        fitted = float(np.polyval(fit, end_column))
        if abs(fitted - rows[LAST_COLUMN]) > 1:
            msg = f"Fig. 2: {method} at 3.0 cycles {fitted:.2f}, at the last column "
            msg += f"{rows[LAST_COLUMN]:.2f}"
            raise ValueError(msg)
        row = rows[LAST_COLUMN]
        found[method] = Reading(row, float(y_axis(row)), fitted, float(y_axis(fitted)))
    values = [found[m].value for m in LEGEND]
    if values != sorted(values, reverse=True):
        msg = (
            f"Fig. 2: the curves at 3.0 cycles are not in the legend's order: {values}"
        )
        raise ValueError(msg)
    return found


def extract() -> dict[str, list[Point]]:
    """Return each method's probability at 3T, one file each."""
    found = readings()
    files: dict[str, list[Point]] = {}
    for method, relative in DATASETS.items():
        r = found[method]
        printed = (
            f"Fig. 2 {method} ({COLOURS[method]}): traced row {r.row:.2f} at column "
            f"{LAST_COLUMN} (2.9939 cycles); a line through the rows of columns "
            f"{LAST_COLUMN - FIT_COLUMNS + 1}-{LAST_COLUMN} gives row "
            f"{r.fit_row:.2f} ({r.fit_value:.4g}) at 3.0 cycles (column "
            f"{X_TICKS[-1][0]})"
        )
        files[relative] = [Point("fig2", INTENSITY, f"{r.value:.4g}", "", printed)]
    return files


def main() -> None:
    """Write the data files and print the calibration and the readings."""
    write(extract())
    x_axis, y_axis = axes()
    print(
        f"Fig. 2: tick residuals x {residuals(x_axis):.3f}, y {residuals(y_axis):.3f} "
        f"pixels; {x_axis.resolution():.5f} cycles and {y_axis.resolution():.4g} per "
        "pixel"
    )
    print(f"TDHF line before the pulse: {check_zero(weights(colour_image())):+.2f} px")
    for method, r in readings().items():
        print(
            f"{method:11s} column {LAST_COLUMN}: {r.value:.4g} (row {r.row:.2f}); "
            f"fit at 3.0 cycles: {r.fit_value:.4g} (row {r.fit_row:.2f})"
        )


if __name__ == "__main__":
    main()
