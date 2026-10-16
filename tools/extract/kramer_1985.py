# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kramer, Payne & Chen (1985): xenon's refractivity per atom at 115-117 nm.

The measured values are in a figure only. **Fig. 3** (p. 3, library figure
``fig_26507279fdf6``) plots the refractive index per xenon atom, n_3(Xe) = (n - 1)/N in
cm^3, negative, on a logarithmic axis (-1e-24 at the bottom to -1e-20 at the top)
against the vacuum wavelength, 115.0-117.2 nm, as open circles, with Mahon's
semiempirical calculation as a line. The page is a scanned image embedded in the PDF
(3168 x 4384 pixels, 400 dpi, one bit), read here from its native pixels (data rule 8,
TC-699):

* the axes: the tick marks, found in the image and checked against the positions
  written below. The scan is slightly skewed and keystoned (the frame's top line falls
  11 pixels across the plot, its bottom line 7.5), so each axis is a bilinear map of
  (column, row), fitted to the ticks of both opposite sides: x to the ten inward ticks
  115.2-117.0 nm on the bottom and on the top line, log10 |n_3| to the decade and
  minor ticks (2-9) on the left and right lines and the frame's top (-1e-20) and
  bottom (-1e-24) lines. Every tick lies on its map within 1.5 pixels.
* the markers: each open circle (13 pixels across, 7-8 inside) encloses a small white
  region; its centroid is the marker's centre. Two circles at 117.02-117.03 nm overlap
  and their insides join (an area twice a single one's); they are split in two by
  k-means on the region's pixels. The text label ("Xe 7s[3/2]1") is excluded.

One pixel is 0.0022 nm in x and 1.3 % in n_3. The ``printed`` cell of each value gives
its pixel (column, row). The per-point uncertainty (the sigma column) is the paper's:
its error limits are 'due almost entirely to the uncertainty in the value for alpha'
(p. 2), alpha = -9.0(+-0.9)e-19 cm^2, which through Eq. (1), n_3 = lambda_3 (alpha R +
beta), gives sigma = 0.1 |n_3 - lambda_3 beta| (beta = 2.6e-18 cm^2); the error bars the
figure draws at 115.23, 115.41 and 115.81 nm agree with it (checked below). It is
combined with 2 % for the reading (1.5 pixels).

The reading is checked against the numbers the text prints (p. 2): the pairs (R,
lambda_3) = (4.00, 115.41 nm) and (40.0, 116.94 nm) with alpha interpolated between
its printed values (-9.9e-19 at 115.24, -9.0e-19 at 116.23, -8.4e-19 at 117.02 nm)
give n_3 within 3 % of the figure's; the third pair, R = 9.11 at '106.23' nm
(116.23 nm; the record's notes), gives 23 % more than the figure's point at 116.23 nm
and is reported, not enforced.
"""

from collections.abc import Sequence

import numpy as np
from _figure import image
from _library import Point, write
from scipy import ndimage
from scipy.cluster.vq import kmeans2

DIRECTORY = "Kramer_1985_VacuumUltravioletRefractive"
PAGE = 3
DARK = 128
FILE = "Xe/kramer_1985__fig3_Xe.csv"
# Tick positions in the native image. x: column centres of the inward ticks, read
# along row 1075 (bottom line) and row 385 (top line), 115.2 to 117.0 nm by 0.2.
X_WAVELENGTHS = tuple(round(115.2 + 0.2 * k, 1) for k in range(10))
X_BOTTOM_ROW = 1075
X_BOTTOM = (438.5, 538.5, 637.5, 738.0, 838.0, 938.5, 1039.5, 1140.0, 1240.5, 1341.5)
X_TOP_ROW = 385
X_TOP = (435.0, 533.5, 633.5, 732.5, 832.5, 932.5, 1032.5, 1132.5, 1232.5, 1333.5)
# y: row centres of the ticks, read along column 346 (left line; from -1e-21 down: the
# top decade's minor ticks are too short there) and column 1430 (right line; from
# -9e-21 down), from top to bottom: 10, 9, 8, ..., 2 x 10^k per decade.
Y_LEFT_COLUMN = 346
Y_LEFT_FIRST = -21  # the exponent of the first (topmost) tick read
Y_LEFT = (
    552.0, 560.0, 569.5, 579.5, 591.5, 605.5, 622.5, 645.0, 676.5,
    730.0, 738.0, 747.0, 757.5, 769.5, 783.5, 801.0, 823.5, 854.5,
    908.0, 916.5, 925.5, 936.0, 948.5, 962.5, 979.5, 1001.5, 1032.5,
)  # fmt: skip
Y_RIGHT_COLUMN = 1430
Y_RIGHT = (
    372.0, 381.0, 391.0, 403.0, 417.0, 434.5, 456.5, 487.5,
    541.5, 549.5, 558.5, 569.5, 581.0, 595.0, 612.5, 635.0, 667.0,
    720.5, 729.0, 738.0, 748.0, 760.0, 774.5, 792.0, 814.5, 846.0,
    900.0, 908.5, 917.5, 928.0, 940.0, 954.5, 972.0, 994.0, 1025.0,
)  # fmt: skip
# The frame's top (-1e-20) and bottom (-1e-24) lines: (column, row of the centre).
FRAME_TOP = ((350, 374.5), (1425, 363.5))
FRAME_BOTTOM = ((350, 1086.5), (1425, 1079.0))
TOLERANCE_PIXELS = 1.5
# The plot's inside, clear of the frame lines, and the region of the text label.
INSIDE = (380, 1080, 342, 1434)  # rows, then columns
LABEL = (780, 840, 1250, 1440)
SINGLE_AREA = (30, 60)  # pixels inside one circle (40-49 found)
DOUBLE_AREA = (61, 120)  # two joined circles (87 found)
EXPECTED = 26
# The paper's constants (p. 2): beta, cm^2, and the relative uncertainty of alpha.
BETA = 2.6e-18
ALPHA_RELATIVE = 0.9 / 9.0
ALPHA = ((115.24, -9.9e-19), (116.23, -9.0e-19), (117.02, -8.4e-19))
READING = 0.02
# The text's pairs (R, lambda_3 in nm) and how closely each must match the figure.
TEXT_PAIRS = ((4.00, 115.41, 0.06), (40.0, 116.94, 0.06), (9.11, 116.23, None))
# Points drawn with an error bar (the nearest marker, nm), and how closely the drawn
# bar's upper end must match sigma = 0.1 |n_3 - lambda_3 beta| (1.09, 0.99, 1.07 found).
BARS = ((115.23, 0.15), (115.41, 0.15), (115.81, 0.15))


def _runs(line: np.ndarray, offset: int) -> list[float]:
    """Return the centres of the runs of dark pixels along a line of the image."""
    index = np.flatnonzero(line)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [(g[0] + g[-1]) / 2 + offset for g in groups]


def _check_ticks(dark: np.ndarray) -> None:
    """Raise unless the ticks are where the constants above say."""
    found = {
        "x bottom": _runs(dark[X_BOTTOM_ROW, 345:1430], 345),
        "x top": _runs(dark[X_TOP_ROW, 345:1430], 345),
        "y left": [r for r in _runs(dark[380:1082, Y_LEFT_COLUMN], 380) if r >= 550],
        "y right": [r for r in _runs(dark[368:1076, Y_RIGHT_COLUMN], 368) if r >= 370],
    }
    expected = {
        "x bottom": X_BOTTOM,
        "x top": X_TOP,
        "y left": Y_LEFT,
        "y right": Y_RIGHT,
    }
    for name, positions in expected.items():
        if tuple(found[name]) != positions:
            msg = f"Fig. 3: {name} ticks at {found[name]}, not {positions}"
            raise ValueError(msg)
    for column, row in (*FRAME_TOP, *FRAME_BOTTOM):
        near = _runs(dark[int(row) - 6 : int(row) + 7, column], int(row) - 6)
        if not any(abs(r - row) < 0.6 for r in near):
            msg = f"Fig. 3: no frame line at ({column}, {row}); found {near}"
            raise ValueError(msg)


def _tick_values(first_exponent: int, count: int, start: int = 10) -> list[float]:
    """Return log10 |value| of log-axis ticks from the top: 10, 9, ..., 2 per decade."""
    values = []
    exponent = first_exponent
    mantissas = list(range(start, 1, -1))
    while len(values) < count:
        values.extend(np.log10(m) + exponent - 1 for m in mantissas)
        mantissas = list(range(10, 1, -1))
        exponent -= 1
    return values[:count]


def _bilinear(points: Sequence[tuple[float, float, float]]) -> np.ndarray:
    """Return the coefficients of v = p0 + p1 c + p2 r + p3 c r fitted to (c, r, v)."""
    design = np.array([[1.0, c, r, c * r] for c, r, _ in points])
    values = np.array([v for _, _, v in points])
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    return coefficients


def _evaluate(p: np.ndarray, column: float, row: float) -> float:
    """Return a bilinear map's value at a pixel."""
    return float(p[0] + p[1] * column + p[2] * row + p[3] * column * row)


def _calibration() -> tuple[np.ndarray, np.ndarray]:
    """Return the x (nm) and y (log10 |n_3|) maps, checked against every tick."""
    x_points = [
        (c, X_BOTTOM_ROW, w) for c, w in zip(X_BOTTOM, X_WAVELENGTHS, strict=True)
    ]
    x_points += [(c, X_TOP_ROW, w) for c, w in zip(X_TOP, X_WAVELENGTHS, strict=True)]
    left = _tick_values(Y_LEFT_FIRST, len(Y_LEFT))
    right = _tick_values(-20, len(Y_RIGHT), start=9)
    y_points = [(Y_LEFT_COLUMN, r, v) for r, v in zip(Y_LEFT, left, strict=True)]
    y_points += [(Y_RIGHT_COLUMN, r, v) for r, v in zip(Y_RIGHT, right, strict=True)]
    y_points += [(c, r, -20.0) for c, r in FRAME_TOP]
    y_points += [(c, r, -24.0) for c, r in FRAME_BOTTOM]
    x_map, y_map = _bilinear(x_points), _bilinear(y_points)
    # The residuals in pixels: 0.2 nm is 100.3 pixels, a decade 178.
    for points, p, per_pixel in (
        (x_points, x_map, 0.2 / 100.3),
        (y_points, y_map, 1 / 178),
    ):
        worst = max(abs(_evaluate(p, c, r) - v) / per_pixel for c, r, v in points)
        if worst > TOLERANCE_PIXELS:
            msg = f"Fig. 3: a tick lies {worst:.2f} pixels off its axis"
            raise ValueError(msg)
    return x_map, y_map


def _markers(dark: np.ndarray) -> list[tuple[float, float]]:
    """Return the centres (column, row) of the open circles, left to right."""
    top, bottom, left, right = INSIDE
    white = ~dark[top:bottom, left:right]
    labels, count = ndimage.label(white)
    centres = []
    for index in range(1, count + 1):
        rows, columns = np.nonzero(labels == index)
        area = rows.size
        row, column = rows.mean() + top, columns.mean() + left
        if LABEL[0] < row < LABEL[1] and LABEL[2] < column < LABEL[3]:
            continue
        if SINGLE_AREA[0] <= area <= SINGLE_AREA[1]:
            centres.append((column, row))
        elif DOUBLE_AREA[0] <= area <= DOUBLE_AREA[1]:
            pixels = np.column_stack([columns + left, rows + top]).astype(float)
            pair, _ = kmeans2(pixels, 2, seed=1, minit="++")
            centres.extend((float(c), float(r)) for c, r in pair)
    if len(centres) != EXPECTED:
        msg = f"Fig. 3: {len(centres)} markers, not {EXPECTED}"
        raise ValueError(msg)
    return sorted(centres)


def _sigma(wavelength_nm: float, value: float) -> float:
    """Return the paper's uncertainty of n_3 (alpha's, by Eq. (1)) with the reading."""
    stated = ALPHA_RELATIVE * abs(value - wavelength_nm * 1e-7 * BETA)
    return float(np.hypot(stated, READING * abs(value)))


def _check_text(points: list[tuple[float, float]]) -> None:
    """Raise unless the text's (R, lambda_3) pairs reproduce the figure's values."""
    wavelengths = np.array([w for w, _ in points])
    log_values = np.log10([-v for _, v in points])
    for ratio, wavelength, tolerance in TEXT_PAIRS:
        alpha = np.interp(wavelength, [w for w, _ in ALPHA], [a for _, a in ALPHA])
        printed = wavelength * 1e-7 * (alpha * ratio + BETA)
        read = -(10 ** np.interp(wavelength, wavelengths, log_values))
        if tolerance is not None and abs(printed / read - 1) > tolerance:
            msg = f"Fig. 3: {read:.3e} at {wavelength} nm, Eq. (1) {printed:.3e}"
            raise ValueError(msg)
        print(
            f"R = {ratio} at {wavelength} nm: Eq. (1) {printed:.3e}, Fig. 3 {read:.3e}"
        )


def _check_bars(
    dark: np.ndarray,
    y_map: np.ndarray,
    markers: list[tuple[float, float, float, float]],
) -> None:
    """Raise unless the drawn error bars' upper ends match the stated uncertainty."""
    for wavelength, tolerance in BARS:
        column, row, value, near = min(markers, key=lambda m: abs(m[3] - wavelength))
        if abs(near - wavelength) > 0.01:
            msg = f"Fig. 3: no marker at {wavelength} nm"
            raise ValueError(msg)
        c, end = round(column), round(row) - 7  # above the circle's outline
        while dark[end - 1, c - 1 : c + 2].any():
            end -= 1
        drawn = 10 ** _evaluate(y_map, column, end) - abs(value)
        stated = ALPHA_RELATIVE * abs(value - near * 1e-7 * BETA)
        if abs(drawn / stated - 1) > tolerance:
            msg = (
                f"Fig. 3: the bar at {near} nm is {drawn:.2e}, alpha gives {stated:.2e}"
            )
            raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return Fig. 3's measured n_3(Xe), with the paper's uncertainty."""
    grey = image(DIRECTORY, PAGE)
    if grey.shape != (4384, 3168):
        msg = f"p. 3: an image of {grey.shape}, not 4384 x 3168"
        raise ValueError(msg)
    dark = grey < DARK
    _check_ticks(dark)
    x_map, y_map = _calibration()
    read = []
    for column, row in _markers(dark):
        wavelength = _evaluate(x_map, column, row)
        value = -(10 ** _evaluate(y_map, column, row))
        read.append((column, row, value, wavelength))
    _check_bars(dark, y_map, read)
    _check_text([(w, v) for _, _, v, w in read])
    points = [
        Point(
            f"{wavelength:.3f}",
            f"{wavelength:.3f}",
            f"{value:.3e}",
            f"{_sigma(wavelength, value):.2e}",
            f"pixel ({column:.1f}, {row:.1f})",
        )
        for column, row, value, wavelength in read
    ]
    return {FILE: points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
