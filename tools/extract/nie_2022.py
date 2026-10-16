# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Nie et al. (2022): measured plasma phase shifts of He, Ar, Kr, Xe and H2.

The paper prints no measured phase shift as a number; they are in figures only
(data rule 8, TC-699). The measured points with their error bars are read from three
figures; the theory curves (ADK, PPT, MO-PPT and the fitted mixtures) are not taken.

**Fig. 3(a)** (p. 6, library figure ``fig_b4bda677c9ff``): He, 800 nm, 55 fs, linear
polarisation, eight orange squares ('Expt.-LP'). Panels (a)-(c) are one image embedded
in the PDF (929 x 918 pixels, lossless), read here from its native pixels:

* the axes: both are logarithmic. The tick marks are the runs of dark grey pixels
  (``grey_ink``: coloured curves are not ink) along a row or column just inside the
  frame, on both opposite sides (``RasterAxis.bands``; a tick hidden on one side by
  a curve counts from the other), less lines that cross the whole plot. The
  intensity's 39 ticks run from 2.2 to 9.8e14 W/cm^2 in steps of 0.2 (long ticks at
  the labelled 3 to 9), the phase's 35 from 2e-3 to 9 rad (long ticks at the
  labelled decades). The count must match, every tick must lie on the fitted scale
  (``TICK_TOLERANCE``), and the long ticks (``RasterAxis.deep``) must be the labelled
  values.
* the markers: the connected regions of orange pixels (``_orange``) in the panel,
  less the legend's (``FIG3_LEGEND``). In each, the marker is the 9 x 9 pixel window
  of the largest orange weight, and its centre the weighted centroid of the colour in
  that window (``_locate``).
* the error bars: from the marker's edge outward along the bar's columns while the
  bar is there, the outermost cap (a short horizontal run of orange either side of
  the bar and no wider, ``_bar_end``); the end is the mean of the cap's rows. A bar
  no longer than the marker is hidden by it, and its ``sigma`` is left empty.
* the pressure printed beside each marker (1500 to 3 Torr), read on the image and
  written in ``FIG3_PRESSURES_TORR`` in the order of the markers; the extremes are
  checked against the text (p. 6).

**Fig. 4** (p. 8, library figure ``fig_5947c22a506f``): Ar, Kr and Xe, 800 nm, 55 fs,
linear (blue squares, 'Expt.-LP') and circular polarisation (red circles, 'Expt.-CP').
One image (1153 x 543 pixels, lossless), three panels:

* the axes: as for Fig. 3, the ticks found on both sides of each panel's frame
  (``FIG4_FRAMES``), the long one under the label '10^14' being 1e14 W/cm^2 and the
  short ones its neighbours k x 10^n (``FIG4_X_VALUES``); the phase axis, 1e-3 to
  1e3 rad, has 53 ticks inside the frame. Lines that run the plot's full height (the
  dashed Keldysh lines) are not ticks and are dropped.
* the markers: blue and red overlap the theory curves of the same colour, so each
  marker is seeded at a pixel found by inspecting the image (``FIG4_SEEDS``) and
  located by rule: the 6 x 6 window of the largest blue (``_blue``) or red (``_red``)
  weight within 2 pixels of the seed, and the weighted centroid of the colour in it.
  The legend's markers are not seeded.
* the error bars: the bars too cross the curves, which a cap-finding rule mistakes
  for caps, so each cap that shows beyond its marker is seeded at a row by
  inspection (``FIG4_CAPS``; at this resolution a cap touching the marker cannot be
  told from the marker's outline, and is not read) and located by rule
  (``_cap_row``): the rows within a pixel of the seed coloured just outside the bar's
  two columns on both sides, weighted by that colour; the bar itself must show
  between the marker and the cap.

**Fig. 6** (p. 9, library figure ``fig_4b9d5a21931f``): H2 at 800, 400 and 267 nm,
linear polarisation ('Exp. data'). A vector drawing, read from the PDF (coordinates of
the figure's form XObject, points):

* the axes: the grey tick segments of each panel's bottom and left frame; the long
  ones are matched to the labels '10' + exponent of the text layer, the short ones are
  the k x 10^n between them (k by rank). Every tick is checked against the fitted
  scale.
* the markers: the blue (0, 114, 189) stroked squares, 1.35 pt a side (the legend's
  are 0.9 pt); their centres.
* the error bars: the blue vertical segments drawn from each marker's centre, one up
  and one down; their far ends.

Every value's ``printed`` cell gives its marker's pixel (column, row) or point (x, y)
and its bar's ends (rows or y). ``sigma`` is half the bar's full length in rad (with
one end hidden in the marker, the visible half's length). On the logarithmic axis a bar
symmetric in rad is longer below the marker than above; :func:`bar_asymmetry` gives the
ratio of the halves in rad, recorded in the notes.

The record's parameters (pulse durations, spot sizes, probe wavelengths) are checked
here against Table 1 (``tbl_01fff8291f09``) by :func:`check_parameters`.
"""

import ctypes
import math
import tomllib
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, pdf
from _library import DATA, LIBRARY, Point, table, write
from pypdfium2 import raw
from scipy import ndimage

DIRECTORY = "Nie_2022_CrossPolarizedCommon"
SOURCE = "nie_2022"
TABLE_1 = "tbl_01fff8291f09"

# The weight of a pixel's colour, 0 to 1: how far it is from grey towards the markers'
# colour. A threshold of ON marks a pixel of the colour; CAP_ON, the faint edge of a
# cap (the antialiased end pixels).
ON = 0.35
CAP_ON = 0.1
INK = 75  # 255 minus the mean of the channels: a dark (tick) pixel
CHROMA = 60  # the spread of the channels above which a pixel is coloured, not grey


def grey_ink(rgb: np.ndarray) -> np.ndarray:
    """Return 255 minus the channels' mean for grey pixels, 0 for coloured ones.

    The axes are drawn in black and grey; a coloured curve crossing a tick band is
    not ink.
    """
    chroma = rgb.max(axis=2) - rgb.min(axis=2)
    return np.where(chroma < CHROMA, 255 - rgb.mean(axis=2), 0.0)


def _orange(rgb: np.ndarray) -> np.ndarray:
    """Fig. 3's markers, RGB (223, 87, 26): red over blue."""
    return np.clip((rgb[..., 0] - rgb[..., 2] - 20) / 150, 0, 1)


def _blue(rgb: np.ndarray) -> np.ndarray:
    """Fig. 4's linear-polarisation markers, RGB about (57, 101, 184)."""
    return np.clip((rgb[..., 2] - rgb[..., 0] - 20) / 100, 0, 1)


def _red(rgb: np.ndarray) -> np.ndarray:
    """Fig. 4's circular-polarisation markers, RGB about (187, 0, 0)."""
    excess = np.minimum(rgb[..., 0] - rgb[..., 1], rgb[..., 0] - rgb[..., 2])
    return np.clip((excess - 20) / 150, 0, 1)


def colour_image(page: int, index: int) -> np.ndarray:
    """Return a page's (1-based) embedded image as RGB levels 0-255, native size."""
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    images = [
        obj
        for obj in document[page - 1].get_objects(max_depth=15)
        if isinstance(obj, pdfium.PdfImage)
    ]
    bitmap = images[index].get_bitmap(render=False)
    if bitmap.mode not in {"BGR", "BGRA", "BGRX"}:
        msg = f"p. {page}: an image in {bitmap.mode}, not BGR"
        raise ValueError(msg)
    return np.asarray(bitmap.to_numpy()[..., 2::-1], dtype=float)


# --- the raster axes -----------------------------------------------------------------


def _centres(line: np.ndarray) -> list[float]:
    """Return the ink-weighted centres of the runs of dark pixels along a line."""
    index = np.flatnonzero(line > INK)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [float((g * line[g]).sum() / line[g].sum()) for g in groups]


@dataclass(frozen=True)
class RasterAxis:
    """Where an axis's ticks are in an image, and what they print.

    Attributes
    ----------
    along : str
        "x" for a horizontal axis (ticks are columns), "y" for a vertical one.
    bands : tuple of int
        The row (x) or column (y) through the short ticks, on each side of the frame.
    deep : tuple of int
        The row or column through only the long ticks, on each side.
    span : tuple of int
        The range of columns (x) or rows (y) searched: inside the frame, clear of it.
    interior : tuple of int
        The rows (x) or columns (y) of the plot's inside, to drop full-length lines.
    values : tuple of float
        The ticks' values, in increasing order of value.
    long : tuple of float
        The values of the long (labelled) ticks.
    """

    along: str
    bands: tuple[int, int]
    deep: tuple[int, int]
    span: tuple[int, int]
    interior: tuple[int, int]
    values: tuple[float, ...]
    long: tuple[float, ...]

    def _line(self, ink: np.ndarray, at: int) -> np.ndarray:
        first, last = self.span
        return ink[at, first:last] if self.along == "x" else ink[first:last, at]

    def axis(self, ink: np.ndarray, name: str) -> Axis:
        """Return the calibrated axis; raise unless every tick is where it should be."""
        first, _ = self.span
        one, two = (
            [first + c for c in _centres(self._line(ink, at))] for at in self.bands
        )
        # A tick on both sides is their mean; on one side only (a curve drawn over the
        # other), that side's.
        found = [(a + b) / 2 for a in one for b in two if abs(a - b) <= 1.0]
        found += [a for a in one + two if all(abs(a - f) > 1.0 for f in found)]
        found.sort()
        low, high = self.interior
        kept = []
        for position in found:
            p = round(position)
            through = ink[low:high, p] if self.along == "x" else ink[p, low:high]
            if (through > INK).mean() < 0.3:  # not a line across the plot
                kept.append(position)
        if len(kept) != len(self.values):
            msg = f"{name}: {len(kept)} ticks found, not {len(self.values)}"
            raise ValueError(msg)
        # Rows grow downwards: the largest value is at the smallest row.
        values = self.values if self.along == "x" else self.values[::-1]
        for position, value in zip(kept, values, strict=True):
            p = round(position)
            deep = [
                ink[at, p - 1 : p + 2] if self.along == "x" else ink[p - 1 : p + 2, at]
                for at in self.deep
            ]
            is_long = all(float(d.max()) > INK for d in deep)
            if is_long != any(math.isclose(value, v) for v in self.long):
                kind = "long" if is_long else "short"
                msg = f"{name}: the tick of {value:g} at {position:.2f} is {kind}"
                raise ValueError(msg)
        calibrated = Axis(tuple(kept), values, log=True, tolerance=TICK_TOLERANCE)
        calibrated.check()
        return calibrated


# The largest residual of a tick from the fitted scale, pixels.
TICK_TOLERANCE = 0.6


def worst_residual(axis: Axis) -> float:
    """Return the largest distance of a tick from the fitted scale (pixels or pt)."""
    slope, intercept = np.polyfit(axis.positions, np.log10(axis.values), 1)
    predicted = (np.log10(axis.values) - intercept) / slope
    return float(np.max(np.abs(predicted - np.asarray(axis.positions))))


# --- the raster markers and their error bars ------------------------------------------


def _window_sums(weight: np.ndarray, size: int) -> np.ndarray:
    """Return the sum of the weight over each size x size window, by its top left."""
    total = np.pad(weight.cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    return (
        total[size:, size:]
        - total[:-size, size:]
        - total[size:, :-size]
        + total[:-size, :-size]
    )


def _locate(
    weight: np.ndarray,
    sums: np.ndarray,
    size: int,
    rows: tuple[int, int],
    columns: tuple[int, int],
) -> tuple[float, float]:
    """Return a marker's centre (column, row) from windows with top left in range.

    The heaviest size x size window (equal maxima: the first); the centre is the
    weighted centroid of the colour inside it.
    """
    r0, r1 = rows
    c0, c1 = columns
    block = sums[r0 : r1 + 1, c0 : c1 + 1]
    rr, cc = np.unravel_index(int(np.argmax(block)), block.shape)
    top, left = r0 + int(rr), c0 + int(cc)
    inside = weight[top : top + size, left : left + size]
    grid_r, grid_c = np.mgrid[top : top + size, left : left + size]
    total = inside.sum()
    return float((grid_c * inside).sum() / total), float(
        (grid_r * inside).sum() / total
    )


@dataclass(frozen=True)
class Bars:
    """How a figure's error bars look, pixels.

    Attributes
    ----------
    half : float
        Half the marker's height: the bar is searched beyond it.
    cap_in : int
        How far beyond the bar's columns a cap reaches at least.
    cap_out : int
        How far beyond them a cap never reaches (a crossing line does).
    longest : int
        The longest half-bar searched.
    """

    half: float
    cap_in: int
    cap_out: int
    longest: int = 120


def _bar_end(
    weight: np.ndarray, centre: tuple[float, float], bars: Bars, direction: int
) -> float | None:
    """Return the row of the end of a marker's half-bar (up -1, down +1), or None.

    From the marker's edge outward along the bar's columns (within 1 pixel of the
    centre) while the bar is there, the outermost
    run of rows with a cap: the colour ``cap_in`` beyond the bar on both sides and not
    ``cap_out`` beyond it. The end is the mean of that run's rows.
    """
    column, row = centre
    cols = [
        c
        for c in range(math.floor(column) - 1, math.ceil(column) + 2)
        if abs(c - column) <= 1
    ]
    start = row + direction * bars.half
    r = math.ceil(start) if direction > 0 else math.floor(start)
    caps: list[int] = []
    while 0 <= r < weight.shape[0] and abs(r - row) <= bars.longest:
        bar = [c for c in cols if weight[r, c] >= ON]
        if not bar:
            break
        left, right = min(bar), max(bar)
        inside = (
            weight[r, left - bars.cap_in] >= CAP_ON
            and weight[r, right + bars.cap_in] >= CAP_ON
        )
        outside = (
            weight[r, left - bars.cap_out] >= CAP_ON
            or weight[r, right + bars.cap_out] >= CAP_ON
        )
        if inside and not outside:
            caps.append(r)
        r += direction
    if not caps:
        return None
    run = [caps[-1]]
    for c in reversed(caps[:-1]):
        if abs(c - run[-1]) != 1:
            break
        run.append(c)
    return float(np.mean(run))


@dataclass(frozen=True)
class Marker:
    """A marker read from an image or a drawing: centre, bar ends (None: hidden)."""

    x: float
    y: float
    top: float | None
    bottom: float | None
    unit: str  # "px" or "pt"


@dataclass(frozen=True)
class Reading:
    """A dataset's markers, the axes that read them, and their extra columns."""

    markers: list[Marker]
    x_axis: Axis
    y_axis: Axis
    extras: list[dict[str, str]] | None = None


def _point(
    row: str,
    marker: Marker,
    x_axis: Axis,
    y_axis: Axis,
    extra: dict[str, str] | None = None,
) -> Point:
    """Return a marker's intensity, phase and sigma (half the bar's length), rad."""
    intensity = float(x_axis(marker.x))
    phase = float(y_axis(marker.y))
    ends = [e for e in (marker.top, marker.bottom) if e is not None]
    digits = 1 if marker.unit == "px" else 3
    if marker.unit == "px":
        where = f"pixel ({marker.x:.1f}, {marker.y:.1f})"
    else:
        where = f"marker at ({marker.x:.3f}, {marker.y:.3f}) pt"
    if len(ends) == 2:
        upper, lower = (float(y_axis(e)) for e in ends)
        sigma = f"{(max(upper, lower) - min(upper, lower)) / 2:.3e}"
        printed = f"{where}; bar ends {ends[0]:.{digits}f} and {ends[1]:.{digits}f}"
    elif len(ends) == 1:
        end = float(y_axis(ends[0]))
        sigma = f"{abs(end - phase):.3e}"
        printed = f"{where}; one bar end, {ends[0]:.{digits}f} (the other hidden)"
    else:
        sigma = ""
        printed = f"{where}; no bar beyond the marker"
    return Point(row, f"{intensity:.4e}", f"{phase:.3e}", sigma, printed, extra or {})


# --- Fig. 3(a) ----------------------------------------------------------------------

FIG3_PAGE, FIG3_INDEX, FIG3_SHAPE = 6, 1, (918, 929)
# The frame of panel (a): columns 100-101 and 897-898, rows 29-30 and 368-369.
FIG3_X = RasterAxis(
    along="x",
    bands=(32, 366),
    deep=(37, 361),
    span=(104, 895),
    interior=(45, 355),
    values=tuple(round(2.2 + 0.2 * k, 1) * 1e14 for k in range(39)),
    long=tuple(k * 1e14 for k in range(3, 10)),
)
FIG3_Y = RasterAxis(
    along="y",
    bands=(103, 895),
    deep=(108, 890),
    span=(33, 366),
    interior=(115, 885),
    values=tuple(k * 10.0**n for n in range(-3, 1) for k in range(1, 10))[1:],
    long=(1e-2, 1e-1, 1.0),
)
# Panel (a)'s inside, and its legend box (rows, columns), whose marker is not a point.
FIG3_PANEL = ((31, 367), (102, 896))
FIG3_LEGEND = ((284, 352), (727, 862))
FIG3_MARKER = 9
FIG3_BARS = Bars(half=5.5, cap_in=1, cap_out=4)
# Printed beside the markers, in their order of intensity (read on the image).
FIG3_PRESSURES_TORR = ("1500", "800", "400", "100", "50", "12.5", "6", "3")


def _fig3() -> Reading:
    """Return Fig. 3(a)'s markers, with the pressures printed beside them."""
    rgb = colour_image(FIG3_PAGE, FIG3_INDEX)
    if rgb.shape[:2] != FIG3_SHAPE:
        msg = f"Fig. 3: an image of {rgb.shape[:2]}, not {FIG3_SHAPE}"
        raise ValueError(msg)
    ink = grey_ink(rgb)
    x_axis = FIG3_X.axis(ink, "Fig. 3(a) x")
    y_axis = FIG3_Y.axis(ink, "Fig. 3(a) y")
    weight = _orange(rgb)
    (r0, r1), (c0, c1) = FIG3_PANEL
    mask = np.zeros(weight.shape, dtype=bool)
    mask[r0:r1, c0:c1] = weight[r0:r1, c0:c1] >= ON
    (l0, l1), (k0, k1) = FIG3_LEGEND
    mask[l0:l1, k0:k1] = False
    labels, count = ndimage.label(mask, structure=np.ones((3, 3)))
    sums = _window_sums(weight, FIG3_MARKER)
    markers = []
    for found in ndimage.find_objects(labels)[:count]:
        if found is None:
            continue
        rows, cols = found
        # Windows that lie within the region, give or take 2 pixels.
        top_lefts_r = (rows.start - 2, max(rows.start, rows.stop - FIG3_MARKER) + 2)
        top_lefts_c = (cols.start - 2, max(cols.start, cols.stop - FIG3_MARKER) + 2)
        centre = _locate(weight, sums, FIG3_MARKER, top_lefts_r, top_lefts_c)
        top = _bar_end(weight, centre, FIG3_BARS, -1)
        bottom = _bar_end(weight, centre, FIG3_BARS, +1)
        markers.append(Marker(*centre, top, bottom, "px"))
    markers.sort(key=lambda m: m.x)
    if len(markers) != len(FIG3_PRESSURES_TORR):
        msg = f"Fig. 3(a): {len(markers)} markers, not {len(FIG3_PRESSURES_TORR)}"
        raise ValueError(msg)
    _check_pressures_in_text()
    extras = [{"pressure": torr} for torr in FIG3_PRESSURES_TORR]
    return Reading(markers, x_axis, y_axis, extras)


def _check_pressures_in_text() -> None:
    """Raise unless p. 6 states the extreme pressures of Fig. 3(a)."""
    text = (LIBRARY / DIRECTORY / "paper.md").read_text(encoding="utf-8")
    for phrase in (
        f"we used a gas pressure of {FIG3_PRESSURES_TORR[-1]} Torr",
        f"we used a gas pressure of {FIG3_PRESSURES_TORR[0]} Torr",
    ):
        if phrase not in text:
            msg = f"Fig. 3(a): the text does not say {phrase!r}"
            raise ValueError(msg)


# --- Fig. 4 -------------------------------------------------------------------------

FIG4_PAGE, FIG4_INDEX, FIG4_SHAPE = 8, 1, (543, 1153)
# Each panel's frame columns (left, right); the frame's rows are 18 and 484 in all.
FIG4_FRAMES = {"Ar": (61, 381), "Kr": (445, 764), "Xe": (828, 1147)}
# The intensities of the ticks inside each frame, W/cm^2 (the frame's edges, 4e13 and
# 6e14 for Ar, 4e13 and 5e14 for Kr, 2e13 and 4e14 for Xe, are not counted).
FIG4_X_VALUES = {
    "Ar": tuple(k * 1e13 for k in range(5, 10)) + tuple(k * 1e14 for k in range(1, 6)),
    "Kr": tuple(k * 1e13 for k in range(5, 10)) + tuple(k * 1e14 for k in range(1, 5)),
    "Xe": tuple(k * 1e13 for k in range(3, 10)) + tuple(k * 1e14 for k in range(1, 4)),
}
FIG4_Y_VALUES = tuple(k * 10.0**n for n in range(-3, 3) for k in range(1, 10))[1:]
FIG4_MARKER = 6
FIG4_SEARCH = 2
# Seeds (column, row): a pixel on each marker, found by inspecting the image, in order
# of intensity. The markers' pixels are then located by rule (the module docstring).
FIG4_SEEDS = {
    ("Ar", "LP"): (
        (98, 394),
        (132, 333),
        (158, 272),
        (180, 232),
        (214, 172),
        (240, 141),
        (262, 110),
        (280, 99),
        (296, 88),
        (322, 68),
    ),
    ("Ar", "CP"): (
        (180, 335),
        (214, 264),
        (240, 223),
        (262, 181),
        (278, 159),
        (294, 134),
        (322, 104),
        (344, 85),
    ),
    ("Kr", "LP"): (
        (485, 276),
        (520, 229),
        (548, 187),
        (572, 158),
        (608, 122),
        (636, 94),
        (659, 77),
        (679, 68),
        (696, 57),
        (711, 49),
    ),
    ("Kr", "CP"): (
        (520, 331),
        (549, 282),
        (572, 245),
        (607, 191),
        (636, 149),
        (660, 123),
        (679, 100),
        (696, 87),
        (711, 77),
        (724, 69),
    ),
    ("Xe", "LP"): (
        (900, 244),
        (943, 171),
        (974, 126),
        (997, 106),
        (1017, 85),
        (1047, 62),
        (1071, 49),
        (1090, 42),
        (1107, 36),
    ),
    ("Xe", "CP"): (
        (900, 345),
        (943, 258),
        (974, 207),
        (998, 166),
        (1017, 133),
        (1048, 98),
        (1072, 77),
        (1091, 62),
        (1107, 51),
        (1122, 46),
    ),
}
WEIGHTS: dict[str, Callable[[np.ndarray], np.ndarray]] = {"LP": _blue, "CP": _red}


def _fig4_axes(ink: np.ndarray, gas: str) -> tuple[Axis, Axis]:
    left, right = FIG4_FRAMES[gas]
    x = RasterAxis(
        along="x",
        bands=(20, 482),
        deep=(22, 480),
        span=(left + 3, right - 2),
        interior=(30, 470),
        values=FIG4_X_VALUES[gas],
        long=(1e14,),
    )
    y = RasterAxis(
        along="y",
        bands=(left + 2, right - 2),
        deep=(left + 4, right - 4),
        span=(21, 482),
        interior=(left + 8, right - 8),
        values=FIG4_Y_VALUES,
        long=(1e-2, 1e-1, 1.0, 10.0, 100.0),
    )
    return x.axis(ink, f"Fig. 4 {gas} x"), y.axis(ink, f"Fig. 4 {gas} y")


# The caps of the error bars that show beyond the markers: a row on each cap, found by
# inspecting the image (top, bottom; None where no bar shows beyond the marker), for
# the markers in the order of FIG4_SEEDS. Markers not listed show no bar. The bars cross
# the theory curves of the same colour, which a cap-finding rule mistakes for caps; and
# a cap that touches the marker (a bar no longer than the marker) is not read.
FIG4_CAPS: dict[tuple[str, str], dict[int, tuple[int | None, int | None]]] = {
    ("Ar", "LP"): {
        1: (388, 402),
        2: (328, 338),
        3: (266, 279),
        4: (None, 236),
        5: (167, 178),
        6: (137, 146),
        7: (105, 115),
    },
    ("Ar", "CP"): {1: (326, 349), 2: (258, 272), 7: (99, 109)},
    ("Kr", "LP"): {1: (268, 285), 2: (224, 235)},
    ("Kr", "CP"): {1: (320, 346)},
    ("Xe", "LP"): {1: (238, 252), 2: (162, 184), 3: (121, 133)},
    ("Xe", "CP"): {1: (338, 354)},
}
# Half the height of a marker, pixels: a bar shows if it is coloured beyond this.
FIG4_HALF = 3.0


def _cap_row(
    weight: np.ndarray, centre: tuple[float, float], seed: int, name: str
) -> float:
    """Return the row of a bar's cap near a seed row, by rule.

    The bar spans the columns floor(x) and ceil(x) of the marker's centre; a cap row is
    one, within a pixel of the seed, coloured (``CAP_ON``) in the columns just outside
    them on both sides. The cap is the mean of its rows weighted by that colour. The bar
    must be coloured (``ON``) in a row between the marker's edge and the cap.
    """
    column, row = centre
    left, right = math.floor(column), math.ceil(column)
    if left == right:
        right += 1
    rows, weights = [], []
    for r in (seed - 1, seed, seed + 1):
        side = (weight[r, left - 1], weight[r, right + 1])
        if min(side) >= CAP_ON:
            rows.append(r)
            weights.append(sum(side))
    if not rows:
        msg = f"{name}: no cap within a pixel of row {seed}"
        raise ValueError(msg)
    cap = float(np.average(rows, weights=weights))
    direction = 1 if seed > row else -1
    edge = row + direction * FIG4_HALF
    between = [
        r
        for r in range(round(row), seed, direction)
        if (r - edge) * direction > 0 and (r - cap) * direction < -0.5
    ]
    if not any(max(weight[r, left], weight[r, right]) >= ON for r in between):
        msg = f"{name}: no bar between the marker and the cap at row {seed}"
        raise ValueError(msg)
    return cap


def _fig4() -> dict[tuple[str, str], Reading]:
    """Return Fig. 4's markers, panel by panel and polarisation by polarisation."""
    rgb = colour_image(FIG4_PAGE, FIG4_INDEX)
    if rgb.shape[:2] != FIG4_SHAPE:
        msg = f"Fig. 4: an image of {rgb.shape[:2]}, not {FIG4_SHAPE}"
        raise ValueError(msg)
    ink = grey_ink(rgb)
    out = {}
    for gas in FIG4_FRAMES:
        x_axis, y_axis = _fig4_axes(ink, gas)
        for polarisation, weigh in WEIGHTS.items():
            weight = weigh(rgb)
            sums = _window_sums(weight, FIG4_MARKER)
            half = (FIG4_MARKER - 1) / 2
            caps = FIG4_CAPS[(gas, polarisation)]
            markers = []
            seeds = FIG4_SEEDS[(gas, polarisation)]
            for k, (column, row) in enumerate(seeds, start=1):
                r0 = math.floor(row - half) - FIG4_SEARCH
                c0 = math.floor(column - half) - FIG4_SEARCH
                rows = (r0, r0 + 2 * FIG4_SEARCH + 1)
                cols = (c0, c0 + 2 * FIG4_SEARCH + 1)
                centre = _locate(weight, sums, FIG4_MARKER, rows, cols)
                name = f"Fig. 4 {gas} {polarisation} p{k:02d}"
                top_seed, bottom_seed = caps.get(k, (None, None))
                top = (
                    None
                    if top_seed is None
                    else _cap_row(weight, centre, top_seed, name)
                )
                bottom = (
                    None
                    if bottom_seed is None
                    else _cap_row(weight, centre, bottom_seed, name)
                )
                markers.append(Marker(*centre, top, bottom, "px"))
            out[(gas, polarisation)] = Reading(markers, x_axis, y_axis)
    return out


# --- Fig. 6 (vector) ----------------------------------------------------------------

FIG6_PAGE = 9
# The figure's form XObject, in its own coordinates (pt): the region of the drawing.
FIG6_REGION = (0.0, 0.0, 380.0, 200.0)
# Each panel's frame (left, right, bottom, top), pt.
FIG6_FRAMES = {
    "800nm": (19.066, 124.199, 21.281, 174.702),
    "400nm": (145.192, 250.224, 21.281, 174.700),
    "267nm": (271.217, 376.350, 21.281, 174.700),
}
GREY = (78, 78, 80)
BLUE = (0, 114, 189)
SQUARE = (1.2, 1.5)  # the side of a data marker, pt (the legend's is 0.9)
MERGE = 0.02  # ticks closer than this, pt, are one tick drawn twice
LONG = 1.2  # a tick longer than this, pt, is a long (labelled) one


@dataclass(frozen=True)
class Path:
    """A vector path inside the figure's form: colour, mode and vertices (pt)."""

    stroke: tuple[int, int, int]
    stroked: bool
    vertices: np.ndarray


@dataclass(frozen=True)
class Label:
    """A tick label '10' with its exponent: the exponent and the label's centre (pt)."""

    exponent: int
    x: float
    y: float


def _text(obj: pdfium.PdfObject, page: pdfium.PdfTextPage) -> str:
    size = raw.FPDFTextObj_GetText(obj.raw, page.raw, None, 0)
    buffer = (ctypes.c_ushort * size)()
    raw.FPDFTextObj_GetText(obj.raw, page.raw, buffer, size * 2)
    return bytes(buffer).decode("utf-16-le").replace("\x00", "").strip()


def _fig6_objects() -> tuple[list[Path], list[Label]]:
    """Return the paths and the tick labels nested in Fig. 6's form."""
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    page = document[FIG6_PAGE - 1]
    textpage = page.get_textpage()
    left, bottom, right, top = FIG6_REGION
    paths, texts = [], []
    for obj in page.get_objects(max_depth=15):
        if obj.level < 1:  # page-level, not in the figure's form
            continue
        west, south, east, north = obj.get_bounds()
        if not (west >= left and east <= right and south >= bottom and north <= top):
            continue
        kind = raw.FPDFPageObj_GetType(obj.raw)
        if kind == raw.FPDF_PAGEOBJ_PATH:
            r, g, b, a = (ctypes.c_uint() for _ in range(4))
            raw.FPDFPageObj_GetStrokeColor(
                obj.raw, *(ctypes.byref(v) for v in (r, g, b, a))
            )
            fill, stroke = ctypes.c_int(), ctypes.c_int()
            raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill), ctypes.byref(stroke))
            paths.append(
                Path((r.value, g.value, b.value), bool(stroke.value), _vertices(obj))
            )
        elif kind == raw.FPDF_PAGEOBJ_TEXT:
            texts.append((_text(obj, textpage), west, south, east, north))
    labels = []
    for text, west, south, east, north in texts:
        if text != "10":
            continue
        exponent = [
            t
            for t in texts
            if abs(t[1] - east) < 1.5
            and south < t[2] < north
            and t[0].lstrip("-").isdigit()
        ]
        if len(exponent) != 1:
            msg = f"Fig. 6: the label '10' at ({west:.1f}, {south:.1f}) has no exponent"
            raise ValueError(msg)
        e = exponent[0]
        labels.append(Label(int(e[0]), (west + e[3]) / 2, (south + e[4]) / 2))
    return paths, labels


def _vertices(obj: pdfium.PdfObject) -> np.ndarray:
    """Return a path object's vertices in the coordinates of its form."""
    matrix = obj.get_matrix()
    points = []
    for index in range(raw.FPDFPath_CountSegments(obj.raw)):
        segment = raw.FPDFPath_GetPathSegment(obj.raw, index)
        x, y = ctypes.c_float(), ctypes.c_float()
        raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
        points.append(matrix.on_point(x.value, y.value))
    return np.array(points, dtype=float).reshape(-1, 2)


def _rank_values(
    positions: list[float], long: list[bool], majors: dict[float, int]
) -> list[float]:
    """Return the ticks' values: long, its label's 10^n; short, k x 10^n by rank.

    A short tick between the long ticks of 10^n and 10^(n+1) is the k-th above 10^n
    (k = 2 for the first); below the lowest long tick, counted down from 9; above the
    highest, counted up from 2.
    """
    values: list[float] = []
    exponents = sorted(majors.values())
    for index, position in enumerate(positions):
        if long[index]:
            values.append(10.0 ** majors[position])
            continue
        below = [p for p in majors if p < position]
        if below:
            base = max(below)
            rank = sum(1 for p in positions if base < p <= position)
            values.append((rank + 1) * 10.0 ** majors[base])
        else:
            base = min(majors)
            rank = sum(1 for p in positions if position <= p < base)
            values.append((10 - rank) * 10.0 ** (exponents[0] - 1))
    return values


def _fig6_ticks(paths: list[Path], panel: str, along: str) -> list[list[float]]:
    """Return a panel's ticks on its bottom (x) or left (y) frame: [position, length].

    A tick is a grey segment across the frame line, shorter than 2 pt: long (1.5 pt)
    or short (0.8 pt). Ticks drawn more than once (within ``MERGE`` pt) are one.
    """
    left, right, bottom, top = FIG6_FRAMES[panel]
    found: list[tuple[float, float]] = []
    for p in paths:
        if p.stroke != GREY or len(p.vertices) != 2:
            continue
        (x0, y0), (x1, y1) = p.vertices
        if along == "x":
            position, across, start, length = x0, x1, min(y0, y1), abs(y1 - y0)
            frame, low, high = bottom, left, right
        else:
            position, across, start, length = y0, y1, min(x0, x1), abs(x1 - x0)
            frame, low, high = left, bottom, top
        if (
            math.isclose(position, across, abs_tol=1e-3)
            and math.isclose(start, frame, abs_tol=0.01)
            and length < 2
            and low - MERGE <= position <= high + MERGE
        ):
            found.append((position, length))
    found.sort()
    ticks: list[list[float]] = []  # [position, longest length]
    for position, length in found:
        if ticks and position - ticks[-1][0] <= MERGE:
            ticks[-1][1] = max(ticks[-1][1], length)
        else:
            ticks.append([position, length])
    return ticks


def _fig6_axis(paths: list[Path], labels: list[Label], panel: str, along: str) -> Axis:
    """Return a panel's axis from its ticks, the long ones matched to their labels."""
    left, _, bottom, _ = FIG6_FRAMES[panel]
    ticks = _fig6_ticks(paths, panel, along)
    positions = [t[0] for t in ticks]
    long = [t[1] > LONG for t in ticks]
    majors: dict[float, int] = {}
    for position, is_long in zip(positions, long, strict=True):
        if not is_long:
            continue
        if along == "x":
            near = [
                lab for lab in labels if lab.y < bottom and abs(lab.x - position) < 3
            ]
        else:
            near = [
                lab
                for lab in labels
                if lab.x < left and left - lab.x < 15 and abs(lab.y - position) < 3
            ]
        if len(near) != 1:
            msg = f"Fig. 6 {panel} {along}: {len(near)} labels at {position}"
            raise ValueError(msg)
        majors[position] = near[0].exponent
    values = _rank_values(positions, long, majors)
    axis = Axis(tuple(positions), tuple(values), log=True, tolerance=0.05)
    axis.check()
    return axis


def _fig6() -> dict[str, Reading]:
    """Return Fig. 6's markers, panel by panel."""
    paths, labels = _fig6_objects()
    squares = []
    for p in paths:
        if p.stroke != BLUE or not p.stroked or len(p.vertices) != 5:
            continue
        corners = p.vertices[:4]
        side = float(np.ptp(corners[:, 0]))
        if SQUARE[0] <= side <= SQUARE[1]:
            squares.append((float(corners[:, 0].mean()), float(corners[:, 1].mean())))
    halves = [
        p.vertices
        for p in paths
        if p.stroke == BLUE
        and p.stroked
        and len(p.vertices) == 2
        and math.isclose(p.vertices[0, 0], p.vertices[1, 0], abs_tol=1e-3)
    ]
    out = {}
    for panel, (left, right, _, _) in FIG6_FRAMES.items():
        x_axis = _fig6_axis(paths, labels, panel, "x")
        y_axis = _fig6_axis(paths, labels, panel, "y")
        markers = []
        for x, y in sorted(s for s in squares if left < s[0] < right):
            ends = [
                v[1, 1]
                for v in halves
                if abs(v[0, 0] - x) < 0.02 and abs(v[0, 1] - y) < 0.02
            ]
            if len(ends) != 2 or not (min(ends) <= y <= max(ends)):
                msg = f"Fig. 6 {panel}: marker ({x:.2f}, {y:.2f}), bar ends {ends}"
                raise ValueError(msg)
            markers.append(Marker(x, y, max(ends), min(ends), "pt"))
        out[panel] = Reading(markers, x_axis, y_axis)
    return out


# --- the files ----------------------------------------------------------------------

FILES = {
    "fig3_He": "He/nie_2022__fig3_He.csv",
    **{
        f"fig4_{gas}_{pol}": f"{gas}/nie_2022__fig4_{gas}_{pol}.csv"
        for gas in FIG4_FRAMES
        for pol in WEIGHTS
    },
    **{
        f"fig6_H2_{panel}": f"H2/nie_2022__fig6_H2_{panel}.csv" for panel in FIG6_FRAMES
    },
}
COUNTS = {
    "fig3_He": 8,
    "fig4_Ar_LP": 10,
    "fig4_Ar_CP": 8,
    "fig4_Kr_LP": 10,
    "fig4_Kr_CP": 10,
    "fig4_Xe_LP": 9,
    "fig4_Xe_CP": 10,
    "fig6_H2_800nm": 11,
    "fig6_H2_400nm": 12,
    "fig6_H2_267nm": 13,
}


def readings() -> dict[str, Reading]:
    """Return every dataset's reading, by dataset id."""
    out = {"fig3_He": _fig3()}
    for (gas, polarisation), reading in _fig4().items():
        out[f"fig4_{gas}_{polarisation}"] = reading
    for panel, reading in _fig6().items():
        out[f"fig6_H2_{panel}"] = reading
    for dataset, count in COUNTS.items():
        if len(out[dataset].markers) != count:
            msg = f"{dataset}: {len(out[dataset].markers)} points, not {count}"
            raise ValueError(msg)
    return out


def extract() -> dict[str, list[Point]]:
    """Return the measured phase shifts of Figs. 3(a), 4 and 6."""
    check_parameters()
    files = {}
    for dataset, reading in readings().items():
        extras = reading.extras or [{} for _ in reading.markers]
        files[FILES[dataset]] = [
            _point(f"p{k:02d}", marker, reading.x_axis, reading.y_axis, extra)
            for k, (marker, extra) in enumerate(
                zip(reading.markers, extras, strict=True), start=1
            )
        ]
    return files


def bar_asymmetry() -> dict[str, tuple[float, float, int]]:
    """Return, per dataset, the range of upper over lower half-bar (rad) and the count.

    A bar symmetric in rad gives 1; on the logarithmic axis it is longer below the
    marker than above, in pixels.
    """
    out = {}
    for dataset, reading in readings().items():
        ratios = []
        for m in reading.markers:
            if m.top is None or m.bottom is None:
                continue
            phase = float(reading.y_axis(m.y))
            upper = abs(float(reading.y_axis(m.top)) - phase)
            lower = abs(phase - float(reading.y_axis(m.bottom)))
            ratios.append(upper / lower)
        if ratios:
            out[dataset] = (min(ratios), max(ratios), len(ratios))
    return out


def check_parameters() -> None:
    """Raise unless the record's beam parameters are those of Table 1 (p. 7)."""
    rows = table(DIRECTORY, TABLE_1)
    if rows[0] != [
        "Pump",
        "Pulse duration (fs)",
        "Spot size (µm)",
        "Probe",
        "Pulse duration (fs)",
        "Spot size (µm)",
    ]:
        msg = f"Table 1: unexpected header {rows[0]}"
        raise ValueError(msg)
    beams = {}
    for pump, pump_fs, pump_um, probe, probe_fs, probe_um in rows[1:]:
        beams[float(pump.removesuffix(" nm"))] = {
            "pulse_fs": float(pump_fs),
            "pump_waist_um": float(pump_um),
            "probe_wavelength_nm": float(probe.removesuffix(" nm")),
            "probe_pulse_fs": float(probe_fs),
            "probe_waist_um": float(probe_um),
        }
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        record = tomllib.load(handle)
    for dataset in record["datasets"]:
        parameters = dataset["parameters"]
        beam = beams[float(parameters["wavelength_nm"])]
        for name, value in beam.items():
            if not math.isclose(float(parameters[name]), value):
                msg = f"{dataset['id']}: {name} = {parameters[name]}, Table 1 {value}"
                raise ValueError(msg)


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
