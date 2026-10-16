# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract L'Huillier, Jönsson and Wendin (1986): two-photon ionisation of He, Fig. 11.

The PDF is a scan (one 5060 x 6898 bilevel image per page, about 612 dpi) with a text
layer; nothing is printed as numbers, so the curves are read from the page image at its
native pixels (data rule 8, TC-699).

* **Fig. 11** (p. 3948, ``fig_6afcab30a614``): the two-photon ionisation rate of He
  divided by I^2, p / I^2 in 10^-16 s^-1 W^-2 cm^4, against the photon energy in Ry:
  the random-phase approximation (solid line) and the one-electron approximation
  (dashed), local-density (LDA) basis. Each curve is drawn in two pieces, 0.62 to 1.02
  Ry and 1.20 Ry onwards (the gap around the LDA threshold, 1.13 Ry). A thin
  long-dashed interpolation joins the solid pieces (not read: the authors' sketch).
* **Fig. 13** (p. 3949, ``fig_4310d9a04e31``): the same quantity; Ritchie's TDHF curve
  (dashed) is read where it stands clear of the authors' solid line and Victor's
  circles, 0.91 to 1.33 Ry, only to check ritchie_1977's Table I against it
  (:func:`ritchie_check`); not data of this record.

How the figures are read:

* the axes: the inward tick marks of the frame, written in :data:`F11` and :data:`F13`
  (the bottom edge's x ticks, both edges' y ticks), each found in the image within
  :data:`FOUND` pixels. x: every 0.1 Ry from 0.7 to 2.4,
  except the ticks under the dotted threshold lines (0.9, 1.8) and under the curves
  near the axis (2.0, 2.1 in Fig. 11; 2.1 in Fig. 13). The frame's left edge, labelled
  '0.57' (Fig. 11) and '0.6' (Fig. 13), lies at 0.62 Ry on this scale: the curves are
  read on the tick scale, and the label is not used. y: every 0.5 from 0 (the frame's
  bottom line) to 6.5. The page is slightly rotated and the drawing's spacing differs
  a little between the two edges (up to 7 pixels): a pixel's value is interpolated
  between the left and the right edge's linear scales by its column.
* the curves: per sampled column (every :data:`EVERY` th), the dark runs (grey below
  :data:`DARK`) inside a row window, :data:`STROKE` pixels long. Each curve piece is
  tracked from a seed: the run whose centre lies nearest the row predicted from the
  last two samples, within :data:`JUMP` pixels; a column where none does (a gap of a
  dashed line) is skipped. The pieces' column ranges (:data:`F11_PIECES`) leave out
  where two lines touch: 1.28-1.55 Ry (the solid line meets the thin interpolation and
  then the dashed line) and beyond the 1.81 dotted line (they run into the axis).
  Columns within :data:`AWAY` pixels of a dotted threshold line or a bottom tick are
  skipped. Two curves never share a run (checked). The value is the run's centre;
  ``printed`` gives the pixel (column, row) and the run's length.
* which line is which beyond 1.55 Ry, where the two lie 12-15 pixels apart: the upper
  one is continuous and the lower one dashed (checked: every column with one run has it
  on the upper line's track), so the upper one is the solid line (RPA).

**Checks** (:func:`checks`; the script fails if one is missed): the text's 'at 1 Ry,
the two-photon ionization rate is lowered by a factor of 1.4' (p. 3949): the dashed
over the solid line at 1.0 Ry; 'This effect is reversed around 1.5 Ry' (p. 3949) and
'increased at higher energies above 1.45 Ry' (p. 3952): the dashed line above the
solid one at 1.20-1.24 Ry, below it beyond 1.54 Ry. Fig. 13's dashed
line against ritchie_1977's Table I, column E, converted to W / I^2 (that record's
definition, sigma_2 = 2 v a_0^2 E_p E_h / I_0): within :data:`RITCHIE` at Table I's
energies inside the range read (three).

**Reported** (:func:`comparisons`): Fig. 11's RPA line as sigma_2 against the curated He
calculations at 0.92-1.00 Ry.
"""

import functools
from dataclasses import dataclass

import numpy as np
import ritchie_1977 as ritchie
from _figure import Axis, image
from _library import Point, write
from scipy import ndimage

DIRECTORY = "LHuillier_1986_ManyElectronEffects"
SOURCE = "lhuillier_1986"
FILES = {
    "rpa": "He/lhuillier_1986__fig11_rpa_He.csv",
    "one_electron": "He/lhuillier_1986__fig11_one_electron_He.csv",
}
DARK = 128
FOUND = 1.5
"""Each written tick is found in the image within this, pixels."""
X_TOLERANCE = 10.0
"""The x ticks lie on the linear scale within this, pixels (the drawing's spacing)."""
Y_TOLERANCE = 6.0
EVERY = 4
STROKE = (4, 24)
"""A curve's run in a column is this many pixels long (longer on the steep left)."""
JUMP = 12.0
"""A run is taken for a curve within this of the predicted row, pixels."""
AWAY = 6
"""Columns within this of a dotted line or a bottom tick are skipped, pixels."""
TICK_HEIGHT = 50
"""The bottom ticks reach this far above the frame's bottom line, pixels (at most)."""


@dataclass(frozen=True)
class Frame:
    """One figure's page, frame and ticks (native pixels: column, row)."""

    page: int
    x_ticks: tuple[tuple[float, float], ...]
    """The bottom edge's ticks: (column, Ry)."""
    left: tuple[tuple[float, float], ...]
    """The left edge's ticks: (row, value); the first is the frame's bottom line."""
    right: tuple[tuple[float, float], ...]
    edges: tuple[float, float]
    """The columns of the frame's left and right edges."""
    bands: tuple[tuple[int, int], tuple[int, int], tuple[int, int]]
    """The rows of the bottom ticks' band, the columns of the left and right bands."""
    lines: tuple[float, ...]
    """The columns of the dotted threshold lines."""


F11 = Frame(
    page=11,
    x_ticks=(
        (730.0, 0.7),
        (813.0, 0.8),
        (990.0, 1.0),
        (1080.0, 1.1),
        (1169.0, 1.2),
        (1250.0, 1.3),
        (1340.0, 1.4),
        (1430.0, 1.5),
        (1513.5, 1.6),
        (1600.0, 1.7),
        (1775.0, 1.9),
        (2036.0, 2.2),
        (2127.0, 2.3),
        (2210.0, 2.4),
    ),
    left=(
        (5977.0, 0.0),
        (5904.5, 0.5),
        (5828.5, 1.0),
        (5752.0, 1.5),
        (5680.5, 2.0),
        (5604.0, 2.5),
        (5524.0, 3.0),
        (5449.5, 3.5),
        (5372.0, 4.0),
        (5295.0, 4.5),
        (5215.0, 5.0),
        (5141.5, 5.5),
        (5064.0, 6.0),
        (4987.5, 6.5),
    ),
    right=(
        (5971.0, 0.0),
        (5897.5, 0.5),
        (5822.0, 1.0),
        (5748.0, 1.5),
        (5677.0, 2.0),
        (5602.5, 2.5),
        (5522.0, 3.0),
        (5449.5, 3.5),
        (5372.5, 4.0),
        (5296.0, 4.5),
        (5216.0, 5.0),
        (5144.0, 5.5),
        (5068.0, 6.0),
        (4991.5, 6.5),
    ),
    edges=(658.0, 2302.0),
    bands=((5962, 5970), (662, 672), (2288, 2298)),
    lines=(904.0, 1095.0, 1687.0),
)
"""Fig. 11, PDF page 11. The bottom ticks at 0.9 and 1.8 lie under the dotted lines
of the thresholds ('0.90' and '1.81'), at 2.0 and 2.1 under the curves."""
F13 = Frame(
    page=12,
    x_ticks=(
        (739.5, 0.7),
        (822.5, 0.8),
        (999.0, 1.0),
        (1089.0, 1.1),
        (1177.5, 1.2),
        (1259.0, 1.3),
        (1347.5, 1.4),
        (1437.5, 1.5),
        (1521.0, 1.6),
        (1609.0, 1.7),
        (1783.0, 1.9),
        (1871.0, 2.0),
        (2043.0, 2.2),
        (2135.5, 2.3),
        (2220.0, 2.4),
    ),
    left=(
        (6067.0, 0.0),
        (5994.5, 0.5),
        (5919.0, 1.0),
        (5843.5, 1.5),
        (5769.0, 2.0),
        (5695.0, 2.5),
        (5616.0, 3.0),
        (5543.5, 3.5),
        (5466.0, 4.0),
        (5391.0, 4.5),
        (5310.5, 5.0),
        (5236.0, 5.5),
        (5158.5, 6.0),
        (5082.0, 6.5),
    ),
    right=(
        (6066.0, 0.0),
        (5992.5, 0.5),
        (5918.5, 1.0),
        (5843.0, 1.5),
        (5770.5, 2.0),
        (5697.0, 2.5),
        (5619.0, 3.0),
        (5545.5, 3.5),
        (5469.5, 4.0),
        (5395.0, 4.5),
        (5316.0, 5.0),
        (5242.0, 5.5),
        (5166.0, 6.0),
        (5089.0, 6.5),
    ),
    edges=(668.0, 2308.0),
    bands=((6052, 6063), (671, 681), (2296, 2306)),
    lines=(),
)
"""Fig. 13, PDF page 12. The bottom ticks at 0.9 and 1.8 lie under the dotted lines,
2.1 under the curve."""


@functools.cache
def _dark(page: int) -> np.ndarray:
    return image(DIRECTORY, page) < DARK


def _marks(band: np.ndarray, axis: int, start: int) -> list[float]:
    """Return the centres of the narrow dark features crossing a band."""
    counts = band.sum(axis=axis)
    labels, count = ndimage.label(counts >= 0.9 * band.shape[axis])
    out = []
    for k in range(1, count + 1):
        index = np.flatnonzero(labels == k)
        out.append(float(np.average(index, weights=counts[index])) + start)
    return out


def _check_ticks(frame: Frame) -> None:
    """Raise unless every written tick (but the frame's bottom line) is in the image."""
    dark = _dark(frame.page)
    (r0, r1), (l0, l1), (q0, q1) = frame.bands
    left_edge, right_edge = (round(e) for e in frame.edges)
    found = {
        "bottom": _marks(dark[r0:r1, left_edge + 4 : right_edge - 4], 0, left_edge + 4),
        "left": _marks(dark[4800:6100, l0:l1], 1, 4800),
        "right": _marks(dark[4800:6100, q0:q1], 1, 4800),
    }
    written = {
        "bottom": [c for c, _ in frame.x_ticks],
        "left": [r for r, _ in frame.left[1:]],
        "right": [r for r, _ in frame.right[1:]],
    }
    for name, positions in written.items():
        for position in positions:
            if not any(abs(m - position) <= FOUND for m in found[name]):
                msg = f"page {frame.page}: no {name} tick at {position}"
                raise ValueError(msg)


@functools.cache
def axes(frame: Frame) -> tuple[Axis, Axis, Axis]:
    """Return a figure's x axis (column, Ry) and its left and right y axes, checked."""
    _check_ticks(frame)
    x_axis = Axis(
        tuple(c for c, _ in frame.x_ticks),
        tuple(v for _, v in frame.x_ticks),
        tolerance=X_TOLERANCE,
    )
    left = Axis(
        tuple(r for r, _ in frame.left),
        tuple(v for _, v in frame.left),
        tolerance=Y_TOLERANCE,
    )
    right = Axis(
        tuple(r for r, _ in frame.right),
        tuple(v for _, v in frame.right),
        tolerance=Y_TOLERANCE,
    )
    for axis in (x_axis, left, right):
        axis.check()
    return x_axis, left, right


def value(frame: Frame, column: float, row: float) -> tuple[float, float]:
    """Return a pixel's photon energy (Ry) and value, the y scale interpolated."""
    x_axis, left, right = axes(frame)
    share = (column - frame.edges[0]) / (frame.edges[1] - frame.edges[0])
    y = (1 - share) * float(left(row)) + share * float(right(row))
    return float(x_axis(column)), y


def residuals(frame: Frame) -> tuple[float, float, float]:
    """Return the worst tick residual of the x, left and right axes, pixels."""
    out = []
    for axis in axes(frame):
        slope, intercept = np.polyfit(axis.positions, axis.values, 1)
        predicted = (np.asarray(axis.values) - intercept) / slope
        out.append(float(np.max(np.abs(predicted - np.asarray(axis.positions)))))
    return out[0], out[1], out[2]


def runs(page: int, column: int, window: tuple[int, int]) -> list[tuple[float, int]]:
    """Return a column's dark runs inside a row window: (centre row, length)."""
    cells = _dark(page)[window[0] : window[1], column].astype(int)
    edges = np.diff(np.concatenate([[0], cells, [0]]))
    starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
    return [
        ((s + e - 1) / 2 + window[0], int(e - s))
        for s, e in zip(starts, ends, strict=True)
        if STROKE[0] <= e - s <= STROKE[1]
    ]


@dataclass(frozen=True)
class Piece:
    """One drawn piece of a curve: its columns, row window and seed row."""

    curve: str
    columns: tuple[int, int]
    window: tuple[int, int]
    seed: float
    pair: str = ""
    """Where two lines run close ('upper' or 'lower'): with two runs in the window
    take that one; with one run, the upper line only (it is the continuous one)."""


F11_PIECES = (
    Piece("one_electron", (702, 1013), (5000, 5560), 5092.0),
    Piece("rpa", (702, 1013), (5150, 5610), 5218.0),
    Piece("one_electron", (1162, 1334), (5600, 5880), 5619.0),
    Piece("rpa", (1162, 1238), (5795, 5830), 5809.5),
    Piece("one_electron", (1470, 1676), (5880, 5960), 5911.0, "lower"),
    Piece("rpa", (1470, 1676), (5880, 5960), 5898.0, "upper"),
)
"""Fig. 11's pieces read. The first two: 0.67-1.02 Ry (the left edge's ticks, which
reach column 698, hide 0.62-0.67); then 1.20 Ry onwards, where the solid line runs
clear of the thin interpolation to 1.28 Ry and the dashed line clear of both to 1.39
Ry; and 1.55-1.78 Ry, the two lines 12-15 pixels apart, up to the 1.81 dotted line."""
F13_PIECE = Piece("ritchie", (924, 1290), (5670, 5945), 5696.0)
"""Fig. 13's dashed line (Ritchie's), 0.91-1.33 Ry, the lowest line there."""


def _skip(frame: Frame, piece: Piece, column: int) -> bool:
    """Return whether a column is near a dotted line, or a tick the piece reaches."""
    near = list(frame.lines)
    if piece.window[1] > frame.left[0][0] - TICK_HEIGHT:
        near += [c for c, _ in frame.x_ticks]
    return any(abs(column - c) <= AWAY for c in near)


def track(frame: Frame, piece: Piece) -> list[tuple[int, float, int]]:
    """Return a piece's samples: (column, run centre, run length), tracked."""
    samples: list[tuple[int, float, int]] = []
    for column in range(piece.columns[0], piece.columns[1] + 1, EVERY):
        if _skip(frame, piece, column):
            continue
        if len(samples) >= 2:
            (c1, r1, _), (c2, r2, _) = samples[-2], samples[-1]
            predicted = r2 + (r2 - r1) / (c2 - c1) * (column - c2)
        elif samples:
            predicted = samples[-1][1]
        else:
            predicted = piece.seed
        found = runs(frame.page, column, piece.window)
        if piece.pair:
            if len(found) == 2:
                found = [found[0] if piece.pair == "upper" else found[1]]
            elif len(found) != 1 or piece.pair != "upper":
                continue
        if not found:
            continue
        centre, length = min(found, key=lambda run: abs(run[0] - predicted))
        if abs(centre - predicted) <= JUMP:
            samples.append((column, centre, length))
    if len(samples) < 5:
        msg = f"{piece.curve} {piece.columns}: {len(samples)} samples"
        raise ValueError(msg)
    return samples


@functools.cache
def fig11() -> dict[str, tuple[tuple[int, float, int], ...]]:
    """Return Fig. 11's samples per curve, checked: no run taken for both curves."""
    out: dict[str, list[tuple[int, float, int]]] = {"rpa": [], "one_electron": []}
    for piece in F11_PIECES:
        out[piece.curve].extend(track(F11, piece))
    shared = {(c, r) for c, r, _ in out["rpa"]} & {
        (c, r) for c, r, _ in out["one_electron"]
    }
    if shared:
        msg = f"Fig. 11: runs taken for both curves: {sorted(shared)}"
        raise ValueError(msg)
    return {name: tuple(sorted(samples)) for name, samples in out.items()}


def _upper_is_continuous() -> int:
    """Return how many one-run columns beyond 1.55 Ry lie on the upper (RPA) track.

    Raises if one lies nearer the lower (dashed) line.
    """
    piece = F11_PIECES[-1]
    rpa = [(c, r) for c, r, _ in fig11()["rpa"] if c >= piece.columns[0]]
    dashed = [(c, r) for c, r, _ in fig11()["one_electron"] if c >= piece.columns[0]]
    upper_c, upper_r = zip(*rpa, strict=True)
    lower_c, lower_r = zip(*dashed, strict=True)
    single = 0
    for column in range(piece.columns[0], piece.columns[1] + 1, EVERY):
        if _skip(F11, piece, column):
            continue
        found = runs(F11.page, column, piece.window)
        if len(found) != 1:
            continue
        row = found[0][0]
        upper = abs(row - float(np.interp(column, upper_c, upper_r)))
        lower = abs(row - float(np.interp(column, lower_c, lower_r)))
        if lower < upper:
            msg = f"Fig. 11: a lone run at column {column} is on the lower line"
            raise ValueError(msg)
        single += 1
    return single


def points(curve: str) -> list[Point]:
    """Return one of Fig. 11's curves as data points, one per sampled column."""
    out = []
    for column, row, length in fig11()[curve]:
        x, y = value(F11, column, row)
        out.append(
            Point(
                f"c{column}",
                f"{x:.4f}",
                f"{y:.3f}",
                "",
                f"pixel ({column}, {row:.1f}), run {length}",
            )
        )
    return out


def _curve(curve: str) -> tuple[np.ndarray, np.ndarray]:
    values = [value(F11, c, r) for c, r, _ in fig11()[curve]]
    x, y = zip(*values, strict=True)
    return np.asarray(x), np.asarray(y)


def at(curve: str, energy: float, half_width: float = 0.03) -> float:
    """Return a Fig. 11 curve at a photon energy (Ry), a local quadratic fit."""
    x, y = _curve(curve)
    near = np.abs(x - energy) <= half_width
    if near.sum() < 4 or x[near].min() > energy or x[near].max() < energy:
        msg = f"{curve}: too few samples around {energy} Ry"
        raise ValueError(msg)
    return float(np.polyval(np.polyfit(x[near], y[near], 2), energy))


FACTOR = 1.4
"""p. 3949: 'at 1 Ry, the two-photon ionization rate is lowered by a factor of 1.4'."""
CROSSING = (1.24, 1.54)
"""p. 3949: 'This effect is reversed around 1.5 Ry'; p. 3952: 'increased at higher
energies above 1.45 Ry'. The dashed line lies above the solid one below the first
energy and below it beyond the second (the pieces read: the lines touch between)."""
RITCHIE = 0.10
"""Fig. 13's dashed line against ritchie_1977's Table I (as W / I^2) at its energies:
the reading's uncertainty there (the x scale's residual, 0.01 Ry, on a slope of about
5 per Ry: 2-4 %; the y scales', 0.05; the dash's centre, 0.01) and the drawing's."""
RYDBERG_PER_HARTREE = 2.0


@functools.cache
def fig13() -> tuple[tuple[int, float, int], ...]:
    """Return Fig. 13's dashed (Ritchie's) line, tracked."""
    return tuple(track(F13, F13_PIECE))


def ritchie_check() -> dict[str, tuple[float, float, float]]:
    """Return Fig. 13's dashed line and ritchie_1977's Table I E at its energies.

    Per Table I energy inside the range read (au): (Ry, drawn W / I^2, Table I's value
    as W / I^2), in 10^-16 s^-1 W^-2 cm^4.
    """
    values = [value(F13, c, r) for c, r, _ in fig13()]
    x, y = (np.asarray(v) for v in zip(*values, strict=True))
    out = {}
    for row in ritchie.rows("tbl1"):
        energy = float(row[0]) * RYDBERG_PER_HARTREE
        near = np.abs(x - energy) <= 0.04
        if near.sum() < 4 or x[near].min() > energy or x[near].max() < energy:
            continue
        drawn = float(np.polyval(np.polyfit(x[near], y[near], 2), energy))
        printed = ritchie.w_over_i2(float(row[4]), float(row[0])) / 1e-16
        out[row[0]] = (energy, drawn, printed)
    return out


def checks() -> dict[str, object]:
    """Return the checks of Figs. 11 and 13; raise if one is missed."""
    one, rpa = at("one_electron", 1.0), at("rpa", 1.0)
    if round(one / rpa, 1) != FACTOR:
        msg = f"Fig. 11 at 1 Ry: dashed / solid = {one / rpa:.3f}, not 1.4"
        raise ValueError(msg)
    x_one, y_one = _curve("one_electron")
    x_rpa, y_rpa = _curve("rpa")
    below = [
        float(np.interp(e, x_one, y_one) - np.interp(e, x_rpa, y_rpa))
        for e in x_rpa
        if 1.15 < e < CROSSING[0]
    ]
    beyond = [
        float(np.interp(e, x_one, y_one) - np.interp(e, x_rpa, y_rpa))
        for e in x_rpa
        if e > CROSSING[1]
    ]
    if min(below) <= 0 or max(beyond) >= 0:
        msg = "Fig. 11: the dashed line does not cross the solid one at 1.24-1.54 Ry"
        raise ValueError(msg)
    single = _upper_is_continuous()
    ritchie_values = ritchie_check()
    if len(ritchie_values) < 3:
        msg = f"Fig. 13: {len(ritchie_values)} of Table I's energies in range"
        raise ValueError(msg)
    for energy, (_, drawn, printed) in ritchie_values.items():
        if abs(drawn / printed - 1) > RITCHIE:
            msg = f"Fig. 13 at {energy} au: {drawn:.3f}, ritchie_1977 {printed:.3f}"
            raise ValueError(msg)
    return {
        "Fig. 11 at 1 Ry, dashed, solid, ratio": (one, rpa, one / rpa),
        "Fig. 11, dashed minus solid, 1.15-1.24 Ry (min) and beyond 1.54 (max)": (
            min(below),
            max(beyond),
        ),
        "Fig. 11, lone runs beyond 1.55 Ry on the solid line's track": single,
        "Fig. 13 against ritchie_1977 (Ry, drawn, Table I E as W/I^2, ratio)": {
            energy: (ry, drawn, printed, drawn / printed)
            for energy, (ry, drawn, printed) in ritchie_values.items()
        },
    }


COMPARE_RY = (0.92, 0.95, 1.0)
"""Fig. 11's RPA line compared where the real two-photon ionisation lies below the
LDA threshold (0.90-1.13 Ry, the authors' region (ii)), in Ry."""


def comparisons() -> dict[float, dict[str, float | None]]:
    """Return Fig. 11's RPA line as sigma_2 over the curated He calculations.

    sigma_2 = (p / I^2) (hbar omega)^2 (W = sigma_2 Phi^2), at :data:`COMPARE_RY`;
    the others as in ritchie_1977's comparison, and OEP 6's He model.
    """
    out: dict[float, dict[str, float | None]] = {}
    for energy in COMPARE_RY:
        photon = energy / RYDBERG_PER_HARTREE * ritchie.AU_ENERGY
        ours = at("rpa", energy) * 1e-16 * photon**2
        ev = energy / RYDBERG_PER_HARTREE * ritchie.HARTREE_EV
        ratios: dict[str, float | None] = {"sigma_2 (cm4 s)": ours, "eV": ev}
        for label, keys in ritchie.OTHERS.items():
            parts = [ritchie._at(key, ev) for key in keys]
            known = [p for p in parts if p is not None]
            ratios[label] = ours / sum(known) if len(known) == len(parts) else None
        ratios["OEP 6 He model"] = ours / ritchie.model_cm4_s(
            energy / RYDBERG_PER_HARTREE
        )
        out[energy] = ratios
    return out


def extract() -> dict[str, list[Point]]:
    """Return Fig. 11's two curves, after the checks."""
    checks()
    return {relative: points(curve) for curve, relative in FILES.items()}


def main() -> None:
    """Write the data files, and print the checks."""
    write(extract())
    for name, frame in (("Fig. 11", F11), ("Fig. 13", F13)):
        x_res, left_res, right_res = residuals(frame)
        print(
            f"{name} ticks within {x_res:.1f} (x), {left_res:.1f} (left) and "
            f"{right_res:.1f} (right) pixels"
        )
    for curve in FILES:
        x, _ = _curve(curve)
        print(f"Fig. 11 {curve}: {len(x)} samples, {x.min():.3f}-{x.max():.3f} Ry")
    for label, result in checks().items():
        print(f"check {label}: {result}")
    for energy, ratios in comparisons().items():
        print(f"Fig. 11 RPA at {energy} Ry: {ratios}")


if __name__ == "__main__":
    main()
