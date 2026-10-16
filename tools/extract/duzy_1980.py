# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Duzy and Hyman (1980): photoionisation of excited Ne, Ar, Kr and Xe.

The library copy is a scan of Phys. Rev. A 22, 1878 (1980) (PDF page = journal page -
1877), each page one embedded image (4807 x 6331 pixels) with a text layer.

* **Table II** (p. 1882, ``tbl_398a22d38504``, JSON only: merged header and atom
  cells): the cross sections, in 10^-18 cm^2, of the first excited p and d states of
  Ar and Kr at 248 nm (KrF) and of Ne and Xe at 308 nm (XeCl), one value each. Each is
  one file of one point; the printed row is kept in ``printed``. Checked against the
  table image.
* **Fig. 3** (p. 1880, ``fig_1e192e00a6b7``): the cross sections (cm^2, log axis,
  1e-20 to 1e-16) of the first excited p states against the photon wavelength (nm,
  linear, 580 at the left frame to 80 at the right), four curves told apart by their
  dashes: Ne 3p short dashes, Ar 4p long dashes, Kr 5p solid, Xe 6p dash-dot. A raster
  scan: read from the page image at native resolution. The Xe 6p curve is written; the
  others are read only at Table II's wavelengths, and the Kr 5p curve at the
  wavelengths of dakka_2018 and payne_1981 (:func:`checks`).

The figure is read as follows:

* the axes: the frame is hand-drafted and the scan is turned by about 0.25 degree, so
  each side has its own calibration from its tick marks (written below: the centre of
  the dark run of each tick ``SCAN`` pixels inside the frame), and a point is
  interpolated between the two sides by its position across the frame. The sigma
  ticks (1, 2, 4, 6, 8 x 10^k) are not evenly spaced decade to decade (559 to 545
  pixels per decade), so each decade has its own logarithmic :class:`Axis` from its
  six ticks (every tick within ``TOLERANCE_Y`` of its fit); the wavelength ticks
  (580 to 80 nm, every 100) give one linear :class:`Axis` per side (within
  ``TOLERANCE_X``). Every tick is checked to be a dark run within ``FOUND`` pixels.
* the runs: in a column, the dark runs between ``MARGIN`` pixels inside the top and
  bottom frame lines, runs closer than ``JOIN`` pixels joined and specks shorter than
  ``SPECK`` left out. A line is ``ONE`` pixels wide or less (8-12 measured).
* the Xe 6p curve, by three rules, each where it identifies the curve without doubt
  (the ranges were found by inspection of the image and are written in ``RULES``):

  - ``isolated`` (514-386 nm): the dash-dot curve is the top one, alone. Traced from
    the column of ``SEED_NM`` both ways, the top run of a column is taken if it is a
    single line and lies within ``GUARD`` pixels of the row extrapolated from the two
    last points taken (a column in a gap of the dashes, where the top run is the Kr
    curve, is skipped).
  - ``coincident`` (335-245 nm): the dash-dot and solid curves are drawn on each other,
    and the top run of a column is both; it is taken where it is at most ``TWO``
    pixels wide (one line, or the two not quite superposed). Where it is wider, the
    curves separate or the Ar curve joins them (386-335 nm and 245-222 nm): not read.
  - ``between`` (222-104 nm): the Kr, Xe and Ar curves are three separate lines, Xe the
    middle one; a column is taken where exactly three runs lie within ``SPAN`` pixels
    of the top one (the Ne curve lies farther below), each at most ``STEEP`` pixels
    tall (a steep line's vertical chord). Where the Xe and Ar curves merge (200-146
    nm) or a dash gap leaves two runs, nothing is taken.

  A value is the centre of the run; one point is written every ``EVERY`` columns of
  those taken; ``printed`` holds the pixel and the rule.

One pixel is 0.35 nm in wavelength and 0.0018 decade (0.41 %) in sigma; the ticks lie
within 2 pixels of their fits. The reading's uncertainty, 8 %, is stated in the record.
"""

import math
import re
from dataclasses import dataclass
from functools import cache

import numpy as np
from _figure import Axis, image
from _json_table import json_table
from _library import Point, write

DIRECTORY = "Duzy_1980_PhotoionizationExcitedRare"
TABLE_II = "tbl_398a22d38504"
FIGURE_3 = "fig_1e192e00a6b7"
PDF_PAGE = 3
"""Fig. 3 is on journal p. 1880, PDF p. 3, the page's one image."""

GASES = {"Ar*": "Ar", "Kr*": "Kr", "Ne*": "Ne", "Xe*": "Xe"}
"""Table II's atom cells, in order."""
STATES = {
    "Ar": ("4p", "3d"),
    "Kr": ("5p", "4d"),
    "Ne": ("3p", "3d"),
    "Xe": ("6p", "5d"),
}
"""Each gas's first excited p and d states, Table II's rows."""
FULL = re.compile(r"^(\d)p\^\{5\}(\d)p$")
"""Table II's p-state cells, ``3p^{5}4p``."""
WAVELENGTHS = (r"\lambda = 248 \text{ nm (KrF)}", r"\lambda = 308 \text{ nm (XeCl)}")
UNIT = r"Cross section (10^{-18} \text{ cm}^{2})"

DARK = 128
SCAN = 15
"""How far inside a frame line, pixels, the ticks are read."""
FOUND = 3.0
"""The farthest, pixels, a written tick may be from the centre of a dark run."""
TOLERANCE_Y = 2.5
"""The largest residual, pixels, of a sigma tick on its decade's fit."""
TOLERANCE_X = 3.5
"""The largest residual, pixels, of a wavelength tick on its side's fit."""
MARGIN = 50
"""Rows this close to the top or bottom frame line are not searched."""
EDGE = 80
"""Columns this close to the right frame line are not searched (its long ticks)."""
JOIN = 3
"""Runs of a column closer than this, pixels, are one line (the scan breaks lines)."""
SPECK = 4
"""Runs shorter than this, pixels, are specks of the scan, not lines."""
ONE = 13
"""The widest run, pixels, that is one line (the lines are 8-12 pixels)."""
TWO = 17
"""The widest top run, pixels, read as the Kr and Xe curves drawn on each other."""
STEEP = 30
"""The tallest run, pixels, of one line where the curves are steep (below 222 nm)."""
SPAN = 300
"""Kr, Xe and Ar lie within this many pixels of the top run below 222 nm; Ne farther."""
GUARD = 6.0
"""The farthest, pixels, an isolated Xe run may be from its extrapolation."""
BASE = 10
"""The extrapolation's slope is taken over at least this many columns."""
GAP = 6
"""The fewest white pixels between the Xe run and the Kr and Ar runs below 222 nm."""
SEED_NM = 500.0
EVERY = 10
"""A point is written every this many columns of those taken."""

COLUMNS = (2875.5, 4298.5)
"""The left and right frame lines' columns (at mid-height)."""
ROWS = (3101.5, 5322.0)
"""The top and bottom frame lines' rows (at mid-width)."""
SIGMA_TICKS = (
    1e-16,
    *(m * 10.0**k for k in range(-17, -21, -1) for m in (8, 6, 4, 2, 1)),
)
"""The sigma ticks from the top: 1e-16, 8e-17, 6e-17, 4e-17, 2e-17, ..., 1e-20."""
LEFT = (
    *(3098.5, 3156.0, 3225.0, 3322.5, 3491.5, 3659.5, 3713.0, 3782.5, 3881.5, 4050.0),
    *(4218.5, 4273.0, 4344.5, 4440.5, 4607.5, 4772.5, 4826.5, 4894.5, 4989.5, 5153.5),
    5318.5,
)
"""Rows of the sigma ticks on the left side (the first and last are the frame)."""
RIGHT = (
    *(3105.0, 3158.5, 3226.5, 3325.0, 3492.5, 3661.5, 3714.5, 3784.5, 3883.5, 4051.5),
    *(4221.0, 4274.5, 4345.5, 4442.5, 4607.5, 4774.5, 4826.5, 4895.5, 4990.5, 5155.5),
    5325.0,
)
WAVELENGTH_TICKS = (580.0, 480.0, 380.0, 280.0, 180.0, 80.0)
BOTTOM = (2874.5, 3156.0, 3442.5, 3728.5, 4013.5, 4297.5)
"""Columns of the bottom side's wavelength ticks (the first and last: the frame)."""
TOP = (2875.5, 3158.0, 3447.0, 3734.0, 4016.0, 4299.0)

RULES = (
    ("isolated", 514.0, 386.0),
    ("coincident", 335.0, 245.0),
    ("between", 222.0, 200.0),
    ("between", 146.0, 0.0),
)
"""Each rule and the wavelengths (nm, long to short) where it is tried."""


def _spans(line: np.ndarray, join: int) -> list[tuple[int, int]]:
    """Return the first and last index of each run of True, joining near runs."""
    index = np.flatnonzero(line)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > join) + 1)
    return [(int(group[0]), int(group[-1])) for group in groups]


@dataclass(frozen=True)
class Calibrated:
    """Fig. 3's side axes (the sigma sides by decade), for reading a pixel."""

    left: tuple[Axis, ...]
    right: tuple[Axis, ...]
    top: Axis
    bottom: Axis

    @staticmethod
    def _sigma(decades: tuple[Axis, ...], row: float) -> float:
        """Return log10 sigma at a row from the decade whose ticks bracket it."""
        for axis in decades:
            if row <= axis.positions[-1]:
                return math.log10(float(axis(row)))
        return math.log10(float(decades[-1](row)))

    def value(self, column: float, row: float) -> tuple[float, float]:
        """Return (wavelength in nm, sigma in cm^2) at a pixel."""
        across = (column - COLUMNS[0]) / (COLUMNS[1] - COLUMNS[0])
        down = (row - ROWS[0]) / (ROWS[1] - ROWS[0])
        wavelength = (1 - down) * float(self.top(column)) + down * float(
            self.bottom(column)
        )
        log_s = (1 - across) * self._sigma(self.left, row) + across * self._sigma(
            self.right, row
        )
        return wavelength, 10**log_s

    def column(self, wavelength: float) -> int:
        """Return the column of a wavelength at mid-height."""
        slope = (WAVELENGTH_TICKS[-1] - WAVELENGTH_TICKS[0]) / (BOTTOM[-1] - BOTTOM[0])
        guess = BOTTOM[0] + (wavelength - WAVELENGTH_TICKS[0]) / slope
        for _ in range(3):  # Newton on the mid-height wavelength
            here = self.value(guess, sum(ROWS) / 2)[0]
            guess += (wavelength - here) / slope
        return round(guess)


def _decades(rows: tuple[float, ...]) -> tuple[Axis, ...]:
    """Return a side's four decade axes, each from its six ticks."""
    return tuple(
        Axis(
            rows[5 * k : 5 * k + 6],
            SIGMA_TICKS[5 * k : 5 * k + 6],
            log=True,
            tolerance=TOLERANCE_Y,
        )
        for k in range(4)
    )


def _check_ticks(dark: np.ndarray) -> None:
    """Raise unless every written tick is a dark run where it is written."""
    left, right = (round(c) for c in COLUMNS)
    top, bottom = (round(r) for r in ROWS)
    sides = {
        "left": (LEFT, dark[:, left + SCAN]),
        "right": (RIGHT, dark[:, right - SCAN]),
        "top": (TOP, dark[top + SCAN, :]),
        "bottom": (BOTTOM, dark[bottom - SCAN, :]),
    }
    for side, (ticks, line) in sides.items():
        centres = np.array([(a + b) / 2 for a, b in _spans(line, 1)])
        for position in ticks:
            if np.min(np.abs(centres - position)) > FOUND:
                msg = f"Fig. 3: no {side} tick at {position}"
                raise ValueError(msg)


@cache
def calibrate() -> tuple[Calibrated, np.ndarray]:
    """Return Fig. 3's calibration (its ticks checked) and its page's dark pixels."""
    dark = image(DIRECTORY, PDF_PAGE) < DARK
    _check_ticks(dark)
    calibrated = Calibrated(
        _decades(LEFT),
        _decades(RIGHT),
        Axis(TOP, WAVELENGTH_TICKS, tolerance=TOLERANCE_X),
        Axis(BOTTOM, WAVELENGTH_TICKS, tolerance=TOLERANCE_X),
    )
    for axis in (
        *calibrated.left,
        *calibrated.right,
        calibrated.top,
        calibrated.bottom,
    ):
        axis.check()
    return calibrated, dark


def _runs(column: int) -> list[tuple[int, int]]:
    """Return the dark runs (first and last row) of a column inside the margins."""
    _, dark = calibrate()
    first = round(ROWS[0]) + MARGIN
    last = round(ROWS[1]) - MARGIN
    return [
        (first + a, first + b)
        for a, b in _spans(dark[first:last, column], JOIN)
        if b - a + 1 >= SPECK
    ]


def _centre(run: tuple[int, int]) -> float:
    """Return the centre row of a run."""
    return (run[0] + run[1]) / 2


def _width(run: tuple[int, int]) -> int:
    """Return the height of a run in pixels."""
    return run[1] - run[0] + 1


def _three(column: int) -> list[tuple[int, int]] | None:
    """Return the Kr, Xe and Ar runs of a column below 222 nm, or None if unresolved."""
    runs = _runs(column)
    if not runs:
        return None
    near = [r for r in runs if _centre(r) - _centre(runs[0]) <= SPAN]
    if len(near) != 3 or any(_width(r) > STEEP for r in near):
        return None
    if near[1][0] - near[0][1] <= GAP or near[2][0] - near[1][1] <= GAP:
        return None
    return near


def _isolated(start: int, stop: int, seed: int) -> dict[int, tuple[int, int, str]]:
    """Return the isolated Xe curve's runs by column, traced from the seed."""
    found: dict[int, tuple[int, int, str]] = {}
    for direction in (1, -1):
        taken: list[tuple[int, float]] = []
        column = seed
        while start <= column <= stop:
            runs = _runs(column)
            if runs and _width(runs[0]) <= ONE:
                row = _centre(runs[0])
                base = [t for t in taken if abs(taken[-1][0] - t[0]) >= BASE]
                if not base:
                    near = not taken or abs(row - taken[-1][1]) <= GUARD
                else:
                    (c1, r1), (c2, r2) = base[-1], taken[-1]
                    predicted = r2 + (r2 - r1) / (c2 - c1) * (column - c2)
                    near = abs(row - predicted) <= GUARD
                if near:
                    found[column] = (*runs[0], "isolated")
                    taken.append((column, row))
            column += direction
    return found


@cache
def xe_6p() -> list[tuple[int, float, float, float, str]]:
    """Return the Xe 6p curve: (column, row, nm, sigma cm^2, rule) per column taken."""
    calibrated, _ = calibrate()
    found: dict[int, tuple[int, int, str]] = {}
    for rule, long, short in RULES:
        start = calibrated.column(long)
        stop = min(calibrated.column(short), round(COLUMNS[1]) - EDGE)
        if rule == "isolated":
            found |= _isolated(start, stop, calibrated.column(SEED_NM))
            continue
        for column in range(start, stop + 1):
            if rule == "coincident":
                runs = _runs(column)
                if runs and _width(runs[0]) <= TWO:
                    found[column] = (*runs[0], rule)
            else:
                three = _three(column)
                if three is not None:
                    found[column] = (*three[1], rule)
    out = []
    for column in sorted(found):
        a, b, rule = found[column]
        wavelength, sigma = calibrated.value(column, (a + b) / 2)
        out.append((column, (a + b) / 2, wavelength, sigma, rule))
    return out


def table_ii() -> dict[tuple[str, str], tuple[str, str, str]]:
    """Return Table II by (gas, state): wavelength (nm), value as printed, row."""
    rows = json_table(DIRECTORY, TABLE_II)
    if rows[0][:3] != ["Atom", "State", UNIT] or list(rows[1][2:]) != list(WAVELENGTHS):
        msg = f"Table II: headers moved {rows[:2]}"
        raise ValueError(msg)
    out = {}
    atom = ""
    for cells in rows[2:]:
        atom = cells[0] or atom
        gas = GASES[atom]
        match = FULL.match(cells[1])
        state = f"{match[2]}p" if match else cells[1]
        if state not in STATES[gas] or (match and state != STATES[gas][0]):
            msg = f"Table II: an unexpected row {cells}"
            raise ValueError(msg)
        filled = [i for i in (2, 3) if cells[i]]
        if len(filled) != 1:
            msg = f"Table II: {gas} {state} has no single wavelength"
            raise ValueError(msg)
        wavelength = "248" if filled[0] == 2 else "308"
        out[(gas, state)] = (
            wavelength,
            cells[filled[0]],
            " | ".join([atom, *cells[1:]]),
        )
    if len(out) != 8 or list(dict.fromkeys(g for g, _ in out)) != list(STATES):
        msg = "Table II: rows moved"
        raise ValueError(msg)
    return out


def _read(wavelength: float, rank: int, count: int) -> tuple[float, float]:
    """Return (nm, sigma cm^2) of the rank-th run where a column has ``count`` runs.

    The column nearest the wavelength (searching outward up to 15 columns) whose runs
    number ``count`` and are each one line (``ONE``), except the top run, which may be
    the Kr and Xe curves drawn on each other (``TWO``).
    """
    calibrated, _ = calibrate()
    centre = calibrated.column(wavelength)
    for offset in (0, *(s * d for d in range(1, 16) for s in (1, -1))):
        runs = _runs(centre + offset)
        if (
            len(runs) == count
            and _width(runs[0]) <= TWO
            and all(_width(r) <= ONE for r in runs[1:])
        ):
            return calibrated.value(centre + offset, _centre(runs[rank]))
    msg = f"no column near {wavelength} nm has {count} single runs"
    raise ValueError(msg)


FIGURE_OVER_TABLE = {
    ("Ne", "3p"): 0.98,
    ("Ar", "4p"): 0.84,
    ("Kr", "5p"): 1.07,
    ("Xe", "6p"): 1.12,
}
"""Fig. 3 read at Table II's wavelength over Table II, as found (the register)."""
READ_TOLERANCE = 0.04
"""How far a ratio may be from the one found: the reading (about 2 % each way)."""
KR_NM = 215.0
"""dakka_2018's wavelength for Kr 5p[3/2]2, where the top of three lines is Kr 5p."""
HYMAN_1977 = {"Ar": 4.3, "Kr": 4.5}
"""hyman_1977's Table I at 2486 A, 1e-18 cm^2: 'the earlier calculation' (p. 1880)."""


def checks() -> dict[str, object]:
    """Return Fig. 3 read at Table II's wavelengths, and the text's comparisons.

    At 248 nm the column holds three runs: the Kr and Xe curves on each other, Ar, Ne;
    at 308 nm the same. Raises unless each ratio of the figure to Table II is the one
    found (``FIGURE_OVER_TABLE``) within ``READ_TOLERANCE``.
    """
    out: dict[str, object] = {}
    printed = table_ii()
    rank = {"Kr": 0, "Xe": 0, "Ar": 1, "Ne": 2}
    for (gas, state), expected in FIGURE_OVER_TABLE.items():
        wavelength, value, _ = printed[(gas, state)]
        nm, sigma = _read(float(wavelength), rank[gas], 3)
        ratio = sigma / (float(value) * 1e-18)
        out[
            f"{gas} {state} at {nm:.1f} nm: Fig. 3 {sigma:.3g} cm2 / Table II {value}"
        ] = round(ratio, 3)
        if abs(ratio - expected) > READ_TOLERANCE:
            msg = f"{gas} {state}: the figure is {ratio:.3f} x Table II"
            raise ValueError(msg)
    for gas, earlier in HYMAN_1977.items():
        value = printed[(gas, STATES[gas][0])][1]
        out[f"{gas} Table II / hyman_1977 {earlier}"] = round(float(value) / earlier, 3)
    for wavelength, quoted in ((249.6, "payne_1981, Xe"), (326.2, "payne_1981, Kr")):
        nm, sigma = _read(wavelength, 0, 3)
        out[f"Kr 5p and Xe 6p (one line) at {nm:.1f} nm ({quoted}), cm2"] = (
            f"{sigma:.3g}"
        )
    calibrated, _ = calibrate()
    column = calibrated.column(KR_NM)
    three = _three(column)
    if three is None:
        msg = f"Kr 5p: the curves are not resolved at {KR_NM} nm"
        raise ValueError(msg)
    nm, sigma = calibrated.value(column, _centre(three[0]))
    out[f"Kr 5p at {nm:.1f} nm (dakka_2018), cm2"] = f"{sigma:.3g}"
    trace = xe_6p()
    out["Xe 6p: columns taken by rule"] = {
        rule: sum(1 for p in trace if p[4] == rule) for rule, _, _ in RULES
    }
    return out


def extract() -> dict[str, list[Point]]:
    """Return Table II's values and the Xe 6p curve of Fig. 3."""
    checks()
    files: dict[str, list[Point]] = {}
    for (gas, state), (wavelength, value, row) in table_ii().items():
        files[f"{gas}/duzy_1980__table2_{state}.csv"] = [
            Point(wavelength, wavelength, value, "", row)
        ]
    points = []
    last = -EVERY
    for column, row, wavelength, sigma, rule in xe_6p():
        if column - last < EVERY:
            continue
        last = column
        points.append(
            Point(
                f"{column}",
                f"{wavelength:.1f}",
                f"{sigma:.3e}",
                "",
                f"pixel ({column}, {row:.1f}); {rule}",
            )
        )
    files["Xe/duzy_1980__fig3_6p.csv"] = points
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
