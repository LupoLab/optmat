# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hummer and Storey (1998): photoionisation of He's singlet S and D levels.

The library copy is the publisher's PDF of Mon. Not. R. Astron. Soc. 297, 1073 (1998)
(PDF page = journal page - 1072); the text and tables are text, the two figures are
embedded images.

* **Table 1** (p. 1074, ``tbl_af91b1b8d77f``, JSON only): threshold photoionisation
  cross sections (Mb) of the n = 1-4 1S, 1P, 1D, 3S, 3P and 3D levels, in four
  columns: 'Extrapolated' (from Drake's bound-bound oscillator strengths, with the
  fitting error), 'Coulomb', 'OP' and 'This Work'. Taken: the singlet S and D rows
  (the levels two photons reach from the ground state), Extrapolated and This Work,
  one file per level and column, one point at the level's threshold. Checked against
  the table image.
* **Table 2** (p. 1075, ``tbl_e11b5c5787d2``, JSON only): the quantum defects of the
  1S and 1P levels; This Work's give the calculation's thresholds, of the 1S levels
  only.
* **Fig. 1** (p. 1076, the first image of the page, 3383 x 2200 pixels): log10 of
  the cross sections (Mb) of the n = 1, 2, 5, 10, 15, 20, 25 1S levels against the
  electron energy (Ryd, 0 to 2.5). Taken: n = 1 and 2, read from the image.

The threshold of a Table 1 point is the calculation's: for the n 1S rows, the binding
energy 1 / (n - delta)^2 Ry with This Work's quantum defect (Table 2: 0.2549, 0.1478,
0.1419, 0.1401 for n = 1-4), in R_infinity hc ('rydberg'), for both columns (TC-934's
correction); the paper prints no D quantum defects and no level energies, so the 1D
rows' threshold is ``not stated``. The NIST ionisation energies of the levels
(nist_asd_2026, ``data/external/nist_asd/He_I_levels.tsv``: the He II limit less the
level) are used only by :func:`checks`, which compares the calculated thresholds with
them.

The photon energy of a Fig. 1 point is its electron energy plus the calculation's
binding energy of the level, 1 / (n - delta)^2 Ry with This Work's quantum defect
(Table 2: n = 1, 0.2549; n = 2, 0.1478), in R_infinity hc ('rydberg'); the quantum
defects are those of the calculated energies.

Fig. 1 is read as follows (:func:`curve`):

* the axes: each side of the frame has its own linear :class:`Axis`, calibrated on
  its labelled ticks (written below, the pixel of each, read from the image; the
  right side's ticks are unlabelled and mirror the left's), and a point is
  interpolated between the two sides by its position across the frame. Every
  labelled tick is checked to be a dark run of its side's scan line (``SCAN`` pixels
  inside the frame) within ``FOUND`` pixels, every fit within ``TOLERANCE``, and every
  minor tick (each 0.5 of log10 sigma, each 0.1 Ryd) predicted by the fit is checked
  to fall on a dark run within ``MINOR`` pixels (on the left side, the n = 5-25
  curves rise along the frame and cover the 1.5 tick).
* the curves: traced from ``SEED`` Ryd, where the first and second dark runs below
  the top margin of the column are the n = 1 and n = 2 curves, column by column both
  ways, each time to the run nearest the previous row within ``STEP`` pixels; the
  trace stops where the run is wider than ``WIDEST`` pixels (the curves meet near
  threshold, or n = 2 turns down into the 2s2p 1P resonance near 2.5 Ryd) and at
  ``LAST`` Ryd. A value is the centre of the run.
* the points written: one every ``EVERY`` columns (0.05 Ryd) from ``FIRST`` Ryd
  (nearer threshold the two curves lie within two line widths of each other; they
  cross at 0.05 Ryd); ``printed`` holds the pixel and the electron energy read.

One pixel is 0.00088 Ryd and 0.0048 in log10 sigma (1.1 %); the lines are 3-6 pixels
wide.
"""

import re
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np
from _figure import Axis, image
from _json_table import json_table
from _library import DATA, Point, write

from optmat.datasets import NOT_STATED

DIRECTORY = "Hummer_1998_RecombinationHeliumLike"
TABLE_1 = "tbl_af91b1b8d77f"
TABLE_2 = "tbl_e11b5c5787d2"
LEVELS = {
    "1^{1}S": ("1s2", "1S", "0", "1s2_1S"),
    "2^{1}S": ("1s.2s", "1S", "0", "1s2s_1S"),
    "3^{1}S": ("1s.3s", "1S", "0", "1s3s_1S"),
    "4^{1}S": ("1s.4s", "1S", "0", "1s4s_1S"),
    "3^{1}D": ("1s.3d", "1D", "2", "1s3d_1D"),
    "4^{1}D": ("1s.4d", "1D", "2", "1s4d_1D"),
}
"""Table 1's rows taken: NIST's configuration, term and J, and the file label."""
ROWS = (
    "1^{1}S",
    "2^{1}S",
    "3^{1}S",
    "4^{1}S",
    "2^{1}P^o",
    "3^{1}P^o",
    "4^{1}P^o",
    "3^{1}D",
    "4^{1}D",
    "2^{3}S",
    "3^{3}S",
    "4^{3}S",
    "2^{3}P^o",
    "3^{3}P^o",
    "4^{3}P^o",
    "3^{3}D",
    "4^{3}D",
)
"""Table 1's rows, in order."""
NIST_LEVELS = DATA / "external" / "nist_asd" / "He_I_levels.tsv"
LIMIT = "He II (1s 2S<1/2>)"
"""NIST's row of the ionisation limit (the configuration column)."""
CELL = re.compile(r"^\d+\.\d+$")
ERROR = re.compile(r"^±(\d\.\d{3})$")

FIGURE_PAGE = 4
FIGURE_INDEX = 0
DARK = 128
SCAN = 8
FOUND = 2.0
"""The farthest, pixels, a written tick may be from the centre of a dark run."""
TOLERANCE = 2.0
"""The largest residual, pixels, of a written tick on its side's fit."""
MINOR = 3.0
"""The farthest, pixels, a minor tick predicted by the fit may be from a dark run."""
MARGIN = 20
"""Rows this close to the top or bottom frame line are not searched."""
EDGE = 20
"""Columns this close to the left or right frame line are not searched."""
SEED = 1.0
"""The electron energy (Ryd) the traces start at."""
LAST = 2.45
"""The electron energy (Ryd) the traces stop at on the right."""
STEP = 12
"""The farthest, pixels, a curve's row may move from one column to the next."""
WIDEST = 14
"""The widest run, pixels, a value is read from (a wider one: the curves meet)."""
EVERY = 57
"""A point is written every this many columns (0.05 Ryd)."""
DEFECTS = {1: "0.2549", 2: "0.1478", 3: "0.1419", 4: "0.1401"}
"""Table 2, This Work, the 1S quantum defects of n = 1-4 (checked by the script)."""
S_ROWS = {"1^{1}S": 1, "2^{1}S": 2, "3^{1}S": 3, "4^{1}S": 4}
"""Table 1's 1S rows and their n, placed at the calculation's thresholds."""


@dataclass(frozen=True)
class Frame:
    """Fig. 1's frame lines and labelled ticks, (pixel, value), per side."""

    columns: tuple[int, int]
    rows: tuple[int, int]
    left: tuple[tuple[float, float], ...]
    right: tuple[tuple[float, float], ...]
    top: tuple[tuple[float, float], ...]
    bottom: tuple[tuple[float, float], ...]


Y_LABELS = (2.0, 0.0, -2.0, -4.0, -6.0)
X_LABELS = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5)
# Labelled ticks read from the image (centres of the dark runs SCAN pixels inside each
# side; the bottom frame line is the -6 tick, the left and right ones 0 and 2.5 Ryd).
FRAME = Frame(
    columns=(437, 3267),
    rows=(40, 1924),
    left=tuple(zip((255.0, 673.5, 1089.5, 1504.0, 1923.0), Y_LABELS, strict=True)),
    right=tuple(zip((248.0, 667.5, 1085.0, 1501.0, 1920.5), Y_LABELS, strict=True)),
    top=tuple(
        zip((434.0, 999.5, 1567.5, 2134.5, 2700.0, 3267.5), X_LABELS, strict=True)
    ),
    bottom=tuple(
        zip((438.0, 1004.0, 1571.5, 2139.0, 2703.0, 3270.0), X_LABELS, strict=True)
    ),
)


def _spans(line: np.ndarray) -> list[tuple[int, int]]:
    """Return the first and last index of each run of True."""
    index = np.flatnonzero(line)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [(int(group[0]), int(group[-1])) for group in groups]


def _axis(ticks: tuple[tuple[float, float], ...]) -> Axis:
    """Return a linear axis calibrated on a side's labelled ticks."""
    return Axis(
        tuple(p for p, _ in ticks), tuple(v for _, v in ticks), False, TOLERANCE
    )


def _check_ticks(dark: np.ndarray) -> None:
    """Raise unless the labelled and minor ticks are dark runs where they should be."""
    left, right = FRAME.columns
    top, bottom = FRAME.rows
    sides = {
        "left": (FRAME.left, dark[:, left + SCAN], 0.5),
        "right": (FRAME.right, dark[:, right - SCAN], 0.5),
        "top": (FRAME.top, dark[top + SCAN, :], 0.1),
        "bottom": (FRAME.bottom, dark[bottom - SCAN, :], 0.1),
    }
    for side, (ticks, line, spacing) in sides.items():
        spans = _spans(line)
        centres = np.array([(a + b) / 2 for a, b in spans])
        _axis(ticks).check()
        for position, value in ticks:
            if np.min(np.abs(centres - position)) > FOUND:
                msg = f"Fig. 1: no {side} tick at {position} ({value})"
                raise ValueError(msg)
        positions = np.array([p for p, _ in ticks])
        values = np.array([v for _, v in ticks])
        slope, intercept = np.polyfit(values, positions, 1)
        low, high = sorted((values[0], values[-1]))
        for k in range(round((high - low) / spacing) + 1):
            guess = slope * (low + k * spacing) + intercept
            if not any(a - MINOR <= guess <= b + MINOR for a, b in spans):
                msg = f"Fig. 1: no {side} minor tick near {guess:.0f}"
                raise ValueError(msg)


@dataclass(frozen=True)
class Calibrated:
    """Fig. 1's four side axes, for reading a pixel."""

    left: Axis
    right: Axis
    top: Axis
    bottom: Axis

    def value(self, column: float, row: float) -> tuple[float, float]:
        """Return (electron energy in Ryd, log10 sigma in Mb) at a pixel."""
        left, right = FRAME.columns
        top, bottom = FRAME.rows
        across = (column - left) / (right - left)
        down = (row - top) / (bottom - top)
        energy = (1 - down) * float(self.top(column)) + down * float(
            self.bottom(column)
        )
        log_sigma = (1 - across) * float(self.left(row)) + across * float(
            self.right(row)
        )
        return energy, log_sigma

    def column(self, energy: float) -> int:
        """Return the column of an electron energy at mid-height."""
        slope = (
            float(self.top(1.0))
            - float(self.top(0.0))
            + float(self.bottom(1.0))
            - float(self.bottom(0.0))
        ) / 2
        intercept = (float(self.top(0.0)) + float(self.bottom(0.0))) / 2
        return round((energy - intercept) / slope)


@cache
def calibrate() -> tuple[Calibrated, np.ndarray]:
    """Return Fig. 1's calibration (its ticks checked) and its dark pixels."""
    dark = image(DIRECTORY, FIGURE_PAGE, FIGURE_INDEX) < DARK
    _check_ticks(dark)
    calibrated = Calibrated(
        _axis(FRAME.left), _axis(FRAME.right), _axis(FRAME.top), _axis(FRAME.bottom)
    )
    return calibrated, dark


@cache
def curve(n: int) -> list[tuple[int, float, float, float]]:
    """Return the n = 1 or 2 curve: (column, row, Ryd, log10 Mb) per column."""
    calibrated, dark = calibrate()
    left, right = FRAME.columns
    top, bottom = FRAME.rows
    first = top + MARGIN
    seed = calibrated.column(SEED)
    runs = _spans(dark[first : bottom - MARGIN, seed])
    a, b = runs[n - 1]
    start = first + (a + b) / 2
    stop = min(calibrated.column(LAST), right - EDGE)
    found: dict[int, tuple[int, int]] = {}
    for direction in (1, -1):
        current, column = start, seed
        while left + EDGE <= column <= stop:
            near = [
                (first + a, first + b)
                for a, b in _spans(dark[first : bottom - MARGIN, column])
                if abs(first + (a + b) / 2 - current) <= STEP
            ]
            if not near:
                break
            a, b = min(near, key=lambda span: abs((span[0] + span[1]) / 2 - current))
            if b - a + 1 > WIDEST:
                break
            found[column] = (a, b)
            current = (a + b) / 2
            column += direction
    out = []
    for column in sorted(found):
        a, b = found[column]
        energy, log_sigma = calibrated.value(column, (a + b) / 2)
        out.append((column, (a + b) / 2, energy, log_sigma))
    return out


def table_1() -> dict[str, list[str]]:
    """Return Table 1's rows by state: Extrapolated, error, Coulomb, OP, This Work."""
    rows = json_table(DIRECTORY, TABLE_1)
    if rows[0] != ["State", "Extrapolated*", "", "Coulomb", "OP", "This Work"]:
        msg = f"Table 1: headers moved {rows[0]}"
        raise ValueError(msg)
    out = {}
    for cells in rows[1:]:
        if (
            not all(CELL.match(cells[i]) for i in (1, 3, 4, 5))
            or ERROR.match(cells[2]) is None
        ):
            msg = f"Table 1: an unexpected row {cells}"
            raise ValueError(msg)
        out[cells[0]] = cells[1:]
    if tuple(out) != ROWS:
        msg = "Table 1: rows moved"
        raise ValueError(msg)
    return out


def table_2() -> dict[int, list[str]]:
    """Return Table 2's rows by n: the 1S and 1P quantum defects, as printed."""
    rows = json_table(DIRECTORY, TABLE_2)
    if rows[1] != ["n", *(["Experiment^†", "OP*", "This Work"] * 2)]:
        msg = f"Table 2: headers moved {rows[1]}"
        raise ValueError(msg)
    return {int(cells[0]): cells[1:] for cells in rows[2:]}


def thresholds() -> dict[str, str]:
    """Return the NIST ionisation energy (cm^-1) of each level taken, as text.

    Used only to compare the calculated thresholds with NIST's (:func:`checks`).
    """
    lines = [
        [cell.strip('"') for cell in line.split("\t")]
        for line in Path(NIST_LEVELS).read_text(encoding="utf-8").splitlines()[1:]
    ]
    limit = [line[4] for line in lines if line[0] == LIMIT and line[1] == "Limit"]
    if len(limit) != 1:
        msg = "NIST: no single ionisation limit"
        raise ValueError(msg)
    out = {}
    for state, (configuration, term, j, _) in LEVELS.items():
        found = [line[4] for line in lines if line[:3] == [configuration, term, j]]
        if len(found) != 1:
            msg = f"NIST: {len(found)} levels {configuration} {term} {j}"
            raise ValueError(msg)
        decimals = len(limit[0].split(".")[1])
        out[state] = f"{float(limit[0]) - float(found[0]):.{decimals}f}"
    return out


def binding(n: int) -> float:
    """Return the calculation's binding energy (Ryd) of n 1S: 1 / (n - delta)^2."""
    delta = table_2()[n][2]
    if delta != DEFECTS[n]:
        msg = f"Table 2: n = {n} This Work is {delta}"
        raise ValueError(msg)
    return 1 / (n - float(delta)) ** 2


KAPRALOVA = 1.727
"""kapralova_zdanska_2025's implied 1s2s cross section, Mb, at 10.31 eV (record)."""
KUMAR = 1.691
"""kumar_2026's implied 1s2s cross section, Mb, at 10.3088 eV (record)."""
KUMAR_EV = 1239.842 / 120.27
"""kumar_2026's photon energy (120.27 nm), eV."""
KAPRALOVA_EV = 1239.842 / 120.3
"""kapralova_zdanska_2025's photon energy (120.3 nm), eV."""
KUMAR_XUV_EV = 1239.842 / 65.05
KUMAR_XUV = 0.459
"""kumar_2026's implied 1s2s cross section, Mb, at 19.06 eV (65.05 nm) (record)."""
RYDBERG_EV = 13.605693122994
"""R_infinity hc in eV (CODATA 2018)."""


def _at(trace: list[tuple[int, float, float, float]], energy: float) -> float:
    """Return a curve's sigma (Mb) at an electron energy (Ryd), interpolating log10."""
    energies = np.array([p[2] for p in trace])
    if not energies.min() <= energy <= energies.max():
        msg = f"{energy:.3g} Ryd is outside the curve read"
        raise ValueError(msg)
    return float(10 ** np.interp(energy, energies, [p[3] for p in trace]))


def checks() -> dict[str, object]:
    """Return the figure against Table 1, and the 1s2s comparisons.

    Raises unless the n = 1 curve, carried to threshold by a straight line fitted to
    log10 sigma from ``FIRST`` to ``FIT_TO`` Ryd,
    is within 3 % of Table 1's This Work (7.394 Mb), unless both curves are read to
    within ``NEAR`` Ryd of threshold (where they cross), or unless the calculated
    thresholds of the n = 1-4 1S levels are within 0.4 % of NIST's.
    """
    out: dict[str, object] = {}
    rows = table_1()
    nist = thresholds()
    for n in (1, 2):
        trace = curve(n)
        if trace[0][2] > NEAR:
            msg = f"n = {n}: the curve read starts at {trace[0][2]:.3f} Ryd"
            raise ValueError(msg)
        out[f"n = {n}: read from {trace[0][2]:.4f} to {trace[-1][2]:.3f} Ryd"] = len(
            trace
        )
    for state, n in S_ROWS.items():
        ev = float(nist[state]) * 1.239841984e-4
        ratio = binding(n) * RYDBERG_EV / ev
        out[f"n = {n}: calculated / NIST threshold"] = round(ratio, 5)
        if abs(ratio - 1) > 0.004:
            msg = f"n = {n}: the calculated threshold is {ratio:.4f} x NIST's"
            raise ValueError(msg)
    head = [p for p in curve(1) if FIRST <= p[2] <= FIT_TO]
    _, intercept = np.polyfit([p[2] for p in head], [p[3] for p in head], 1)
    printed = float(rows["1^{1}S"][4])
    out["n = 1: curve carried to threshold / Table 1 This Work"] = round(
        10**intercept / printed, 4
    )
    if abs(10**intercept / printed - 1) > 0.03:
        msg = f"n = 1: the curve gives {10**intercept:.2f} Mb at threshold"
        raise ValueError(msg)
    trace = curve(2)
    for label, photon, quoted in (
        ("kumar_2026 at 120.27 nm", KUMAR_EV, KUMAR),
        ("kapralova_zdanska_2025 at 120.3 nm", KAPRALOVA_EV, KAPRALOVA),
        ("kumar_2026 at 65.05 nm", KUMAR_XUV_EV, KUMAR_XUV),
    ):
        electron = photon / RYDBERG_EV - binding(2)
        sigma = _at(trace, electron)
        out[f"1s2s 1S at {photon:.3f} eV: Fig. 1 (Mb)"] = round(sigma, 3)
        out[f"1s2s 1S: {label} {quoted} Mb / Fig. 1"] = round(quoted / sigma, 3)
    return out


NEAR = 0.06
"""The curves must be read to within this electron energy (Ryd) of threshold."""
FIRST = 0.08
"""The electron energy (Ryd) of the first point written: below it the n = 1 and 2
curves lie within two line widths of each other (they cross at 0.05 Ryd)."""
FIT_TO = 0.3
"""The n = 1 curve from ``FIRST`` to this (Ryd) is fitted to carry it to threshold."""


def extract() -> dict[str, list[Point]]:
    """Return Table 1's singlet S and D values and Fig. 1's n = 1 and 2 curves."""
    checks()
    files: dict[str, list[Point]] = {}
    for state, cells in table_1().items():
        if state not in LEVELS:
            continue
        label = LEVELS[state][3]
        row = " | ".join([state, *cells])
        if state in S_ROWS:
            n = S_ROWS[state]
            threshold = f"{binding(n):.6f}"
            row += f" || Table 2 n = {n} This Work delta {DEFECTS[n]}"
        else:
            threshold = NOT_STATED
        extrapolated = Point(state, threshold, cells[0], cells[1][1:], row)
        this_work = Point(state, threshold, cells[4], "", row)
        files[f"He/hummer_1998__table1_extrapolated_{label}.csv"] = [extrapolated]
        files[f"He/hummer_1998__table1_{label}.csv"] = [this_work]
    for n, label in ((1, "1s2_1S"), (2, "1s2s_1S")):
        threshold = binding(n)
        points = []
        for column, row, energy, log_sigma in [p for p in curve(n) if p[2] >= FIRST][
            ::EVERY
        ]:
            points.append(
                Point(
                    f"{column}",
                    f"{energy + threshold:.5f}",
                    f"{10**log_sigma:.4g}",
                    "",
                    f"pixel ({column}, {row:.1f}); electron {energy:.4f} Ryd",
                )
            )
        files[f"He/hummer_1998__fig1_{label}.csv"] = points
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
