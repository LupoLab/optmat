# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Petrov et al. (2011): photoionisation of excited np levels of Ar, Kr, Xe.

Eur. Phys. J. D 62, 347-359 (2011); journal p. 347 is PDF page 1, so the journal page
is the PDF page plus 346. CIPFCP calculations (length gauge) of the one-photon
photoionisation cross sections of the ten levels mp5 (m+1)p, J = 0-3 (Paschen 2p1-2p10)
of Ar, Kr and Xe above the mp5 2P1/2 threshold. Taken:

* **Table 2** (PDF p. 6, ``tbl_63001e98ed74``; read from the library's JSON cells, every
  cell checked against the table image, 2026-10-12): for each gas and level, the total
  cross section at the 2P1/2 threshold (sigma_tot^0, a), the partial into the 2P3/2
  ion at that threshold (sigma_3/2^0, b) and the partial into the 2P1/2 ion at the
  relative photon energy omega_1/2 = Delta_FS (sigma_1/2^Delta_FS, c). No unit is
  printed; it is the figures' Mb (the totals agree with Figs. 3-5, :func:`checks`).
  The rows are the Paschen labels of argon with the jK label beside them; footnote e:
  for Kr 2p2/2p3 and Xe 2p8/2p9 the gas's own Paschen labels are swapped, so a row is
  identified by its jK label (as the figures' panels confirm, :func:`checks`). The
  ratio column r (d) is derived and not taken (checked).
* **Figs. 3, 4, 5** (PDF p. 5, journal p. 351; ``fig_cdcdef4cd679``,
  ``fig_c5f6d8609180``, ``fig_f7565eddc4c9``): the total cross section sigma_tot
  (the solid black line) of each level against the relative photon energy
  omega_1/2, 0-6 eV, read from the vector drawing (data rule 8). The partial curves
  (the blue sigma_3/2 and the dashed red sigma_1/2, several drawn x10) are not read.

How the drawing is read: each figure is ten panels (five rows, two columns). The frames
and ticks are stroked black lines 0.61 pt wide, drawn as move-line pairs; a frame edge
is a segment longer than 40 pt, a tick a segment shorter than 4 pt with one end on an
edge: long ticks (3.25 pt) at the labelled values, short ones (1.65 pt) halfway. The y
labels are text objects beside the outer edge of each column (0, 4, 8, 12 [, 16, 20]
Mb); each long tick takes the label within :data:`LABEL_DY` of it, and every tick, long
and short, lies on the fitted scale (:meth:`_figure.Axis.check`). The x labels are one
text object '0246' under each column (the digits are not separate objects); the long
ticks are 0, 2, 4, 6 eV as printed on the rendered page, the short ones 1, 3, 5. The
total is the one black stroke 0.85 pt wide with more than two vertices whose middle
lies in the panel (a two-vertex one in the first panel is the legend); it is a single
polyline of line segments, drawn left to right and clipped by the frame: the vertices
outside the frame (beyond 6 eV) are dropped. Each panel's level is read from its own
labels ('2p' with the Paschen subscript, the jK label with its J subscript), checked
against :data:`PANELS`.

The abscissa written is the photon energy, converted from the printed omega_1/2 (a
decision recorded in the record's notes): E_photon = omega_1/2 + E(2P1/2 limit) -
E(level), with NIST's levels (``data/external/nist_asd/<X>_I_levels.tsv``: the level,
the 2P3/2 limit; ``<X>_II_levels.tsv``: the ion's 2P1/2 level). The paper's own
energies are experimental (Fig. 1's caption, ref. 55). sigma_1/2^Delta_FS is placed at
omega_1/2 = Delta_FS, NIST's ion fine-structure splitting (the paper prints 0.178,
0.666 and 1.306 eV, p. 349; 177.5 meV for Ar, p. 353). The ``printed`` column keeps the
row's cells, or the vertex in page points with its omega_1/2.

**Checks** (:func:`checks`; the script fails if one is missed): every axis's ticks; the
panel labels; each curve's value at omega_1/2 = 0 against Table 2's sigma_tot^0 of
the level with the same jK label (all 30, within :data:`AT_ZERO`); r = sigma_1/2^DFS /
sigma_3/2^0 against the printed ratio (all 30, to the rounding of the three numbers);
and the paper's three Delta_FS against NIST's.
"""

import csv
import functools
import re
from dataclasses import dataclass

import numpy as np
from _figure import Axis, painted, texts
from _json_table import json_table
from _library import DATA, Point, write

from optmat.constants import E_CHARGE, HBAR, C

DIRECTORY = "Petrov_2011_NearThresholdPhotoionization"
SOURCE = "petrov_2011"
TABLE_2 = "tbl_63001e98ed74"
PAGE = 5
GASES = ("Ar", "Kr", "Xe")
FIGURES = {
    "Ar": (60.0, 450.0, 290.0, 760.0),
    "Kr": (330.0, 460.0, 580.0, 760.0),
    "Xe": (330.0, 110.0, 580.0, 415.0),
}
"""Each figure's region on PDF page 5 (left, bottom, right, top), pt."""
FRAME_WIDTH = (0.60, 0.62)
CURVE_WIDTH = (0.84, 0.86)
EDGE = 0.02
LABEL_DY = 1.5
START = 1.0
"""A curve's first vertex lies within this of the frame's left edge (0.07 eV), pt."""
AT_ZERO = 0.1
"""Largest difference, Mb, of a curve at omega_1/2 = 0 from Table 2's sigma_tot^0."""
CM_PER_EV = E_CHARGE / (2 * np.pi * HBAR * C * 100)
"""Wavenumbers (cm^-1) per electronvolt."""
DELTA_FS = {"Ar": 0.178, "Kr": 0.666, "Xe": 1.306}
"""The paper's ion fine-structure splittings, eV (p. 349)."""

# The panels as printed (rows top to bottom, left then right): the gas's Paschen label
# and the jK label, the core written ' for 2P1/2.
_ORDER = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
PANELS = {
    "Ar": ("p'[1/2]0", "p'[1/2]1", "p'[3/2]2", "p'[3/2]1", "p[1/2]0", "p[3/2]2",
           "p[3/2]1", "p[5/2]2", "p[5/2]3", "p[1/2]1"),
    "Kr": ("p'[1/2]0", "p'[3/2]2", "p'[1/2]1", "p'[3/2]1", "p[1/2]0", "p[3/2]2",
           "p[3/2]1", "p[5/2]2", "p[5/2]3", "p[1/2]1"),
    "Xe": ("p'[1/2]0", "p'[1/2]1", "p'[3/2]2", "p'[3/2]1", "p[1/2]0", "p[3/2]2",
           "p[3/2]1", "p[5/2]3", "p[5/2]2", "p[1/2]1"),
}  # fmt: skip
SHELL = {"Ar": "4", "Kr": "5", "Xe": "6"}
_NIST = {"Ar": "3s2.3p5", "Kr": "4s2.4p5", "Xe": "5p5"}
_ION = {"Ar": "3s2.3p5", "Kr": "4s2.4p5", "Xe": "5s2.5p5"}
_JK = re.compile(r"^p('?)\[(\d)/2\](\d)$")


def tag(jk: str) -> str:
    """Return a dataset-id fragment of a jK label: p'[1/2]0 -> pp12_0."""
    match = _JK.match(jk)
    if match is None:
        msg = f"unexpected jK label {jk!r}"
        raise ValueError(msg)
    prime, k, j = match.groups()
    return f"{'pp' if prime else 'p'}{k}2_{j}"


def _table_jk(cell: str) -> str:
    """Return Table 2's jK cell, [1/2]_0', as the panel's p'[1/2]0."""
    match = re.match(r"^\[(\d)/2\]_(\d)('?)$", cell.strip())
    if match is None:
        msg = f"unexpected jK cell {cell!r}"
        raise ValueError(msg)
    k, j, prime = match.groups()
    return f"p{prime}[{k}/2]{j}"


@functools.cache
def _levels(gas: str) -> dict[str, float]:
    """Return NIST's level energies (cm^-1) of the gas's np levels by jK, and limits."""
    out: dict[str, float] = {}
    base = DATA / "external" / "nist_asd"
    with (base / f"{gas}_I_levels.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            conf = row["Configuration"].strip('"')
            term, j = row["Term"].strip('"'), row["J"].strip('"')
            energy = row["Level (cm-1)"].strip('"')
            if term == "Limit" and "2P*<3/2>" in conf:
                out["limit"] = float(energy)
            for core, prime in (("2P*<3/2>", ""), ("2P*<1/2>", "'")):
                if conf == f"{_NIST[gas]}.({core}).{SHELL[gas]}p":
                    out[f"p{prime}{term.removeprefix('2')}{j}"] = float(energy)
    with (base / f"{gas}_II_levels.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            conf = row["Configuration"].strip('"')
            term = (row["Term"].strip('"'), row["J"].strip('"'))
            if conf == _ION[gas] and term == ("2P*", "1/2"):
                out["core"] = float(row["Level (cm-1)"].strip('"'))
    if len(out) != 12:
        msg = f"{gas}: NIST levels {sorted(out)}"
        raise ValueError(msg)
    return out


def threshold(gas: str, jk: str) -> float:
    """Return the photon energy (eV) of the level's 2P1/2 threshold, omega_1/2 = 0."""
    levels = _levels(gas)
    return (levels["limit"] + levels["core"] - levels[jk]) / CM_PER_EV


def splitting(gas: str) -> float:
    """Return NIST's ion fine-structure splitting Delta_FS, eV."""
    return _levels(gas)["core"] / CM_PER_EV


@dataclass(frozen=True)
class Panel:
    """One panel: its frame, axes, labels and the total's vertices (page points)."""

    gas: str
    paschen: str
    jk: str
    x: Axis
    y: Axis
    curve: np.ndarray


def _segments(gas: str) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Return the figure's frame-and-tick segments and its total curves."""
    lines, curves = [], []
    for path in painted(DIRECTORY, PAGE, FIGURES[gas]):
        if path.stroke != (0, 0, 0):
            continue
        if FRAME_WIDTH[0] < path.width < FRAME_WIDTH[1]:
            lines.extend(path.vertices.reshape(-1, 2, 2))
        elif CURVE_WIDTH[0] < path.width < CURVE_WIDTH[1] and len(path.vertices) > 2:
            curves.append(path.vertices)
    return lines, curves


def _clusters(values: list[float], gap: float = 0.5) -> list[float]:
    """Return the distinct values of a list (within a gap), sorted."""
    out: list[float] = []
    for value in sorted(values):
        if not out or value - out[-1] > gap:
            out.append(value)
    return out


def _label(found: list[tuple[str, float, float]], key: str) -> str:
    """Return a panel's label text with its subscript (the next text down-right)."""
    heads = [(s, x, y) for s, x, y in found if s.startswith(key)]
    if len(heads) != 1:
        msg = f"{key!r}: {heads}"
        raise ValueError(msg)
    head, x, y = heads[0]
    subs = [
        s
        for s, sx, sy in found
        if s.isdigit() and x < sx < x + 25 and y - 5 < sy < y - 1.5
    ]
    if len(subs) != 1:
        msg = f"{head!r}: subscripts {subs}"
        raise ValueError(msg)
    return head + subs[0]


def _panels(gas: str) -> list[Panel]:
    """Return the ten panels of a gas's figure, rows top to bottom, left then right."""
    lines, curves = _segments(gas)
    long = [s for s in lines if np.hypot(*(s[1] - s[0])) > 40]
    xs = _clusters([p[0] for s in long if abs(s[0, 0] - s[1, 0]) < EDGE for p in s])
    ys = _clusters([p[1] for s in long if abs(s[0, 1] - s[1, 1]) < EDGE for p in s])
    if len(xs) != 4 or len(ys) != 6:
        msg = f"{gas}: frame edges x {xs}, y {ys}"
        raise ValueError(msg)
    words = texts(DIRECTORY, PAGE, FIGURES[gas])
    panels = []
    for row in range(5):
        bottom, top = ys[4 - row], ys[5 - row]
        for column in range(2):
            left, right = xs[2 * column], xs[2 * column + 1]
            index = 2 * row + column
            panels.append(
                _panel(gas, index, (left, bottom, right, top), lines, curves, words)
            )
    return panels


def _y_axis(
    name: str,
    index: int,
    frame: tuple[float, float, float, float],
    ticks: list[np.ndarray],
    words: list[tuple[str, float, float]],
) -> Axis:
    """Return a panel's y axis, calibrated from its outer edge's ticks and labels.

    The labels are the digits outside the frame; the ticks have one end on the outer
    edge (left column: left; right: right). The top edge carries no tick of its own
    (it is above the last label), only the panel above's zero.
    """
    left, bottom, right, top = frame
    edge = left if index % 2 == 0 else right
    on_y = [
        s
        for s in ticks
        if abs(s[0, 1] - s[1, 1]) < EDGE
        and min(abs(s[:, 0] - edge)) < EDGE
        and bottom - EDGE <= s[0, 1] < top - 0.5
    ]
    labels = [
        (float(w), wy)
        for w, wx, wy in words
        if w.isdigit()
        and 0 < (edge - wx if index % 2 == 0 else wx - edge) < 9
        and bottom - 2 < wy < top + 2
    ]
    # a tick drawn twice (the edge shared with the panel above or below) counts once
    lengths: dict[float, float] = {}
    for s in on_y:
        position = float(s[0, 1])
        key = next((k for k in lengths if abs(k - position) < 0.05), position)
        lengths[key] = max(lengths.get(key, 0.0), float(np.hypot(*(s[1] - s[0]))))
    positions, values = [], []
    for position, length in sorted(lengths.items()):
        if length > 2.5:
            near = [v for v, wy in labels if abs(wy - position) < LABEL_DY]
            if len(near) != 1:
                msg = f"{name} panel {index}: y tick {position:.2f} labels {near}"
                raise ValueError(msg)
            positions.append(position)
            values.append(near[0])
    if len(positions) < 4:
        msg = f"{name} panel {index}: {len(positions)} labelled y ticks"
        raise ValueError(msg)
    step = (values[1] - values[0]) / 2
    major = Axis(tuple(positions), tuple(values))
    for position, length in sorted(lengths.items()):  # the short ticks halfway
        if length <= 2.5:
            positions.append(position)
            values.append(round(float(major(position)) / step) * step)
    y = Axis(tuple(positions), tuple(values))
    y.check()
    return y


def _panel(
    gas: str,
    index: int,
    frame: tuple[float, float, float, float],
    lines: list[np.ndarray],
    curves: list[np.ndarray],
    words: list[tuple[str, float, float]],
) -> Panel:
    """Return one panel, its axes calibrated and its labels checked."""
    left, bottom, right, top = frame
    ticks = [s for s in lines if 0.5 < np.hypot(*(s[1] - s[0])) < 4]
    y = _y_axis(f"{gas}", index, frame, ticks, words)
    # x: ticks with one end on the bottom edge; the panel's own point down out of it,
    # and the panel below draws its top ticks at the same places (the same values).
    on_x = _clusters(
        [
            float(s[0, 0])
            for s in ticks
            if abs(s[0, 0] - s[1, 0]) < EDGE
            and min(abs(s[:, 1] - bottom)) < EDGE
            and left - EDGE <= s[0, 0] <= right + EDGE
        ],
        gap=0.05,
    )
    if len(on_x) != 7:
        msg = f"{gas} panel {index}: x ticks {on_x}"
        raise ValueError(msg)
    x = Axis(tuple(on_x), (0, 1, 2, 3, 4, 5, 6))
    x.check()
    # the labels of the panel
    inside = [
        (s, wx, wy) for s, wx, wy in words if left < wx < right and bottom < wy < top
    ]
    paschen = _label(inside, "2p").removeprefix("2p")
    jk = _label(inside, f"{SHELL[gas]}p").removeprefix(SHELL[gas])
    if paschen != _ORDER[index] or jk != PANELS[gas][index]:
        msg = f"{gas} panel {index}: 2p{paschen} {jk}"
        raise ValueError(msg)
    mine = [
        c
        for c in curves
        if left < c[len(c) // 2, 0] < right and bottom < c[len(c) // 2, 1] < top
    ]
    if len(mine) != 1:
        msg = f"{gas} panel {index}: {len(mine)} curves"
        raise ValueError(msg)
    curve = mine[0]
    if abs(curve[0, 0] - left) > START:
        msg = f"{gas} panel {index}: the curve starts at {curve[0, 0]:.2f} pt"
        raise ValueError(msg)
    curve = curve[curve[:, 0] <= right + EDGE]
    if np.any(np.diff(curve[:, 0]) <= 0):
        msg = f"{gas} panel {index}: the curve steps back"
        raise ValueError(msg)
    return Panel(gas, paschen, jk, x, y, curve)


def file(kind: str, gas: str, jk: str) -> str:
    """Return a data file's name: kind fig, tot, s32 or s12."""
    return f"{gas}/{SOURCE}__{kind}_{gas}_{tag(jk)}.csv"


def _figures() -> dict[str, list[Point]]:
    """Return the totals of Figs. 3-5, one file per gas and level."""
    out = {}
    for gas in GASES:
        for panel in _panels(gas):
            base = threshold(gas, panel.jk)
            points = []
            for number, (px, py) in enumerate(panel.curve, start=1):
                omega = float(panel.x(px))
                points.append(
                    Point(
                        f"v{number:03d}",
                        f"{base + omega:.4f}",
                        f"{float(panel.y(py)):.3f}",
                        "",
                        f"vertex {number} at ({px:.3f}, {py:.3f}) pt, "
                        f"omega_1/2 = {omega:.4f} eV",
                    )
                )
            out[file("fig", gas, panel.jk)] = points
    return out


@functools.cache
def _table() -> list[list[str]]:
    """Return Table 2's ten rows of cells."""
    grid = json_table(DIRECTORY, TABLE_2)
    if len(grid) != 12:
        msg = f"Table 2 has {len(grid)} rows"
        raise ValueError(msg)
    return grid[2:]


def _value(cell: str) -> str:
    """Return a Table 2 cell's number, without footnote e."""
    text = cell.replace(r"^{\mathrm{e}}", "").strip()
    float(text)
    return text


def _tables() -> dict[str, list[Point]]:
    """Return Table 2's three cross sections, one file per gas, level and column."""
    out = {}
    for cells in _table():
        jk = _table_jk(cells[1])
        row = cells[0].replace("_{", "").replace("_", "").replace("}", "")
        for g, gas in enumerate(GASES):
            tot, s32, s12, _ = cells[2 + 4 * g : 6 + 4 * g]
            printed = " | ".join(cells[:2] + cells[2 + 4 * g : 6 + 4 * g])
            zero = f"{threshold(gas, jk):.4f}"
            fs = f"{threshold(gas, jk) + splitting(gas):.4f}"
            out[file("tot", gas, jk)] = [Point(row, zero, _value(tot), "", printed)]
            out[file("s32", gas, jk)] = [Point(row, zero, _value(s32), "", printed)]
            out[file("s12", gas, jk)] = [Point(row, fs, _value(s12), "", printed)]
    return out


def extract() -> dict[str, list[Point]]:
    """Return every file of the record."""
    return _tables() | _figures()


def checks() -> None:
    """Check the curves against Table 2, the ratio column and Delta_FS."""
    table = {(_table_jk(cells[1]), g): cells for cells in _table() for g in range(3)}
    worst = 0.0
    for g, gas in enumerate(GASES):
        if abs(splitting(gas) - DELTA_FS[gas]) > 0.0006:
            msg = f"{gas}: Delta_FS {splitting(gas):.4f} against {DELTA_FS[gas]}"
            raise ValueError(msg)
        for panel in _panels(gas):
            cells = table[panel.jk, g]
            drawn = float(panel.y(panel.curve[0, 1]))
            printed = float(_value(cells[2 + 4 * g]))
            worst = max(worst, abs(drawn - printed))
            if abs(drawn - printed) > AT_ZERO:
                msg = f"{gas} {panel.jk}: drawn {drawn:.3f} against {printed}"
                raise ValueError(msg)
    for cells in _table():
        for g in range(3):
            s32, s12, ratio = (cells[3 + 4 * g], cells[4 + 4 * g], cells[5 + 4 * g])
            low = (float(s12) - _half(s12)) / (float(s32) + _half(s32))
            high = (float(s12) + _half(s12)) / (float(s32) - _half(s32))
            if not low - _half(ratio) <= float(ratio) <= high + _half(ratio):
                msg = f"{cells[0]} gas {g}: r {ratio} against {s12}/{s32}"
                raise ValueError(msg)
    print(f"curves at omega_1/2 = 0 against Table 2: worst {worst:.3f} Mb")


def _half(cell: str) -> float:
    """Return half a unit of the last printed digit."""
    text = cell.strip()
    decimals = len(text.split(".")[1]) if "." in text else 0
    return 0.5 * 10.0**-decimals


def main() -> None:
    """Check the figures and the table, and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
