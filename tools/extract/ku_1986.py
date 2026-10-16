# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ku and Setser (1986): lifetimes and quenching of Xe 6p in Xe and Ar.

J. Chem. Phys. 84, 4304-4316 (1986); journal p. 4304 is PDF page 2 (PDF page 1 is the
AIP download cover), so the journal page is the PDF page plus 4302. A scanned (raster)
PDF: every cell of Tables I-III was checked against the table images and the rendered
pages (2026-10-09); the cells come from the library's ``.json``. Rows are the levels in
the paper's Paschen notation (``2p5``, ``1s2``). Taken, 'This work' only:

* **Table I** (PDF p. 5, journal p. 4307, ``tbl_e76f9b89f7ed``), in pure xenon: the
  radiative lifetimes (ns) and total quenching rate constants k_Q (10^-11 cm^3/s) of
  the six 6p levels; and the slope of the 'Slow component' after 2p10 excitation
  (footnote g), a coupled-pair rate, in a file of its own. The parenthesised values in
  the k_Q column (Refs. 24 and 7) and the row 1s3 (Ref. 23) are other papers'.
* **Table II** (PDF p. 7, journal p. 4309, ``tbl_c5ee4460f517``), in argon: k_Ar for
  2p6 and 2p9; for 2p5 the fast and slow decay components (footnotes b, c), whose
  cell the library runs together, '14.0 \pm 1.0^b(3.0 \pm 0.5)^c', split into the rows
  ``2p5 fast`` and ``2p5 slow`` of a separate file; and the radiative lifetimes from the
  argon intercepts (footnote a).
* **Table III** (PDF p. 8, journal p. 4310, ``tbl_0fc2938f26a2``), the state-to-state
  rate constants k_ij by xenon. Each reactant's processes and values are stacked lines
  in one cell, which the library runs together ('1.3 \pm 0.22.1 \pm 0.2...'); the
  script holds the lines as the table image shows them (:data:`LINES`) and checks that
  they join to the library's cells. One file per product (the channel), the row the
  reactant; the two upper limits '<0.2' (2p6 -> 2p10, 2p9 -> 1s2) are held as the
  interval 0-0.2 (value 0.1, bound half-width 0.1) in two more files. The Gamma column
  (k_ij / k_Q) and the tau column (the lifetimes used in the simulations) are not taken.

**Checks** (:func:`checks`; the script fails if one is missed): the counts; each
reactant's k_ij summed against Table I's k_Q (Table III's footnote d, 'Assigned by
difference', has no mark in the table; the '1s3 or Xe(6s)' rates close the sums
exactly); every Gamma against k_ij / k_Q (within 0.012); detailed balance at the stated
300 K with NIST's level energies: 1s2 -> 2p10 = (2p10 -> 1s2) / exp(84.10 cm^-1 / kT)
(6.68 against 6.7), the record's inline 3d5 -> 2p5 = (1/3) exp(-132.34 / kT) (2p5 ->
3d5) (1.77 against 1.8), and the 2p10/1s3 equilibrium constant '(1/57)' (1/57.1); and
Fig. 1's level energies (read from the rendered page, :data:`FIGURE_1`) against NIST's,
a uniform offset of 0.50-0.54 cm^-1 (the evidence of the level mapping).
"""

import csv
import math
import re
import tomllib

from _json_table import json_table
from _library import DATA, Point, number, write

from optmat.constants import HBAR, K_B, C

DIRECTORY = "Ku_1986_CollisionalDeactivationXe"
SOURCE = "ku_1986"
TABLE_I = "tbl_e76f9b89f7ed"
TABLE_II = "tbl_c5ee4460f517"
TABLE_III = "tbl_0fc2938f26a2"
_PASCHEN = re.compile(r"^([12][ps])_\{?(\d+)\}?$")
_VALUE = re.compile(r"^([\d.]+) \\pm ([\d.]+)(?:\^\{?[a-z]\}?)?$")
OTHER = "1s3 or Xe(6s)"
# Table III as the table image shows it: reactant -> (product, k_ij, Gamma) lines.
LINES = {
    "2p6": [
        ("2p7", "1.3 \\pm 0.2", "0.16 \\pm 0.03"),
        ("2p8", "2.1 \\pm 0.2", "0.26 \\pm 0.03"),
        ("2p9", "2.8 \\pm 0.2", "0.34 \\pm 0.03"),
        ("2p10", "<0.2", "<0.02"),
        (OTHER, "1.9 \\pm 0.2", "0.23 \\pm 0.03"),
    ],
    "2p7": [
        ("2p6", "0.6 \\pm 0.2", "0.05 \\pm 0.02"),
        ("2p8", "6.0 \\pm 2.0", "0.47 \\pm 0.16"),
        (OTHER, "6.2 \\pm 2.0", "0.48 \\pm 0.16"),
    ],
    "2p8": [
        ("2p9", "3.1 \\pm 0.2", "0.58 \\pm 0.04"),
        ("2p10", "0.4 \\pm 0.2", "0.08 \\pm 0.03"),
        ("1s2", "0.4 \\pm 0.2", "0.08 \\pm 0.03"),
        (OTHER, "1.4 \\pm 0.2", "0.26 \\pm 0.04"),
    ],
    "2p9": [
        ("2p8", "1.1 \\pm 0.2", "0.11 \\pm 0.02"),
        ("2p10", "4.5 \\pm 0.2", "0.47 \\pm 0.02"),
        ("1s2", "<0.2", "<0.02"),
        (OTHER, "4.0 \\pm 0.2", "0.42 \\pm 0.02"),
    ],
    "2p10": [
        ("1s2", "10.0 \\pm 2.0", "0.75 \\pm 0.15"),
        (OTHER, "3.3 \\pm 2.0", "0.25 \\pm 0.15"),
    ],
    "1s2": [
        ("2p10", "6.7 \\pm 0.7^{e}", "0.90 \\pm 0.10"),
        (OTHER, "0.70 \\pm 0.3^{f}", "0.10 \\pm 0.10"),
    ],
}
# Table III's product -> file stem.
CHANNELS = {
    "2p6": "k_2p6",
    "2p7": "k_2p7",
    "2p8": "k_2p8",
    "2p9": "k_2p9",
    "2p10": "k_2p10",
    "1s2": "k_1s2",
    OTHER: "k_1s3_6s",
}
# Fig. 1 (PDF p. 2), the level energies (cm^-1) as printed, and NIST's.
FIGURE_1 = {
    "2p5": ("80119.47", 80118.962),
    "2p6": ("79212.97", 79212.465),
    "2p7": ("78956.54", 78956.031),
    "2p8": ("78403.56", 78403.061),
    "2p9": ("78120.30", 78119.798),
    "2p10": ("77269.65", 77269.145),
    "1s2": ("77185.56", 77185.041),
    "1s3": ("76197.29", 76196.767),
    "3d3": ("80323.28", 80322.746),
    "3d4": ("80197.16", 80196.629),
    "3d5": ("79987.16", 79986.618),
    "3d6": ("79771.80", 79771.267),
}


def _level(cell: str) -> str:
    match = _PASCHEN.match(cell.strip())
    if match is None:
        msg = f"unexpected level {cell!r}"
        raise ValueError(msg)
    return "".join(match.groups())


def _value(cell: str) -> tuple[str, str]:
    match = _VALUE.match(cell.strip())
    if match is None:
        msg = f"unexpected value {cell!r}"
        raise ValueError(msg)
    return number(match.group(1)), match.group(2)


def _file(stem: str) -> str:
    return f"Xe/{SOURCE}__{stem}_Xe.csv"


def _table_i() -> dict[str, list[Point]]:
    """Return Table I's lifetimes, k_Q and the slow component's slope."""
    tau, kq, slow = [], [], []
    for cells in json_table(DIRECTORY, TABLE_I)[2:]:
        if cells[0].startswith("Slow component"):
            row = "slow component"
            slow.append(Point(row, row, *_value(cells[4]), " | ".join(cells)))
            continue
        if cells[0] == "1s_3":  # Ref. 23's values
            continue
        level = _level(cells[0])
        printed = " | ".join(cells)
        tau.append(Point(level, level, *_value(cells[2]), printed))
        own = cells[4].split("(")[0]  # the parenthesised value is another paper's
        kq.append(Point(level, level, *_value(own), printed))
    return {_file("tau"): tau, _file("kq_self"): kq, _file("kq_slow"): slow}


def _table_ii() -> dict[str, list[Point]]:
    """Return Table II's argon rates (single, coupled) and lifetimes."""
    single, coupled, tau = [], [], []
    for cells in json_table(DIRECTORY, TABLE_II)[2:]:
        level = _level(cells[0])
        printed = " | ".join(cells)
        tau.append(Point(level, level, *_value(cells[4]), printed))
        if level == "2p5":
            fast, slow = cells[2].split("(")
            if not slow.endswith(")^c"):
                msg = f"unexpected 2p5 cell {cells[2]!r}"
                raise ValueError(msg)
            fast = fast.removesuffix("^b")
            for row, cell in (("2p5 fast", fast), ("2p5 slow", slow[:-3])):
                coupled.append(Point(row, row, *_value(cell), printed))
        else:
            single.append(Point(level, level, *_value(cells[2]), printed))
    return {
        _file("kq_Ar"): single,
        _file("kq_Ar_coupled"): coupled,
        _file("tau_Ar"): tau,
    }


def _lines() -> dict[str, list[tuple[str, str, str]]]:
    """Return Table III's lines, checked against the library's joined cells."""
    for cells in json_table(DIRECTORY, TABLE_III)[2:]:
        if not cells[3]:  # 1s3: only Ref. 15's lifetime
            continue
        reactant = _level(cells[0])
        steps = [s.strip() for s in cells[3].split("\\rightarrow")]
        products = [OTHER if p == "1s_{3} or Xe(6s)" else _level(p) for p in steps[1:]]
        lines = LINES[reactant]
        if (
            _level(steps[0]) != reactant
            or products != [line[0] for line in lines]
            or "".join(line[1] for line in lines) != cells[4]
            or "".join(line[2] for line in lines) != cells[5]
        ):
            msg = f"Table III {reactant}: {cells}"
            raise ValueError(msg)
    return LINES


def _table_iii() -> dict[str, list[Point]]:
    """Return Table III's k_ij, one file per product, and the two bounds."""
    files: dict[str, list[Point]] = {_file(s): [] for s in CHANNELS.values()}
    files[_file("k_2p10_bound")] = []
    files[_file("k_1s2_bound")] = []
    for reactant, lines in _lines().items():
        for product, cell, _ in lines:
            printed = f"{reactant} -> {product}: {cell}"
            if cell == "<0.2":
                point = Point(reactant, reactant, "0.1", "0.1", printed)
                files[_file(f"{CHANNELS[product]}_bound")].append(point)
                continue
            point = Point(reactant, reactant, *_value(cell), printed)
            files[_file(CHANNELS[product])].append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return every file of the record."""
    return _table_i() | _table_ii() | _table_iii()


def _inline() -> dict[str, float]:
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        record = tomllib.load(handle)
    return {
        d["id"]: float(d["points"][0][2])
        for d in record["datasets"]
        if d.get("table") == "text"
    }


def _nist_energy(level: str) -> float:
    return FIGURE_1[level][1]


def _check_nist_file() -> None:
    """Check FIGURE_1's NIST energies against the held NIST level list."""
    held = set()
    path = DATA / "external" / "nist_asd" / "Xe_I_levels.tsv"
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            held.add(row["Level (cm-1)"].strip('"'))
    for level, (_, energy) in FIGURE_1.items():
        if f"{energy:.3f}" not in held:
            msg = f"NIST energy of {level}: {energy}"
            raise ValueError(msg)


def checks() -> None:
    """Check the counts, Table III against Table I, detailed balance and Fig. 1."""
    files = extract()
    counts = [len(points) for points in files.values()]
    if counts != [6, 6, 1, 2, 2, 3, 1, 1, 3, 2, 3, 2, 6, 1, 1]:
        msg = f"counts {counts}"
        raise ValueError(msg)
    kq = {p.row: float(p.value) for p in files[_file("kq_self")]}
    for reactant, lines in LINES.items():
        rates = [float(_value(c)[0]) for _, c, _ in lines if c != "<0.2"]
        total = kq.get(reactant, sum(rates))
        if abs(sum(rates) - total) > 0.15:
            msg = f"sum for {reactant}: {sum(rates)} against {total}"
            raise ValueError(msg)
        for _, cell, gamma in lines:
            if cell == "<0.2":
                continue
            ratio = float(_value(cell)[0]) / total
            if abs(ratio - float(_value(gamma)[0])) > 0.012:
                msg = f"Gamma of {reactant}: {ratio:.3f} against {gamma}"
                raise ValueError(msg)
    kt = K_B * 300.0 / (2 * math.pi * HBAR * C * 100)  # cm^-1
    rates = {p.row: float(p.value) for p in files[_file("k_2p10")]}
    up = {p.row: float(p.value) for p in files[_file("k_1s2")]}
    gap = _nist_energy("2p10") - _nist_energy("1s2")
    inline = _inline()
    gap_3d5 = _nist_energy("2p5") - _nist_energy("3d5")
    gap_1s3 = _nist_energy("2p10") - _nist_energy("1s3")
    tests = {
        "1s2 -> 2p10": (up["2p10"] / math.exp(gap / kt), rates["1s2"], 0.05),
        "3d5 -> 2p5": (
            inline["k1a_Ar"] * math.exp(-gap_3d5 / kt) / 3,
            inline["k2a_Ar"],
            0.05,
        ),
        "1/57": (1 / (3 * math.exp(-gap_1s3 / kt)), 57.0, 0.5),
    }
    for name, (computed, printed, tolerance) in tests.items():
        if abs(computed - printed) > tolerance:
            msg = f"{name}: {computed:.3f} against {printed}"
            raise ValueError(msg)
    _check_nist_file()
    for level, (printed, nist) in FIGURE_1.items():
        if abs(float(printed) - nist - 0.52) > 0.03:
            msg = f"Fig. 1 {level}: {printed} against NIST {nist}"
            raise ValueError(msg)


def main() -> None:
    """Check the tables and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
