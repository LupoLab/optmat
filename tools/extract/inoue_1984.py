# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Inoue, Ku and Setser (1984): lifetimes and quenching of Xe 6p, 6p', 7p, 6d.

J. Chem. Phys. 81, 5760-5774 (1984); journal p. 5760 is PDF page 2 (PDF page 1 is the
AIP download cover), so the journal page is the PDF page plus 5758. A scanned (raster)
PDF with an OCR text layer: every taken cell of Tables I and II was checked against the
table images and the rendered pages (2026-10-09); the cells come from the library's
``.json``. Rows are the levels in the paper's Paschen notation (``2p1``, ``4d1'``,
``4d1''``), as both tables print them; Table I's Racah column is checked against the
NIST mapping below (:data:`LEVELS`). Taken, 'This work' only:

* **Table I** (PDF p. 3, journal p. 5761, ``tbl_19e693573e85``), radiative lifetimes
  (ns), the zero-pressure intercepts of the Stern-Volmer plots in neon and in argon
  buffer gas. One file per buffer; the three values the authors call the less
  reliable of a pair (2p4 and 3p9 in neon, 3p10 in argon; pp. 5765, 5769-5770) in two
  more files.
* **Table II** (PDF p. 5, journal p. 5763, ``tbl_e8107e3f6246``), total quenching rate
  constants (10^-11 cm^3 atom^-1 s^-1) by neon and by argon, one file each; and by
  xenon, 2p10's '13 \pm 2' (the fast component, Fig. 10b) in a file of its own. The
  parenthesised values (Refs. 4, 8 and 10) and 2p5's and 2p6's xenon values (footnote
  d, the note added in proof: Setser and Ku's two-photon measurement, ku_1986's
  Table I) are other work, not taken.

**Checks** (:func:`checks`; the script fails if one is missed): the counts; Table I's
Racah labels and both tables' rows against :data:`LEVELS`, whose NIST energies are
checked in the held level list; Fig. 2's printed energies (read from the rendered page,
:data:`FIGURE_2`) against NIST's, a uniform offset of 0.43-0.57 cm^-1; the energy gaps
the text prints against NIST's (:data:`GAPS`); the equilibrium constants K(4d1'/2p1) =
33, K(4d1'/2p2) = 1.5 and K(2p10/1s3) = 57 at 300 K; the air wavelengths the text
prints (828.01, 820.63, 823.16, 450.1, 452.5 nm) against NIST's vacuum lines; every
value the text repeats (3p10's lifetimes and rates, 2p10's, 4d1''s '88 ns', the
conclusions' ranges and ratios); the average discrepancy '1.7 ns' of the 6p and 6p'
lifetimes from Aymar and Coulombe's tau_v (Table I's theoretical column); and the
record's NIST lifetimes, 1/sum A over the NIST lines from each level (:data:`NIST_TAU`).
"""

import csv
import math
import re
import statistics

from _json_table import json_table
from _library import DATA, Point, number, write

from optmat.constants import HBAR, K_B, C

DIRECTORY = "Inoue_1984_LaserInducedFluorescence"
SOURCE = "inoue_1984"
TABLE_I = "tbl_19e693573e85"
TABLE_II = "tbl_e8107e3f6246"
_PASCHEN = re.compile(r"^(\d[spd])('*)_\{?(\d+)\}?('*)$")
_VALUE = re.compile(r"^([\d.]+) (?:±|\\pm) ([\d.]+)$")
_MISSING = "..."
P32 = "5p5.(2P*<3/2>)."
P12 = "5p5.(2P*<1/2>)."
# Paschen -> (Table I's Racah, NIST configuration, term, J, NIST energy in cm^-1).
LEVELS = {
    "2p1": ("6p'[1/2]_0", P12 + "6p", "2[1/2]", "0", "89860.015"),
    "2p2": ("6p'[1/2]_1", P12 + "6p", "2[1/2]", "1", "89278.706"),
    "2p3": ("6p'[3/2]_2", P12 + "6p", "2[3/2]", "2", "89162.356"),
    "2p4": ("6p'[3/2]_1", P12 + "6p", "2[3/2]", "1", "88379.126"),
    "2p5": ("6p[1/2]_0", P32 + "6p", "2[1/2]", "0", "80118.962"),
    "2p6": ("6p[3/2]_2", P32 + "6p", "2[3/2]", "2", "79212.465"),
    "2p7": ("6p[3/2]_1", P32 + "6p", "2[3/2]", "1", "78956.031"),
    "2p8": ("6p[5/2]_3", P32 + "6p", "2[5/2]", "3", "78403.061"),
    "2p9": ("6p[5/2]_2", P32 + "6p", "2[5/2]", "2", "78119.798"),
    "2p10": ("6p[1/2]_1", P32 + "6p", "2[1/2]", "1", "77269.145"),
    "3p5": ("7p[1/2]_0", P32 + "7p", "2[1/2]", "0", "88842.256"),
    "3p7": ("7p[3/2]_1", P32 + "7p", "2[3/2]", "1", "88744.559"),
    "3p6": ("7p[3/2]_2", P32 + "7p", "2[3/2]", "2", "88686.500"),
    "3p8": ("7p[5/2]_3", P32 + "7p", "2[5/2]", "3", "88469.213"),
    "3p9": ("7p[5/2]_2", P32 + "7p", "2[5/2]", "2", "88351.681"),
    "3p10": ("7p[1/2]_1", P32 + "7p", "2[1/2]", "1", "87927.131"),
    "4d2": ("6d[3/2]_1", P32 + "6d", "2[3/2]*", "1", "90032.155"),
    "4d1'": ("6d[5/2]_3", P32 + "6d", "2[5/2]*", "3", "89534.568"),
    "4d1''": ("6d[5/2]_2", P32 + "6d", "2[5/2]*", "2", "89243.258"),
    "4d4": ("6d[7/2]_3", P32 + "6d", "2[7/2]*", "3", "89024.890"),
    "4d4'": ("6d[7/2]_4", P32 + "6d", "2[7/2]*", "4", "88911.692"),
    "4d3": ("6d[3/2]_2", P32 + "6d", "2[3/2]*", "2", "88708.466"),
    "4d5": ("6d[1/2]_1", P32 + "6d", "2[1/2]*", "1", "88549.775"),
    "4d6": ("6d[1/2]_0", P32 + "6d", "2[1/2]*", "0", "88491.020"),
}
# Lower levels the text names (NIST energies, cm^-1).
LOWER = {
    "1s5": "67067.547",  # 6s[3/2]2
    "1s4": "68045.156",  # 6s[3/2]1
    "1s3": "76196.767",  # 6s'[1/2]0
    "1s2": "77185.041",  # 6s'[1/2]1
    "2s4": "85440.017",  # 7s[3/2]1
}
# Fig. 2 (PDF p. 7, journal p. 5765), the energies (cm^-1) as printed.
FIGURE_2 = {
    "2p1": "89860.5",
    "2p2": "89279.2",
    "2p3": "89162.9",
    "2p4": "88379.6",
    "3p5": "88842.8",
    "3p7": "88745.1",
    "3p6": "88687.0",
    "3p8": "88469.7",
    "3p9": "88352.2",
    "3p10": "87927.7",
    "4d2": "90032.7",
    "4d1'": "89535.0",
    "4d1''": "89243.7",
    "4d4": "89025.4",
    "4d4'": "88912.2",
    "4d3": "88709.0",
    "4d5": "88550.3",
    "4d6": "88491.5",
}
# Gaps the text prints (cm^-1): (upper, lower, printed, tolerance).
GAPS = [
    ("2p1", "2p2", "581", 0.5),  # p. 5765
    ("4d2", "2p1", "172", 0.5),  # p. 5765
    ("2p1", "4d1'", "325", 0.5),  # p. 5765
    ("3p10", "3p9", "-425", 0.5),  # p. 5765
    ("3p10", "2s4", "2487", 0.5),  # p. 5765
    ("4d4'", "3p5", "69.4", 0.05),  # p. 5768
    ("4d6", "3p8", "22", 0.5),  # p. 5768
    ("3p7", "3p6", "58", 0.1),  # p. 5768, 'less than 58'
    ("2p10", "1s3", "1073", 0.7),  # p. 5771
]
# Air wavelengths the text prints (nm): (upper, lower, printed, NIST vacuum line).
AIR = [
    ("2p5", "1s4", "828.01", "828.23921"),  # pp. 5770, 5771
    ("2p4", "1s3", "820.63", "820.8592"),  # p. 5771, Fig. 11
    ("2p6", "1s5", "823.16", "823.38964"),  # p. 5771, Fig. 11
    ("2p2", "1s5", "450.1", "450.2241"),  # p. 5771
    ("2p3", "1s5", "452.5", "452.59493"),  # p. 5771
]
# NIST 1/sum A (ns) over the NIST lines from each level (both line files), as the
# record's notes quote them; None where NIST gives no transition probability.
NIST_TAU = {
    "2p1": 27.5,
    "2p2": 25.9,
    "2p3": 22.8,
    "2p4": 22.8,
    "2p5": 27.1,
    "2p6": 27.2,
    "2p7": 37.0,
    "2p8": 33.3,
    "2p9": 39.4,
    "2p10": 30.4,
    "3p5": 69.0,
    "3p7": 100.8,
    "3p6": 106.4,
    "3p8": 123.7,
    "3p9": 129.2,
    "3p10": 109.1,
    "4d2": 6.7,
    "4d1'": None,
    "4d1''": None,
    "4d4": None,
    "4d4'": None,
    "4d3": None,
    "4d5": 227.3,
    "4d6": None,
}
# Table I's lifetimes the authors call the less reliable of the pair (pp. 5765, 5770).
LESS_RELIABLE = {"Ne": ("2p4", "3p9"), "Ar": ("3p10",)}


def _level(cell: str) -> str:
    match = _PASCHEN.match(cell.strip())
    if match is None:
        msg = f"unexpected level {cell!r}"
        raise ValueError(msg)
    base, before, digits, after = match.groups()
    return f"{base}{digits}{before}{after}"


def _value(cell: str) -> tuple[str, str]:
    match = _VALUE.match(cell.strip())
    if match is None:
        msg = f"unexpected value {cell!r}"
        raise ValueError(msg)
    return number(match.group(1)), number(match.group(2))


def _file(stem: str) -> str:
    return f"Xe/{SOURCE}__{stem}_Xe.csv"


def _table_i() -> dict[str, list[Point]]:
    """Return Table I's lifetimes in neon and argon, the less reliable apart."""
    files: dict[str, list[Point]] = {
        _file("tau_Ne"): [],
        _file("tau_Ar"): [],
        _file("tau_Ne_check"): [],
        _file("tau_Ar_check"): [],
    }
    for cells in json_table(DIRECTORY, TABLE_I)[2:]:
        level = _level(cells[0])
        if cells[1] != LEVELS[level][0]:
            msg = f"Table I Racah of {level}: {cells[1]!r}"
            raise ValueError(msg)
        printed = " | ".join(cells[:4])
        for buffer, cell in (("Ne", cells[2]), ("Ar", cells[3])):
            if cell == _MISSING:
                continue
            stem = f"tau_{buffer}"
            if level in LESS_RELIABLE[buffer]:
                stem += "_check"
            files[_file(stem)].append(Point(level, level, *_value(cell), printed))
    return files


def _table_ii() -> dict[str, list[Point]]:
    """Return Table II's total quenching by neon, argon and (2p10 only) xenon."""
    ne, ar, xe = [], [], []
    for cells in json_table(DIRECTORY, TABLE_II)[1:]:
        level = _level(cells[0])
        printed = " | ".join(cells)
        for target, cell in ((ne, cells[1]), (ar, cells[2])):
            if cell != _MISSING:
                target.append(Point(level, level, *_value(cell), printed))
        if level == "2p10":
            own, other = cells[4].split("(")
            if other != "41 \\pm 5)^{c}":
                msg = f"unexpected 2p10 cell {cells[4]!r}"
                raise ValueError(msg)
            xe.append(Point(level, level, *_value(own), printed))
    return {_file("kq_Ne"): ne, _file("kq_Ar"): ar, _file("kq_self"): xe}


def extract() -> dict[str, list[Point]]:
    """Return every file of the record."""
    return _table_i() | _table_ii()


def _held_levels() -> set[str]:
    held = set()
    path = DATA / "external" / "nist_asd" / "Xe_I_levels.tsv"
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["Configuration"], row["Term"], row["J"], row["Level (cm-1)"])
            held.add("|".join(cell.strip('"') for cell in key))
    return held


def _nist_lines() -> list[dict[str, str]]:
    rows = []
    for name in ("Xe_I_lines.tsv", "Xe_I_lines_nir.tsv"):
        path = DATA / "external" / "nist_asd" / name
        with path.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows.append({k: (v or "").strip('"') for k, v in row.items() if k})
    return rows


def _energy(level: str) -> float:
    return float(LEVELS[level][4] if level in LEVELS else LOWER[level])


def _check_levels() -> None:
    """Check the mapping, Fig. 2, the gaps and the equilibrium constants."""
    held = _held_levels()
    for level, (_, configuration, term, j, energy) in LEVELS.items():
        if f"{configuration}|{term}|{j}|{energy}" not in held:
            msg = f"NIST level of {level}"
            raise ValueError(msg)
    for level, printed in FIGURE_2.items():
        if abs(float(printed) - _energy(level) - 0.50) > 0.075:
            msg = f"Fig. 2 {level}: {printed} against NIST {_energy(level)}"
            raise ValueError(msg)
    for upper, lower, printed, tolerance in GAPS:
        gap = _energy(upper) - _energy(lower)
        if abs(gap - float(printed)) > tolerance:
            msg = f"gap {upper} - {lower}: {gap:.2f} against {printed}"
            raise ValueError(msg)
    kt = K_B * 300.0 / (2 * math.pi * HBAR * C * 100)  # cm^-1
    constants = {
        # (g(4d1') / g(2p1)) exp(gap / kT)
        "K(4d1'/2p1)": (7 * math.exp((_energy("2p1") - _energy("4d1'")) / kt), 33.0),
        # printed as the inverse: (g(2p2) / g(4d1')) exp(gap / kT)
        "K(4d1'/2p2)": (3 / 7 * math.exp((_energy("4d1'") - _energy("2p2")) / kt), 1.5),
        "K(2p10/1s3)": (math.exp((_energy("2p10") - _energy("1s3")) / kt) / 3, 57.0),
    }
    for name, (computed, printed) in constants.items():
        if abs(computed / printed - 1) > 0.04:
            msg = f"{name}: {computed:.2f} against {printed}"
            raise ValueError(msg)


def _check_wavelengths(lines: list[dict[str, str]]) -> None:
    """Check the printed air wavelengths against NIST's vacuum lines."""
    for upper, lower, printed, vacuum in AIR:
        found = [
            r
            for r in lines
            if r["Ek(cm-1)"] == LEVELS[upper][4] and r["Ei(cm-1)"] == LOWER[lower]
        ]
        if len(found) != 1 or found[0]["obs_wl_vac(nm)"] != vacuum:
            msg = f"NIST line {upper} -> {lower}"
            raise ValueError(msg)
        ratio = float(vacuum) / float(printed) - 1
        if not 2.0e-4 < ratio < 3.4e-4:  # the air-vacuum difference, to the digits
            msg = f"{upper} -> {lower}: {printed} against vacuum {vacuum}"
            raise ValueError(msg)
    # p. 5769 prints '3p8 -> 1s5 transition at 452.7 nm' (record doubt): NIST's line
    # from 3p8 to 1s5 is 467.25 nm (vacuum); no line from 3p8 lies near 452.7 nm.
    from_3p8 = [r for r in lines if r["Ek(cm-1)"] == LEVELS["3p8"][4]]
    to_1s5 = [r for r in from_3p8 if r["Ei(cm-1)"] == LOWER["1s5"]]
    near = [r for r in from_3p8 if abs(float(r["obs_wl_vac(nm)"] or 0) - 452.8) < 2]
    if [r["obs_wl_vac(nm)"] for r in to_1s5] != ["467.25334"] or near:
        msg = "the 3p8 -> 1s5 line"
        raise ValueError(msg)


def _check_nist_lifetimes(lines: list[dict[str, str]]) -> None:
    """Check the record's NIST 1/sum A against the held line lists."""
    for level, quoted in NIST_TAU.items():
        rates = [
            float(r["Aki(s^-1)"])
            for r in lines
            if r["Ek(cm-1)"] == LEVELS[level][4] and r["Aki(s^-1)"]
        ]
        if not rates:
            if quoted is not None:
                msg = f"NIST lifetime of {level}: no A, quoted {quoted}"
                raise ValueError(msg)
            continue
        tau = 1e9 / sum(rates)
        if quoted is None or abs(tau - quoted) > 0.05:
            msg = f"NIST lifetime of {level}: {tau:.2f} against {quoted}"
            raise ValueError(msg)


def _check_text(files: dict[str, list[Point]]) -> None:
    """Check the values the text repeats against the tables."""

    def value(stem: str, level: str) -> float:
        for point in files[_file(stem)]:
            if point.row == level:
                return float(point.value)
        msg = f"{stem} has no {level}"
        raise KeyError(msg)

    tests = {
        "tau(3p10) Ne 123": (value("tau_Ne", "3p10"), 123.0),  # p. 5765
        "tau(3p10) Ar 131": (value("tau_Ar_check", "3p10"), 131.0),
        "k_Q(3p10) Ne 7e-12": (value("kq_Ne", "3p10"), 0.7),
        "k_Q(3p10) Ar 2.1e-10": (value("kq_Ar", "3p10"), 21.0),
        "tau(4d1') '88 ns'": (
            (value("tau_Ne", "4d1'") + value("tau_Ar", "4d1'")) / 2,
            88.0,
        ),  # p. 5765
        "tau(2p2) '38 ns'": (value("tau_Ar", "2p2"), 38.0),  # p. 5765
        "tau(3p8) '154 ns'": (value("tau_Ar", "3p8"), 154.0),  # p. 5769
        "k_Q(2p10) Ar 0.3e-11": (value("kq_Ar", "2p10"), 0.3),  # p. 5770
        "k_Q(2p10) Ne 0.3e-11": (value("kq_Ne", "2p10"), 0.3),
        "k(2p10) Xe 1.3e-10": (value("kq_self", "2p10"), 13.0),
        "k_Q(4d6) Ar 52": (value("kq_Ar", "4d6"), 52.0),  # p. 5773
        "k_Q(4d6) Ne 71": (value("kq_Ne", "4d6"), 71.0),
        "3p10 Ar/Ne '30 times'": (
            value("kq_Ar", "3p10") / value("kq_Ne", "3p10"),
            30.0,
        ),  # p. 5765
        "4d4 Ar/Ne 'four times'": (
            round(value("kq_Ar", "4d4") / value("kq_Ne", "4d4")),
            4.0,
        ),  # p. 5767
    }
    for name, (computed, printed) in tests.items():
        if abs(computed - printed) > 1e-9:
            msg = f"{name}: {computed}"
            raise ValueError(msg)
    # p. 5767: 'the present k_Q(2p2) is 25% lower' than Ref. 4's (15 +- 1).
    if abs(1 - value("kq_Ar", "2p2") / 15.0 - 0.25) > 0.03:
        msg = "k_Q(2p2) against Ref. 4"
        raise ValueError(msg)
    # p. 5773: in argon 'from 30 to 0.3' apart from 4d6.
    argon = [float(p.value) for p in files[_file("kq_Ar")] if p.row != "4d6"]
    if min(argon) != 0.3 or not 25 <= max(argon) <= 30:
        msg = f"argon range {min(argon)}-{max(argon)}"
        raise ValueError(msg)


def _check_tau_v() -> None:
    """Check the average discrepancy of 1.7 ns from tau_v for the 6p and 6p' levels.

    p. 5770: 'The average discrepancies are 1.7 and 3 ns for the Xe(6p and 6p') and
    Xe(7p) states'. All the printed 6p and 6p' values (neon and argon) give 1.69 ns;
    the 7p '3 ns' is not reproduced (3.6-5.5 ns, whichever values are taken).
    """
    differences = []
    for cells in json_table(DIRECTORY, TABLE_I)[2:]:
        if not _level(cells[0]).startswith("2p"):
            continue
        tau_v = float(cells[5])
        differences += [
            abs(float(_value(c)[0]) - tau_v) for c in cells[2:4] if c != _MISSING
        ]
    if abs(statistics.mean(differences) - 1.7) > 0.05:
        msg = f"6p, 6p' discrepancy {statistics.mean(differences):.2f}"
        raise ValueError(msg)


def checks() -> None:
    """Check the counts, the level mapping, the text, tau_v and NIST's lifetimes."""
    files = extract()
    counts = [len(points) for points in files.values()]
    if counts != [11, 19, 2, 1, 22, 20, 1]:
        msg = f"counts {counts}"
        raise ValueError(msg)
    rows_i = [_level(c[0]) for c in json_table(DIRECTORY, TABLE_I)[2:]]
    rows_ii = [_level(c[0]) for c in json_table(DIRECTORY, TABLE_II)[1:]]
    if rows_i != list(LEVELS) or rows_ii != list(LEVELS):
        msg = "the tables' rows"
        raise ValueError(msg)
    _check_levels()
    lines = _nist_lines()
    _check_wavelengths(lines)
    _check_nist_lifetimes(lines)
    _check_text(files)
    _check_tau_v()


def main() -> None:
    """Check the tables and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
