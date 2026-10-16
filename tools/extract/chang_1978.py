# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chang and Setser (1978): lifetimes and quenching by argon of Ar 4p, 4p'.

J. Chem. Phys. 69, 3885-3897 (1978); journal p. 3885 is PDF page 2 (PDF page 1 is the
AIP download cover), so the journal page is the PDF page plus 3883. The tables are read
from the library's cells (``.json``; merged headers), every taken cell checked against
the table images and the rendered pages (2026-10-12). Rows are the levels in the
paper's Paschen notation (2p2 ... 2p10; the record's notes give the NIST levels), and
``x`` the same label. Taken, 'Present work' / 'CS' only:

* **Table II** (PDF p. 6, journal p. 3889, ``tbl_56408e2b767a``), radiative lifetimes
  (ns) of the eight levels laser-excited from 1s5 (2p1 and 2p5 were not studied).
* **Table III** (PDF p. 6, ``tbl_dda0b9a86083``), deactivation by argon, 10^-11
  cm^3 atom^-1 s^-1: the total k_Q with its +- (seven levels), the upper limit of 2p6
  ('<=1.8', footnote b; held as the interval 0-1.8, value 0.9 with half-width 0.9,
  ``printed`` keeps it), the state-to-state k_nn' one file per product level (rows the
  initial level), the two values 'assigned using detailed balance' (footnote c, k43 and
  k98) apart, and the intermultiplet sums to the 4s, 4s' manifold that are numbers
  (k_Q minus the k_nn', so derived). Not taken: k_{9,10} '> 1.8' (a lower limit), the
  intermultiplet '<3.0' and '<3.5' and 2p6's 0.0 (differences with a bound), the 0.0
  of k_nn' for 2p2 and 2p10 (no product observed), sigma = k_Q / <v> (derived).
* **Table IV** (PDF p. 9, journal p. 3892, ``tbl_ea55ed0bd3cf``), column CS: the
  transition probabilities A_ki (10^6 s^-1) of the 26 lines of the eight levels (the
  product of the measured branching ratio and the inverse lifetime); the parenthesised
  entries of 2p1 and 2p5 are branching ratios, not A, and are not taken; the f column
  is derived from A (footnote a) and is not taken.

**Checks** (:func:`checks`; the script fails if one is missed): the counts; every
row's Paschen label against its printed Racah label (NIST's levels,
``data/external/nist_asd/Ar_I_levels.tsv``); Table III's sums k_Q = sum k_nn' + sum
k_nm' and sigma = k_Q / <v> with <v> the mean relative Ar-Ar speed at 300 K; the
detailed-balance values k43 and k98 from k34 and k89 with NIST's energies at 300 K;
Table IV's wavelengths against NIST's levels (standard air, within the 1 Angstrom of
the print), g_k and g_i against the levels' J, sum A x tau = 1 for each level, and the
CS f against f = 1.4992e-16 lambda^2 (g_k / g_i) A, which holds within the rounding
and 0.6 % for every line but 2p10-1s2 (printed 0.028 for 0.0028) and 2p10-1s5 (0.134
for 0.1318, 1.7 %); f is not taken.
"""

import csv
import math
import re

from _json_table import json_table
from _library import DATA, Point, number, write

from optmat.conditions import standard_air_refractivity
from optmat.constants import ATOMIC_MASS, HBAR, K_B, C

DIRECTORY = "Chang_1978_RadiativeLifetimesTwo"
SOURCE = "chang_1978"
TABLE_II = "tbl_56408e2b767a"
TABLE_III = "tbl_dda0b9a86083"
TABLE_IV = "tbl_ea55ed0bd3cf"
_CORE = {"3/2": "3s2.3p5.(2P*<3/2>)", "1/2": "3s2.3p5.(2P*<1/2>)"}
PASCHEN = {
    "1s5": ("4s", "3/2", "3/2", "2"),
    "1s4": ("4s", "3/2", "3/2", "1"),
    "1s3": ("4s", "1/2", "1/2", "0"),
    "1s2": ("4s", "1/2", "1/2", "1"),
    "2p1": ("4p", "1/2", "1/2", "0"),
    "2p2": ("4p", "1/2", "1/2", "1"),
    "2p3": ("4p", "1/2", "3/2", "2"),
    "2p4": ("4p", "1/2", "3/2", "1"),
    "2p5": ("4p", "3/2", "1/2", "0"),
    "2p6": ("4p", "3/2", "3/2", "2"),
    "2p7": ("4p", "3/2", "3/2", "1"),
    "2p8": ("4p", "3/2", "5/2", "2"),
    "2p9": ("4p", "3/2", "5/2", "3"),
    "2p10": ("4p", "3/2", "1/2", "1"),
}
"""Paschen label -> (orbital, core j, K, J): NIST's (2P*<core>).orbital 2[K] J."""
_PASCHEN = re.compile(r"^(\d)p_\{?(\d+)\}?$")
_RACAH = re.compile(r"^4p ?(')? ?\[(\d/2)\]_\{(\d)\}$")
_PM = re.compile(r"^([\d.]+)\((?:\\pm|±)([\d.]+)\)$")
_KNN = re.compile(r"k_\{(\d)(\d+)\} ([=>]) ([\d.]+)(\^\{\\text\{c\}\})?")
_LINE = re.compile(r"^2p_\{?(\d+)\}?-1s_(\d)$")
ANGSTROM2 = 1e-16
"""1 Angstrom^2 in cm^2."""
KT_300K_CM = K_B * 300.0 / (2 * math.pi * HBAR * C * 100)
"""kT at 300 K in cm^-1."""


def _level(cell: str) -> str:
    """Return a Paschen label, '2p_{10}' -> '2p10'."""
    match = _PASCHEN.match(cell.strip())
    if match is None:
        msg = f"unexpected level {cell!r}"
        raise ValueError(msg)
    return f"{match.group(1)}p{match.group(2)}"


def _racah_matches(level: str, cell: str) -> bool:
    """Return whether a printed Racah label, "4p' [1/2]_{0}", is the Paschen level."""
    match = _RACAH.match(cell.strip())
    if match is None:
        return False
    core = "1/2" if match.group(1) else "3/2"
    return PASCHEN[level] == ("4p", core, match.group(2), match.group(3))


def _pm(cell: str) -> tuple[str, str]:
    r"""Return the value and the +- of a cell '26.5(\pm0.4)'."""
    match = _PM.match(cell.strip())
    if match is None:
        msg = f"unexpected value {cell!r}"
        raise ValueError(msg)
    return match.group(1), match.group(2)


def _lifetimes() -> dict[str, list[Point]]:
    """Return Table II's 'Present work' lifetimes."""
    grid = json_table(DIRECTORY, TABLE_II)
    points = []
    for cells in grid[1:]:
        if not cells[2]:
            continue
        level = _level(cells[0])
        value, sigma = _pm(cells[2])
        points.append(Point(level, level, value, sigma, " | ".join(cells[:3])))
    return {f"Ar/{SOURCE}__tau_Ar.csv": points}


def _state_to_state(level: str, cell: str, printed: str) -> dict[str, list[Point]]:
    """Return the k_nn' of a Table III cell, by file (product, detailed balance)."""
    out: dict[str, list[Point]] = {}
    for match in _KNN.finditer(cell):
        initial, final, relation, value, balance = match.groups()
        if f"2p{initial}" != level:
            msg = f"{level}: {match.group(0)!r}"
            raise ValueError(msg)
        if relation == ">":
            continue  # k_{9,10} > 1.8: a lower limit, not taken
        name = f"k_2p{final}" + ("_db" if balance else "")
        point = Point(level, level, number(value), "", printed)
        out.setdefault(f"Ar/{SOURCE}__{name}.csv", []).append(point)
    return out


def _quenching() -> dict[str, list[Point]]:
    """Return Table III: k_Q, its upper limit, the k_nn' and the 4s sums."""
    grid = json_table(DIRECTORY, TABLE_III)
    out: dict[str, list[Point]] = {}

    def add(name: str, point: Point) -> None:
        out.setdefault(f"Ar/{SOURCE}__{name}.csv", []).append(point)

    for cells in grid[1:]:
        if not cells[4]:
            continue  # 2p1 and 2p5: no data
        level = _level(cells[0])
        printed = " | ".join(cells)
        if cells[4] == "≤1.8b":
            add("kq_self_bound", Point(level, level, "0.9", "0.9", printed))
        else:
            value, sigma = _pm(cells[4])
            add("kq_self", Point(level, level, value, sigma, printed))
        for name, points in _state_to_state(level, cells[2], printed).items():
            out.setdefault(name, []).extend(points)
        intermultiplet = cells[3].removesuffix("e")
        if re.fullmatch(r"[\d.]+", intermultiplet) and float(intermultiplet) > 0:
            add("k_4s", Point(level, level, intermultiplet, "", printed))
    return dict(sorted(out.items()))


def _transition_probabilities() -> dict[str, list[Point]]:
    """Return Table IV's CS transition probabilities (the numbers, not the ratios)."""
    grid = json_table(DIRECTORY, TABLE_IV)
    points = []
    for cells in grid[2:]:
        if cells[4].startswith("("):
            continue  # 2p1, 2p5: branching ratios
        upper, lower = _line(cells[0])
        label = f"{upper}-{lower}"
        points.append(Point(label, label, number(cells[4]), "", " | ".join(cells)))
    return {f"Ar/{SOURCE}__A_4p_4s.csv": points}


def _line(cell: str) -> tuple[str, str]:
    """Return the upper and lower Paschen levels of '2p_{10}-1s_2'."""
    match = _LINE.match(cell.strip())
    if match is None:
        msg = f"unexpected transition {cell!r}"
        raise ValueError(msg)
    return f"2p{match.group(1)}", f"1s{match.group(2)}"


def extract() -> dict[str, list[Point]]:
    """Return every file of the record."""
    return _lifetimes() | _quenching() | _transition_probabilities()


def nist_levels() -> dict[str, float]:
    """Return NIST's energies (cm^-1) of the Paschen levels 1s2-1s5, 2p1-2p10."""
    path = DATA / "external" / "nist_asd" / "Ar_I_levels.tsv"
    keys = {
        (f"{_CORE[core]}.{orbital}", f"2[{k}]" + ("*" if orbital == "4s" else ""), j): (
            level
        )
        for level, (orbital, core, k, j) in PASCHEN.items()
    }
    found = {}
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (
                row["Configuration"].strip('"'),
                row["Term"].strip('"'),
                row["J"].strip('"'),
            )
            if key in keys:
                found[keys[key]] = float(row["Level (cm-1)"].strip('"'))
    if len(found) != len(PASCHEN):
        msg = f"NIST levels {found}"
        raise ValueError(msg)
    return found


def balance(initial: str, final: str, energy: dict[str, float]) -> float:
    """Return k(initial -> final) / k(final -> initial) at 300 K (detailed balance)."""
    g_initial = 2 * int(PASCHEN[initial][3]) + 1
    g_final = 2 * int(PASCHEN[final][3]) + 1
    gap = energy[initial] - energy[final]
    return g_final / g_initial * math.exp(gap / KT_300K_CM)


def _check_labels() -> None:
    """Check each Paschen label of Tables II and III against its Racah label."""
    for table_id in (TABLE_II, TABLE_III):
        for cells in json_table(DIRECTORY, table_id)[1:]:
            if not _racah_matches(_level(cells[0]), cells[1]):
                msg = f"{table_id}: {cells[0]!r} is not {cells[1]!r}"
                raise ValueError(msg)


def _check_table_iii(files: dict[str, list[Point]]) -> None:
    """Check k_Q = sum k_nn' + sum k_nm', sigma = k_Q / <v>, and detailed balance."""
    grid = json_table(DIRECTORY, TABLE_III)
    mass = 39.948 / 2 * ATOMIC_MASS
    speed = math.sqrt(8 * K_B * 300.0 / (math.pi * mass)) * 100  # cm/s
    knn: dict[str, float] = {}
    for name, points in files.items():
        if "__k_2p" in name:
            for p in points:
                knn[p.row] = knn.get(p.row, 0.0) + float(p.value)
    knn["2p9"] = knn["2p9"] + 1.8  # k_{9,10} > 1.8 enters the sum at its bound
    for cells in grid[1:]:
        if not cells[4]:
            continue
        level = _level(cells[0])
        total = float(re.sub(r"[^\d.]", "", cells[4].split("(")[0].replace("b", "")))
        intermultiplet = float(re.sub(r"[^\d.]", "", cells[3].removesuffix("e")))
        if abs(knn.get(level, 0.0) + intermultiplet - total) > 0.051:
            msg = f"{level}: sum {knn.get(level, 0.0)} + {intermultiplet} != {total}"
            raise ValueError(msg)
        sigma = float(cells[5].lstrip("≤"))
        model = total * 1e-11 / speed / ANGSTROM2
        if abs(model - sigma) > 0.05 + 0.05 * 1e-11 / speed / ANGSTROM2:
            msg = f"{level}: sigma {sigma} against k_Q / <v> = {model:.2f}"
            raise ValueError(msg)
    energy = nist_levels()
    value = {
        (name.split("__")[1].removesuffix(".csv"), p.row): float(p.value)
        for name, points in files.items()
        for p in points
    }
    tests = {
        "k43": (value["k_2p4", "2p3"] / balance("2p3", "2p4", energy), "k_2p3_db"),
        "k98": (value["k_2p9", "2p8"] / balance("2p8", "2p9", energy), "k_2p8_db"),
    }
    for name, (model, dataset) in tests.items():
        printed = next(v for (d, _), v in value.items() if d == dataset)
        if abs(model - printed) > 0.05:
            msg = f"{name}: detailed balance {model:.3f} against {printed}"
            raise ValueError(msg)


def _check_table_iv(files: dict[str, list[Point]]) -> None:
    """Check Table IV's wavelengths, weights, sum A tau = 1 and the f column."""
    energy = nist_levels()
    tau = {p.row: float(p.value) for p in files[f"Ar/{SOURCE}__tau_Ar.csv"]}
    total: dict[str, float] = {}
    misprinted = []
    g_upper = ""
    for cells in json_table(DIRECTORY, TABLE_IV)[2:]:
        upper, lower = _line(cells[0])
        g_upper = cells[2] or g_upper
        weights = (2 * int(PASCHEN[upper][3]) + 1, 2 * int(PASCHEN[lower][3]) + 1)
        if (int(g_upper), int(cells[3])) != weights:
            msg = f"{cells[0]}: g {g_upper}, {cells[3]} against {weights}"
            raise ValueError(msg)
        wavenumber = (energy[upper] - energy[lower]) * 100  # m^-1
        air = 1e10 / wavenumber / (1 + float(standard_air_refractivity(wavenumber)))
        if abs(air - float(cells[1])) > 1.0:
            msg = f"{cells[0]}: {cells[1]} Angstrom against NIST {air:.2f}"
            raise ValueError(msg)
        if cells[4].startswith("("):
            continue
        a = float(cells[4])
        total[upper] = total.get(upper, 0.0) + a
        f_model = 1.4992e-16 * float(cells[1]) ** 2 * weights[0] / weights[1] * a * 1e6
        f_printed = cells[7]
        digits = len(f_printed.split(".")[1]) if "." in f_printed else 0
        if abs(f_model - float(f_printed)) > 0.5 * 10**-digits + 0.006 * f_model:
            misprinted.append(cells[0])
    if misprinted != ["2p_{10}-1s_2", "2p_{10}-1s_5"]:
        msg = f"f column disagrees with A for {misprinted}"
        raise ValueError(msg)
    for level, a in total.items():
        if abs(a * tau[level] * 1e-3 - 1) > 0.01:
            msg = f"{level}: sum A tau = {a * tau[level] * 1e-3:.4f}"
            raise ValueError(msg)


def checks() -> None:
    """Check the counts, labels and the tables' internal relations (module notes)."""
    files = extract()
    counts = {
        name.split("__")[1].removesuffix(".csv"): len(points)
        for name, points in files.items()
    }
    expected = {
        "tau_Ar": 8,
        "k_2p3_db": 1,
        "k_2p4": 1,
        "k_2p6": 1,
        "k_2p7": 1,
        "k_2p8": 2,
        "k_2p8_db": 1,
        "k_2p9": 2,
        "k_4s": 5,
        "kq_self": 7,
        "kq_self_bound": 1,
        "A_4p_4s": 26,
    }
    if counts != expected:
        msg = f"counts {counts}"
        raise ValueError(msg)
    _check_labels()
    _check_table_iii(files)
    _check_table_iv(files)


def main() -> None:
    """Check the tables and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
