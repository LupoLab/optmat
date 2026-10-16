# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bruce et al. (1990): lifetimes and quenching of two-photon excited Xe.

J. Chem. Phys. 92, 2917-2926 (1990); journal p. 2917 is PDF page 2 (PDF page 1 is the
AIP download cover), so the journal page is the PDF page plus 2915. The tables are read
from the library's cells (``.json``; merged headers), every cell checked against the
table images (2026-10-09); the row is the level in the paper's (modified Racah)
notation, ``6p[5/2]2``, and ``x`` the same label. Errors are 'two standard deviations'
throughout (p. 2919). Taken, the column 'This work' only:

* **Table III** (PDF p. 7, journal p. 2922, ``tbl_3cb99347e778``), radiative lifetimes
  (ns), the zero-pressure intercepts in pure xenon: the eight 6p, 7p and 6p' levels,
  and 5d[1/2]1 (radiation trapped; a separate file).
* **Table II** (PDF p. 4, journal p. 2919, ``tbl_3514df8be433``), total quenching by
  xenon: the bimolecular rates (10^-12 cm^3/s) of nine levels; the termolecular rate
  (10^-31 cm^6/s) of 6p[1/2]0; and the bound '<4' of 6p[3/2]2, held as the interval
  0-4, value 2 with half-width 2 (``printed`` keeps '<4').
* **Table IV** (PDF p. 8, journal p. 2923, ``tbl_b2aa0277aa53``), total quenching by
  argon (10^-12 cm^3/s): the seven levels with a single decay, and, in a separate file,
  the three rates of the coupled pair 6p[1/2]0 / 5d[1/2]1 (footnotes b, c, d: the
  slopes of the fast and slow decay components), which are not the levels' removal
  coefficients (the appendix's are; :func:`checks`).

The other columns (other papers) are not taken. The appendix's state-to-state rates
are inline points of the record.

**Checks** (:func:`checks`; the script fails if one is missed): the table counts; the
appendix's scheme (record inline datasets) against the two-state model (Eqs. A9-A12)
with chi = 0.1752: k21 = chi k12 (29.8 against 30 +- 1), k1 = k12 + k10 (183), k2 =
k21 + k20 (34.6 against 34 +- 2), Table IV's fast-component slope k12 (1 + chi) (199.8
against 200 +- 4.0) and slow-component slope (k10 chi + k20) / (1 + chi) (5.85 against
5.9 +- 0.4); and chi = (1/3) exp(-132.3 cm^-1 / kT) = 0.1752 at 296 K (0.1767 at the
printed 300 K).
"""

import math
import re
import tomllib

from _json_table import json_table
from _library import DATA, Point, number, write

from optmat.constants import HBAR, K_B, C

DIRECTORY = "Bruce_1990_RadiativeLifetimesCollisional"
SOURCE = "bruce_1990"
TABLE_II = "tbl_3514df8be433"
TABLE_III = "tbl_3cb99347e778"
TABLE_IV = "tbl_b2aa0277aa53"
_LEVEL = re.compile(r"^(\d[a-z]'?\[\d/2\])_\{?(\d)\}?$")
_VALUE = re.compile(r"^([\d.]+) (?:\\pm|±) ([\d.]+)(?:\^\{?([a-z])\}?)?$")
CHI = 0.1752
# Table IV's rows of the coupled pair: footnote -> row id.
COUPLED = {"b": "6p[1/2]0", "c": "5d[1/2]1_slow", "d": "5d[1/2]1_vuv"}


def _level(cell: str) -> str:
    """Return a level label as the paper prints it, without the LaTeX subscript."""
    match = _LEVEL.match(cell.strip())
    if match is None:
        msg = f"unexpected level {cell!r}"
        raise ValueError(msg)
    return "".join(match.groups())


def _value(cell: str) -> tuple[str, str, str]:
    """Return the value, the 2-sigma error and the footnote mark of a cell."""
    match = _VALUE.match(cell.strip())
    if match is None:
        msg = f"unexpected value {cell!r}"
        raise ValueError(msg)
    return number(match.group(1)), match.group(2), match.group(3) or ""


def _lifetimes() -> dict[str, list[Point]]:
    """Return Table III's 'This work' lifetimes: the 6p, 7p, 6p' levels and 5d[1/2]1."""
    grid = json_table(DIRECTORY, TABLE_III)
    row = next(cells for cells in grid if cells[0] == "This work")
    printed = " | ".join(row)
    main, trapped = [], []
    for header, cell in zip(grid[0][1:], row[1:], strict=True):
        level = _level(header)
        value, sigma, _ = _value(cell)
        point = Point(level, level, value, sigma, f"{header}: {printed}")
        (trapped if level == "5d[1/2]1" else main).append(point)
    return {
        f"Xe/{SOURCE}__tau_Xe.csv": main,
        f"Xe/{SOURCE}__tau_5d_Xe.csv": trapped,
    }


def _xenon() -> dict[str, list[Point]]:
    """Return Table II's 'This work': bimolecular, termolecular, and the bound."""
    grid = json_table(DIRECTORY, TABLE_II)
    section = ""
    two, three, bound = [], [], []
    for cells in grid[1:]:
        if not cells[1]:
            section = cells[0]
            continue
        level = _level(cells[0])
        printed = f"{section}: " + " | ".join(cells)
        if cells[1] == "<4":
            bound.append(Point(level, level, "2", "2", printed))
            continue
        value, sigma, _ = _value(cells[1])
        target = two if section.startswith("Bimolecular") else three
        target.append(Point(level, level, value, sigma, printed))
    return {
        f"Xe/{SOURCE}__kq_self_Xe.csv": two,
        f"Xe/{SOURCE}__kq3_self_Xe.csv": three,
        f"Xe/{SOURCE}__kq3_bound_Xe.csv": bound,
    }


def _argon() -> dict[str, list[Point]]:
    """Return Table IV's 'This work': single levels, and the coupled pair."""
    grid = json_table(DIRECTORY, TABLE_IV)
    single, coupled = [], []
    level = ""
    for cells in grid[2:]:
        level = _level(cells[0]) if cells[0] else level
        value, sigma, mark = _value(cells[1])
        printed = " | ".join(cells)
        if mark:
            coupled.append(Point(COUPLED[mark], level, value, sigma, printed))
        else:
            single.append(Point(level, level, value, sigma, printed))
    return {
        f"Xe/{SOURCE}__kq_Ar_Xe.csv": single,
        f"Xe/{SOURCE}__kq_Ar_coupled_Xe.csv": coupled,
    }


def extract() -> dict[str, list[Point]]:
    """Return every file of the record."""
    return _lifetimes() | _xenon() | _argon()


def _inline() -> dict[str, float]:
    """Return the record's appendix rates (inline points) by dataset id."""
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        record = tomllib.load(handle)
    return {
        entry["id"]: float(entry["points"][0][2])
        for entry in record["datasets"]
        if entry.get("table") == "text"
    }


def checks() -> None:
    """Check the counts and the appendix against the two-state model (module notes)."""
    files = extract()
    counts = [len(points) for points in files.values()]
    if counts != [8, 1, 9, 1, 1, 7, 3]:
        msg = f"counts {counts}"
        raise ValueError(msg)
    k = _inline()
    table = {p.row: float(p.value) for p in files[f"Xe/{SOURCE}__kq_Ar_coupled_Xe.csv"]}
    tests = {
        "k21 = chi k12": (CHI * k["k12_Ar"], k["k21_Ar"], 1.0),
        "k1 = k12 + k10": (k["k12_Ar"] + k["k10_Ar"], k["k1_Ar"], 0.5),
        "k2 = k21 + k20": (k["k21_Ar"] + k["k20_Ar"], k["k2_Ar"], 1.0),
        "fast": (k["k12_Ar"] * (1 + CHI), table["6p[1/2]0"], 0.5),
        "slow": (
            (k["k10_Ar"] * CHI + k["k20_Ar"]) / (1 + CHI),
            table["5d[1/2]1_slow"],
            0.1,
        ),
    }
    for name, (model, printed, tolerance) in tests.items():
        if abs(model - printed) > tolerance:
            msg = f"{name}: {model:.3f} against {printed}"
            raise ValueError(msg)
    # 1 cm^-1 is 100 h c J; chi at the printed 300 K and at 296 K.
    hc = 2 * math.pi * HBAR * C * 100
    chi = {t: math.exp(-132.3 * hc / (K_B * t)) / 3 for t in (296.0, 300.0)}
    if round(chi[296.0], 4) != CHI or round(chi[300.0], 4) != 0.1767:
        msg = f"chi {chi}"
        raise ValueError(msg)


def main() -> None:
    """Check the tables and the appendix, and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
