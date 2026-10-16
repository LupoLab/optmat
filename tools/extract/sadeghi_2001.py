# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Sadeghi et al. (2001): quenching of two-photon excited Ar 4p by 22 gases.

J. Chem. Phys. 115, 3144-3154 (2001); journal p. 3144 is PDF page 2 (PDF page 1 is the
AIP download cover), so the journal page is the PDF page plus 3142. The tables are read
from the library's cells (``.json``; merged headers), every cell checked against the
table images (2026-10-12). Rows are the levels in the paper's Paschen notation, 2p1,
2p5, 2p6, 2p8 (4p'[1/2]0, 4p[1/2]0, 4p[3/2]2, 4p[5/2]2; the record's notes give the
NIST levels), and ``x`` the same label. Taken, one file per reagent gas (the
perturber), the k_Q columns in 10^-10 cm^3 s^-1 (the sigma_Q columns are k_Q / <v>,
derived, and not taken):

* **Table II** (PDF p. 6, journal p. 3148, ``tbl_4c7e7e30c054``): He, Ne, Ar, Kr, Xe,
  H2, D2, N2, CO. Kr has two values for 2p6 on two lines: 0.47 (footnotes c, d: a lower
  limit, from all the data of Fig. 5) and 0.75 (d, f: the fit below 15 mbar, which the
  text prefers); the first goes to its own file, the second is the row of the Kr file.
* **Table III** (PDF p. 8, journal p. 3150, ``tbl_f1b459c74236``): O2, NO, CO2, CF4,
  CHF3, CH4, C2H2, C2H4, C2H6, SF6, NO2, F2, Cl2. Its last column (Ar(4s,4s'), Ref.
  1's averages for the metastable levels) is another paper's and is not taken.

Footnote marks are dropped from the value and kept in ``printed``; the record gives the
footnotes' uncertainties (d: 30 %, Table III's f: 20 %, otherwise 10 %).

**Checks** (:func:`checks`; the script fails if one is missed): the table counts; Table
I's two-photon energies, 2 / lambda_VUV, against its printed level energies (within the
rounding of lambda_VUV to 0.1 nm) and the printed energies against NIST's levels
(``data/external/nist_asd/Ar_I_levels.tsv``; within 1 cm^-1); the radiative lifetimes
(the record's inline points) against the printed radiative decay rates (p. 3147); the
text's values of k_Ar (0.16 +- 0.02 for 2p1, 0.12, 0.13 and 0.24) and of k_Kr for 2p8
(1.6 +- 0.3) against the tables.
"""

import csv
import re
import tomllib

from _json_table import json_table
from _library import DATA, Point, write

DIRECTORY = "Sadeghi_2001_QuenchingRateConstants"
SOURCE = "sadeghi_2001"
TABLE_I = "tbl_db6e23a95906"
TABLE_II = "tbl_4c7e7e30c054"
TABLE_III = "tbl_f1b459c74236"
LEVELS = ("2p1", "2p5", "2p6", "2p8")
NIST = {
    "2p1": ("3s2.3p5.(2P*<1/2>).4p", "2[1/2]", "0"),
    "2p5": ("3s2.3p5.(2P*<3/2>).4p", "2[1/2]", "0"),
    "2p6": ("3s2.3p5.(2P*<3/2>).4p", "2[3/2]", "2"),
    "2p8": ("3s2.3p5.(2P*<3/2>).4p", "2[5/2]", "2"),
}
"""The NIST configuration, term and J of each Paschen label."""
_REAGENT = re.compile(r"^([A-Za-z_0-9]+) \(")
_VALUE = re.compile(r"^([\d.]+)(?:\^\{?([a-z,]+)\}?)?$")
RATES = {"2p1": 4.2, "2p5": 4.2, "2p6": 3.6, "2p8": 3.4}
"""The radiative decay rates printed on p. 3147, 10^7 s^-1."""


def perturber(cell: str) -> str:
    """Return a reagent cell's gas, 'C_2H_6 (11.5)^a' -> 'C2H6'; argon is 'self'."""
    match = _REAGENT.match(cell.strip())
    if match is None:
        msg = f"unexpected reagent {cell!r}"
        raise ValueError(msg)
    gas = match.group(1).replace("_", "")
    return "self" if gas == "Ar" else gas


def _value(cell: str) -> tuple[str, str]:
    """Return a k_Q cell's number and its footnote marks."""
    match = _VALUE.match(cell.strip())
    if match is None:
        msg = f"unexpected value {cell!r}"
        raise ValueError(msg)
    return match.group(1), match.group(2) or ""


def file(gas: str) -> str:
    """Return the data file of a perturber."""
    return f"Ar/{SOURCE}__kq_{gas}.csv"


def _rates(table_id: str, label: str) -> dict[str, list[Point]]:
    """Return a table's k_Q columns of 2p1, 2p5, 2p6, 2p8, one file per reagent."""
    grid = json_table(DIRECTORY, table_id)
    out: dict[str, list[Point]] = {}
    gas = ""
    for cells in grid[2:]:
        if cells[0]:
            gas = perturber(cells[0])
        printed = f"{label}" + ("" if cells[0] else f" (second line of {gas})")
        printed += ": " + " | ".join(cells[:9])
        for index, level in enumerate(LEVELS):
            cell = cells[1 + 2 * index]
            if not cell:
                continue
            value, marks = _value(cell)
            # Kr, 2p6: the second line (0.75, d, f) is the row; the first (0.47, c, d)
            # the lower limit from all the data, in its own file.
            name = "Kr_full" if gas == "Kr" and marks == "c,d" else gas
            out.setdefault(file(name), []).append(
                Point(level, level, value, "", printed)
            )
    return out


def extract() -> dict[str, list[Point]]:
    """Return every file of the record."""
    return _rates(TABLE_II, "Table II") | _rates(TABLE_III, "Table III")


def _nist() -> dict[str, float]:
    """Return NIST's energies (cm^-1) of the four levels."""
    path = DATA / "external" / "nist_asd" / "Ar_I_levels.tsv"
    found = {}
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = tuple(row[name].strip('"') for name in ("Configuration", "Term", "J"))
            for level, nist in NIST.items():
                if key == nist:
                    found[level] = float(row["Level (cm-1)"].strip('"'))
    return found


def checks() -> None:
    """Check counts, Table I, the lifetimes and the text's rates (module notes)."""
    files = extract()
    counts = {name: len(points) for name, points in files.items()}
    if len(files) != 23 or counts[file("Kr_full")] != 1 or counts[file("Kr")] != 4:
        msg = f"files {counts}"
        raise ValueError(msg)
    if any(n != 4 for name, n in counts.items() if "Kr_full" not in name):
        msg = f"counts {counts}"
        raise ValueError(msg)
    nist = _nist()
    grid = json_table(DIRECTORY, TABLE_I)
    for cells, level in zip(grid[1:], LEVELS, strict=True):
        if not cells[0].endswith(level.replace("p", "p_")):
            msg = f"Table I row {cells[0]!r} is not {level}"
            raise ValueError(msg)
        energy = float(cells[1].replace(" ", ""))
        vuv = float(cells[2])
        two_photon = 2e7 / vuv
        rounding = 2e7 / (vuv - 0.05) - two_photon
        if abs(two_photon - energy) > rounding or abs(nist[level] - energy) > 1:
            msg = f"{level}: {energy} against 2/lambda {two_photon:.0f}, NIST {nist}"
            raise ValueError(msg)
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        record = tomllib.load(handle)
    tau = next(d for d in record["datasets"] if d["quantity"] == "radiative_lifetime")
    for row, _, value, _ in tau["points"]:
        if round(100 / float(value), 1) != RATES[row]:
            msg = f"{row}: 1/{value} ns against {RATES[row]}e7 s^-1"
            raise ValueError(msg)
    table = {(name, p.row): float(p.value) for name, ps in files.items() for p in ps}
    text = {
        (file("self"), "2p1"): 0.16,
        (file("self"), "2p5"): 0.12,
        (file("self"), "2p6"): 0.13,
        (file("self"), "2p8"): 0.24,
        (file("Kr"), "2p8"): 1.6,
    }
    for key, value in text.items():
        if table[key] != value:
            msg = f"{key}: {table[key]} against the text's {value}"
            raise ValueError(msg)


def main() -> None:
    """Check the tables and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
