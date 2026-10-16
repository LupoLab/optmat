# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kumar & Meath (1994): CO isotropic and anisotropic dipole properties.

The tables print powers of ten in parentheses (``1.308(1)`` = 13.08); numbers are
written in e-notation with the digits unchanged, the cell kept in ``printed``. All the
tables are unlabelled in the extraction; the paper's numbers are used here. Images
checked against the cells (2026-09-27).

- Table 1 (``tbl_a5187b988bd4``, cells file): S_k, k = 2 ... -8: the isotropic results
  of this work and of Ref. [22] (Jhanwar and Meath 1982) are extracted, rows by k
  (``-0.5`` for -1/2); the perpendicular and parallel components have no quantity
  (record notes). The mark ``^b`` (a constraint) sits before the power of ten
  (``1.380^b(4)``) and is stripped.
- Table 3 (``tbl_77c53d70381b``, cells file): the dynamic polarizability alpha and its
  anisotropy gamma against wavelength (A; the cell ``∞`` is the static limit), this
  work, Ref. [49], Ref. [50] and "experiment" (note b: Ref. [41] at infinite
  wavelength, Ref. [16] elsewhere, split into two files here).
- Table 4 (``tbl_d2eb2911fa57``, CSV): the total depolarization ratio (printed as 100
  rho) and the Rayleigh cross section (10^-27 cm^2) against wavelength.
- Table 9 (``tbl_7df843159942``): the perpendicular and parallel pseudo-DOSDs are
  formulas in the record; this script asserts that the record's terms equal the cells.
"""

import re
import tomllib
from pathlib import Path

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Kumar_1994_ReliableIsotropicAnisotropic"
TABLE_1 = "tbl_a5187b988bd4"
TABLE_3 = "tbl_77c53d70381b"
TABLE_4 = "tbl_d2eb2911fa57"
TABLE_9 = "tbl_7df843159942"
RECORD = Path(__file__).resolve().parents[2] / (
    "src/optmat/data/sources/kumar_1994.toml"
)
PREFIX = "CO/kumar_1994__"
_NUMBER = re.compile(r"(?P<mantissa>-?\d+\.\d+)(?:\((?P<exponent>-?\d+)\))?")
TABLE_3_COLUMNS = {
    1: "alpha_present",
    2: "alpha_oddershede_1982",
    3: "alpha_jorgensen_1978",
    4: "alpha_experiment",
    5: "gamma_present",
    6: "gamma_oddershede_1982",
    7: "gamma_jorgensen_1978",
}


def _value(cell: str) -> str:
    """Return a cell's number in e-notation, the digits unchanged."""
    found = _NUMBER.fullmatch(cell.strip().replace("^b", "").replace(" ", ""))
    if found is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    mantissa, exponent = found["mantissa"], found["exponent"]
    return number(f"{mantissa}e{exponent}" if exponent else mantissa)


def _order(text: str) -> str:
    text = text.strip()
    return str(int(text.removesuffix("/2")) / 2) if text.endswith("/2") else text


def _table_1() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_1)[2:]:
        order, printed = _order(cells[0]), " | ".join(cells)
        for column, stem in ((3, "moments"), (4, "moments_jhanwar_1982")):
            point = Point(order, order, _value(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}{stem}.csv", []).append(point)
    return files


def _table_3() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_3)[2:]:
        wavelength = "inf" if cells[0] == "\N{INFINITY}" else number(cells[0])
        printed = " | ".join(cells)
        for column, column_stem in TABLE_3_COLUMNS.items():
            if not cells[column]:
                continue
            stem = column_stem
            if column_stem == "alpha_experiment":
                # Note b: Watson and Ramaswamy at infinite wavelength, Smith et al.
                # at the others.
                stem = (
                    "alpha_watson_1936" if wavelength == "inf" else "alpha_smith_1976"
                )
            point = Point(wavelength, wavelength, _value(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}{stem}.csv", []).append(point)
    return files


def _table_4() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_4)[1:]:
        wavelength, printed = number(cells[0]), " | ".join(cells)
        for column, stem in ((1, "depolarization"), (2, "rayleigh")):
            point = Point(wavelength, wavelength, _value(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}{stem}.csv", []).append(point)
    return files


def _check_pseudospectra() -> None:
    """Assert that the record's Table 9 terms are the table's cells."""
    record = tomllib.loads(RECORD.read_text(encoding="utf-8"))
    formulas = {f["id"]: f for f in record["formulas"]}
    rows = table(DIRECTORY, TABLE_9)[1:]
    for name, (e_column, f_column) in {
        "perpendicular": (0, 1),
        "parallel": (2, 3),
    }.items():
        terms = formulas[f"pseudospectrum_{name}"]["coefficients"]
        for i, cells in enumerate(rows, start=1):
            assert terms[f"E{i}"] == _value(cells[e_column]), (name, i)
            assert terms[f"f{i}"] == _value(cells[f_column]), (name, i)


def extract() -> dict[str, list[Point]]:
    """Return the sums (Table 1), alpha and gamma (Table 3), rho and Q (Table 4)."""
    _check_pseudospectra()
    return _table_1() | _table_3() | _table_4()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
