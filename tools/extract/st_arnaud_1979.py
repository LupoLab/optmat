# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract St-Arnaud & Bose (1979): refractivity and dielectric virials of CO2 and SF6.

Three tables (cells files only; images checked, cells agree). In each, the gas is
named on the first row of its block, and the last column names the reference ("This
work" with a footnote mark, "Theory", or a reference number); "..." is an empty cell;
values print as "value +- sigma" (sigma three standard deviations for this work) or
with a footnote letter, stripped here.

- Table I (``tbl_42f95d223f73``): A_R (6328 A) of this work (323 K) and of Buckingham
  and Graham (Ref. 11, 298 K); the dielectric A_eps of Refs. 4 (CO2) and 5 (SF6). The
  rows without a reference (Landolt-Bornstein A_R, footnote a) and Refs. 7, 8 (no
  wavelength) are not extracted, nor the derived alpha columns.
- Table II (``tbl_4533c14b419a``): B_R of this work (323 K), the theory (classical
  dipole-induced dipole, 323 K) and Buckingham and Graham (298 K). B_P (a density
  virial) and the B_R rows of Refs. 7, 8 and 9 (no single wavelength) are not
  extracted.
- Table III (``tbl_4e565abf97f8``): the dielectric B_eps of Refs. 5 (SF6) and 4 (CO2)
  at 323 K; the B_R column repeats Table II; B_FIR and B_IR have no quantity.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "StArnaud_1979_DirectDeterminationIntermolecular"
TABLE_I = "tbl_42f95d223f73"
TABLE_II = "tbl_4533c14b419a"
TABLE_III = "tbl_4e565abf97f8"
GASES = {"SF_6": "SF6", "CO_2": "CO2"}
_VALUE = re.compile(r"^(?P<value>-?[\d.]+)\s*(?:±\s*(?P<sigma>[\d.]+))?[a-e]?$")


def _split(cell: str) -> tuple[str, str] | None:
    """Return (value, sigma) of a cell, or None for "..." or an empty cell."""
    text = cell.strip()
    if text in {"", "..."}:
        return None
    match = _VALUE.match(text)
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    return number(match["value"]), match["sigma"] or ""


def _rows(table_id: str) -> list[tuple[str, list[str]]]:
    """Return (gas, cells) for each data row, the gas carried down its block."""
    out, gas = [], ""
    for cells in json_table(DIRECTORY, table_id)[1:]:
        gas = GASES.get(cells[0].strip(), gas)
        out.append((gas, cells))
    return out


def _add(
    files: dict[str, list[Point]],
    name: str,
    gas: str,
    wavelength: str,
    cell: str,
    printed: str,
) -> None:
    parsed = _split(cell)
    assert parsed is not None, (name, gas)
    value, sigma = parsed
    files[f"{gas}/st_arnaud_1979__{name}_{gas}.csv"] = [
        Point(wavelength, wavelength, value, sigma, printed)
    ]


def extract() -> dict[str, list[Point]]:
    """Return the A_R, A_eps, B_R and B_eps datasets."""
    files: dict[str, list[Point]] = {}
    for gas, cells in _rows(TABLE_I):
        reference, printed = cells[7].strip().rstrip("*"), " | ".join(cells)
        if reference == "This work":
            _add(files, "A_R_this_work", gas, "6328", cells[2], printed)
        elif reference == "11":
            _add(files, "A_R_buckingham_graham_1974", gas, "6328", cells[2], printed)
        elif reference in {"4", "5"}:
            _add(files, "A_eps_dielectric", gas, "inf", cells[3], printed)
    for gas, cells in _rows(TABLE_II):
        reference, printed = cells[4].strip(), " | ".join(cells)
        if reference == "This workc":
            _add(files, "B_R_this_work", gas, "6328", cells[2], printed)
        elif reference == "Theory":
            _add(files, "B_R_theory", gas, "inf", cells[2], printed)
        elif reference == "11":
            _add(files, "B_R_buckingham_graham_1974", gas, "6328", cells[2], printed)
    for gas, cells in _rows(TABLE_III):
        if _split(cells[2]) is not None:
            _add(files, "B_eps_dielectric", gas, "inf", cells[2], " | ".join(cells))
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
