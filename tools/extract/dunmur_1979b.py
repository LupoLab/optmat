# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Dunmur, Hunt & Jessup (1979), Part II: Kerr virial coefficients.

All images checked against the cells (no differences; decimal points printed as raised
dots in Tables 3, 4 and 8, turned into points here). Only Ar, Kr, Xe, CH4 and SF6 are
taken (CF4 and neopentane, Tables 6-7, are out of scope).

* Tables 1-3 (Ar, Kr, Xe; CSV; lambda = 458 nm): per temperature the first Kerr virial
  coefficient A_K (10^29 C^2 m^5 J^-2; one temperature-independent fit coefficient,
  printed on every row: taken once, ``A_K_<gas>``, against the wavelength 458 nm, from
  the first row) and the second, B_K (10^33 C^2 m^8 J^-2; ``B_K_<gas>``, against
  temperature). Each is printed ``value ± error``. The B_K^expt/B_K^DID column is not
  extracted (no quantity).
* Table 4 (SF6, CSV; lambda = 633 nm except the 297 K row, "297†", measured at 458 nm
  at low density): A_K (the 297 K row only), B_K (``B_K_SF6`` at 633 nm;
  ``B_K_458_SF6``, the 297 K row) and C_K (10^36 C^2 m^11 J^-2, ``C_K_SF6``); "—" is
  no value.
* Table 5 (CH4, CSV; 458 nm): A_K (the 298 K low-density row only), B_K and C_K
  (10^38 C^2 m^11 J^-2).
* Table 8 (``tbl_623370bf3c99``, JSON cells): the second hyperpolarizability
  10^62 gamma-bar (C^4 m^4 J^-3), this work (``gamma_<gas>``), and the published
  values beside it (reproduced, one file per gas and process, row = the reference
  number as printed): dc Kerr [3], [5]; three-wave mixing [22], [23]; ESHG [22],
  [24], [25]; THG [26]. A value printed without an error ("20 [22]") has sigma "".
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Dunmur_1979_InfluenceIntermolecularInteractions_0bb7bb"
DOT = "\N{MIDDLE DOT}"
DAGGER = "\N{DAGGER}"
DASH = "\N{EM DASH}"
KERR_TABLES = {
    "Ar": "tbl_3848d596ace9",
    "Kr": "tbl_6fa5749b14cf",
    "Xe": "tbl_c76c4b08dfc1",
    "SF6": "tbl_08c025d2f30b",
    "CH4": "tbl_82697e3efe49",
}
TABLE_8 = "tbl_623370bf3c99"
GASES_8 = {"Ar": "Ar", "Kr": "Kr", "Xe": "Xe", "CH_4": "CH4", "SF_6": "SF6"}
# Table 8's columns after "This work": dc Kerr (value | reference), three-wave mixing,
# ESHG, THG; the reference is printed in brackets inside the last three.
_CITED = re.compile(r"^([\d.]+)(?: ± ([\d.]+))? \[(\d+)\]$")


def _clean(cell: str) -> str:
    return cell.replace(DOT, ".").strip()


def _split(cell: str) -> tuple[str, str]:
    """Return the value and the error of "6.2 ± 0.5"."""
    value, sign, sigma = _clean(cell).partition("±")
    if not sign:
        msg = f"no error in {cell!r}"
        raise ValueError(msg)
    return number(value), number(sigma)


def _kerr(gas: str) -> dict[str, list[Point]]:
    grid = table(DIRECTORY, KERR_TABLES[gas])
    if grid[0][0] != "T/K" or "A_{\\text{K}}" not in grid[0][2]:
        msg = f"unexpected header in {KERR_TABLES[gas]}: {grid[0]!r}"
        raise ValueError(msg)
    has_c = len(grid[0]) == 6
    prefix = f"{gas}/dunmur_1979b__"
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        printed = " | ".join(cells)
        low_density = cells[0].endswith(DAGGER)  # measured at 458 nm (SF6)
        temperature = number(cells[0].removesuffix(DAGGER))
        a_k = _clean(cells[2])
        if a_k != DASH and f"{prefix}A_K_{gas}.csv" not in files:
            files[f"{prefix}A_K_{gas}.csv"] = [
                Point(temperature, "458", *_split(a_k), printed)
            ]
        b_k = _clean(cells[3])
        if b_k != DASH:
            name = f"B_K_458_{gas}" if low_density else f"B_K_{gas}"
            files.setdefault(f"{prefix}{name}.csv", []).append(
                Point(temperature, temperature, *_split(b_k), printed)
            )
        if has_c and _clean(cells[5]) != DASH:
            files.setdefault(f"{prefix}C_K_{gas}.csv", []).append(
                Point(temperature, temperature, *_split(cells[5]), printed)
            )
    return files


def _gamma() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_8)
    if grid[1][1:3] != ["This work", "Published values"]:
        msg = f"unexpected Table 8 header {grid[1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in grid[2:]:
        if cells[0]:
            gas = GASES_8.get(cells[0], "")
        if not gas:
            continue
        printed = " | ".join(cells)
        prefix = f"{gas}/dunmur_1979b__"
        if cells[1]:
            files[f"{prefix}gamma_{gas}.csv"] = [
                Point("this work", "458", *_split(cells[1]), printed)
            ]
        if cells[2]:
            value, sigma = _split(cells[2])
            reference = cells[3].strip("[]")
            files.setdefault(f"{prefix}gamma_dc_kerr_lit_{gas}.csv", []).append(
                Point(reference, "not stated", value, sigma, printed)
            )
        for column, process in ((4, "twm"), (5, "eshg"), (6, "thg")):
            cell = _clean(cells[column])
            if not cell:
                continue
            match = _CITED.match(re.sub(r"\s+", " ", cell.replace("± ", "± ")))
            if match is None:
                msg = f"unexpected Table 8 cell {cells[column]!r}"
                raise ValueError(msg)
            value, sigma, reference = match[1], match[2] or "", match[3]
            files.setdefault(f"{prefix}gamma_{process}_lit_{gas}.csv", []).append(
                Point(reference, "not stated", number(value), sigma, printed)
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables 1-5 and 8."""
    files: dict[str, list[Point]] = {}
    for gas in KERR_TABLES:
        files |= _kerr(gas)
    return files | _gamma()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
