# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bishop, Pipin & Rérat (1990): the averaged gamma of H2 and D2.

No table has a CSV rendering except Table X; all are read from their JSON cells. Every
image was checked against the cells (no differences in the extracted columns; decimal
points and minus signs agree; thin spaces group the digits of the wavelengths).

* Tables II-VII (H2 dc Kerr, ESHG, THG; D2 the same): the averaged (295 K) electronic
  gamma_par and gamma_perp, the vibrational contributions and their totals (THG:
  gamma_par only). Each table has two blocks with different abscissae: photon energies
  hbar omega = 0.00-0.05 a.u. (hartree) and six laser wavelengths in angstrom. Both are
  kept as printed, so each column gives two files:
  ``<process>_<component>_<part>_<gas>`` (omega) and ``..._lasers_<gas>`` (lambda; the
  thin space in "10 640" dropped).
* Table IX (``tbl_cf07a9edd99b``): the ratios S = gamma_par^ESHG(total, H2)/(total, D2)
  and T = gamma_par^ESHG(total, H2)/gamma_par^ESHG(He), theory and experiment
  (Refs. 10), against lambda in angstrom; "(4)" on the last digits is written as the
  sigma; "—" (no theoretical T at 5900 A) is skipped.
* Table X (``tbl_ba8c25793416``, CSV): gamma^K = (3/2)(gamma_par^Kerr(total) -
  gamma_perp^Kerr(total)) of H2 and D2 against lambda.

Table I (molecular-frame components at R = 1.4 a0), Table VIII (the ratio R =
gamma_par/gamma_perp), Table XI (dispersion coefficients: formulas in the record) and
Tables XII-XIV (field-gradient polarizabilities, comparisons) are not extracted.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Bishop_1990_NonlinearOpticalProperties"
PROCESS_TABLES = {  # table id -> (gas, process, has gamma_perp)
    "tbl_3542c447427e": ("H2", "dc_kerr", True),
    "tbl_998201c95e2d": ("H2", "eshg", True),
    "tbl_147e4bcd04b8": ("H2", "thg", False),
    "tbl_50515faab394": ("D2", "dc_kerr", True),
    "tbl_76a1d5aaa9d0": ("D2", "eshg", True),
    "tbl_13d5a3420da6": ("D2", "thg", False),
}
PARTS = ("elec", "vib", "total")
TABLE_IX = "tbl_cf07a9edd99b"
TABLE_X = "tbl_ba8c25793416"
_UNCERTAIN = re.compile(r"^([\d.]+)\((\d+)\)$")


def _processes() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for table_id, (gas, process, perp) in PROCESS_TABLES.items():
        grid = json_table(DIRECTORY, table_id)
        if grid[0][2] != "\\bar{\\gamma}_{\\parallel} (elec)":
            msg = f"unexpected header in {table_id}: {grid[0]!r}"
            raise ValueError(msg)
        columns = [("par", part) for part in PARTS]
        if perp:
            columns += [("perp", part) for part in PARTS]
        block = ""
        for cells in grid[1:]:
            if cells[0]:
                block = {"\\omega (a.u.)": "", "\\lambda (Å)": "_lasers"}[cells[0]]
            x = number(cells[1].replace(" ", ""))
            printed = " | ".join(cells)
            for (component, part), value in zip(columns, cells[2:], strict=True):
                name = f"{process}_{component}_{part}{block}_{gas}"
                files.setdefault(f"{gas}/bishop_1990b__{name}.csv", []).append(
                    Point(x, x, number(value), "", printed)
                )
    return files


def _parenthesised(cell: str) -> tuple[str, str]:
    """Split "1.023(4)" into ("1.023", "0.004"): the sigma on the last digits."""
    match = _UNCERTAIN.match(cell)
    if match is None:
        return number(cell), ""
    value, digits = match[1], match[2]
    decimals = len(value.partition(".")[2])
    return number(value), f"{int(digits) / 10**decimals:.{decimals}f}"


def _ratios() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_IX)
    if grid[1] != ["", "Theorya", "Expt.b", "Thoeryc", "Expt.b"]:
        msg = f"unexpected Table IX header {grid[1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        lam = number(cells[0].replace(" ", ""))
        printed = " | ".join(cells)
        for name, cell in (
            ("s_theory", cells[1]),
            ("s_expt", cells[2]),
            ("t_theory", cells[3]),
            ("t_expt", cells[4]),
        ):
            if cell == "\N{EM DASH}":
                continue
            files.setdefault(f"H2/bishop_1990b__ratio_{name}_H2.csv", []).append(
                Point(lam, lam, *_parenthesised(cell), printed)
            )
    return files


def _kerr() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_X)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        lam = number(cells[0].replace(" ", ""))
        for gas, value in (("H2", cells[1]), ("D2", cells[2])):
            files.setdefault(f"{gas}/bishop_1990b__gamma_k_total_{gas}.csv", []).append(
                Point(lam, lam, number(value), "", " | ".join(cells))
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables II-VII, IX and X."""
    return _processes() | _ratios() | _kerr()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
