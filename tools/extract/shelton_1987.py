# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton (1987), Mol. Phys. 60, 65: the vibrational gamma^v of diatomics.

All tables checked against their images. Decimal points are printed as raised dots
("59·32"), written as "." in the numbers; the printed cells are kept.

* Table 2 (``tbl_6c3a7b517193``, CSV): per molecule the vibrational frequency nu_01
  (``table2_nu01_<gas>``, a Raman shift) and B_0/2 pi c (``table2_B0_<gas>``, a
  rotational constant), both in cm^-1. The nuclear-spin weights and the transition
  polarizabilities are in the record's notes (no quantity holds a Raman transition
  polarizability). The row marks (dagger, double dagger, section sign) are dropped; the
  extraction reads D2's double dagger as a dagger ("D_{2\dagger}"): the image shows the
  double dagger (J = 2 state).
* Table 3 (``tbl_d9a1dad537c6``, CSV; its caption, "Table 3. Calculated values of the
  static gamma^v at 293 K" and the unit 10^-62 C^4 m^4 J^-3, is detached from it in the
  extraction): the static gamma^v_xxxx from Eq. (17) (``table3_eq17_<gas>``) and Eq.
  (18) (``table3_eq18_<gas>``). The note marks "||" and the pilcrow are dropped. The
  ratio columns are in the record's notes.
* Table 4 (``tbl_a48923e63b94``, JSON cells only): gamma^v at nu_1 = 15800 cm^-1 for
  the dc Kerr effect, DFWM, CARS, ESHG and THG (``table4_<gas>_<process>``), the
  asterisks of the dc Kerr column dropped. Printed values; Shelton's corrections (Phys.
  Rev. A 36, 3461, ref. 10) are errata in the record. The nu_01 column (D2 printed 2998
  for Table 2's 2988) is not taken.
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Shelton_1987_VibrationalContributionsHyperpolarizabilities"
TABLE_2 = "tbl_6c3a7b517193"
TABLE_3 = "tbl_d9a1dad537c6"
TABLE_4 = "tbl_a48923e63b94"
DOT = "\N{MIDDLE DOT}"
GASES = {"H": "H2", "D": "D2", "N": "N2", "O": "O2"}
NU_1 = "15800"  # cm^-1, Table 4's fundamental (the caption's "lambda = 633 nm")
TABLE_4_PROCESSES = ("dc_kerr", "dfwm", "cars", "eshg", "thg")


def _gas(cell: str) -> str:
    r"""Return the molecule of a row cell ("H_2\dagger", "D_{2\dagger}", "N_{2}")."""
    return GASES[cell.strip()[0]]


def _value(cell: str) -> str:
    """Return the number of a cell, raised dots as points and note marks dropped."""
    return number(cell.replace(DOT, ".").rstrip("|¶*").strip())


def _table2() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        gas = _gas(cells[0])
        printed = " | ".join(cells)
        prefix = f"{gas}/shelton_1987__table2_"
        files[f"{prefix}nu01_{gas}.csv"] = [
            Point("v=0-1", "v=0-1", _value(cells[1]), "", printed)
        ]
        files[f"{prefix}B0_{gas}.csv"] = [
            Point("B_0", "B_0", _value(cells[2]), "", printed)
        ]
    return files


def _table3() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = table(DIRECTORY, TABLE_3)
    if grid[0][1:3] != ["Equation (17)", "Equation (18)"]:
        msg = f"unexpected Table 3 header {grid[0]!r}"
        raise ValueError(msg)
    for cells in grid[1:]:
        gas = _gas(cells[0])
        printed = " | ".join(cells)
        for column, name in ((1, "eq17"), (2, "eq18")):
            files[f"{gas}/shelton_1987__table3_{name}_{gas}.csv"] = [
                Point("static", "inf", _value(cells[column]), "", printed)
            ]
    return files


def _table4() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = json_table(DIRECTORY, TABLE_4)
    if [cell.rstrip("†‡§|¶") for cell in grid[1][2:]] != [
        "dc Kerr",
        "DFWM",
        "CARS",
        "ESHG",
        "THG",
    ]:
        msg = f"unexpected Table 4 header {grid[1]!r}"
        raise ValueError(msg)
    for cells in grid[2:]:
        gas = _gas(cells[0])
        printed = " | ".join(cells)
        for process, cell in zip(TABLE_4_PROCESSES, cells[2:], strict=True):
            files[f"{gas}/shelton_1987__table4_{gas}_{process}.csv"] = [
                Point(process, NU_1, _value(cell), "", printed)
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table 2, 3 and 4 files."""
    return _table2() | _table3() | _table4()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
