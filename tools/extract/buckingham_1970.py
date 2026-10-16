# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Buckingham, Bogaard, Dunmur, Hobbs & Orr (1970): dc Kerr effect of N2, CO2.

All tables checked against the page images (pp. 1550-1552).

* Table 1 is extracted as eight fragments without headers (the blocks of six gases,
  with the polarizability lines between them as separate one-cell tables). Only the
  first two blocks are taken: CARBON DIOXIDE (``tbl_2089ba2fdb4d``, JSON cells only,
  which keeps the column header and the gas label) and NITROGEN
  (``tbl_c5d4dfda3344``, JSON cells only, headerless: the gas is identified by its
  label cell "NITROGEN" and by its position after the CO2 block on p. 1550 of the
  image). Per temperature: 10^27 A_K (``table1_A_K_<gas>``) and 10^32 B_K +- its
  standard deviation (``table1_B_K_<gas>``).
* Table 2 (``tbl_53d38067698d``, CSV): the gamma column, 10^61 gamma +- error, for
  carbon dioxide and nitrogen (``table2_gamma_<gas>``); the intercept and slope are
  the record's formulas; the other gases are out of scope.
* Table 3 (``tbl_a0c4c55e9fd4``, CSV): the optical anisotropy (a_par - a_perp) at 632.8
  nm (``table3_anisotropy_optical_<gas>``) and the static (alpha_par - alpha_perp)
  derived with the infrared atomic contributions (``table3_anisotropy_static_<gas>``).
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Buckingham_1970_KerrEffectSome"
CO2_BLOCK = "tbl_2089ba2fdb4d"
N2_BLOCK = "tbl_c5d4dfda3344"
TABLE_2 = "tbl_53d38067698d"
TABLE_3 = "tbl_a0c4c55e9fd4"
PLUS_MINUS = "\N{PLUS-MINUS SIGN}"
GASES = {"carbon dioxide": "CO2", "nitrogen": "N2"}
WAVELENGTH = "632.8"  # nm, the He-Ne line (p. 1549)


def _split(cell: str, separator: str) -> tuple[str, str]:
    value, error = cell.split(separator)
    return number(value), number(error)


def _block(table_id: str, label: str, gas: str) -> dict[str, list[Point]]:
    """Return the A_K and B_K files of one gas's block of Table 1."""
    grid = json_table(DIRECTORY, table_id)
    start = next(i for i, row in enumerate(grid) if label in row) + 1
    files: dict[str, list[Point]] = {}
    for cells in grid[start:]:
        temperature = number(cells[0])
        printed = " | ".join(cells)
        files.setdefault(f"{gas}/buckingham_1970__table1_A_K_{gas}.csv", []).append(
            Point(temperature, temperature, number(cells[1]), "", printed)
        )
        value, error = _split(cells[2], PLUS_MINUS)
        files.setdefault(f"{gas}/buckingham_1970__table1_B_K_{gas}.csv", []).append(
            Point(temperature, temperature, value, error, printed)
        )
    return files


def _tables_2_3() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = table(DIRECTORY, TABLE_2)
    if not grid[0][3].startswith("10^{61}\\gamma"):
        msg = f"unexpected Table 2 header {grid[0]!r}"
        raise ValueError(msg)
    for cells in grid[1:]:
        gas = GASES.get(cells[0])
        if gas is None:
            continue  # ethane, ethylene, acetylene, cyclopropane: out of scope
        value, error = _split(cells[3], "\\pm")
        files[f"{gas}/buckingham_1970__table2_gamma_{gas}.csv"] = [
            Point(WAVELENGTH, WAVELENGTH, value, error, " | ".join(cells))
        ]
    grid = table(DIRECTORY, TABLE_3)
    if not (
        grid[0][3].startswith("10^{40}(a_{")
        and grid[0][4].startswith("10^{40}(\\alpha")
    ):
        msg = f"unexpected Table 3 header {grid[0]!r}"
        raise ValueError(msg)
    for cells in grid[1:]:
        gas = GASES.get(cells[0])
        if gas is None:
            continue
        printed = " | ".join(cells)
        stem = f"{gas}/buckingham_1970__table3_anisotropy"
        files[f"{stem}_optical_{gas}.csv"] = [
            Point(WAVELENGTH, WAVELENGTH, number(cells[3]), "", printed)
        ]
        files[f"{stem}_static_{gas}.csv"] = [
            Point("inf", "inf", number(cells[4]), "", printed)
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the CO2 and N2 files."""
    return (
        _block(CO2_BLOCK, "CARBON DIOXIDE", "CO2")
        | _block(N2_BLOCK, "NITROGEN", "N2")
        | _tables_2_3()
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
