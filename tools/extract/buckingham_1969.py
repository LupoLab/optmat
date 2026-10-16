# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Buckingham & Orr (1969): the dc Kerr effect of methane (Table 1).

Table 1 (``tbl_d2bf39485ac8``, JSON cells only; image checked) gives, per gas and
temperature, the pressure range, the number of measurements, 10^27 A_K and 10^31 B_K.
Only METHANE (the first block, three temperatures) is taken: A_K
(``table1_A_K_CH4``) and B_K (``table1_B_K_CH4``), the temperature as printed as the
row id and the abscissa. CF4, CH3F, CH2F2 and CHF3 are out of scope (the record's
notes). Decimal points are printed as raised dots ("0·150"), written as ".". The A_K
header is extracted as "C^2m^3J^{-2}"; the image shows C^2 m^5 J^-2 (structural: the
unit is the record's, and the script checks only the start of the header).
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Buckingham_1969_KerrEffectMethane"
TABLE_1 = "tbl_d2bf39485ac8"
DOT = "\N{MIDDLE DOT}"


def extract() -> dict[str, list[Point]]:
    """Return the methane A_K and B_K files."""
    grid = json_table(DIRECTORY, TABLE_1)
    a_k, b_k = grid[0][3], grid[0][4]
    if not (a_k.startswith("10^{27}A_K") and b_k.startswith("10^{31}B_K")):
        msg = f"unexpected Table 1 header {grid[0]!r}"
        raise ValueError(msg)
    if grid[1][0] != "METHANE":
        msg = f"unexpected first block {grid[1][0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        if not cells[1]:
            break  # the next gas's heading
        temperature = number(cells[0].replace(DOT, "."))
        printed = " | ".join(cells)
        for column, name in ((3, "A_K"), (4, "B_K")):
            value = number(cells[column].replace(DOT, "."))
            files.setdefault(f"CH4/buckingham_1969__table1_{name}_CH4.csv", []).append(
                Point(temperature, temperature, value, "", printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
