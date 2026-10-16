# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Buckingham & Orr (1968), Proc. R. Soc. A 305, 259: the dc Kerr effect of H2.

Both tables checked against their images (no differences). Decimal points are printed
as raised dots ("5·54"), written as ".".

* Table 1 (``tbl_4157737602db``, CSV): 17 sets of measurements, each at one
  temperature (245.4-332.4 K, most temperatures twice), with 10^14 A_K and 10^12 B_K in
  e.s.u. Rows are the sets in printed order, "set 1" ... "set 17"; the temperature is
  the abscissa. Files ``table1_A_K_H2`` and ``table1_B_K_H2``.
* Table 2 (``tbl_a2e5fa68c008``, CSV): 10^36 gamma (e.s.u.) per set, printed as two
  side-by-side lists (sets 1-9 on the left, 10-17 on the right, in Table 1's order;
  the temperatures agree set by set, which the script checks). File
  ``table2_gamma_H2``, rows "set n" with the temperature as an extra column.
"""

from _library import Point, number, table, write

DIRECTORY = "Buckingham_1968_ElectricBirefringenceMolecular"
TABLE_1 = "tbl_4157737602db"
TABLE_2 = "tbl_a2e5fa68c008"
DOT = "\N{MIDDLE DOT}"
WAVELENGTH = "6328"  # angstrom, the He-Ne line as printed (in air)


def _number(cell: str) -> str:
    return number(cell.replace(DOT, "."))


def extract() -> dict[str, list[Point]]:
    """Return the Table 1 and Table 2 files."""
    files: dict[str, list[Point]] = {}
    temperatures = []
    for index, cells in enumerate(table(DIRECTORY, TABLE_1)[1:], start=1):
        temperature = _number(cells[0])
        temperatures.append(temperature)
        printed = " | ".join(cells)
        for column, name in ((1, "A_K"), (2, "B_K")):
            files.setdefault(f"H2/buckingham_1968b__table1_{name}_H2.csv", []).append(
                Point(f"set {index}", temperature, _number(cells[column]), "", printed)
            )
    rows = table(DIRECTORY, TABLE_2)[1:]
    pairs = [(cells[0], cells[1], cells) for cells in rows]
    pairs += [(cells[2], cells[3], cells) for cells in rows if cells[2]]
    if len(pairs) != len(temperatures):
        msg = "Table 2 should give one gamma per set of Table 1"
        raise ValueError(msg)
    for index, (temperature, gamma, cells) in enumerate(pairs, start=1):
        if _number(temperature) != temperatures[index - 1]:
            msg = f"Table 2 set {index} is at {temperature}, not Table 1's"
            raise ValueError(msg)
        files.setdefault("H2/buckingham_1968b__table2_gamma_H2.csv", []).append(
            Point(
                f"set {index}",
                WAVELENGTH,
                _number(gamma),
                "",
                " | ".join(cells),
                {"temperature": _number(temperature)},
            )
        )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
