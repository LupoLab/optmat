# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Tarazkar, Romanov & Levis (2014, Phys. Rev. A): the gamma of Xe (Table II).

Table II (``tbl_6edfad81b32e``, CSV; image checked, every cell agrees): CCSD static, dc
Kerr, ESHG and DFWM gamma of Xe at 800 nm in four x-aug-cc-pVQZ-pp basis sets, and two
quoted literature values (rows 'a' and 'b', static only). One file per column; row ids
are the basis names. The fourth-order (gamma^(4)) tables are not taken (no quantity; see
the record).
"""

from _library import Point, table, write

DIRECTORY = "Tarazkar_2014_HighOrderNonlinear"
TABLE_II = "tbl_6edfad81b32e"
COLUMNS = {1: "static", 2: "dc_kerr", 3: "eshg", 4: "kerr"}  # DFWM is the kerr process
PREFIX = "Xe/tarazkar_2014__"


def extract() -> dict[str, list[Point]]:
    """Return Table II by column, and the quoted values."""
    grid = table(DIRECTORY, TABLE_II)
    if grid[0] != ["", "Static", "dc Kerr", "ESHG", "DFWM"]:
        msg = f"unexpected Table II header: {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        printed = " | ".join(cells)
        if cells[0] in {"a", "b"}:  # literature values, static column only
            value, _, sigma = cells[1].partition(" ± ")
            float(value)
            files.setdefault(f"{PREFIX}quoted_static.csv", []).append(
                Point(cells[0], "inf", value, sigma, printed)
            )
            continue
        for column, name in COLUMNS.items():
            value = cells[column]
            float(value)
            wavelength = "inf" if name == "static" else "800"
            files.setdefault(f"{PREFIX}{name}.csv", []).append(
                Point(cells[0], wavelength, value, "", printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
