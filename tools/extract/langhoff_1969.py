# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Langhoff & Karplus (1969): Cauchy moments, Pade values (Tables I, III-V).

Table I (``tbl_2b7dd49ddcb9``): the four Cauchy coefficients alpha_0..alpha_3 per
species. The paper identifies alpha_k with the sum S(-2k-2) (footnote 19), so each row
is extracted as Cauchy moments S(-2), S(-4), S(-6), S(-8) in atomic units, the moment
order in the wavelength column. Atomic hydrogen (not a species of the package) is
skipped; the H_2/N_2/O_2 labels are mapped to the package's names.

Table III (merged header; JSON cells): predicted n - 1 (x 10^4, NTP) for Ne-Xe (the
paper's) and Kingston's previous calculation (reproduced); ``...`` is no value.
Table IV (argon at Heddle's wavelengths): each cell ``value+-spread``; the 1216 A row of
the measured and previous-calculation columns holds two values each (``5.65, 5.57e``:
Gill and Heddle, and Marmo marked e; ``5.30, 4.98f``: Kingston, and Liggett and
Levinger marked f), split here by the footnote letters.
Table V (merged header; JSON cells): Wilkinson's measured nitrogen values and the five
Pade approximants.
"""

import json

from _library import LIBRARY, Point, number, table, write

DIRECTORY = "Langhoff_1969_PadeSummationCauchy"
TABLE_I = "tbl_2b7dd49ddcb9"
TABLE_III = "tbl_af4050e393d8"
TABLE_IV = "tbl_c1a45b3fe5fa"
TABLE_V = "tbl_6e6971f10d5e"
SPECIES = {"He": "He", "Ne": "Ne", "Ar": "Ar", "Kr": "Kr", "Xe": "Xe"}
SPECIES |= {"H_2": "H2", "N_2": "N2", "O_2": "O2"}
ORDERS = ("-2", "-4", "-6", "-8")  # alpha_0..alpha_3 = S(-2)..S(-8)
APPROXIMANTS = ("pade11f", "pade21f", "pade22f", "pade11k", "pade10k")


def grid(table_id: str) -> list[list[str]]:
    """Return a library table's cell grid from its JSON file (merged cells once)."""
    path = LIBRARY / DIRECTORY / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))
    rows = [[""] * cells["columns"] for _ in range(cells["rows"])]
    for cell in cells["cells"]:
        rows[cell["row"]][cell["column"]] = cell["text"].strip()
    return rows


def _split(cell: str) -> tuple[str, str]:
    """Return the value and the spread of a ``3.25±0.03`` cell."""
    value, _, spread = cell.partition("±")
    return number(value), spread


def extract() -> dict[str, list[Point]]:
    """Return the moments per species and the prediction columns per gas."""
    files: dict[str, list[Point]] = {}

    def add(gas: str, name: str, point: Point) -> None:
        files.setdefault(f"{gas}/langhoff_1969__{name}.csv", []).append(point)

    for cells in table(DIRECTORY, TABLE_I)[1:]:
        if cells[0] not in SPECIES:
            continue  # atomic hydrogen
        gas = SPECIES[cells[0]]
        for order, cell in zip(ORDERS, cells[1:5], strict=True):
            add(
                gas,
                f"cauchy_{gas}",
                Point(order, order, number(cell), "", " | ".join(cells)),
            )

    columns = [("theory", "Ne"), ("theory", "Ar"), ("theory", "Kr"), ("theory", "Xe")]
    columns += [
        ("kingston_1964", "Ar"),
        ("kingston_1964", "Kr"),
        ("kingston_1964", "Xe"),
    ]
    for cells in grid(TABLE_III)[2:]:
        for (name, gas), cell in zip(columns, cells[1:], strict=True):
            if cell != "...":
                point = Point(cells[0], cells[0], number(cell), "", " | ".join(cells))
                add(gas, f"{name}_{gas}", point)

    for cells in table(DIRECTORY, TABLE_IV)[1:]:
        printed = " | ".join(cells)
        wavelength = cells[0].split("±")[0]
        if wavelength == "1216":
            gill, marmo = (c.strip() for c in cells[1].split(","))
            kingston, liggett = (c.strip() for c in cells[2].split(","))
            add("Ar", "gill_1963_Ar", Point("1216", "1216", number(gill), "", printed))
            marmo_point = Point(
                "1216", "1216", number(marmo.removesuffix("e")), "", printed
            )
            add("Ar", "marmo_Ar", marmo_point)
            point = Point("1216", "1216", number(kingston), "", printed)
            add("Ar", "kingston_1964_table4_Ar", point)
            point = Point(
                "1216", "1216", number(liggett.removesuffix("f")), "", printed
            )
            add("Ar", "liggett_1968_Ar", point)
        else:
            value, spread = _split(cells[1])
            add(
                "Ar",
                "heddle_1963_Ar",
                Point(wavelength, wavelength, value, spread, printed),
            )
            value, spread = _split(cells[2])
            point = Point(wavelength, wavelength, value, spread, printed)
            add("Ar", "kingston_1964_table4_Ar", point)
        for name, cell in zip(APPROXIMANTS, cells[3:], strict=True):
            value, spread = _split(cell)
            add(
                "Ar",
                f"{name}_Ar",
                Point(wavelength, wavelength, value, spread, printed),
            )

    for cells in grid(TABLE_V)[2:]:
        printed = " | ".join(cells)
        point = Point(cells[0], cells[0], number(cells[1]), "", printed)
        add("N2", "wilkinson_1960_N2", point)
        for name, cell in zip(APPROXIMANTS, cells[2:], strict=True):
            add(
                "N2", f"{name}_N2", Point(cells[0], cells[0], number(cell), "", printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
