# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bishop (1990, review): the H2 gamma values of Tables VI and X.

Table VI (``tbl_9c908a3a7ef7``, CSV; image checked): total (electronic + vibrational)
gamma_zzzz of H2 (v = 0) by four methods, all quoted from other papers. Table X
(``tbl_b41995b8424a``, CSV; image checked): electronic and vibrational gamma_par and
gamma_perp of H2 at infinite wavelength and at 4880 A for the dc Kerr effect and ESHG,
'From Bishop et al. (1990)' (Bishop, Pipin and Rerat). Table X's rows mix processes
("4880-Kerr", "4880-ESHG"); each (process, component, part) is its own file, the static
row going to both processes.
"""

from _library import Point, table, write

DIRECTORY = "Bishop_1990_MolecularVibrationalRotational"
TABLE_VI = "tbl_9c908a3a7ef7"
TABLE_X = "tbl_b41995b8424a"
PREFIX = "H2/bishop_1990__"
COLUMNS = {1: ("par", "e"), 2: ("par", "v"), 3: ("perp", "e"), 4: ("perp", "v")}


def extract() -> dict[str, list[Point]]:
    """Return Table VI's values and Table X's columns by process."""
    files: dict[str, list[Point]] = {}
    grid = table(DIRECTORY, TABLE_VI)
    for cells in grid[1:]:
        value = cells[3].removesuffix("^b").removesuffix("^{b}")
        float(value)
        files.setdefault(f"{PREFIX}table6.csv", []).append(
            Point(cells[1], "inf", value, "", " | ".join(cells))
        )
    grid = table(DIRECTORY, TABLE_X)
    for cells in grid[1:]:
        label = cells[0]
        if label == "\\infty":
            processes, wavelength = ("dc_kerr", "eshg"), "inf"
        elif label == "4880-Kerr":
            processes, wavelength = ("dc_kerr",), "4880"
        elif label == "4880-ESHG":
            processes, wavelength = ("eshg",), "4880"
        else:
            msg = f"unexpected Table X row {cells!r}"
            raise ValueError(msg)
        for process in processes:
            for column, (component, part) in COLUMNS.items():
                value = cells[column]
                float(value)
                files.setdefault(
                    f"{PREFIX}{process}_{component}_{part}.csv", []
                ).append(Point(label, wavelength, value, "", " | ".join(cells)))
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
