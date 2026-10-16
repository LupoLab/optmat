# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Shaw, Hooker & Wilson (1993): n2 relative to air at 248 nm (Table 1).

Table 1 (``tbl_d1f5f8c27975``; image checked, cells equal) prints the relative n2 of
air (1.00, the reference), N2, O2, Ar, CH4, Ne and He with "the standard errors in the
slopes of the fitted lines". Air, the reference, is not a curated species and is not
extracted; every other row is a ratio to air.
"""

from _library import Point, number, table, write

DIRECTORY = "Shaw_1993_MeasurementNonlinearRefractive"
TABLE = "tbl_d1f5f8c27975"
GASES = {
    "N_{2}": "N2",
    "O_{2}": "O2",
    "Ar": "Ar",
    "CH_{4}": "CH4",
    "Ne": "Ne",
    "He": "He",
}


def extract() -> dict[str, list[Point]]:
    """Return one file per gas with its ratio to air at 248.6 nm."""
    header, *rows = table(DIRECTORY, TABLE)
    assert header == ["Gas", "Relative n_{2}"], header
    files: dict[str, list[Point]] = {}
    for cells in rows:
        if cells[0] == "Air":
            continue
        gas = GASES[cells[0]]
        value, sigma = cells[1].split("±")
        files[f"{gas}/shaw_1993__ratio_{gas}.csv"] = [
            Point("248", "248.6", number(value), number(sigma), " | ".join(cells))
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
