# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Ordal et al. (1988): n, k, -epsilon_1 and epsilon_2 of Al, 50-15000 cm-1.

Table I (``tbl_131e600939fe``, p. 1204) prints omega (cm-1), n, k, -epsilon_1, epsilon_2
and R for 51 wavenumbers. It is the aluminium table on the page image, not only by its
caption: the caption "Table I. Optical Constants of Aluminum" sits directly above it on
PDF page 2, in the aluminium section, and the table has aluminium's signature, the
interband dip of R to 0.868 at 12000 cm-1 (1.49 eV) and a Drude wing matching the
section's omega_p ~ 9.5e4 cm-1 (the other five tables' low-frequency n and k match their
own metals' Drude wings). The CSV agrees with the image cell for cell.

The fourth column is printed "-epsilon_1": its values are extracted as printed, and the
dataset's scale of -1 makes them epsilon_1. R has no quantity and stays in ``printed``.

The epsilon columns carry a misprint the record corrects by errata (not here): every
|epsilon| of 1000 or more is printed ten times too small (the decimal point one place
left), e.g. -epsilon_1 "4445.8912" at 50 cm-1 where k^2 - n^2 = 44458.9 from the same
row. The n, k and R columns are mutually consistent throughout (R from n and k to 1e-6).
"""

from _library import Point, number, table, write

DIRECTORY = "Ordal_1988_OpticalPropertiesAl"
TABLE_I = "tbl_131e600939fe"
HEADER = ["ω", "n", "k", "-\\epsilon_1", "\\epsilon_2", "R"]
# Dataset suffix -> column.
COLUMNS = {"n": 1, "k": 2, "eps1": 3, "eps2": 4}


def extract() -> dict[str, list[Point]]:
    """Return aluminium's n, k, -epsilon_1 and epsilon_2 against the wavenumber."""
    rows = table(DIRECTORY, TABLE_I)
    if rows[0] != HEADER:
        msg = f"{TABLE_I}: unexpected header {rows[0]}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {
        f"Al/ordal_1988__{name}_Al.csv": [] for name in COLUMNS
    }
    for cells in rows[1:]:
        wavenumber = number(cells[0])
        printed = " | ".join(cells)
        for name, column in COLUMNS.items():
            files[f"Al/ordal_1988__{name}_Al.csv"].append(
                Point(wavenumber, wavenumber, number(cells[column]), "", printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
