# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Borzsonyi et al. (2010): Table 1, n2 per bar of Ar, Ne, Xe, N2 (p. 25853).

Row 1 prints the measured "Pressure dependent n2 [x 10^-19 cm2/W at 1bar]" as ``value
\pm error``; row 2 the literature values at 1 bar, several per gas, each with its
reference number as a superscript, e.g. ``0.98^{[16]}``. Air (a column) is not curated.
The CSV matches the table image (checked).
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Borzsonyi_2010_MeasurementPressureDependent"
TABLE = "tbl_7c2c2d4e525c"
GASES = {"Argon": "Ar", "Neon": "Ne", "Xenon": "Xe", "Nitrogen": "N2"}
LITERATURE = re.compile(r"([\d.]+)\^\{\[([\d, ]+)\]\}")


def extract() -> dict[str, list[Point]]:
    """Return the measured value and the literature values of each curated gas."""
    header, measured, literature = table(DIRECTORY, TABLE)
    files = {}
    for column, name in enumerate(header):
        if name not in GASES:
            continue  # the label column and air
        gas = GASES[name]
        printed = f"{name} | {measured[column]} | {literature[column]}"
        value, error = measured[column].split(r"\pm")
        files[f"{gas}/borzsonyi_2010__n2_{gas}.csv"] = [
            Point("800", "800", number(value), number(error), printed)
        ]
        # One row per literature value; the row id names the reference(s).
        files[f"{gas}/borzsonyi_2010__literature_{gas}.csv"] = [
            Point(
                "ref" + refs.replace(" ", "").replace(",", "_"),
                "not stated",
                value_,
                "",
                printed,
            )
            for value_, refs in LITERATURE.findall(literature[column])
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
