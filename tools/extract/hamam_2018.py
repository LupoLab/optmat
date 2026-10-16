# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Hamam & Gamal (2018): Liu et al.'s silica breakdown thresholds (reproduced).

Table 2 (``tbl_5037c0738f9f``, JSON cells only: merged headers) lists, for nine pulse
durations (s) and two focusing conditions, the threshold intensity I_th (W/cm2)
"measured ... by Liu et al., (2005)" and the laser power P (W) the authors calculate
from it. Checked against the table image: the cells agree. Only I_th is extracted (there
is no quantity for the power at threshold; P is quoted in the record, including the cell
``1,13E04``, which the image prints with a decimal comma). Table 1
(``tbl_9e3c689c1c99``) repeats the same I_th beside the Keldysh parameter; the script
checks that the two tables' I_th agree cell for cell.

Both tables head the second focusing condition "NA = 0.225" (the image confirms it),
while the abstract, the text and every figure caption give 0.255: see the record.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Hamam_2018_NumericalAnalysisBreakdown"
TABLE1 = "tbl_9e3c689c1c99"
TABLE2 = "tbl_5037c0738f9f"


def extract() -> dict[str, list[Point]]:
    """Return the threshold intensities for the two numerical apertures."""
    rows2 = json_table(DIRECTORY, TABLE2)[2:]
    rows1 = json_table(DIRECTORY, TABLE1)[2:]
    if [(r[0], r[1], r[3]) for r in rows1] != [(r[0], r[1], r[3]) for r in rows2]:
        msg = "Tables 1 and 2 disagree on the durations or I_th"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for stem, column in (("NA0126", 1), ("NA0255", 3)):
        files[f"SiO2/hamam_2018__threshold_{stem}.csv"] = [
            Point(
                cells[0],
                number(cells[0]),
                number(cells[column]),
                "",
                " | ".join(cells),
            )
            for cells in rows2
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
