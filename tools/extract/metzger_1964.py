# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Metzger & Cook (1964): O2 absorption coefficients of bands and windows.

Table 2 (``tbl_cbf80ad2566e``, CSV; image checked, the cells agree): designation, the
wavelength (A) and the absorption coefficient k (cm^-1) of the diffuse bands A-F, the
seven "windows" 1-7 and the maximum of the Schumann-Runge continuum. The paper prints a
raised decimal point ("1215·7", "0·28"); the value columns keep the digits
with a full stop, the printed column the raw cells. The last row's designation is
"S-R maximum" on two lines in the image; the cell lost the space. The bands and windows
go to one file, the Schumann-Runge maximum (the only value with a stated uncertainty) to
another.
"""

from _library import Point, table, write

DIRECTORY = "Metzger_1964_ReinvestigationAbsorptionCross"
TABLE_2 = "tbl_cbf80ad2566e"
PREFIX = "O2/metzger_1964__"


def _number(cell: str) -> str:
    text = cell.strip().replace("\N{MIDDLE DOT}", ".")
    float(text)
    return text


def extract() -> dict[str, list[Point]]:
    """Return Table 2's bands and windows, and its Schumann-Runge maximum."""
    points, maximum = [], []
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        designation, wavelength, k = (c.strip() for c in cells)
        row = "S-R maximum" if designation == "S-Rmaximum" else designation
        (maximum if row == "S-R maximum" else points).append(
            Point(
                row,
                _number(wavelength),
                _number(k),
                "",
                " | ".join((designation, wavelength, k)),
            )
        )
    return {
        f"{PREFIX}bands_windows.csv": points,
        f"{PREFIX}sr_maximum.csv": maximum,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
