# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Brion et al. (1998): O3 cross sections at six visible lines, 295 K.

Table II (``tbl_cc4f708ad271``, CSV; image checked, every cell agrees): the cross
section (sigma x 10^21 cm^2, the caption) at six discrete wavelengths (543.516-632.816
nm) from Anderson and Mauersberger, Hearn and this work, each printed ``value +-
sigma``, and the deviation (derived, not extracted). The paper prints decimal commas
("3,075 +- 0,027", "543,516") and puts a space after the comma in the cross sections
("3, 075"); one cell mixes a comma and a point ("4,634 +- 0.070", 611.971 nm; the image
shows the point too). Commas become points in the value and sigma columns; ``printed``
keeps the raw cells.
"""

from _library import Point, table, write

DIRECTORY = "Brion_1998_AbsorptionSpectraMeasurements"
TABLE_II = "tbl_cc4f708ad271"
PREFIX = "O3/brion_1998__"
COLUMNS = {1: "anderson_mauersberger_1992", 2: "hearn_1961", 3: "present"}


def _decimal(text: str) -> str:
    number = text.strip().replace(" ", "").replace(",", ".")
    float(number)
    return number


def extract() -> dict[str, list[Point]]:
    """Return Table II's three columns of cross sections with their uncertainties."""
    files: dict[str, list[Point]] = {f"{PREFIX}{n}.csv": [] for n in COLUMNS.values()}
    for cells in table(DIRECTORY, TABLE_II)[1:]:
        printed = " | ".join(c.strip() for c in cells)
        wavelength = _decimal(cells[0])
        for column, name in COLUMNS.items():
            cell = cells[column].strip()
            if not cell:
                continue
            value, sigma = (_decimal(p) for p in cell.split("\N{PLUS-MINUS SIGN}"))
            files[f"{PREFIX}{name}.csv"].append(
                Point(wavelength, wavelength, value, sigma, printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
