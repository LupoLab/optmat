# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Griesmann & Burnett (1999): N2 refractivity at 149-279 nm (Table 1).

Table 1 (``tbl_0e295af4655d``, checked against the image) gives, for 31 Fe-Ne
hollow-cathode lines, the vacuum wavelength (nm), the vacuum wave number (um^-1) and
the refractivity x 10^6 at 273.15 K and 101 325 Pa. The library garbles the unit row
('(mm21)' is um^-1, '(310^6)' is x 10^6); that row is skipped. The wave-number column
repeats 1000/wavelength (to 5e-7 um^-1) and is not extracted.
"""

from _library import Point, number, table, write

DIRECTORY = "Griesmann_1999_RefractivityNitrogenGas"
TABLE = "tbl_0e295af4655d"


def extract() -> dict[str, list[Point]]:
    """Return the 31 refractivities."""
    points = [
        Point(number(c[0]), number(c[0]), number(c[2]), "", " | ".join(c))
        for c in table(DIRECTORY, TABLE)[2:]
    ]
    return {"N2/griesmann_1999__experiment.csv": points}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
