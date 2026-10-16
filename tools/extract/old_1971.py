# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Old, Gentili & Peck (1971): CO2 refractivity at 0.48-1.82 um (Table I).

The library holds this paper as ``Old_1971_DispersionCarbonDioxide`` (identity now
validated; it was unverified when this was written); its text is
"Dispersion of Carbon Dioxide" by J. G. Old, K. L. Gentili and E. R. Peck (University
of Idaho, received 29 June 1970).

Table I (``tbl_02a196e116e3``, checked against the image) gives 19 vacuum wavelengths
(um), the experimental and the predicted (Sellmeier) (n - 1) x 10^7 and their
difference. The experimental and the predicted columns are written to separate files;
two predicted values are misprinted (errata in the source record).
"""

from _library import Point, number, table, write

DIRECTORY = "Old_1971_DispersionCarbonDioxide"
TABLE = "tbl_02a196e116e3"


def extract() -> dict[str, list[Point]]:
    """Return the experimental and the predicted columns."""
    experiment = []
    predicted = []
    for cells in table(DIRECTORY, TABLE)[1:]:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        experiment.append(Point(wavelength, wavelength, number(cells[1]), "", printed))
        predicted.append(Point(wavelength, wavelength, number(cells[2]), "", printed))
    return {
        "CO2/old_1971__experiment.csv": experiment,
        "CO2/old_1971__predicted.csv": predicted,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
