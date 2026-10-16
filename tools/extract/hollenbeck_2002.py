# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Hollenbeck & Cantrell (2002): the 13-mode Raman response of silica fibre.

Table 1 (p. 2889) gives, per mode i, the component position (cm^-1), the "peak
intensity" A_i, and the "Gaussian FWHM" and "Lorentzian FWHM" (cm^-1) of the
intermediate-broadening model, Eq. (9). The extracted cells agree with the table image
(checked cell by cell); the header's "cm21" is the extraction's reading of cm^-1.
"""

from _library import Point, number, table, write

DIRECTORY = "Hollenbeck_2002_MultipleVibrationalMode"
TABLE = "tbl_3554e41bb9a2"


def extract() -> dict[str, list[Point]]:
    """Return the 13 modes: amplitude, with position and the two widths as columns."""
    rows = table(DIRECTORY, TABLE)[1:]
    modes = []
    for cells in rows:
        mode, position, amplitude, gaussian, lorentzian = cells
        modes.append(
            Point(
                mode,
                mode,
                number(amplitude),
                "",
                " | ".join(cells),
                {
                    "shift": number(position),
                    "gaussian_width": number(gaussian),
                    "lorentzian_width": number(lorentzian),
                },
            )
        )
    return {"SiO2/hollenbeck_2002__modes.csv": modes}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
