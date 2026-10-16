# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton and Mizrahi (1985): the SF6/N2 polarizability-dispersion ratio.

The only table (``tbl_72ec8b2e57c2``; image checked) gives, at the fundamental
wavelength (nm; dye- and Ar+-laser lines, air wavelengths), the ratio
``Delta alpha_SF6 / Delta alpha_N2`` (Delta alpha = alpha(2 omega) - alpha(omega)),
which equals the phase-match density ratio rho_N2 / rho_SF6, with its +-, and the
hyperpolarizability ratios (not extracted). The 694.3 nm row (a gamma value from Ward
and Miller only, the ratio printed as a dash) is skipped.

Third-order data: the same table's gamma_SF6/gamma_N2
(``gamma_SF6_over_N2``, filed under SF6 with reference N2) and gamma_SF6 in 1e-62 C^4
m^4 J^-3 (``gamma_SF6``, the six rows measured here; ``gamma_SF6_ward_miller``, the
694.3 nm row marked c) "From ref. [12]", Ward & Miller 1979).
"""

from _library import Point, number, table, write

DIRECTORY = "Shelton_1985_FrequencyDependenceHyperpolarizability"
TABLE = "tbl_72ec8b2e57c2"


def extract() -> dict[str, list[Point]]:
    """Return the dispersion ratios and the gamma data as data files."""
    points = []
    for cells in table(DIRECTORY, TABLE)[1:]:
        if cells[1].strip() == "-":
            continue
        value, error = cells[1].split("\N{PLUS-MINUS SIGN}")
        wavelength = number(cells[0])
        points.append(
            Point(
                wavelength, wavelength, number(value), number(error), " | ".join(cells)
            )
        )
    files = {"N2/shelton_1985_sf6__density_N2_over_SF6.csv": points}
    ratios, gammas, ward_miller = [], [], []
    for cells in table(DIRECTORY, TABLE)[1:]:
        printed = " | ".join(cells)
        wavelength = number(cells[0])
        gamma, gamma_error = cells[3].removesuffix("c)").split("\N{PLUS-MINUS SIGN}")
        point = Point(
            wavelength, wavelength, number(gamma), number(gamma_error), printed
        )
        if cells[3].endswith("c)"):
            ward_miller.append(point)
            continue
        gammas.append(point)
        ratio, ratio_error = cells[2].split("\N{PLUS-MINUS SIGN}")
        ratios.append(
            Point(wavelength, wavelength, number(ratio), number(ratio_error), printed)
        )
    files["SF6/shelton_1985_sf6__gamma_SF6_over_N2.csv"] = ratios
    files["SF6/shelton_1985_sf6__gamma_SF6.csv"] = gammas
    files["SF6/shelton_1985_sf6__gamma_SF6_ward_miller.csv"] = ward_miller
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
