# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Peck & Khanna (1966): N2 refractivity at 0.47-2.06 um (Tables I and II).

The library holds this paper unverified (``Peck_1966_DispersionNitrogen``); its text is
"Dispersion of Nitrogen" by E. R. Peck and B. N. Khanna, J. Opt. Soc. Am. 56, 1059.

Table I (``tbl_c2970086b7ad``, checked against the image) gives 19 vacuum wavelengths
(A) with a source label (Cd, Hg, Ne, A = argon, Rb, He; the library's 'Hg^{193}' is
the image's 198Hg), the measured refractivity x 10^8 at 15 C standardised to
28 352.40 at 5462.258 A, the value of Eq. (3) and the difference. Digit groups are
separated by a thin space, removed here. Three measured values are printed in
parentheses (with their differences): they were omitted from the fit (p. 1060) and are
written to a separate file.

Table II (``tbl_029c38d597ef``, image checked) gives the Verdet constant at air
wavelengths (A) relative to its value at 5780 A, computed from Eq. (3) and as measured
by Ingersoll & Liebenberg. Its row '5460-61' (the mercury green line as listed by two
sources) is given the wavelength 5460.5 A, the middle of the printed range; the half
angstrom is immaterial for a ratio known to 1e-3.
"""

from _library import Point, number, table, write

DIRECTORY = "Peck_1966_DispersionNitrogen"
TABLE_I = "tbl_c2970086b7ad"
TABLE_II = "tbl_029c38d597ef"


def _number(cell: str) -> str:
    return number(cell.strip().strip("()").replace(" ", ""))


def _table_i() -> dict[str, list[Point]]:
    measured: list[Point] = []
    excluded: list[Point] = []
    fit: list[Point] = []
    for cells in table(DIRECTORY, TABLE_I)[1:]:
        wavelength = number(cells[1])
        printed = " | ".join(cells)
        point = Point(wavelength, wavelength, _number(cells[2]), "", printed)
        (excluded if cells[2].strip().startswith("(") else measured).append(point)
        fit.append(Point(wavelength, wavelength, _number(cells[3]), "", printed))
    return {
        "N2/peck_1966__experiment.csv": measured,
        "N2/peck_1966__experiment_excluded.csv": excluded,
        "N2/peck_1966__eq3.csv": fit,
    }


def _table_ii() -> dict[str, list[Point]]:
    computed: list[Point] = []
    measured: list[Point] = []
    for cells in table(DIRECTORY, TABLE_II)[1:]:
        label = cells[0].strip()
        wavelength = "5460.5" if label == "5460\N{EN DASH}61" else number(label)
        printed = " | ".join(cells)
        computed.append(Point(label, wavelength, number(cells[1]), "", printed))
        measured.append(Point(label, wavelength, number(cells[2]), "", printed))
    return {
        "N2/peck_1966__verdet_ratio_eq3.csv": computed,
        "N2/peck_1966__verdet_ratio_ingersoll_liebenberg.csv": measured,
    }


def extract() -> dict[str, list[Point]]:
    """Return the Table I and Table II data files."""
    return _table_i() | _table_ii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
