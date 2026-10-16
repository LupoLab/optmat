# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Gaiser & Fellmuth (2010): the static molar polarizability of neon by DCGT.

Table 2 (``tbl_7249e8b90cbe``; image checked, the cells agree) gives, per isotherm
temperature T_i (26-36 K), the molar polarizability A_eps = N_A alpha_0 / (3 eps0) in
``x 10^7 / (m^3/mol)`` (printed values times 1e-7 m^3/mol), the atomic polarizability
alpha_0 in atomic units derived from it with CODATA 2006, and the relative standard
uncertainty u_r in ppm; the last row is the weighted mean, the paper's result.

Files: the weighted mean of A_eps (row "mean"), the seven per-temperature A_eps (rows by
temperature), and the weighted mean of alpha_0 (a conversion of the same number). The
u_r column is not a data column: the record carries it as relative uncertainties. All
values are static (wavelength ``inf``).
"""

from _library import Point, number, table, write

DIRECTORY = "Gaiser_2010_ExperimentalBenchmarkValue"
TABLE_2 = "tbl_7249e8b90cbe"
PREFIX = "Ne/gaiser_2010__"


def extract() -> dict[str, list[Point]]:
    """Return the weighted means and the per-temperature values of Table 2."""
    per_temperature = []
    mean: list[Point] = []
    alpha: list[Point] = []
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        printed = " | ".join(cells)
        if cells[0].strip() == "Weighted mean":
            mean.append(Point("mean", "inf", number(cells[1]), "", printed))
            alpha.append(Point("mean", "inf", number(cells[2]), "", printed))
        else:
            row = number(cells[0])  # the temperature, K
            per_temperature.append(Point(row, "inf", number(cells[1]), "", printed))
    return {
        f"{PREFIX}molar_polarizability.csv": mean,
        f"{PREFIX}molar_polarizability_by_temperature.csv": per_temperature,
        f"{PREFIX}alpha_au.csv": alpha,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
