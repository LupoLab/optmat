# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Donley & Shelton (1993): ESHG ratios against density, 514.5 nm.

Table 2 (``tbl_6498821f87d2``, JSON cells; image checked) gives gamma_A/gamma_B and
rho_B/rho_A for the pairs H2/He, Ar/H2 and N2/H2, measured with four electrode arrays,
that is at four phase-match densities (Table 1, ``tbl_81eec7b9dfdc``: rho(H2) = 80.5,
152.4, 202.3 and 323.5 mol m^-3 for arrays 1-4). Table 3 (``tbl_d50d62a52000``, CSV;
image checked) gives the zero-density intercepts (gamma_A/gamma_B)_0 of the linear fits
and the ratios calculated from ab initio gamma.

Files: the gamma ratios per array (``gamma_<A>_over_<B>_arrays``, gas A, reference B),
the density ratios per array (``density_<B>_over_<A>``, filed under gas B with reference
A, so that the printed rho_B/rho_A is the package's rho_gas/rho_reference), and the
Table 3 intercepts and calculated ratios. Row ids name the array and its H2
phase-match density; every row is at 514.5 nm (19430 cm^-1, the table caption).
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Donley_1993_HyperpolarizabilitiesMeasuredInteracting"
TABLE1 = "tbl_81eec7b9dfdc"
TABLE2 = "tbl_6498821f87d2"
TABLE3 = "tbl_d50d62a52000"
SPECIES = {"H_2": "H2", "He": "He", "Ar": "Ar", "N_2": "N2"}
WAVENUMBER = "19430"


def _split(cell: str) -> tuple[str, str]:
    value, error = cell.split("\N{PLUS-MINUS SIGN}")
    return number(value), number(error)


def extract() -> dict[str, list[Point]]:
    """Return the per-array ratios and the Table 3 ratios."""
    densities = {row[0]: row[5] for row in table(DIRECTORY, TABLE1)[1:]}
    labels = [f"array {n} (rho_H2 = {densities[str(n)]} mol/m3)" for n in (1, 2, 3, 4)]
    grid = json_table(DIRECTORY, TABLE2)
    if not grid[1][0].startswith("\\gamma") or not grid[4][0].startswith("\\rho"):
        msg = f"Table 2 layout changed: {grid[1][0]!r}, {grid[4][0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for index, cells in enumerate(grid[1:], 1):
        a, b = SPECIES[cells[1]], SPECIES[cells[2]]
        if index <= 3:
            path = f"{a}/donley_1993__gamma_{a}_over_{b}_arrays.csv"
        else:
            path = f"{b}/donley_1993__density_{b}_over_{a}.csv"
        files[path] = [
            Point(label, WAVENUMBER, *_split(cell), " | ".join([*cells[1:3], cell]))
            for label, cell in zip(labels, cells[3:7], strict=True)
        ]
    for cells in table(DIRECTORY, TABLE3)[1:]:
        a, b = SPECIES[cells[0]], SPECIES[cells[1]]
        printed = " | ".join(cells)
        files[f"{a}/donley_1993__gamma_{a}_over_{b}_zero_density.csv"] = [
            Point("zero density", WAVENUMBER, *_split(cells[2]), printed)
        ]
        calc = cells[3].split(",")[0].rstrip("abcd")  # "19.550a,b)" -> "19.550"
        files[f"{a}/donley_1993__gamma_{a}_over_{b}_calc.csv"] = [
            Point("ab initio", WAVENUMBER, number(calc), "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
