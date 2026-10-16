# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Wang et al. (2015): Table 2, the N2 constants of the simulations (p. 2242).

Table 2 lists the values used as the "true" parameters of the numerical experiments:
n2,instant at 800 and 400 nm (um2/(TW bar)), B0, D0 (cm^-1; not extracted), Delta alpha
(angstrom3) and the spin weights (not extracted). All are quoted or derived from other
works (footnotes a-c). The CSV matches the table image (checked).
"""

from _library import Point, number, table, write

DIRECTORY = "Wang_2015_MeasuringElectronicKerr"
TABLE = "tbl_1bb9ec2e9554"


def extract() -> dict[str, list[Point]]:
    """Return n2,instant (two wavelengths) and Delta alpha of N2."""
    rows = {cells[0]: cells for cells in table(DIRECTORY, TABLE)[1:]}
    n2 = []
    for wavelength, footnote in (("800", "a"), ("400", "b")):
        key = (
            rf"n_{{2,\text{{instant}}}}^{{{footnote}}} at {wavelength} nm"
            r"(\mu\text{m}^2/(TWbar))"
        )
        cells = rows[key]
        n2.append(
            Point(wavelength, wavelength, number(cells[1]), "", " | ".join(cells))
        )
    dalpha = rows[r"\Delta\alpha^{c} (A3)"]
    return {
        "N2/wang_2015__n2_instant_N2.csv": n2,
        "N2/wang_2015__dalpha_herzberg_N2.csv": [
            Point("inf", "not stated", number(dalpha[1]), "", " | ".join(dalpha))
        ],
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
