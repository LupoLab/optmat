# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Franta et al. (2019): the band gap of the current float-zone silicon model.

Table 2 (p. 12) lists the band gap and critical-point energies of the authors' current
temperature-dependent model of c-Si at 0 K and at the reference temperature 300 K,
with the average phonon temperature Theta of Eq. (42) and the equation of each point's
temperature dependence. Only E_g has a quantity: its 0 K and 300 K values are two
datasets. The other energies are the record's formula coefficients. The CSV agrees with
the table image (checked); the row is found by its label and its model equation.
"""

from _library import Point, number, table, write

DIRECTORY = "Franta_2019_TemperatureDependentDispersion"
TABLE_2 = "tbl_09a36f80965c"
HEADER = [
    "",
    "standard notation",
    "E_{j}^{0\\text{K}} (eV)",
    "E_{j}^{300\\text{K}} (eV)",
]


def extract() -> dict[str, list[Point]]:
    """Return E_g at 0 K and at 300 K as two datasets."""
    rows = table(DIRECTORY, TABLE_2)
    if rows[0][:4] != HEADER:
        msg = f"unexpected Table 2 header {rows[0]!r}"
        raise ValueError(msg)
    (gap,) = [cells for cells in rows if cells[0] == "E_{g}"]
    if gap[5] != "(42)":
        msg = f"unexpected E_g row {gap!r}"
        raise ValueError(msg)
    printed = " | ".join(gap)
    return {
        "Si/franta_2019__band_gap_0K.csv": [
            Point("1", "", number(gap[2]), "", printed)
        ],
        "Si/franta_2019__band_gap_300K.csv": [
            Point("1", "", number(gap[3]), "", printed)
        ],
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
