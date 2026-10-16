# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bejot et al. (2010): the n2 of N2, O2 and Ar used in their simulations.

Table I (``tbl_994ad741fed5``) prints n2 ... n8 of N2, O2 and air, Table II
(``tbl_c5fc8d94f60e``, captioned in the library with the text of Fig. 3; the table's
own caption became its note) n2 ... n10 of Ar, all "as used in the present work [17]",
i.e. Loriot et al. (2009) with the erratum's division by (j + 1). Only the n2 columns
are third-order values; they are extracted as reproduced datasets. Cells are
``1.1 \pm 0.2`` (LaTeX); the table images match the CSVs. Air is skipped.

The first library copy (HAL author version) is used; its tables equal those of the
published copy ``Bejot_2010_HigherOrderKerr_c718b7`` (images compared), whose CSVs are
garbled by the extraction.
"""

from _library import Point, number, table, write

DIRECTORY = "Bejot_2010_HigherOrderKerr"
TABLE_I = "tbl_994ad741fed5"
TABLE_II = "tbl_c5fc8d94f60e"


def _plus_minus(cell: str) -> tuple[str, str]:
    r"""Split a printed ``a \pm b`` cell into its two numbers, as printed."""
    value, sigma = cell.split(r"\pm")
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return the n2 column of Tables I (N2, O2) and II (Ar)."""
    files: dict[str, list[Point]] = {}
    header, *rows = table(DIRECTORY, TABLE_I)
    assert header[1].startswith("n_{2}"), header
    for cells in rows:
        gas = {"N_{2}": "N2", "O_{2}": "O2"}.get(cells[0])
        if gas is None:  # air
            continue
        value, sigma = _plus_minus(cells[1])
        files[f"{gas}/bejot_2010__n2_{gas}.csv"] = [
            Point("n2", "800", value, sigma, " | ".join(cells))
        ]
    header, cells = table(DIRECTORY, TABLE_II)
    assert header[0].startswith("n_{2}"), header
    value, sigma = _plus_minus(cells[0])
    files["Ar/bejot_2010__n2_Ar.csv"] = [
        Point("n2", "800", value, sigma, " | ".join(cells))
    ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
