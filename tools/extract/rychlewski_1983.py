# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Rychlewski (1983): ab initio dynamic polarizabilities of D2 and H2.

All tables list 16 vacuum wavelengths in angstroms (``inf`` for the static limit,
printed as the infinity sign) with polarizabilities in atomic units (a0^3), checked
against the table images.

- Table V (``tbl_01d4df32ff6f``): thermally averaged (293 K) isotropic polarizability
  alpha and anisotropy gamma of HD and D2. Only alpha(D2) is extracted (HD is not a
  curated species; the package has no anisotropy quantity).
- Table II (``tbl_5741c5d89712``): alpha of D2 in the rotational states J = 0-6 (v = 0),
  one file per J.
- Table VI (``tbl_895c7be743f8``): H2, thermally averaged; only the author's own column
  (a) is extracted. Column b (Bishop & Cheung 1980) and the 'experimental' column c
  (Victor & Dalgarno 1969) are the same values as in Bishop & Cheung's own tables,
  extracted from that source.
Tables I and III (HD), IV (D2 anisotropies) and VII (non-diagonal matrix elements) are
not extracted.
"""

from _library import Point, number, table, write

DIRECTORY = "Rychlewski_1983_FrequencyDependentPolarizabilities"
TABLE_II = "tbl_5741c5d89712"
TABLE_V = "tbl_01d4df32ff6f"
TABLE_VI = "tbl_895c7be743f8"


def _wavelength(cell: str) -> str:
    cell = cell.strip()
    return "inf" if cell in {"\N{INFINITY}", "$\\infty$", "\\infty"} else number(cell)


def _column(table_id: str, column: int) -> list[Point]:
    points = []
    for cells in table(DIRECTORY, table_id)[1:]:
        wavelength = _wavelength(cells[0])
        points.append(
            Point(wavelength, wavelength, number(cells[column]), "", " | ".join(cells))
        )
    return points


def extract() -> dict[str, list[Point]]:
    """Return alpha(D2) thermal and per J, and alpha(H2) thermal."""
    files = {
        "D2/rychlewski_1983__alpha_293K.csv": _column(TABLE_V, 3),
        "H2/rychlewski_1983__alpha_293K.csv": _column(TABLE_VI, 1),
    }
    for j in range(7):
        files[f"D2/rychlewski_1983__alpha_J{j}.csv"] = _column(TABLE_II, j + 1)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
