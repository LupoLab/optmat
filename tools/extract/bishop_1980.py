# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Bishop & Cheung (1980): isotropic dynamic polarizability of H2 (Table V).

Table V (``tbl_7b56aafde20c``, checked against the image) lists 22 vacuum wavelengths
(A; the static limit printed as an infinity sign with a footnote mark 'c') with the
isotropic polarizability (a0^3) of H2 for J = 0-4 (v = 0), its thermal average at
293 K, and an 'Experimental' column from Victor & Dalgarno (1969) whose values in
parentheses are that work's semi-empirical values. The parenthesised values are
written to a separate file. Tables VI (anisotropies), VII-VIII (HeH+) and IX (static
components) are not extracted.
"""

from _library import Point, number, table, write

DIRECTORY = "Bishop_1980_DynamicDipolePolarizability"
TABLE_V = "tbl_7b56aafde20c"


def _wavelength(cell: str) -> str:
    cell = cell.strip()
    return "inf" if cell.startswith("\\infty") else number(cell)


def extract() -> dict[str, list[Point]]:
    """Return the per-J, thermal and Victor-Dalgarno columns."""
    rows = table(DIRECTORY, TABLE_V)[1:]
    files: dict[str, list[Point]] = {}
    for j in range(5):
        files[f"H2/bishop_1980__alpha_J{j}.csv"] = [
            Point(
                _wavelength(c[0]),
                _wavelength(c[0]),
                number(c[j + 1]),
                "",
                " | ".join(c),
            )
            for c in rows
        ]
    files["H2/bishop_1980__alpha_293K.csv"] = [
        Point(_wavelength(c[0]), _wavelength(c[0]), number(c[6]), "", " | ".join(c))
        for c in rows
    ]
    experimental: list[Point] = []
    semiempirical: list[Point] = []
    for c in rows:
        cell = c[7].strip()
        point = Point(
            _wavelength(c[0]),
            _wavelength(c[0]),
            number(cell.strip("()")),
            "",
            " | ".join(c),
        )
        (semiempirical if cell.startswith("(") else experimental).append(point)
    files["H2/bishop_1980__victor_dalgarno.csv"] = experimental
    files["H2/bishop_1980__victor_dalgarno_semiempirical.csv"] = semiempirical
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
