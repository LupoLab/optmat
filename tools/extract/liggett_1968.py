# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Liggett & Levinger (1968): mu_-4 and n - 1 at Lyman alpha (Tables II, III).

Table II (``tbl_699c13f7b3b4``, no caption in the library; the caption "TABLE II.
Calculations and measurements of mu_-4" is printed in the text above it) gives the
fourth inverse moment mu_-4 of Ne and Ar in Rydberg units from four sources; Table III
gives delta(E) = n - 1 (x 10^4) at the Lyman alpha line from five sources (an em dash
is no value). Each column becomes one dataset per gas.

Tables IV and V (neon and argon in the x-ray region, delta(E) x 10^8, all negative) are
not extracted: a negative refractivity cannot be expressed as a (positive) dataset in
the package; they are described in the source record.
"""

from _library import Point, number, table, write

DIRECTORY = "Liggett_1968_CalculationIndexRefraction"
TABLE_II = "tbl_699c13f7b3b4"
TABLE_III = "tbl_7b2a63bc8b39"
GASES = {"Neon": "Ne", "Argon": "Ar"}
MU_COLUMNS = ("cuthbertson", "ingersoll", "dalgarno", "theory")
LYMAN_COLUMNS = ("marmo", "dalgarno", "heddle_1963", "kingston_1964", "theory")
# The Lyman alpha line of atomic hydrogen, which the paper names without a wavelength
# ("the Lyman alpha line (0.75 Ry)", p. 1); the value is Mahon et al.'s (1979, p. 4).
LYMAN_ALPHA = "1215.67"


def extract() -> dict[str, list[Point]]:
    """Return one file per column and gas of Tables II and III."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_II)[1:]:
        gas = GASES[cells[0].removeprefix(r"\mu_{-4}")]
        for name, cell in zip(MU_COLUMNS, cells[1:], strict=True):
            files[f"{gas}/liggett_1968__mu4_{name}_{gas}.csv"] = [
                Point("-4", "-4", number(cell), "", " | ".join(cells))
            ]
    for cells in table(DIRECTORY, TABLE_III)[1:]:
        gas = GASES[cells[0].removeprefix("δ(E)\N{MULTIPLICATION SIGN}104")]
        for name, cell in zip(LYMAN_COLUMNS, cells[1:], strict=True):
            if cell == "—":
                continue
            files[f"{gas}/liggett_1968__lyman_{name}_{gas}.csv"] = [
                Point(LYMAN_ALPHA, LYMAN_ALPHA, number(cell), "", " | ".join(cells))
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
