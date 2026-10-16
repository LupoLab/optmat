# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Dawes (1968): calculated third-harmonic coefficients at 6943 A (Table I).

Table I (``tbl_f296ada000f8``) lists per rare gas Cuthbertson's oscillator parameters
N_0 and nu_0 (reproduced inputs, not extracted), (A) the calculated alpha_zzzz(3nu) in
10^-39 esu/atom (the caption prints "10^39", a misprint: the abstract gives 0.9 x
10^-39 esu/atom for He), (B) the same relative to He, and (C), (D) New and Ward's two
experimental sets of ratios as ranges (not single values: not extracted). The cells
agree with the table image.
"""

from _library import Point, number, table, write

DIRECTORY = "Dawes_1968_OpticalThirdHarmonic"
TABLE = "tbl_f296ada000f8"
RUBY = "694.3"  # nm, lambda = 6943 A (Table I caption)


def extract() -> dict[str, list[Point]]:
    """Return one file per gas for columns (A) and (B)."""
    files = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        gas = cells[0]
        printed = " | ".join(cells)
        files[f"{gas}/dawes_1968__alpha_3nu_{gas}.csv"] = [
            Point("calc", RUBY, number(cells[3]), "", printed)
        ]
        if gas != "He":  # the helium row of (B) is the reference, 1
            files[f"{gas}/dawes_1968__ratio_{gas}.csv"] = [
                Point("calc", RUBY, number(cells[4]), "", printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
