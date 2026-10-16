# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Lehmeier, Leupacher & Penzkofer (1985): Tables 1 and 2.

Table 1 (``tbl_3d2f3ba324ee``) lists per gas the standard (1 atm, 0 degrees C) n - 1 at
the fundamental lambda_1 = 1.055 um and at its third harmonic, both x 10^5, each with a
footnote letter giving its source (``3.472a)``; a) Dalgarno and Kingston, b)
Landolt-Boernstein), which is stripped. The third-harmonic wavelength is not printed; it
is lambda_1 / 3 = 0.3516667 um. The same table gives the paper's third-harmonic results:
chi_0(3) / chi_0(3)(He) (column 5), xi/f^2 relative to He (column 6, an anharmonic
coupling constant, not extracted) and gamma_1111(3)(-w3; w1, w1, w1) in
10^-64 C m^4/V^3 (column 7) and 10^-39 esu (column 8). Delta k_0 and the coherence
length (columns 3, 4) follow from the refractivities and are not extracted.

Table 2 (``tbl_905211afbf03``, published only as JSON cells) compares gamma(3) at the
ruby wavelength 694.3 nm, in 10^-39 esu: (a) this work, transferred from 1.055 um by
the paper's Eq. (18); (b) Ward and New's ranges (not single values: not extracted);
(c) Ward and New normalised to gamma(He) = 1 x 10^-39 esu; (d)-(i) theories; and the
last column, the dc Kerr gamma(-w; w, 0, 0) of Buckingham and Dunmur at 632.8 nm. The
cells were compared with the table image: they agree, including the blank (f) cells
below helium, which the image also leaves blank, and the dashes.
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Lehmeier_1985_NonresonantThirdOrder"
TABLE = "tbl_3d2f3ba324ee"
TABLE_2 = "tbl_905211afbf03"
GASES = {
    "He": "He",
    "Ne": "Ne",
    "Ar": "Ar",
    "Kr": "Kr",
    "Xe": "Xe",
    "N_{2}": "N2",
    "N_2": "N2",
}
FUNDAMENTAL = "1.055"  # um, Table 1 caption
HARMONIC = "0.3516667"  # um, lambda_1 / 3 (not printed)
RUBY = "694.3"  # nm, Table 2 caption
HE_NE = "632.8"  # nm, p. 72: the dc Kerr measurements' wavelength

# Table 2 column index -> dataset prefix (column 3, (b), holds ranges: not extracted).
TABLE_2_COLUMNS = {
    2: "gamma_694",
    4: "ward_new_normalized",
    5: "sitz_yaris",
    6: "mizuno",
    7: "klingbeil",
    8: "dawes",
    9: "tripathi",
    10: "langhoff",
    11: "dc_kerr_buckingham",
}


def extract() -> dict[str, list[Point]]:
    """Return the refractivity and third-order data files."""
    files = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        files[f"{gas}/lehmeier_1985__standard_{gas}.csv"] = [
            Point("n1", FUNDAMENTAL, number(cells[1][:-2]), "", printed),
            Point("n3", HARMONIC, number(cells[2][:-2]), "", printed),
        ]
        files[f"{gas}/lehmeier_1985__gamma_thg_{gas}.csv"] = [
            Point("gamma", FUNDAMENTAL, number(cells[7]), "", printed)
        ]
        files[f"{gas}/lehmeier_1985__gamma_thg_esu_{gas}.csv"] = [
            Point("gamma", FUNDAMENTAL, number(cells[8]), "", printed)
        ]
        if gas != "He":  # the helium row is the reference, 1 by definition
            files[f"{gas}/lehmeier_1985__chi_ratio_{gas}.csv"] = [
                Point("ratio", FUNDAMENTAL, number(cells[5]), "", printed)
            ]
    for cells in json_table(DIRECTORY, TABLE_2)[2:]:
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        for column, prefix in TABLE_2_COLUMNS.items():
            cell = cells[column]
            if cell in {"", "-"}:  # blank or a dash: no value for this gas
                continue
            wavelength = HE_NE if prefix == "dc_kerr_buckingham" else RUBY
            files[f"{gas}/lehmeier_1985__{prefix}_{gas}.csv"] = [
                Point("gamma", wavelength, number(cell), "", printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
