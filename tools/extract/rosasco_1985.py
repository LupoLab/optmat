# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Rosasco & Hurst (1985): electronic 3 sigma/24 at 1 amagat (Table IV).

Table IV (``tbl_61c221923f79``, published only as JSON cells; the cells agree with the
table image) lists 3 sigma/24 in 10^-18 esu at 1 amagat and 296 K: this work (Ar), the
paper's reevaluations of the three-wave-mixing (CARS) data of Rado (Ref. 5), De Martini
et al. (Ref. 36) and Hauchecorne et al. (Ref. 37), Lundeen et al.'s CARS (Ref. 38), Huo
and Jaffe's theory (Ref. 2), a reevaluation of Hauchecorne et al.'s FISHG (Ref. 37),
Shelton and Buckingham's FISHG (Ref. 6), dc Kerr (Refs. 7, 39, 9, 8 by footnote) and THG
(Ref. 8). Values in parentheses are the normalization references of their column
(note a); they are not extracted. Footnote letters are stripped. The paper's own ac
Kerr results (Table III) are text values in the record.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Rosasco_1985_MeasurementResonantNonresonant"
TABLE = "tbl_61c221923f79"
GASES = {
    "Ar": "Ar",
    "N_2": "N2",
    "He": "He",
    "D_2": "D2",
    "H_2": "H2",
    "O_2": "O2",
    "CH_4": "CH4",
    "SF_6": "SF6",
}
# column -> (dataset prefix, fundamental or pump wavelength in nm): the wavelengths
# are those of the original experiments (see the record's notes), not printed here.
COLUMNS = {
    1: ("this_work", "676.4"),
    2: ("twm_rado", "694.3"),
    3: ("twm_demartini", "694.3"),
    4: ("twm_hauchecorne", "694.3"),
    5: ("twm_lundeen", "532"),
    6: ("twm_theory_huo_jaffe", "694.3"),
    7: ("fishg_hauchecorne", "694.3"),
    8: ("fishg_shelton_buckingham", "514.5"),
    9: ("dc_kerr", "632.8"),
    10: ("thg", "694.3"),
}
_CELL = re.compile(r"^(?P<value>[\d.]+)(?:\^(?P<mark>[a-z]))?$")


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and column (reference values in parentheses skipped)."""
    files = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        for column, (prefix, wavelength) in COLUMNS.items():
            cell = cells[column].strip()
            if not cell or cell.startswith("("):  # none, or a normalization reference
                continue
            match = _CELL.match(cell)
            if match is None:
                msg = f"unexpected cell {cell!r}"
                raise ValueError(msg)
            files[f"{gas}/rosasco_1985__{prefix}_{gas}.csv"] = [
                Point("3sigma/24", wavelength, number(match["value"]), "", printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
