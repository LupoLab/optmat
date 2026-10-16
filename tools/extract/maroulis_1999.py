# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Maroulis (1999): mean static polarizability of ozone (Table III).

The library holds this paper as ``Maroulis_1999_ExtremeProblemRevisited`` (identity now
validated; it was unverified when this was written); its text is
"An extreme problem revisited: The hyperpolarizability of the open and cyclic forms of
ozone" by G. Maroulis. Table III (``tbl_20d61f051169``, JSON cells only; checked
against the image) compares the open and cyclic forms. Only the mean polarizability of
the open (ground-state) form is extracted: this paper's CCSD(T) value and the
'Experiment' value (footnote e, an extrapolation of old refractivity data to infinite
wavelength, U. Hohm, private communication). The CCD+ST(CCD) column repeats Maroulis
(1994) and the cyclic form is not the gas's ground state.

Third-order data: the gamma rows of Table I (``tbl_f41afb24c0cd``, open form,
methods SCF to CCSD(T)); one file per row (six components and the mean gamma-bar), row
ids the methods. STRUCTURAL FIX: the library's cells render the minus sign as a leading
"2" in some cells (and not others); the image shows: gamma_xxxx -2962 (cell "22962"),
-14832 ("214832"), -1217 ("21217"), -294 ("2294"); gamma_zzxx (MP3) -47 ("247");
gamma-bar (SDQ-MP4) -1282 ("21282"). Those cells are mapped below, each asserting the
garbled text; every mean is consistent with its components (checked). Table II (the
cyclic D3h form, not the gas's ground state) is not taken.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Maroulis_1999_ExtremeProblemRevisited"
TABLE_III = "tbl_20d61f051169"
COLUMNS = {2: "ccsd_t", 3: "hohm_refractivity"}
TABLE_I = "tbl_f41afb24c0cd"
GAMMA_ROWS = {
    "gxxxx": "xxxx",
    "gyyyy": "yyyy",
    "gzzzz": "zzzz",
    "gxxyy": "xxyy",
    "gyyzz": "yyzz",
    "gzzxx": "zzxx",
    "g\u00af": "mean",
}
METHODS = ("SCF", "MP2", "MP3", "SDQ-MP4", "MP4", "CCSD", "CCSD(T)")
# (row, method) -> (library cell, value on the image): lost minus signs.
MINUS = {
    ("xxxx", "SCF"): ("22962", "-2962"),
    ("xxxx", "SDQ-MP4"): ("214832", "-14832"),
    ("xxxx", "MP4"): ("21217", "-1217"),
    ("xxxx", "CCSD"): ("2294", "-294"),
    ("zzxx", "MP3"): ("247", "-47"),
    ("mean", "SDQ-MP4"): ("21282", "-1282"),
}


def _gamma() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_I)[1:]:
        name = GAMMA_ROWS.get(cells[0])
        if name is None:
            continue
        for method, cell in zip(METHODS, cells[1:], strict=True):
            value = cell
            if (name, method) in MINUS:
                garbled, value = MINUS[(name, method)]
                if cell != garbled:
                    msg = f"Table I {name}/{method}: expected {garbled!r}, got {cell!r}"
                    raise ValueError(msg)
            files.setdefault(f"O3/maroulis_1999__gamma_{name}.csv", []).append(
                Point(method, "inf", number(value), "", " | ".join(cells))
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the two mean-polarizability values, one file each, and the gammas."""
    (row,) = (c for c in json_table(DIRECTORY, TABLE_III) if c[0] == "\\bar{\\alpha}")
    printed = " | ".join(row)
    return _gamma() | {
        f"O3/maroulis_1999__{name}.csv": [
            Point(
                "alpha_bar", "inf", number(re.sub(r"\^.*$", "", row[col])), "", printed
            )
        ]
        for col, name in COLUMNS.items()
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
