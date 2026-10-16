# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Gerardo & Johnson (1974): Table II, oscillator strengths of Kr and Xe lines.

Table II (``tbl_5fe8d04e80c5``, cells file only; image checked, cells agree) lists, per
resonance line (gas, state, wavelength in angstrom as printed, rounded to 1 A, g), two
literature values of the absorption oscillator strength f_k1 with their reference
number: 22 (Wilkinson 1966), 23 (Anderson 1965), 24 (Lewis 1967). Each line's second
value is on a row with only f and the reference; the line is carried down. Every value
is another paper's: the files are reproductions, one per gas.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Gerardo_1974_PhotoattenuationExtremeRed"
TABLE_II = "tbl_5fe8d04e80c5"
HEADER = ["Gas", "State", "\\lambda (Å)", "g", "f_{k1}", "Reference"]
FILES = {"Xe": "Xe/gerardo_1974__table2_Xe.csv", "Kr": "Kr/gerardo_1974__table2_Kr.csv"}
REFERENCES = {"22", "23", "24"}


def extract() -> dict[str, list[Point]]:
    """Return one file per gas: row <state>_ref<n>, x the printed wavelength (A)."""
    rows = json_table(DIRECTORY, TABLE_II)
    if rows[0] != HEADER:
        msg = f"Table II: unexpected header {rows[0]}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {name: [] for name in FILES.values()}
    gas = state = wavelength = ""
    for row in rows[1:]:
        if row[0]:
            gas, state, wavelength = row[0], row[1], number(row[2])
        elif any(row[:4]):
            msg = f"Table II: a continuation row with a line of its own: {row}"
            raise ValueError(msg)
        if row[5] not in REFERENCES:
            msg = f"Table II: unexpected reference {row[5]!r}"
            raise ValueError(msg)
        label = state.replace("^", "").replace("_", "")  # '^3P_1' -> '3P1'
        files[FILES[gas]].append(
            Point(
                f"{label}_ref{row[5]}",
                wavelength,
                number(row[4]),
                "",
                " | ".join(row),
            )
        )
    for name, points in files.items():
        if len(points) != 4:
            msg = f"Table II: {len(points)} values in {name}, not 4"
            raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
