# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hohm (2013): experimental static polarizabilities (au) of the curated gases.

The compilation's tables (cells files only; every in-scope row checked against its
image) print, per molecule, alpha_el "value ± error" with the reference numbers in
brackets and a dagger on the recommended value, and, for dielectric (DP) rows, alpha_el
+ alpha_vib (the measured static value, with its reference) and alpha_vib (with its
reference); alpha_el of a DP row is Hohm's difference of the two. The tables used:

- Table 1 (``tbl_3b3edcf5d4b1``, homonuclear diatomics): H2, D2, N2, O2;
- Table 2 (``tbl_56b3b73b65bf``, heteronuclear diatomics): CO;
- Table 3 (``tbl_0e38258e2324``, triatomics; the body was re-read by OCR, "all 229
  numbers agree with the text layer", and the in-scope rows agree with the image):
  CO2, N2O (the O3 cell is a range, "18.68 - 19.67", and is not extracted);
- Table 5 (``tbl_38210b618411``, penta-atomics): CH4;
- Table 7 (``tbl_a47f36ad8e7f``, hepta-atomics): SF6.

The captions are detached in the extraction (Table 1's caption sits in the text of
p. 2; Tables 6 and 7 carry the next table's caption as a footnote); the tables are
identified by their molecules and the text's table numbers. Each dataset is a set of
rows named by the reference of the value (an electronic value derived from a DP row is
named by that row's reference). The rows are addressed by index and checked against
their molecule and reference, so a change in the library's cells fails loudly.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Hohm_2013_ExperimentalStaticDipole"
_VALUE = re.compile(
    r"^(?P<value>\d+\.\d+) ± (?P<error>\d+\.\d+)(?P<dagger>†)?\s*(?P<refs>\[[^\]]+\])?$"
)
_CORRECTION = re.compile(r"^(?P<value>\d+\.\d+) (?P<refs>\[[^\]]+\])$")

TABLE_1 = "tbl_3b3edcf5d4b1"
TABLE_2 = "tbl_56b3b73b65bf"
TABLE_3 = "tbl_0e38258e2324"
TABLE_5 = "tbl_38210b618411"
TABLE_7 = "tbl_a47f36ad8e7f"

# Table 1 has the columns Molecule, alpha_el, Method, Remarks; the others Molecule,
# alpha_el, alpha_el + alpha_vib, alpha_vib, Method, Remarks.
EL, TOTAL, VIB = 1, 2, 3

# dataset file -> [(table, row index, molecule cell of the block, column, row id)]
SPEC: dict[str, list[tuple[str, int, str, int, str]]] = {
    "H2/hohm_2013__dp_H2.csv": [(TABLE_1, 1, "H_2", EL, "[43]")],
    "D2/hohm_2013__ri_D2.csv": [
        (TABLE_1, 2, "D_2", EL, "[77]"),
        (TABLE_1, 3, "D_2", EL, "[78]"),
    ],
    "D2/hohm_2013__dp_essen_D2.csv": [(TABLE_1, 4, "D_2", EL, "[79]")],
    "N2/hohm_2013__dp_N2.csv": [(TABLE_1, 10, "N_2", EL, "[44]")],
    "N2/hohm_2013__dosd_N2.csv": [(TABLE_1, 11, "N_2", EL, "[84]")],
    "O2/hohm_2013__dp_O2.csv": [(TABLE_1, 12, "O_2", EL, "[44,45]")],
    "O2/hohm_2013__dosd_O2.csv": [(TABLE_1, 13, "O_2", EL, "[85]")],
    "CO/hohm_2013__dosd_CO.csv": [(TABLE_2, 4, "CO", EL, "[97]")],
    "CO/hohm_2013__ri_mean_CO.csv": [(TABLE_2, 5, "CO", EL, "[98-102]")],
    "CO/hohm_2013__dp_total_CO.csv": [
        (TABLE_2, 6, "CO", TOTAL, "[103]"),
        (TABLE_2, 7, "CO", TOTAL, "[104]"),
        (TABLE_2, 8, "CO", TOTAL, "[105]"),
        (TABLE_2, 9, "CO", TOTAL, "[46]"),
        (TABLE_2, 10, "CO", TOTAL, "[106]"),
    ],
    "CO/hohm_2013__dp_electronic_CO.csv": [
        (TABLE_2, 6, "CO", EL, "[103]-[16]"),
        (TABLE_2, 7, "CO", EL, "[104]-[16]"),
        (TABLE_2, 8, "CO", EL, "[105]-[16]"),
        (TABLE_2, 9, "CO", EL, "[46]-[16]"),
        (TABLE_2, 10, "CO", EL, "[106]-[16]"),
    ],
    # The five CO rows print the same alpha_vib, 0.12 [16]; it is taken once.
    "CO/hohm_2013__alpha_vib_CO.csv": [(TABLE_2, 6, "CO", VIB, "[16]")],
    "CO2/hohm_2013__dosd_CO2.csv": [(TABLE_3, 15, "CO_2", EL, "[97]")],
    # Its reference list continues on the next row ("Ref. cont. [...]"), added to the
    # printed cells by extract().
    "CO2/hohm_2013__ri_mean_CO2.csv": [(TABLE_3, 16, "CO_2", EL, "[115,123-126]")],
    "CO2/hohm_2013__dp_total_schmidt_CO2.csv": [(TABLE_3, 18, "CO_2", TOTAL, "[44]")],
    # Rows 19 and 20 print the same mean, 19.67 ± 0.10 [130-132]; it is taken once.
    "CO2/hohm_2013__dp_total_mean_CO2.csv": [(TABLE_3, 19, "CO_2", TOTAL, "[130-132]")],
    "CO2/hohm_2013__dp_electronic_CO2.csv": [
        (TABLE_3, 18, "CO_2", EL, "[44]-[129]"),
        (TABLE_3, 19, "CO_2", EL, "[130-132]-[129]"),
        (TABLE_3, 20, "CO_2", EL, "[130-132]-[16]"),
    ],
    "CO2/hohm_2013__alpha_vib_thomas_CO2.csv": [(TABLE_3, 18, "CO_2", VIB, "[129]")],
    "CO2/hohm_2013__alpha_vib_bishop_CO2.csv": [(TABLE_3, 20, "CO_2", VIB, "[16]")],
    "N2O/hohm_2013__dosd_N2O.csv": [(TABLE_3, 21, "N_2O", EL, "[84]")],
    "N2O/hohm_2013__ri_hohm_N2O.csv": [(TABLE_3, 22, "N_2O", EL, "[125]")],
    "N2O/hohm_2013__ri_cuthbertson_N2O.csv": [(TABLE_3, 23, "N_2O", EL, "[50]")],
    "N2O/hohm_2013__ri_koch_N2O.csv": [(TABLE_3, 24, "N_2O", EL, "[107]")],
    "N2O/hohm_2013__dp_total_N2O.csv": [(TABLE_3, 25, "N_2O", TOTAL, "[133]")],
    "N2O/hohm_2013__dp_electronic_N2O.csv": [(TABLE_3, 25, "N_2O", EL, "[133]-[16]")],
    "N2O/hohm_2013__alpha_vib_N2O.csv": [(TABLE_3, 25, "N_2O", VIB, "[16]")],
    "CH4/hohm_2013__dosd_CH4.csv": [(TABLE_5, 1, "CH_4", EL, "[169]")],
    "CH4/hohm_2013__ri_mean_CH4.csv": [
        (TABLE_5, 2, "CH_4", EL, "[100,115,124,146,170-172]")
    ],
    "CH4/hohm_2013__dp_total_mean_CH4.csv": [
        (TABLE_5, 3, "CH_4", TOTAL, "[44,173,14,132,105]")
    ],
    "CH4/hohm_2013__dp_electronic_CH4.csv": [
        (TABLE_5, 3, "CH_4", EL, "[44,173,14,132,105]-[16]")
    ],
    "CH4/hohm_2013__alpha_vib_CH4.csv": [(TABLE_5, 3, "CH_4", VIB, "[16]")],
    "SF6/hohm_2013__dosd_SF6.csv": [(TABLE_7, 12, "SF_6", EL, "[227]")],
    "SF6/hohm_2013__ri_watson_SF6.csv": [(TABLE_7, 13, "SF_6", EL, "[115]")],
    "SF6/hohm_2013__ri_burns_SF6.csv": [(TABLE_7, 14, "SF_6", EL, "[172]")],
    # Rows 15 and 16 print the same mean, 44.21 ± 0.50 [228,105,192,132,229].
    "SF6/hohm_2013__dp_total_mean_SF6.csv": [
        (TABLE_7, 15, "SF_6", TOTAL, "[228,105,192,132,229]")
    ],
    "SF6/hohm_2013__dp_electronic_SF6.csv": [
        (TABLE_7, 15, "SF_6", EL, "[228,105,192,132,229]-[230]"),
        (TABLE_7, 16, "SF_6", EL, "[228,105,192,132,229]-[155]"),
    ],
    "SF6/hohm_2013__alpha_vib_pack_SF6.csv": [(TABLE_7, 15, "SF_6", VIB, "[230]")],
    "SF6/hohm_2013__alpha_vib_bulanin_SF6.csv": [(TABLE_7, 16, "SF_6", VIB, "[155]")],
}


def molecule(rows: list[list[str]], index: int) -> str:
    """Return the molecule of a row: its own first cell, or the block's above it."""
    while not rows[index][0]:
        index -= 1
    return rows[index][0]


def point(rows: list[list[str]], index: int, column: int, row_id: str) -> Point:
    """Return one value of a row, checked against the reference its row id names."""
    cells = rows[index]
    cell = cells[column].strip()
    printed = " | ".join(cells)
    if cells[1].startswith("17.50 ± 0.02†"):
        printed += " || " + " | ".join(rows[index + 1])  # the continued reference list
    if column == VIB:
        match = _CORRECTION.match(cell)
        if match is None or match["refs"] != row_id:
            msg = f"row {index}: alpha_vib cell {cell!r} is not {row_id}"
            raise ValueError(msg)
        return Point(row_id, "inf", number(match["value"]), "", printed)
    match = _VALUE.match(cell)
    if match is None:
        msg = f"row {index}: unreadable cell {cell!r}"
        raise ValueError(msg)
    # A value cell carries its reference; a derived alpha_el (DP rows) carries none and
    # is named by the references of its row's total and alpha_vib.
    expected = None if column == EL and "]-[" in row_id else row_id
    if match["refs"] != expected:
        msg = f"row {index}: cell {cell!r} does not carry {row_id}"
        raise ValueError(msg)
    return Point(row_id, "inf", number(match["value"]), number(match["error"]), printed)


def extract() -> dict[str, list[Point]]:
    """Return one file per in-scope gas and kind of value."""
    tables = {
        table_id: json_table(DIRECTORY, table_id)
        for table_id in (TABLE_1, TABLE_2, TABLE_3, TABLE_5, TABLE_7)
    }
    files: dict[str, list[Point]] = {}
    for path, entries in SPEC.items():
        points = []
        for table_id, index, gas, column, row_id in entries:
            rows = tables[table_id]
            if molecule(rows, index) != gas:
                msg = f"{table_id} row {index} is {molecule(rows, index)!r}, not {gas}"
                raise ValueError(msg)
            points.append(point(rows, index, column, row_id))
        files[path] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
