# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Majety & Scrinzi (2015): haCC static-field rates and ionic-state energies.

Every table was checked against its image; the cells agree digit for digit. Rates are
printed ``a \times 10^{-k}`` per cell (in the CSV "a", the multiplication sign, "10-k")
and are written as ``ae-k``; a cell printed ``<10^{-k}`` is an upper bound, written as
``1e-k`` for the record and excluded from its dataset (the record's ``excluded_rows``).

* Table 1 (``tbl_f4022cb09486``): He, rates (au) against F (au) for haCC(1), haCC(5),
  the literature columns [11], [17], [32] and "ADK theory"; one file per column
  (``table1_<column>_He``). Blank cells (no literature value) are skipped.
* Tables 2 and 3 (``tbl_111da5b2e5e8``, ``tbl_73f51612f76d``): H2, field parallel and
  perpendicular to the molecular axis, haCC(1), (2), (4), (5) and (Table 2 only) [18].
* Tables 4 and 5 (``tbl_c324d8092e8b``, ``tbl_80132d1acc1b``): vertical energies (eV)
  of the ionic states of N2 and CO, one file per column, the rows labelled by the ionic
  state; a dash (no value) is skipped. The state labels are extracted as "X2S+g",
  "A2Pu", ...; the images show X ^2Sigma_g^+, A ^2Pi_u, ...: the labels are written in
  ASCII from the images (``_STATES``), and a mismatch with the extraction stops the
  script.
* Table A1 (``tbl_8edc6fedf30a``): N2, haCC(5), one file per angle (0-90 deg).
* Table A2 (``tbl_4299943f807c``, 0-90 deg, and its continuation ``tbl_141c6c5af7dd``,
  105-180 deg, unlabelled in the library): CO, haCC(6), one file per angle.
  STRUCTURAL FIX: in ``tbl_4299943f807c`` the last row (F = 0.090) has the
  continuation's header merged into its cells (the 0 deg cell ends "10-3105°", ..., the
  75 deg cell "10-3180°"); the image shows that the row ends at the rate and that the
  angles 105°-180° head the next block. The script strips exactly the expected angle
  from each of those cells.
"""

import re

from _library import Point, table, write

DIRECTORY = "Majety_2015_StaticFieldIonization"
TABLE_1 = "tbl_f4022cb09486"
TABLE_2 = "tbl_111da5b2e5e8"
TABLE_3 = "tbl_73f51612f76d"
TABLE_4 = "tbl_c324d8092e8b"
TABLE_5 = "tbl_80132d1acc1b"
TABLE_A1 = "tbl_8edc6fedf30a"
TABLE_A2 = "tbl_4299943f807c"
TABLE_A2_CONTINUED = "tbl_141c6c5af7dd"

_RATE = re.compile(r"^(\d+\.\d+) \N{MULTIPLICATION SIGN} 10-(\d+)$")
_BOUND = re.compile(r"^<10-(\d+)$")

# Column headers as extracted, and the dataset stem each becomes.
TABLE_1_COLUMNS = {
    "haCC(1)": "hacc1",
    "haCC(5)": "hacc5",
    "[11]": "parker_2009",
    "[17]": "scrinzi_1999",
    "[32]": "themelis_1999",
    "ADK theory": "adk",
}
HACC_COLUMNS = {
    "haCC(1)": "hacc1",
    "haCC(2)": "hacc2",
    "haCC(4)": "hacc4",
    "haCC(5)": "hacc5",
    "[18]": "saenz_2000",
}
# The ionic states as extracted, and as the images print them (ASCII).
_STATES = {
    "X2S+g": "X2Sigma_g+",
    "A2Pu": "A2Pi_u",
    "B2S+u": "B2Sigma_u+",
    "C2S+u": "C2Sigma_u+",
    "D2Pg": "D2Pi_g",
    "X2S+": "X2Sigma+",
    "A2P": "A2Pi",
    "B2S+": "B2Sigma+",
    "D2Pu": "D2Pi_u",
    "32S+": "3 2Sigma+",
}
TABLE_4_COLUMNS = {
    "m-aug-cc-pvtz": "tz",
    "m-aug-cc-pvqz": "qz",
    "Svensson et al [35]": "svensson_1991",
}
TABLE_5_COLUMNS = {
    "m-aug-cc-pvtz": "tz",
    "Baltzer et al [36]": "baltzer_1994",
    "Svensson et al [35]": "svensson_1991",
}
ANGLES_A = ("0", "15", "30", "45", "60", "75", "90")
ANGLES_A2_CONTINUED = ("105", "120", "135", "150", "165", "180")
DEGREE = "\N{DEGREE SIGN}"


def _rate(cell: str) -> str:
    """Return a printed rate "a x 10-k" (or the bound "<10-k") as "ae-k" ("1e-k")."""
    text = cell.strip()
    match = _RATE.match(text)
    if match:
        return f"{match.group(1)}e-{match.group(2)}"
    match = _BOUND.match(text)
    if match:
        return f"1e-{match.group(1)}"
    msg = f"not a rate: {cell!r}"
    raise ValueError(msg)


def _columns(
    grid: list[list[str]], names: dict[str, str], material: str, stem: str
) -> dict[str, list[Point]]:
    """Return one file per named column of a rate table (field in column 0)."""
    header = grid[0]
    if header[0] != "F (a.u)" or header[1:] != list(names):
        msg = f"unexpected header {header!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for column, name in enumerate(names.values(), start=1):
        points = []
        for cells in grid[1:]:
            if not cells[column].strip():
                continue  # no value printed (a literature column)
            printed = " | ".join(cells)
            points.append(Point(cells[0], cells[0], _rate(cells[column]), "", printed))
        files[f"{material}/majety_2015__{stem}_{name}_{material}.csv"] = points
    return files


def _energies(
    table_id: str, names: dict[str, str], material: str, stem: str
) -> dict[str, list[Point]]:
    """Return one file per column of a table of vertical energies (rows: ion states)."""
    grid = table(DIRECTORY, table_id)
    if grid[0][0] != "State" or grid[0][1:] != list(names):
        msg = f"unexpected header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for column, name in enumerate(names.values(), start=1):
        points = []
        for cells in grid[1:]:
            if cells[column].strip() == "-":
                continue  # printed as a dash: no value
            state = _STATES[cells[0]]
            printed = " | ".join(cells)
            points.append(Point(state, state, cells[column].strip(), "", printed))
        files[f"{material}/majety_2015__{stem}_{name}_{material}.csv"] = points
    return files


def _angles(
    grid: list[list[str]], angles: tuple[str, ...], material: str, stem: str
) -> dict[str, list[Point]]:
    """Return one file per angle of an appendix table (field in column 0)."""
    files: dict[str, list[Point]] = {}
    for column, angle in enumerate(angles, start=1):
        points = []
        for cells in grid:
            printed = " | ".join(cells)
            points.append(Point(cells[0], cells[0], _rate(cells[column]), "", printed))
        files[f"{material}/majety_2015__{stem}_{material}_{angle}deg.csv"] = points
    return files


def _table_a1() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_A1)
    if grid[1] != ["(a.u)", *(f"{angle}{DEGREE}" for angle in ANGLES_A)]:
        msg = f"unexpected Table A1 header {grid[1]!r}"
        raise ValueError(msg)
    return _angles(grid[2:], ANGLES_A, "N2", "tableA1")


def _table_a2() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_A2)
    if grid[1] != ["(a.u)", *(f"{angle}{DEGREE}" for angle in ANGLES_A)]:
        msg = f"unexpected Table A2 header {grid[1]!r}"
        raise ValueError(msg)
    rows = [list(cells) for cells in grid[2:]]
    # Structural fix (see the module docstring): strip the continuation's header from
    # the last row, cell by cell, checking that exactly the expected angle is there.
    last = rows[-1]
    for column, angle in enumerate(ANGLES_A2_CONTINUED, start=1):
        suffix = f"{angle}{DEGREE}"
        if not last[column].endswith(suffix):
            msg = f"Table A2 last row: {last[column]!r} does not end with {suffix!r}"
            raise ValueError(msg)
        last[column] = last[column].removesuffix(suffix)
    continued = table(DIRECTORY, TABLE_A2_CONTINUED)
    if continued[0] != ["", *(f"{angle}{DEGREE}" for angle in ANGLES_A2_CONTINUED)]:
        msg = f"unexpected Table A2 continuation header {continued[0]!r}"
        raise ValueError(msg)
    if [cells[0] for cells in continued[1:]] != [cells[0] for cells in rows]:
        msg = "the Table A2 continuation's fields differ from the first part's"
        raise ValueError(msg)
    return _angles(rows, ANGLES_A, "CO", "tableA2") | _angles(
        continued[1:], ANGLES_A2_CONTINUED, "CO", "tableA2"
    )


def extract() -> dict[str, list[Point]]:
    """Return the rate and energy files of Tables 1-5, A1 and A2."""
    return (
        _columns(table(DIRECTORY, TABLE_1), TABLE_1_COLUMNS, "He", "table1")
        | _columns(table(DIRECTORY, TABLE_2), HACC_COLUMNS, "H2", "table2_parallel")
        | _columns(
            table(DIRECTORY, TABLE_3),
            {k: v for k, v in HACC_COLUMNS.items() if k != "[18]"},
            "H2",
            "table3_perpendicular",
        )
        | _energies(TABLE_4, TABLE_4_COLUMNS, "N2", "table4")
        | _energies(TABLE_5, TABLE_5_COLUMNS, "CO", "table5")
        | _table_a1()
        | _table_a2()
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
