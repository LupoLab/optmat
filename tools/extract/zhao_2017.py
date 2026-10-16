# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Zhao, Li, Wang, Li & Zhou (2017): MO-ADK structure parameters, LB-alpha.

Only the curated species are taken, of the 123 molecules of the paper: H2+ (an ion),
H2, N2, O2, CO, CO2 and N2O. Every table is read from the library's cells
(``.json``), which exist for every fragment; the rows taken were checked against the
table images (2026-09-28).

* Table 1 (pp. 292-293), four page fragments (``tbl_3e9e6c6c77bf``,
  ``tbl_5c04cdea4543``, ``tbl_8aa77f541a54``, ``tbl_eeca443820ae``): R (angstrom), the
  LB-alpha ionisation energy, the "experimental vertical" I_p (eV) and the alpha of
  the LB-alpha potential. The two energies are written, one file each per molecule;
  R and alpha stay in the ``printed`` cells, and R is the per-row column
  ``internuclear_distance`` of the molecule's C_lm files. N2O has two distances, on
  two rows of the table ("1.128 (R_{NN})", "1.184 (R_{NO})"): the columns
  ``internuclear_distance_nn`` and ``internuclear_distance_no``.
* Table 2 (pp. 294-297), four page fragments (``tbl_5cfe79479e68``,
  ``tbl_a21cdf242a86``, ``tbl_8caed4c2ffc4``, ``tbl_3ac1855eb657``): C_lm of the
  HOMO, this work's row first (the molecule in its first cell), then the rows the
  table reprints, each marked by its reference in a last, unlabelled column ([30] Tong
  et al. 2002, [33] Kjeldsen and Madsen 2005, [53] Abu-samha and Madsen 2010, [34]
  Zhao et al. 2009, [60] Le et al. 2007). Blank cells are not values.

Structural fixes, all checked on the images:

* The rows of O2 run over a page break: the first fragment ends with O2's own row and
  its [30] row, and the second begins with two rows without a molecule ([53], [33]),
  which are O2's. The fragments are read in page order and the molecule is carried
  over.
* The cells of the first and last fragments of Table 2 have eleven columns, and the
  reference mark is in the last one, under the header C_9m (the page prints it to the
  right of C_9m); no row taken has an l = 9 value.
* The cells of the second and third fragments merge the columns C_5m ... C_9m and the
  reference mark into one cell, in which the l of a number cannot be told. The rows
  taken from them (O2's two reprinted rows, CO's three rows) print nothing beyond
  C_4m, and the script raises if such a cell holds more than a reference mark.

Table 3 (the z coordinates of the atoms of the polar molecules) has no quantity and
is quoted in the source record's notes.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Zhao_2017_AccurateStructureParameters"
ALPHA = "\N{GREEK SMALL LETTER ALPHA}"
SIGMA = "\N{GREEK SMALL LETTER SIGMA}"
PI = "\N{GREEK SMALL LETTER PI}"
RING = "\N{RING ABOVE}A"
ANGSTROM = "\N{LATIN CAPITAL LETTER A WITH RING ABOVE}"

# Table 1: fragment, its header and the cells of the molecules taken from it.
PLAIN_HEADER = ["Molecule", f"R/{RING}", f"LB{ALPHA}/eV", "Ip/eV", ALPHA]
TABLE_1 = (
    (
        "tbl_3e9e6c6c77bf",
        PLAIN_HEADER,
        {"H+2": "H2+", "H2": "H2", "N2": "N2", "O2": "O2"},
    ),
    ("tbl_5c04cdea4543", PLAIN_HEADER, {}),
    ("tbl_8aa77f541a54", PLAIN_HEADER, {"CO": "CO"}),
    (
        "tbl_eeca443820ae",
        ["Molecule", f"R/{ANGSTROM}", f"LB{ALPHA}/eV", "I_p/eV", ALPHA],
        {"CO_2": "CO2", "N_2O": "N2O"},
    ),
)
# Table 2: fragment (in page order), its layout and the cells of the molecules taken.
WIDE_HEADER = ["Molecule", *(f"C_{{{wave}m}}" for wave in range(10))]
MERGED_HEADER = [
    "Molecule",
    *(f"C{wave}m" for wave in range(5)),
    "C5m C6m C7m C8m C9m",
]
TABLE_2 = (
    (
        "tbl_5cfe79479e68",
        WIDE_HEADER,
        {
            f"H_2^+({SIGMA}_g)": "H2+",
            f"H_2({SIGMA}_g)": "H2",
            f"N_2({SIGMA}_g)": "N2",
            f"O_2({PI}_g)": "O2",
        },
    ),
    ("tbl_a21cdf242a86", MERGED_HEADER, {f"CO({SIGMA})": "CO"}),
    ("tbl_8caed4c2ffc4", MERGED_HEADER, {}),
    (
        "tbl_3ac1855eb657",
        WIDE_HEADER,
        {
            "\\text{CO}_{2}(\\pi_{g})": "CO2",
            "\\text{N}_{2}\\text{O}(\\pi)": "N2O",
        },
    ),
)
GASES = ("H2+", "H2", "N2", "O2", "CO", "CO2", "N2O")
# Table 2's reference marks -> dataset name.
REFERENCES = {
    "[30]": "tong",
    "[33]": "kjeldsen",
    "[53]": "abusamha",
    "[34]": "zhao2009",
    "[60]": "le",
}
# The datasets Table 2 must give: this work's row and the reprinted rows.
EXPECTED = {
    "H2+": ("", "tong"),
    "H2": ("", "tong", "kjeldsen"),
    "N2": ("", "tong", "abusamha", "kjeldsen"),
    "O2": ("", "tong", "abusamha", "kjeldsen"),
    "CO": ("", "tong", "kjeldsen"),
    "CO2": ("", "abusamha", "zhao2009", "le"),
    "N2O": ("",),
}
LOWEST_MERGED_WAVE = 5
_MARK = re.compile(r"^\[\d+\]$")
_BOND = re.compile(r"^(\d\.\d{3}) \(R_\{(NN|NO)\}\)$")


def _table_1() -> tuple[dict[str, list[Point]], dict[str, dict[str, str]]]:
    """Return the energy files and each molecule's distance column(s)."""
    files: dict[str, list[Point]] = {}
    distance: dict[str, dict[str, str]] = {}
    for table_id, header, wanted in TABLE_1:
        grid = json_table(DIRECTORY, table_id)
        if grid[0] != header:
            msg = f"unexpected header of {table_id}: {grid[0]!r}"
            raise ValueError(msg)
        for index, cells in enumerate(grid[1:], start=1):
            gas = wanted.get(cells[0])
            if gas is None:
                continue
            if gas in distance:
                msg = f"Table 1 has {gas} twice"
                raise ValueError(msg)
            printed = " | ".join(cells)
            if gas == "N2O":
                below = grid[index + 1]
                first, second = _BOND.match(cells[1]), _BOND.match(below[1])
                if (
                    first is None
                    or second is None
                    or (first.group(2), second.group(2)) != ("NN", "NO")
                    or any(cell for i, cell in enumerate(below) if i != 1)
                ):
                    msg = f"unexpected N2O rows of Table 1: {cells!r}, {below!r}"
                    raise ValueError(msg)
                printed += " || " + " | ".join(below)
                distance[gas] = {
                    "internuclear_distance_nn": number(first.group(1)),
                    "internuclear_distance_no": number(second.group(1)),
                }
            else:
                distance[gas] = {"internuclear_distance": number(cells[1])}
            files[f"{gas}/zhao_2017__table1_lbalpha_{gas}.csv"] = [
                Point("I_p", "", number(cells[2]), "", printed)
            ]
            files[f"{gas}/zhao_2017__table1_ip_{gas}.csv"] = [
                Point("I_p", "", number(cells[3]), "", printed)
            ]
    if set(distance) != set(GASES):
        msg = f"Table 1 gave {sorted(distance)}"
        raise ValueError(msg)
    return files, distance


def _reference(cells: list[str], header: list[str]) -> str:
    """Return a row's reference mark ("" for this work's row), checking the tail."""
    tail = cells[-1]
    if header is MERGED_HEADER and tail and not _MARK.match(tail):
        msg = f"the merged cell C_5m ... C_9m holds numbers: {cells!r}"
        raise ValueError(msg)
    if tail and not _MARK.match(tail):
        msg = f"a value under C_9m: {cells!r}"
        raise ValueError(msg)
    return tail


def _coefficients(
    cells: list[str], header: list[str], distance: dict[str, str]
) -> list[Point]:
    """Return a row's C_lm, one point per partial wave printed."""
    printed = " | ".join(cells)
    # Columns 1 ... 5 are l = 0 ... 4 in both layouts; the eleven-column layout
    # continues to l = 8 in column 9 (column 10 is the reference mark).
    last = LOWEST_MERGED_WAVE if header is MERGED_HEADER else len(header) - 2
    return [
        Point(str(wave), str(wave), number(cells[wave + 1]), "", printed, distance)
        for wave in range(last)
        if cells[wave + 1]
    ]


def _table_2(distance: dict[str, dict[str, str]]) -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    found: dict[str, list[str]] = {}
    gas = ""
    for table_id, header, wanted in TABLE_2:
        grid = json_table(DIRECTORY, table_id)
        if grid[0] != header:
            msg = f"unexpected header of {table_id}: {grid[0]!r}"
            raise ValueError(msg)
        for cells in grid[1:]:
            if cells[0]:
                gas = wanted.get(cells[0], "")
            if not gas:
                continue
            mark = _reference(cells, header)
            if bool(cells[0]) == bool(mark):
                msg = f"a row of {gas} is neither this work's nor a reprint: {cells!r}"
                raise ValueError(msg)
            name = REFERENCES[mark] if mark else ""
            found.setdefault(gas, []).append(name)
            stem = f"{name}_{gas}" if name else gas
            files[f"{gas}/zhao_2017__table2_{stem}.csv"] = _coefficients(
                cells, header, distance[gas]
            )
    if found != {gas: list(names) for gas, names in EXPECTED.items()}:
        msg = f"Table 2 gave the rows {found!r}"
        raise ValueError(msg)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the datasets of Tables 1 and 2."""
    table_1, distance = _table_1()
    return table_1 | _table_2(distance)


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
