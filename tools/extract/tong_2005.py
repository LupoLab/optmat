# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Tong & Lin (2005): the barrier-suppression constant alpha, model potentials.

Both tables are CSV renderings, checked against their images (p. 2596).

* Table 2 (``tbl_87ba350a98ce``): per atom or ion (columns Rb, H, He, Ar+, Ne+, Ar, Ne)
  the fitted alpha of Eq. (2), the critical field F_b of Eq. (1) and I_p, both in
  atomic units. alpha (``alpha_<target>``) and I_p (``ip_<target>``) are written for
  the curated targets; F_b is derived from I_p (record's notes).
* Table 1 (``tbl_3240f8c18932``): the model-potential parameters Z_c, a_1 ... a_6 of
  Eq. (4), one file per curated target (``potential_<target>``), rows labelled by the
  parameter.

H and Rb are not curated and are skipped. "Ar^+" and "Ne^+" are the ions Ar+ and Ne+.
"""

from _library import Point, number, table, write

DIRECTORY = "Tong_2005_EmpiricalFormulaStatic"
TABLE_1 = "tbl_3240f8c18932"
TABLE_2 = "tbl_87ba350a98ce"
TARGETS = {"H": "H", "He": "He", "Ne": "Ne", "Ar": "Ar", "Ne+": "Ne+", "Ar+": "Ar+"}
TABLE_1_TARGETS = {"He": "He", "Ne": "Ne", "Ar": "Ar", "Ne^+": "Ne+", "Ar^+": "Ar+"}
SKIPPED = {"Rb"}
# Hydrogen's alpha and I_p are curated since atomic hydrogen became a material (TC-833);
# its potential row (Z_c = 1, every a_i = 0) is the bare Coulomb potential, not curated.
TABLE_1_SKIPPED = {"H", "Rb"}
ALPHA = "\N{GREEK SMALL LETTER ALPHA}"  # the row label of alpha
PARAMETERS = ("Z_c", "a_1", "a_2", "a_3", "a_4", "a_5", "a_6")


def _table_2() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = table(DIRECTORY, TABLE_2)
    rows = {cells[0].strip(): cells for cells in grid[1:]}
    if set(rows) != {ALPHA, "Fb", "Ip"}:
        msg = f"unexpected Table 2 rows {sorted(rows)}"
        raise ValueError(msg)
    for column, heading in enumerate(grid[0][1:], start=1):
        if heading in SKIPPED:
            continue
        target = TARGETS[heading]
        printed = f"{heading} | " + " | ".join(
            f"{name}: {rows[name][column]}" for name in (ALPHA, "Fb", "Ip")
        )
        files[f"{target}/tong_2005__alpha_{target}.csv"] = [
            Point("alpha", "", number(rows[ALPHA][column]), "", printed)
        ]
        files[f"{target}/tong_2005__ip_{target}.csv"] = [
            Point("I_p", "", number(rows["Ip"][column]), "", printed)
        ]
    return files


def _table_1() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = table(DIRECTORY, TABLE_1)
    if tuple(grid[0][1:]) != PARAMETERS:
        msg = f"unexpected Table 1 header {grid[0]!r}"
        raise ValueError(msg)
    for cells in grid[1:]:
        label = cells[0].strip()
        if label in TABLE_1_SKIPPED:
            continue
        target = TABLE_1_TARGETS[label]
        printed = " | ".join(cells)
        files[f"{target}/tong_2005__potential_{target}.csv"] = [
            Point(name, name, number(cell), "", printed)
            for name, cell in zip(PARAMETERS, cells[1:], strict=True)
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return alpha, I_p and the model-potential parameters of the curated targets."""
    return _table_2() | _table_1()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
