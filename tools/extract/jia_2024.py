# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Jia, Yuen, Jing, Zhou, Lin & Zhao (2024): N2 SDI against R (TC-957).

Every cell was checked against the table images (2026-10-09). Pages are those of the
accepted manuscript.

* Table I (p. 4, ``tbl_99caa8ce4558``, cells ``.json`` only, the I_p header spans
  the five distances): the vertical I_p (eV) of N2 to the eleven N2^2+ states at R =
  1.00 ... 1.20 angstrom (Pandey et al. 2014), with each state's limit. One file; row
  ``<index>@<R>``, x the state, the column ``internuclear_distance`` (angstrom).
* Table II (p. 5, ``tbl_d1566ca9a202``, cells ``.json`` only): per orbital (3sigma_g,
  1pi_u, 2sigma_u) and R, the LB-alpha binding energy, the CASSCF+MRCI ("C+M")
  ionisation energy of Bhattacharya et al. 2021, three C_l and the alpha of the
  LB-alpha potential. The two energies: one file each, row ``<orbital>@<R>``. The C_l:
  one file per orbital and R (the C_l reader takes one orbital at one geometry, rows
  l). The sub-header rows ("C_{0m} C_{2m} C_{4m}", "C_{1m} C_{3m} C_{5m}") give the l
  of the columns and are checked. alpha stays in the printed cells.
* Table III (p. 5, ``tbl_250cc6975ec7``): the N2+ X-A (along x) and X-B (along z)
  transition dipoles (au) against R (bohr), from Langhoff et al. 1987 and Langhoff and
  Bauschlicher 1988. One file per pair, row R, the column ``internuclear_distance``
  (bohr); the unit-vector suffixes are checked and dropped.

* Table IV (p. 6, ``tbl_34e26f76cb3c``): the transition dipoles (au) between N2^2+
  states m and n (Table I's indices) at R = 1.00 ... 1.20 angstrom, this paper's
  calculation. One file per pair, as Table III: row R, the column
  ``internuclear_distance`` (angstrom); each pair's unit vector (x-hat or z-hat) is
  checked to be the same at every R and dropped. Material N2++ (registered in TC-954).
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Jia_2024_ImprovedModelDissociative"
TABLE_1 = "tbl_99caa8ce4558"
TABLE_2 = "tbl_d1566ca9a202"
TABLE_3 = "tbl_250cc6975ec7"
TABLE_4 = "tbl_34e26f76cb3c"
# Table IV: the pairs (m, n) of its columns, in order, and each one's axis.
PAIRS_4 = (
    (1, 4, "x"),
    (1, 10, "z"),
    (2, 5, "x"),
    (2, 6, "z"),
    (3, 6, "x"),
    (4, 7, "x"),
    (4, 8, "x"),
    (4, 9, "z"),
    (8, 10, "z"),
    (9, 10, "x"),
)
ANGSTROM = "\N{LATIN CAPITAL LETTER A WITH RING ABOVE}"
DISTANCES = ("1.00", "1.05", "1.10", "1.15", "1.20")
# Table II: the orbital's cell -> (file name, the l of the three C columns).
ORBITALS = {
    "3\\sigma_g": ("3sigmag", (0, 2, 4)),
    "1\\pi_u": ("1piu", (1, 3, 5)),
    "2\\sigma_u": ("2sigmau", (1, 3, 5)),
}


def _table_1() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_1)
    distances = [f"R=1.00{ANGSTROM}"] + [f"{r}{ANGSTROM}" for r in DISTANCES[1:]]
    if (
        grid[0][:5] != ["Index", "State", "Config.", "Limit", "I_p(eV)"]
        or grid[1][4:] != distances
    ):
        msg = f"unexpected header of {TABLE_1}: {grid[:2]!r}"
        raise ValueError(msg)
    points: list[Point] = []
    for cells in grid[2:]:
        index, state = cells[0], cells[1]
        printed = " | ".join(cells)
        for r, cell in zip(DISTANCES, cells[4:], strict=True):
            points.append(
                Point(
                    f"{index}@{r}",
                    state,
                    number(cell),
                    "",
                    printed,
                    {"internuclear_distance": r},
                )
            )
    if len(points) != 11 * len(DISTANCES):
        msg = f"Table I gave {len(points)} values"
        raise ValueError(msg)
    return {"N2/jia_2024__table1_ip_N2.csv": points}


def _table_2() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_2)
    if grid[0] != [
        "Orbitals",
        "R(\\text{\\AA})",
        "",
        "I_p(\\text{eV})",
        "",
        "C_l",
        "",
        "\\alpha",
    ] or grid[1][1:7] != ["", "LB\\alpha", "C+M", "C_{0m}", "C_{2m}", "C_{4m}"]:
        msg = f"unexpected header of {TABLE_2}: {grid[:2]!r}"
        raise ValueError(msg)
    lbalpha: list[Point] = []
    mrci: list[Point] = []
    files: dict[str, list[Point]] = {}
    orbital = ""
    waves: tuple[int, ...] = ()
    for index, cells in enumerate(grid[1:], start=1):
        if not cells[1]:
            # A sub-header row: the l of the C_l columns of the block it opens. The
            # first one also names the orbital; the others leave it to the next row.
            heads = tuple(
                int(c.removeprefix("C_{").removesuffix("m}")) for c in cells[4:7]
            )
            opens = cells[0] or grid[index + 1][0]
            if heads != ORBITALS[opens][1]:
                msg = f"unexpected C_l header {cells!r}"
                raise ValueError(msg)
            waves = heads
            orbital = cells[0] or orbital
            continue
        if cells[0]:
            orbital = cells[0]
        name, expected = ORBITALS[orbital]
        if waves != expected:
            msg = f"the C_l of {orbital} under the header {waves!r}"
            raise ValueError(msg)
        r = cells[1]
        printed = " | ".join(cells)
        extra = {"internuclear_distance": number(r)}
        row = f"{name}@{r}"
        lbalpha.append(Point(row, name, number(cells[2]), "", printed, extra))
        mrci.append(Point(row, name, number(cells[3]), "", printed, extra))
        files[f"N2/jia_2024__table2_cl_{name}_{r.replace('.', '')}_N2.csv"] = [
            Point(str(l), str(l), number(c), "", printed, extra)
            for l, c in zip(waves, cells[4:7], strict=True)  # noqa: E741 - l is the partial wave
        ]
    if len(lbalpha) != len(ORBITALS) * len(DISTANCES):
        msg = f"Table II gave {len(lbalpha)} rows"
        raise ValueError(msg)
    return {
        "N2/jia_2024__table2_lbalpha_N2.csv": lbalpha,
        "N2/jia_2024__table2_cm_N2.csv": mrci,
    } | files


def _table_3() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_3)
    if rows[0] != [
        "R(a_0)\\backslash State 1 \N{RIGHTWARDS ARROW} 2",
        "X^2\\Sigma_g^+ \\to \\mathrm{A}^2\\Pi_u",
        "X^2\\Sigma_g^+ \\to \\mathrm{B}^2\\Sigma_u^+",
    ]:
        msg = f"unexpected header of {TABLE_3}: {rows[0]!r}"
        raise ValueError(msg)
    pairs = {"xa": (1, "\\hat{x}"), "xb": (2, "\\hat{z}")}
    files: dict[str, list[Point]] = {}
    for name, (column, axis) in pairs.items():
        points = []
        for cells in rows[1:]:
            cell = cells[column]
            if not cell.endswith(axis):
                msg = f"the X-{name[1].upper()} dipole is not along {axis}: {cells!r}"
                raise ValueError(msg)
            r = number(cells[0])
            points.append(
                Point(
                    r,
                    r,
                    number(cell.removesuffix(axis)),
                    "",
                    " | ".join(cells),
                    {"internuclear_distance": r},
                )
            )
        files[f"N2+/jia_2024__table3_{name}_N2+.csv"] = points
    return files


def _table_4() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_4)
    header = ["R(\\text{\\AA})\\backslash d_{\\text{m}\\rightarrow\\text{n}}"] + [
        f"d_{{{m}\\rightarrow {n}}}" for m, n, _ in PAIRS_4
    ]
    if rows[0] != header:
        msg = f"unexpected header of {TABLE_4}: {rows[0]!r}"
        raise ValueError(msg)
    if [cells[0] for cells in rows[1:]] != list(DISTANCES):
        msg = f"unexpected distances of {TABLE_4}: {rows[1:]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for column, (m, n, axis) in enumerate(PAIRS_4, start=1):
        suffix = f"\\hat{{{axis}}}"
        points = []
        for cells in rows[1:]:
            cell = cells[column]
            if not cell.endswith(suffix):
                msg = f"d({m}-{n}) is not along {axis}: {cells!r}"
                raise ValueError(msg)
            r = number(cells[0])
            points.append(
                Point(
                    r,
                    r,
                    number(cell.removesuffix(suffix)),
                    "",
                    " | ".join(cells),
                    {"internuclear_distance": r},
                )
            )
        files[f"N2++/jia_2024__table4_d{m}_{n}_N2++.csv"] = points
    return files


def extract() -> dict[str, list[Point]]:
    """Return the datasets of Tables I-IV."""
    return _table_1() | _table_2() | _table_3() | _table_4()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
