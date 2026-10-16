# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Zhao, Xu, Jin, Le & Lin (2011): MO-ADK C_lm of nonlinear molecules, HF.

CH4 and SF6 are taken, of the 17 molecules of Table 1. CH4 is on p. 7
(``tbl_6768a1c91416``, read from its ``.json`` cells: the table has merged headers and
no ``.csv``); its rows were checked against the table image and the page
(2026-09-28).

SF6, on the table's continuation (p. 8, ``tbl_4df3e935bc02``), is NOT taken from the
library's cells, which are not a rendering of the page: the page prints the table
rotated by 90 degrees, and the library's reading mixes CCl4's and SF6's coefficients,
loses SF6's ionisation energy and six of its columns, and prints every paired sign as
"+-". By an exception to the extraction rule (2026-09-28, TC-684) SF6's cells are those
read on the table image turned upright, written below as the page prints them
(``SF6_PAGE``). The script checks that the library's table is still the reading the
exception was made for (it fails if the table is re-extracted, and SF6 is then to be
taken from the new cells), that the cells the library does hold agree with the page
reading in header and magnitude (``CORROBORATED``), and that the three degenerate
orbitals have the same sum of |C_lm|^2 over m for each l, as orbitals that are rotations
of one another must (to the table's rounding).


CH4's block in the cells, which the script checks cell by cell:

* a row with the molecule, "1b_2(HOMO1)" and the ionisation energy, followed by the
  headers ``C_{1±1}`` ... ``C_{4±3}``; the next row has HOMO1's values; then a row
  with "1b_2(HOMO2)" and HOMO2's values; a row of HOMO3's headers ``C_{10}`` ...
  ``C_{4±2}``; a row with "1b_2(HOMO3)" and its values. The page prints the orbitals
  "1t_2(HOMO1)" ... and each label on its row of values; the library's "1b_2" is a
  misreading, and the record gives the printed label.

Structural fix: the library misreads two of CH4's paired signs, which the page prints
"\mp" (the upper sign minus) and the cells "\pm": HOMO1's C_{3±1} (page "-+0.13") and
HOMO3's C_{4±2} (page "-+0.11i"). Both were read on the page at 800 dpi and agree with
the PDF's own text layer (2026-09-28). ``MISREAD`` restores them, after checking that
the cell is the library's reading; ``printed`` keeps the page's cell. The other 18
cells of CH4 agree with the image.

A header ``C_{l±m}`` holds C_{l,+m} and C_{l,-m}; a value cell reads "\mp a" (C_{l,+m}
= -a, C_{l,-m} = +a: the upper sign for +m, the record's reading, not stated by the
paper), "\pm a" (+a, -a) or a plain number, the same for both. A value ending in "i"
is imaginary. Each file is one orbital and one signed m, the rows the partial waves l,
``value`` the real part and the extra column ``imaginary`` the imaginary part.
"""

import math
import re
from collections import defaultdict

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Zhao_2011_EffectOrbitalSymmetry"
TABLE = "tbl_6768a1c91416"
PLUS_MINUS = "\N{PLUS-MINUS SIGN}"
HEADER = ["Molecule", "Orbitals", "I_p (eV)", *[""] * 5, "C_{lm}", *[""] * 5]
FIRST_COEFFICIENT = 3
# The library's orbital cell -> the orbital as printed on the page, and its file name.
ORBITALS = {
    "1b_2(HOMO1)": ("1t2 HOMO1", "homo1"),
    "1b_2(HOMO2)": ("1t2 HOMO2", "homo2"),
    "1b_2(HOMO3)": ("1t2 HOMO3", "homo3"),
}
# The datasets the table must give: (orbital, m) -> the partial waves l.
EXPECTED = {
    ("homo1", 1): (1, 2, 3, 4),
    ("homo1", -1): (1, 2, 3, 4),
    ("homo1", 3): (3, 4),
    ("homo1", -3): (3, 4),
    ("homo2", 1): (1, 2, 3, 4),
    ("homo2", -1): (1, 2, 3, 4),
    ("homo2", 3): (3, 4),
    ("homo2", -3): (3, 4),
    ("homo3", 0): (1, 3),
    ("homo3", 2): (2, 4),
    ("homo3", -2): (2, 4),
}
# (orbital, header) -> (the library's cell, the cell on the page).
MISREAD = {
    ("homo1", f"C_{{3{PLUS_MINUS}1}}"): ("\\pm 0.13", "\\mp 0.13"),
    ("homo3", f"C_{{4{PLUS_MINUS}2}}"): ("\\pm 0.11i", "\\mp 0.11i"),
}
_HEADER = re.compile(rf"^C_\{{(?:(\d+){PLUS_MINUS}(\d+)|(\d)(0))\}}$")
_VALUE = re.compile(r"^(\\mp |\\pm )?(-?\d+\.\d+)(i?)$")

SF6_TABLE = "tbl_4df3e935bc02"
SF6_ENERGY = "15.69"
# SF6's cells as the page prints them (p. 8, Table 1 continued, the SF6 block), read
# on the library's table image turned upright (2026-09-28; TC-684): per orbital, the
# header and the value cell in the notation of the library's cells ("\mp" the upper
# sign minus). Not the library's cells: see the module docstring.
SF6_PAGE = {
    "homo1": (
        ("C_{4±1}", "\\mp 7.98"),
        ("C_{4±3}", "\\pm 3.01"),
        ("C_{6±1}", "\\mp 0.55"),
        ("C_{6±3}", "\\mp 0.87"),
        ("C_{6±5}", "\\pm 0.75"),
        ("C_{8±1}", "\\mp 0.37"),
        ("C_{8±3}", "\\pm 0.13"),
        ("C_{8±5}", "\\mp 0.12"),
        ("C_{8±7}", "\\pm 0.10"),
    ),
    "homo2": (
        ("C_{4±1}", "7.98i"),
        ("C_{4±3}", "3.01i"),
        ("C_{6±1}", "0.55i"),
        ("C_{6±3}", "-0.87i"),
        ("C_{6±5}", "-0.75i"),
        ("C_{8±1}", "0.37i"),
        ("C_{8±3}", "0.13i"),
        ("C_{8±5}", "0.12i"),
        ("C_{8±7}", "0.10i"),
    ),
    "homo3": (
        ("C_{4±4}", "\\pm 8.53i"),
        ("C_{6±4}", "\\mp 1.28i"),
        ("C_{8±4}", "\\pm 0.13i"),
        ("C_{8±8}", "\\pm 0.40i"),
        ("C_{10±4}", "\\mp 0.01i"),
        ("C_{10±8}", "\\mp 0.02i"),
    ),
}
# The page's cells that the library's garbled cells also hold, with the same header
# and magnitude (its paired signs are all "±"): 11 of the 24.
CORROBORATED = {
    ("homo1", "C_{4±1}"),
    ("homo1", "C_{4±3}"),
    ("homo1", "C_{6±1}"),
    ("homo1", "C_{8±5}"),
    ("homo1", "C_{8±7}"),
    ("homo2", "C_{4±1}"),
    ("homo2", "C_{6±1}"),
    ("homo2", "C_{8±7}"),
    ("homo3", "C_{4±4}"),
    ("homo3", "C_{8±4}"),
    ("homo3", "C_{8±8}"),
}
# The first cells of the library's SF6 row: the reading the exception was made for.
SF6_LIBRARY_ROW = ["SF_6", "2t_1(HOMO3)", "-"]
_LIBRARY_CELL = re.compile(rf"C_\{{(\d+{PLUS_MINUS}\d+)\}}{PLUS_MINUS}?-?(\d+\.\d+i?)")
# The largest difference between the orbitals' sums of |C_lm|^2 over m at one l that
# rounding to 0.01 allows: 2 |C| 0.005 per cell, summed.
ROUNDING = 0.05


def _signed(m: int) -> str:
    """Return m as the dataset id writes it: p1, m1, 0."""
    return "0" if m == 0 else f"{'p' if m > 0 else 'm'}{abs(m)}"


def _entries(header: str, cell: str) -> list[tuple[int, int, str, str]]:
    """Return (l, m, real part, imaginary part) of a header and its value cell."""
    heading = _HEADER.match(header)
    value = _VALUE.match(cell)
    if heading is None or value is None:
        msg = f"unexpected header or value: {header!r}, {cell!r}"
        raise ValueError(msg)
    paired = heading.group(1) is not None
    wave = int(heading.group(1) or heading.group(3))
    order = int(heading.group(2) or heading.group(4))
    sign, magnitude, imaginary = value.groups()
    magnitude = number(magnitude)
    if paired == (order == 0) or (sign and not paired):
        msg = f"a paired sign or m does not fit: {header!r}, {cell!r}"
        raise ValueError(msg)
    if not paired:
        signs = {0: magnitude}
    elif sign is None:
        signs = {order: magnitude, -order: magnitude}
    else:
        # The upper sign for +m (the record's reading); "\mp" is minus for +m.
        upper, lower = ("-", "") if sign == "\\mp " else ("", "-")
        signs = {order: upper + magnitude, -order: lower + magnitude}
    return [
        (wave, m, "0", text) if imaginary else (wave, m, text, "0")
        for m, text in signs.items()
    ]


def _block(
    headers: list[str], values: list[str], orbital: str
) -> dict[tuple[str, int], list[Point]]:
    """Return one orbital's points by (orbital, m), from its header and value rows."""
    points: dict[tuple[str, int], list[Point]] = {}
    for column in range(FIRST_COEFFICIENT, len(headers)):
        header, cell = headers[column], values[column]
        if not header and not cell:
            continue
        if (orbital, header) in MISREAD:
            library, page = MISREAD[orbital, header]
            if cell != library:
                msg = f"{orbital} {header} is {cell!r}, not the misread {library!r}"
                raise ValueError(msg)
            cell = page
        printed = f"{header} | {cell}"
        for wave, m, real, imaginary in _entries(header, cell):
            points.setdefault((orbital, m), []).append(
                Point(str(wave), str(wave), real, "", printed, {"imaginary": imaginary})
            )
    return points


def _magnitude(cell: str) -> str:
    r"""Return a value cell without its paired or minus sign: '\mp 7.98' -> '7.98'."""
    return cell.removeprefix("\\mp ").removeprefix("\\pm ").removeprefix("-")


def _check_sf6_library() -> None:
    """Check the library's SF6 cells: still the garbled reading, and corroborating."""
    grid = json_table(DIRECTORY, SF6_TABLE)
    starts = [index for index, cells in enumerate(grid) if cells[0] == "SF_6"]
    if len(starts) != 1 or grid[starts[0]][:3] != SF6_LIBRARY_ROW:
        msg = (
            f"the library's {SF6_TABLE} is no longer the reading of TC-684: take SF6 "
            "from its cells and retire SF6_PAGE"
        )
        raise ValueError(msg)
    held = {
        (match.group(1), match.group(2))
        for cells in grid[starts[0] :]
        for cell in cells[3:]
        for match in _LIBRARY_CELL.finditer(cell)
    }
    found = {
        (orbital, header)
        for orbital, cells in SF6_PAGE.items()
        for header, cell in cells
        if (header.removeprefix("C_{").removesuffix("}"), _magnitude(cell)) in held
    }
    if found != CORROBORATED:
        msg = f"the library's SF6 cells corroborate {sorted(found)}"
        raise ValueError(msg)


def _check_degenerate(points: dict[tuple[str, int], list[Point]]) -> None:
    """Check that the three orbitals have one sum of |C_lm|^2 over m per l."""
    sums: dict[int, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for (orbital, _m), rows in points.items():
        for point in rows:
            square = float(point.value) ** 2 + float(point.extra["imaginary"]) ** 2
            sums[int(point.row)][orbital] += square
    for wave, by_orbital in sums.items():
        if len(by_orbital) < len(SF6_PAGE):
            continue  # l = 10 is printed for HOMO3 only (0.001: below the rounding)
        values = list(by_orbital.values())
        if max(values) - min(values) > ROUNDING and not math.isclose(
            max(values), min(values), rel_tol=0.02
        ):
            msg = f"SF6's orbitals differ at l = {wave}: {dict(by_orbital)}"
            raise ValueError(msg)


def _sf6() -> dict[str, list[Point]]:
    """Return SF6's ionisation energy and the C_lm of its three HOMOs (TC-684)."""
    _check_sf6_library()
    points: dict[tuple[str, int], list[Point]] = {}
    for orbital, cells in SF6_PAGE.items():
        for header, cell in cells:
            for wave, m, real, imaginary in _entries(header, cell):
                points.setdefault((orbital, m), []).append(
                    Point(
                        str(wave),
                        str(wave),
                        real,
                        "",
                        f"{header} | {cell}",
                        {"imaginary": imaginary},
                    )
                )
    _check_degenerate(points)
    files = {
        "SF6/zhao_2011__table1_ip_SF6.csv": [
            Point("I_p", "", SF6_ENERGY, "", f"SF_6 | 1t_1g(HOMO1) | {SF6_ENERGY}")
        ]
    }
    for (orbital, m), rows in sorted(points.items()):
        files[f"SF6/zhao_2011__table1_SF6_{orbital}_m{_signed(m)}.csv"] = rows
    return files


def extract() -> dict[str, list[Point]]:
    """Return the ionisation energies and the C_lm of the HOMOs of CH4 and SF6."""
    return _ch4() | _sf6()


def _ch4() -> dict[str, list[Point]]:
    """Return CH4's ionisation energy and the C_lm of its three HOMOs."""
    grid = json_table(DIRECTORY, TABLE)
    if grid[0] != HEADER:
        msg = f"unexpected header of {TABLE}: {grid[0]!r}"
        raise ValueError(msg)
    starts = [index for index, cells in enumerate(grid) if cells[0] == "CH_4"]
    if len(starts) != 1:
        msg = f"CH_4 is on the rows {starts} of {TABLE}"
        raise ValueError(msg)
    start = starts[0]
    first, homo1, homo2, third, homo3 = grid[start : start + 5]
    labels = (first[:3], homo1[:3], homo2[:3], third[:3], homo3[:3])
    if labels != (
        ["CH_4", "1b_2(HOMO1)", "12.60"],
        ["", "", ""],
        ["", "1b_2(HOMO2)", ""],
        ["", "", ""],
        ["", "1b_2(HOMO3)", ""],
    ):
        msg = f"unexpected layout of CH4's rows: {labels!r}"
        raise ValueError(msg)
    if grid[start + 5][0] != "C_2H_6":
        msg = f"CH4's block does not end where expected: {grid[start + 5]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {
        "CH4/zhao_2011__table1_ip_CH4.csv": [
            Point("I_p", "", number(first[2]), "", " | ".join(first[:3]))
        ]
    }
    points = (
        _block(first, homo1, ORBITALS[first[1]][1])
        | _block(first, homo2, ORBITALS[homo2[1]][1])
        | _block(third, homo3, ORBITALS[homo3[1]][1])
    )
    found = {key: tuple(int(p.row) for p in rows) for key, rows in points.items()}
    if found != EXPECTED:
        msg = f"CH4's coefficients gave {found!r}"
        raise ValueError(msg)
    for (orbital, m), rows in points.items():
        name = f"table1_CH4_{orbital}_m{_signed(m)}"
        files[f"CH4/zhao_2011__{name}.csv"] = rows
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
