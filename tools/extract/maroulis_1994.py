# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Maroulis (1994): mean static polarizability of ozone (Table IV).

Table IV (``tbl_590eb102ed40``, JSON cells only; checked against the image) compares
electric properties of O3 from experiment and theory. Only the mean polarizability row
is extracted: the 'Experiment' value (footnote c, 'Referenced by Mack and Muenter'),
the CASPT value of Andersson et al. (footnote d) and this paper's CCD+ST(CCD) value.
Footnote letters (``^c``) are stripped from the cells.

Third-order data: the gamma rows of Tables I (``tbl_1c2bf21aff33``, SCF by
basis set), II (``tbl_357dbeb4985d``, MBPT) and III (``tbl_9fd97a4e3ee8``, coupled
cluster), images checked, every cell agrees. One file per row (the six independent
components and the mean gamma-bar); the rows of Table I are the four basis sets, those
of Tables II and III together the methods SCF, MP2, MP3, DQ-MP4, CCD, CCD+S(CCD),
CCD+ST(CCD) (the increment columns D2, D3, DQR4, Delta-CCD, S(CCD), T(CCD) are not
taken; Table III's SCF column repeats Table II's). Cells printed "-285(1)" mean -285 x
10^1 (footnote a of Table I): the power is applied here, exactly, by moving the
decimal point. Table IV's gamma-bar "325x10^1" repeats Table III's CCD+ST(CCD) value.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Maroulis_1994_AccurateElectricDipole"
TABLE_IV = "tbl_590eb102ed40"
COLUMNS = {1: "mack_muenter", 3: "caspt_andersson", 4: "ccd_st"}
TABLE_I = "tbl_1c2bf21aff33"
TABLE_II = "tbl_357dbeb4985d"
TABLE_III = "tbl_9fd97a4e3ee8"
GAMMA_ROWS = {
    "\\gamma_{xxxx}": "xxxx",
    "\\gamma_{yyyy}": "yyyy",
    "\\gamma_{zzzz}": "zzzz",
    "\\gamma_{xxyy}": "xxyy",
    "\\gamma_{yyzz}": "yyzz",
    "\\gamma_{zzxx}": "zzxx",
    "\\bar{\\gamma}": "mean",
}
BASES = ("6s4p3d", "6s4p3d1f", "7s5p4d1f", "7s5p4d2f")
_POWER = re.compile(r"^(-?)(\d+)\((\d)\)$")


def _value(cell: str) -> str:
    """Return a cell's number, applying a printed power of ten: "-285(1)" -> "-2850"."""
    match = _POWER.match(cell.strip())
    if match is None:
        return number(cell)
    sign, digits, power = match.groups()
    return f"{sign}{digits}{'0' * int(power)}"


def _gamma() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    prefix = "O3/maroulis_1994__gamma_"
    for cells in table(DIRECTORY, TABLE_I)[1:]:
        name = GAMMA_ROWS.get(cells[0])
        if name is None:
            continue
        for basis, cell in zip(BASES, cells[1:], strict=True):
            files.setdefault(f"{prefix}scf_{name}.csv", []).append(
                Point(basis, "inf", _value(cell), "", " | ".join(cells))
            )
    for table_id, columns in (
        (TABLE_II, {1: "SCF", 5: "MP2", 6: "MP3", 7: "DQ-MP4"}),
        (TABLE_III, {5: "CCD", 6: "CCD+S(CCD)", 7: "CCD+ST(CCD)"}),
    ):
        for cells in table(DIRECTORY, table_id)[1:]:
            name = GAMMA_ROWS.get(cells[0])
            if name is None:
                continue
            for column, method in columns.items():
                files.setdefault(f"{prefix}corr_{name}.csv", []).append(
                    Point(method, "inf", _value(cells[column]), "", " | ".join(cells))
                )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the three mean-polarizability values, one file each, and the gammas."""
    (row,) = (c for c in json_table(DIRECTORY, TABLE_IV) if c[0] == "\\bar{\\alpha}")
    printed = " | ".join(row)
    return _gamma() | {
        f"O3/maroulis_1994__{name}.csv": [
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
