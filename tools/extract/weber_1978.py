# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Weber, Milam & Smith (1978): gamma of glasses and crystals (Tables 1-3).

- Table 1 (``tbl_b9afd7663851``, p. 466; JSON only): glasses, with n_D, gamma in 10^-20
  m^2/W with its +- and the reference letter (a-f, listed under the table). Group
  headings ("Nd Laser Glasses", "Faraday Rotator Glass") span the row and are skipped.
  ED-2 has a second row with empty glass, type and n_D cells and its own gamma and
  reference (d); it is ED-2's (the table image shows the value under ED-2's). An
  asterisk after the +- ("21.2 ± 2.1*") marks the footnote "Measured relative to ED-2
  with gamma = 4.05 x 10^-20 m^2/W"; it is stripped and kept in ``printed``.
- Table 2 (``tbl_279cf82b90e4``, p. 466; JSON only): crystals, the same columns (its own
  reference list a-e). The "≈" printed before the n_D of La2Be2O5:Nd and CeF3 is
  extracted at the end of the structure-type cell (the image shows it before n_D);
  n_D is not extracted as a value, so only ``printed`` carries it.
- Table 3 (``tbl_b4032b6777cc``, p. 467): the nuclear parameters A0 and B0 (Heiman and
  Hellwarth) and n2(nuclear) and n2(total) in 10^-14 esu, and f_e. n2(nuclear) and
  n2(total) are extracted; A0, B0 and f_e have no quantity and stay in ``printed``.

All cells were compared with the table images: no misread cell.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Weber_1978_NonlinearRefractiveIndex"
TABLE_1 = "tbl_b9afd7663851"
TABLE_2 = "tbl_279cf82b90e4"
TABLE_3 = "tbl_b4032b6777cc"
PUMP_NM = "1064"

# Table 1 "Type" cell -> registered material id.
GLASSES = {
    "4000": "SiO2",
    "BeF_2": "BeF2",
    "SF-6": "SF-6",
    "BK-7": "BK7",
    "FK-5": "FK-5",
    "FK-51": "FK-51",
    "C-835": "C-835",
    "C-2828": "C-2828",
    "ED-2": "ED-2",
    "EV-1": "EV-1",
    "LHG-5": "LHG-5",
    "LHG-6": "LHG-6",
    "Q-88": "Q-88",
    "A86-82": "A86-82",
    "LG-812": "LG-812",
    "B-101": "B-101",
    "FR-4": "FR-4",
    "EY-1": "EY-1",
    "FR-5": "FR-5",
}
TABLE_1_GROUPS = ("Nd Laser Glasses", "Faraday Rotator Glass")

# Table 2 "Material" cell -> registered material id.
CRYSTALS = {
    "Y_{3}Al_{5}O_{12}": "YAG",
    "Al_{2}O_{3}": "Al2O3",
    "Al_{2}O_{3}:Cr": "Al2O3-Cr",
    "BeAl_{2}O_{4}": "BeAl2O4",
    "La_{2}Be_{2}O_{5}:Nd": "La2Be2O5",
    "KH_{2}PO_{4}": "KDP",
    "NaCl": "NaCl",
    "LiF": "LiF",
    "NaF": "NaF",
    "CaF_{2}": "CaF2",
    "SrF_{2}": "SrF2",
    "BaF_{2}": "BaF2",
    "CdF_{2}": "CdF2",
    "PbF_{2}": "PbF2",
    "MgF_{2}": "MgF2",
    "LaF_{3}": "LaF3",
    "CeF_{3}": "CeF3",
    "LiYF_{4}": "LiYF4",
}
TABLE_2_GROUPS = ("Oxide:", "Halide:")

# Table 3 "Type" cell -> registered material id.
TABLE_3_GLASSES = {
    "SF-6": "SF-6",
    "ED-4": "ED-4",
    "SiO_{2}": "SiO2",
    "BK-7": "BK7",
    "LHG-5": "LHG-5",
    "LHG-6": "LHG-6",
    "FK-5": "FK-5",
    "FK-51": "FK-51",
    "B-101": "B-101",
    "BeF_{2}": "BeF2",
}

_GAMMA = re.compile(r"^([\d.]+)\s*±\s*([\d.]+)(\*?)$")


def _gamma(cell: str) -> tuple[str, str]:
    """Return the value and the +- of a cell ``v ± u`` or ``v ± u*``, as printed."""
    match = _GAMMA.match(cell.strip())
    if match is None:
        msg = f"unexpected gamma cell {cell!r}"
        raise ValueError(msg)
    return number(match.group(1)), number(match.group(2))


def _table_1() -> dict[str, list[Point]]:
    """Return one file per glass and row of Table 1."""
    rows = json_table(DIRECTORY, TABLE_1)
    if rows[0][3] != r"\gamma (10^{-20} m^2/W)":
        msg = f"unexpected header {rows[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    material = ""
    for cells in rows[1:]:
        if cells[0] in TABLE_1_GROUPS:
            continue
        if cells[1]:
            material = GLASSES[cells[1]]
            name = f"{material}/weber_1978__gamma_{material}.csv"
        else:
            # ED-2's second row (reference d): empty glass, type and n_D cells.
            if material != "ED-2":
                msg = f"unexpected continuation row {cells!r}"
                raise ValueError(msg)
            name = f"{material}/weber_1978__gamma_{material}_{cells[4]}.csv"
        value, sigma = _gamma(cells[3])
        files[name] = [Point("1", PUMP_NM, value, sigma, " | ".join(cells))]
    if len(files) != len(GLASSES) + 1:
        msg = f"unexpected Table 1 rows: {sorted(files)!r}"
        raise ValueError(msg)
    return files


def _table_2() -> dict[str, list[Point]]:
    """Return one file per crystal of Table 2."""
    rows = json_table(DIRECTORY, TABLE_2)
    if rows[0][3] != r"\gamma (10^{-20} m^{2}/W)":
        msg = f"unexpected header {rows[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        if cells[0] in TABLE_2_GROUPS:
            continue
        material = CRYSTALS[cells[0]]
        value, sigma = _gamma(cells[3])
        files[f"{material}/weber_1978__gamma_{material}.csv"] = [
            Point("1", PUMP_NM, value, sigma, " | ".join(cells))
        ]
    if len(files) != len(CRYSTALS):
        msg = f"unexpected Table 2 rows: {sorted(files)!r}"
        raise ValueError(msg)
    return files


def _table_3() -> dict[str, list[Point]]:
    """Return n2(nuclear) and n2(total) of each glass of Table 3."""
    rows = table(DIRECTORY, TABLE_3)
    header = rows[0]
    if header[4] != "10^{14} n_{2}(nuclear)" or header[5] != "10^{14} n_{2}(total)":
        msg = f"unexpected header {header!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        material = TABLE_3_GLASSES[cells[1]]
        printed = " | ".join(cells)
        files[f"{material}/weber_1978__n2_nuclear_{material}.csv"] = [
            Point("1", "not stated", number(cells[4]), "", printed)
        ]
        files[f"{material}/weber_1978__n2_total_{material}.csv"] = [
            Point("1", PUMP_NM, number(cells[5]), "", printed)
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the data files of Tables 1-3."""
    return _table_1() | _table_2() | _table_3()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
