# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Zhao, Jin, Le, Jiang & Lin (2010): MO-ADK structure parameters, LDA+LB.

All tables were checked against their images (pp. 4-6, 10). Only the curated species
are taken: the rare gases, H2+ (an ion), the hydrogen molecule (printed "D_2/H_2",
recorded as H2), N2, O2, CO and CO2; F2, S2, NO, SO, C2H2, HF, HCl and HCN are not.

* Table I (``tbl_fba864ad07cc``, cells only, p. 4): rare-gas ionisation energies in
  atomic units: this work's exchange-only LDA + LB value, a second value marked ^b
  (Banerjee and Harbola, ref. 26) on the next row for He, Ne and Ar, and the
  experimental I_p marked ^a (Emsley, ref. 27). The marks are stripped.
* Table II (``tbl_b4bf810303ee``, cells only, p. 4): R (angstrom), the LDA + LB
  ionisation energy and the experimental vertical I_p (eV) of each molecule. R is
  written as the per-row column ``internuclear_distance`` of every C_l file of the
  molecule (Tables III and V give none).
* Table III (``tbl_7f5181e18a32``, cells only, p. 4): C_lm (columns C_0m ... C_6m) per
  molecule, this work's row first (the molecule label in its first cell), then the
  earlier values, each row marked in the last cell by its reference ([9] Tong et al.
  2002, [30] Awasthi et al. 2008, [17] Le et al. 2007). Blank cells are not values.
* Table IV (``tbl_2bc9b27f5a3b``, cells only, p. 6): binding energies (eV) of HOMO,
  HOMO-1 and HOMO-2: LDA + LB (this work), LB-alpha (^a ref. 29, ^d ref. 6) and the
  experimental vertical I_p (^b ref. 43, ^c ref. 44, ^e ref. 45); one file per molecule
  and column, rows labelled by the orbital. The marks are stripped (the note ^d,
  "Reference [6]", is lost from the library's notes and was read on the page).
* Table V (``tbl_825e661a911b``, cells only, p. 6): C_l of HOMO-1 and HOMO-2, each
  under its own row of column headers C_lm. The HOMO rows repeat Table III's first rows
  and are checked equal to them, not written again. The cells file drops the fifth
  column of the CO and HCl blocks (the image shows C_4m 0.05 for CO's 5 sigma, C_5m
  0.02 for CO's 1 pi, and HCl's fifth values): CO's 1 pi C_5m = 0.02 is added below
  from the image, guarded by the exact extracted row, and CO's 5 sigma C_4m is Table
  III's 0.05.
* Table VI (``tbl_56d1c04a1790``, CSV, p. 10): H2+ at five R. Only the 1s sigma_g row at
  R = 2.0 a.u. (the ground state at equilibrium) concerns a curated datum, and it is
  checked equal to Table III's H2+ row; nothing is written from it.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Zhao_2010_DeterminationStructureParameters"
TABLE_I = "tbl_fba864ad07cc"
TABLE_II = "tbl_b4bf810303ee"
TABLE_III = "tbl_7f5181e18a32"
TABLE_IV = "tbl_2bc9b27f5a3b"
TABLE_V = "tbl_825e661a911b"
TABLE_VI = "tbl_56d1c04a1790"
ATOMS = ("He", "Ne", "Ar", "Kr", "Xe")
SIGMA = "\N{GREEK SMALL LETTER SIGMA}"
# Table II and III molecule cells -> curated material.
MOLECULES = {
    "H_2^+": "H2+",
    "D_2/H_2": "H2",
    "N_2": "N2",
    "O_2": "O2",
    "CO": "CO",
    "CO_2": "CO2",
}
TABLE_III_LABELS = {
    "H_2^+(\\sigma_g)": "H_2^+",
    "D_2/H_2(\\sigma_g)": "D_2/H_2",
    "N_2(\\sigma_g)": "N_2",
    "O_2(\\pi_g)": "O_2",
    f"CO({SIGMA})": "CO",
    "CO_2(\\pi_g)": "CO_2",
}
TABLE_III_SKIPPED = {
    "F_2(\\pi_g)",
    "S_2(\\pi_g)",
    "NO(π)",
    "SO(π)",
    "C_2H_2(\\pi_u)",
    "HF(π)",
    "HCl(π)",
    "HCN(π)",
}
# Table III reference marks -> dataset name prefix.
REFERENCES = {"[9]": "tong", "[30]": "awasthi", "[17]": "le"}
TABLE_IV_MOLECULES = {"N_2": "N2", "O_2": "O2", "CO_2": "CO2", "CO": "CO"}
TABLE_IV_COLUMNS = {2: "ldalb", 3: "lbalpha", 4: "expt"}
_MARK = re.compile(r"\^\{?[a-f]\}?$")
_WAVE = re.compile(r"^C_\{(\d)m\}$")
_ORBITAL = re.compile(r"^(.+)\((HOMO(?:-\d)?)\)$")
# Table V, CO 1 pi (HOMO-1): the cells file keeps four of the five printed values; the
# image (p. 6) shows C_1m ... C_5m = 1.73, 0.14, 0.21, 0.02, 0.02.
CO_HOMO1_ROW = ["", "1\\pi(HOMO-1)", "1.73", "0.14", "0.21", "0.02"]
CO_HOMO1_LOST = ("5", "0.02")


def _strip(cell: str) -> str:
    """Return a cell's number without its footnote mark."""
    return number(_MARK.sub("", cell.strip()))


def _orbital(cell: str) -> str:
    r"""Return an orbital label as ASCII ("3\sigma_g(HOMO)" -> "3sigma_g (HOMO)")."""
    match = _ORBITAL.match(cell.strip())
    if match is None:
        msg = f"unexpected orbital {cell!r}"
        raise ValueError(msg)
    name = match.group(1).replace("\\sigma", "sigma").replace("\\pi", "pi")
    return f"{name} ({match.group(2)})"


def _table_i() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    atom = ""
    for cells in json_table(DIRECTORY, TABLE_I)[1:]:
        printed = " | ".join(cells)
        if cells[0].strip():
            atom = cells[0].strip()
            files[f"{atom}/zhao_2010__table1_ldalb_{atom}.csv"] = [
                Point("I_p", "", number(cells[1]), "", printed)
            ]
            files[f"{atom}/zhao_2010__table1_expt_{atom}.csv"] = [
                Point("I_p", "", _strip(cells[2]), "", printed)
            ]
        else:  # the ^b row (Banerjee and Harbola) under He, Ne and Ar
            if not cells[1].endswith("^b"):
                msg = f"unexpected Table I row {cells!r}"
                raise ValueError(msg)
            files[f"{atom}/zhao_2010__table1_banerjee_{atom}.csv"] = [
                Point("I_p", "", _strip(cells[1]), "", printed)
            ]
    if {name.split("/")[0] for name in files} != set(ATOMS):
        msg = "Table I lost an atom"
        raise ValueError(msg)
    return files


def _table_ii() -> tuple[dict[str, list[Point]], dict[str, str]]:
    files: dict[str, list[Point]] = {}
    distance: dict[str, str] = {}
    for cells in json_table(DIRECTORY, TABLE_II)[1:]:
        gas = MOLECULES.get(cells[0].strip())
        if gas is None:
            continue
        printed = " | ".join(cells)
        distance[gas] = number(cells[1])
        files[f"{gas}/zhao_2010__table2_ldalb_{gas}.csv"] = [
            Point("I_p", "", number(cells[2]), "", printed)
        ]
        files[f"{gas}/zhao_2010__table2_ip_{gas}.csv"] = [
            Point("I_p", "", number(cells[3]), "", printed)
        ]
    if set(distance) != set(MOLECULES.values()):
        msg = "Table II lost a molecule"
        raise ValueError(msg)
    return files, distance


def _waves(header: list[str]) -> dict[int, str]:
    """Return the columns of a C_lm header row and their l."""
    return {
        column: match.group(1)
        for column, cell in enumerate(header)
        if (match := _WAVE.match(cell.strip()))
    }


def _coefficients(
    cells: list[str], waves: dict[int, str], distance: str
) -> list[Point]:
    printed = " | ".join(cells)
    return [
        Point(
            wave,
            wave,
            number(cells[column]),
            "",
            printed,
            {"internuclear_distance": distance},
        )
        for column, wave in waves.items()
        if cells[column].strip()
    ]


def _table_iii(distance: dict[str, str]) -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = json_table(DIRECTORY, TABLE_III)
    waves = _waves(grid[0])
    gas = ""
    for cells in grid[1:]:
        label = cells[0].strip()
        if label in TABLE_III_SKIPPED:
            gas = ""
            continue
        if label:
            gas = MOLECULES[TABLE_III_LABELS[label]]
            name = gas
        elif not gas:  # a reference row under a skipped molecule
            continue
        else:
            name = f"{REFERENCES[cells[-1].strip()]}_{gas}"
        if label and cells[-1].strip():
            msg = f"this work's row carries a reference: {cells!r}"
            raise ValueError(msg)
        files[f"{gas}/zhao_2010__table3_{name}.csv"] = _coefficients(
            cells, waves, distance[gas]
        )
    return files


def _table_iv() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in json_table(DIRECTORY, TABLE_IV)[1:]:
        if cells[0].strip():
            gas = TABLE_IV_MOLECULES.get(cells[0].strip(), "")
        if not gas:
            continue
        printed = " | ".join(cells)
        orbital = _orbital(cells[1])
        for column, name in TABLE_IV_COLUMNS.items():
            if cells[column].strip():
                path = f"{gas}/zhao_2010__table4_{name}_{gas}.csv"
                files.setdefault(path, []).append(
                    Point(orbital, orbital, _strip(cells[column]), "", printed)
                )
    return files


def _table_v(
    distance: dict[str, str], table_iii: dict[str, list[Point]]
) -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    gas = ""
    waves: dict[int, str] = {}
    for cells in json_table(DIRECTORY, TABLE_V)[1:]:
        if cells[0].strip():
            gas = TABLE_IV_MOLECULES.get(cells[0].strip(), "")
        if not gas:
            continue
        if not cells[1].strip():  # a row of column headers
            waves = _waves(cells)
            continue
        orbital = _orbital(cells[1])
        points = _coefficients(cells, waves, distance[gas])
        if orbital.endswith("(HOMO)"):
            # Table III's row, in the columns the cells file keeps.
            homo = table_iii[f"{gas}/zhao_2010__table3_{gas}.csv"]
            full = {p.row: p.value for p in homo}
            if any(full.get(p.row) != p.value for p in points):
                msg = f"Table V's {gas} HOMO differs from Table III: {cells!r}"
                raise ValueError(msg)
            continue
        if gas == "CO" and orbital == "1pi (HOMO-1)":
            if cells != CO_HOMO1_ROW:
                msg = f"the CO 1 pi row changed in the library: {cells!r}"
                raise ValueError(msg)
            wave, value = CO_HOMO1_LOST
            printed = " | ".join(cells) + f" | {value} (from the image)"
            points = [
                Point(p.row, p.x, p.value, "", printed, p.extra) for p in points
            ] + [
                Point(
                    wave,
                    wave,
                    value,
                    "",
                    printed,
                    {"internuclear_distance": distance[gas]},
                )
            ]
        level = orbital.rsplit("(", 1)[1].rstrip(")").replace("HOMO-", "homo")
        files[f"{gas}/zhao_2010__table5_{level}_{gas}.csv"] = points
    return files


def _check_table_vi(table_iii: dict[str, list[Point]]) -> None:
    """Check Table VI's 1s sigma_g row at R = 2.0 a.u. against Table III's H2+."""
    grid = table(DIRECTORY, TABLE_VI)
    if grid[1][0] != f"1s{SIGMA}g" or grid[2][0] != "2.0":
        msg = f"unexpected Table VI layout {grid[:3]!r}"
        raise ValueError(msg)
    c0, c2 = grid[2][3].split()
    printed = {"0": c0, "2": c2, "4": grid[2][4]}
    ours = {p.row: p.value for p in table_iii["H2+/zhao_2010__table3_H2+.csv"]}
    if printed != ours:
        msg = f"Table VI's H2+ ground state {printed} differs from Table III {ours}"
        raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return the datasets of Tables I-V."""
    table_ii, distance = _table_ii()
    table_iii = _table_iii(distance)
    _check_table_vi(table_iii)
    return (
        _table_i() | table_ii | table_iii | _table_iv() | _table_v(distance, table_iii)
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
