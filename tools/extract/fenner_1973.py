# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Fenner, Hyatt, Kellam and Porto (1973): Raman cross sections of gases.

All four tables were checked against their images (the CSV cells match them):

- Table I (``tbl_0791ca776a3c``, p. 2), 'Pertinent properties of molecules studied': the
  vibrational shift (cm^-1) and the rotational constant (cm^-1) of each molecule,
  printed without a source. The curated molecules' shifts and constants are written.
- Table II (``tbl_422afc34de03``, p. 3): the ratio of each vibrational line's
  differential cross section to that of the N2 Q branch at 488.0 nm, this paper's column
  and four other workers' columns (adjusted by the authors to 488.0 nm with a nu_R^4
  law).
  ``...`` marks an empty cell. N2's own row (1.0 in every column, the reference) is not
  written.
- Table III (``tbl_81ab5bdc357b``, p. 3): the ratio of the rotational lines (the
  strongest Stokes line, and the sum of all lines) to the Q branch of the same molecule,
  and for N2 and O2 of the vibration-rotation (O and S branch) lines. The peak line is
  printed ``0(J)``: kept as printed in the row label. The superscripts a and b of H2 and
  CO2 (their reference is Q(1), resp. nu1) are dropped from the gas name; the record
  says what they mean.
- Table IV (``tbl_8d6a0ed08ffd``, p. 3): the absolute N2 Q-branch cross section, this
  paper's two determinations and four others' (adjusted to 488.0 nm). The column factor
  ``\times 10^{-31}`` is the dataset scale.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Fenner_1973_RamanCrossSection"
TABLE_I = "tbl_0791ca776a3c"
TABLE_II = "tbl_422afc34de03"
TABLE_III = "tbl_81ab5bdc357b"
TABLE_IV = "tbl_8d6a0ed08ffd"
EXCITATION = "488.0"  # nm, the Ar+ line of every value (standard air)
EMPTY = "..."

# Table I: (printed molecule, printed designation) -> (material, row label).
TABLE_I_ROWS = {
    ("N_{2}", ""): ("N2", "v=0-1"),
    ("O_{2}", ""): ("O2", "v=0-1"),
    ("H_{2}", ""): ("H2", "v=0-1"),
    ("CO", ""): ("CO", "v=0-1"),
    ("CO_{2}", "v_{1}"): ("CO2", "nu1"),
    ("CO_{2}", "2v_{2}"): ("CO2", "2nu2"),
    ("N_{2}O", "v_{1}"): ("N2O", "nu1"),
    ("N_{2}O", "v_{3}"): ("N2O", "nu3"),
    ("CH_{4}", "v_{1}"): ("CH4", "nu1"),
}

# Table II: printed gas -> (material, dataset suffix).
TABLE_II_ROWS = {
    "O_{2}": ("O2", "O2"),
    "H_{2} (sum)": ("H2", "H2_sum"),
    "H_{2} (Q(1))": ("H2", "H2_Q1"),
    "CO": ("CO", "CO"),
    "CO_{2}(v_{1})": ("CO2", "CO2_nu1"),
    "CO_{2}(2v_{2})": ("CO2", "CO2_2nu2"),
    "N_{2}O(v_{1})": ("N2O", "N2O_nu1"),
    "N_{2}O(v_{3})": ("N2O", "N2O_nu3"),
    "CH_{4}(v_{1})": ("CH4", "CH4_nu1"),
}
# Table II's columns after 'This paper', with the references of the table notes.
OTHERS = (
    "F & C (refs. 3, 4)",
    "M, H, B (ref. 2)",
    "S, C, W (ref. 1)",
    "P, G, L (ref. 5)",
)

# Table III: printed gas (superscript dropped) -> material.
TABLE_III_GASES = {
    "N_{2}": "N2",
    "O_{2}": "O2",
    "H_{2}": "H2",
    "CO": "CO",
    "CO_{2}": "CO2",
}

# Table IV: printed row -> row label; the first two are this paper's.
TABLE_IV_ROWS = {
    "This paper—calibrated spectrometer": "calibrated spectrometer",
    "This paper—comparison with v_{2} of C_{6}H_{6}^{a}": "against benzene 992 cm-1",
    "P, G, & Lb": "P, G, & L (ref. 5)",
    "F & Cc": "F & C (ref. 4)",
    "M, H, & Bd": "M, H, & B (ref. 2)",
    "S, C, & We": "S, C, & W (ref. 1)",
}
_ABSOLUTE = re.compile(
    r"^(?P<value>[\d.]+)(?: \\pm (?P<sigma>[\d.]+))? \\times 10\^\{-31\}"
    r"(?: cm2/ster-mol)?$"
)


def _table_i() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_I)[1:]:
        key = (cells[0], cells[2])
        if key not in TABLE_I_ROWS:
            continue
        gas, label = TABLE_I_ROWS[key]
        printed = " | ".join(cells)
        files.setdefault(f"{gas}/fenner_1973__table1_shift_{gas}.csv", []).append(
            Point(label, label, number(cells[1]), "", printed)
        )
        if cells[6]:
            files[f"{gas}/fenner_1973__table1_B_{gas}.csv"] = [
                Point("B", "B", number(cells[6]), "", printed)
            ]
    return files


def _table_ii() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_II)
    if rows[0][2] != "This paper":
        msg = f"Table II header {rows[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        if cells[0] not in TABLE_II_ROWS:
            continue
        gas, suffix = TABLE_II_ROWS[cells[0]]
        printed = " | ".join(cells)
        files[f"{gas}/fenner_1973__ratio_{suffix}.csv"] = [
            Point("this paper", EXCITATION, number(cells[2]), "", printed)
        ]
        others = [
            Point(label, EXCITATION, number(cell), "", printed)
            for label, cell in zip(OTHERS, cells[3:], strict=True)
            if cell != EMPTY
        ]
        if others:
            files[f"{gas}/fenner_1973__ratio_{suffix}_others.csv"] = others
    return files


def _table_iii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_III)[1:]:
        gas = TABLE_III_GASES[re.sub(r"\^\{[ab]\}$", "", cells[0])]
        printed = " | ".join(cells)
        if cells[1] == "vib-rot":
            if cells[2] != EMPTY:
                msg = f"Table III vib-rot peak cell {cells[2]!r}"
                raise ValueError(msg)
            files[f"{gas}/fenner_1973__vibrot_{gas}.csv"] = [
                Point("all vib-rot lines", EXCITATION, number(cells[3]), "", printed)
            ]
            continue
        files[f"{gas}/fenner_1973__rotational_{gas}.csv"] = [
            Point(f"peak line {cells[1]}", EXCITATION, number(cells[2]), "", printed),
            Point("all lines", EXCITATION, number(cells[3]), "", printed),
        ]
    return files


def _table_iv() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for index, (name, cell) in enumerate(table(DIRECTORY, TABLE_IV)):
        match = _ABSOLUTE.match(cell)
        if name not in TABLE_IV_ROWS or match is None:
            msg = f"Table IV row {name!r} | {cell!r}"
            raise ValueError(msg)
        point = Point(
            TABLE_IV_ROWS[name],
            EXCITATION,
            number(match["value"]),
            match["sigma"] or "",
            f"{name} | {cell}",
        )
        dataset = ("absolute_N2", "absolute_N2_benzene")[index] if index < 2 else None
        stem = f"N2/fenner_1973__{dataset or 'absolute_N2_others'}.csv"
        files.setdefault(stem, []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables I-IV as data files."""
    return _table_i() | _table_ii() | _table_iii() | _table_iv()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
