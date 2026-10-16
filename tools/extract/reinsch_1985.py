# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Reinsch (1985): MC-SCF dynamic polarizabilities of He, Ne and CH4.

All tables are cells files only (images checked); omega is the photon energy in
hartree (Table I note a), the polarizabilities are in atomic units; rows are the
printed omega. HF, H2O and NH3 are not curated species and are skipped.

- Table I (``tbl_26de7c6c1f20``): He, this work SCF and MC-SCF, Chung, Starkschall and
  Gordon, experiment (Pade fit); one file per column (blank cells skipped).
- Table III (``tbl_e741f654ef3f``): Ne, this work in the standard (Sbs.) and extended
  (Ebs.) basis, experiment; the HF columns are skipped.
- Table V (``tbl_eecb6e4f851f``): CH4, this work and experiment (NH3 skipped). The
  footnote mark ``^{f}`` on the static value (``{15.616}^{f}``) and the braces are
  stripped.
- Table VI (``tbl_c3d0b718917c``, CSV): the Ne row, static SCF, MC-SCF, PNO-CI and CEPA
  (Werner and Meyer except this work's MC-SCF).

H2 (Table II) is not extracted here: the upper block holds tensor components at fixed
internuclear distances, and the lower block (the mean and vibrationally averaged
values) is missing from the library's cells; the record enters the latter inline from
the table image.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Reinsch_1985_CalculationDynamicPolarizabilities"
TABLE_I = "tbl_26de7c6c1f20"
TABLE_III = "tbl_e741f654ef3f"
TABLE_V = "tbl_eecb6e4f851f"
TABLE_VI = "tbl_c3d0b718917c"
HE = {1: "scf", 2: "mcscf", 3: "chung", 4: "starkschall_gordon", 5: "expt"}
NE = {1: "mcscf_standard", 2: "mcscf_extended", 3: "expt"}
CH4 = {5: "mcscf", 6: "expt"}
NE_STATIC = {
    1: "werner_meyer_scf",
    2: "mcscf_static",
    3: "werner_meyer_pno_ci",
    4: "werner_meyer_cepa",
}
_MARKED = re.compile(r"^\{(?P<value>[\d.]+)\}\^\{\w\}$")


def _value(cell: str) -> str:
    """Return a cell's number, without the braces of a footnote-marked value."""
    match = _MARKED.match(cell.strip())
    return number(match["value"] if match else cell)


def _columns(
    rows: list[list[str]], columns: dict[int, str], gas: str
) -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in rows:
        omega = cells[0].strip()
        for column, name in columns.items():
            if not cells[column].strip():
                continue
            path = f"{gas}/reinsch_1985__{name}_{gas}.csv"
            files.setdefault(path, []).append(
                Point(omega, omega, _value(cells[column]), "", " | ".join(cells))
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the He, Ne and CH4 columns."""
    files = _columns(json_table(DIRECTORY, TABLE_I)[2:], HE, "He")
    files |= _columns(json_table(DIRECTORY, TABLE_III)[3:], NE, "Ne")
    files |= _columns(json_table(DIRECTORY, TABLE_V)[2:], CH4, "CH4")
    (neon,) = (r for r in table(DIRECTORY, TABLE_VI) if r[0].strip() == "Ne")
    for column, name in NE_STATIC.items():
        path = f"Ne/reinsch_1985__{name}_Ne.csv"
        files[path] = [Point("0", "0", number(neon[column]), "", " | ".join(neon))]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
