# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Franta et al. (2013), a-Si:H: the band gaps of the six dispersion models.

Table 2 (``tbl_d0ad1c1043c4``, cells only) lists, for the six models fitted to the same
data (PJDOS-EX, PDOS, PJDOS-TL, PJDOS-CC, PJDOS-CL, PJDOS-2EX), the fitted
parameters. The only row that is a value of a vocabulary quantity is E_g, the band
gap (``band_gap``, a model parameter): the authors' most reliable model, PJDOS-2EX
(p. 240), is one file, the other five another, rows labelled by the model.

The fitted parameters themselves are the record's ``[[formulas]]`` (one per model, and
one for the phonon terms of Table 3, ``tbl_4125a753d140``), transcribed from the
images. This script checks each recorded coefficient against the table cells, so that
the record and the library cannot drift apart:

- a cell "x^a" is the value x, fixed in the fit (footnote a); the record lists it
  under ``fixed``;
- a cell "-" is a parameter the model does not have;
- structural fix: Table 3's cells merge the two columns under "v_p (cm^-1)", "Ref.
  [21]" and "This work", into one cell ("2000 1996.6"), and read the footnote mark of
  a fixed value as a trailing "a" ("2120a"); the image shows two columns and a
  superscript a. The script splits them.

All cells agree with the images (checked 2026-09-26).
"""

import tomllib

from _json_table import json_table
from _library import DATA, Point, number, write

DIRECTORY = "Franta_2013_ApplicationSumRule"
TABLE_2 = "tbl_d0ad1c1043c4"
TABLE_3 = "tbl_4125a753d140"
PREFIX = "a-Si_H/franta_2013b__"
MODELS = ("PJDOS-EX", "PDOS", "PJDOS-TL", "PJDOS-CC", "PJDOS-CL", "PJDOS-2EX")
RECOMMENDED = "PJDOS-2EX"

# Table 2 row label -> coefficient name in the record (the exciton rows are named
# E_c, B_c for the one-Lorentzian models and E_c1, B_c1 for PDOS and PJDOS-2EX, as the
# row label "E_c|E_c1" says).
ROWS = {
    "N_a (eV^2)": "N_a",
    "f_H": "f_H",
    "n_K": "n_K",
    "E_K (eV)": "E_K",
    "n_L": "n_L",
    "E_L (eV)": "E_L",
    "E_g (eV)": "E_g",
    "E_p (eV)": "E_p",
    "E_h (eV)": "E_h",
    "A_{c1}": "A_c1",
    "E_c|E_{c1} (eV)": "E_c",
    "B_c|B_{c1} (eV)": "B_c",
    "A_{c2}": "A_c2",
    "E_{c2} (eV)": "E_c2",
    "B_{c2} (eV)": "B_c2",
    "\\alpha_x": "alpha_x",
    "E_x (eV)": "E_x",
    "\\alpha_u": "alpha_u",
    "E_u (eV)": "E_u",
    "\\alpha_{loc} \\times 10^6": "alpha_loc_1e6",
    "E_{loc} (eV)": "E_loc",
    "B_{loc} (eV)": "B_loc",
}
TWO_EXCITON = {"PDOS", "PJDOS-2EX"}


def _record() -> dict[str, dict[str, object]]:
    """Return the record's formulas by id."""
    path = DATA / "sources" / "franta_2013b.toml"
    with path.open("rb") as handle:
        return {entry["id"]: entry for entry in tomllib.load(handle)["formulas"]}


def _split(cell: str) -> tuple[str, bool]:
    """Return a cell's number as printed and whether it is marked fixed."""
    for mark in ("^a", "a"):
        if cell.endswith(mark):
            return number(cell.removesuffix(mark)), True
    return number(cell), False


def _check(
    formula: dict[str, object], expected: dict[str, str], fixed: set[str]
) -> None:
    """Check a recorded formula's coefficients and fixed list against the cells."""
    recorded = set(formula.get("fixed", []))  # pyright: ignore[reportArgumentType]
    if formula["coefficients"] != expected or recorded != fixed:
        msg = f"{formula['id']}: the record differs from the table, {expected}, {fixed}"
        raise ValueError(msg)


def _table_2() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_2)
    if tuple(rows[0][1:7]) != MODELS:
        msg = f"unexpected Table 2 header {rows[0]!r}"
        raise ValueError(msg)
    by_label = {cells[0]: cells for cells in rows if cells[0] in ROWS}
    if set(by_label) != set(ROWS):
        msg = f"Table 2 lacks the rows {set(ROWS) - set(by_label)}"
        raise ValueError(msg)
    formulas = _record()
    files: dict[str, list[Point]] = {}
    for column, model in enumerate(MODELS, start=1):
        expected: dict[str, str] = {}
        fixed: set[str] = set()
        for label, name in ROWS.items():
            cell = by_label[label][column]
            if cell == "-":
                continue
            two = model in TWO_EXCITON and name in {"E_c", "B_c"}
            key = f"{name}1" if two else name
            value, is_fixed = _split(cell)
            expected[key] = value
            if is_fixed:
                fixed.add(key)
        _check(formulas[f"udm_{model.lower().replace('-', '_')}"], expected, fixed)
        cells = by_label["E_g (eV)"]
        printed = f"{cells[0]} | {model} | {cells[column]}"
        point = Point(model, "", expected["E_g"], "", printed)
        name = "band_gap_pjdos_2ex" if model == RECOMMENDED else "band_gap_other_models"
        files.setdefault(f"{PREFIX}{name}.csv", []).append(point)
    return files


def _check_table_3() -> None:
    rows = json_table(DIRECTORY, TABLE_3)
    if rows[1][1] != "Ref. [21] This work" or len(rows) != 10:
        msg = f"unexpected Table 3 layout {rows[:2]!r}"
        raise ValueError(msg)
    expected: dict[str, str] = {}
    fixed: set[str] = set()
    for cells in rows[2:]:
        peak = cells[0]
        # Structural fix: the two position columns are one cell (see the docstring).
        reference, this_work = cells[1].split()
        expected[f"nu_{peak}_ref21"] = number(reference)
        expected[f"nu_{peak}"], is_fixed = _split(this_work)
        expected[f"beta_{peak}"] = number(cells[3])
        expected[f"alpha_{peak}"] = number(cells[4])
        if is_fixed:
            fixed.add(f"nu_{peak}")
    _check(_record()["phonons"], expected, fixed)


def extract() -> dict[str, list[Point]]:
    """Return the band gaps of Table 2, after checking the recorded formulas."""
    _check_table_3()
    return _table_2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
