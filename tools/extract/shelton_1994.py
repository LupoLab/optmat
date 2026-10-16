# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton & Rice (1994), Chem. Rev. 94, 3: the gamma tables for curated gases.

Every table is read from the library's JSON cells (``_json_table``); each was compared
with its image (``tables/<id>.jpg``). Each extracted cell carries footnote marks
(``#``, ``*``, ``@``, letters) glued to the number, printed as ``1870*,w``,
``(682.7)^{\#,c}`` or ``3257 ± 16#,x``; the script splits the number (and the ``±``
uncertainty) from the marks and checks the marks against the ones the image shows, so
that a change in the library's extraction cannot silently move a value to another
footnote (and so to another original source). Parentheses (vibrationally averaged
calculations, Table IX) are dropped from the number and noted in the record.

One data file per dataset, ``<gas>/shelton_1994__<id>.csv``. The row id
names the method, reference or process of the value; the wavelength is the laser line
as printed (nm in air, or the printed vacuum wavenumber of Table IX), ``inf`` for a
static value. Only curated species are taken (the tables also hold H, ions, HF, HCl,
NO, NH3, H2O, CF4, the halomethanes and hydrocarbons).

Structural notes from the images:

- Table V (``tbl_38363e48379b``) prints two side-by-side method lists; "CCSD+ST(CCD)"
  wraps onto a second line, beside which the SCF/IVa row (70.8, ref 176) sits (JSON row
  6, left half empty): the rows are matched as printed. The VPT cell holds two values,
  "8, 65" (ref 179): split into rows ``VPT a`` and ``VPT b``.
- Table VI (``tbl_678f6d31104f``) repeats gamma_calc = 44.211 on the three 632.8 nm dc
  Kerr rows: the calculated dc Kerr dataset keeps one 632.8 nm row (the first) and the
  514.5 nm row.
- Table XII (``tbl_f5d00dad6453``): part a (beta, not extracted) is misaligned in the
  JSON (the SO2, CH3OH, (CH3)2O rows are shifted), but the gamma rows of CH4, SF6 and
  CO2 (JSON rows 31-32, 37-38, 43-44) match the image cell for cell.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Shelton_1994_MeasurementsCalculationsHyperpolarizabilities"
T4 = "tbl_fbf8e009a36f"
T5 = "tbl_38363e48379b"
T6 = "tbl_678f6d31104f"
T8 = "tbl_a56d72648e32"
T9 = "tbl_a4a5bfa3d452"
T11 = "tbl_c6b7f0da4741"
T12 = "tbl_f5d00dad6453"

_CELL = re.compile(r"^\(?\s*(-?[\d.]+)\s*\)?\s*(?:±\s*([\d.]+))?\s*(.*)$")


def split(text: str) -> tuple[str, str, str]:
    """Return a cell's number, its ``±`` uncertainty (or "") and its footnote marks."""
    cleaned = text.replace("\\pm", "±").replace("\N{MINUS SIGN}", "-").strip()
    match = _CELL.match(cleaned)
    if match is None:
        msg = f"not a value cell: {text!r}"
        raise ValueError(msg)
    value, sigma, rest = match.groups()
    marks = re.sub(r"[\^{}\\\s]", "", rest)
    return number(value), number(sigma) if sigma else "", marks


# (dataset file stem, gas, table, [(json row, json column, row id, wavelength, marks)])
Cell = tuple[int, int, str, str, str]
SPECS: list[tuple[str, str, str, list[Cell]]] = []


def _add(stem: str, gas: str, table: str, cells: list[Cell]) -> None:
    SPECS.append((stem, gas, table, cells))


# Table IV: static gamma_0 of the atoms, correlated (col 3) and SCF (col 5); the
# footnote letters stand in their own "ref" columns (4 and 6), checked in extract().
T4_REFS = {
    3: ("d#", "e#"),
    22: ("k#", "e#"),
    26: ("n#", "n#"),
    28: ("n#", "n#"),
    29: ("n", "n"),
}
for _gas, _row, _method in [
    ("He", 3, "CI-Hylleraas"),
    ("Ne", 22, "CCSD(T)"),
    ("Ar", 26, "CCSD(T)"),
    ("Kr", 28, "CCSD(T)"),
    ("Xe", 29, "CCSD(T)"),
]:
    _add(
        f"table4_{_gas}",
        _gas,
        T4,
        [
            (_row, 3, f"{_method} (note {T4_REFS[_row][0]})", "inf", ""),
            (_row, 5, f"SCF (note {T4_REFS[_row][1]})", "inf", ""),
        ],
    )

# Table V: Ne static gamma_0 by method and basis set (left half cols 0-3, right 4-7).
_T5_LEFT = [(1, ""), (2, ""), (3, ""), (4, ""), (5, ""), (7, ""), (8, ""), (9, "*")]
_T5_RIGHT = [
    (1, ""),
    (2, ""),
    (3, ""),
    (4, ""),
    (5, ""),
    (6, ""),
    (7, ""),
    (8, "*"),
    (9, ""),
    (11, ""),
    (12, ""),
]

# Table VI: He, absolute experiments (rows 1-5) and Bishop & Pipin's values.
_add(
    "table6_He_dc_kerr_expt",
    "He",
    T6,
    [
        (1, 3, "ref 11", "632.8", ""),
        (2, 3, "ref 19", "514.5", ""),
        (3, 3, "ref 30", "632.8", ""),
        (4, 3, "ref 31", "632.8", ""),
    ],
)
_add("table6_He_thg_expt", "He", T6, [(5, 3, "ref 60", "1055", "")])
_add(
    "table6_He_dc_kerr_calc",
    "He",
    T6,
    [(1, 4, "632.8 nm", "632.8", ""), (2, 4, "514.5 nm", "514.5", "")],
)
_add("table6_He_thg_calc", "He", T6, [(5, 4, "1055 nm", "1055", "")])

# Table VIII: gamma^v and gamma^R by process (T = 295 K).
_T8_COLUMNS = [
    (1, "static", "inf"),
    (2, "dc_kerr", "632.8"),
    (3, "kerr", "632.8"),
    (4, "eshg", "694.3"),
    (5, "thg", "694.3"),
]
_T8_ROWS = [
    ("H2", "v", 5, ["e", "e", "f", "e", "e"]),
    ("D2", "v", 6, ["e", "e", "f", "e", "e"]),
    ("N2", "v", 8, ["f", "f", "f", "f", "f"]),
    ("O2", "v", 9, ["f", "f", "f", "f", "f"]),
    ("H2", "R", 12, ["c", "g", "h", "h", "h"]),
    ("N2", "R", 13, ["i", "i", "i", "i", "i"]),
]
for _gas, _part, _row, _marks in _T8_ROWS:
    for (_column, _process, _wavelength), _mark in zip(
        _T8_COLUMNS, _marks, strict=True
    ):
        _add(
            f"table8_{_gas}_gamma{_part}_{_process}",
            _gas,
            T8,
            [(_row, _column, f"gamma^{_part} {_process}", _wavelength, _mark)],
        )

# Table IX: electronic gamma of the diatomics (calc col 4, expt col 5).
for _gas, _row, _calc, _expt, _wavenumber in [
    ("H2", 3, "#,c", "#,d", "0"),
    ("D2", 4, "#,c", "#,d", "0"),
    ("N2", 6, "#,f", "#,d", "0"),
    ("O2", 7, None, "#,d", "0"),
    ("CO", 20, "#,p", "q", "14399"),
]:
    _wavelength = "inf" if _wavenumber == "0" else _wavenumber
    if _calc is not None:
        _add(
            f"table9_{_gas}_calc",
            _gas,
            T9,
            [(_row, 4, "calc", _wavelength, _calc)],
        )
    _add(f"table9_{_gas}_expt", _gas, T9, [(_row, 5, "expt", _wavelength, _expt)])

# Table XI: CO2 gamma^v (MP2) by process at omega = 0.07 au (15363 cm^-1).
for _column, _process in enumerate(["static", "dc_kerr", "kerr", "eshg", "thg"], 1):
    _add(
        f"table11_CO2_gammav_{_process}",
        "CO2",
        T11,
        [(5, _column, f"gamma^v {_process}", "inf" if _column == 1 else "15363", "")],
    )

# Table XII part b: CH4, SF6, CO2 (JSON rows 31-32, 37-38, 43-44).
_add("table12_CH4_calc_static", "CH4", T12, [(31, 1, "SCF (ref 227)", "inf", "*,w")])
_add(
    "table12_CH4_dc_kerr",
    "CH4",
    T12,
    [(31, 3, "ref 28", "632.8", "h"), (32, 3, "ref 22", "458", "@,aa")],
)
_add(
    "table12_CH4_eshg",
    "CH4",
    T12,
    [(31, 4, "ref 34", "694.3", "#,x"), (32, 4, "ref 55", "694.3", "m")],
)
_add("table12_CH4_thg", "CH4", T12, [(31, 5, "ref 62", "694.3", "y")])
_add("table12_CH4_cars", "CH4", T12, [(31, 6, "ref 73", "532", "z")])
_add(
    "table12_SF6_dc_kerr",
    "SF6",
    T12,
    [(37, 3, "ref 30", "632.8", "bb"), (38, 3, "ref 22", "632.8", "aa")],
)
_add(
    "table12_SF6_eshg",
    "SF6",
    T12,
    [(37, 4, "ref 34", "694.3", "#,x"), (38, 4, "ref 53", "694.3", "d")],
)
_add("table12_SF6_thg", "SF6", T12, [(37, 5, "ref 62", "694.3", "y")])
_add("table12_SF6_cars", "SF6", T12, [(37, 6, "ref 73", "532", "z")])
_add("table12_CO2_calc_static", "CO2", T12, [(43, 1, "CCSD(T) (ref 123)", "inf", "ff")])
_add("table12_CO2_calc_eshg", "CO2", T12, [(43, 2, "CCSD(T) (ref 123)", "694.3", "ff")])
_add("table12_CO2_dc_kerr", "CO2", T12, [(43, 3, "ref 18", "632.8", "gg")])
_add(
    "table12_CO2_eshg",
    "CO2",
    T12,
    [(43, 4, "ref 34", "694.3", "#,x"), (44, 4, "ref 53", "694.3", "d")],
)
_add("table12_CO2_thg", "CO2", T12, [(43, 5, "ref 64", "694.3", "hh")])
_add("table12_CO2_cars", "CO2", T12, [(43, 6, "ref 73", "532", "z")])


def _table5() -> list[Point]:
    """Return Table V's Ne values, left list then right list, as printed."""
    grid = json_table(DIRECTORY, T5)
    points = []
    for halves in ((0, _T5_LEFT), (4, _T5_RIGHT)):
        offset, rows = halves
        for row, mark in rows:
            cells = grid[row]
            method, basis, value, ref = cells[offset : offset + 4]
            number_, sigma, marks = split(value)
            assert marks == mark, (row, value)
            label = f"{method} {basis} (ref {ref})".replace("  ", " ")
            points.append(Point(label, "inf", number_, sigma, " | ".join(cells)))
    # The VPT cell "8, 65" (ref 179) holds two values: one row each.
    cells = grid[10]
    assert cells[4:8] == ["VPT", "", "8, 65", "179"], cells
    for suffix, value in (("a", "8"), ("b", "65")):
        points.append(
            Point(f"VPT (ref 179) {suffix}", "inf", value, "", " | ".join(cells))
        )
    return points


def extract() -> dict[str, list[Point]]:
    """Return one file per dataset."""
    files: dict[str, list[Point]] = {}
    grids: dict[str, list[list[str]]] = {}
    for stem, gas, table, cells in SPECS:
        grid = grids.setdefault(table, json_table(DIRECTORY, table))
        points = []
        for row, column, label, wavelength, expected in cells:
            value, sigma, marks = split(grid[row][column])
            if marks != expected:
                msg = (
                    f"{table} row {row} col {column}: marks {marks!r}, not {expected!r}"
                )
                raise ValueError(msg)
            points.append(Point(label, wavelength, value, sigma, " | ".join(grid[row])))
        files[f"{gas}/shelton_1994__{stem}.csv"] = points
    for row, refs in T4_REFS.items():
        if (grids[T4][row][4], grids[T4][row][6]) != refs:
            msg = f"Table IV row {row}: footnotes {grids[T4][row]}, not {refs}"
            raise ValueError(msg)
    files["Ne/shelton_1994__table5_Ne.csv"] = _table5()
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
