# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bose & Cole (1970): dielectric and pressure virials of CO2, Ar and mixtures.

Table II (``tbl_d868b7eaf734``, CSV; image checked, no differences) gives A_eps
(cm^3 mol^-1), B_eps (cm^6 mol^-2) and C_eps (cm^9 mol^-3) of the pure gases, each
printed ``value±error``, with the temperature in degrees Celsius and the source ("This
work" or a reference number). Rows of this work: A_eps to one file per gas (row = the
temperature in degrees C as printed), B_eps and C_eps to one file per gas and
temperature. The literature rows (Refs. 1, 3, 4, 5) are reproduced, one file per
quantity and reference. "..." marks a value not given (skipped); Ref. 4 prints its CO2
values without an error.

Table III (``tbl_a3d065b6813b``, cells file only; image checked, no differences): the
pressure second virial coefficients B_p (cm^3 mol^-1) of this work against the mole
fraction X of CO2, at 49.7 degrees C except the X = 0.313 mixture (footnote b, 47.7
degrees C; the footnote markers ``a``, ``b``, ``e`` are stripped from the cells). X = 0
and 1.000 are pure Ar and CO2 (``B_p_Ar``, ``B_p_CO2``); each mixture is a file of its
own, filed under Ar (the first of the pair in sorted order) with Ar's mole fraction
1 - X in the dataset id (``B_p_Ar_0.874_CO2`` for X = 0.126). The rows are the
temperatures in degrees C. The literature and calculated columns are not extracted
(the record's notes).

Table IV (``tbl_c9ddaa569d7c``, cells file only; image checked): only the experimental
cross coefficient B_p12 of the CO2-Ar row (``B_p12_Ar_CO2``), at 49.7 degrees C (the
caption). Its experimental B_p of the pure gases repeat Table III's.

Table I (the mixtures' dielectric coefficients) is not extracted: the vocabulary has no
quantity for a mixture's dielectric virials (the record's notes).
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Bose_1970_DielectricPressureVirial"
TABLE_II = "tbl_d868b7eaf734"
TABLE_III = "tbl_a3d065b6813b"
TABLE_IV = "tbl_c9ddaa569d7c"
CELSIUS = "49.7"
_PM = "\N{PLUS-MINUS SIGN}"
GASES = {"CO2": "CO2", "Ar": "Ar"}
REFERENCES = {
    "1": "orcutt_1967",
    "3": "orcutt_cole_1965",
    "4": "johnston_cole_1962",
    "5": "johnston_1960",
}
QUANTITIES = {2: "A_eps", 3: "B_eps", 4: "C_eps"}


def _split(cell: str) -> tuple[str, str]:
    """Return the value and the error of ``7.358±0.002`` (the error "" if none)."""
    value, _, sigma = cell.partition(_PM)
    return number(value), number(sigma) if sigma else ""


def _dielectric() -> dict[str, list[Point]]:
    """Return Table II."""
    grid = table(DIRECTORY, TABLE_II)
    if grid[0][0] != "Gas" or grid[0][5] != "Reference":
        msg = f"unexpected Table II header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        gas = GASES[cells[0]]
        celsius = number(cells[1])
        printed = " | ".join(cells)
        source = cells[5]
        for column, name in QUANTITIES.items():
            cell = cells[column]
            if cell == "...":
                continue
            point = Point(celsius, "inf", *_split(cell), printed)
            if source == "This work":
                stem = f"A_eps_{gas}" if name == "A_eps" else f"{name}_{celsius}_{gas}"
                files.setdefault(f"{gas}/bose_1970__{stem}.csv", []).append(point)
            else:
                stem = f"{name}_{REFERENCES[source]}_{gas}"
                files[f"{gas}/bose_1970__{stem}.csv"] = [point]
    return files


def _pressure() -> dict[str, list[Point]]:
    """Return Table III's B_p of this work, one file per gas or mixture."""
    grid = json_table(DIRECTORY, TABLE_III)
    if not grid[0][0].startswith("Mole fraction X of CO2") or grid[1][1] != "This work":
        msg = f"unexpected Table III header {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        marked = cells[0].rstrip("ab")
        # footnote b: the X = 0.313 mixture was measured at 47.7 degrees C
        celsius = "47.7" if cells[0].endswith("b") else CELSIUS
        printed = " | ".join(cells)
        value, sigma = _split(cells[1])
        point = Point(celsius, celsius, value, sigma, printed)
        if marked == "0":
            stem = "B_p_Ar"
        elif marked == "1.000":
            stem = "B_p_CO2"
        else:
            fraction = f"{1 - float(marked):.3f}"  # Ar's mole fraction, 1 - X(CO2)
            stem = f"B_p_Ar_{fraction}_CO2"
        gas = "CO2" if stem == "B_p_CO2" else "Ar"
        files[f"{gas}/bose_1970__{stem}.csv"] = [point]
    return files


def _cross() -> dict[str, list[Point]]:
    """Return Table IV's experimental B_p12 of CO2-Ar."""
    grid = json_table(DIRECTORY, TABLE_IV)
    if grid[1][4:6] != ["Calc", "Exptl"]:
        msg = f"unexpected Table IV header {grid[:2]!r}"
        raise ValueError(msg)
    (cells,) = [cells for cells in grid[2:] if cells[0] == "CO_2-Ar"]
    point = Point(CELSIUS, CELSIUS, number(cells[5]), "", " | ".join(cells))
    return {"Ar/bose_1970__B_p12_Ar_CO2.csv": [point]}


def extract() -> dict[str, list[Point]]:
    """Return Tables II, III and IV."""
    return _dielectric() | _pressure() | _cross()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
