# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kingston (1964): calculated n - 1 and Verdet constants of Ar, Kr, Xe.

Table I: n - 1 (x 10^4, NTP) from adjusted oscillator strengths; ``...`` marks no value.
Table II: Verdet constants in micro-min/(Oe cm atm), each printed as a mantissa with
the power of ten in parentheses (``2.872(1)`` = 28.72; footnote a). The value column
keeps the mantissa and exponent as ``2.872e1``; ``printed`` keeps the cell.
Table III (merged header, so read from the JSON cells): Heddle's and Gill's measured
argon values beside the calculated ones, with the wavelength band after ``+-``.
"""

import json
import re

from _library import LIBRARY, Point, number, table, write

DIRECTORY = "Kingston_1964_RefractiveIndicesVerdet"
TABLE_I = "tbl_368145603687"
TABLE_II = "tbl_2ffb8b888ca1"
TABLE_III = "tbl_ff52cd5ab774"
GASES = ("Ar", "Kr", "Xe")
_POWER = re.compile(r"^([\d.]+)\((\d+)\)(\^a)?$")


def grid(table_id: str) -> list[list[str]]:
    """Return a library table's cell grid from its JSON file (merged cells once)."""
    path = LIBRARY / DIRECTORY / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))
    rows = [[""] * cells["columns"] for _ in range(cells["rows"])]
    for cell in cells["cells"]:
        rows[cell["row"]][cell["column"]] = cell["text"].strip()
    return rows


def _powered(cell: str) -> str:
    """Return ``2.872(1)`` as ``2.872e1`` (the footnote's power-of-ten notation)."""
    match = _POWER.match(cell)
    if match is None:
        msg = f"not a mantissa(power) cell: {cell!r}"
        raise ValueError(msg)
    return f"{match.group(1)}e{match.group(2)}"


def extract() -> dict[str, list[Point]]:
    """Return Tables I and II per gas and Table III's two argon columns."""
    files: dict[str, list[Point]] = {}
    for table_id, name, convert in (
        (TABLE_I, "theory", number),
        (TABLE_II, "verdet_theory", _powered),
    ):
        for cells in table(DIRECTORY, table_id)[1:]:
            for gas, cell in zip(GASES, cells[1:], strict=True):
                if cell == "...":
                    continue
                relative = f"{gas}/kingston_1964__{name}_{gas}.csv"
                files.setdefault(relative, []).append(
                    Point(cells[0], cells[0], convert(cell), "", " | ".join(cells))
                )
    heddle, gill, calculated = [], [], []
    for cells in grid(TABLE_III)[2:]:
        printed = " | ".join(cells)
        wavelength = cells[0].split("±")[0]
        measured, _, sigma = cells[1].partition("±")
        target = gill if wavelength == "1216" else heddle
        target.append(Point(wavelength, wavelength, number(measured), sigma, printed))
        value, _, sigma = cells[2].partition("±")
        calculated.append(Point(wavelength, wavelength, number(value), sigma, printed))
    files["Ar/kingston_1964__heddle_1963_Ar.csv"] = heddle
    files["Ar/kingston_1964__gill_1963_Ar.csv"] = gill
    files["Ar/kingston_1964__theory_table3_Ar.csv"] = calculated
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
