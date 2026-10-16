# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Burns, Graham & Weller (1986): A_R and B_R at 457.9-632.8 nm.

Tables 2, 3, 6 and 12.

The tables have no CSV rendering (merged headings) and are read from their JSON cells,
checked against the images. Numbers print a raised decimal point (``0·5202``), turned
into ``.`` here; the ``printed`` column keeps the raw cells.

Table 2, ``10^6 A_R / m^3 mol^-1`` (= cm^3/mol): each cell is this work's value, then
in parentheses a literature value and its reference numbers, ``6·362(6·349)[22]``.
Table 3, ``10^12 B_R / m^6 mol^-2`` (= cm^6/mol^2): ``value ± uncertainty``, then any
number of literature values ``(value ± uncertainty)[ref]``; "—" is no measurement.
Tables 6 and 12: the calculated B_R ("Calc.") at 632.8 nm and 298.2 K beside the
observed one (which repeats Table 3). The fluorinated methanes (not curated species)
are skipped.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Burns_1986_DirectMeasurementCalculation"
TABLE_2 = "tbl_f1e1a8af1bc0"
TABLE_3 = "tbl_04b95e2c3c8c"
TABLE_6 = "tbl_e61576675901"
TABLE_12 = "tbl_b91da9afcf21"
GASES = {
    "He": "He",
    "Ne": "Ne",
    "Ar": "Ar",
    "Kr": "Kr",
    "Xe": "Xe",
    "CH_4": "CH4",
    "CO_2": "CO2",
    "N_2": "N2",
    "SF_6": "SF6",
}
_A_CELL = re.compile(r"^([-\d.]+)(?:\(([-\d.]+)\)\[([\d, ]+)\])?$")
_B_OWN = re.compile(r"^([-\d.]+) ± ([\d.]+)")
_B_LITERATURE = re.compile(r"\(([-\d.]+) ± ([\d.]+)\)\[(\d+)\]")


def _text(cell: str) -> str:
    """Return a cell with its raised decimal point as ``.``, stripped."""
    return cell.replace("·", ".").strip()


def extract() -> dict[str, list[Point]]:
    """Return A_R and B_R per gas, the literature values and the calculated B_R."""
    files: dict[str, list[Point]] = {}

    def add(gas: str, name: str, point: Point) -> None:
        files.setdefault(f"{gas}/burns_1986__{name}.csv", []).append(point)

    rows = json_table(DIRECTORY, TABLE_2)
    heads = rows[1]
    for cells in rows[2:]:
        wavelength, printed = _text(cells[0]), " | ".join(cells)
        for head, cell in zip(heads[1:], cells[1:], strict=True):
            if head not in GASES:
                continue
            match = _A_CELL.match(_text(cell))
            if match is None:
                msg = f"unexpected Table 2 cell {cell!r}"
                raise ValueError(msg)
            gas = GASES[head]
            add(
                gas,
                f"A_R_{gas}",
                Point(wavelength, wavelength, number(match[1]), "", printed),
            )
            add(
                gas,
                f"A_R_literature_{gas}",
                Point(wavelength, wavelength, number(match[2]), "", printed),
            )

    rows = json_table(DIRECTORY, TABLE_3)
    heads = rows[1]
    for cells in rows[2:]:
        wavelength, printed = _text(cells[0]), " | ".join(cells)
        for head, cell in zip(heads[1:], cells[1:], strict=True):
            text = _text(cell)
            if head not in GASES or text == "—":
                continue
            gas = GASES[head]
            own = _B_OWN.match(text)
            if own is None:
                msg = f"unexpected Table 3 cell {cell!r}"
                raise ValueError(msg)
            add(
                gas,
                f"B_R_{gas}",
                Point(wavelength, wavelength, number(own[1]), own[2], printed),
            )
            for value, sigma, ref in _B_LITERATURE.findall(text):
                add(
                    gas,
                    f"B_R_ref{ref}_{gas}",
                    Point(wavelength, wavelength, number(value), sigma, printed),
                )

    for table_id, column in ((TABLE_6, 4), (TABLE_12, 5)):
        for cells in json_table(DIRECTORY, table_id)[2:]:
            if cells[0] not in GASES:
                continue
            gas = GASES[cells[0]]
            add(
                gas,
                f"B_R_theory_{gas}",
                Point(
                    "632.8",
                    "632.8",
                    number(_text(cells[column])),
                    "",
                    " | ".join(cells),
                ),
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
