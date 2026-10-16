# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Laporte, Damany & Damany (1975): Table I, oscillator strengths of Kr and Xe.

Table I (``tbl_8c662e0c2319``, cells file only; image checked, cells agree) compiles
the oscillator strengths of the resonance lines Kr 123.58 nm and Xe 146.96, 129.56,
125.02 and 119.22 nm (the column headers; the text calls the last 119.20 nm) from the
literature: an "Experimental" block (Refs. 32-42) and a "Theoretical" block (Refs.
43-45). Every value is another paper's: the files are reproductions. "..." is an empty
cell. Refs. 43 and 45 print two values per line, joined by a brace (``{`` in the
cells) on two rows, the second with no label: they become rows ``ref43a`` and
``ref43b``. The column's wavelength (nm, as printed in the header) is the abscissa.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Laporte_1975_HighDensityResonant"
TABLE_I = "tbl_8c662e0c2319"
HEADER = re.compile(r"^(?P<gas>Kr|Xe)`?(?P<nm>\d+\.\d+) nm$")
FILES = {
    ("Experimental", "Kr"): "Kr/laporte_1975__experimental_Kr.csv",
    ("Experimental", "Xe"): "Xe/laporte_1975__experimental_Xe.csv",
    ("Theoretical", "Kr"): "Kr/laporte_1975__theoretical_Kr.csv",
    ("Theoretical", "Xe"): "Xe/laporte_1975__theoretical_Xe.csv",
}
EXPECTED = {
    ("Experimental", "Kr"): 7,
    ("Experimental", "Xe"): 17,
    ("Theoretical", "Kr"): 5,
    ("Theoretical", "Xe"): 10,
}


def _reference(label: str) -> str:
    """Return a row label's reference key: 'Ref. 32' -> 'ref32', 'Refs. 39 and 40'."""
    numbers = re.findall(r"\d+", label)
    if not numbers:
        msg = f"Table I: no reference number in {label!r}"
        raise ValueError(msg)
    return "ref" + "_".join(numbers)


def extract() -> dict[str, list[Point]]:
    """Return one file per block and gas, each value with its reference and line."""
    rows = json_table(DIRECTORY, TABLE_I)
    columns = []
    for cell in rows[0][1:]:
        match = HEADER.match(cell)
        if match is None:
            msg = f"Table I: unexpected header {cell!r}"
            raise ValueError(msg)
        columns.append((match["gas"], match["nm"]))
    files: dict[str, list[Point]] = {name: [] for name in FILES.values()}
    block, reference, suffix = "", "", ""
    for row in rows[1:]:
        if row[0] in {"Experimental", "Theoretical"} and not any(row[1:]):
            block = row[0]
            continue
        braced = row[1].startswith("{")
        if row[0]:
            reference = _reference(row[0])
            suffix = "a" if braced else ""
        elif braced and suffix == "a":
            suffix = "b"  # the second row of a brace
        else:
            msg = f"Table I: an unlabelled row that is not a brace's second: {row}"
            raise ValueError(msg)
        for (gas, nm), cell in zip(columns, row[1:], strict=True):
            text = cell.removeprefix("{").strip()
            if text == "...":
                continue
            key = f"{reference}{suffix}" if gas == "Kr" else f"{reference}{suffix}_{nm}"
            files[FILES[block, gas]].append(
                Point(key, nm, number(text), "", " | ".join(row))
            )
    for name, count in EXPECTED.items():
        if len(files[FILES[name]]) != count:
            msg = f"Table I: {len(files[FILES[name]])} values for {name}, not {count}"
            raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
