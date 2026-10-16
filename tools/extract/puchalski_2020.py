# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Puchalski, Szalewicz, Lesiuk & Jeziorski (2020): the polarizability of 4He.

Table VIII (``tbl_d93800f29266``; image checked, the cells agree) lists the
contributions to the static dipole polarizability of 4He in atomic units (a0^3), their
total, the molar polarizability (4 pi/3) alpha_d N_A in cm^3/mol and Gaiser and
Fellmuth's (2018) measured value. Only the three rows that are values of the whole
quantity are extracted: "Total" (a0^3), the molar polarizability and the experiment
(cm^3/mol). The component rows (nonrelativistic, alpha^2, alpha^2/M, alpha^3 QED in
two parts, alpha^3/M, alpha^4, finite nuclear size) are corrections, not totals, and
are recorded in the source record's notes. The cells print the digits in groups of
three separated by a LaTeX thin space (``\ ``) and carry footnote markers
(``^{\text{f,g}}``): both are removed before the number is read; the ``printed``
column keeps the raw cells. The parenthesised last-digit uncertainty goes to the sigma
column. Tables I-VII (convergence of the electric-field derivative of the Bethe
logarithm) hold no optical data.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Puchalski_2020_QEDCalculationDipole"
TABLE_VIII = "tbl_d93800f29266"
# Row label as printed (its start) -> data file.
ROWS = {
    "Total": "He/puchalski_2020__alpha_static.csv",
    "Molar polarizability": "He/puchalski_2020__molar_polarizability.csv",
    "Experiment, Ref. [58]": "He/puchalski_2020__gaiser_2018.csv",
}
_FOOTNOTE = re.compile(r"\^\{\\text\{[a-z,]+\}\}$")
_VALUE = re.compile(r"^([\d.]+)\((\d+)\)$")


def _split(cell: str) -> tuple[str, str]:
    """Return a value and its last-digit uncertainty in the value's unit."""
    text = _FOOTNOTE.sub("", cell.strip()).replace("\\ ", "").replace(" ", "")
    match = _VALUE.match(text)
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    value, digits = match.groups()
    decimals = len(value.split(".")[1])
    return number(value), f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return the total, the molar polarizability and the quoted experiment."""
    rows = table(DIRECTORY, TABLE_VIII)
    assert rows[0] == ["Contribution", "Value"], rows[0]
    files: dict[str, list[Point]] = {}
    for label, cell in rows[1:]:
        (file,) = [f for start, f in ROWS.items() if label.startswith(start)] or [None]
        if file is None:
            continue  # a component, recorded in the notes
        value, sigma = _split(cell)
        files[file] = [Point("static", "inf", value, sigma, f"{label} | {cell}")]
    assert len(files) == len(ROWS), files
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
