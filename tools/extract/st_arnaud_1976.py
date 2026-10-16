# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract St-Arnaud & Bose (1976): the refractivity virial coefficients of CH4.

The one table (``tbl_e86bfa5ea90a``, CSV; image checked, cells agree) lists, per
author, the temperature (K), the wavelength (A; "..." for the theory row), A_R
(cm^3/mol), B_P (cm^3/mol) and B_R (cm^6/mol^2), each as "value +- sigma". A_R goes to
``molar_refractivity`` files, B_R to ``refractivity_virial_coefficient`` files, one per
author; B_P (a density virial coefficient) has no quantity and is not extracted. The
theory row (Kirkwood) has no wavelength and its B_R carries the footnote mark "a"
(stripped); it is filed at "inf".
"""

import re

from _library import Point, number, table, write

DIRECTORY = "StArnaud_1976_DirectDeterminationSecond"
TABLE = "tbl_e86bfa5ea90a"
AUTHORS = {
    "Olson": "olson_1975",
    "Buckingham and Graham": "buckingham_graham_1974",
    "Our work": "this_work",
    "Kirkwood (theory)": "kirkwood_theory",
}
_VALUE = re.compile(r"^\s*(?P<value>-?[\d.]+)\s*(?:±\s*(?P<sigma>[\d.]+))?\s*a?\s*$")


def _split(cell: str) -> tuple[str, str] | None:
    """Return (value, sigma) of a "value +- sigma" cell, or None for "..."."""
    if cell.strip() == "...":
        return None
    match = _VALUE.match(cell)
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    return number(match["value"]), match["sigma"] or ""


def extract() -> dict[str, list[Point]]:
    """Return A_R and B_R of each author as data files."""
    rows = table(DIRECTORY, TABLE)
    assert rows[0][3].startswith("A_R") and rows[0][5].startswith("B_R"), rows[0]
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        name = AUTHORS[cells[0].strip()]
        wavelength = "inf" if cells[2].strip() == "..." else number(cells[2])
        printed = " | ".join(cells)
        for column, quantity in ((3, "A_R"), (5, "B_R")):
            parsed = _split(cells[column])
            if parsed is None:
                continue
            value, sigma = parsed
            files[f"CH4/st_arnaud_1976__{quantity}_{name}_CH4.csv"] = [
                Point(wavelength, wavelength, value, sigma, printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
