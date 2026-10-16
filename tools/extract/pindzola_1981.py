# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Pindzola (1981): photoionisation cross sections of excited Xe 6p (J = 0).

The PDF is a scan with a text layer (Phys. Rev. A 23, 201: journal page = PDF page +
200). **Table I** (p. 204) is printed in four parts, each a library table without a
merged header (the parts' headings are outside the crops; their order was checked on
the rendered page):

* A, the ``2P3/2 es (J = 1)`` partial cross section (``tbl_933b3111b7ca``, JSON only);
* B, ``2P3/2 ed-bar (J = 1)`` (``tbl_e18e3ec44547``);
* C, ``2P3/2 ed (J = 1)`` (``tbl_6fcf7de022b4``);
* D, the total (``tbl_a31586f75c71``).

Each part has seven photon energies (2.70 to 18.20 eV) and four columns: Hartree-Fock
length and velocity (HFL, HFV) and Dirac-Fock length and velocity (DFL, DFV), in cm^2,
printed as ``1.367^{-18}`` ('Denotes the power of ten by which the entry is to be
multiplied', note b). Two cells carry a footnote mark (``1.367^{-18}b`` and
``1.283^{-18}a``, row 2.70 of part A), dropped from the value and kept in ``printed``.
Every cell was checked against the table images and the rendered page.

:func:`checks` adds the three partial cross sections and compares the sum with the
printed total, within the rounding of the four printed numbers. One cell does not
close (DFV at 6.20 eV: the partials add to 8.7386e-19, the total prints 8.735e-19,
3.6 units of its last figure where rounding allows 1.1); it is listed in
``NOT_CLOSING`` and the check fails if any other cell does not close, or this one does.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Pindzola_1981_PhotoionizationExcitedState"
PARTS = {
    "es": "tbl_933b3111b7ca",
    "dbar": "tbl_e18e3ec44547",
    "d": "tbl_6fcf7de022b4",
    "total": "tbl_a31586f75c71",
}
"""Table I's parts A, B, C and D, by final channel (the total last)."""
METHODS = ("hfl", "hfv", "dfl", "dfv")
"""The columns after the photon energy, in order (headers HFL, HFV, DFL, DFV)."""
ENERGIES = ("2.70", "3.20", "4.20", "4.97", "6.20", "10.20", "18.20")
"""The photon energies (eV) printed in every part."""
CELL = re.compile(r"^(\d\.\d{3})\^\{(-\d{2})\}([ab]?)$")
"""A cell: mantissa, power of ten (note b) and an optional footnote mark."""
NOT_CLOSING = {("6.20", "dfv")}
"""The one cell whose partials do not add to the printed total within rounding."""


def _number(cell: str) -> tuple[str, int]:
    """Return a cell's value in cm^2 as text, and the exponent of its last figure."""
    match = CELL.match(cell)
    if match is None:
        msg = f"Table I: an unexpected cell {cell!r}"
        raise ValueError(msg)
    mantissa, power, _ = match.groups()
    return f"{mantissa}e{power}", int(power) - 3


def tables() -> dict[str, list[list[str]]]:
    """Return each part's data rows (photon energy, then the four columns)."""
    out = {}
    for part, table_id in PARTS.items():
        rows = [row for row in json_table(DIRECTORY, table_id) if row[0] in ENERGIES]
        if [row[0] for row in rows] != list(ENERGIES) or any(
            len(row) != 1 + len(METHODS) for row in rows
        ):
            msg = f"Table I, part {part}: rows or columns moved"
            raise ValueError(msg)
        out[part] = rows
    return out


def checks() -> dict[str, object]:
    """Return the closure of the partial cross sections on the totals.

    Raises unless every printed total is the sum of its three partials within the
    rounding of the four numbers, except the cells of ``NOT_CLOSING``, which must not
    close.
    """
    parts = tables()
    out: dict[str, object] = {}
    failures = []
    for i, energy in enumerate(ENERGIES):
        for j, method in enumerate(METHODS, start=1):
            total, last = _number(parts["total"][i][j])
            terms = [_number(parts[p][i][j]) for p in ("es", "dbar", "d")]
            bound = 0.5 * (10.0**last + sum(10.0**e for _, e in terms))
            excess = sum(float(v) for v, _ in terms) - float(total)
            units = excess / 10.0**last
            if (abs(excess) > bound) != ((energy, method) in NOT_CLOSING):
                failures.append(f"{energy} eV {method}: off by {units:.2f} units")
            out[f"{energy} eV {method}: partials - total (last-figure units)"] = round(
                units, 2
            )
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return Table I's four parts, one file per part and method."""
    checks()
    files: dict[str, list[Point]] = {}
    for part, rows in tables().items():
        for j, method in enumerate(METHODS, start=1):
            files[f"Xe/pindzola_1981__{part}_{method}.csv"] = [
                Point(row[0], row[0], _number(row[j])[0], "", " | ".join(row))
                for row in rows
            ]
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
