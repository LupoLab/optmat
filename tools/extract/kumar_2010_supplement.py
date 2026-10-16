# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kumar & Thakkar (2010), supplement: integrated DOS of Ne, Ar, Kr and Xe.

Table I (``tbl_8fedf8f12d10``, CSV; image checked, every cell) lists, per gas, the
interval edges E_j (eV) and the integrated dipole oscillator strength (IDOS) of the
adopted constrained DOSD over [E_j, E_{j+1}] (the last interval ends at infinity), and
in its first row the summed discrete oscillator strength. "The numbers in parentheses
are powers of ten": ``2.156(1)`` is written 2.156e1, the printed cell kept. Per gas,
one file of intervals (rows j, abscissa E_j, the lower edge) and one of the discrete
sum (index abscissa). The gas headers carry footnote letters (``Nea``); the column
pairs are fixed by position. The ten-term pseudospectra (page 3) are recorded as
formulas in the source record, not extracted here.
"""

import dataclasses
import re

from _library import Point, table, write

DIRECTORY = "Kumar_2010_DipoleOscillatorStrength/supplement_01"
TABLE = "tbl_8fedf8f12d10"
GASES = {1: "Ne", 3: "Ar", 5: "Kr", 7: "Xe"}  # column of E_j; IDOS is the next one
_CELL = re.compile(r"^(?P<mantissa>\d+\.\d+)(?:\((?P<exponent>-?\d+)\))?$")


def _value(cell: str) -> str:
    """Return a ``mantissa(power of ten)`` cell in e-notation."""
    match = _CELL.match(cell.strip().replace("\N{MINUS SIGN}", "-"))
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    mantissa, exponent = match["mantissa"], match["exponent"]
    return f"{mantissa}e{exponent}" if exponent else mantissa


def extract() -> dict[str, list[Point]]:
    """Return the IDOS intervals and the discrete sums per gas."""
    rows = table(DIRECTORY, TABLE)
    if rows[1][:3] != ["j", "Ej", "IDOS"]:
        msg = f"unexpected header {rows[1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[2:]:
        label = cells[0].strip()
        for column, gas in GASES.items():
            energy, idos = cells[column].strip(), cells[column + 1].strip()
            if not idos:
                continue
            printed = f"{gas} | {label} | {energy} | {idos}"
            if label.startswith("Discrete"):
                path = f"{gas}/kumar_2010_supplement__discrete_{gas}.csv"
                files[path] = [Point("discrete", "1", _value(idos), "", printed)]
                continue
            path = f"{gas}/kumar_2010_supplement__idos_{gas}.csv"
            files.setdefault(path, []).append(
                Point(label, _value(energy), _value(idos), "", printed)
            )
    return {
        path: _with_upper_edges(points) if "__idos_" in path else points
        for path, points in files.items()
    }


def _with_upper_edges(points: list[Point]) -> list[Point]:
    """Close each interval at the next row's lower edge, the last at infinity.

    E_j is the lower edge of interval j (it matches kumar_1985's explicit intervals:
    Ar 15.76 eV, f = 0.0494 against 15.76-15.94 eV, 0.0490), and the last interval
    runs to infinity as kumar_1985's 10^5-infinity row does.
    """
    edges = [p.x for p in points[1:]] + ["inf"]
    return [
        dataclasses.replace(p, extra={"x_upper": edge})
        for p, edge in zip(points, edges, strict=True)
    ]


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
