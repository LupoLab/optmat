# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Gornik, Kindt, Matthias & Schmidt (1981): Xe two-photon pressure shifts.

Table II (p. 72, library table ``tbl_2928c3b46e82``), 'Measured line shift parameters
for two-photon transitions into the 6p configuration', prints the shift parameter
beta/N in 10^-9 s^-1 cm^3, 'beta is the line shift in angular frequencies and N the
number of atoms per cm^3' (p. 72), as 'value (uncertainty)', for the lines from 1S0 to
6p[1/2]0, 6p[3/2]2 and 6p[5/2]2 (the paper writes 6p|1/2|0; the second and third rows
print only '-6p|3/2|2', '-6p|5/2|2' under the first row's '1S0'). One file, the
``pressure_line_coefficient`` with effect shift (approved in TC-1040). The
cells were compared with the table image (2026-10-09): no misread cell.

**Checks** (:func:`checks`): the three values and uncertainties must be the abstract's
'-(2.5 +- 0.2) x 10^-9, -(2.7 +- 0.2) x 10^-9, and -(2.8 +- 0.2) x 10^-9 s^-1 cm^3'
(p. 68), and the 6p[1/2]0 row's Delta C6 must be the text's '2.0 x 10^7 cm^-1
Angstrom^6' (p. 72). Reported: each shift in ordinary frequency per bar and per torr
at 295 K (the scale of Fig. 5's 'Line Shift / GHz' against 'Gas Pressure / bar').
"""

import math
import re

from _library import Point, table, write

from optmat.constants import K_B

DIRECTORY = "Gornik_1981_TwoPhotonExcitation"
TABLE_II = "tbl_2928c3b46e82"
FILE = "Xe/gornik_1981__beta_shift_Xe.csv"
HEADER_SHIFT = r"Shift \beta/N(10^{-9} \text{ s}^{-1} \text{ cm}^{3})"
ROWS = {
    r"^{1}S_{0} - 6p \mid 1/2 \mid_{0}": "6p[1/2]0",
    r"-6p \mid 3/2 \mid_{2}": "6p[3/2]2",
    r"-6p \mid 5/2 \mid_{2}": "6p[5/2]2",
}
"""Table II's transition cells and the record's row labels (jl notation)."""
_CELL = re.compile(r"^(-?\d+\.\d+) \((\d+\.\d+)\)$")
ABSTRACT = {"6p[1/2]0": (-2.5, 0.2), "6p[3/2]2": (-2.7, 0.2), "6p[5/2]2": (-2.8, 0.2)}
"""p. 68, the abstract's beta/N, 10^-9 s^-1 cm^3."""
TEXT_DELTA_C6 = "2.0"
"""p. 72: 'Delta C6 = 2.0 x 10^7 cm^-1 Angstrom^6' for 1S0 -> 6p[1/2]0."""
TEMPERATURE_K = 295.0
"""Room temperature for the reported conversions (the paper states none)."""


def _rows() -> list[list[str]]:
    """Return Table II's three body rows, the header checked."""
    grid = table(DIRECTORY, TABLE_II)
    if grid[0][1] != HEADER_SHIFT or [row[0] for row in grid[1:]] != list(ROWS):
        msg = f"Table II has changed: {grid}"
        raise ValueError(msg)
    return grid[1:]


def extract() -> dict[str, list[Point]]:
    """Return Table II's shift parameters with their printed uncertainties."""
    points = []
    for row in _rows():
        match = _CELL.match(row[1])
        if match is None:
            msg = f"Table II: unexpected shift cell {row[1]!r}"
            raise ValueError(msg)
        label = ROWS[row[0]]
        points.append(
            Point(label, label, match.group(1), match.group(2), " | ".join(row[:2]))
        )
    return {FILE: points}


def checks() -> dict[str, str]:
    """Raise unless Table II agrees with the abstract and the text; report more."""
    points = extract()[FILE]
    failures = [
        f"{p.row}: {p.value} ({p.sigma}) against the abstract's {ABSTRACT[p.row]}"
        for p in points
        if (float(p.value), float(p.sigma)) != ABSTRACT[p.row]
    ]
    if _rows()[0][2] != TEXT_DELTA_C6:
        failures.append(f"Delta C6 {_rows()[0][2]!r} is not the text's 2.0")
    if failures:
        raise ValueError("; ".join(failures))
    per_bar = 1e5 / (K_B * TEMPERATURE_K) * 1e-6
    per_torr = per_bar * 101325 / 760 / 1e5
    return {
        p.row: (
            f"{float(p.value) * 1e-9 * per_bar / (2 * math.pi) / 1e9:.2f} GHz/bar, "
            f"{float(p.value) * 1e-9 * per_torr / (2 * math.pi) / 1e6:.1f} MHz/torr "
            f"(ordinary frequency, {TEMPERATURE_K:.0f} K)"
        )
        for p in points
    }


def main() -> None:
    """Write the data file and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
