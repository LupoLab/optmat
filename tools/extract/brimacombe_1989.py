# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Brimacombe, Taylor & Leopold (1989): TPA of fused-silica fibres (Table II).

Table II (published only as JSON: merged header cells) gives the two-photon absorption
coefficient, called alpha_1 (cm/MW), of fused silica against wavelength: "This work"
(fibres, 193-351 nm), "Other work" and its "Reference". Each cell prints its own
uncertainty and power of ten, "2.0 \pm 0.8 \times 10^{-3}": the value and sigma are
stored in e-notation (``2.0e-3``, ``0.8e-3``; scale 1). "..." is an empty cell.

Checked against the table image (tbl_b42460a7da40.jpg): every cell agrees. The image
prints the "< 4.5 \times 10^{-5}" of Ref. 21 on its own line under the 266 nm row, with
no wavelength of its own: the extraction leaves the wavelength cell empty, and the
script gives the row the wavelength of the row above (266 nm), as the image and the
paper's Fig. 9 ("266 and 355 nm, Ref. 21") show. "<" cells are upper bounds: kept,
and excluded by the record. The footnote marker "a" (substrates rather than fibres) is
extracted as "^{a}" in the reference cell.

Table I (fibre types, core, cladding, OH content) holds no optical data; it is
described in the record's notes.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Brimacombe_1989_DependenceNonlinearTransmission"
TABLE = "tbl_b42460a7da40"

# The data file of each reference in the "Other work" column (the paper's numbering).
REFERENCES = {
    "9^{a}": "SiO2/brimacombe_1989__taylor_1988_248.csv",
    "6^{a}": "SiO2/brimacombe_1989__liu_1979_266.csv",
    "21^{a}": "SiO2/brimacombe_1989__liu_1978.csv",
    "8": "SiO2/brimacombe_1989__mizunami_1988_282.csv",
    "7": "SiO2/brimacombe_1989__mizunami_1987_308.csv",
}
THIS_WORK = "SiO2/brimacombe_1989__alpha1_fibre.csv"

_CELL = re.compile(
    r"^(?P<bound><)?\s*(?P<mantissa>[\d.]+)(?:\s*\\pm\s*(?P<sigma>[\d.]+))?"
    r"\s*\\times\s*10\^\{(?P<power>[-\d]+)\}$"
)


def parse(cell: str) -> tuple[str, str]:
    r"""Return the value and sigma of ``[<] m [\pm s] \times 10^{k}``, in e-notation."""
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"not a value cell: {cell!r}"
        raise ValueError(msg)
    power = match["power"]
    value = f"{match['mantissa']}e{power}"
    sigma = f"{match['sigma']}e{power}" if match["sigma"] else ""
    float(value)
    return value, sigma


def extract() -> dict[str, list[Point]]:
    """Return this work's alpha_1 and the other works' values, by reference."""
    grid = json_table(DIRECTORY, TABLE)
    if grid[1][1:] != ["This work", "Other work", "Reference"]:
        msg = f"unexpected sub-header {grid[1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    wavelength = ""
    for cells in grid[2:]:
        # A row without a wavelength continues the row above (266 nm, Ref. 21).
        wavelength = cells[0] or wavelength
        printed = " | ".join(cells)
        if cells[1] not in {"", "..."}:
            value, sigma = parse(cells[1])
            files.setdefault(THIS_WORK, []).append(
                Point(wavelength, wavelength, value, sigma, printed)
            )
        if cells[2] not in {"", "..."}:
            value, sigma = parse(cells[2])
            files.setdefault(REFERENCES[cells[3]], []).append(
                Point(wavelength, wavelength, value, sigma, printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
