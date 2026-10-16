# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Lehmberg et al. (1995): the nonlinear refractive indices of Table 1.

Table 1 (unlabelled in the extraction, ``tbl_30b21d3fcab2``, CSV; image checked, p. 87)
prints n2 in 10^-14 esu at 760 Torr and 300 K: xenon at 248.4 nm by the ab initio
calculation (Eqs. 1 or 3), by the self-defocusing experiment and by the CO2
compensation experiment, and CO2 at 248.4 nm by the self-focusing experiment, each as
'value \pm error'. The rows for air (Ref. [2]) and fused silica (Ref. [1]) are other
papers' values and are not taken.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Lehmberg_1995_TwoPhotonResonantly"
TABLE_1 = "tbl_30b21d3fcab2"
HEADER = ["Medium (source)", "n_{2} (10^{-14} esu)"]
# Table 1's row -> (material, dataset stem).
ROWS = {
    "Xe (Eqs. (1) or (3))": ("Xe", "calculated_Xe"),
    "Xe (self-defocusing exp.)": ("Xe", "defocusing_Xe"),
    "Xe (CO_{2} compensation exp.)": ("Xe", "compensation_Xe"),
    "CO_{2} (self-focusing exp.)": ("CO2", "focusing_CO2"),
}
SKIPPED = {"Air (Ref. [2])", "Fused SiO_{2} (Ref. [1])"}
WAVELENGTH_NM = "248.4"
_CELL = re.compile(r"^(?P<value>[-+]?[\d.]+) \\pm (?P<error>[\d.]+)$")


def extract() -> dict[str, list[Point]]:
    """Return the four n2 values of Table 1 with their errors."""
    grid = table(DIRECTORY, TABLE_1)
    if grid[0] != HEADER:
        msg = f"unexpected Table 1 header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        label = cells[0].strip()
        if label in SKIPPED:
            continue
        gas, stem = ROWS[label]
        match = _CELL.match(cells[1].strip())
        if match is None:
            msg = f"unexpected Table 1 cell {cells[1]!r}"
            raise ValueError(msg)
        files[f"{gas}/lehmberg_1995__{stem}.csv"] = [
            Point(
                "1",
                WAVELENGTH_NM,
                number(match["value"]),
                number(match["error"]),
                " | ".join(cells),
            )
        ]
    if len(files) != len(ROWS):
        msg = f"Table 1 rows missing: {sorted(files)}"
        raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
