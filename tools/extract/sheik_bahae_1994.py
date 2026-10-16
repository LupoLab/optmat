# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Sheik-Bahae, Wang and Van Stryland (1994), the unlabelled Table III.

Table III (``tbl_e5d825fda3b1``, p. 252; published only as JSON, merged cells) gives for
polycrystalline ZnSe and ZnS the linear index at 0.532 and 1.06 um, the gap, the
degenerate n2 and beta at 0.532 um and at 1.06 um, and the nondegenerate n2(2 omega;
omega) and beta(2 omega; omega) (probe 0.532 um, pump 1.06 um). Units from the caption:
n2 in 1e-14 cm2/W, beta in 1e-9 cm/W. The cells agree with the table image.

Layout (image and JSON): rows 0-2 are headers (the wavelength header of the
nondegenerate columns reads "(0.532 um, 1.06 um)"), each material spans two rows, the
second holding only the index at 1.06 um ("2.5 (1.06 um)"). Each index cell carries its
wavelength in parentheses, which becomes the abscissa. "N.A." (ZnS n2 at 0.532 um),
"<0.073" (ZnS nondegenerate beta, an upper limit) and the "0" of beta at 1.06 um (below
the two-photon edge, hbar omega < E_g / 2) are not measurements of a value and are left
out (the record's notes quote them).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "SheikBahae_1994_NondegenerateOpticalKerr"
TABLE_III = "tbl_e5d825fda3b1"

# Grid rows of each material: the main row, then the row with its second index.
MATERIALS = {"ZnSe": (3, 4), "ZnS": (5, 6)}

# The degenerate columns (index -> wavelength in um) and the nondegenerate ones, whose
# abscissa is the probe wavelength (0.532 um; the pump, 1.06 um, is a parameter).
DEGENERATE_N2 = {3: "0.532", 5: "1.06"}
DEGENERATE_BETA = {4: "0.532", 6: "1.06"}
NONDEGENERATE_N2, NONDEGENERATE_BETA, PROBE = 7, 8, "0.532"

_INDEX = re.compile(r"^\s*([\d.]+)\s*\(([\d.]+)\s*μm\)\s*$")


def _value(cell: str) -> str | None:
    """Return the number of a cell, or None for "N.A.", an upper limit or a zero."""
    text = cell.strip()
    if text in {"N.A.", "0", ""} or text.startswith("<"):
        return None
    return number(text)


def extract() -> dict[str, list[Point]]:
    """Return the indices, gaps, degenerate and nondegenerate n2 and beta."""
    grid = json_table(DIRECTORY, TABLE_III)
    files: dict[str, list[Point]] = {}
    for material, (main, second) in MATERIALS.items():
        cells = grid[main]
        printed = " | ".join(cells)

        def put(name: str, points: list[Point], material: str = material) -> None:
            if points:
                files[f"{material}/sheik_bahae_1994__{name}_{material}.csv"] = points

        indices = []
        for row in (main, second):
            match = _INDEX.match(grid[row][1])
            if match is None:
                msg = f"{material}: index cell {grid[row][1]!r}"
                raise ValueError(msg)
            value, wavelength = match.groups()
            indices.append(
                Point(wavelength, wavelength, value, "", " | ".join(grid[row]))
            )
        put("n0", indices)
        put("eg", [Point("E_g", "", number(cells[2]), "", printed)])
        for name, columns in (("n2_deg", DEGENERATE_N2), ("beta_deg", DEGENERATE_BETA)):
            points = []
            for column, wavelength in columns.items():
                value = _value(cells[column])
                if value is not None:
                    points.append(Point(wavelength, wavelength, value, "", printed))
            put(name, points)
        for name, column in (
            ("n2_nondeg", NONDEGENERATE_N2),
            ("beta_nondeg", NONDEGENERATE_BETA),
        ):
            value = _value(cells[column])
            if value is not None:
                put(name, [Point(PROBE, PROBE, value, "", printed)])
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
