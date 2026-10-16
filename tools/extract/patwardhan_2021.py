# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Patwardhan et al. (2021): mid-infrared n2 and the two-band theory (Table 1).

Table 1 (published only as cells, read with ``_json_table``) has, per material, the band
gap used by the theory and, at 2.3, 3.5, 4.5 and 8.3 um, the measured n2 ("Exp.",
value ± sigma) and the two-band-model n2 ("Theory"), in 10^-20 m^2/W. The wavelength
sits in the merged header cell over each Exp./Theory pair (the image,
``tables/tbl_2daca1ccb18c.jpg``, agrees with the cells); a dash marks a wavelength that
was not measured and is not extracted.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Patwardhan_2021_NonlinearRefractiveIndex"
TABLE = "tbl_2daca1ccb18c"

# The first cell as printed -> the registered material id.
MATERIALS = {
    "ZnSe": "ZnSe",
    "ZnS (polycrystalline)": "ZnS",
    "Sapphire": "Al2O3",
    "CaF_2": "CaF2",
    "MgF_2": "MgF2",
    "Fused silica": "SiO2",
}
DASH = "-"


def extract() -> dict[str, list[Point]]:
    """Return the measured and the theoretical n2, and the band gap, per material."""
    rows = json_table(DIRECTORY, TABLE)
    header, kinds = rows[0], rows[1]
    columns = []  # (column index, wavelength in um as printed, "exp" or "theory")
    for index in range(2, len(header)):
        wavelength = header[index] or header[index - 1]
        kind = {"Exp.": "exp", "Theory": "theory"}[kinds[index]]
        columns.append((index, number(wavelength.removesuffix("μm")), kind))
    files: dict[str, list[Point]] = {}
    for row, cells in enumerate(rows[2:], start=1):
        material = MATERIALS[cells[0]]
        printed = " | ".join(cells)
        files[f"{material}/patwardhan_2021__band_gap_{material}.csv"] = [
            Point(str(row), "", number(cells[1]), "", printed)
        ]
        for index, wavelength, kind in columns:
            cell = cells[index].strip()
            if cell == DASH:
                continue
            value, _, sigma = cell.partition("±")
            relative = f"{material}/patwardhan_2021__{kind}_{material}.csv"
            files.setdefault(relative, []).append(
                Point(
                    f"{row}-{wavelength}",
                    wavelength,
                    number(value),
                    number(sigma) if sigma else "",
                    printed,
                )
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
