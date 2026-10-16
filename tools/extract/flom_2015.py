# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Flom et al. (2015): n2 of six window materials at 772, 1030 and 1550 nm.

Table 1 (``tbl_5dfb2dc9eeb0``, published only as JSON) gives, per material and
wavelength, n2 in "10^-7 cm2/GW" (= 1e-16 cm2/W) with its +-, and the range of peak
intensities (GW/cm2) of the Z-scans averaged. The cells agree with the table image.
The intensity range is split into the columns ``intensity_min`` and ``intensity_max``
(one printed intensity gives both); its footnote letters ("15-144a", "16-228b": n2
increases, respectively decreases, with I) stay in ``printed``. "—" (Spinel at 772
nm, no n2 extracted) is left out. The pulse duration at each wavelength, from the
text (p. F125), is the per-row column ``pulse_fs``.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Flom_2015_UltrafastScanMeasurements"
TABLE = "tbl_5dfb2dc9eeb0"

MATERIALS = {
    "Fused silica": "SiO2",
    "BK7": "BK7",
    "AlON": "AlON",
    "Spinel": "MgAl2O4",
    "BGG": "BGG",
    "CLEARTRAN": "ZnS",
}

# FWHM by SHG FROG (p. F125): 275 +- 20 fs at 772 nm, 140 +- 10 fs at 1030 nm,
# 97 +- 7 fs at 1550 nm.
PULSE_FS = {"772": "275", "1030": "140", "1550": "97"}

# Column of each wavelength's n2; the intensity range is the next column.
N2_COLUMNS = {1: "772", 3: "1030", 5: "1550"}


def _intensities(cell: str) -> tuple[str, str]:
    """Return the lowest and highest intensity of a "lo-hi" cell (en dash) or "v"."""
    text = cell.strip().rstrip("ab")
    low, _, high = text.partition("\N{EN DASH}")
    return number(low), number(high or low)


def extract() -> dict[str, list[Point]]:
    """Return one data file per material."""
    files = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        material = MATERIALS[cells[0]]
        printed = " | ".join(cells)
        points = []
        for column, wavelength in N2_COLUMNS.items():
            if cells[column].strip() == "—":
                continue
            value, _, sigma = cells[column].partition("±")
            low, high = _intensities(cells[column + 1])
            extra = {
                "pulse_fs": PULSE_FS[wavelength],
                "intensity_min": low,
                "intensity_max": high,
            }
            points.append(
                Point(
                    wavelength,
                    wavelength,
                    number(value),
                    number(sigma) if sigma else "",
                    printed,
                    extra,
                )
            )
        files[f"{material}/flom_2015__n2_{material}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
