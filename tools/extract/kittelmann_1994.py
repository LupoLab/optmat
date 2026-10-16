# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kittelmann & Ringling (1994): two-photon absorption at 193 nm (Table 1).

Table 1 (published only as JSON: its header has merged cells) gives, per material,
the sample ("probe") thickness, the band gap, beta at 193 nm fitted for a flat-topped
and for a Gaussian beam profile, and beta at 248 nm "taken from Ref. 3" (Taylor et
al. 1988), which is reproduced. Each beta cell prints its own power of ten; it is
stored in e-notation (``1.9e-10``, scale 1).

Checked against the table image (tbl_27e00c033f1a.jpg): the two header rows (the
"beta at 193 nm" span over "Flat Topped" and "Gaussian") and every cell agree. The
image prints "CaF2" and "BaF2" without subscripts, as extracted. The "\leq" of MgF2 at
248 nm is Taylor's upper bound; it is kept (the record excludes it).
"""

import re

from _json_table import json_table
from _library import MINUS, Point, write

DIRECTORY = "Kittelmann_1994_IntensityDependentTransmission"
TABLE = "tbl_27e00c033f1a"

MATERIALS = {"MgF_2": "MgF2", "CaF2": "CaF2", "BaF2": "BaF2", "Fused silica": "SiO2"}

_CELL = re.compile(
    r"^(?P<bound>\\leq)?\s*(?P<mantissa>[\d.]+)\s*\\times\s*10\^\{(?P<power>[-\d]+)\}$"
)


def scientific(cell: str) -> str:
    r"""Return a cell ``[\leq] m \times 10^{k}`` as ``mek`` (bound in printed)."""
    match = _CELL.match(cell.strip().replace(MINUS, "-"))
    if match is None:
        msg = f"not a power-of-ten cell: {cell!r}"
        raise ValueError(msg)
    value = f"{match['mantissa']}e{match['power']}"
    float(value)
    return value


def extract() -> dict[str, list[Point]]:
    """Return the band gaps, the two 193 nm fits and Taylor's 248 nm betas."""
    grid = json_table(DIRECTORY, TABLE)
    if grid[1][3:5] != ["Flat Topped", "Gaussian"]:
        msg = f"unexpected sub-header {grid[1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        material = MATERIALS[cells[0]]
        printed = " | ".join(cells)
        stem = f"{material}/kittelmann_1994__"
        files[f"{stem}band_gap_{material}.csv"] = [
            Point(material, "", cells[2], "", printed)
        ]
        files[f"{stem}beta_193_flat_top_{material}.csv"] = [
            Point("193", "193", scientific(cells[3]), "", printed)
        ]
        files[f"{stem}beta_193_gaussian_{material}.csv"] = [
            Point("193", "193", scientific(cells[4]), "", printed)
        ]
        files[f"{stem}taylor_248_{material}.csv"] = [
            Point("248", "248", scientific(cells[5]), "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
