# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Taylor, Gibson & Roberts (1988): two-photon absorption at 248 nm (Table 1).

Table 1 gives, per material, the band gap, the two-photon absorption coefficient
beta measured at 248 nm (0.7 ps) and beta at 266 nm from Liu et al. (the two Liu
papers of footnote b), which is reproduced. Each beta cell prints its own power of
ten ("4.5 \times 10^{-11}"); it is stored in e-notation (``4.5e-11``, scale 1).

Checked against the table image (tbl_d491af0a4326.jpg): every cell agrees. The
footnote marker c of LiF and MgF2 ("nonlinear absorption fitted as two-photon
absorption") is extracted as a trailing "c" in the material cell; the script strips
it. Cells printed "<" or "\leq" are upper bounds: they are kept (the bound's number
in ``value``, the raw cell in ``printed``) and the record marks them excluded.
"""

import re

from _library import MINUS, Point, table, write

DIRECTORY = "Taylor_1988_TwoPhotonAbsorption"
TABLE = "tbl_d491af0a4326"

MATERIALS = {
    "Fused silica": "SiO2",
    "BaF_{2}": "BaF2",
    "SrF_{2}": "SrF2",
    "CaF_{2}": "CaF2",
    # The trailing "c" is footnote c, lost as a superscript by the extraction.
    "LiFc": "LiF",
    "MgF_{2}c": "MgF2",
}

_CELL = re.compile(
    r"^(?P<bound><|\\leq)?\s*(?P<mantissa>[\d.]+)\s*\\times\s*10\^\{(?P<power>[-\d]+)\}$"
)


def scientific(cell: str) -> str:
    r"""Return a cell ``[<|\leq] m \times 10^{k}`` as ``mek`` (bound in printed)."""
    match = _CELL.match(cell.strip().replace(MINUS, "-"))
    if match is None:
        msg = f"not a power-of-ten cell: {cell!r}"
        raise ValueError(msg)
    value = f"{match['mantissa']}e{match['power']}"
    float(value)
    return value


def extract() -> dict[str, list[Point]]:
    """Return the band gaps, the 248 nm betas and Liu et al.'s 266 nm betas."""
    rows = table(DIRECTORY, TABLE)
    header, body = rows[0], rows[1:]
    if not header[2].startswith(r"\beta (cm/W) at 2h\nu = 10"):
        msg = f"unexpected header {header!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in body:
        material = MATERIALS[cells[0]]
        printed = " | ".join(cells)
        files[f"{material}/taylor_1988__band_gap_{material}.csv"] = [
            Point(material, "", cells[1], "", printed)
        ]
        files[f"{material}/taylor_1988__beta_248_{material}.csv"] = [
            Point("248", "248", scientific(cells[2]), "", printed)
        ]
        files[f"{material}/taylor_1988__liu_266_{material}.csv"] = [
            Point("266", "266", scientific(cells[3]), "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
