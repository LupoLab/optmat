# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Liu, Yen & Bloembergen (1979): two-photon absorption at 355 and 266 nm.

Table I (``tbl_14fb6f23a7e6``, p. 1017) gives per material the band gap ("Band gaps are
listed when available", footnote a) and beta in cm/MW at 355 and 266 nm. The
extraction loses the superscripts of the powers of ten ("<6.2 x 10-6" is the image's
"<6.2 x 10^-6"), which the script reads as exponents. beta cells printed "<" are upper
bounds: kept (the number in ``value``, the raw cell in ``printed``), and the record
marks them excluded. The relative uncertainty ("+- 15%") stays in ``printed``; the
footnote mark "b" of fused silica is stripped. Band gaps printed "~6" and "~11" (the
image shows the tilde) are approximate: the number is kept and the tilde stays in
``printed``; SrF2's gap is printed "?" and is not extracted.

The table image (``tables/tbl_14fb6f23a7e6.jpg``) matches every extracted cell.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Liu_1979_TwoPhotonAbsorption"
TABLE = "tbl_14fb6f23a7e6"
MATERIALS = {
    "MgF_2": "MgF2",
    "SrF_2": "SrF2",
    "BaF_2": "BaF2",
    "CdF_2": "CdF2",
    "LiYF_4": "LiYF4",
    "Fused silica": "SiO2",
}
COLUMNS = {2: "355", 3: "266"}

TIMES = "\N{MULTIPLICATION SIGN}"
PLUS_MINUS = "\N{PLUS-MINUS SIGN}"
_BETA = re.compile(
    r"^(?P<bound><)?(?P<mantissa>[\d.]+)\s*" + TIMES + r"\s*10(?P<power>-\d+)"
    r"(?:\s*" + PLUS_MINUS + r"\s*(?P<percent>\d+)%)?(?:b)?$"
)


def _beta(cell: str) -> str:
    """Return a beta cell ``[<] m x 10-k [+- p%]`` as ``me-k`` (the rest is printed)."""
    match = _BETA.match(cell.strip())
    if match is None:
        msg = f"unexpected beta cell {cell!r}"
        raise ValueError(msg)
    value = f"{match['mantissa']}e{match['power']}"
    float(value)
    return value


def extract() -> dict[str, list[Point]]:
    """Return beta at each wavelength and the band gap, per material."""
    rows = table(DIRECTORY, TABLE)
    header, body = rows[0], rows[1:]
    if header[2] != "β(cm/MW)355 nm" or header[3] != "β(cm/MW)266 nm":
        msg = f"unexpected header {header!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in body:
        material = MATERIALS[cells[0]]
        printed = " | ".join(cells)
        for column, wavelength in COLUMNS.items():
            if cells[column]:
                files[f"{material}/liu_1979__beta_{wavelength}_{material}.csv"] = [
                    Point("1", wavelength, _beta(cells[column]), "", printed)
                ]
        if cells[1] != "?":
            files[f"{material}/liu_1979__band_gap_{material}.csv"] = [
                Point("1", "", number(cells[1].removeprefix("~")), "", printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
