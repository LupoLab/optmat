# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract White, Smith & Milam (1984): gamma of fused silica and BK-10 (Table 1).

Table 1 has one row per glass and two columns, gamma at 1064 nm ("Previous Work") and at
355 nm ("This Work"), in units of 10^-16 cm^2/W (footnote a). The wavelength of each
column is in its header, not in a cell. The cells are LaTeX, "2.7 \pm 0.5^{b}": the
value, the uncertainty and a footnote mark (b: Ref. 2; c: Ref. 11), which the record
explains. The table image (``tables/tbl_985eae28abf9.jpg``) matches the cells.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "White_1984_DirectMeasurementNonlinear"
TABLE = "tbl_985eae28abf9"
ROWS = {"Fused Silica": "SiO2", "BK-10": "BK-10"}
COLUMNS = {1: ("1064", "previous_1064"), 2: ("355", "this_work_355")}
_CELL = re.compile(r"^\s*([\d.]+)\s*\\pm\s*([\d.]+)\s*(\^\{[a-z]\})?\s*$")


def extract() -> dict[str, list[Point]]:
    """Return one file per glass and column."""
    files = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        material = ROWS[cells[0].strip()]
        for column, (wavelength, name) in COLUMNS.items():
            match = _CELL.match(cells[column])
            if match is None:
                msg = f"unexpected cell {cells[column]!r}"
                raise ValueError(msg)
            point = Point(
                "1",
                wavelength,
                number(match.group(1)),
                number(match.group(2)),
                " | ".join(cells),
            )
            files[f"{material}/white_1984__{name}_{material}.csv"] = [point]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
