# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract C. and M. Cuthbertson (1910, neon): doubled n - 1 of neon at three lines.

The unlabelled table ``tbl_9757a187853d`` (``lambda x 10^8``, ``(mu - 1) 10^6``) has
its three rows collapsed into one cell per column in the library
(``643854614800`` and ``134·02134·30134·63``). The table image shows three rows,
6438, 5461 and 4800 A with 134.02, 134.30 and 134.63; the script splits the cells into
four-digit wavelengths and ``ddd·dd`` values and checks that three of each come out.
The values are doubled refractivities, 2 (mu - 1) (the text doubles the 5461 A index
"for comparison with diatomic gases"); the raised decimal point is read as a point.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Cuthbertson_1910_RefractionDispersionNeon"
TABLE = "tbl_9757a187853d"


def extract() -> dict[str, list[Point]]:
    """Return the three neon values."""
    (cells,) = table(DIRECTORY, TABLE)[1:]
    wavelengths = re.findall(r"\d{4}", cells[0].strip())
    values = re.findall(r"\d{3}\N{MIDDLE DOT}\d{2}", cells[1].strip())
    if len(wavelengths) != 3 or len(values) != 3:
        msg = f"cannot split the collapsed row {cells!r}"
        raise ValueError(msg)
    printed = " | ".join(cells)
    points = [
        Point(w, w, number(v.replace("\N{MIDDLE DOT}", ".")), "", printed)
        for w, v in zip(wavelengths, values, strict=True)
    ]
    return {"Ne/cuthbertson_1910b__experiment_Ne.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
