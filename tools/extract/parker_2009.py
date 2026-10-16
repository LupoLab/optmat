# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Parker, Armstrong, Boca & Taylor (2009): static-field ionisation rates, He.

Table 1 (``tbl_f4122d1f56cb``) gives Gamma (au) against the peak field E (au) from this
work (full-dimensional time-dependent integration), Scrinzi et al [3], Themelis et al
[4] and "ADK theory [9]"; one file per column. Each value is printed
``a \times 10^{-k}`` and written as ``ae-k``; a blank cell (no value at that field) is
skipped.

The table is scrambled in the library: its one data row holds, per column, the
characters of every row interleaved (e.g. the E column
``000...0...0000000000001111...``). STRUCTURAL FIX: the table is transcribed below
from the table image (``tables/tbl_f4122d1f56cb.jpg``), and the script checks that each
scrambled cell is exactly the multiset of the characters of the transcribed column
(spaces dropped, the minus sign extracted as "-"); a mismatch stops the extraction. For
the E, Scrinzi and ADK columns the scrambled cell is also exactly the transcription
interleaved by character position; for the other two the interleaving follows the
printed alignment instead, and only the multiset check applies.
"""

from collections import Counter

from _library import Point, table, write

DIRECTORY = "Parker_2009_UVStaticField"
TABLE = "tbl_f4122d1f56cb"
TIMES = "\N{MULTIPLICATION SIGN}"

# A rate as printed: (mantissa, k) for "mantissa x 10^-k"; None where the cell is empty.
Rate = tuple[str, str] | None
# Table 1 as printed (the image): E (au), this work, Scrinzi et al [3], Themelis et al
# [4], ADK theory [9].
TABLE_1_IMAGE: list[tuple[str, Rate, Rate, Rate, Rate]] = [
    ("0.045", ("9.4463", "15"), None, None, ("1.75", "14")),
    ("0.048", ("9.0006", "14"), None, None, ("1.61", "13")),
    ("0.05", ("3.3536", "13"), None, None, ("6.10", "13")),
    ("0.055", ("6.7817", "12"), None, None, ("1.11", "11")),
    ("0.06", ("7.8936", "11"), None, None, ("1.23", "10")),
    ("0.067", ("1.4045", "9"), None, ("1.6", "9"), ("1.96", "9")),
    ("0.068", ("1.8365", "9"), None, ("2.0", "9"), ("2.78", "9")),
    ("0.07", ("3.6360", "9"), None, ("3.6", "9"), ("5.42", "9")),
    ("0.075", ("1.6945", "8"), None, ("1.42", "8"), ("2.45", "8")),
    ("0.08", ("6.2093", "8"), ("4.63", "8"), ("5.57", "8"), ("9.15", "8")),
    ("0.09", ("5.3387", "7"), ("5.09", "7"), ("5.23", "7"), ("8.20", "7")),
    ("0.095", ("1.3413", "6"), None, ("1.30", "6"), ("2.06", "6")),
    ("0.10", ("2.9391", "6"), ("2.88", "6"), ("2.92", "6"), ("4.71", "6")),
    ("0.11", ("1.1755", "5"), ("1.15", "5"), ("1.17", "5"), ("1.96", "5")),
    ("0.12", ("3.6829", "5"), ("3.62", "5"), ("3.66", "5"), ("6.42", "5")),
    ("0.13", ("9.5748", "5"), ("9.43", "5"), ("9.50", "5"), ("1.74", "4")),
    ("0.14", ("2.1495", "4"), ("2.12", "4"), ("2.13", "4"), ("4.09", "4")),
    ("0.15", ("4.2913", "4"), ("4.23", "4"), ("4.25", "4"), ("8.56", "4")),
    ("0.16", ("7.7875", "4"), ("7.68", "4"), ("7.70", "4"), ("1.63", "3")),
    ("0.18", ("2.0578", "3"), ("2.03", "3"), ("2.03", "3"), ("4.73", "3")),
    ("0.20", ("4.3347", "3"), ("4.31", "3"), ("4.30", "3"), ("1.07", "2")),
    ("0.225", ("9.0115", "3"), None, None, ("2.57", "2")),
    ("0.25", ("1.5793", "2"), ("1.57", "2"), ("1.56", "2"), ("5.02", "2")),
    ("0.30", ("3.5857", "2"), ("3.56", "2"), ("3.52", "2"), ("1.35", "1")),
    ("0.34", ("5.7545", "2"), ("5.72", "2"), None, ("2.40", "1")),
    ("0.3773", ("8.1210", "2"), None, None, ("3.65", "1")),
    ("0.40", None, ("9.77", "2"), ("9.64", "2"), ("4.53", "1")),
]
COLUMNS = ("present", "scrinzi_1999", "themelis_1999", "adk")
HEADER = ["E (au)", "Present work", "Scrinzi et al [3]", "Themelis et al [4]"]


def _printed(rate: Rate) -> str:
    """Return a rate as the characters the print shows, "a", the times sign, "10-k"."""
    return "" if rate is None else f"{rate[0]}{TIMES}10-{rate[1]}"


def _check() -> None:
    """Check the transcription against the scrambled library cells."""
    grid = table(DIRECTORY, TABLE)
    if grid[1] != [*HEADER, "ADK theory [9]"]:
        msg = f"unexpected Table 1 header {grid[1]!r}"
        raise ValueError(msg)
    transcribed = ["".join(row[0] for row in TABLE_1_IMAGE)] + [
        "".join(_printed(row[1:][k]) for row in TABLE_1_IMAGE) for k in range(4)
    ]
    for index, (cell, text) in enumerate(zip(grid[2], transcribed, strict=True)):
        if Counter(cell.replace(" ", "")) != Counter(text):
            msg = f"Table 1 column {index}: {cell!r} does not scramble {text!r}"
            raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return the four rate columns of Table 1."""
    _check()
    files: dict[str, list[Point]] = {}
    for column, name in enumerate(COLUMNS):
        points = []
        for row in TABLE_1_IMAGE:
            rate = row[1:][column]
            if rate is None:
                continue
            printed = " | ".join([row[0], *(_printed(r) for r in row[1:])])
            value = f"{rate[0]}e-{rate[1]}"
            printed = f"{printed} (from the image)"
            points.append(Point(row[0], row[0], value, "", printed))
        files[f"He/parker_2009__table1_{name}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
