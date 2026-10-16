# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hajigeorgiou (2013), J. Chem. Phys. 138, 014309: 16O2, X 3Sigma_g-.

All tables checked against their images.

* Table I (``tbl_4bb1a4a3289b``, CSV): the G_v, B_v, D_v (x 1e-6) and H_v (x 1e-12)
  "employed in calculating term values", v = 0-31 (the N_max column stays in
  ``printed``). One file per constant; rows are v.
* Table III (``tbl_a6413e37f185``): the quantum-mechanical G_v, B_v, D_v and H_v of the
  fitted potential, v = 0-38. **Structural fix:** the library's cells hold each column
  as one string of its characters in glyph order (column-major by position on the
  page), not as rows. The table is therefore transcribed below from the PDF's own text
  layer (``pdftotext -layout``, page 6 of ``original/source_01``), checked row by row
  against the table image; the script checks the transcription against the library
  cells (each column's characters, spaces dropped, as a multiset) and the G_v column
  against Table V's cells (identical). D_v and H_v are printed as "m x 10^-k" with
  exponents that change down the column; they are written as "me-k". Table III's G_v
  is not a separate file: Table V holds the same numbers.
* Table V (``tbl_02f68c84e070``, CSV): G_v, v = 0-38 (the turning points r_inner and
  r_outer stay in ``printed``: no quantity).
* Table IV (``tbl_9bf6778409d5``, JSON cells): G_v - G_0 and B_v for v = 29-38 of Yang
  and Wodtke (YW) and of Bytautas et al. (Ref. 7, in parentheses beside the present
  work's values, which repeat Table III rounded and are not taken). The YW note marks
  b, c, d are dropped from the numbers: b "Experimental value with uncertainty +-40
  cm^-1" gives the sigma column (40); d "Extrapolated value with unknown uncertainty"
  leaves it empty; c "Value calculated from hybrid potential with uncertainty +-0.01
  cm^-1" is the B dataset's absolute uncertainty in the record.
* Table II (the potential parameters) is a formula in the record, not a data file.
"""

import re
from collections import Counter

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Hajigeorgiou_2013_PotentialEnergyFunction"
TABLE_1 = "tbl_4bb1a4a3289b"
TABLE_3 = "tbl_a6413e37f185"
TABLE_4 = "tbl_9bf6778409d5"
TABLE_5 = "tbl_02f68c84e070"
PREFIX = "O2/hajigeorgiou_2013__"
UPSILON = "\N{GREEK SMALL LETTER UPSILON}"
TIMES = "\N{MULTIPLICATION SIGN}"

# Table III as printed (the PDF's text layer, checked against the image): v, G_v, B_v,
# D_v, H_v; the digit-group spaces kept, the printed minus sign (U+2212) and
# multiplication sign written as "-" and "x".
TABLE_3_PRINTED = (
    ("0", "787.3223", "1.437 674 521", "4.839 482 x 10-6", "-4.220 376 x 10-14"),
    ("1", "2343.7129", "1.421 859 938", "4.837 533 x 10-6", "-6.110 061 x 10-14"),
    ("2", "3876.6372", "1.406 114 625", "4.837 444 x 10-6", "-1.310 273 x 10-13"),
    ("3", "5386.1599", "1.390 428 695", "4.837 576 x 10-6", "-3.438 115 x 10-13"),
    ("4", "6872.3983", "1.374 786 275", "4.836 856 x 10-6", "-6.889 675 x 10-13"),
    ("5", "8335.5004", "1.359 173 600", "4.835 212 x 10-6", "-1.092 499 x 10-12"),
    ("6", "9775.6206", "1.343 581 809", "4.833 536 x 10-6", "-1.485 936 x 10-12"),
    ("7", "11 192.9010", "1.328 006 430", "4.833 271 x 10-6", "-1.803 383 x 10-12"),
    ("8", "12 587.4579", "1.312 446 017", "4.836 187 x 10-6", "-2.006 540 x 10-12"),
    ("9", "13 959.3743", "1.296 900 135", "4.844 118 x 10-6", "-2.092 867 x 10-12"),
    ("10", "15 308.6965", "1.281 367 233", "4.858 815 x 10-6", "-2.081 087 x 10-12"),
    ("11", "16 635.4312", "1.265 842 907", "4.881 828 x 10-6", "-2.027 810 x 10-12"),
    ("12", "17 939.5457", "1.250 318 462", "4.914 388 x 10-6", "-1.995 512 x 10-12"),
    ("13", "19 220.9671", "1.234 780 076", "4.957 515 x 10-6", "-2.066 210 x 10-12"),
    ("14", "20 479.5826", "1.219 208 311", "5.011 832 x 10-6", "-2.346 373 x 10-12"),
    ("15", "21 715.2399", "1.203 577 827", "5.077 664 x 10-6", "-2.947 441 x 10-12"),
    ("16", "22 927.7465", "1.187 857 180", "5.155 102 x 10-6", "-3.997 457 x 10-12"),
    ("17", "24 116.8688", "1.172 008 537", "5.244 020 x 10-6", "-5.658 037 x 10-12"),
    ("18", "25 282.3274", "1.155 987 038", "5.344 206 x 10-6", "-8.093 603 x 10-12"),
    ("19", "26 423.7907", "1.139 740 021", "5.455 480 x 10-6", "-1.149 868 x 10-11"),
    ("20", "27 540.8630", "1.123 205 847", "5.577 855 x 10-6", "-1.605 080 x 10-11"),
    ("21", "28 633.0686", "1.106 312 739", "5.711 620 x 10-6", "-2.192 116 x 10-11"),
    ("22", "29 699.8298", "1.088 977 478", "5.857 560 x 10-6", "-2.922 764 x 10-11"),
    ("23", "30 740.4395", "1.071 104 214", "6.017 021 x 10-6", "-3.802 036 x 10-11"),
    ("24", "31 754.0273", "1.052 583 444", "6.191 946 x 10-6", "-4.825 588 x 10-11"),
    ("25", "32 739.5196", "1.033 291 086", "6.384884 x 10-6", "-5.972 574 x 10-11"),
    ("26", "33 695.5917", "1.013 087 819", "6.598 778 x 10-6", "-7.201 281 x 10-11"),
    ("27", "34 620.6125", "0.991 818 644", "6.836 596 x 10-6", "-8.447 021 x 10-11"),
    ("28", "35 512.5781", "0.969 312 472", "7.100 642 x 10-6", "-9.611 354 x 10-11"),
    ("29", "36 369.0317", "0.945 382 084", "7.391 422 x 10-6", "-1.055 474 x 10-10"),
    ("30", "37 186.9666", "0.919 824 531", "7.706 355 x 10-6", "-1.109 109 x 10-10"),
    ("31", "37 962.7048", "0.892 422 226", "8.038 158 x 10-6", "-1.100 649 x 10-10"),
    ("32", "38 691.7440", "0.862 944 204", "8.374 440 x 10-6", "-1.012 869 x 10-10"),
    ("33", "39 368.5673", "0.831 144 654", "8.701 664 x 10-6", "-8.559 929 x 10-11"),
    ("34", "39 986.4128", "0.796 747 157", "9.022 935 x 10-6", "-7.290 530 x 10-11"),
    ("35", "40 537.0281", "0.759 376 380", "9.415 421 x 10-6", "-9.768 058 x 10-11"),
    ("36", "41 010.5181", "0.718 317 371", "1.020 175 x 10-5", "-2.645 162 x 10-10"),
    ("37", "41 395.6545", "0.671 740 031", "1.249 039 x 10-5", "-9.419 619 x 10-10"),
    ("38", "41 681.7298", "0.614 420 969", "2.025 242 x 10-5", "-4.162 905 x 10-9"),
)
_POWER = re.compile(r"^(-?[\d. ]+?) x 10-(\d+)$")


def _characters(text: str) -> Counter[str]:
    """Return the characters of a cell, spaces dropped, minus and times as ASCII."""
    ascii_text = text.replace("\N{MINUS SIGN}", "-").replace(
        "\N{MULTIPLICATION SIGN}", "x"
    )
    return Counter(ascii_text.replace(" ", ""))


def _check_table_3() -> None:
    """Check the transcription against the library's (column-scrambled) cells."""
    cells = json_table(DIRECTORY, TABLE_3)
    if cells[0] != [UPSILON] + [f"{c}{UPSILON}" for c in "GBDH"] or len(cells) != 2:
        msg = f"unexpected Table III layout {cells[0]!r}, {len(cells)} rows"
        raise ValueError(msg)
    for column, scrambled in enumerate(cells[1]):
        printed = "".join(row[column] for row in TABLE_3_PRINTED)
        if _characters(printed) != _characters(scrambled):
            msg = f"Table III column {cells[0][column]!r} does not match the cells"
            raise ValueError(msg)
    table_5 = [(row[0], row[1]) for row in table(DIRECTORY, TABLE_5)[1:]]
    printed_g = [(row[0], row[1].replace(" ", "")) for row in TABLE_3_PRINTED]
    if printed_g != table_5:
        msg = "Table III's G_v differs from Table V's"
        raise ValueError(msg)


def _power(cell: str) -> str:
    """Return "m x 10-k" as "me-k" (the digit-group spaces dropped)."""
    match = _POWER.match(cell)
    if match is None:
        msg = f"not a power of ten: {cell!r}"
        raise ValueError(msg)
    mantissa = number(match.group(1).replace(" ", ""))
    return f"{mantissa}e-{match.group(2)}"


def _table1() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_1)
    header = [UPSILON, "Nmax", f"G{UPSILON}", f"B{UPSILON}"]
    header += [f"D{UPSILON} {TIMES} 10-6", f"H{UPSILON} {TIMES} 10-12"]
    if grid[0] != header:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        v = str(int(cells[0]))
        printed = " | ".join(cells)
        for column, name in ((2, "G"), (3, "B"), (4, "D"), (5, "H")):
            point = Point(v, v, number(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}table1_{name}_O2.csv", []).append(point)
    return files


def _table3() -> dict[str, list[Point]]:
    _check_table_3()
    files: dict[str, list[Point]] = {}
    for row in TABLE_3_PRINTED:
        v = row[0]
        printed = " | ".join(row)
        values = {
            "B": number(row[2].replace(" ", "")),
            "D": _power(row[3]),
            "H": _power(row[4]),
        }
        for name, value in values.items():
            point = Point(v, v, value, "", printed)
            files.setdefault(f"{PREFIX}table3_{name}_O2.csv", []).append(point)
    return files


def _table5() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_5)
    if grid[0][:2] != [UPSILON, f"G{UPSILON} (cm-1)"]:
        msg = f"unexpected Table V header {grid[0]!r}"
        raise ValueError(msg)
    points = [
        Point(cells[0], cells[0], number(cells[1]), "", " | ".join(cells))
        for cells in grid[1:]
    ]
    return {f"{PREFIX}table5_G_O2.csv": points}


_MARK = re.compile(r"^([\d.]+)\^\{([bcd])\}$")
_PAIR = re.compile(r"^([\d.]+) \(([\d.]+)\)$")


def _table4() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_4)
    if grid[1] != ["", "YW^{a}", "PW^{a} (Ref. 7)", "YW^{a}", "PW^{a} (Ref. 7)"]:
        msg = f"unexpected Table IV header {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        v = cells[0]
        printed = " | ".join(cells)
        for column, name in ((1, "G"), (3, "B")):
            match = _MARK.match(cells[column])
            if match is None:
                msg = f"unexpected Table IV cell {cells[column]!r}"
                raise ValueError(msg)
            # Note b: "Experimental value with uncertainty +-40 cm^-1".
            sigma = "40" if match.group(2) == "b" else ""
            point = Point(v, v, number(match.group(1)), sigma, printed)
            files.setdefault(f"{PREFIX}table4_{name}_yang_wodtke_O2.csv", []).append(
                point
            )
            pair = _PAIR.match(cells[column + 1])
            if pair is None:
                msg = f"unexpected Table IV cell {cells[column + 1]!r}"
                raise ValueError(msg)
            point = Point(v, v, number(pair.group(2)), "", printed)
            files.setdefault(f"{PREFIX}table4_{name}_bytautas_O2.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the files of Tables I, III, IV and V."""
    return _table1() | _table3() | _table4() | _table5()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
