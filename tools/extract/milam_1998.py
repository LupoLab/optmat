# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Milam (1998): the measured n2 (gamma) of fused silica, Tables 1 and 2.

Table 1 lists 24 absolute measurements, each with its reference number; every
measurement is another paper's value, so the rows are grouped into one reproduced
dataset per cited reference. Table 2 is the author's assembly of eleven values used for
the weighted average "at 1053 nm", three of them (804, 1319 and 1550 nm) translated to
1064 nm with the PERT equation.

Structural fix, read from the table images (``tables/tbl_ee3c0e192f6c.jpg`` and
``tables/tbl_cda7eab4b25e.jpg``): the extraction renders the printed "±" as "6" (the
paper's font maps ± to the glyph of 6, as its text maps - to 2 and ( ) to ~ !), so
"5.6 6 0.8" is "5.6 ± 0.8". Every cell of the first column has exactly this form.

Table 1's misprints (row 3's ±0.17 for ±1.7; row 17's 532 nm for 1064 nm) are errata in
the record, not fixed here.
"""

from _library import Point, number, table, write

DIRECTORY = "Milam_1998_ReviewAssessmentMeasured"
TABLE_1 = "tbl_ee3c0e192f6c"
TABLE_2 = "tbl_cda7eab4b25e"

# Reference number (Table 1's third column) -> dataset id, as the record names them.
REFERENCES = {
    "4": "desalvo_1996",
    "5": "taylor_1996",
    "6": "shimada_1995",
    "7": "kim_1994",
    "9": "ross_1990",
    "10": "kato_1995",
    "11": "kim_hutchinson_1989",
    "12": "vasileva_1985",
    "13": "white_1984",
    "14": "williams_1985",
    "16": "altshuler_1977",
    "17": "milam_weber_1976",
    "18": "smith_1975",
    "19": "stolen_ashkin_1973",
    "21": "veduta_1968",
    "22": "perry",
}


def _value_sigma(cell: str) -> tuple[str, str]:
    """Split "value ± sigma", extracted as "value 6 sigma" (module docstring)."""
    value, plus_minus, sigma = cell.split()
    if plus_minus != "6":
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return one file per cited reference of Table 1, and Table 2's assembly."""
    files: dict[str, list[Point]] = {
        f"SiO2/milam_1998__{name}.csv": [] for name in REFERENCES.values()
    }
    for index, cells in enumerate(table(DIRECTORY, TABLE_1)[1:], start=1):
        value, sigma = _value_sigma(cells[0])
        name = REFERENCES[cells[2].strip()]
        files[f"SiO2/milam_1998__{name}.csv"].append(
            Point(str(index), number(cells[1]), value, sigma, " | ".join(cells))
        )
    # Table 2: the caption and the text (p. 548) call these "eleven 1064-nm values"
    # for the weighted average; the printed wavelength is where each was measured,
    # kept as the column measured_nm, and x is the 1064 nm the values stand for.
    assembly = []
    for index, cells in enumerate(table(DIRECTORY, TABLE_2)[1:], start=1):
        value, sigma = _value_sigma(cells[0])
        assembly.append(
            Point(
                str(index),
                "1064",
                value,
                sigma,
                " | ".join(cells),
                {"measured_nm": number(cells[1])},
            )
        )
    files["SiO2/milam_1998__table2_assembly.csv"] = assembly
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
