# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rakic (1995): n and k of aluminium, 6.2 meV to 10 keV (Table 2).

Table 2 ("Values of the Intrinsic Optical Constants n(omega) and k(omega) and
Reflectance R(omega) of Aluminum Determined in this Study") prints hbar omega (eV),
lambda (um), n, k and R in scientific notation on four pages: ``tbl_a8fb03d43d0c`` (p.
4762, 6.1993E-03 to 2.4000E+00 eV, 62 rows), ``tbl_956d16c6a12e`` (p. 4763, 2.6 to 73.6
eV, 64 rows), ``tbl_7415e4ccca18`` (p. 4764, 73.7 eV to 1.58 keV, 63 rows) and
``tbl_784e13d5a0b4`` (p. 4765, 1.595 to 10 keV, 17 rows): 206 rows. The page order is
the images' (each page repeats the header; the energies continue); the script checks it
by the wavelength column, which falls in every row.

STRUCTURAL FIX. The journal's PDF text layer encodes the parentheses and the signs of
the exponents with the digits 1 and 2 (the header reads '"v 1eV2' for "hbar omega
(eV)"), and the library's cells keep that encoding: "6.1993E203" is the printed
"6.1993E-03" and "2.0000E102" the printed "2.0000E+02" (checked on all four images). The
script decodes every cell, "E2" to "E-" and "E1" to "E+", and keeps the decoded strings
(the page's notation, with an ASCII minus) as x, value and ``printed``: x and value are
the numbers as the page prints them. Two energies print a one-digit exponent
("2.1697E-2", "4.0000E-2"), as the page does.

The row id is the row's position in the table (1-206), because the energy column is
not unique as printed (see the errata: row 74, 10 eV, prints "1.0000E+00").

CHECKS (2026-09-28): every cell was compared with the images. R has no quantity and
stays in ``printed``; the printed R agrees with R = ((n - 1)^2 + k^2)/((n + 1)^2 + k^2)
from the same row to 4.4e-4 (relative) in every row up to 66 eV (the worst at 62 eV,
where R = 2.3e-5); above 68 eV, where R < 4e-4 and depends on 1 - n, it agrees within
the rounding of the printed n and k except at 700 eV (3e-5 outside), 6 keV (4 %) and 10
keV (2 %), and at 1.1 keV, where n is misprinted. hbar omega times lambda is 1.23984
eV um to 5e-5 in every row but four misprinted energies (errata) and 1.9 eV, whose
wavelength is misprinted (below).

The misprints the checks found, all confirmed on the images, are corrected by errata in
the record, not here: the energies of rows 74 (10 eV printed "1.0000E+00"), 177, 179
and 187 (1450, 1520 and 1562 eV printed with "E-03"), and n at 1.1 keV (row 173,
"9.9961E-05" for 9.9961E-01). The wavelength at 1.9 eV (row 59) is printed
"6.5225E-01" where 1.23984/1.9 = 0.65255 um (the 19 and 190 eV rows print 6.5255E-02 and
6.5255E-03); lambda is not a stored column, so this one is only noted.
"""

import re

from _library import Point, table, write

DIRECTORY = "Rakic_1995_AlgorithmDeterminationIntrinsic"
PARTS = ("tbl_a8fb03d43d0c", "tbl_956d16c6a12e", "tbl_7415e4ccca18", "tbl_784e13d5a0b4")
HEADER = ['"v 1eV2', "l 1µm2", "n", "k", "R"]
ROWS = 206
# A library cell: mantissa, the encoded exponent sign (1 for +, 2 for -), exponent.
_ENCODED = re.compile(r"^(\d\.\d+)E([12])(\d+)$")


def decode(cell: str) -> str:
    """Return a cell as the page prints it: "6.1993E203" becomes "6.1993E-03"."""
    match = _ENCODED.match(cell.strip())
    if match is None:
        msg = f"not an encoded number: {cell!r}"
        raise ValueError(msg)
    mantissa, sign, exponent = match.groups()
    return f"{mantissa}E{'-' if sign == '2' else '+'}{exponent}"


def extract() -> dict[str, list[Point]]:
    """Return aluminium's n and k against the photon energy."""
    files: dict[str, list[Point]] = {
        "Al/rakic_1995__n_Al.csv": [],
        "Al/rakic_1995__k_Al.csv": [],
    }
    previous = float("inf")
    row = 0
    for table_id in PARTS:
        rows = table(DIRECTORY, table_id)
        if rows[0] != HEADER:
            msg = f"{table_id}: unexpected header {rows[0]}"
            raise ValueError(msg)
        for cells in rows[1:]:
            energy, wavelength, n, k, reflectance = (decode(cell) for cell in cells)
            if float(wavelength) >= previous:
                msg = f"{table_id}: wavelengths out of order at {cells}"
                raise ValueError(msg)
            previous = float(wavelength)
            row += 1
            printed = " | ".join((energy, wavelength, n, k, reflectance))
            for name, value in (("n", n), ("k", k)):
                files[f"Al/rakic_1995__{name}_Al.csv"].append(
                    Point(str(row), energy, value, "", printed)
                )
    if row != ROWS:
        msg = f"Table 2 has {row} rows, not {ROWS}"
        raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
