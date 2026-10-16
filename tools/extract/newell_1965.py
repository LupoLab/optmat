# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Newell & Baird (1965): microwave refractivities of gases at 47.7 GHz.

Table III (``tbl_601dafe76d7c``, cells file only): (n - 1) x 1e6 at 0 degrees C and 760
Torr, per gas (a heading row), for Newell and Baird's two resonators and the earlier
microwave values they compare. One file per gas and author (and cavity). Each value is
printed "value+-uncertainty"; Essen and Froome's CO2 value carries the footnote mark
"a" (corrected by the authors with Berthelot's law), stripped. The image agrees with
the cells ("Cylondrical" on the image's CO2 row is a misprint the cells correct).

Wavelengths: the table gives the frequency (nominal: 9, 24, 48, 72 GHz); the
wavelength column holds the vacuum wavelength c/f in m, computed here from the
paper's precise 47.736 GHz for Newell and Baird (p. 3751) and from the nominal
frequency for the others; Maryott and Buckley's means (rf to optical) have no
frequency and are stored as static ("inf"). The printed frequency stays in the
printed cells. Dry CO2-free air is not a curated species; the abstract's table
(``tbl_c32ef3abe0f0``) repeats the cylindrical-cavity values and is not extracted.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Newell_1965_AbsoluteDeterminationRefractive"
TABLE_III = "tbl_601dafe76d7c"
C = 299_792_458.0  # speed of light, m/s
GASES = {
    "Dry CO2-free air": None,
    "Argon": "Ar",
    "Carbon dioxide": "CO2",
    "Helium": "He",
    "Hydrogen": "H2",
    "Nitrogen": "N2",
    "Oxygen": "O2",
}
AUTHORS = {
    ("Newell-Baird", "Cylindrical"): ("newell_baird_cylindrical", 47.736e9),
    ("Newell-Baird", "Fabry-Perot"): ("newell_baird_fabry_perot", 47.736e9),
    ("Maryott-Buckley", ""): ("maryott_buckley_1953", None),
    ("Essen", ""): ("essen_1953", 9e9),
    ("Essen-Froome", ""): ("essen_froome_1951", 24e9),
    ("Froome", ""): ("froome_1955", 72e9),
}


def extract() -> dict[str, list[Point]]:
    """Return the refractivities per gas and author."""
    files: dict[str, list[Point]] = {}
    gas: str | None = None
    for cells in json_table(DIRECTORY, TABLE_III)[1:]:
        if not any(cells[1:]):
            gas = GASES[cells[0]]
            continue
        if gas is None:
            continue
        name, frequency = AUTHORS[cells[0], cells[1]]
        value, _, sigma = cells[2].partition("\N{PLUS-MINUS SIGN}")
        wavelength = "inf" if frequency is None else f"{C / frequency:.5g}"
        point = Point(
            cells[3] or "unstated",
            wavelength,
            number(value),
            sigma.strip().removesuffix("a"),
            " | ".join(cells),
        )
        path = f"{gas}/newell_1965__{name}_{gas}.csv"
        assert path not in files, path
        files[path] = [point]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
