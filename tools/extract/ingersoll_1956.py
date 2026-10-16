# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ingersoll and Liebenberg (1956): Verdet constants of 13 curated gases.

Tables I (noble gases, ``tbl_67fd0411dc50``), II (hydrocarbons, ``tbl_321e7af15c60``)
and III (common gases, ``tbl_de2fd916f64f``), images checked, each print in an upper
part the Verdet constant (microminutes per oersted-centimeter, gas at 76 cm and 0 C)
at 4360, 5460, 5780 and 5893 A, and in a lower part the Faraday dispersion at 17
wavelengths, 3635-9875 A, 'with all values reduced to 1.000 at 5780 A'. Both parts are
extracted per gas: ``absolute_<gas>`` (the four upper rows) and ``dispersion_<gas>``
(the relative values, which the source record scales by the gas's printed V(5780)).
The upper part is the first four data rows. Helium has no value at 9875 A (``...``),
skipped. Table I's 1/lambda^2 column and the non-curated gases (hydrocarbons but
methane, freon-12) are not extracted.
"""

from _library import Point, number, table, write

DIRECTORY = "Ingersoll_1956_FaradayEffectGases"
# Table -> {column: gas}.
TABLES = {
    "tbl_67fd0411dc50": {2: "He", 3: "Ne", 4: "Ar", 5: "Kr", 6: "Xe"},
    "tbl_321e7af15c60": {1: "CH4"},
    "tbl_de2fd916f64f": {
        1: "H2",
        2: "D2",
        3: "N2",
        4: "O2",
        5: "N2O",
        6: "CO2",
        7: "CO",
        8: "SF6",
    },
}
UPPER_ROWS = 4


def extract() -> dict[str, list[Point]]:
    """Return the absolute and relative Verdet constants of Tables I-III."""
    files: dict[str, list[Point]] = {}
    for table_id, gases in TABLES.items():
        rows = table(DIRECTORY, table_id)[1:]
        for index, cells in enumerate(rows):
            wavelength = number(cells[0])
            part = "absolute" if index < UPPER_ROWS else "dispersion"
            for column, gas in gases.items():
                if cells[column].strip() == "...":
                    continue
                files.setdefault(f"{gas}/ingersoll_1956__{part}_{gas}.csv", []).append(
                    Point(
                        wavelength,
                        wavelength,
                        number(cells[column]),
                        "",
                        " | ".join(cells),
                    )
                )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
