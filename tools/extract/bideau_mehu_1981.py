# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Bideau-Mehu et al. (1981): Ne, Ar, Kr, Xe n - 1 at 140-254 nm (Tables 2-5).

The library's table captions are wrong for three of the four tables (checked on the
table images): ``tbl_804a29eb8dd3`` (uncaptioned) is Table 2, neon;
``tbl_724798d9c01d`` is Table 3, argon; ``tbl_7062fe0f03a9`` (captioned "Table 5
Xenon") is Table 4, krypton, as its Ref. 22 (Smith et al., krypton) column and its
values show; ``tbl_7ab7e916e7b3`` (uncaptioned) is Table 5, xenon, whose values run
through the 146.9 nm resonance.

Each table's first column is the paper's measurement and its last the paper's own
Sellmeier fit. The columns between reproduce other authors' values; the printed tables
mark those that are *extrapolations* of the other authors' formulas with a raised dot,
which the extracted CSVs lose. The marks, read from the table images, are listed in
``EXTRAPOLATED``; marked cells are not extracted, because an extrapolated formula value
is not a measurement.
"""

from _library import Point, number, table, write

DIRECTORY = "BideauMehu_1981_MeasurementRefractiveIndices"
TABLES = {
    "Ne": "tbl_804a29eb8dd3",
    "Ar": "tbl_724798d9c01d",
    "Kr": "tbl_7062fe0f03a9",
    "Xe": "tbl_7ab7e916e7b3",
}

# gas -> column index -> dataset id, for every column after the wavelength.
COLUMNS = {
    "Ne": {
        1: "experiment_Ne",
        2: "cuthbertson_1932_Ne",
        3: "dalgarno_kingston_1960_Ne",
        4: "langhoff_karplus_1969_Ne",
        5: "kaveeshwar_1968_Ne",
        6: "fit_Ne",
    },
    "Ar": {
        1: "experiment_Ar",
        2: "chashchina_1968_Ar",
        3: "larsen_1934_Ar",
        4: "heddle_1963_Ar",
        5: "peck_fisher_1964_Ar",
        6: "kingston_1964_Ar",
        7: "langhoff_karplus_1969_Ar",
        8: "fit_Ar",
    },
    "Kr": {
        1: "experiment_Kr",
        2: "smith_1975_Kr",
        3: "kronjager_1936_Kr",
        4: "koch_1949_Kr",
        5: "chashchina_shreider_1969_Kr",
        6: "langhoff_karplus_1969_Kr",
        7: "kingston_1964_Kr",
        8: "fit_Kr",
    },
    "Xe": {
        1: "experiment_Xe",
        2: "kronjager_1936_Xe",
        3: "koch_1949_Xe",
        4: "chashchina_shreider_1969_Xe",
        5: "langhoff_karplus_1969_Xe",
        6: "kingston_1964_Xe",
        7: "fit_Xe",
    },
}

# (gas, column index) -> wavelengths (as printed) whose cell carries the raised dot of
# an extrapolated value, read from the table images.
EXTRAPOLATED = {
    ("Ne", 2): "all",
    ("Ar", 3): {
        "228.8",
        "213.9",
        "194.2",
        "184.9",
        "180.7",
        "170.2",
        "164.1",
        "152.5",
        "140.4",
    },
    ("Ar", 5): "all",
    ("Kr", 3): {"213.9", "206.2", "184.4", "170.2", "164.1", "152.5", "140.4"},
    ("Kr", 4): "all",
    ("Xe", 2): {
        "213.9",
        "194.2",
        "184.9",
        "180.7",
        "170.2",
        "164.1",
        "152.5",
        "149.5",
        "148.3",
        "140.4",
    },
    ("Xe", 3): {
        "228.8",
        "213.9",
        "194.2",
        "184.9",
        "180.7",
        "170.2",
        "164.1",
        "152.5",
        "149.5",
        "148.3",
        "140.4",
    },
    # Not dotted, but the text (p. 395) says Chashchina and Shreider "give values down
    # to 135 nm but they stopped their measurements at 155 nm": below 155 nm their
    # xenon values are extrapolations of their formula.
    ("Xe", 4): {"152.5", "149.5", "148.3", "140.4"},
}


def _extrapolated(gas: str, column: int, wavelength: str) -> bool:
    marks = EXTRAPOLATED.get((gas, column), set())
    return marks == "all" or wavelength in marks


def extract() -> dict[str, list[Point]]:
    """Return one data file per column and gas, without extrapolated or empty cells."""
    files: dict[str, list[Point]] = {}
    for gas, table_id in TABLES.items():
        rows = table(DIRECTORY, table_id)[1:]
        for column, dataset in COLUMNS[gas].items():
            points = []
            for cells in rows:
                wavelength, cell = cells[0].strip(), cells[column].strip()
                if cell in {"", "-"} or _extrapolated(gas, column, wavelength):
                    continue
                points.append(
                    Point(wavelength, wavelength, number(cell), "", " | ".join(cells))
                )
            if points:
                files[f"{gas}/bideau_mehu_1981__{dataset}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
