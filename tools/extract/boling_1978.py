# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Boling, Glass & Owyoung (1978): n, n_d and n2 of Tables I and III.

Table I lists, per material, n, the n2 "measured" (from Owyoung's ellipse rotation)
and the n2 calculated with Eqs. (25)-(26); Table III lists n_d, the Abbe number, the n2
calculated with Eq. (38) and the n2 measured by time-resolved interferometry (Milam and
Weber), with its +-. The n2 are in esu of n = n0 + n2 <E^2> (nonlinear_index_field),
printed as "10^13 n2", i.e. in units of 1e-13 esu. The linear indices are Table I's n
(the wavelength is not stated) and Table III's n_d (the d line, 0.58756 um, Eq. 32).
One file per material and column.

The cells match the table images (``tables/tbl_1104b7cc486c.jpg``,
``tables/tbl_6d8d3f9d443b.jpg``); in Table III the image's "1:769" for the silicate
glass is a speck on the decimal point, extracted as 1.769. Table III's header prints
"x 10^23" for the calculated column (the image shows it): a misprint of the paper,
corrected by a scale erratum in the record, not here.
"""

from _library import Point, number, table, write

DIRECTORY = "Boling_1978_EmpiricalRelationshipsPredicting"
TABLE_I = "tbl_1104b7cc486c"
TABLE_III = "tbl_6d8d3f9d443b"
D_LINE = "0.58756"
NOT_STATED = "not stated"

# Row names as printed -> registered material ids.
MATERIALS = {
    "FK-6": "FK-6",
    "fused silica": "SiO2",
    "BK-7": "BK7",
    "LSO": "LSO",
    "ED-4": "ED-4",
    "EY-1": "EY-1",
    "SF-7": "SF-7",
    "LaK-3": "LaK-3",
    "YAG": "YAG",
    "LaSF-7": "LaSF-7",
    "silicate": "silicate_boling_1978",
    "FR-5": "FR-5",
    "ED-2": "ED-2",
    "LHG-5": "LHG-5",
    "LHG-6": "LHG-6",
    "EV-1": "EV-1",
    "FK-5": "FK-5",
    "fluoro-phosphate": "fluorophosphate_boling_1978",
    "FK-51": "FK-51",
}

# (table, column index) -> (dataset name, wavelength cell).
COLUMNS = {
    (TABLE_I, 1): ("n_table1", NOT_STATED),
    (TABLE_I, 2): ("n2_measured_table1", NOT_STATED),
    (TABLE_I, 3): ("n2_calculated_table1", NOT_STATED),
    (TABLE_III, 1): ("nd_table3", D_LINE),
    (TABLE_III, 3): ("n2_calculated_table3", NOT_STATED),
    (TABLE_III, 4): ("n2_measured_table3", NOT_STATED),
}


def _value_sigma(cell: str) -> tuple[str, str]:
    """Split "v ± s" (Table III's measured column); a bare number has no sigma."""
    value, _, sigma = cell.partition("±")
    return number(value), number(sigma) if sigma else ""


def extract() -> dict[str, list[Point]]:
    """Return one file per material and column."""
    files = {}
    for (table_id, column), (name, wavelength) in COLUMNS.items():
        for cells in table(DIRECTORY, table_id)[1:]:
            material = MATERIALS[cells[0]]
            value, sigma = _value_sigma(cells[column])
            point = Point("1", wavelength, value, sigma, " | ".join(cells))
            files[f"{material}/boling_1978__{name}_{material}.csv"] = [point]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
