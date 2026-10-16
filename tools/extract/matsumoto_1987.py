# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Matsumoto (1987): CO2 refractivity at the 10.6 um CO2-laser line.

Table 2 (``tbl_a8b02d951b7a``, image checked) gives ``(n - 1)_m x 10^6`` of CO2 at
101325 Pa and 20 degrees C on the 10.4 um band's P(20) line, at its centre frequency
v0 and 25 MHz either side, with ``(n_m - n_c) x 10^6`` against Koch's formula. The
measured column is extracted, one point per row; the Koch column (a calculation, its
last row's sign misprinted; see the source record) is not.

The table prints no wavelength: the text calls the laser the "10.6-um CO2-laser"
(p. 148) and Fig. 2 labels its fringes "lambda = 10.6 um", so every row's abscissa is
the printed 10.6 um. The three rows differ by 25 MHz, 1e-6 of the frequency.
"""

from _library import Point, number, table, write

DIRECTORY = "Matsumoto_1987_RefractivityCO2Gas"
TABLE = "tbl_a8b02d951b7a"
FILE = "CO2/matsumoto_1987__measured_CO2.csv"
WAVELENGTH_UM = "10.6"  # the text's (p. 148) and Fig. 2's label of the CO2 laser
# The first column's LaTeX labels -> row ids.
ROWS = {
    r"v_0 - \Delta\nu": "v0-dv",
    "v_0": "v0",
    r"v_0 + \Delta\nu": "v0+dv",
}


def extract() -> dict[str, list[Point]]:
    """Return Table 2's measured refractivities as one data file."""
    header, *rows = table(DIRECTORY, TABLE)
    if not header[1].startswith("(n-1)_m"):
        msg = f"unexpected header {header!r}"
        raise ValueError(msg)
    points = []
    for cells in rows:
        row = ROWS.get(cells[0].strip())
        if row is None:
            msg = f"unexpected row label {cells[0]!r}"
            raise ValueError(msg)
        value = number(cells[1])
        points.append(Point(row, WAVELENGTH_UM, value, "", " | ".join(cells)))
    found = [p.row for p in points]
    if found != list(ROWS.values()):
        msg = f"expected the rows {list(ROWS.values())}, found {found}"
        raise ValueError(msg)
    return {FILE: points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
