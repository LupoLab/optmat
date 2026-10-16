# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Pigeon et al. (2016): n2,eff of N2 and O2 at 10.6 um (Table 1).

Table 1 (``tbl_a7c0c0dc1804``, "Measured Values of n2,eff of Major Air Constituents,
Scaled to Atmospheric Pressure") has three columns of n2,eff in 10^-19 cm^2/W: 0.8 um
(Nibbering et al. 1997, ref. 11), 2.4 um ("Calculated from measurements in Zahedpour et
al. using n2,eff = n2 + n2,rot", ref. 13; footnote b) and 10.6 um (this work). The CSV
equals the image except that the footnote marks are kept as ``^{b}`` (stripped here) and
the "+" of footnote b is lost in the markdown (read from the page image). Air is
skipped.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Pigeon_2016_MeasurementsNonlinearRefractive"
TABLE = "tbl_a7c0c0dc1804"
COLUMNS = {
    1: ("nibbering_1997", "0.8"),
    2: ("zahedpour_2015", "2.4"),
    3: ("fwm", "10.6"),
}
_MARK = re.compile(r"\^\{[ab]\}$")


def _plus_minus(cell: str) -> tuple[str, str]:
    r"""Split ``a \pm b`` (with an optional footnote mark) into its two numbers."""
    value, sigma = _MARK.sub("", cell.strip()).split(r"\pm")
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return the three columns for N2 and O2."""
    header, _units, *rows = table(DIRECTORY, TABLE)
    assert "0.8" in header[1] and "2.4" in header[2] and "10.6" in header[3], header
    files: dict[str, list[Point]] = {}
    for cells in rows:
        gas = {"N_{2}": "N2", "O_{2}": "O2"}.get(cells[0])
        if gas is None:  # air
            continue
        for column, (name, wavelength) in COLUMNS.items():
            value, sigma = _plus_minus(cells[column])
            files[f"{gas}/pigeon_2016__{name}_{gas}.csv"] = [
                Point(wavelength, wavelength, value, sigma, " | ".join(cells))
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
