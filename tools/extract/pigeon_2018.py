# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Pigeon et al. (2018): n2 of Kr and Xe at 10.6 um (Table I).

Table I (``tbl_d5296ec423fe``, "Nonlinear refractive indices for noble gases scaled to
1 atm"; image checked) has three columns of n2 in 10^-19 cm^2/W: 0.8 um (Wahlstrand et
al. 2012, ref. 27), 10.6 um "calculated" from the 0.8 um values with Miller's formula
(footnote a) and 10.6 um "measured" (this work). The measured Ar cell, ``<1^b``
("Not measurable in this experiment"), is a bound, not a value, and is not extracted.

Table II (``tbl_2537ec41e92d``, N2, O2, air) is garbled in the library (rows and digits
interleaved); its image shows exactly pigeon_2016's Table 1, so it is not extracted.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Pigeon_2018_ExperimentalStudyThird"
TABLE_I = "tbl_d5296ec423fe"
COLUMNS = {1: ("wahlstrand_2012b", "0.8"), 2: ("miller", "10.6"), 3: ("fwm", "10.6")}
_MARK = re.compile(r"\^[ab]$")


def _plus_minus(cell: str) -> tuple[str, str] | None:
    r"""Split ``a \pm b`` into its two numbers; None for a bound such as ``<1^b``."""
    text = _MARK.sub("", cell.strip())
    if text.startswith("<"):
        return None
    value, sigma = text.split(r"\pm")
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return the three columns for Ar (two), Kr and Xe."""
    header, *rows = table(DIRECTORY, TABLE_I)
    assert "0.8" in header[1] and "calculated" in header[2], header
    assert "measured" in header[3], header
    files: dict[str, list[Point]] = {}
    for cells in rows:
        gas = cells[0]
        assert gas in {"Ar", "Kr", "Xe"}, cells
        for column, (name, wavelength) in COLUMNS.items():
            parsed = _plus_minus(cells[column])
            if parsed is None:  # Ar at 10.6 um: "<1", not measurable
                continue
            value, sigma = parsed
            files[f"{gas}/pigeon_2018__{name}_{gas}.csv"] = [
                Point(wavelength, wavelength, value, sigma, " | ".join(cells))
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
