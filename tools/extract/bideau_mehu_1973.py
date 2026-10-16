# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Bideau-Mehu et al. (1973): CO2 n - 1 at 180.7-253.7 nm (Tables 1 and 2).

Both tables print ``(n - 1) x 10^5`` (checked on the images). Table 1: the measured
column ("present work") beside the extrapolated formulas of Cuthbertson and
Cuthbertson, Koch and Old et al., which are not extracted (extrapolations, not
measurements). Its 2062.300 A cell prints 53.07 where Table 2 prints 52.07: an erratum
in the source record. Table 2: the paper's five-term formula ("Calculated values")
against the measured values it was fitted to; each wavelength cell carries a mark,
``*`` for this work and ``[4]`` (Koch) or ``[5]`` (Old et al.) for the reproduced
values, which the script strips from the row id and uses to split the measured column.
The ``*`` rows of the measured column repeat Table 1 and are not extracted again.
"""

from _library import Point, number, table, write

DIRECTORY = "BideauMehu_1973_InterferometricDeterminationRefractive"
TABLE_1 = "tbl_9f56a49e9efd"
TABLE_2 = "tbl_9bf7f917ab3a"
MEASURED = {"[4]": "koch_1914_CO2", "[5]": "old_1971_CO2"}


def _wavelength(cell: str) -> tuple[str, str]:
    """Return a Table 2 wavelength cell's mark ("*", "[4]", "[5]") and number."""
    text = cell.strip()
    if text.startswith("*"):
        return "*", text[1:].strip()
    mark, _, wavelength = text.partition(" ")
    return mark, wavelength.strip()


def extract() -> dict[str, list[Point]]:
    """Return the measured column, the formula column and the reproduced values."""
    files: dict[str, list[Point]] = {
        "CO2/bideau_mehu_1973__experiment_CO2.csv": [
            Point(c[0], c[0], number(c[1]), "", " | ".join(c))
            for c in table(DIRECTORY, TABLE_1)[1:]
        ]
    }
    fit, reproduced = [], {name: [] for name in MEASURED.values()}
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        mark, wavelength = _wavelength(cells[0])
        printed = " | ".join(cells)
        fit.append(Point(wavelength, wavelength, number(cells[1]), "", printed))
        if mark in MEASURED:
            reproduced[MEASURED[mark]].append(
                Point(wavelength, wavelength, number(cells[2]), "", printed)
            )
        elif mark != "*":
            msg = f"unexpected Table 2 mark {mark!r}"
            raise ValueError(msg)
    files["CO2/bideau_mehu_1973__fit_CO2.csv"] = fit
    for name, points in reproduced.items():
        files[f"CO2/bideau_mehu_1973__{name}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
