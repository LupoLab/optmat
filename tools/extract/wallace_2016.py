# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Wallace et al. (2016): the yield curve's parameters, the Kr+ yields.

Table I (``tbl_03b3d0db2cd3``, p. 4): alpha, gamma and I_c of Eq. (4) for the Ar+, Kr+
and Xe+ yields. Per gas, one file of the shape parameters alpha and gamma
(``tableI_shape_<gas>``, dimensionless, printed "arb. units") and one of the
characteristic intensity I_c (``tableI_Ic_<gas>``, in 10^14 W/cm^2, with its
uncertainty), rows labelled by the parameter.

The column headers name the ion produced; the library's text layer gives the
superscript plus as a thorn ("Arþ"). The material is the neutral gas that is
ionised (Ar, Kr, Xe).

STRUCTURAL FIX: the library prints the I_c cells as "3.86 - 0.05". The table image
shows "3.86 ± 0.05", "2.06 ± 0.03" and "1.18 ± 0.03": the plus-minus sign
is lost in the text layer, as the plus of the headers is (checked against the image,
2026-09-28). The second number is the uncertainty, not a subtrahend, and goes in the
sigma column. The script stops if a cell is not of that form, or if the row labels or
headers differ from the ones the image was checked with.

Table I of the Supplemental Material (``supplement_01``, ``tbl_2cab2b41f725``, its p.
7): the 37 measured Kr+ yields, with the columns I_est and I_0 (printed in 10^14
W/cm^2), Y (arbitrary units) and its absolute error. One file,
``supplement_tableI_yield_Kr``: the yield and its error against the calibrated
intensity I_0, with the estimated intensity I_est as an extra column. An intensity
printed "1.31" under "x10^14 W/cm^2" is written "1.31e14". Every cell was checked
against the table image (2026-09-28): they agree. The header's symbols are garbled in
the text layer (no subscripts or superscripts) and are checked as the library has them.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Wallace_2016_PreciseAccurateMeasurements"
TABLE = "tbl_03b3d0db2cd3"
THORN = "\N{LATIN SMALL LETTER THORN}"
GASES = ("Ar", "Kr", "Xe")
HEADER = ["Fit parameter", *(f"{gas}{THORN}" for gas in GASES)]
ROWS = (
    "\N{GREEK SMALL LETTER ALPHA} (arb. units)",
    "\N{GREEK SMALL LETTER GAMMA} (arb. units)",
    "Ic (\N{MULTIPLICATION SIGN}1014 W=cm2)",
)
# "value - uncertainty", the text layer's rendering of "value +- uncertainty".
_PLUS_MINUS = re.compile(r"^(?P<value>\d+\.\d+) - (?P<sigma>\d+\.\d+)$")
SUPPLEMENT = f"{DIRECTORY}/supplement_01"
SUPPLEMENT_TABLE = "tbl_2cab2b41f725"
SUPPLEMENT_HEADER = [
    "Iest (x1014 W/cm2)",
    "YKr+ (arb. units)",
    "\N{GREEK SMALL LETTER SIGMA}YKr+ (abs. error)",
    "I0 (\N{MULTIPLICATION SIGN}1014 W/cm2)",
]
SUPPLEMENT_ROWS = 37


def _krypton_yields() -> dict[str, list[Point]]:
    """Return the measured Kr+ yields of the supplement's Table I."""
    grid = table(SUPPLEMENT, SUPPLEMENT_TABLE)
    if [cell.strip() for cell in grid[0]] != SUPPLEMENT_HEADER or (
        len(grid) != SUPPLEMENT_ROWS + 1
    ):
        msg = f"unexpected layout of the supplement's Table I: {grid[0]!r}"
        raise ValueError(msg)
    points = []
    for row, cells in enumerate(grid[1:], start=1):
        estimated, value, sigma, calibrated = (number(cell) for cell in cells)
        points.append(
            Point(
                f"{row:02d}",
                f"{calibrated}e14",
                value,
                sigma,
                " | ".join(cell.strip() for cell in cells),
                {"intensity_estimated": f"{estimated}e14"},
            )
        )
    return {"Kr/wallace_2016__supplement_tableI_yield_Kr.csv": points}


def extract() -> dict[str, list[Point]]:
    """Return the fit parameters per gas and the measured yields of krypton."""
    grid = table(DIRECTORY, TABLE)
    labels = tuple(row[0].strip() for row in grid[1:])
    if [cell.strip() for cell in grid[0]] != HEADER or labels != ROWS:
        msg = f"unexpected Table I layout {grid!r}"
        raise ValueError(msg)
    alphas, gammas, intensities = grid[1:]
    files: dict[str, list[Point]] = {}
    for column, gas in enumerate(GASES, start=1):
        cell = intensities[column].strip()
        match = _PLUS_MINUS.match(cell)
        if match is None:
            msg = f"Table I, I_c of {gas}: {cell!r} is not 'value - uncertainty'"
            raise ValueError(msg)
        heading = f"{gas}+ (printed {grid[0][column].strip()})"
        files[f"{gas}/wallace_2016__tableI_shape_{gas}.csv"] = [
            Point(
                "alpha",
                "alpha",
                number(alphas[column]),
                "",
                f"{heading} | {alphas[0]} | {alphas[column]}",
            ),
            Point(
                "gamma",
                "gamma",
                number(gammas[column]),
                "",
                f"{heading} | {gammas[0]} | {gammas[column]}",
            ),
        ]
        files[f"{gas}/wallace_2016__tableI_Ic_{gas}.csv"] = [
            Point(
                "I_c",
                "I_c",
                match["value"],
                match["sigma"],
                f"{heading} | {intensities[0]} | {cell} (the image has "
                f"{match['value']} \N{PLUS-MINUS SIGN} {match['sigma']})",
            )
        ]
    return files | _krypton_yields()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
