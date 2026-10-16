# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rollefson & Havens (1940): CH4 refractivity at 49 wavelengths, 1.68-14.8 um.

Table I (``tbl_d3b491ff89fa``, p. 713; image checked, and every number agreeing with
the PDF's text layer) has seven columns: the wavelength in microns; ``(n - 1)_EXP x
10^6``, the average of several runs reduced to 0 degrees C and 760 mm with the
International Critical Tables' expansion coefficient and pv data (p. 712);
``(n - 1)_EXT x 10^6``, Friberg's visible and ultraviolet dispersion extrapolated into
the infrared and adjusted vertically (p. 713); their difference ``Delta(n - 1)``; the
two band constants ``Delta_3 (nu_3^2 - nu^2)`` and ``Delta_4 (nu_4^2 - nu^2)``, some in
parentheses (near the bands and beyond 11 um, left out of the authors' averages); and
``(n - 1)_CALC x 10^6``, the authors' two-band model with the averages K3 = 16.50 and
K4 = 7.50.

The experimental column is the measurement and the calculated column the authors' fit
(a check). The extrapolated column (a model of other data), the difference and the
band constants (derived) are not extracted; the difference is checked here on every
run against |EXP - EXT| at the printed precision, which catches a misread or shifted
cell in either column. Row ids are the wavelengths as printed.
"""

from _library import Point, number, table, write

DIRECTORY = "Rollefson_1940_IndexRefractionMethane"
TABLE = "tbl_d3b491ff89fa"
ROWS = 49  # "measured at 49 different wave-lengths between 1.6 and 14.8 u" (p. 712)


def _decimals(text: str) -> int:
    """Return the number of decimals printed in a number."""
    return len(text.split(".")[1]) if "." in text else 0


def _check_difference(cells: list[str]) -> None:
    """Check the printed Delta(n - 1) against |EXP - EXT| at its printed precision."""
    experimental, extrapolated, difference = (number(cell) for cell in cells[1:4])
    computed = abs(float(experimental) - float(extrapolated))
    if f"{computed:.{_decimals(difference)}f}" != difference:
        msg = f"row {cells[0]}: |EXP - EXT| = {computed} against Delta {difference}"
        raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return the experimental and calculated columns of Table I."""
    rows = table(DIRECTORY, TABLE)[1:]
    if len(rows) != ROWS:
        msg = f"Table I has {len(rows)} rows, not {ROWS}"
        raise ValueError(msg)
    experiment, calculated = [], []
    for cells in rows:
        _check_difference(cells)
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        experiment.append(Point(wavelength, wavelength, number(cells[1]), "", printed))
        calculated.append(Point(wavelength, wavelength, number(cells[6]), "", printed))
    return {
        "CH4/rollefson_1940__experiment_CH4.csv": experiment,
        "CH4/rollefson_1940__calc_CH4.csv": calculated,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
