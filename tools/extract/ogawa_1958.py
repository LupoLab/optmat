# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Ogawa & Cook (1958): O3 absorption coefficients, 526-1305 A.

Table I (``tbl_987cd6327574``, CSV; image checked row by row): the approximate
wavelength (A), the relative absorption coefficient in arbitrary units printed as
``value +- sigma`` (Lyman alpha set to 500), its percent error, and the absolute
absorption coefficient mu (cm^-1, at 0 C and 760 mm Hg), scaled to Tanaka, Inn and
Watanabe's 614 cm^-1 at 1216 A.

What the image shows and the cells lose (OCR misreadings, fixed here):

- three absolute values carry a subscript last digit, the old convention for a digit
  that is not significant: ``119_1`` (718 A), ``111_2`` (703 A) and ``100_8`` (645 A).
  The cells read ``1191``, ``llh`` and ``100.``; the digits are ``1191``, ``1112`` and
  ``1008`` (each reproduces the relative value times 614/500 to 0.2 %);
- the letter O for a zero in two relative cells (``41O+-14.2`` at 879 A and
  ``704+-1O.2`` at 539 A): the image prints 410 and 10.2.

Two files: the absolute coefficients (the values), and the relative coefficients with
their printed +- as the per-point sigma (a check of the per-point precision).
"""

from _library import Point, table, write

DIRECTORY = "Ogawa_1958_AbsorptionCoefficientsO3"
TABLE_I = "tbl_987cd6327574"
PREFIX = "O3/ogawa_1958__"

# Absolute cells the OCR misread: the cell, and the digits the image prints (the last
# one as a subscript).
_ABSOLUTE_OCR = {"llh": "1112", "100.": "1008"}


def _digits(cell: str) -> str:
    text = cell.strip().replace("O", "0")  # the letter O read for a zero
    float(text)
    return text


def extract() -> dict[str, list[Point]]:
    """Return Table I's absolute and relative absorption coefficients."""
    absolute, relative = [], []
    for cells in table(DIRECTORY, TABLE_I)[3:]:
        wavelength, rel, percent, mu = (c.strip() for c in cells)
        printed = " | ".join((wavelength, rel, percent, mu))
        value, sigma = (_digits(part) for part in rel.split("\N{PLUS-MINUS SIGN}"))
        float(percent)
        mu_value = _ABSOLUTE_OCR.get(mu, mu)
        float(mu_value)
        absolute.append(Point(wavelength, wavelength, mu_value, "", printed))
        relative.append(Point(wavelength, wavelength, value, sigma, printed))
    return {f"{PREFIX}absolute.csv": absolute, f"{PREFIX}relative.csv": relative}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
