# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Daumont et al. (1992): O3 cross sections at the mercury lines, 295 K.

Table III (``tbl_d67497d2c571``, CSV; image checked, every cell agrees): the absorption
cross section at five Hg I lines (253.65-334.15 nm, air wavelengths) from this work,
Molina and Molina, Bass and Paur (renormalised to Mauersberger's 113.6 at 253.65 nm,
note a), Yoshino et al., and the mean +- RMS (with the RMS in per cent in parentheses).
The table prints no unit: the text gives the 253.65 nm value as (113.05 +- 1.1) x
10^-19 cm^2, so the columns are in 10^-19 cm^2 (the datasets' scale). The last row's
cells carry their own factor ``\times 10^{-2}`` (in the mean-value cell after the RMS:
``4.54 +- 0.14(3%) \times 10^{-2}``, the factor applying to both numbers, as the image
shows); these are written in e-notation. Bass and Paur have no value at 253.65 nm (their
reference line) and the mean is not given there: empty cells, no point.
"""

import re

from _library import Point, table, write

DIRECTORY = "Daumont_1992_OzoneUVSpectroscopy"
TABLE_III = "tbl_d67497d2c571"
PREFIX = "O3/daumont_1992__"
COLUMNS = ("present", "molina_1986", "bass_paur_1984", "yoshino_1988")

_FACTOR = re.compile(r"^([\d.]+)\\times 10\^\{(-\d+)\}$")
_MEAN = re.compile(
    r"^([\d.]+) \N{PLUS-MINUS SIGN} ([\d.]+)\(([\d.]+)%\)(?: \\times 10\^\{(-\d+)\})?$"
)


def _value(cell: str) -> str:
    match = _FACTOR.match(cell)
    text = f"{match.group(1)}e{match.group(2)}" if match else cell
    float(text)
    return text


def extract() -> dict[str, list[Point]]:
    """Return Table III's columns, one file each, and the mean +- RMS."""
    files: dict[str, list[Point]] = {f"{PREFIX}{c}.csv": [] for c in (*COLUMNS, "mean")}
    for cells in table(DIRECTORY, TABLE_III)[1:]:
        wavelength, *values, mean = (c.strip() for c in cells)
        printed = " | ".join((wavelength, *values, mean))
        float(wavelength)
        for column, cell in zip(COLUMNS, values, strict=True):
            if cell:
                files[f"{PREFIX}{column}.csv"].append(
                    Point(wavelength, wavelength, _value(cell), "", printed)
                )
        if mean:
            match = _MEAN.match(mean)
            if match is None:
                msg = f"unexpected mean cell {mean!r}"
                raise ValueError(msg)
            value, rms, _, power = match.groups()
            suffix = f"e{power}" if power else ""
            files[f"{PREFIX}mean.csv"].append(
                Point(wavelength, wavelength, value + suffix, rms + suffix, printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
