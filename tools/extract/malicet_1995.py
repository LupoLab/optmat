# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Malicet et al. (1995): O3 cross sections, temperature dependence.

Table III (``tbl_b1fba2648b6c``, cells file only; image checked, every cell agrees):
absolute cross sections at 253.65 nm against temperature, this work (295, 243, 228,
218 K, printed ``113.7 \pm 1.5``; "-" at 273 K) and Barnes and Mauersberger (297, 273,
253, 221 K; "-" at 228 K), with their relative difference (derived, not extracted). The
unit is the header's ``sigma x 10^19 cm^2``. The temperature differs from row to row
and a dataset has one temperature, so each temperature is its own file,
``t<T>_253`` (this work, the +- as sigma) and ``barnes_mauersberger_<T>``.

Table IV (``tbl_418f10b3b9e9``, CSV; image checked, every cell agrees): the cross
sections at 228 K at 13 fixed wavelengths (230-345 nm), this work and Yoshino et al.
(1988), in sigma x 10^19 cm^2; the relative-difference column is derived and not
extracted. Table IV's 253.65 nm value (113.5) is Table III's 228 K value.
"""

import re

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Malicet_1995_OzoneUVSpectroscopy"
TABLE_III = "tbl_b1fba2648b6c"
TABLE_IV = "tbl_418f10b3b9e9"
PREFIX = "O3/malicet_1995__"
WAVELENGTH_III = "253.65"  # Table III's caption: "Absolute cross-sections at 253.65 nm"

_PLUS_MINUS = re.compile(r"^([\d.]+) \\pm ([\d.]+)$")


def _table_iii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_III)[2:]:
        t_present, present, t_barnes, barnes, difference = cells
        printed = " | ".join(cells)
        if present != "-":
            match = _PLUS_MINUS.match(present)
            if match is None:
                msg = f"unexpected cell {present!r}"
                raise ValueError(msg)
            value, sigma = match.groups()
            files[f"{PREFIX}t{int(t_present)}_253.csv"] = [
                Point(WAVELENGTH_III, WAVELENGTH_III, value, sigma, printed)
            ]
        if barnes != "-":
            float(barnes)
            files[f"{PREFIX}barnes_mauersberger_{int(t_barnes)}.csv"] = [
                Point(WAVELENGTH_III, WAVELENGTH_III, barnes, "", printed)
            ]
        del difference  # derived from the two columns: not extracted
    return files


def _table_iv() -> dict[str, list[Point]]:
    present, yoshino = [], []
    for cells in table(DIRECTORY, TABLE_IV)[1:]:
        wavelength, ours, theirs, difference = (c.strip() for c in cells)
        printed = " | ".join((wavelength, ours, theirs, difference))
        for text in (wavelength, ours, theirs):
            float(text)
        present.append(Point(wavelength, wavelength, ours, "", printed))
        yoshino.append(Point(wavelength, wavelength, theirs, "", printed))
    return {f"{PREFIX}t228.csv": present, f"{PREFIX}yoshino_1988_228.csv": yoshino}


def extract() -> dict[str, list[Point]]:
    """Return Table III per temperature and Table IV's two columns."""
    return _table_iii() | _table_iv()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
