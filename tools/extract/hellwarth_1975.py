# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hellwarth, Cherlow & Yang (1975): the nonlinear parameters of five glasses.

Table I (``tbl_aa2b556d00dc``, p. 967; JSON only), one row per glass: n; 10^16 B0 and
10^16 (A0 + B0) (this work, from Raman spectra); 10^15 (sigma + 2 B0) (Owyoung's IIPC,
694 nm); 10^15 sigma; 10^14 n2 of Eq. (9) (the header prints "[Eq. (4)]": the caption
and the text make it Eq. (9)); 10^14 n2 measured; f_e (%). Every cell is ``v ± u``.

- sigma, sigma + 2 B0 and n2 [Eq. (9)] are extracted for all five glasses, the +- as
  sigma; the IIPC columns at 694 nm, n2 at "not stated" (a long-wavelength prediction).
- n2 (Meas.) is printed for fused quartz (``18 ± 3^c``) and, as three values under a
  brace, for ED-4 (``\begin{cases} 26 ± 3^c \\ 15 ± 2^d \\ 15 ± 3^e \end{cases}``); the
  rows are named by the footnote letter. The other glasses' cells are empty.
- f_e, the electronic fraction of the long-pulse n2 of Eq. (9) (linear polarization),
  is recorded as its complement, the nuclear share 1 - f_e, a ``raman_fraction``
  (``nuclear_fraction_*``): the percent cell's integer is subtracted from 100, exactly,
  its +- kept (TC-694). B0 and A0 + B0 have no quantity and stay in ``printed``.

All cells were compared with the table image (pdf p. 4): no misread cell.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Hellwarth_1975_OriginFrequencyDependence"
TABLE = "tbl_aa2b556d00dc"
IIPC_NM = "694"

HEADER = [
    "Material",
    "n",
    "10^{16}B_0^a",
    "10^{16}(A_0 + B_0)^a",
    r"10^{15}(\sigma + 2B_0)^b",
    r"10^{15}\sigma",
    "10^{14}n_2[Eq. (4)]",
    "10^{14}n_2(Meas.)",
    "f_e(%)",
]

# Table I "Material" cell -> registered material id.
GLASSES = {
    "Fused quartz": "SiO2",
    "LSO glass": "LSO",
    "ED-4 glass": "ED-4",
    "SF-7 glass": "SF-7",
    "LASF-7 glass": "LaSF-7",
}

_PLUS_MINUS = re.compile(r"^([\d.]+)\s*±\s*([\d.]+)$")
_MEASURED = re.compile(r"([\d.]+)\s*±\s*([\d.]+)\^([a-e])")


def _value(cell: str) -> tuple[str, str]:
    """Return the value and the +- of a cell ``v ± u``, as printed."""
    match = _PLUS_MINUS.match(cell.strip())
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    return number(match.group(1)), number(match.group(2))


def _measured(cell: str) -> list[tuple[str, str, str]]:
    """Return (footnote letter, value, +-) of each measured n2 in a cell."""
    found = _MEASURED.findall(cell)
    return [(letter, number(value), number(sigma)) for value, sigma, letter in found]


def extract() -> dict[str, list[Point]]:
    """Return the data files of Table I."""
    rows = json_table(DIRECTORY, TABLE)
    if rows[0] != HEADER:
        msg = f"unexpected header {rows[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    body = [cells for cells in rows[1:] if cells[0] in GLASSES]
    if [cells[0] for cells in body] != list(GLASSES):
        msg = f"unexpected rows {[cells[0] for cells in rows]!r}"
        raise ValueError(msg)
    for cells in body:
        material = GLASSES[cells[0]]
        printed = " | ".join(cells)
        stem = f"{material}/hellwarth_1975__"
        value, sigma = _value(cells[5])
        files[f"{stem}sigma_{material}.csv"] = [
            Point("1", IIPC_NM, value, sigma, printed)
        ]
        value, sigma = _value(cells[4])
        files[f"{stem}iipc_{material}.csv"] = [
            Point("1", IIPC_NM, value, sigma, printed)
        ]
        value, sigma = _value(cells[6])
        files[f"{stem}n2_{material}.csv"] = [
            Point("1", "not stated", value, sigma, printed)
        ]
        electronic, spread = _value(cells[8])
        # The complement of the printed percent, exactly (both are whole percents).
        files[f"{stem}nuclear_fraction_{material}.csv"] = [
            Point(
                "1-f_e",
                "1-f_e",
                f"{(100 - int(electronic)) / 100:.2f}",
                f"{int(spread) / 100:.2f}",
                printed,
            )
        ]
        measured = _measured(cells[7])
        if bool(measured) != (material in {"SiO2", "ED-4"}):
            msg = f"unexpected measured n2 cell {cells[7]!r}"
            raise ValueError(msg)
        if measured:
            # One row per value, named by its footnote letter (c, d, e); fused quartz
            # has one (c), ED-4 three under a brace.
            files[f"{stem}n2_measured_{material}.csv"] = [
                Point(letter, "not stated", value, sigma, printed)
                for letter, value, sigma in measured
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
