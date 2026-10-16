# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Nguyen and Sadeghi (1978): collisional transfer between Ar 3p5 4p levels.

Phys. Rev. A 18, 1388-1395 (1978); journal p. 1388 is PDF page 1, so the journal page is
the PDF page plus 1387. Taken: **Table II** (PDF p. 5, journal p. 1392,
``tbl_c96ebe695f0d``), the experimental rate coefficients k_ij (10^-12 cm^3 s^-1) of
collisional transfer from the i to the j level of Ar 3p5 4p (Paschen 2p_i) in
collisions with ground-state argon at 300 K, all 28 printed cells, every one checked
against the table image (2026-10-12). One file per product level j (the
``channel``), rows the initial level i (``2p3``), ``x`` the same label. Table III's
ratios k_ij / k_ji (from Table II, and the detailed-balance theory of Eq. 6) are
derived and are not taken; they are checked.

**Checks** (:func:`checks`; the script fails if one is missed): the count; Table III's
experimental ratios against Table II (each within the rounding of Table II's cells
and of the ratio: k68 / k86, printed 21.3 where 6 / 0.3 = 20.0, needs k86 near 0.28);
its theoretical
ratios against Eq. (6) with NIST's level energies
(``data/external/nist_asd/Ar_I_levels.tsv``) at 300 K (within 3 %); the text's k3 =
sum_j k3j = 0.76 x 10^-10 cm^3 s^-1 (p. 1395); the energy gap '>= 1361 cm-1' between
2p10 and the levels above it (p. 1391).
"""

import math

from _json_table import json_table
from _library import Point, number, write
from chang_1978 import balance, nist_levels

DIRECTORY = "Nguyen_1978_RateCoefficientsCollisional"
SOURCE = "nguyen_1978"
TABLE_II = "tbl_c96ebe695f0d"
TABLE_III = "tbl_24eb8d79a965"


def file(final: str) -> str:
    """Return the data file of the transfers into a product level."""
    return f"Ar/{SOURCE}__k_to_{final}.csv"


def _matrix() -> dict[tuple[str, str], str]:
    """Return Table II's printed cells by (initial, final) level."""
    grid = json_table(DIRECTORY, TABLE_II)
    header = grid[0]
    if header[1:] != [str(j) for j in range(1, 11)]:
        msg = f"Table II header {header}"
        raise ValueError(msg)
    cells = {}
    for row in grid[1:]:
        for j, cell in zip(header[1:], row[1:], strict=True):
            if cell:
                cells[f"2p{row[0]}", f"2p{j}"] = cell
    return cells


def extract() -> dict[str, list[Point]]:
    """Return one file per product level, rows the initial levels."""
    grid = json_table(DIRECTORY, TABLE_II)
    out: dict[str, list[Point]] = {}
    for (initial, final), cell in sorted(
        _matrix().items(), key=lambda item: (int(item[0][1][2:]), int(item[0][0][2:]))
    ):
        row = next(r for r in grid[1:] if f"2p{r[0]}" == initial)
        printed = f"i = {initial[2:]}, j = {final[2:]}: " + " | ".join(row)
        point = Point(initial, initial, number(cell), "", printed)
        out.setdefault(file(final), []).append(point)
    return out


def checks() -> None:
    """Check the count, Table III, the text's k3 and 2p10's gap (module notes)."""
    k = {key: float(value) for key, value in _matrix().items()}
    if len(k) != 28 or any(i == j for i, j in k):
        msg = f"Table II: {len(k)} cells"
        raise ValueError(msg)
    energy = nist_levels()
    for cells in json_table(DIRECTORY, TABLE_III)[1:]:
        digits = cells[0].removeprefix("k_{").split("}")[0]
        i, j = f"2p{digits[0]}", f"2p{digits[1]}"
        ratio = k[i, j] / k[j, i]
        # Table II's cells are rounded to their last printed digit.
        low = (k[i, j] - 0.05) / (k[j, i] + 0.05)
        high = (k[i, j] + 0.05) / (k[j, i] - 0.05)
        experimental = float(cells[1])
        half = 0.5 * 10 ** -len(cells[1].partition(".")[2])  # the ratio's rounding
        inside = low - half <= experimental <= high + half
        if not inside:
            msg = f"{cells[0]}: {experimental} against {ratio:.3f}"
            raise ValueError(msg)
        theory = balance(i, j, energy)
        if abs(theory / float(cells[2]) - 1) > 0.03:
            msg = f"{cells[0]}: theory {cells[2]} against Eq. (6) {theory:.3f}"
            raise ValueError(msg)
    k3 = sum(value for (i, _), value in k.items() if i == "2p3")
    if round(k3 / 100, 2) != 0.76:
        msg = f"k3 = {k3}e-12"
        raise ValueError(msg)
    gap = min(energy[f"2p{n}"] for n in range(1, 10)) - energy["2p10"]
    if not math.isclose(round(gap), 1361):
        msg = f"2p10 gap {gap}"
        raise ValueError(msg)


def main() -> None:
    """Check the tables and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
