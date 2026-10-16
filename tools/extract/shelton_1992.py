# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton & Donley (1992): the neon ESHG measurements of Table 1.

Table 1 (``tbl_9856decaf635``, CSV cells; checked against the image) lists, per laser
line, the phase-match density ratio rho_He/rho_Ne, the ratio gamma_Ne/gamma_He and
gamma_Ne in atomic units, with footnote marks glued to the wavelength: c) present
work, d) Refs. [1, 2] (shelton_1990, shelton_1989), e) Ref. [11] (Finn & Ward 1971),
f) Ref. [9] (Shelton & Lu 1988, shown here to be spurious). One file per source and
quantity; the row id is the printed air wavelength, the abscissa the printed vacuum
wavenumber. The density ratio is taken for the present work only (the other rows'
density ratios are in ``shelton_1989``).
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Shelton_1992_HyperpolarizabilityDispersionNeon"
TABLE = "tbl_9856decaf635"
GROUPS = {"c": "present", "d": "infrared", "e": "finn_ward", "f": "shelton_lu"}


def extract() -> dict[str, list[Point]]:
    """Return the gamma_Ne values by source, and the present ratios."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        match = re.fullmatch(r"([\d.]+)([cdef])\)", cells[0])
        if match is None:
            msg = f"unexpected wavelength cell {cells[0]!r}"
            raise ValueError(msg)
        label, mark = match.groups()
        group = GROUPS[mark]
        printed = " | ".join(cells)
        gamma, gamma_error = cells[4].split("\N{PLUS-MINUS SIGN}")
        files.setdefault(f"Ne/shelton_1992__gamma_{group}.csv", []).append(
            Point(label, cells[1], number(gamma), number(gamma_error), printed)
        )
        if group != "present":
            continue
        ratio, ratio_error = cells[3].split("\N{PLUS-MINUS SIGN}")
        files.setdefault("Ne/shelton_1992__ratio_present.csv", []).append(
            Point(label, cells[1], number(ratio), number(ratio_error), printed)
        )
        density, density_error = cells[2].split("\N{PLUS-MINUS SIGN}")
        files.setdefault("He/shelton_1992__density_He_over_Ne.csv", []).append(
            Point(label, cells[1], number(density), number(density_error), printed)
        )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
