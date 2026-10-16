# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Decleva et al. (2016): vertical ionisation energies of O3, theory.

Table 1 (``tbl_873f20c6cf53``, CSV; image checked, PDF p. 18): per cation state j, its
energy above the neutral's X state and above its B (Hartley) state at the Franck-Condon
point (MRCI-SD(Q), eV), and the squares of the Dyson norms I_ik from X and from B
(CASSCF/aug-cc-pVQZ). One file per initial state, rows labelled by the cation's C2v
state in ASCII ("1^2A_1" -> "1 2A1"). The Dyson norms are not a curated quantity (the
record's notes give them). The two columns differ by the X-B gap, 5.78 eV in the
caption: checked to lie within 5.77-5.79 eV for every row (the printed rounding).

Register: TC-978.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Decleva_2016_AttosecondElectronicNuclear"
TABLE = "tbl_873f20c6cf53"
_STATE = re.compile(r"^\((\d)\^2([AB])_(\d)\)$")


def extract() -> dict[str, list[Point]]:
    """Return the vertical ionisation energies from X and from B."""
    from_x: list[Point] = []
    from_b: list[Point] = []
    for cells in table(DIRECTORY, TABLE)[1:]:
        match = _STATE.match(cells[1].strip())
        if match is None:
            msg = f"unreadable state {cells[1]!r}"
            raise ValueError(msg)
        state = f"{match.group(1)} 2{match.group(2)}{match.group(3)}"
        printed = " | ".join(cells)
        gap = float(cells[2]) - float(cells[4])
        if not 5.765 <= gap <= 5.795:  # 5.78 eV +- 0.01, both energies rounded to 0.01
            msg = f"{state}: X-B gap {gap:.3f} eV, not 5.78"
            raise ValueError(msg)
        from_x.append(Point(state, state, number(cells[2]), "", printed))
        from_b.append(Point(state, state, number(cells[4]), "", printed))
    return {
        "O3/decleva_2016__vertical_from_X.csv": from_x,
        "O3/decleva_2016__vertical_from_B.csv": from_b,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
