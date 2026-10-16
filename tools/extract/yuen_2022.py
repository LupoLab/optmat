# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Yuen & Lin (2022): the N2^2+ states of the DM-SDI model (TC-955).

Table I (PDF p. 5 of the accepted manuscript), read from the library's cells
(``tbl_076ba83089ee.json``; the table has no CSV rendering, its first row is the
neutral's configuration spanning the table). Every cell was checked against the table
image (2026-10-09). Columns: index, state, configuration, the vertical I_p (eV) to the
N2^2+ state from the neutral (Bhattacharya et al. 2021), the dissociation limit
(Iwayama et al. 2016) and the KER (eV).

* ``table1_ip``: the eleven I_p, row = the table's index, x = the state.
* ``table1_ker_3p3p`` and ``table1_ker_1d1d``: the KER of the states that dissociate,
  one file per limit (states 5-9 to N+(3P) + N+(3P); 10 and 11 to N+(1D) + N+(1D)).
  States 1-4 are metastable ("m.s.", KER "-") and give no KER.

The script checks that I_p - KER is one number per limit (the limit above the neutral,
which the paper does not print): 38.96 eV for 3P + 3P, 42.76 eV for 1D + 1D.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Yuen_2022_DensityMatrixApproach"
TABLE = "tbl_076ba83089ee"
HEADER = ["Index", "States", "Config.", "I_{p} (eV)", "Limit", "KER (eV)"]
LIMITS = {
    "N^{+}(^{3}P)+N^{+}(^{3}P)": "3p3p",
    "N^{+}(^{1}D)+N^{+}(^{1}D)": "1d1d",
    # State 10: metastable, or dissociating to 1D + 1D (the table prints both).
    "m.s. / N^{+}(^{1}D)+N^{+}(^{1}D)": "1d1d",
}
# I_p - KER, eV: the limit above the neutral implied by the table (not printed).
IMPLIED = {"3p3p": 38.96, "1d1d": 42.76}
ROUNDING = 0.005  # half the last printed digit of I_p and KER


def extract() -> dict[str, list[Point]]:
    """Return the I_p file and the KER files of Table I."""
    grid = json_table(DIRECTORY, TABLE)
    if not grid[0][0].startswith("N_{2}") or grid[1] != HEADER:
        msg = f"unexpected header of {TABLE}: {grid[:2]!r}"
        raise ValueError(msg)
    energies: list[Point] = []
    ker: dict[str, list[Point]] = {"3p3p": [], "1d1d": []}
    for cells in grid[2:]:
        index, state, _, energy, limit, release = cells
        printed = " | ".join(cells)
        energies.append(Point(index, state, number(energy), "", printed))
        if limit == "m.s.":
            if release != "-":
                msg = f"a metastable state with a KER: {cells!r}"
                raise ValueError(msg)
            continue
        name = LIMITS[limit]
        if abs(float(energy) - float(release) - IMPLIED[name]) > 2 * ROUNDING:
            msg = f"I_p - KER is not the limit's {IMPLIED[name]} eV: {cells!r}"
            raise ValueError(msg)
        ker[name].append(Point(index, state, number(release), "", printed))
    if [p.row for p in energies] != [str(i) for i in range(1, 12)]:
        msg = f"Table I rows {[p.row for p in energies]!r}"
        raise ValueError(msg)
    return {
        "N2/yuen_2022__table1_ip_N2.csv": energies,
        "N2/yuen_2022__table1_ker_3p3p_N2.csv": ker["3p3p"],
        "N2/yuen_2022__table1_ker_1d1d_N2.csv": ker["1d1d"],
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
