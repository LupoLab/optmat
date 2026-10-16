# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Spackman (1991): the DOSD Cauchy moments of Meath and co-workers (Table I).

Table I (``tbl_6f890f395ef6``, p. 1290; no CSV rendering, read from its JSON cells and
checked against the image) compares TDHF values of alpha(0) = S(-2), S(-4) and S(-6)
with the "Experiment" columns (note b: "Experimental DOSD alpha(0), S(-4), and S(-6)
from Meath et al. (Refs. 47-60)"). Only the experiment columns are taken, for the
curated species; the TDHF columns (and the ELP-basis values in their parentheses) are
theory that the record describes in its notes.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Spackman_1991_TimeDependentHartree"
TABLE_I = "tbl_6f890f395ef6"
GASES = {"H_2": "H2", "CH_4": "CH4", "N_2": "N2", "CO": "CO", "N_2O": "N2O"}
GASES |= {"CO_2": "CO2", "SF_6": "SF6"}
EXPERIMENT = {"-2": 2, "-4": 4, "-6": 6}  # the order and its "Experiment" column


def extract() -> dict[str, list[Point]]:
    """Return the Meath moments S(-2), S(-4), S(-6) per curated gas."""
    rows = json_table(DIRECTORY, TABLE_I)
    if rows[1][2] != "Experimentb":
        msg = "Table I's second header row has moved"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[2:]:
        if cells[0] not in GASES:
            continue  # NH3, H2O, HF, HCl, H2S, OCS, SO2, CS2 and the organics
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        files[f"{gas}/spackman_1991__meath_{gas}.csv"] = [
            Point(order, order, number(cells[column]), "", printed)
            for order, column in EXPERIMENT.items()
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
