# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Grebenshchikov et al. (2007): O3 thresholds and transition moments.

Table 1 (``tbl_6aebd02f5e84``, CSV; image checked, p. 2047): per electronic state of
O3, its dissociation channel (the numbers of eqns (1)-(5)), the calculated threshold
D0 of that channel (eV), and the transition dipole moment mu (D) from the ground state
at the Franck-Condon point, "4.6 (-4)" meaning 4.6 x 10^-4 (footnote o).

* **D0:** the row repeats one value per channel across the states that correlate with
  it; the repeats are checked equal and one value per channel is written.
* **mu:** written in E notation (the bracketed power of ten moved into the number,
  nothing else changed). The triplets' and 1 1A'''s values are this work's (footnote
  o); those of 2 1A'' and 3 1A' are the authors' earlier papers' (refs. 43, 49), a
  separate, reproduced file.

Register: TC-973.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Grebenshchikov_2007_NewTheoreticalInvestigations"
TABLE = "tbl_6aebd02f5e84"
# The table's channel cells (eqn number in brackets) and the data files they go to.
CHANNELS = {"1": "channel_1", "2": "channel_2", "5": "channel_5"}
# The C_s labels of the header, in ASCII, by column.
STATES = ("1 1A' (X)", "1 3A'", "1 3A''", "2 3A''", "1 1A''", "2 1A''", "3 1A'")
REPRODUCED = ("2 1A''", "3 1A'")  # footnote o: from refs. 43 and 49
_CHANNEL = re.compile(r"\((\d)\)")
_POWER = re.compile(r"^([\d.]+) \((-\d+)\)$")


def _moment(cell: str) -> str:
    """Return a printed moment, "4.6 (-4)", in E notation, "4.6e-4"."""
    match = _POWER.match(cell.strip())
    if match is None:
        return number(cell)
    return f"{match.group(1)}e{match.group(2)}"


def extract() -> dict[str, list[Point]]:
    """Return the calculated thresholds D0 and the transition moments."""
    rows = {cells[0]: cells[1:] for cells in table(DIRECTORY, TABLE)}
    channels, d0, c2v = rows["diss. channel"], rows["D0 b"], rows["C2v"]
    files: dict[str, list[Point]] = {}
    for column, (cell, value) in enumerate(zip(channels, d0, strict=True)):
        match = _CHANNEL.search(cell)
        if match is None:
            msg = f"no channel number in {cell!r}"
            raise ValueError(msg)
        path = f"O3/grebenshchikov_2007__d0_{CHANNELS[match.group(1)]}.csv"
        printed = f"{STATES[column]} ({c2v[column]}) | {cell} | D0 = {value}"
        point = Point(match.group(1), match.group(1), number(value), "", printed)
        if path in files:
            if files[path][0].value != point.value:
                msg = f"{path}: {value} differs from {files[path][0].value}"
                raise ValueError(msg)
            continue
        files[path] = [point]
    for column, cell in enumerate(rows["mo"]):
        if not cell:
            continue  # the ground state has no moment
        state = STATES[column]
        dataset = "tdm_earlier" if state in REPRODUCED else "tdm"
        printed = f"{state} ({c2v[column]}) | mu = {cell}"
        files.setdefault(f"O3/grebenshchikov_2007__{dataset}.csv", []).append(
            Point(state, state, _moment(cell), "", printed)
        )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
