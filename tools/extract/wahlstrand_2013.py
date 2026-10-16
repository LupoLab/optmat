# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Wahlstrand et al. (2013): Table I, n2 of Ar, N2, O2 compared (p. 8).

Table I compares the n2 of wahlstrand_2012 (Ref. [17]) with Loriot et al.'s erratum
values (Ref. [10]) and with Loriot's values corrected here for the rotational grating
(N2 and O2 only; Ar has no rotational grating). Each cell prints ``value \pm error`` in
1e-20 cm2/W. The library has no CSV for this table; its JSON cells match the table
image (checked). Each column and gas becomes one single-row data file.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Wahlstrand_2013_EffectTwoBeam"
TABLE = "tbl_2cd9818d136f"
GASES = {"Ar": "Ar", "N_{2}": "N2", "O_{2}": "O2"}
COLUMNS = {1: "ref17", 2: "ref10", 3: "corrected"}


def _value_error(cell: str) -> tuple[str, str]:
    r"""Split a ``value \pm error`` cell into the two printed numbers."""
    value, error = cell.split(r"\pm")
    return number(value), number(error)


def extract() -> dict[str, list]:
    """Return one data file per column and gas (empty cells skipped)."""
    files: dict[str, list] = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        for column, name in COLUMNS.items():
            if not cells[column]:
                continue  # Ar has no rotational-grating correction
            value, error = _value_error(cells[column])
            path = f"{gas}/wahlstrand_2013__{name}_{gas}.csv"
            files[path] = [Point("800", "800", value, error, printed)]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
