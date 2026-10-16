# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Piszczatowski et al. (2015): helium polarizability (Table I totals).

Table I (``tbl_df068b216aa3``, cells file only) gives the components of the 4He
polarizability at 632.9908 nm; the static "Total" (row 5) and the dynamic "Total,
Present" (row 18) are extracted, with the last-digit uncertainty in parentheses. The
cells print digit groups with spaces (``1.383 760 77(14)``), removed here. The table
image agrees with the cells. Table II (virial expansion of n) is not extracted: its
a_r include the magnetic term and its b_n, c_n are not the package's quantities (see
the source record).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Piszczatowski_2015_FrequencyDependentPolarizability"
TABLE_I = "tbl_df068b216aa3"
_VALUE = re.compile(r"^(-?[\d.]+)(?:\((\d+)\))?$")


def split(cell: str) -> tuple[str, str]:
    """Return the value and the last-digits uncertainty of ``1.234 5(6)``."""
    match = _VALUE.match(cell.replace(" ", ""))
    if match is None:
        msg = f"cannot read {cell!r}"
        raise ValueError(msg)
    value, digits = match.group(1), match.group(2)
    if digits is None:
        return number(value), ""
    decimals = len(value.split(".")[1]) if "." in value else 0
    return number(value), f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return the static and the 632.9908 nm totals."""
    rows = json_table(DIRECTORY, TABLE_I)
    assert rows[5][1] == "Total", rows[5]
    assert rows[18][:2] == ["Total", "Present"], rows[18]
    static = split(rows[5][2])
    dynamic = split(rows[18][2])
    return {
        "He/piszczatowski_2015__alpha.csv": [
            Point("static", "inf", *static, " | ".join(rows[5])),
            Point("632.9908", "632.9908", *dynamic, " | ".join(rows[18])),
        ]
    }


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
