# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Kumar & Thakkar (2011): dipole sums S_k of O2 and O3 (Tables I and III).

Both tables (checked against their images) list the order k in the first column.

- Table I (``tbl_46bad8286a3e``): S_k of O2, k = 2 ... -8, from the new constrained DOSD
  ('This work') and two earlier constrained DOSDs (Zeiss et al. 1977; Kumar et al.
  1996); '...' cells are skipped.
- Table III (``tbl_ce305bfd8b15``): S_k of O3, k = 2 ... -6, for the five
  distributions A, B1, B2, C and D, and an 'Other' column of literature values with
  footnote letters. Only the k = -2 'Other' cell (static polarizabilities) is
  extracted, one file per footnote letter: '18.83,^{f} 18.65,^{g} 19.27,^{h}
  17.12,^{i} 19.7,^{j} 20.2^k'.
Tables II and IV (mean excitation energies) and V (C6 coefficients) are not
extracted. The trailing '.' of '3134.' etc. is kept as printed.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Kumar_2011_OzoneUnresolvedDiscrepancies"
TABLE_I = "tbl_46bad8286a3e"
TABLE_III = "tbl_ce305bfd8b15"

O2_COLUMNS = {1: "o2_zeiss_1977", 2: "o2_kumar_1996", 3: "o2"}
O3_COLUMNS = {1: "o3_A", 2: "o3_B1", 3: "o3_B2", 4: "o3_C", 5: "o3_D"}
_OTHER = re.compile(r"([\d.]+),?\^\{?([a-z])\}?")


def _k(cell: str) -> str:
    return number(cell)


def extract() -> dict[str, list[Point]]:
    """Return the S_k columns and the static-polarizability literature values."""
    files: dict[str, list[Point]] = {}
    rows = table(DIRECTORY, TABLE_I)[1:]
    for column, name in O2_COLUMNS.items():
        files[f"O2/kumar_2011__{name}.csv"] = [
            Point(_k(c[0]), _k(c[0]), number(c[column]), "", " | ".join(c))
            for c in rows
            if c[column].strip() not in {"...", "\N{HORIZONTAL ELLIPSIS}", ""}
        ]
    rows = table(DIRECTORY, TABLE_III)[1:]
    for column, name in O3_COLUMNS.items():
        files[f"O3/kumar_2011__{name}.csv"] = [
            Point(_k(c[0]), _k(c[0]), number(c[column]), "", " | ".join(c[:6]))
            for c in rows
        ]
    (static,) = (c for c in rows if _k(c[0]) == "-2")
    for value, letter in _OTHER.findall(static[6]):
        files[f"O3/kumar_2011__static_{letter}.csv"] = [
            Point(f"S-2 {letter}", "inf", number(value), "", static[6])
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
