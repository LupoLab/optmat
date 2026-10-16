# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kumar & Thakkar (2010): dipole sums S_k of Ne, Ar, Kr and Xe.

"In all tables and their footnotes, the numbers in parentheses are powers of ten"
(Table I caption): ``2.586(4)`` is 2.586e4. Values are written in e-notation, the
printed cell kept in ``printed``. All tables were checked against their images (Table I
image; Tables II-IV on the PDF page).

- Table I (``tbl_2e89f4a65bc0``): S_k of the recommended constrained DOSD, k = 2, 1, 0,
  -0.5, -1, -1.5, -2, ..., -10, one file per gas, rows by k.
- Tables II (S_2, ``tbl_4b977c5b7875``), III (S_1, ``tbl_9235bcb67ec8``) and IV
  (S_-1, ``tbl_f451ae7d305a``): comparisons with other work, one file per method row and
  gas (the row of k = 2, 1 or -1). The rows repeating this work (Table I) and Kumar and
  Meath's 1985 values (extracted in kumar_1985, identical) are skipped. Table IV's
  "Correlated" row carries a footnote letter per cell (f: Saxon's CI for Ne, g: Naon et
  al.'s Bethe-Goldstone for Ar), which names the file.
Tables V (mean excitation energies), VI (C6) and VII (C9) are not extracted.
"""

import re

from _library import Point, table, write

DIRECTORY = "Kumar_2010_DipoleOscillatorStrength"
TABLE_I = "tbl_2e89f4a65bc0"
# Comparison table -> (moment order, {method row label: dataset stem}).
COMPARISONS = {
    "tbl_4b977c5b7875": (
        "2",
        {
            "HF^a": "hf_koga",
            "CI^d": "ci_saxon",
            "CI^e": "ci_bunge",
            "CI^f": "ci_bundgen",
        },
    ),
    "tbl_9235bcb67ec8": (
        "1",
        {"HFa": "hf", "CId": "ci_saxon", "MRCIe": "mrci_bundgen"},
    ),
    "tbl_f451ae7d305a": (
        "-1",
        {
            "HFa": "hf",
            "DOSDb": "dosd_olney",
            "TDMP2e": "tdmp2_hattig",
            "MRCIh": "mrci_bundgen",
        },
    ),
}
# Table IV "Correlated" cells: footnote letter -> dataset stem.
CORRELATED = {"f": "ci_saxon", "g": "bethe_goldstone_naon"}
GASES = {1: "Ne", 2: "Ar", 3: "Kr", 4: "Xe"}
_CELL = re.compile(r"^(?P<mantissa>-?\d+(?:\.\d*)?)(?:\((?P<exponent>-?\d+)\))?$")
_LETTER = re.compile(r"^(?P<number>[\d.]+)(?P<letter>[a-z])$")


def _value(cell: str) -> str:
    """Return a printed ``mantissa(power of ten)`` cell's number in e-notation."""
    match = _CELL.match(cell.strip().replace("\N{MINUS SIGN}", "-"))
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    mantissa, exponent = match["mantissa"], match["exponent"]
    return f"{mantissa}e{exponent}" if exponent else mantissa


def _moments() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_I)[1:]:
        k = _value(cells[0])
        printed = " | ".join(cells)
        for column, gas in GASES.items():
            point = Point(k, k, _value(cells[column]), "", printed)
            files.setdefault(f"{gas}/kumar_2010__moments_{gas}.csv", [])
            files[f"{gas}/kumar_2010__moments_{gas}.csv"].append(point)
    return files


def _comparisons() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for table_id, (order, methods) in COMPARISONS.items():
        for cells in table(DIRECTORY, table_id)[1:]:
            label = cells[0].strip()
            printed = " | ".join(cells)
            for column, gas in GASES.items():
                cell = cells[column].strip()
                if not cell:
                    continue
                if label == "Correlated":
                    found = _LETTER.match(cell)
                    if found is None:
                        msg = f"unreadable cell {cell!r}"
                        raise ValueError(msg)
                    stem, cell = CORRELATED[found["letter"]], found["number"]
                elif label in methods:
                    stem = methods[label]
                else:
                    continue  # the constrained-DOSD rows (this work, 1985)
                path = f"{gas}/kumar_2010__s{order}_{stem}_{gas}.csv"
                files[path] = [Point(order, order, _value(cell), "", printed)]
    return files


def extract() -> dict[str, list[Point]]:
    """Return Table I's moments and the literature values of Tables II-IV."""
    return _moments() | _comparisons()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
