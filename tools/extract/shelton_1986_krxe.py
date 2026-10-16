# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton (1986): the ESHG hyperpolarizabilities of Kr and Xe (Tables I, II).

Table I (``tbl_c9d92892bdc1``, CSV; image checked) gives gamma_Kr/gamma_Ar,
gamma_Xe/gamma_Ar and the absolute gamma_Kr, gamma_Xe (1e-61 C^4 m^4 J^-3) at four
laser lines, with 'the uncertainty of the last digits ... in parentheses': "2.359(6)" is
2.359 +- 0.006. The script turns the parenthesised digits into the sigma column (same
decimal places as the value). Table II (``tbl_e36736d86d09``, CSV; image checked) gives
gamma_Kr and gamma_Xe from other processes and authors (dc Kerr, FWM, ESHG, THG), one
file per gas and process; Xe's FWM cell is "..." (none). The abscissa of FWM is the
pump of the printed pair "532,683". Row ids are the printed air wavelengths.
"""

import re
from decimal import Decimal

from _library import Point, number, table, write

DIRECTORY = "Shelton_1986_HyperpolarizabilityDispersionMeasured_5d296f"
TABLE1 = "tbl_c9d92892bdc1"
TABLE2 = "tbl_e36736d86d09"
PROCESSES = {
    "dc Kerr^{\\text{b}}": ("dc_kerr", "ref 19 Buckingham and Dunmur"),
    "FWM^{\\text{c}}": ("cars", "ref 10 Lundeen et al."),
    "ESHG^{\\text{d}}": ("eshg", "ref 20 Finn and Ward"),
    "THG^{\\text{e}}": ("thg", "ref 21 Ward and New"),
}
_PARENTHESES = re.compile(r"^(\d+\.(\d+))\((\d+)\)$")


def parenthesised(cell: str) -> tuple[str, str]:
    """Return the value and sigma of a cell printed as "2.359(6)"."""
    match = _PARENTHESES.match(cell.strip())
    if match is None:
        msg = f"not a value(uncertainty) cell: {cell!r}"
        raise ValueError(msg)
    value, decimals, digits = match.groups()
    sigma = Decimal(digits).scaleb(-len(decimals))
    return value, f"{sigma:f}"


def extract() -> dict[str, list[Point]]:
    """Return the Table I ratios and gammas and the Table II values."""
    columns = {
        1: "Kr/shelton_1986_krxe__gamma_Kr_over_Ar",
        2: "Xe/shelton_1986_krxe__gamma_Xe_over_Ar",
        3: "Kr/shelton_1986_krxe__gamma_Kr",
        4: "Xe/shelton_1986_krxe__gamma_Xe",
    }
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE1)[1:]:
        label = number(cells[0])
        for column, stem in columns.items():
            value, sigma = parenthesised(cells[column])
            files.setdefault(f"{stem}.csv", []).append(
                Point(label, label, value, sigma, " | ".join(cells))
            )
    for cells in table(DIRECTORY, TABLE2)[1:]:
        process, row = PROCESSES[cells[0]]
        wavelength = number(cells[1].split(",")[0])
        for column, gas in ((2, "Kr"), (3, "Xe")):
            if cells[column] == "...":
                continue
            value, error = cells[column].split("\\pm")
            path = f"{gas}/shelton_1986_krxe__table2_{gas}_{process}.csv"
            files[path] = [
                Point(row, wavelength, number(value), number(error), " | ".join(cells))
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
