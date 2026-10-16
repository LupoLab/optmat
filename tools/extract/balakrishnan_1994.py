# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Balakrishnan, Smith and Stoicheff (1994): dissociation energies of H2, D2.

Phys. Rev. A 49, 2460-2469 (1994); journal p. 2460 is PDF page 1. The second
dissociation limit, H2(B' 1Sigma_u+) -> H(1s) + H(2s), measured by XUV
fluorescence-excitation spectroscopy near 84 nm with delayed quenching of the H(2s)
atoms; D0 follows by subtracting the atomic 2s - 1s interval.

Table III (``tbl_26f7348398d0``, CSV; image checked, p. 2468) prints D0(H2), D0(D2)
and, below a second header, D0(H2+) and D0(D2+) in cm^-1: one row per determination,
the authors and their reference numbers in the last column. Cells print digit-group
spaces and the uncertainty after '±' ("36 118.11±0.08"); the spaces are dropped. The
'Present' rows are this paper's values; the others are reproduced, one file per row
and molecule, named by the paper's reference list (p. 2469).

:func:`checks` tests the paper's arithmetic against the numbers it prints (pp.
2465-2468): D_L(H2) = 118 258.57 + F''(1) 118.49; D0 = D_L - E(2s - 1s) for H2 and D2;
De = D0 + the zero-point energy; D0(ion) = D0 + I(atom) - I(molecule); the D2 J'' = 2
limit; the uncertainty budget of Table II (quadrature, plus the hyperfine term
linearly for H2); the average of the two I(H2) values; and the Kr light shift against
the two values of 2 nu1 it explains.
"""

import math
import re
from typing import Final

from _library import Point, number, table, write

DIRECTORY: Final = "Balakrishnan_1994_DissociationEnergiesHydrogen"
TABLE_III: Final = "tbl_26f7348398d0"
_CELL: Final = re.compile(r"^([\d .]+?)(?:±([\d.]+))?$")

# (authors cell as printed, column) -> (file stem, material). The references are the
# paper's list (p. 2469): [2] Wolniewicz 1993, 'see also' Kolos, Szalewicz and
# Monkhorst 1986; [3] McCormack and Eyler 1991; [4] Balakrishnan, Smith and Stoicheff
# 1992; [5] Herzberg 1970; [6] Stwalley 1970; [7] LeRoy and Barwell 1975; [9] Eyler
# and Melikechi 1993; [10] Balakrishnan and Stoicheff 1992; [31] Wolniewicz and
# Orlikowski 1991.
ROWS: Final = {
    ("Herzberg [5]", "H2"): ("herzberg_1970_H2", "H2"),
    ("Herzberg [5]", "D2"): ("herzberg_1970_D2", "D2"),
    ("Stwalley [6]", "H2"): ("stwalley_1970_H2", "H2"),
    ("LeRoy and Barwell [7]", "D2"): ("leroy_1975_D2", "D2"),
    ("McCormack and Eyler [3]", "H2"): ("mccormack_1991_H2", "H2"),
    ("Eyler and Melikechi [9]", "H2"): ("eyler_1993_H2", "H2"),
    ("Eyler and Melikechi [9]", "D2"): ("eyler_1993_D2", "D2"),
    ("Present [4,10]", "H2"): ("d0_H2", "H2"),
    ("Present [4,10]", "D2"): ("d0_D2", "D2"),
    ("Wolniewicz [2]", "H2"): ("wolniewicz_1993_H2", "H2"),
    ("Wolniewicz [2]", "D2"): ("wolniewicz_1993_D2", "D2"),
    ("Kolos, Szalewicz,and Monkhorst [2]", "H2"): ("kolos_1986_H2", "H2"),
    ("Kolos, Szalewicz,and Monkhorst [2]", "D2"): ("kolos_1986_D2", "D2"),
    ("Present [10]", "H2+"): ("d0_H2+", "H2+"),
    ("Present [10]", "D2+"): ("d0_D2+", "D2+"),
    ("Wolniewicz andOrlikowski [31]", "H2+"): ("wolniewicz_1991_H2+", "H2+"),
    ("Wolniewicz andOrlikowski [31]", "D2+"): ("wolniewicz_1991_D2+", "D2+"),
}
"""Table III's cells taken, by the authors cell and the column's molecule."""
HEADERS: Final = (
    ["", "H_{2}", "D_{2}", "Authors"],
    ["", "H_{2}^{+}", "D_{2}^{+}", "Authors"],
)

# Numbers printed in the text, used only by checks() (journal pages).
DL_J1_H2: Final = 118258.57
"""D_L(H2) above X(v = 0, J'' = 1), +-0.04 (p. 2465)."""
F1_H2: Final = 118.49
"""The first rotational term value F''(1) of H2 (p. 2466, their Ref. [18])."""
DL_H2: Final = 118377.06
"""D_L(H2) above X(v = 0, J'' = 0), +-0.04 (abstract; p. 2466)."""
DL_D2: Final = 119029.72
"""D_L(D2) above X(v = 0, J'' = 0), +-0.04 (abstract; p. 2466)."""
DL_J2_D2: Final = 118850.68
"""D_L(D2) above J'' = 2, +-0.09 (p. 2466)."""
F2_D2: Final = 179.029
"""The rotational term value F''(2) of D2 (p. 2466, their Ref. [21])."""
DL_J0_FROM_J2_D2: Final = 119029.71
"""D_L(D2) from the J'' = 2 spectrum, +-0.09 (p. 2466)."""
INTERVAL_H: Final = 82258.954
"""E(2s) - E(1s) of H (p. 2466, their Ref. [19])."""
INTERVAL_D: Final = 82281.336
"""The 2s - 1s interval of D (p. 2466, their Ref. [22])."""
ZERO_POINT: Final = {"H2": 2179.27, "D2": 1546.49}
"""The zero-point vibrational energies, +-0.05 (pp. 2465-2466, their Ref. [20])."""
DE: Final = {"H2": 38297.38, "D2": 38294.87}
"""De(H2) +-0.10 and De(D2) +-0.09 (pp. 2465-2466)."""
I_H2: Final = (124417.507, 124417.484)
"""I(H2) of Gilligan and Eyler [27] and of Jungen et al. [28] (p. 2468)."""
I_H2_MEAN: Final = 124417.500
"""'The average of these two values is I(H2) = 124 417.500+-0.01' (p. 2468)."""
I_ATOM: Final = {"H2": 109678.764, "D2": 109708.617}
"""I(H) [29] and I(D) [28, 30] (p. 2468)."""
I_D2: Final = 124745.353
"""I(D2) (p. 2468)."""
BUDGET: Final = {
    "H2": ((0.04, 0.01, 0.04, 0.03, 0.01), 0.013, 0.08),
    "D2": ((0.04, 0.015, 0.04, 0.03, 0.008), 0.0, 0.07),
}
"""Table II (p. 2467): extrapolation, nu2, 2 nu1, ac Stark, Doppler (in quadrature);
the hyperfine term (added linearly); the quoted total."""
KR_SHIFTED: Final = 94093.00
"""2 nu1 at about 8e8 W/cm2 (pp. 2464, 2466)."""
KR_ZERO: Final = 94092.93
"""2 nu1 extrapolated to zero intensity, +-0.05 (p. 2467)."""
KR_INTENSITY: Final = 8e8
"""The nu1 intensity of the dissociation experiments, '~8 x 10^8 W/cm2' (p. 2466)."""
KR_SLOPE: Final = (8.6e-11, 3.8e-11)
"""The Stark shift, '(8.6+-3.8) x 10^-11 cm-1/W/cm2' (p. 2467)."""


def value_and_sigma(cell: str) -> tuple[str, str]:
    """Return the value and absolute uncertainty ("" if none) of a printed cell."""
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"Table III: not a 'value±uncertainty' cell: {cell!r}"
        raise ValueError(msg)
    value = number(match.group(1).replace(" ", ""))
    return value, number(match.group(2)) if match.group(2) else ""


def _cells() -> dict[tuple[str, str], tuple[str, list[str]]]:
    """Return Table III's non-empty cells by (authors, molecule), with their row."""
    grid = table(DIRECTORY, TABLE_III)
    if grid[0] != HEADERS[0] or grid[9] != HEADERS[1]:
        msg = f"unexpected Table III headers {grid[0]!r}, {grid[9]!r}"
        raise ValueError(msg)
    cells = {}
    for index, row in enumerate(grid):
        if index in {0, 9}:
            continue
        molecules = ("H2", "D2") if index < 9 else ("H2+", "D2+")
        for column, molecule in zip((1, 2), molecules, strict=True):
            if row[column].strip():
                cells[row[3].strip(), molecule] = (row[column], row)
    if set(cells) != set(ROWS):
        msg = f"Table III cells {sorted(set(cells) ^ set(ROWS))} unexpected or missing"
        raise ValueError(msg)
    return cells


def extract() -> dict[str, list[Point]]:
    """Return Table III's dissociation energies, one file per row and molecule."""
    files = {}
    for key, (cell, row) in _cells().items():
        stem, material = ROWS[key]
        value, sigma = value_and_sigma(cell)
        label = f"{row[0]}, {key[0]}"
        point = Point(label, material, value, sigma, " | ".join(row))
        files[f"{material}/balakrishnan_1994__{stem}.csv"] = [point]
    return files


def _close(name: str, got: float, printed: float, tolerance: float) -> None:
    """Raise if a recomputed number is not the printed one within ``tolerance``."""
    if abs(got - printed) > tolerance:
        msg = f"{name}: recomputed {got:.4f}, printed {printed:.4f}"
        raise ValueError(msg)


def checks() -> None:
    """Check the paper's arithmetic and Table III against its text; raise if not."""
    cells = {key: value_and_sigma(cell) for key, (cell, _) in _cells().items()}
    d0 = {m: float(cells["Present [4,10]", m][0]) for m in ("H2", "D2")}
    ion = {m: float(cells["Present [10]", f"{m}+"][0]) for m in ("H2", "D2")}
    _close("D_L(H2) from J'' = 1", DL_J1_H2 + F1_H2, DL_H2, 0.005)
    _close("D0(H2)", DL_H2 - INTERVAL_H, d0["H2"], 0.005)
    _close("D0(D2)", DL_D2 - INTERVAL_D, d0["D2"], 0.005)
    _close("D_L(D2) from J'' = 2", DL_J2_D2 + F2_D2, DL_J0_FROM_J2_D2, 0.005)
    for molecule in ("H2", "D2"):
        _close(
            f"De({molecule})", d0[molecule] + ZERO_POINT[molecule], DE[molecule], 0.005
        )
        squares, hyperfine, total = BUDGET[molecule]
        budget = math.sqrt(sum(term**2 for term in squares)) + hyperfine
        _close(f"Table II total ({molecule})", budget, total, 0.005)
        sigma = float(cells["Present [4,10]", molecule][1])
        _close(f"Table III sigma ({molecule})", sigma, total, 1e-9)
    # The printed average of I(H2) is 124 417.500; the two values average 124 417.4955.
    _close("I(H2) average", sum(I_H2) / 2, I_H2_MEAN, 0.005)
    _close("D0(H2+)", d0["H2"] + I_ATOM["H2"] - I_H2_MEAN, ion["H2"], 0.005)
    _close("D0(D2+)", d0["D2"] + I_ATOM["D2"] - I_D2, ion["D2"], 0.005)
    # The Kr light shift: the 0.07 cm^-1 between the two 2 nu1 values at ~8e8 W/cm2.
    implied = (KR_SHIFTED - KR_ZERO) / KR_INTENSITY
    if abs(implied - KR_SLOPE[0]) > KR_SLOPE[1]:
        msg = f"Kr shift: implied {implied:.2e} against printed {KR_SLOPE[0]:.2e}"
        raise ValueError(msg)


def main() -> None:
    """Check the paper, then write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
