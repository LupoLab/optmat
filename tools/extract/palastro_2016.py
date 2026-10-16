# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Palastro et al. (2016): Table I, the N2 and O2 parameters (p. 023816-2).

Table I lists the simulation parameters with their reference numbers: fraction, n2
(m2/W), Omega_v (s^-1), I_M (J s^2), Delta a (m3), mu (kg), d a/dQ (m2). Only n2 and the
polarizability anisotropy (Delta a) are extracted; the rest have no quantity. Cells
print ``mantissa \times 10^{k}``; the mantissa is kept and the power goes into each
dataset's scale (1e-24 for n2, 1e-31 m3 = 1e-25 cm3 for Delta a; checked here). The CSV
matches the table image (checked).
"""

from _library import Point, number, table, write

DIRECTORY = "Palastro_2016_TwoPhotonVibrational"
TABLE = "tbl_8f0115c294cd"
ROWS = {"n_{2}(m^{2}/W)": ("n2", "-24"), r"\Delta a(m^{3})": ("dalpha", "-31")}


def extract() -> dict[str, list[Point]]:
    """Return n2 and Delta a of N2 and O2."""
    _, *rows = table(DIRECTORY, TABLE)
    files = {}
    for cells in rows:
        if cells[0] not in ROWS:
            continue
        name, power = ROWS[cells[0]]
        for column, gas in ((2, "N2"), (3, "O2")):
            if not cells[column].endswith(f"10^{{{power}}}"):
                raise ValueError(
                    cells[column]
                )  # the dataset's scale assumes this power
            path = f"{gas}/palastro_2016__{name}_{gas}.csv"
            # Both come from 800-nm measurements (wahlstrand_2012, zahedpour_2015),
            # used by the paper as frequency independent: recorded at 800 nm.
            files[path] = [
                Point("800", "800", number(cells[column]), "", " | ".join(cells))
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
