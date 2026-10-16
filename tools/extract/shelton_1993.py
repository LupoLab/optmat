# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton & Rugar (1993): the dc Kerr effect of H2 at 632.8 nm.

Table 1 (``tbl_c405d96f0062``, CSV; image checked) gives gamma_eff (Eq. 1) in 1e-62
C^4 m^4 J^-3 at six gas densities and its zero-density extrapolation (the row "0.0a)",
the molecular gamma^evR). The zero-density row goes to ``gamma_evR``, the six measured
rows to ``gamma_eff`` (row id: the density in mol/m^3).

Table 2 (``tbl_84bb0efc5489``, CSV; image checked) gives the H2 polarizability
anisotropies per rotational level J, Delta alpha_JJ and Delta alpha_J,J+2, static and at
632.8 nm, in au: one file per column (row id "J=<J>"). Two cells carry the note mark
"e)" ("2.217 e)", "2.266 e)": interpolated from ref. [14]); the mark is dropped from the
number and stated in the record.
"""

from _library import Point, number, table, write

DIRECTORY = "Shelton_1993_KerrEffectHydrogen"
TABLE1 = "tbl_c405d96f0062"
TABLE2 = "tbl_84bb0efc5489"
COLUMNS = {1: "JJ_static", 2: "JJ_632nm", 3: "JJ2_static", 4: "JJ2_632nm"}


def _split(cell: str) -> tuple[str, str]:
    value, error = cell.replace("a)", "").split("\\pm")
    return number(value), number(error)


def extract() -> dict[str, list[Point]]:
    """Return gamma_evR, the gamma_eff series and the anisotropy columns."""
    files: dict[str, list[Point]] = {}
    measured = []
    for cells in table(DIRECTORY, TABLE1)[1:]:
        value, error = _split(cells[1])
        printed = " | ".join(cells)
        if cells[0] == "0.0a)":
            files["H2/shelton_1993__gamma_evR.csv"] = [
                Point("zero density", "632.8", value, error, printed)
            ]
        else:
            row = f"rho = {number(cells[0])} mol/m3"
            measured.append(Point(row, "632.8", value, error, printed))
    files["H2/shelton_1993__gamma_eff.csv"] = measured
    rows = table(DIRECTORY, TABLE2)[1:]
    for column, name in COLUMNS.items():
        wavelength = "inf" if name.endswith("static") else "632.8"
        files[f"H2/shelton_1993__anisotropy_{name}.csv"] = [
            Point(
                f"J={cells[0]}",
                wavelength,
                number(cells[column].removesuffix(" e)")),
                "",
                " | ".join(cells),
            )
            for cells in rows
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
