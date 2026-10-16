# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Trinh, Tolstikhin, Madsen & Morishita (2013): WFAT coefficients of Ne-Xe.

Table I (``tbl_44397cd5a3f8``, CSV; image checked, p. 043426-9) prints, per noble-gas
atom (the active p electron, m = 0, in a one-electron model potential), in atomic
units: N, u, v (the potential of Eqs. (63)-(64)), E0 and alpha (the model state's
energy and polarizability), g00 and g10 (the asymptotic coefficients of the parabolic
channels (0, 0) and (1, 0), Eq. (37); the structure factor G_nu = g_nu since mu_z = 0),
a00 (Eq. (38)), and the first-order correction coefficients A00 and B00 (Eqs. (40)).

Per atom four files:
- ``wfat_<atom>``: g00 and g10, rows "(0,0)" and "(1,0)"
  (tunnelling_structure_coefficient);
- ``correction_<atom>``: a00, A00, B00 (ionization_model_parameter, model wfat);
- ``potential_<atom>``: N, u, v (ionization_model_parameter, the model potential);
- ``state_<atom>``: E0 and alpha of the model state (ionization_model_parameter).
``printed`` keeps the atom's whole row.

Structural fix: the B00 cells of Ne and Ar are extracted as "- 2.6", "- 7.7"; the image
prints a minus sign set apart by the column alignment, so the space is dropped.
"""

from _library import Point, number, table, write

DIRECTORY = "Trinh_2013_FirstOrderCorrection"
TABLE_I = "tbl_44397cd5a3f8"
ALPHA = "\N{GREEK SMALL LETTER ALPHA}"
HEADER = ["Atom", "N", "u", "v", "E0", ALPHA, "g00", "g10", "a00", "A00", "B00"]
ATOMS = {"Ne(2p)": "Ne", "Ar(3p)": "Ar", "Kr(4p)": "Kr", "Xe(5p)": "Xe"}
FILES = {
    "wfat": (("(0,0)", "g00"), ("(1,0)", "g10")),
    "correction": (("a00", "a00"), ("A00", "A00"), ("B00", "B00")),
    "potential": (("N", "N"), ("u", "u"), ("v", "v")),
    "state": (("E0", "E0"), ("alpha", ALPHA)),
}


def extract() -> dict[str, list[Point]]:
    """Return the four files of each atom."""
    grid = table(DIRECTORY, TABLE_I)
    if grid[0] != HEADER:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        row = dict(zip(HEADER, (cell.strip() for cell in cells), strict=True))
        gas = ATOMS[row["Atom"]]
        printed = " | ".join(cells)
        for stem, entries in FILES.items():
            files[f"{gas}/trinh_2013__{stem}_{gas}.csv"] = [
                Point(label, label, number(row[column].replace(" ", "")), "", printed)
                for label, column in entries
            ]
    if len(files) != 4 * len(ATOMS):
        msg = "Table I atoms missing"
        raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
