# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Christiansen, Hattig & Gauss (1998): mean polarizabilities of CO, N2 and Ne.

Only Table X (``tbl_08a3d26f2c43``, cells file only; image checked, cells agree) goes to
data files: it is the one table that prints the isotropic (mean) polarizability alpha,
the quantity of the package; the other tables give the tensor components alpha_xx and
alpha_zz (recorded in the source record's notes) or Cauchy moments and vibrational
oscillators (recorded as formulas). Its columns: label, wavelength (nm, "∞" for the
static limit), then alpha and Delta alpha for CO, N2 and HF, and alpha for Ne. Rows:

- "Electronic^a" (the static row and, as the second "Electronic^a" block, 457.9-632.8
  nm): the vibrationally averaged electronic polarizability, this study.
- "Electronic+Vibrational^a" (static only): the same plus the pure vibrational term.
- "Vibrational^a": the pure vibrational term at the four wavelengths (a contribution,
  not a polarizability: not extracted, quoted in the record).
- "Experiment": literature values marked b (Ref. 36), c (Ref. 58), d (Ref. 55) or f
  (Ref. 57); a blank label continues the block above it. The marker goes to the
  dataset, so that the two 632.8 nm rows (d and c) stay apart.

HF is not extracted. The Delta alpha (anisotropy) columns of CO and N2 were added for
the third-order curation: the same row blocks, files
``anisotropy_<block>_<gas>``; the experimental N2 anisotropies carry the markers d
(Alms et al.), g and c, the CO one c.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Christiansen_1998_PolarizabilitiesCON2"
TABLE_X = "tbl_08a3d26f2c43"
# alpha column of each curated gas.
ALPHA_COLUMNS = {"CO": 2, "N2": 4, "Ne": 8}
# Delta alpha column of each curated gas.
ANISOTROPY_COLUMNS = {"CO": 3, "N2": 5}
# (row label, footnote marker) -> dataset id prefix.
DATASETS = {
    ("Electronica", ""): "electronic",
    ("Electronic+Vibrationala", ""): "total_static",
    ("Experiment", "b"): "werner_meyer",
    ("Experiment", "c"): "bridge_buckingham",
    ("Experiment", "d"): "alms",
    ("Experiment", "f"): "experiment",
    ("Experiment", "g"): "experiment_g",
}
_MARKED = re.compile(r"^(?P<value>[-\d.]+)(?:\^(?P<mark>[a-z]))?$")


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and row block of Table X."""
    rows = json_table(DIRECTORY, TABLE_X)
    alpha = "\N{GREEK SMALL LETTER ALPHA}"
    assert rows[1][2] == alpha and rows[1][8] == alpha, rows[1]
    files: dict[str, list[Point]] = {}
    label = ""
    for cells in rows[2:]:
        label = cells[0].strip() or label
        if label == "Vibrationala":
            continue  # contributions, not polarizabilities
        wavelength = "inf" if cells[1].strip() == "∞" else number(cells[1])
        for gas, column in ALPHA_COLUMNS.items():
            cell = cells[column].strip()
            if not cell:
                continue
            match = _MARKED.match(cell)
            if match is None:
                msg = f"unexpected Table X cell {cell!r}"
                raise ValueError(msg)
            name = DATASETS[label, match["mark"] or ""]
            if name == "total_static" and gas != "CO":
                continue  # N2 has no pure vibrational term: the electronic value again
            files.setdefault(f"{gas}/christiansen_1998__{name}_{gas}.csv", []).append(
                Point(
                    wavelength,
                    wavelength,
                    number(match["value"]),
                    "",
                    " | ".join(cells),
                )
            )
    return files | _anisotropy(rows)


def _anisotropy(rows: list[list[str]]) -> dict[str, list[Point]]:
    """Return the Delta alpha columns of CO and N2, by row block."""
    delta = "\N{GREEK CAPITAL LETTER DELTA}\N{GREEK SMALL LETTER ALPHA}"
    assert rows[1][3] == delta and rows[1][5] == delta, rows[1]
    files: dict[str, list[Point]] = {}
    label = ""
    for cells in rows[2:]:
        label = cells[0].strip() or label
        if label == "Vibrationala":
            continue  # contributions, not anisotropies
        wavelength = "inf" if cells[1].strip() == "\u221e" else number(cells[1])
        for gas, column in ANISOTROPY_COLUMNS.items():
            cell = cells[column].strip()
            if not cell:
                continue
            match = _MARKED.match(cell)
            if match is None:
                msg = f"unexpected Table X cell {cell!r}"
                raise ValueError(msg)
            name = DATASETS[label, match["mark"] or ""]
            if name == "total_static" and gas != "CO":
                continue  # N2 has no pure vibrational term: the electronic value again
            files.setdefault(
                f"{gas}/christiansen_1998__anisotropy_{name}_{gas}.csv", []
            ).append(
                Point(
                    wavelength,
                    wavelength,
                    number(match["value"]),
                    "",
                    " | ".join(cells),
                )
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
