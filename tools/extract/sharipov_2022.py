# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Sharipov, Loukhovitski & Loukhovitskaya (2022): electronic polarizabilities.

Table 3.4 (``tbl_268cb2f58193``, CSV; image checked, every cell of the curated rows
agrees) reprints, per molecule, the results of the authors' 2017 paper (their Ref. 41,
``sharipov_2017``): alpha_el at the B97-2 and MP4(SDQ) levels, the zero-point
correction Delta alpha^ZPV, the pure vibrational polarizability Delta alpha^PV (all in
Angstrom^3) and the coefficients b and c of Delta alpha(T) = bT + cT^2 for the
unperturbed (UP) and perturbed-potential (PP) averages. Only the curated species CO2,
O3, CH4 and SF6 are taken, and only alpha_el(MP4) goes to a data file (a reproduction
of ``sharipov_2017``); the vibrational and thermal terms are corrections without a
quantity of their own, quoted in the source notes. The superscript mark ``a`` of O3's
MP4 cell (triple substitutions included, table note) is a trailing letter in the CSV
("2.801a") and is stripped here.
"""

from _library import Point, number, table, write

DIRECTORY = "Sharipov_2022_InfluenceInternalDegrees"
TABLE = "tbl_268cb2f58193"
GASES = ("CO2", "O3", "CH4", "SF6")
HEADER = [
    "Species",
    "\N{GREEK SMALL LETTER ALPHA}el(B97-2)",
    "\N{GREEK SMALL LETTER ALPHA}el(MP4)",
]


def extract() -> dict[str, list[Point]]:
    """Return alpha_el(MP4) of each curated species, one file each."""
    rows = table(DIRECTORY, TABLE)
    assert rows[0][:3] == HEADER, rows[0]
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        gas = cells[0].strip()
        if gas not in GASES:
            continue
        value = number(cells[2].strip().removesuffix("a"))
        files[f"{gas}/sharipov_2022__alpha_el_mp4_{gas}.csv"] = [
            Point("static", "inf", value, "", " | ".join(cells))
        ]
    assert len(files) == len(GASES), sorted(files)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
