# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Sharipov, Loukhovitski & Starik (2017): static electronic polarizabilities.

Table A1 (``tbl_6536153eecc8`` and its continuation ``tbl_5c6a04a7a846``, CSV; images
checked, cells agree) lists, per molecule, alpha_e at the B97-2 and MP4 levels, the
zero-point vibrational correction Delta alpha^ZPV, the pure vibrational
polarizability Delta alpha^PV (all in Angstrom^3) and the coefficients of
Delta alpha^UP(T) = bT + cT^2.
Only the curated species CO2, O3, CH4 and SF6 are taken, and only alpha_e(MP4) (the
authors' electronic value, at the equilibrium geometry) goes to a data file: the
vibrational terms are corrections, not polarizabilities, and have no quantity of their
own (they are quoted in the source notes). The superscript mark ``a`` of O3's MP4 cell
(triple substitutions included, table note) is lost as a trailing letter in the CSV
("2.801a") and stripped here.
"""

from _library import Point, number, table, write

DIRECTORY = "Sharipov_2017_InfluenceVibrationsPolyatomic"
TABLES = ("tbl_6536153eecc8", "tbl_5c6a04a7a846")
GASES = {"CO2": "CO2", "O3": "O3", "CH4": "CH4", "SF6": "SF6"}
HEADER = ["Species", "ae(B97-2)", "ae(MP4)", "DaZPV", "DaPV", "b · 106", "c · 109"]


def extract() -> dict[str, list[Point]]:
    """Return alpha_e(MP4) of each curated species, one file each."""
    files: dict[str, list[Point]] = {}
    for table_id in TABLES:
        rows = table(DIRECTORY, table_id)
        assert rows[0] == HEADER, rows[0]
        for cells in rows[1:]:
            gas = GASES.get(cells[0].strip())
            if gas is None:
                continue
            value = number(cells[2].strip().removesuffix("a"))
            files[f"{gas}/sharipov_2017__alpha_e_mp4_{gas}.csv"] = [
                Point("static", "inf", value, "", " | ".join(cells))
            ]
    assert len(files) == len(GASES), sorted(files)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
