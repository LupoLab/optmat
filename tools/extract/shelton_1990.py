# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton (1990): ESHG phase-match density ratios, 1319 and 1064 nm (Table I).

Table I (``tbl_5295cecf2a78``, JSON cells only; checked against the image) gives, for
pairs of a sample gas A and a reference gas B, the ratio of the densities at which
each periodically phase-matches ESHG, printed as :math:`\rho_B/\rho_A` (reference over
sample), and hyperpolarizability ratios (not extracted: third-order susceptibility, not
refractive index). Because the package's ``phase_match_density_ratio`` is
:math:`\rho_\mathrm{gas}/\rho_\mathrm{reference}`, each printed :math:`\rho_B/\rho_A` is
filed under gas B with reference gas A, so that the printed number is the quantity
without inversion (the source record says so per dataset). Rows whose gases are not
curated species (CF4, C2H6, C3H8, n-C4H10, C4H6) are skipped. The second copy of the
paper in the library (``..._9e0b1a``) prints the same table and is not used.

Third-order data. Table I's hyperpolarizability ratios
:math:`\langle\gamma_A\rangle/\langle\gamma_B\rangle` are filed under the sample
gas A with reference gas B (the quantity's own orientation), ``gamma_A_over_B``.
Table III (``tbl_33fb8f18316f`` part a, ``tbl_9bb3e64084db`` parts b and c; JSON cells,
images checked) gives :math:`\langle\gamma\rangle` in 1e-63 C^4 m^4 J^-3 at 34
wavelengths: the 1319 and 1064 nm rows (marked "a", present work) go to
``table3_<gas>_present``, the visible rows (literature values rescaled by the author)
to ``table3_<gas>_literature``, and the calculated
:math:`\langle\gamma^{vR}\rangle` columns of H2 and D2 to ``table3_<gas>_gammavR``.
The row id is the printed air wavelength (without the "a" mark), the abscissa the
printed vacuum wavenumber. The image of part (c) shows the CH4, CF4 and CO2 columns
under their headings exactly as the JSON cells (columns 3, 5 and 7).
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Shelton_1990_NonlinearOpticalSusceptibilities"
TABLE = "tbl_5295cecf2a78"

SPECIES = {
    "He": "He",
    "Ne": "Ne",
    "Ar": "Ar",
    "Kr": "Kr",
    "Xe": "Xe",
    "H_{2}": "H2",
    "D_{2}": "D2",
    "N_{2}": "N2",
    "O_{2}": "O2",
    "CO_{2}": "CO2",
    "SF_{6}": "SF6",
    "CH_{4}": "CH4",
}
WAVELENGTHS = {2: "1319", 3: "1064"}  # columns of rho_B/rho_A


TABLE3 = {"a": "tbl_33fb8f18316f", "bc": "tbl_9bb3e64084db"}
# (gas, table, first JSON row, last JSON row, value column, gammavR column or None)
TABLE3_COLUMNS = [
    ("Ar", "a", 2, 22, 2, None),
    ("H2", "a", 2, 22, 3, 4),
    ("D2", "a", 2, 22, 5, 6),
    ("N2", "a", 2, 22, 7, None),
    ("O2", "a", 2, 22, 8, None),
    ("Ne", "bc", 2, 12, 2, None),
    ("Kr", "bc", 2, 12, 3, None),
    ("Xe", "bc", 2, 12, 4, None),
    ("SF6", "bc", 2, 12, 5, None),
    ("CH4", "bc", 15, 33, 3, None),
    ("CO2", "bc", 15, 33, 7, None),
]
# The part (c) headings, checked so that a shifted column cannot go unnoticed.
TABLE3_HEADINGS = {
    ("a", 1, 2): "Ar(Ref. 28)",
    ("a", 1, 3): "H_2(Ref. 13)",
    ("a", 1, 4): "H_2^{vR}",
    ("a", 1, 5): "D_2(Ref. 13)",
    ("a", 1, 6): "D_2^{vR}",
    ("a", 1, 7): "N_2(Ref. 28)",
    ("a", 1, 8): "O_2(Ref. 28)",
    ("bc", 1, 2): "Ne(Ref. 30)",
    ("bc", 1, 3): "Kr(Ref. 29)",
    ("bc", 1, 4): "Xe(Ref. 29)",
    ("bc", 1, 5): "\\mathrm{SF}_6(Ref. 31)",
    ("bc", 14, 3): "\\mathrm{CH_4}(Ref. 32)",
    ("bc", 14, 7): "\\mathrm{CO_2}(Ref. 22)",
}


def gamma_ratios() -> dict[str, list[Point]]:
    """Return Table I's gamma ratios, one file per curated pair, filed under gas A."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        sample, reference = cells[0], cells[1]
        if sample not in SPECIES or reference not in SPECIES:
            continue
        gas, other = SPECIES[sample], SPECIES[reference]
        points = []
        for column, wavelength in {4: "1319", 5: "1064"}.items():
            value, error = cells[column].split("\N{PLUS-MINUS SIGN}")
            points.append(
                Point(
                    wavelength,
                    wavelength,
                    number(value),
                    number(error),
                    " | ".join([*cells[:2], *cells[4:6]]),
                )
            )
        files[f"{gas}/shelton_1990__gamma_{gas}_over_{other}.csv"] = points
    return files


def table3() -> dict[str, list[Point]]:
    """Return Table III's <gamma> values: present work, literature, gamma^vR."""
    grids = {key: json_table(DIRECTORY, table) for key, table in TABLE3.items()}
    for (key, row, column), heading in TABLE3_HEADINGS.items():
        if grids[key][row][column] != heading:
            msg = f"Table III heading {grids[key][row][column]!r}, not {heading!r}"
            raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for gas, key, first, last, column, vr_column in TABLE3_COLUMNS:
        present, literature, vr = [], [], []
        for cells in grids[key][first : last + 1]:
            label = cells[0].removesuffix("^a")
            if cells[column]:
                value, error = cells[column].split("\N{PLUS-MINUS SIGN}")
                point = Point(
                    label, cells[1], number(value), number(error), " | ".join(cells)
                )
                (present if cells[0].endswith("^a") else literature).append(point)
            if vr_column is not None and cells[vr_column]:
                vr.append(
                    Point(
                        label, cells[1], number(cells[vr_column]), "", " | ".join(cells)
                    )
                )
        stem = f"{gas}/shelton_1990__table3_{gas}"
        files[f"{stem}_present.csv"] = present
        files[f"{stem}_literature.csv"] = literature
        if vr:
            files[f"{stem}_gammavR.csv"] = vr
    return files


def extract() -> dict[str, list[Point]]:
    """Return every file: density ratios, gamma ratios and Table III's gammas."""
    return {**density_ratios(), **gamma_ratios(), **table3()}


def density_ratios() -> dict[str, list[Point]]:
    """Return one file per curated pair, with the 1319 and 1064 nm ratios."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        sample, reference = cells[0], cells[1]
        if sample not in SPECIES or reference not in SPECIES:
            continue
        gas, other = SPECIES[reference], SPECIES[sample]
        points = []
        for column, wavelength in WAVELENGTHS.items():
            if not cells[column]:
                continue
            value, error = cells[column].split("\N{PLUS-MINUS SIGN}")
            points.append(
                Point(
                    wavelength,
                    wavelength,
                    number(value),
                    number(error),
                    " | ".join(cells[:4]),
                )
            )
        files[f"{gas}/shelton_1990__density_{gas}_over_{other}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
