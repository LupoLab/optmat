# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Böwering, Bruce and Keto (1986), part I: Xe 6p state-to-state rates in Xe.

J. Chem. Phys. 84, 709-714 (1986); journal p. 709 is PDF page 2 (PDF page 1 is the AIP
download cover), so the journal page is the PDF page plus 707. A scanned (raster) PDF:
every cell of Tables I-III was checked against the table images (2026-10-09). The
cells come from the library's ``.json`` (merged headers).

Taken: **Table II** (PDF p. 4, journal p. 711, ``tbl_9b11ef553745``), the state-to-state
rate constants k_ij of the three two-photon excited levels (columns, the reactants
6p[1/2]0, 6p[3/2]2, 6p[5/2]2) into each product level (rows), perturber xenon. One file
per product (the dataset's ``channel``), the row the reactant in the paper's Racah
label, ``6p[1/2]0``; the row Sigma_j nu_ij, the sum of the observed channels, is a file
of its own. '...' (the reactant itself) and 'b' (no fluorescence detected) are not
values. The caption prints the unit as '10^-2 cm^3 s^-1', a misprint for 10^-12
(record erratum); the files keep the printed numbers.

**Checks** (:func:`checks`; the script fails if one is missed): the cell counts; each
column's Sigma_j nu_ij against the sum of its printed rates (6.57, 87.3, 77.0 against
6.6, 87, 77); Table III's measured intramultiplet fractions (k_ij over the sum of the
6p products) against Table II for the reactants 6p[1/2]0 and 6p[3/2]2 (within 0.006;
the 6p[5/2]2 column does not reproduce, record notes); and Table I's wavelengths, which
are air wavelengths, against NIST's vacuum wavelengths (TC-985) converted to air with
Edlén's formula (all ten within 0.05 nm), the evidence of the level mapping.
"""

import csv
import re

from _json_table import json_table
from _library import DATA, Point, number, write

DIRECTORY = "Bowering_1986_CollisionalDeactivationTwo"
SOURCE = "bowering_1986"
TABLE_I = "tbl_069bc173c0f8"
TABLE_II = "tbl_9b11ef553745"
TABLE_III = "tbl_55ad1bdc7ac3"
_LEVEL = re.compile(r"^(\d[a-z]'?\[\d/2\])_\{?(\d)\}?(?:\^a)?$")
_VALUE = re.compile(r"^([\d.]+) ± ([\d.]+)$")
# Table II's product rows -> the file's id (Paschen notation, Table I).
PRODUCTS = {
    "6p[3/2]2": "2p6",
    "6p[3/2]1": "2p7",
    "6p[5/2]3": "2p8",
    "6p[5/2]2": "2p9",
    "6p[1/2]1": "2p10",
    "6s'[1/2]1": "1s2",
}
# Table I's levels (Racah) -> NIST energies (cm^-1) of the upper and lower levels.
NIST_LEVELS = {
    "6p[1/2]0": "80118.962",
    "6p[3/2]2": "79212.465",
    "6p[3/2]1": "78956.031",
    "6p[5/2]3": "78403.061",
    "6p[5/2]2": "78119.798",
    "6p[1/2]1": "77269.145",
    "6s[3/2]2": "67067.547",
    "6s[3/2]1": "68045.156",
}


def _level(cell: str) -> str:
    """Return a level label as the paper prints it, without the LaTeX subscript."""
    match = _LEVEL.match(cell.strip())
    if match is None:
        msg = f"unexpected level {cell!r}"
        raise ValueError(msg)
    return "".join(match.groups())


def _file(product: str) -> str:
    """Return the data file of one product channel."""
    return f"Xe/{SOURCE}__k_{product}_Xe.csv"


def extract() -> dict[str, list[Point]]:
    """Return Table II's rates, one file per product level, and the sums."""
    grid = json_table(DIRECTORY, TABLE_II)
    reactants = [_level(cell) for cell in grid[1][1:]]
    files: dict[str, list[Point]] = {_file(p): [] for p in PRODUCTS.values()}
    files[_file("sum")] = []
    for cells in grid[2:]:
        label = cells[0]
        product = "sum" if label.startswith("\\Sigma") else PRODUCTS.get(_level(label))
        if product is None:  # the row 6p[1/2]0: '...', 'b', 'b'
            if any(cell not in {"...", "b"} for cell in cells[1:]):
                msg = f"unexpected values in {cells}"
                raise ValueError(msg)
            continue
        for reactant, cell in zip(reactants, cells[1:], strict=True):
            if cell in {"...", "b"}:
                continue
            match = _VALUE.match(cell)
            if match is None:
                msg = f"unexpected value {cell!r}"
                raise ValueError(msg)
            printed = f"{reactant} -> {label}: {cell}"
            point = Point(
                reactant, reactant, number(match.group(1)), match.group(2), printed
            )
            files[_file(product)].append(point)
    return files


def _air(vacuum_nm: float) -> float:
    """Return the air wavelength (nm) of a vacuum wavelength (Edlén 1966)."""
    s2 = (1e3 / vacuum_nm) ** 2
    n = 1 + 8.34213e-5 + 2.406030e-2 / (130 - s2) + 1.5997e-4 / (38.9 - s2)
    return vacuum_nm / n


def _nist() -> dict[tuple[str, str], float]:
    """Return NIST's vacuum wavelengths (nm) by (upper, lower) energy string."""
    path = DATA / "external" / "nist_asd" / "Xe_I_lines_nir.tsv"
    lines = {}
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            clean = {k: (v or "").strip('"') for k, v in row.items() if k}
            wavelength = clean["obs_wl_vac(nm)"] or clean["ritz_wl_vac(nm)"]
            key = (clean["Ek(cm-1)"], clean["Ei(cm-1)"])
            lines[key] = float(wavelength.rstrip("+*"))
    return lines


def checks() -> None:
    """Check the counts, the sums, Table III's fractions and Table I's wavelengths."""
    files = extract()
    counts = [len(points) for points in files.values()]
    if counts != [2, 3, 3, 2, 1, 3, 3]:
        msg = f"counts {counts}"
        raise ValueError(msg)
    rates: dict[tuple[str, str], float] = {}
    for name, points in files.items():
        product = name.split("__k_")[1].removesuffix("_Xe.csv")
        for p in points:
            rates[p.row, product] = float(p.value)
    for reactant, tolerance in (
        ("6p[1/2]0", 0.05),
        ("6p[3/2]2", 0.5),
        ("6p[5/2]2", 0.5),
    ):
        total = sum(v for (r, q), v in rates.items() if r == reactant and q != "sum")
        if abs(total - rates[reactant, "sum"]) > tolerance:
            msg = f"sum for {reactant}: {total} against {rates[reactant, 'sum']}"
            raise ValueError(msg)
    grid = json_table(DIRECTORY, TABLE_III)
    for column, reactant in ((1, "6p[1/2]0"), (3, "6p[3/2]2")):
        six = {q: v for (r, q), v in rates.items() if r == reactant and q[0] == "2"}
        for cells in grid[3:]:
            product = PRODUCTS.get(_level(cells[0]))
            if product is None or cells[column] == "...":
                continue
            fraction = six.get(product, 0.0) / sum(six.values())
            if abs(fraction - float(cells[column])) > 0.006:
                msg = f"Table III {reactant} -> {product}: {fraction:.3f}"
                raise ValueError(msg)
    nist = _nist()
    upper = ""
    for cells in json_table(DIRECTORY, TABLE_I)[1:]:
        upper = _level(cells[1]) if cells[1] else upper
        lower = _level(cells[2])
        printed = float(cells[4].removesuffix("^{b}"))
        vacuum = nist[NIST_LEVELS[upper], NIST_LEVELS[lower]]
        if abs(_air(vacuum) - printed) > 0.05:
            msg = f"Table I {upper} -> {lower}: {printed} against NIST {vacuum}"
            raise ValueError(msg)


def main() -> None:
    """Check the tables and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
