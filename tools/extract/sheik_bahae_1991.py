# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Sheik-Bahae et al. (1991), Table III: n2 of 40 rows of solids.

Table III (``tbl_1be7b7bfe319``, p. 1305) gives, per row, the material, the wavelength
(um), the band gap (eV) and the linear index used in the calculation (from the CRC
Handbook and Palik's Handbook, Refs. [41], [42]), the measured n2 and two predicted n2
(the two-band theory with K' = 1.50e-8 and, marked *, with K' = 0.86e-8 fitted to the
wide-gap solids), all in units of 1e-13 esu. The cells agree with the table image.

Structural facts the extraction loses, from the image (p. 1305): the table is printed in
four blocks separated by blank lines, and the blocks are the provenance of the measured
column (text p. 1304 and Fig. 5):

1. rows 1-14 (Ge to BaF2 at 0.53 um): the authors' own Z-scan measurements (the open
   squares of Fig. 5(a); Ge at 10.6 um with a picosecond CO2 laser, p. 1305);
2. rows 15-21 (AlGaAs, 0.850-0.810 um): LaGasse et al. [40] (the triangles of Fig.
   5(b));
3. rows 22-39 (CdS at 1.06 um to LiF): Adair et al. [6] (``adair_1989``; the solid
   circles of Fig. 5(a)); each value is Adair's Table IV value rounded to two figures
   (checked against ``adair_1989``'s files): KDP 0.7 and KH2PO4 0.8 are Adair's KDP (o)
   0.72 and (e) 0.78, and "SiO2" 1.1 is Adair's crystal quartz (o) 1.12, not his fused
   silica 0.85, so that row is ``SiO2_quartz``;
4. row 40 (SiO2 at 0.249 um): Ross et al. [39] (the diamond of Fig. 5(a)).

Rows of materials that are not registered (CdTe, CdSe, CdS0.5Se0.5, SBN of unstated
composition, AlGaAs of unstated composition) are not extracted; their values are in the
record's notes. The band gap cell "0.87^{\dagger}" (Ge) carries the dagger "direct gap
rather than the fundamental absorption edge": the dagger is dropped from the number.
"CaCo_3" is printed so (CaCO3). The paper prints no uncertainty for any cell.
"""

from _library import Point, number, table, write

DIRECTORY = "SheikBahae_1991_DispersionBoundElectron"
TABLE_III = "tbl_1be7b7bfe319"

# The printed material cell -> (registered id, dataset name). None: not registered.
# KDP and KH_2PO_4 are the same crystal (Adair's o and e rows), kept apart by name;
# the 1.06-um "SiO_2" row is crystal quartz (see the module docstring).
MATERIALS: dict[str, tuple[str, str] | None] = {
    "Ge": ("Ge", "Ge"),
    "GaAs": ("GaAs", "GaAs"),
    "CdTe": None,
    "CdSe": None,
    "CdS_{0.5}Se_{0.5}": None,
    "ZnTe": ("ZnTe", "ZnTe"),
    "CdS": ("CdS", "CdS"),
    "ZnSe": ("ZnSe", "ZnSe"),
    "SBN": None,
    "ZnS": ("ZnS", "ZnS"),
    "KTP": ("KTP", "KTP"),
    "BaF_2": ("BaF2", "BaF2"),
    "AlGaAs": None,
    "AgCl": ("AgCl", "AgCl"),
    "ZnO": ("ZnO", "ZnO"),
    "NaBr": ("NaBr", "NaBr"),
    "CaCo_3": ("CaCO3", "CaCO3"),
    "KBr": ("KBr", "KBr"),
    "KCl": ("KCl", "KCl"),
    "KDP": ("KDP", "KDP"),
    "KH_2PO_4": ("KDP", "KH2PO4"),
    "NaCl": ("NaCl", "NaCl"),
    "Al_2O_3": ("Al2O3", "Al2O3"),
    "KF": ("KF", "KF"),
    "MgO": ("MgO", "MgO"),
    "SrF_2": ("SrF2", "SrF2"),
    "CaF_2": ("CaF2", "CaF2"),
    "MgF_2": ("MgF2", "MgF2"),
    "LiF": ("LiF", "LiF"),
}

# The measured column's provenance by row number (1-based, the module docstring):
# the dataset-name prefix of each block.
OWN, LAGASSE, ADAIR, ROSS = "n2", "n2_lagasse", "n2_adair", "n2_ross"


def _block(row: int) -> str:
    """Return the dataset prefix of a row's measured n2 (the printed blocks)."""
    if row <= 14:
        return OWN
    if row <= 21:
        return LAGASSE
    if row <= 39:
        return ADAIR
    return ROSS


def _identity(row: int, cell: str) -> tuple[str, str] | None:
    """Return (material id, dataset name) of a row, or None if not registered."""
    if cell == "SiO_2":
        return ("SiO2_quartz", "SiO2_quartz") if row < 40 else ("SiO2", "SiO2")
    return MATERIALS[cell]


def _add(files: dict[str, list[Point]], name: str, material: str, point: Point) -> None:
    files.setdefault(f"{material}/sheik_bahae_1991__{name}.csv", []).append(point)


def extract() -> dict[str, list[Point]]:
    """Return the measured and predicted n2, band gaps and indices of Table III."""
    files: dict[str, list[Point]] = {}
    rows = table(DIRECTORY, TABLE_III)
    seen_gap: set[str] = set()
    for index, cells in enumerate(rows[1:], start=1):
        identity = _identity(index, cells[0])
        if identity is None:
            continue
        material, name = identity
        printed = " | ".join(cells)
        row, wavelength = str(index), number(cells[1])
        _add(
            files,
            f"{_block(index)}_{name}",
            material,
            Point(row, wavelength, number(cells[4]), "", printed),
        )
        _add(
            files,
            f"n2_pred_{name}",
            material,
            Point(row, wavelength, number(cells[5]), "", printed),
        )
        if cells[6].strip():
            _add(
                files,
                f"n2_pred_wide_{name}",
                material,
                Point(row, wavelength, number(cells[6]), "", printed),
            )
        _add(
            files,
            f"n0_{name}",
            material,
            Point(row, wavelength, number(cells[3]), "", printed),
        )
        # One gap per dataset name: repeated rows (ZnSe, BaF2, CdS) print the same.
        if name not in seen_gap:
            seen_gap.add(name)
            gap = cells[2].replace("^{\\dagger}", "")
            point = Point(row, "", number(gap), "", printed)
            _add(files, f"eg_{name}", material, point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
