# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract DeSalvo et al. (1996): band gaps, n2, 2PA coefficients and linear indices.

Table I (``tbl_9572f7715e3f``, published only as JSON) gives, per material, the band
gap E_g, the effective gap E'_g, n2 in units of 1e-14 esu at 1064/532/355/266 nm and
the 2PA coefficient beta (cm/GW) at 532/355/266 nm. Table II (``tbl_87c9f8dcbfd0``)
gives the linear index used at each wavelength. The cells were checked against both
table images: they agree, including the quoted cells "6.8" and "8.0" of E'_g.

The n2 columns (n = n0 + n2<E^2>, in units of 1e-14 esu) are nonlinear_index_field,
one file per material; KTP is split by axis, since its 1064- and 532-nm values come
from two samples (E || y and E || z). "NM" (not measured) and "≈0" (no measurable 2PA)
cells are not numbers and are left out.

The pulse duration of each harmonic (text, p. 1326) is attached to each n2 and beta as
the per-row column ``pulse_fs``.
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "DeSalvo_1996_InfraredUltravioletMeasurements"
TABLE_I = "tbl_9572f7715e3f"
TABLE_II = "tbl_87c9f8dcbfd0"

MATERIALS = {
    "LiF": "LiF",
    "MgF_2": "MgF2",
    "BaF_2": "BaF2",
    "SiO_2": "SiO2",
    "Al_2O_3": "Al2O3",
    "BBO": "BBO",
    "KBr": "KBr",
    "CaCO_3": "CaCO3",
    "LiNbO_3": "LiNbO3",
    "KTP": "KTP",
    "Te Glass": "tellurite_75TeO2-20ZnO-5Na2O",
}

# FWHM pulse durations at each harmonic (p. 1326): 30, 22, 17 and 15 ps.
PULSE_FS = {"1064": "30000", "532": "22000", "355": "17000", "266": "15000"}

# Table I's n2 columns (index -> wavelength in nm).
N2_COLUMNS = {3: "1064", 4: "532", 5: "355", 6: "266"}

# Table I's beta columns (index -> wavelength in nm).
BETA_COLUMNS = {7: "532", 8: "355", 9: "266"}

# KTP: two samples (footnote 2, p. 1325), x cut with E || z at 532 nm and z cut with
# E || y at 1064 nm, so its n2 and Table II indices split by axis.
KTP_AXIS = {"1064": "y", "532": "z"}


def _per_material(
    material: str, name: str, points: list[Point]
) -> dict[str, list[Point]]:
    """Return one file for the material, or one per KTP sample axis."""
    if material != "KTP":
        return {f"{material}/desalvo_1996__{name}_{material}.csv": points}
    return {
        f"KTP/desalvo_1996__{name}_KTP_{axis}.csv": [p for p in points if p.x == x]
        for x, axis in KTP_AXIS.items()
    }


def _harmonics(cells: list[str], columns: dict[int, str]) -> list[Point]:
    """Return the "v±s" cells of the given columns, with their pulse durations."""
    printed = " | ".join(cells)
    points = []
    for column, wavelength in columns.items():
        parsed = _value(cells[column])
        if parsed is not None:
            extra = {"pulse_fs": PULSE_FS[wavelength]}
            points.append(Point(wavelength, wavelength, *parsed, printed, extra))
    return points


def _value(cell: str) -> tuple[str, str] | None:
    """Return (value, sigma) of a "v±s" cell, or None for "NM" and "≈0"."""
    text = cell.strip()
    if text in {"NM", "≈0", ""}:
        return None
    value, _, sigma = text.partition("±")
    return number(value), number(sigma) if sigma else ""


def extract() -> dict[str, list[Point]]:
    """Return the band gaps, n2, the 2PA coefficients and the linear indices."""
    files: dict[str, list[Point]] = {}
    rows = json_table(DIRECTORY, TABLE_I)
    for cells in rows[2:]:
        material = MATERIALS[cells[0]]
        printed = " | ".join(cells)
        files[f"{material}/desalvo_1996__eg_{material}.csv"] = [
            Point("1", "E_g", number(cells[1]), "", printed)
        ]
        # E'_g of BBO and calcite are printed in quotation marks ("6.8", "8.0"):
        # the further-adjusted gaps of the text (pp. 1329-1330), kept as values.
        prime = cells[2].strip().strip('“”"')
        if prime != "NM":
            files[f"{material}/desalvo_1996__eg_prime_{material}.csv"] = [
                Point("1", "E'_g", number(prime), "", printed)
            ]
        files.update(_per_material(material, "n2", _harmonics(cells, N2_COLUMNS)))
        betas = _harmonics(cells, BETA_COLUMNS)
        if betas:
            files[f"{material}/desalvo_1996__beta_{material}.csv"] = betas
    indices = table(DIRECTORY, TABLE_II)
    wavelengths = [cell.removesuffix(" nm") for cell in indices[0][1:]]
    for cells in indices[1:]:
        material = MATERIALS[cells[0]]
        printed = " | ".join(cells)
        points = [
            Point(wavelength, wavelength, number(cell), "", printed)
            for wavelength, cell in zip(wavelengths, cells[1:], strict=True)
            if cell.strip()
        ]
        files.update(_per_material(material, "n", points))
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
