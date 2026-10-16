# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bergé, Mauger & Skupin (2010): the silica parameters of their simulations.

Table I (``tbl_3463a4cb91f2``), "Physical parameters for silica at 355 nm and 1550 nm",
checked against the table image: the cells agree. Every entry is a model input. The
column factors of the row labels (``n_2 (10-16cm2/W)``, ``\sigma (10-19cm2)``) are the
datasets' scales. The row ``k''' (fs2/cm)`` prints the unit of k'' again (the image
confirms "fs2/cm"): a third-order dispersion is in fs^3/cm, so it is recorded in fs3/m
with scale 100, as the record explains.
"""

from _library import Point, number, table, write

DIRECTORY = "Berge_2010_MultifilamentationPowerfulOptical"
TABLE = "tbl_3463a4cb91f2"
WAVELENGTHS = ("355", "1550")

# Row label as extracted -> dataset stem; one file holds both wavelengths.
BOTH = {
    "k'' (fs2/cm)": "gvd",
    "k''' (fs2/cm)": "tod",
    "n_0": "refractive_index",
    "n_2 (10-16cm2/W)": "n2",
    "P_{\\text{cr}} (MW)": "critical_power",
    "\\sigma (10-19cm2)": "avalanche_cross_section",
}
# Rows of quantities with no wavelength abscissa: one file per wavelength.
EACH = {
    "g": "raman_fraction",
    "\\tau_R (fs)": "raman_tau_R",
    "\\omega_R^{-1} (fs)": "raman_inverse_omega_R",
}


def extract() -> dict[str, list[Point]]:
    """Return Table I's rows as data files."""
    rows = table(DIRECTORY, TABLE)
    if rows[0][1:] != [f"{w} nm" for w in WAVELENGTHS]:
        msg = f"unexpected column heads {rows[0]}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    seen = set()
    for cells in rows[1:]:
        label = cells[0]
        printed = " | ".join(cells)
        if label in BOTH:
            files[f"SiO2/berge_2010__{BOTH[label]}.csv"] = [
                Point(w, w, number(cell), "", printed)
                for w, cell in zip(WAVELENGTHS, cells[1:], strict=True)
            ]
        elif label in EACH:
            for w, cell in zip(WAVELENGTHS, cells[1:], strict=True):
                files[f"SiO2/berge_2010__{EACH[label]}_{w}.csv"] = [
                    Point("1", "", number(cell), "", printed)
                ]
        else:
            msg = f"unexpected row {label!r}"
            raise ValueError(msg)
        seen.add(label)
    if seen != set(BOTH) | set(EACH):
        msg = f"rows missing: {sorted(set(BOTH) | set(EACH) - seen)}"
        raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
