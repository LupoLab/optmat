# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Patankar et al. (2017): 2PA coefficients and linear losses at 213 nm.

Table 1 (``tbl_c8b5daab8ce4``, p. 8310) lists the samples (manufacturer, thickness), the
band gap with its reference ("13.6 [10]") and n0 at 213 nm; Table 2
(``tbl_8426527ad90f``, p. 8311) the linear loss alpha (cm^-1) and the 2PA coefficient
beta (cm/GW, "v ± s"). Both agree with their images, row for row; the fused-silica rows
("FS—1AK" ...) are grades of one material. LiF's beta is printed "<0.001", an upper
limit, and is not extracted (the record's notes quote it). The band gap and n0 are the
same on the four fused-silica rows ("Band gap is assumed to be the same for different
grades", Table 1 note); they are taken once, from the FS—1AK row, which carries the
gap's reference.
"""

from _library import Point, number, table, write

DIRECTORY = "Patankar_2017_TwoPhotonAbsorption"
TABLE_1 = "tbl_c8b5daab8ce4"
TABLE_2 = "tbl_8426527ad90f"

# The printed material cell -> (material id, dataset suffix).
SAMPLES = {
    "LiF": ("LiF", "LiF"),
    "MgF_2": ("MgF2", "MgF2"),
    "CaF_2": ("CaF2", "CaF2"),
    "BaF_2": ("BaF2", "BaF2"),
    "Sapphire": ("Al2O3", "Al2O3"),
    "FS—1AK": ("SiO2", "SiO2_1AK"),
    "FS—7980 (ArF)": ("SiO2", "SiO2_7980"),
    "FS—8655": ("SiO2", "SiO2_8655"),
    "FS—312": ("SiO2", "SiO2_312"),
}

# Fused-silica rows after the first repeat its gap and index (Table 1 note).
REPEATED = {"FS—7980 (ArF)", "FS—8655", "FS—312"}

WAVELENGTH = "213"


def extract() -> dict[str, list[Point]]:
    """Return beta, alpha, the band gaps and the indices at 213 nm."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        material, name = SAMPLES[cells[0]]
        printed = " | ".join(cells)
        stem = f"{material}/patankar_2017__"
        files[f"{stem}alpha_{name}.csv"] = [
            Point(WAVELENGTH, WAVELENGTH, number(cells[3]), "", printed)
        ]
        beta = cells[4].strip()
        if not beta.startswith("<"):
            value, _, sigma = beta.partition("±")
            files[f"{stem}beta_{name}.csv"] = [
                Point(WAVELENGTH, WAVELENGTH, number(value), number(sigma), printed)
            ]
    for cells in table(DIRECTORY, TABLE_1)[1:]:
        if cells[0] in REPEATED:
            continue
        material, _ = SAMPLES[cells[0]]
        printed = " | ".join(cells)
        stem = f"{material}/patankar_2017__"
        gap = cells[3].split("[")[0]  # "13.6 [10]": the reference is in the record
        files[f"{stem}eg_{material}.csv"] = [
            Point(material, "", number(gap), "", printed)
        ]
        files[f"{stem}n0_{material}.csv"] = [
            Point(WAVELENGTH, WAVELENGTH, number(cells[4]), "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
