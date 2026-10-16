# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Dragomir, McInerney & Nikogosyan (2002): TPA at 264 nm (Tables 1-7).

The paper prints the first author as "Dragonmir" (the library's directory name).

- Table 1: literature beta of fused silica (x 1e-11 cm/W) with wavelength, pulse
  length tau_p (ps), uncertainties and reference. The wavelength cell spans the rows of
  its wavelength (merged cells); the extraction gives it on the first row only and
  the script carries it down.
- Table 2: the six fused-silica samples (length, n, 1 - R, spectrophotometric alpha).
- Table 3: beta and alpha of each recorded curve of the fused-silica samples, and
  their averages (the glass and average cells span the glass's rows).
- Table 4: the crystal samples (quartz, sapphire, BBO in three configurations).
- Table 5: beta and alpha of each curve for quartz and sapphire, and the averages.
- Table 6: literature beta of crystals, in sections headed by a material row.
- Table 7: beta of each curve of BBO, per configuration, and the averages.
- Tables 8-10 are liquids and are not extracted.

Checked against the table images, all ten tables. Structural handling:

- Spanned cells (the wavelength of Table 1, the sample names and averages of Tables
  3, 5 and 7) are carried down to the rows they span; Table 6's material rows set the
  material of the rows below them.
- Literature values are split into one data file per reference and pulse length (a
  dataset has one pulse length). "<" cells are upper bounds and "~" an approximate
  value: the number is kept, the raw cell stays in ``printed``, and the record says what
  it is. Table 1's "35-52" (Ref. 16, 213 nm) is a range, not a value: it is not
  extracted (the record's notes give it). Thousands separators in tau_p ("25,000")
  are in ``printed`` only.
- Averages printed "2.4 \pm 0.3" are stored as value and sigma; "(5 \pm 9) \times
  10^{-3}" in e-notation (``5e-3``, ``9e-3``; scale 1). Table 4's alpha_sp column
  factor (x 1e-2) is the dataset scale.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Dragonmir_2002_FemtosecondMeasurementsTwo"
T1 = "tbl_208bbc898a3d"
T2 = "tbl_7dda4df6e106"
T3 = "tbl_2c36f1fa9b1c"
T4 = "tbl_1acb0895da28"
T5 = "tbl_83771c0b842b"
T6 = "tbl_bd4d7e6dd4aa"
T7 = "tbl_b99ed48cd500"
WAVELENGTH = "264"

# Table 1 (fused silica): (reference, tau_p) -> dataset stem.
T1_STEMS = {
    ("4", "15"): "gurzadyan_1992",
    ("5", "30"): "liu_1979",
    ("5, 6", "30"): "liu_1978_1979",
    ("7", "15"): "desalvo_1996",
    ("8", "0.28"): "hata_1990",
    ("9", "0.4"): "tomie_1989_0p4ps",
    ("9", "8"): "tomie_1989_8ps",
    ("10", "0.7"): "taylor_1988",
    ("11", "0.45"): "simon_1989",
    ("12", "4"): "ross_1990",
    ("13", "25,000"): "eva_1996",
    ("14", "9000"): "brimacombe_1989_9ns",
    ("14", "5000"): "brimacombe_1989_5ns",
    ("15", "9000"): "taylor_1988b_9ns",
    ("15", "10,000"): "taylor_1988b_10ns",
    ("17", "26"): "repeyev_1992",
    ("18", "0.8"): "kittelmann_1994",
}
# Table 6 (crystals): material row -> material id, and (reference, tau_p) -> stem.
T6_MATERIALS = {"Crystalline Quartz": "SiO2_quartz", "Sapphire": "Al2O3", "BBO": "BBO"}
T6_STEMS = {
    ("4", "15"): "gurzadyan_1992",
    ("5, 6", "30"): "liu_1978_1979",
    ("13", "25,000"): "eva_1996",
    ("16", "?"): "schenker_1995",
    ("7", "15"): "desalvo_1996",
    ("6", "30"): "liu_1978",
    ("22", "0.2"): "veitas_1997",
    ("34", "0.8"): "dubietis_2000",
}
# The fused-silica samples (Tables 2 and 3): printed name -> grade slug.
GLASSES = {
    "Fused silica KU-1": "ku1",
    "Corning 7940": "corning7940",
    "Fused silica SQ": "sq",
    "Suprasil": "suprasil",
    "Herasil": "herasil",
    "Infrasil": "infrasil",
}
# The crystal samples (Tables 4, 5 and 7): printed name -> (material, slug).
CRYSTALS = {
    "Crystalline quartz": ("SiO2_quartz", "SiO2_quartz"),
    "Sapphire": ("Al2O3", "Al2O3"),
    "BBO (\\parallel c)": ("BBO", "BBO_z_o"),
    "BBO (\\perp c, ordinary beam)": ("BBO", "BBO_x_o"),
    "BBO (\\perp c, extraordinary beam)": ("BBO", "BBO_x_e"),
    "\\parallel c": ("BBO", "BBO_z_o"),
    "\\perp c, ordinary beam": ("BBO", "BBO_x_o"),
    "\\perp c, extraordinary beam": ("BBO", "BBO_x_e"),
}

_LITERATURE = re.compile(r"^(?:<|~)?(?P<value>[\d.]+)$")
_PLUS_MINUS = re.compile(
    r"^\(?(?P<value>[\d.]+)\s*(?:±|\\pm)\s*(?P<sigma>[\d.]+)\)?"
    r"(?:\s*\\times\s*10\^\{(?P<power>-\d+)\})?$"
)


def literature(cell: str) -> str:
    """Return the number of a literature cell ("<" bound, "~" approximate)."""
    match = _LITERATURE.match(cell.strip())
    if match is None:
        msg = f"not a literature value: {cell!r}"
        raise ValueError(msg)
    return match["value"]


def plus_minus(cell: str) -> tuple[str, str]:
    r"""Return value, sigma of ``v ± s``, ``(v ± s)`` or ``(v \pm s) \times 10^{k}``."""
    match = _PLUS_MINUS.match(cell.strip())
    if match is None:
        msg = f"not a value with uncertainty: {cell!r}"
        raise ValueError(msg)
    power = f"e{match['power']}" if match["power"] else ""
    return match["value"] + power, match["sigma"] + power


def _literature_tables(files: dict[str, list[Point]]) -> None:
    """Add the literature values of Tables 1 (fused silica) and 6 (crystals)."""
    wavelength = ""
    for cells in json_table(DIRECTORY, T1)[1:]:
        wavelength = cells[0] or wavelength  # the spanned wavelength cell
        if cells[6] == "16":
            continue  # "35-52": a range, not a value (the record's notes)
        stem = T1_STEMS[cells[6], cells[1]]
        files.setdefault(f"SiO2/dragomir_2002__{stem}_SiO2.csv", []).append(
            Point(wavelength, wavelength, literature(cells[2]), "", " | ".join(cells))
        )
    material = ""
    for cells in json_table(DIRECTORY, T6)[1:]:
        if cells[0] in T6_MATERIALS:  # a section row: the material of the rows below
            material = T6_MATERIALS[cells[0]]
            continue
        stem = T6_STEMS[cells[6], cells[1]]
        files.setdefault(f"{material}/dragomir_2002__{stem}_{material}.csv", []).append(
            Point(cells[0], cells[0], literature(cells[2]), "", " | ".join(cells))
        )


def _runs(
    files: dict[str, list[Point]],
    grid: list[list[str]],
    names: dict[str, tuple[str, str]],
    alpha: bool,
) -> None:
    """Add the per-curve values and averages of Tables 3, 5 and 7."""
    material = slug = ""
    run = 0
    for cells in grid[1:]:
        printed = " | ".join(cells)
        if cells[0]:  # the first row of a sample: its name and averages
            material, slug = names[cells[0]]
            run = 0
            value, sigma = plus_minus(cells[2])
            files[f"{material}/dragomir_2002__beta_{slug}.csv"] = [
                Point("mean", WAVELENGTH, value, sigma, printed)
            ]
            if alpha and cells[4]:
                value, sigma = plus_minus(cells[4])
                files[f"{material}/dragomir_2002__alpha_fit_{slug}.csv"] = [
                    Point("mean", WAVELENGTH, value, sigma, printed)
                ]
        run += 1
        files.setdefault(f"{material}/dragomir_2002__beta_runs_{slug}.csv", []).append(
            Point(str(run), WAVELENGTH, number(cells[1]), "", printed)
        )
        if alpha and cells[3]:
            files.setdefault(
                f"{material}/dragomir_2002__alpha_fit_runs_{slug}.csv", []
            ).append(Point(str(run), WAVELENGTH, number(cells[3]), "", printed))


def extract() -> dict[str, list[Point]]:
    """Return every solid's values of Tables 1-7."""
    files: dict[str, list[Point]] = {}
    _literature_tables(files)

    # Table 2: the fused-silica samples.
    for cells in table(DIRECTORY, T2)[1:]:
        slug = GLASSES[cells[0]]
        printed = " | ".join(cells)
        files[f"SiO2/dragomir_2002__index_{slug}.csv"] = [
            Point(slug, WAVELENGTH, cells[2], "", printed)
        ]
        alpha_sp = cells[4].removeprefix("≤")  # "≤0.01" is a bound (the record)
        files[f"SiO2/dragomir_2002__alpha_sp_{slug}.csv"] = [
            Point(slug, WAVELENGTH, number(alpha_sp), "", printed)
        ]

    # Table 4: the crystal samples; BBO's two ordinary-beam rows share a dataset.
    for cells in table(DIRECTORY, T4)[1:]:
        material, slug = CRYSTALS[cells[0]]
        axis_slug = {"BBO_z_o": "BBO_o", "BBO_x_o": "BBO_o", "BBO_x_e": "BBO_e"}
        stem = axis_slug.get(slug, slug)
        row = {"BBO_z_o": "z-cut", "BBO_x_o": "x-cut", "BBO_x_e": "x-cut"}.get(
            slug, slug
        )
        printed = " | ".join(cells)
        files.setdefault(f"{material}/dragomir_2002__index_{stem}.csv", []).append(
            Point(row, WAVELENGTH, cells[2], "", printed)
        )
        value, sigma = plus_minus(cells[4])
        files.setdefault(f"{material}/dragomir_2002__alpha_sp_{stem}.csv", []).append(
            Point(row, WAVELENGTH, value, sigma, printed)
        )

    glasses = {name: ("SiO2", slug) for name, slug in GLASSES.items()}
    _runs(files, json_table(DIRECTORY, T3), glasses, alpha=True)
    _runs(files, json_table(DIRECTORY, T5), CRYSTALS, alpha=True)
    _runs(files, json_table(DIRECTORY, T7), CRYSTALS, alpha=False)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
