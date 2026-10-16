# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Slattery & Nikogosyan (2003): two-photon absorption at 211 nm (Tables 1, 2).

Table 1 ("Parameters of samples at lambda = 211 nm") gives, per sample, the length,
the refractive index (calculated from a handbook's dispersion relations), the Fresnel
transmission (1 - R)^2 in per cent and the linear absorption coefficient alpha with
its uncertainty ("0.041 \pm 0.001" is stored as value and sigma). Table 2 gives, per
material, the band gap with its citation, and beta (x 1e-11 cm/W) with its uncertainty
at several wavelengths: this work at 211 nm, and literature values (reproduced), each
with the pulse length tau_p and the reference.

Checked against the table images. Table 1 (tbl_43ed8a53704c.jpg): the material names
wrap onto two lines in print ("Suprasil / (sample#1)", "Crystalline / quartz (||c)"),
which the extraction joins without a space; every number agrees. Table 2
(tbl_0a42752ea230.jpg, published only as JSON): the material and E_g cells span the
rows of their material; the extraction gives them on the first row only, so the
script carries them down. The superscript "a" (only a relative uncertainty is given
in the original) is extracted as a trailing "a" on the uncertainty and stripped; the
record says which values carry it.
"""

import re

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Slattery_2003_TwoPhotonAbsorption"
SAMPLES = "tbl_43ed8a53704c"
BETAS = "tbl_0a42752ea230"

# Table 1 rows: (material id, row id).
SAMPLE_ROWS = {
    "Suprasil(sample#1)": ("SiO2", "sample 1"),
    "Suprasil(sample#2)": ("SiO2", "sample 2"),
    "Crystallinequartz (\\|c)": ("SiO2_quartz", "SiO2_quartz"),
    "Barium fluoride": ("BaF2", "BaF2"),
    "Calcium fluoride": ("CaF2", "CaF2"),
}
MATERIALS = {
    "Fused silica": "SiO2",
    "Crystalline quartz (\\|c)": "SiO2_quartz",
    "Barium fluoride": "BaF2",
    "Calcium fluoride": "CaF2",
}
# Table 2's reference cells and the dataset stem of each.
REFERENCES = {
    "This work": "beta_211",
    "[17]": "kittelmann_1994",
    "[18]": "schenker_1995",
    "[19]": "gurzadyan_1992",
    "[20]": "taylor_1988",
    "[1]": "dragomir_2002",
}

_PLUS_MINUS = re.compile(r"^(?P<value>[\d.]+)\s*±\s*(?P<sigma>[\d.]+)a?$")
_GAP = re.compile(r"^(?:≈)?(?P<value>[\d.]+) \[\d+\]$")


def plus_minus(cell: str) -> tuple[str, str]:
    """Return the value and uncertainty of ``v ± s`` (footnote marker "a" dropped)."""
    match = _PLUS_MINUS.match(cell.strip())
    if match is None:
        msg = f"not a value with uncertainty: {cell!r}"
        raise ValueError(msg)
    return match["value"], match["sigma"]


def extract() -> dict[str, list[Point]]:
    """Return the sample indices and alphas, the band gaps and the betas."""
    files: dict[str, list[Point]] = {}

    def add(name: str, point: Point) -> None:
        files.setdefault(name, []).append(point)

    samples = table(DIRECTORY, SAMPLES)
    if samples[0] != ["Material", "l (cm)", "n", "(1 - R)^2", "\\alpha (cm-1)"]:
        msg = f"unexpected Table 1 header {samples[0]!r}"
        raise ValueError(msg)
    for cells in samples[1:]:
        material, row = SAMPLE_ROWS[cells[0]]
        printed = " | ".join(cells)
        add(
            f"{material}/slattery_2003__index_{material}.csv",
            Point(row, "211", cells[2], "", printed),
        )
        value, sigma = plus_minus(cells[4])
        add(
            f"{material}/slattery_2003__alpha_{material}.csv",
            Point(row, "211", value, sigma, printed),
        )

    grid = json_table(DIRECTORY, BETAS)
    if grid[0][4] != "β \N{MULTIPLICATION SIGN} 10^{-11} (cm/W)":
        msg = f"unexpected Table 2 header {grid[0]!r}"
        raise ValueError(msg)
    material = ""
    for cells in grid[1:]:
        # The material and E_g cells span their material's rows (see the docstring).
        if cells[0]:
            material = MATERIALS[cells[0]]
            gap = _GAP.match(cells[1])
            if gap is None:
                msg = f"not a band-gap cell: {cells[1]!r}"
                raise ValueError(msg)
            add(
                f"{material}/slattery_2003__band_gap_{material}.csv",
                Point(material, "", gap["value"], "", " | ".join(cells)),
            )
        value, sigma = plus_minus(cells[4])
        stem = REFERENCES[cells[5]]
        add(
            f"{material}/slattery_2003__{stem}_{material}.csv",
            Point(cells[2], cells[2], value, sigma, " | ".join(cells)),
        )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
