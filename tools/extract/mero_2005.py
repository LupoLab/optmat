# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Mero et al. (2005): breakdown figures of merit of oxide films.

- Table I (``tbl_2091d43cd8a8``, p. 115109-2): per film, the index n0 at 800 nm, the
  thickness D (nm; no quantity, kept in ``printed``) and the band gap E_g (eV).
- Table II (``tbl_e839fc341559``, p. 115109-3): per film, the multiphoton order m, the
  avalanche coefficient alpha (cm^2/J), the multiphoton coefficient beta_m (cm^(2m-3)
  fs^(m-1)/J^m), the effective decay time T (fs) and the reduced mass m_r/m0 obtained
  from beta_m with Eq. (2).

Structural fix (Table II, image checked): the extraction reads each beta_m cell
"1.5 x 10^23" as "1.531023": the multiplication sign became "3" and the exponent lost
its superscript. The script splits the cell into mantissa (one decimal), the "3" and
"10", and the exponent (two digits), and checks the pattern; the values then agree
with the image (1.5e23, 6.7e24, 2.9e25, 2.3e24, 9.9e25). Table I's header units "snmd"
and "seVd" are the extraction's rendering of "(nm)" and "(eV)" (image). Table III
(kappa of the power law) has no quantity; the record gives it in the formula.

Every other cell agrees with the table images.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Mero_2005_ScalingLawsFemtosecond"
TABLE_I = "tbl_2091d43cd8a8"
TABLE_II = "tbl_e839fc341559"
WAVELENGTH_NM = "800"

MATERIALS = {
    "TiO2": "TiO2_amorphous",
    "Ta2O5": "Ta2O5",
    "HfO2": "HfO2",
    "Al2O3": "Al2O3_amorphous",
    "SiO2": "SiO2",
}

# "1.531023" = 1.5 x 10^23: mantissa, the multiplication sign read as "3", 10, exponent.
_GARBLED_POWER = re.compile(r"^(\d\.\d)310(\d\d)$")


def _power(cell: str) -> str:
    """Return a garbled ``m x 10^k`` cell of Table II as ``mek``."""
    match = _GARBLED_POWER.match(cell.strip())
    if match is None:
        msg = f"unexpected beta_m cell {cell!r}"
        raise ValueError(msg)
    return f"{match.group(1)}e{match.group(2)}"


def extract() -> dict[str, list[Point]]:
    """Return n0 and E_g (Table I) and alpha, beta_m, T and m_r (Table II) per film."""
    files: dict[str, list[Point]] = {}
    rows = table(DIRECTORY, TABLE_I)
    if rows[0] != ["Material", "n0", "D snmd", "Eg seVd"]:
        msg = f"unexpected Table I header {rows[0]!r}"
        raise ValueError(msg)
    for cells in rows[1:]:
        material = MATERIALS[cells[0]]
        printed = " | ".join(cells)
        files[f"{material}/mero_2005__n_{material}.csv"] = [
            Point("1", WAVELENGTH_NM, number(cells[1]), "", printed)
        ]
        files[f"{material}/mero_2005__band_gap_{material}.csv"] = [
            Point("1", "", number(cells[3]), "", printed)
        ]
    rows = table(DIRECTORY, TABLE_II)
    if rows[0] != ["Material", "m", "a", "bmfit", "T", "mrMPA/m0"]:
        msg = f"unexpected Table II header {rows[0]!r}"
        raise ValueError(msg)
    for cells in rows[1:]:
        material = MATERIALS[cells[0]]
        printed = " | ".join(cells)
        files[f"{material}/mero_2005__avalanche_{material}.csv"] = [
            Point("1", WAVELENGTH_NM, number(cells[2]), "", printed)
        ]
        files[f"{material}/mero_2005__beta_m_{material}.csv"] = [
            Point("1", WAVELENGTH_NM, _power(cells[3]), "", printed)
        ]
        files[f"{material}/mero_2005__decay_time_{material}.csv"] = [
            Point("1", "", number(cells[4]), "", printed)
        ]
        files[f"{material}/mero_2005__reduced_mass_{material}.csv"] = [
            Point("1", "", number(cells[5]), "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
