# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bishop, Kirtman, Kurtz & Rice (1993): mean vibrational gamma of CO2.

Table VI (``tbl_bcd9303581c8``, CSV; image checked, every cell agrees): the isotropic
pure vibrational gamma^v of CO2 for the static case, the dc Kerr effect, ESHG, THG and
the intensity-dependent refractive index (IDRI), and the Kerr combination gamma_K^v, at
hbar omega = 0.07 E_h, by SCF and MP2 without the second derivatives of beta (columns
"SCF^a", "MP2^a") and by complete SCF ("SCF"). One file per row (process); the row ids
are the three columns. H2O and NH3 (Tables II-IV, VII-X) are not curated.
"""

from _library import Point, table, write

DIRECTORY = "Bishop_1993_CalculationVibrationalDynamic"
TABLE_VI = "tbl_bcd9303581c8"
ROWS = {  # printed row label -> (file name, wavelength in hartree)
    "\\overline{\\gamma}^{v}(0;0,0,0)": ("static", "0"),
    "\\overline{\\gamma}^{v}(-\\omega;\\omega,0,0)": ("dc_kerr", "0.07"),
    "\\overline{\\gamma}^{v}(-2\\omega;\\omega,\\omega,0)": ("eshg", "0.07"),
    "\\overline{\\gamma}^{v}(-3\\omega;\\omega,\\omega,\\omega)": ("thg", "0.07"),
    "\\overline{\\gamma}^{v}(-\\omega;\\omega,\\omega,-\\omega)": ("kerr", "0.07"),
    "\\overline{\\gamma}_{K}^{v}(-\\omega;\\omega,0,0)": ("dc_kerr_k", "0.07"),
}
METHODS = ("scf_partial", "mp2_partial", "scf")


def extract() -> dict[str, list[Point]]:
    """Return Table VI, one file per process."""
    grid = table(DIRECTORY, TABLE_VI)
    if grid[0][1:] != ["SCFa", "MP2a", "SCF"]:
        msg = f"unexpected Table VI header: {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        name, omega = ROWS[cells[0]]
        points = []
        for method, value in zip(METHODS, cells[1:], strict=True):
            float(value)
            points.append(Point(method, omega, value, "", " | ".join(cells)))
        files[f"CO2/bishop_1993__{name}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
