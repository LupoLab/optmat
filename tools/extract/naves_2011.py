# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Naves et al. (2011): O3 electronic mean polarizability, static and dynamic.

Table III (``tbl_a5631d2a2b1a``, cells file only; image checked, every cell agrees):
electronic contributions in atomic units at omega = 0, 0.0239, 0.0428 and 0.0656
hartree, CCSD at the CCSD-optimised and at the experimental geometry and CCSD(T) at the
experimental geometry (the dynamic CCSD(T) values scaled by Eq. (29)). Only the mean
polarizability rows (the first four, under "\alpha") are taken; the beta and gamma
rows (hyperpolarizabilities) are for a later phase.

Table II (``tbl_6e45fd86f96d``, cells file only; image checked): static components and
means at the experimental geometry by four methods. Only the orbital-relaxed CCSD mean
(the row of the mean, "alpha" with a macron, the third value column) is taken: the
analytical and orbital-unrelaxed finite-field means repeat Table III's CCSD value
(19.77) and the CCSD(T) mean is Table III's static CCSD(T) value (17.12).

Third-order data:
- Table III's gamma rows (static, THG, dc-SHG, IDRI, dc-K at omega = 0, 0.0239,
  0.0428, 0.0656 hartree), one file per process and method column: ``gamma_<process>_
  <method>``, rows omega.
- Tables I (``tbl_4dfc86410542``, HF and CCSD in aug-, d-aug-, t-aug-cc-pVTZ) and II
  (four methods), the static gamma components and mean: one file per row
  (``gamma_basis_<component>``, ``gamma_method_<component>``), rows the columns.
- Table VI (``tbl_cb8015ae109a``), the pure vibrational terms. STRUCTURAL FIX: the
  library merged each block's four omegas (0.0239, 0.0428, 0.0656, infinity) and their
  values into single cells ('-2-100' for -2, -1, 0, 0). The image's values are listed
  below; each block asserts that their concatenation equals the library cell. The
  omega = infinity rows (the high-frequency limit) have the wavelength "inf" in hartree.
All three images checked, every value agrees.
"""

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Naves_2011_DynamicHyperPolarizabilities"
TABLE_II = "tbl_6e45fd86f96d"
TABLE_III = "tbl_a5631d2a2b1a"
PREFIX = "O3/naves_2011__"
COLUMNS_III = ("ccsd_opt", "ccsd_expt", "ccsd_t_expt")
TABLE_I = "tbl_4dfc86410542"
TABLE_VI = "tbl_cb8015ae109a"
PROCESSES = {"Static": "static", "THG": "thg", "dc-SHG": "eshg", "IDRI": "kerr"}
PROCESSES["dc-K"] = "dc_kerr"
TERMS = ("alpha2_00", "mubeta_00", "mu2alpha_10", "mu2alpha_01")
OMEGAS = ("0.0239", "0.0428", "0.0656", "inf")
# Table VI as the image shows it: process -> the four terms' values at OMEGAS.
TABLE_VI_IMAGE = {
    "THG": (
        ("-2", "-1", "0", "0"),
        ("-33", "-10", "-4", "0"),
        ("1", "0", "0", "0"),
        ("0", "0", "0", "0"),
    ),
    "dc-SHG": (
        ("-5", "-2", "-1", "0"),
        ("148", "164", "169", "172"),
        ("-13", "-4", "-2", "0"),
        ("0", "0", "0", "0"),
    ),
    "IDRI": (
        ("144", "145", "145", "145"),
        ("-43", "-13", "-5", "0"),
        ("0", "0", "0", "0"),
        ("0", "0", "0", "0"),
    ),
    "dc-K": (
        ("75", "77", "78", "79"),
        ("322", "337", "340", "343"),
        ("78", "95", "100", "103"),
        ("-3", "-3", "-3", "-3"),
    ),
}
GAMMA = {
    "\\gamma_{xxxx}": "xxxx",
    "\u03b3xxxx": "xxxx",
    "\\gamma_{xxyy}": "xxyy",
    "\u03b3xxyy": "xxyy",
    "\\gamma_{xxzz}": "xxzz",
    "\u03b3xxzz": "xxzz",
    "\\gamma_{yyyy}": "yyyy",
    "\u03b3yyyy": "yyyy",
    "\\gamma_{yyzz}": "yyzz",
    "\u03b3yyzz": "yyzz",
    "\\gamma_{zzzz}": "zzzz",
    "\u03b3zzzz": "zzzz",
    "\\bar{\\gamma}": "mean",
    "\u03b3\u00af": "mean",
}
PREFIX_G = "O3/naves_2011__gamma_"


def _gamma_table_iii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = json_table(DIRECTORY, TABLE_III)
    process = ""
    for cells in grid[2:]:
        if cells[0] == "\\gamma":
            process = "static"
        elif not process:
            continue  # alpha and beta rows
        if cells[1]:
            process = PROCESSES[cells[1]]
        omega = cells[2]
        float(omega)
        for name, value in zip(COLUMNS_III, cells[3:], strict=True):
            float(value)
            files.setdefault(f"{PREFIX_G}{process}_{name}.csv", []).append(
                Point(omega, omega, value, "", " | ".join(cells))
            )
    return files


def _gamma_static() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = json_table(DIRECTORY, TABLE_I)
    columns = [f"{b}_{m}" for b in ("aug", "d-aug", "t-aug") for m in ("HF", "CCSD")]
    for cells in grid[2:]:
        name = GAMMA.get(cells[0])
        if name is None:
            continue
        for column, value in zip(columns, cells[1:], strict=True):
            float(value)
            files.setdefault(f"{PREFIX_G}basis_{name}.csv", []).append(
                Point(column, "0", value, "", " | ".join(cells))
            )
    grid = json_table(DIRECTORY, TABLE_II)
    columns = [
        "CCSD_analytical",
        "CCSD_ff_unrelaxed",
        "CCSD_ff_relaxed",
        "CCSD(T)_ff_relaxed",
    ]
    for cells in grid[2:]:
        name = GAMMA.get(cells[0])
        if name is None:
            continue
        for column, value in zip(columns, cells[1:], strict=True):
            float(value)
            files.setdefault(f"{PREFIX_G}method_{name}.csv", []).append(
                Point(column, "0", value, "", " | ".join(cells))
            )
    return files


def _gamma_vibrational() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = json_table(DIRECTORY, TABLE_VI)
    for cells in grid[1:]:
        label = cells[0]
        printed = " | ".join(cells)
        process = PROCESSES[label]
        if label == "Static":
            for term, value in zip(TERMS, cells[2:], strict=True):
                float(value)
                files.setdefault(f"{PREFIX_G}pv_{process}_{term}.csv", []).append(
                    Point("0", "0", value, "", printed)
                )
            continue
        if cells[1] != "0.02390.04280.0656\u221e":
            msg = f"unexpected Table VI omega cell {cells[1]!r}"
            raise ValueError(msg)
        for term, cell, values in zip(
            TERMS, cells[2:], TABLE_VI_IMAGE[label], strict=True
        ):
            if "".join(values) != cell:
                msg = f"Table VI {label}/{term}: {values!r} do not make {cell!r}"
                raise ValueError(msg)
            for omega, value in zip(OMEGAS, values, strict=True):
                files.setdefault(f"{PREFIX_G}pv_{process}_{term}.csv", []).append(
                    Point(omega, omega, value, "", printed)
                )
    return files


def extract() -> dict[str, list[Point]]:
    """Return Table III's mean polarizability per method and Table II's relaxed CCSD."""
    grid = json_table(DIRECTORY, TABLE_III)
    if grid[2][0] != "\\alpha":
        msg = f"unexpected Table III layout: {grid[2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {f"{PREFIX}{c}.csv": [] for c in COLUMNS_III}
    for cells in grid[2:6]:  # the four alpha rows (the gamma rows: _gamma_table_iii)
        omega = cells[2]
        float(omega)
        printed = " | ".join(cells)
        for name, value in zip(COLUMNS_III, cells[3:], strict=True):
            float(value)
            files[f"{PREFIX}{name}.csv"].append(Point(omega, omega, value, "", printed))
    (mean,) = [
        r
        for r in json_table(DIRECTORY, TABLE_II)
        if r[0] == "\N{GREEK SMALL LETTER ALPHA}\N{MACRON}"
    ]
    float(mean[3])
    files[f"{PREFIX}ccsd_relaxed_static.csv"] = [
        Point("0", "0", mean[3], "", " | ".join(mean))
    ]
    return files | _gamma_table_iii() | _gamma_static() | _gamma_vibrational()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
