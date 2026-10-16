# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kestin, Korfali, Sengers & Kamgar-Parsi (1981): Ne-Ar mixture densities.

A viscosity paper; only its equation-of-state data are extracted here.

- Table III (p. 421, ``tbl_68ea087a83f3``, cells file only; image checked): for two
  Ne-Ar mixtures at 25 °C, x_Ar = 0.2727 and 0.7722, the pressure P (MPa), the
  measured density rho x 10^3 (mol/cm^3), Z_exp = PV/RT and (Z_exp - Z_calc) x 10^3,
  in two column groups. Written: the density against pressure
  (``Ar/kestin_1981__rho_xAr_<x>.csv``) and Z_exp against pressure
  (``Ar/kestin_1981__Z_xAr_<x>.csv``), the row id the printed pressure; the residual
  column stays in ``printed``.
- Tables IIA and IIB (p. 419, ``tbl_b62f6e4e3136`` and ``tbl_6a5c1d106b1d``): B
  (cm^3/mol) and C ((cm^3/mol)^2) of Ne-Ar and He-Ne mixtures at 25 °C. Written: B of
  the mixtures (the rows 0 < x < 1), one single-point file per composition; the pure
  gases' rows (refits of Michels' isotherms) and C are not.

Pairs are filed under the first species in sorted order: Ne-Ar under ``Ar`` (partner
Ne, x_Ar as printed), He-Ne under ``He`` (partner Ne). Table IIB prints x_Ne; the mole
fraction of the file name (and of the record) is helium's, 1 - x_Ne (0.8021, 0.5965,
0.2172), the printed x_Ne kept in ``printed``.
"""

from decimal import Decimal

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Kestin_1981_DensityDependenceComposition"
TABLE_III = "tbl_68ea087a83f3"
TABLE_IIA = "tbl_b62f6e4e3136"
TABLE_IIB = "tbl_6a5c1d106b1d"
TEMPERATURE_C = "25"


def table_iii() -> dict[str, list[Point]]:
    """Return Table III's densities and Z_exp, per mixture, against pressure."""
    grid = json_table(DIRECTORY, TABLE_III)
    if not grid[0][0].startswith("Experimental pressures and densities") or [
        grid[1][0],
        grid[1][4],
    ] != ["P(MPa)", "P(MPa)"]:
        msg = f"unexpected Table III header {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for start in (0, 4):
        label = grid[2][start]
        x_ar = label.removeprefix("x_{Ar} = ")
        if x_ar == label:
            msg = f"Table III: no x_Ar in {label!r}"
            raise ValueError(msg)
        density, compressibility = [], []
        for cells in grid[3:]:
            pressure, rho, z, residual = cells[start : start + 4]
            if not pressure:
                continue
            pressure = number(pressure)
            printed = f"x_Ar = {x_ar} | {pressure} | {rho} | {z} | {residual}"
            density.append(Point(pressure, pressure, number(rho), "", printed))
            compressibility.append(Point(pressure, pressure, number(z), "", printed))
        files[f"Ar/kestin_1981__rho_xAr_{x_ar}.csv"] = density
        files[f"Ar/kestin_1981__Z_xAr_{x_ar}.csv"] = compressibility
    return files


def table_ii() -> dict[str, list[Point]]:
    """Return B of each mixture of Tables IIA (Ne-Ar) and IIB (He-Ne)."""
    files: dict[str, list[Point]] = {}
    for table_id, header, material, complement in (
        (TABLE_IIA, "x_{\\text{Ar}}", "Ar", False),  # x_Ar, argon's own
        (TABLE_IIB, "x_{\\text{Ne}}", "He", True),  # x_Ne, the partner's
    ):
        grid = table(DIRECTORY, table_id)
        if grid[0][:2] != [header, "B(cm3/mol)"]:
            msg = f"unexpected Table II header {grid[0]!r}"
            raise ValueError(msg)
        for x, b, c in grid[1:]:
            if x in {"0.0000", "1.0000"}:
                continue  # a pure gas: Michels' isotherms refitted
            own = str(Decimal(1) - Decimal(x)) if complement else x
            name = f"{material}/kestin_1981__B_x{material}_{own}.csv"
            printed = f"x_{'Ne' if complement else 'Ar'} = {x} | {b} | {c}"
            files[name] = [Point(TEMPERATURE_C, TEMPERATURE_C, number(b), "", printed)]
    return files


def extract() -> dict[str, list[Point]]:
    """Return Table III's densities and Z, and Table II's mixture B."""
    return table_iii() | table_ii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
