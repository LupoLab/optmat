# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Zhao, Nie, Wang & Hu (2026): CH4, CO, O2 molar polarizabilities at 1566 nm.

The paper prints only the *static* molar polarizabilities A_eps it derives; the 1566 nm
values it measured are not printed. Three tables (CSV; images checked, cells agree):

- Table 3 (``tbl_5032bd285fb3``): A_eps (cm^3/mol) of "This work" and the literature
  (RI, DP, DOSD) for CH4, CO, O2 (columns in that order). "—" cells are
  empty. The parenthesised last-digit standard uncertainty goes to the sigma column.
- Table 1 (``tbl_632ce5af228b``): the dipole sums S(-4), S(-6), S(-8) (a.u.), A_mu x
  10^6 (cm^3/mol), B_T and B_eps for CH4, O2, CO (columns in that order).
- Table 2 (``tbl_72bc6d72934f``): the relative standard uncertainty budget (ppm).

THE 1566 NM VALUE IS RECONSTRUCTED here, by inverting the authors' own conversion
(Eqs. 2 and 3): the refractometer measures, to first order in density, the refractive
molar polarization A_R(omega) = A_eps'(omega) + A_mu (Eq. 2: p/RT = 2(n - 1)/(3(A_eps +
A_mu)) + ...), and Eq. (3) gives A_eps'(omega) = A_eps + (4 pi/3) N_A a0^3 (1 +
m_e/m_w)^3 [S(-4) omega^2 + S(-6) omega^4 + S(-8) omega^6] with omega = 45.56397823/
lambda (lambda the vacuum wavelength in nm, as printed). So

    A_R(1566 nm) = A_eps(Table 3)
                   + (4 pi/3) N_A a0^3 [S(-4) w^2 + S(-6) w^4 + S(-8) w^6]
                   + A_mu(Table 1),

with N_A and a0 from optmat.constants (CODATA via SciPy). The reduced-mass factor (1 +
m_e/m_w)^3, whose "nuclear mass" m_w the paper does not define for a molecule, is left
out: it changes the dispersion term by less than 1.4e-4 of itself (m_w >= 12 u), under
0.5 ppm of A_R. The laser wavelength is given only as "1566 nm"; +-0.5 nm moves A_R by
+-2 ppm (CH4). The uncertainty of A_R is Table 2's budget without the S(i) and A_mu
lines (which the reconstruction removes again): sqrt(T^2 + P^2 + B^2 + statistical^2 +
purity^2 + outgassing^2 + thermal expansion^2 + frequency^2) ppm, the "<0.1" entries
taken as 0.1.
"""

import math
import re

from _library import Point, number, table, write

from optmat.constants import BOHR_RADIUS, N_A

DIRECTORY = "Zhao_2026_ElectricPolarizabilitiesCH"
TABLE_1 = "tbl_632ce5af228b"
TABLE_2 = "tbl_72bc6d72934f"
TABLE_3 = "tbl_5032bd285fb3"
WAVELENGTH_NM = "1566"
# Eq. (3): omega (hartree) = 45.56397823 / lambda (nm), as printed (p. 4).
OMEGA_NUMERATOR = 45.56397823
# (4 pi / 3) N_A a0^3 in cm^3/mol: one atomic unit of polarizability as a molar value.
MOLAR_AU = 4 * math.pi / 3 * N_A * (BOHR_RADIUS * 100) ** 3
_VALUE = re.compile(r"^(?P<number>-?\d+(?:\.\d+)?)\((?P<digits>\d+)\)$")
# Table 2 rows that stay in the uncertainty of the measured A_R(1566 nm).
MEASUREMENT_ROWS = (
    "T",
    "P",
    "B",
    "Statistical",
    "Gas purity",
    "Out gassing",
    "Cavity thermal expansion",
    "Optical frequency fluctuation",
)


def _split(cell: str) -> tuple[str, str]:
    """Return a value and its last-digit uncertainty in the value's unit."""
    match = _VALUE.match(number_like(cell))
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    value, digits = match["number"], match["digits"]
    decimals = len(value.split(".")[1]) if "." in value else 0
    sigma = int(digits) * 10.0**-decimals
    return value, f"{sigma:.{decimals}f}" if decimals else f"{sigma:.0f}"


def number_like(cell: str) -> str:
    """Strip a cell and turn a Unicode minus into ``-``."""
    return cell.strip().replace("\N{MINUS SIGN}", "-")


def _columns(rows: list[list[str]]) -> dict[str, int]:
    """Map a header row's gas labels to column indices."""
    return {label.strip().replace("_", ""): i for i, label in enumerate(rows[0])}


def _table_1() -> dict[str, dict[str, list[str]]]:
    """Return Table 1 by gas: {row label: [cell, ...]} with the printed row."""
    rows = table(DIRECTORY, TABLE_1)
    columns = _columns(rows)
    assert list(columns)[1:] == ["CH4", "O2", "CO"], rows[0]
    labels = [r[0] for r in rows[1:]]
    assert labels[:4] == [
        "S(-4) (a.u.)",
        "S(-6) (a.u.)",
        "S(-8) (a.u.)",
        r"A_\mu \times 10^6 (cm^3·mol^{-1})",
    ], labels
    return {
        gas: {r[0]: [r[columns[gas]], " | ".join(r)] for r in rows[1:]}
        for gas in ("CH4", "O2", "CO")
    }


def _budget(gas: str) -> float:
    """Return the relative standard uncertainty of A_R(1566 nm), from Table 2."""
    rows = table(DIRECTORY, TABLE_2)
    column = _columns(rows)[gas]
    parts = {r[0].strip(): r[column].strip() for r in rows[1:]}
    assert set(MEASUREMENT_ROWS) | {"S(i)", "Aµ", "Total"} == set(parts), parts
    ppm = [float(parts[name].removeprefix("<")) for name in MEASUREMENT_ROWS]
    return math.sqrt(sum(p * p for p in ppm)) * 1e-6


def _reconstructed(gas: str, static: str, table_1: dict[str, list[str]]) -> Point:
    """Return A_R(1566 nm) = A_eps + electronic dispersion + A_mu."""
    omega = OMEGA_NUMERATOR / float(WAVELENGTH_NM)
    s4, _ = _split(table_1["S(-4) (a.u.)"][0])
    s6, _ = _split(table_1["S(-6) (a.u.)"][0])
    s8, _ = _split(table_1["S(-8) (a.u.)"][0])
    mu, _ = _split(table_1[r"A_\mu \times 10^6 (cm^3·mol^{-1})"][0])
    dispersion = (
        float(s4) * omega**2 + float(s6) * omega**4 + float(s8) * omega**6
    ) * MOLAR_AU
    value = float(static) + dispersion + float(mu) * 1e-6
    sigma = _budget(gas) * value
    printed = " || ".join(
        [
            f"Table 3 This work {gas}: {static}",
            *(table_1[k][1] for k in list(table_1)[:4]),
        ]
    )
    return Point(WAVELENGTH_NM, WAVELENGTH_NM, f"{value:.6f}", f"{sigma:.6f}", printed)


def extract() -> dict[str, list[Point]]:
    """Return the static values, the literature values, the sums and A_R(1566 nm)."""
    rows = table(DIRECTORY, TABLE_3)
    columns = _columns(rows)
    assert list(columns)[1:] == ["CH4", "CO", "O2"], rows[0]
    by_row = {r[0].strip(): r for r in rows[1:]}
    assert list(by_row) == ["This work", "RI", "DP", "DOSD"], list(by_row)
    table_1 = _table_1()
    files: dict[str, list[Point]] = {}
    names = {
        ("CH4", "This work"): "static",
        ("CH4", "RI"): "hohm_2013_ri",
        ("CH4", "DP"): "hohm_2013_dp",
        ("O2", "This work"): "static",
        ("O2", "DP"): "may_2008_297K",
        ("O2", "DOSD"): "kumar_2011_dosd",
        ("CO", "This work"): "static",
        ("CO", "RI"): "hohm_2013_ri",
        ("CO", "DP"): "tsankova_2018_dp",
        ("CO", "DOSD"): "kumar_meath_1994_dosd",
    }
    for (gas, label), name in names.items():
        cells = by_row[label]
        value, sigma = _split(cells[columns[gas]])
        files[f"{gas}/zhao_2026__{name}_{gas}.csv"] = [
            Point("static", "inf", value, sigma, " | ".join(cells))
        ]
    for gas, name in (
        ("CH4", "olney_1997_moments"),
        ("O2", "kumar_2011_moments"),
        ("CO", "kumar_meath_1994_moments"),
    ):
        points = []
        for k in ("-4", "-6", "-8"):
            cell, printed = table_1[gas][f"S({k}) (a.u.)"]
            value, sigma = _split(cell)
            points.append(Point(k, k, value, sigma, printed))
        files[f"{gas}/zhao_2026__{name}_{gas}.csv"] = points
    for gas in ("CH4", "O2", "CO"):
        static = files[f"{gas}/zhao_2026__static_{gas}.csv"][0].value
        files[f"{gas}/zhao_2026__A_R_1566_{gas}.csv"] = [
            _reconstructed(gas, static, table_1[gas])
        ]
    assert all(number(p.value) for points in files.values() for p in points)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
