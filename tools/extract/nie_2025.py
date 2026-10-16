# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Nie et al. (2025): N2, Ne and H2 molar polarizabilities at 1566 nm.

The paper prints only the *static* molar polarizabilities A_eps it derives from its
1566 nm refractometer measurements (the method of ``zhao_2026``, by the same group).
Two tables (images checked, cells agree):

- Table I (``tbl_c8f73c6e439f``, CSV): A_mu x 10^6 (cm^3/mol), B_T, B_eps and the dipole
  sums S(-4), S(-6), S(-8) (a.u.) for Ar, Ne, N2, H2 (columns in that order).
- Table II (``tbl_2c3bf20e3331``, cells file only): the relative standard uncertainty
  budget (ppm) for N2, Ne, H2 (columns in that order), then, under the row
  "A_eps (cm^3 mol^-1)", the static A_eps of "This work", "Expt.^a" and "Theor.^b".
  The parenthesised last-digit standard uncertainty goes to the sigma column.

THE 1566 NM VALUE IS RECONSTRUCTED here by inverting the authors' Eq. (3), as for
``zhao_2026``: A_R(1566 nm) = A_eps + (4 pi/3) N_A a0^3 (1 + m_e/m_N)^3 [S(-4) w^2 +
S(-6) w^4 + S(-8) w^6] + A_mu, with w = 45.56397823/1566 (hartree; lambda the vacuum
wavelength in nm, as printed) and N_A, a0, m_e, m_p from optmat.constants and SciPy
(CODATA). The factor (1 + m_e/m_N)^3 (m_N "the nuclear mass") is applied to H2 with
m_N the proton mass (it raises the dispersion term by 1.6e-3 of itself, 5 ppm of
A_R); for N2 and Ne it is left out (under 0.1 ppm of A_R for any nuclear mass above 14
u). The uncertainty of A_R is Table II's budget without the S(i) and A_mu lines: the
root sum of squares of T, p, K, gas purity, outgas, frequency, L_cav, statistics,
B_eps and B_T, the "<0.1" entries taken as 0.1.
"""

import math
import re

from _json_table import json_table
from _library import Point, table, write
from scipy import constants

from optmat.constants import BOHR_RADIUS, M_E, N_A

DIRECTORY = "Nie_2025_PolarizabilityMolecularHydrogen"
TABLE_I = "tbl_c8f73c6e439f"
TABLE_II = "tbl_2c3bf20e3331"
WAVELENGTH_NM = "1566"
OMEGA_NUMERATOR = 45.56397823  # Eq. (3): omega (hartree) = 45.563 978 23 / lambda (nm)
MOLAR_AU = 4 * math.pi / 3 * N_A * (BOHR_RADIUS * 100) ** 3  # cm^3/mol per a.u.
# The reduced-mass factor of Eq. (3), applied to H2 only (see the module docstring).
MASS_FACTOR = {"H2": (1 + M_E / constants.m_p) ** 3, "N2": 1.0, "Ne": 1.0}
GASES = ("N2", "Ne", "H2")
MEASUREMENT_ROWS = (
    "T",
    "p",
    "K",
    "Gas purity",
    "Outgas",
    "Frequency",
    "L_{cav}",
    "Statistics",
    "B_{\\varepsilon}",
    "B_{T}",
)
# Table II row label -> dataset name, per gas.
LITERATURE = {
    ("This work", "N2"): "static",
    ("This work", "Ne"): "static",
    ("This work", "H2"): "static",
    ("Expt.^{a}", "N2"): "egan_2019",
    ("Expt.^{a}", "Ne"): "gaiser_2018",
    ("Expt.^{a}", "H2"): "hohm_kerl_1990",
    ("Theor.^{b}", "N2"): "olney_1997",
    ("Theor.^{b}", "Ne"): "lesiuk_2020",
    ("Theor.^{b}", "H2"): "sun_2021",
}
_VALUE = re.compile(r"^(?P<number>-?\d+(?:\.\d+)?)\((?P<digits>\d+)\)$")


def _split(cell: str) -> tuple[str, str]:
    """Return a value and its last-digit uncertainty in the value's unit."""
    match = _VALUE.match(cell.strip().replace("\N{MINUS SIGN}", "-"))
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    value, digits = match["number"], match["digits"]
    decimals = len(value.split(".")[1]) if "." in value else 0
    sigma = int(digits) * 10.0**-decimals
    return value, f"{sigma:.{decimals}f}" if decimals else f"{sigma:.0f}"


def _table_i() -> dict[str, dict[str, list[str]]]:
    """Return Table I by gas: {row label: [cell, printed row]}."""
    rows = table(DIRECTORY, TABLE_I)
    header = [c.replace("_{", "").replace("}", "").strip() for c in rows[0]]
    assert header[1:] == ["Ar", "Ne", "N2", "H2"], header
    labels = [r[0].strip() for r in rows[1:]]
    assert labels[0].startswith("A_{\\mu}") and labels[3:] == [
        "S(-4) (a.u.)",
        "S(-6) (a.u.)",
        "S(-8) (a.u.)",
    ], labels
    return {
        gas: {r[0].strip(): [r[header.index(gas)], " | ".join(r)] for r in rows[1:]}
        for gas in GASES
    }


def _table_ii() -> tuple[dict[str, dict[str, str]], dict[str, list[str]]]:
    """Return the budget by gas {row: ppm} and the A_eps rows {label: cells}."""
    rows = json_table(DIRECTORY, TABLE_II)
    header = [c.replace("_{", "").replace("}", "").strip() for c in rows[0]]
    assert header[1:] == ["N2", "Ne", "H2"], header
    split = [r[0] for r in rows].index("A_{\\varepsilon} (cm^{3} mol^{-1})")
    budget = {
        gas: {r[0].strip(): r[header.index(gas)].strip() for r in rows[1:split]}
        for gas in GASES
    }
    values = {r[0].strip(): r for r in rows[split + 1 :]}
    assert list(values) == ["This work", "Expt.^{a}", "Theor.^{b}"], list(values)
    return budget, values


def _relative_uncertainty(parts: dict[str, str]) -> float:
    """Return the root sum of squares of the measurement rows of the budget (ppm)."""
    assert set(MEASUREMENT_ROWS) | {"A_{\\mu}", "S(i)", "Total"} == set(parts), parts
    ppm = [float(parts[name].removeprefix("<")) for name in MEASUREMENT_ROWS]
    return math.sqrt(sum(p * p for p in ppm)) * 1e-6


def extract() -> dict[str, list[Point]]:
    """Return the static and literature values, the sums and A_R(1566 nm)."""
    table_i = _table_i()
    budget, values = _table_ii()
    files: dict[str, list[Point]] = {}
    for (label, gas), name in LITERATURE.items():
        cells = values[label]
        value, sigma = _split(cells[1 + GASES.index(gas)])
        files[f"{gas}/nie_2025__{name}_{gas}.csv"] = [
            Point("static", "inf", value, sigma, " | ".join(cells))
        ]
    omega = OMEGA_NUMERATOR / float(WAVELENGTH_NM)
    for gas in GASES:
        moments = []
        for k in ("-4", "-6", "-8"):
            cell, printed = table_i[gas][f"S({k}) (a.u.)"]
            value, sigma = _split(cell)
            moments.append(float(value))
            files.setdefault(f"{gas}/nie_2025__olney_moments_{gas}.csv", [])
            files[f"{gas}/nie_2025__olney_moments_{gas}.csv"].append(
                Point(k, k, value, sigma, printed)
            )
        mu_label = next(k for k in table_i[gas] if k.startswith("A_{\\mu}"))
        mu, _ = _split(table_i[gas][mu_label][0])
        static = files[f"{gas}/nie_2025__static_{gas}.csv"][0].value
        dispersion = (
            MASS_FACTOR[gas]
            * MOLAR_AU
            * (moments[0] * omega**2 + moments[1] * omega**4 + moments[2] * omega**6)
        )
        a_r = float(static) + dispersion + float(mu) * 1e-6
        sigma = _relative_uncertainty(budget[gas]) * a_r
        printed = " || ".join(
            [
                f"Table II This work {gas}: {static}",
                table_i[gas][mu_label][1],
                *(table_i[gas][f"S({k}) (a.u.)"][1] for k in ("-4", "-6", "-8")),
            ]
        )
        files[f"{gas}/nie_2025__A_R_1566_{gas}.csv"] = [
            Point(WAVELENGTH_NM, WAVELENGTH_NM, f"{a_r:.6f}", f"{sigma:.6f}", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
