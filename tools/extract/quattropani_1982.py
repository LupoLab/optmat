# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Quattropani, Bassani and Carillo (1982): hydrogen 1s-ns amplitudes.

Phys. Rev. A 25, 3079-3089 (1982); journal p. 3079 is PDF page 1, so the journal page is
the PDF page plus 3078.

* **Table I** (PDF pp. 4-5, journal pp. 3082-3083; library tables ``tbl_bd158100b034``,
  parts (a) 1s-3s and (b) 1s-6s, and ``tbl_84d4824d7da0``, part (c) 1s-20s; JSON cells
  only): the dimensionless amplitude ``D_1^n`` of Eq. (3), against ``nu_1 = omega_1 /
  R``, in the length gauge ``J0`` and the velocity gauge ``J``, each as a discrete part,
  a continuum part and their total. Taken: the totals, one file per final level and
  gauge, ``row`` and ``x`` the printed ``nu_1``. The parts are gauge-dependent and not
  taken. Every cell was checked against the table images (2026-10-09).

Two totals are misprinted (the record's errata; the files keep the printed strings):
(a) ``nu_1 = 0.6875`` prints ``0.6963`` in both gauges for -0.6963, and (b) ``nu_1 =
0.9400`` prints ``-115.1149`` in the velocity gauge for -151.1149.

**Checks** (:func:`checks`; the script fails if one is missed):

* each total against the sum of its printed parts, within 1.5e-4: every row in the
  length gauge; in the velocity gauge every row but three whose parts are themselves
  inconsistent (listed in ``J_PARTS_INCONSISTENT``; asserted to stay so);
* the two gauges against each other, within 0.011;
* every total against the exact amplitude of non-relativistic hydrogen, within 0.006
  (:func:`amplitude`): Eq. (3a) is ``D = [T(E_1 + omega_1) + T(E_1 + omega_2)] / 4``
  with ``T(E) = <ns| r G_1(E) r |1s>`` in hartree (radial functions, ``G_1`` the l = 1
  Coulomb Green's function), solved by Dalgarno and Lewis's method in the basis
  ``r^(j+2) exp(-a r)``, ``a = sqrt(-2E)``, in Decimal arithmetic. The worst difference
  is 0.0052 (n = 6 and 20 near the poles); each erratum's correction is within 2e-4 of
  the exact value and its printed string is not.
"""

import functools
from decimal import Decimal, getcontext
from math import comb, factorial

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Quattropani_1982_TwoPhotonTransitions"
SOURCE = "quattropani_1982"
# Table, the parts it holds, and the final level of each part.
TABLES = (
    ("tbl_bd158100b034", {"(a)": 3, "(b)": 6}),
    ("tbl_84d4824d7da0", {"(c)": 20}),
)
# A row: nu_1, then J0 discrete, continuum, total, then J discrete, continuum, total.
GAUGES = (("length", 1), ("velocity", 4))
ERRATA = {
    (3, "length", "0.6875"): ("0.6963", "-0.6963"),
    (3, "velocity", "0.6875"): ("0.6963", "-0.6963"),
    (6, "velocity", "0.9400"): ("-115.1149", "-151.1149"),
}
"""Misprinted totals, (n, gauge, nu_1) -> (printed, corrected): the record's errata."""
J_PARTS_INCONSISTENT = {(3, "0.8860"), (6, "0.9000"), (6, "0.9625")}
"""Velocity-gauge rows whose printed parts do not sum to the printed total (the total is
right: it agrees with the length gauge and the exact amplitude)."""
PARTS_TOLERANCE = 1.5e-4
GAUGE_TOLERANCE = 0.011
EXACT_TOLERANCE = 0.006
BASIS = {3: (36, 60), 6: (36, 60), 20: (70, 80)}
"""Basis size and Decimal digits for each final level (converged to below 1e-5)."""


def _rows() -> list[tuple[int, list[str]]]:
    """Return (n, cells) for every data row of Table I, the cells without spaces."""
    rows = []
    for table_id, parts in TABLES:
        n = 0
        for cells in json_table(DIRECTORY, table_id):
            if cells[0] in parts:
                n = parts[cells[0]]
            elif n and cells[0][:1].isdigit():
                rows.append((n, [c.replace(" ", "") for c in cells]))
    return rows


def _poly_ns(n: int) -> dict[int, Decimal]:
    """Return c_p of r R_n0(r) = sum c_p r^p exp(-r/n), in atomic units."""
    two_n = Decimal(2) / n
    norm = (two_n**3 * factorial(n - 1) / (2 * n * factorial(n))).sqrt()
    return {
        m + 1: norm * (-1) ** m * comb(n, n - 1 - m) / factorial(m) * two_n**m
        for m in range(n)
    }


def _solve(matrix: list[list[Decimal]], rhs: list[Decimal]) -> list[Decimal]:
    """Return the solution of a linear system by Gaussian elimination with pivoting."""
    size = len(rhs)
    rows = [[*row, rhs[i]] for i, row in enumerate(matrix)]
    for col in range(size):
        pivot = max(range(col, size), key=lambda r: abs(rows[r][col]))
        rows[col], rows[pivot] = rows[pivot], rows[col]
        for r in range(col + 1, size):
            factor = rows[r][col] / rows[col][col]
            for k in range(col, size + 1):
                rows[r][k] -= factor * rows[col][k]
    solution = [Decimal(0)] * size
    for r in range(size - 1, -1, -1):
        rest = sum(rows[r][k] * solution[k] for k in range(r + 1, size))
        solution[r] = (rows[r][size] - rest) / rows[r][r]
    return solution


def _t(n: int, energy: Decimal, size: int) -> Decimal:
    """Return <ns| r G_1(E) r |1s> (radial, hartree) for E below the threshold."""
    a = (-2 * energy).sqrt()

    def moment(k: int) -> Decimal:
        return factorial(k) / (2 * a) ** (k + 1)

    half = Decimal(1) / 2
    matrix = [
        [
            half
            * (
                (i + 2) * (j + 2) * moment(i + j + 2)
                - a * (i + j + 4) * moment(i + j + 3)
                + a * a * moment(i + j + 4)
            )
            + moment(i + j + 2)
            - moment(i + j + 3)
            - energy * moment(i + j + 4)
            for j in range(size)
        ]
        for i in range(size)
    ]
    rhs = [2 * factorial(i + 4) / (a + 1) ** (i + 5) for i in range(size)]
    coefficients = _solve(matrix, rhs)
    decay = a + Decimal(1) / n
    return sum(
        (
            c * cp * factorial(p + j + 3) / decay ** (p + j + 4)
            for j, c in enumerate(coefficients)
            for p, cp in _poly_ns(n).items()
        ),
        Decimal(0),
    )


@functools.cache
def amplitude(n: int, nu_1: str, size: int = 36, digits: int = 60) -> float:
    """Return the exact D_1^n of Eq. (3a) at nu_1 (rydberg), both photons below 1 Ry."""
    getcontext().prec = digits
    nu1 = Decimal(nu_1)
    nu2 = 1 - Decimal(1) / (n * n) - nu1
    e1 = Decimal(-1) / 2
    return float((_t(n, e1 + nu1 / 2, size) + _t(n, e1 + nu2 / 2, size)) / 4)


def _total(n: int, gauge: str, nu_1: str, printed: str) -> float:
    """Return a printed total, corrected where the record has an erratum."""
    erratum = ERRATA.get((n, gauge, nu_1))
    if erratum is None:
        return float(printed)
    if printed != erratum[0]:
        msg = f"{n}s {gauge} {nu_1}: erratum expects {erratum[0]!r}, not {printed!r}"
        raise ValueError(msg)
    return float(erratum[1])


def checks() -> dict[str, object]:
    """Return the internal and exact checks of Table I; raise on a miss."""
    out: dict[str, object] = {}
    failures: list[str] = []
    worst = 0.0
    for n, cells in _rows():
        nu = cells[0]
        totals = {}
        for gauge, first in GAUGES:
            d, c, t = (float(number(x)) for x in cells[first : first + 3])
            total = _total(n, gauge, nu, number(cells[first + 2]))
            totals[gauge] = total
            inconsistent = gauge == "velocity" and (n, nu) in J_PARTS_INCONSISTENT
            parts_ok = abs(d + c - total) <= PARTS_TOLERANCE
            if parts_ok == inconsistent:
                failures.append(f"{n}s {gauge} {nu}: parts {d + c:.4f}, total {t}")
        if abs(totals["length"] - totals["velocity"]) > GAUGE_TOLERANCE:
            failures.append(f"{n}s {nu}: gauges differ")
        size, digits = BASIS[n]
        exact = amplitude(n, nu, size, digits)
        for gauge, total in totals.items():
            worst = max(worst, abs(total - exact))
            if abs(total - exact) > EXACT_TOLERANCE:
                failures.append(f"{n}s {gauge} {nu}: {total} against exact {exact:.5f}")
            erratum = ERRATA.get((n, gauge, nu))
            if erratum and abs(float(erratum[0]) - exact) <= EXACT_TOLERANCE:
                failures.append(f"{n}s {gauge} {nu}: the printed value is not wrong")
        out[f"{n}s nu_1 = {nu}"] = (*totals.values(), round(exact, 5))
    out["worst |total - exact|"] = round(worst, 5)
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return the totals, one file per final level and gauge, after the checks."""
    checks()
    files: dict[str, list[Point]] = {}
    for n, cells in _rows():
        for gauge, first in GAUGES:
            name = f"H/{SOURCE}__D_{n}s_{gauge}.csv"
            printed = " | ".join(cells)
            value = number(cells[first + 2])
            point = Point(cells[0], cells[0], value, "", printed)
            files.setdefault(name, []).append(point)
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
