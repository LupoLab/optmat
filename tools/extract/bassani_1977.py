# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bassani, Forney and Quattropani (1977): hydrogen 1s-2s amplitudes.

Phys. Rev. Lett. 39, 1070-1073 (1977); journal p. 1070 is PDF page 1, so the journal
page is the PDF page plus 1069.

* **Table I** (PDF p. 2, journal p. 1071; library table ``tbl_76d993024c8f``) and
  **Table II** (PDF p. 3, journal p. 1072; ``tbl_bb20405731ed``): the dimensionless
  amplitude ``D`` of Eq. (3) for 1s-2s in the length gauge ``J0`` (E.r, Table I) and the
  velocity gauge ``J`` (A.p, Table II), at nine photon frequencies ``nu_1 = omega_1 /
  (2 pi R)`` (columns, 0.3750-0.7475; the other photon ``nu_2 = 3/4 - nu_1``), each as
  the discrete sum up to ``N = 2, ..., 13`` and ``N -> infinity`` (Eqs. 4a, 4b), the
  continuum part (Eqs. 5a, 5b) and the total ``D``. Taken: the totals ``D[J0]`` and
  ``D[J]``, one file per gauge, ``row`` and ``x`` the printed ``nu_1``, ``printed`` the
  column's last three cells (``D^d`` at infinity, ``D^c``, ``D``). The partial sums and
  the continuum parts are gauge-dependent, not observables, and not taken; they enter
  :func:`checks`. Every cell was checked against the table images (2026-10-09).

**Checks** (:func:`checks`; the script fails if one is missed). Every one of the 270
printed cells (15 in each of nine columns of both tables) must equal the value computed
here, rounded or truncated to the printed digits (22 cells are truncated, e.g. the
velocity-gauge total at 0.6875, -49.687778 printed ``-49.6877``; no cell is off):

* the discrete partial sums (``N = 2, ..., 13``) of Eqs. (4a) and (4b), with the
  closed-form radial integrals of Bethe and Salpeter (Eq. 63.4 there; ``R_10^n1`` and
  ``R_20^n1``, positive except ``R_20^21 = -sqrt(27)``, as the paper says). Eq. (4b) is
  evaluated with ``[nu_1 (3/4 - nu_1)]^-1``, as in Eq. (5b) and the text's factor
  ``(E_2s - E_n)(E_n - E_1s) / hbar^2 omega_1 omega_2``; with the ``(1/4 - nu_1)`` that
  Eq. (4b) prints, no ``N = 3`` sum reproduces (asserted);
* the ``N -> infinity`` sums: the same series to ``n = 4000`` plus its ``n^-3`` tail;
* the totals ``D[J0]`` and ``D[J]``: the exact amplitude of quattropani_1982's script
  (:func:`quattropani_1982.amplitude` with ``n = 2``);
* the continuum parts: the exact amplitude less the exact discrete sum.

Also each total against the sum of its printed parts (within 2 units of its last
digit), and the column "Ref. 1" of drachman_1990's Table IV (``tbl_a9b929131d75``),
which reproduces this paper: equal as strings to ``D[J0]`` at all nine ``nu_1``
(``D[J]`` differs at 0.6875 only, ``-49.6877`` against ``-49.6878``).
"""

import math
from decimal import ROUND_DOWN, ROUND_HALF_EVEN, Decimal

from _json_table import json_table
from _library import Point, number, table, write
from quattropani_1982 import amplitude

DIRECTORY = "Bassani_1977_ChoiceGaugeTwo"
SOURCE = "bassani_1977"
TABLES = (("length", "tbl_76d993024c8f"), ("velocity", "tbl_bb20405731ed"))
DRACHMAN = ("Drachman_1990_TwoPhotonTransitions", "tbl_a9b929131d75")
PARTIAL_N = tuple(range(2, 14))
SERIES_N = 4000
PARTS_UNITS = 2


def _exponent(text: str) -> int:
    """Return the decimal exponent of the last printed digit of a number."""
    return -len(text.split(".")[1]) if "." in text else 0


def _printing(value: float, text: str) -> str | None:
    """Return how text prints value to its digits: "rounded", "truncated" or None."""
    printed = Decimal(text)
    digits = Decimal(1).scaleb(_exponent(text))
    exact = Decimal(repr(value))
    if exact.quantize(digits, ROUND_HALF_EVEN) == printed:
        return "rounded"
    if exact.quantize(digits, ROUND_DOWN) == printed:
        return "truncated"
    return None


def _radial_1s(n: int) -> float:
    """Return R_10^n1 = <np| r |1s> (radial, atomic units), positive."""
    log = 8 * math.log(2) + 7 * math.log(n) + (2 * n - 5) * math.log(n - 1)
    return math.exp(0.5 * (log - (2 * n + 5) * math.log(n + 1)))


def _radial_2s(n: int) -> float:
    """Return R_20^n1 = <np| r |2s> (radial, atomic units); -sqrt(27) for n = 2."""
    if n == 2:  # the 2p term, the one negative integral
        return -math.sqrt(27)
    log = 17 * math.log(2) + 7 * math.log(n) + math.log(n * n - 1)
    log += (2 * n - 6) * math.log(n - 2) - (2 * n + 6) * math.log(n + 2)
    return math.exp(0.5 * log)


def _term(n: int, nu_1: float, gauge: str, second: float = 0.75) -> float:
    """Return the np term of Eq. (4a) or (4b), ``second`` the 3/4 of Eq. (4b)."""
    level = 1 / n**2
    term = 0.5 * (1 / (1 - level - nu_1) + 1 / (0.25 - level + nu_1))
    term *= _radial_1s(n) * _radial_2s(n)
    if gauge == "velocity":
        term *= -(0.25 - level) * (1 - level) / (nu_1 * (second - nu_1))
    return term


def _discrete(big_n: int, nu_1: float, gauge: str, second: float = 0.75) -> float:
    """Return the discrete sum of Eq. (4a) or (4b) over n = 2, ..., big_n."""
    return sum(_term(n, nu_1, gauge, second) for n in range(2, big_n + 1))


def _discrete_infinite(nu_1: float, gauge: str) -> float:
    """Return the discrete sum to infinity: to SERIES_N, plus its n^-3 tail."""
    tail = _term(SERIES_N, nu_1, gauge) * SERIES_N**3 / (2 * SERIES_N**2)
    return _discrete(SERIES_N, nu_1, gauge) + tail


def _columns(table_id: str) -> list[tuple[str, dict[str, str]]]:
    """Return (nu_1, {row label: cell}) for each column of a table."""
    grid = table(DIRECTORY, table_id)
    labels = ["N=2", *(str(n) for n in PARTIAL_N[1:]), "inf", "c", "total"]
    if len(grid) != len(labels) + 1:
        msg = f"{table_id}: {len(grid)} rows"
        raise ValueError(msg)
    out = []
    for j, nu in enumerate(grid[0][1:]):
        column = [row[j + 1].replace(" ", "") for row in grid[1:]]
        out.append((nu, dict(zip(labels, column, strict=True))))
    return out


def _drachman_ref1() -> dict[str, str]:
    """Return drachman_1990's Table IV column "Ref. 1" (first line of each group)."""
    grid = json_table(*DRACHMAN)
    return {grid[i][0]: grid[i][4] for i in range(1, len(grid), 3)}


def _unit(text: str) -> float:
    """Return one unit of the last printed digit of a number."""
    return 10.0 ** _exponent(text)


def _cells(gauge: str, nu: str, cells: dict[str, str]) -> list[tuple[str, str, float]]:
    """Return (row label, printed cell, computed value) for every cell of a column."""
    x = float(nu)
    labels = ["N=2", *map(str, PARTIAL_N[1:])]
    out = [
        (label, number(cells[label]), _discrete(big_n, x, gauge))
        for big_n, label in zip(PARTIAL_N, labels, strict=True)
    ]
    exact = amplitude(2, nu)
    infinite = _discrete_infinite(x, gauge)
    out += [
        ("inf", number(cells["inf"]), infinite),
        ("c", number(cells["c"]), exact - infinite),
        ("total", number(cells["total"]), exact),
    ]
    return out


def checks() -> dict[str, object]:
    """Return the checks of both tables; raise on a miss."""
    out: dict[str, object] = {}
    failures: list[str] = []
    totals: dict[str, dict[str, str]] = {}
    truncated: list[str] = []
    for gauge, table_id in TABLES:
        for nu, cells in _columns(table_id):
            for label, text, value in _cells(gauge, nu, cells):
                how = _printing(value, text)
                if how is None:
                    failures.append(f"{gauge} {nu} {label}: {text} against {value:.7f}")
                elif how == "truncated":
                    truncated.append(f"{gauge} {nu} {label}")
            total = number(cells["total"])
            summed = float(number(cells["inf"])) + float(number(cells["c"]))
            if abs(summed - float(total)) > PARTS_UNITS * _unit(total):
                failures.append(f"{gauge} {nu}: parts sum to {summed}, total {total}")
            totals.setdefault(nu, {})[gauge] = total
            if gauge == "velocity":
                exact = amplitude(2, nu)
                out[f"nu_1 = {nu}"] = (
                    totals[nu]["length"],
                    total,
                    round(exact, 6),
                    f"continuum share {float(number(cells['c'])) / exact:.3f} (A.p)",
                )
    # Eq. (4b) as printed, [nu_1 (1/4 - nu_1)]^-1, does not give the N = 3 sums.
    for nu, cells in _columns(TABLES[1][1]):
        text = number(cells["3"])
        if _printing(_discrete(3, float(nu), "velocity", 0.25), text) is not None:
            failures.append(f"Eq. (4b) as printed reproduces N = 3 at {nu}")
    ref1 = _drachman_ref1()
    reproduced = [ref1.get(nu, "").replace(" ", "") for nu in totals]
    if reproduced != [t["length"] for t in totals.values()]:
        failures.append(f"drachman_1990 Ref. 1 {ref1} is not D[J0]")
    differ = [nu for nu, t in totals.items() if t["length"] != t["velocity"]]
    if differ != ["0.6875"]:
        failures.append(f"the gauges' printed totals differ at {differ}")
    out["cells printed truncated, not rounded"] = (len(truncated), truncated)
    out["drachman_1990 Ref. 1"] = "equals D[J0] at all nine nu_1"
    out["printed D[J0] != D[J]"] = differ
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return the totals D[J0] and D[J], one file per gauge, after the checks."""
    checks()
    files: dict[str, list[Point]] = {}
    for gauge, table_id in TABLES:
        for nu, cells in _columns(table_id):
            printed = " | ".join([nu, cells["inf"], cells["c"], cells["total"]])
            point = Point(nu, nu, number(cells["total"]), "", printed)
            files.setdefault(f"H/{SOURCE}__D_2s_{gauge}.csv", []).append(point)
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
