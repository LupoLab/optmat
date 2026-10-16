# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Drachman, Bhatia and Shabazz (1990): hydrogen 1S-2S two-photon amplitude.

Phys. Rev. A 42, 6333-6335 (1990); journal p. 6333 is PDF page 1, so the journal page is
the PDF page plus 6332.

* **Table IV** (PDF p. 3, journal p. 6335; library table ``tbl_a9b929131d75``, JSON
  cells only): the amplitude in its length form ``D1`` (Eq. 5), its momentum form ``D2``
  (Eq. 6) and the momentum form rewritten by the commutator relation, ``D3`` (Eq. 8), at
  nine photon energies ``omega`` (rydberg; the other photon is ``3/4 - omega``), each
  for ``N0 = 5``, 10 and 15 pseudostates (``a = 1``), and the column "Ref. 1" (Bassani,
  Forney and Quattropani 1977). Taken: the ``N0 = 15`` line of ``D1``, ``D2`` and
  ``D3``, and the "Ref. 1" column, one file each, ``row`` and ``x`` the printed
  ``omega``. Every cell was checked against the table image (2026-10-09). The cells
  print digit groups ("-11.780 483"); ``value`` drops the spaces, ``printed`` keeps the
  line as printed, its ``omega`` filled in.

Not taken: the ``N0 = 5`` and 10 lines (convergence), Tables II and III (cumulative
sums at ``omega = 0.375``; used by :func:`checks`) and Table I (pseudostate energies).

**Checks** (:func:`checks`; the script fails if one is missed), against the exact
amplitude of quattropani_1982's script (:func:`quattropani_1982.amplitude` with ``n =
2``; Eq. 5 is its Eq. 3a for the 2s level):

* the ``N0 = 15`` ``D1`` and ``D2`` within 1.5 units of their last printed digit;
* ``D3`` within 5e-6 relative (it converges more slowly, as the text says);
* "Ref. 1" within 1.5 units of its last printed digit;
* Tables II and III: the ``N = 20`` sums of ``D1`` and ``D2`` (both scales ``a``) within
  1.5e-5 of the exact value at ``omega = 0.375``, and their first term, -17.87845, the
  2p term alone, ``(1/2)(8/3 + 8/3) R_1s^2p R_2s^2p`` with the closed-form radial
  integrals ``2^7 sqrt(6) / 3^5`` and ``-3 sqrt(3)``.
"""

import math

from _json_table import json_table
from _library import Point, number, write
from quattropani_1982 import amplitude

DIRECTORY = "Drachman_1990_TwoPhotonTransitions"
SOURCE = "drachman_1990"
TABLE_IV = "tbl_a9b929131d75"
TABLES_II_III = ("tbl_91240e2dabf7", "tbl_d35f7eb262b9")
HEADER = ["\\omega", "D_{1}", "D_{2}", "D_{3}", "Ref. 1"]
# Dataset, column of Table IV, line of each omega's group (0: N0 = 5, 2: N0 = 15).
TAKEN = (
    ("D1_N15", 1, 2),
    ("D2_N15", 2, 2),
    ("D3_N15", 3, 2),
    ("ref1_bassani_1977", 4, 0),
)
D3_RELATIVE = 5e-6
SUM_TOLERANCE = 1.5e-5


def _groups() -> list[tuple[str, list[list[str]]]]:
    """Return (omega, its three lines) for every omega of Table IV."""
    grid = json_table(DIRECTORY, TABLE_IV)
    if grid[0] != HEADER:
        msg = f"Table IV heading {grid[0]!r}"
        raise ValueError(msg)
    body = grid[1:]
    groups = [(body[i][0], body[i : i + 3]) for i in range(0, len(body), 3)]
    if len(groups) != 9 or any(lines[1][0] or lines[2][0] for _, lines in groups):
        msg = "Table IV is not nine groups of three lines"
        raise ValueError(msg)
    return groups


def _plain(cell: str) -> str:
    """Return a cell's number without its digit-group spaces."""
    return number(cell.replace(" ", ""))


def _unit(text: str) -> float:
    """Return one unit of the last printed digit of a number."""
    decimals = text.split(".")[1] if "." in text else ""
    return 10.0 ** -len(decimals)


def checks() -> dict[str, object]:
    """Return the checks against the exact amplitude; raise on a miss."""
    out: dict[str, object] = {}
    failures: list[str] = []
    for omega, lines in _groups():
        exact = amplitude(2, omega)
        for name, column, line in TAKEN:
            text = _plain(lines[line][column])
            value = float(text)
            if name.startswith("D3"):
                ok = abs(value - exact) <= D3_RELATIVE * abs(exact)
            else:
                ok = abs(value - exact) <= 1.5 * _unit(text)
            if not ok:
                failures.append(f"{name} {omega}: {text} against exact {exact:.7f}")
        taken = (_plain(lines[2][1]), _plain(lines[2][2]))
        out[f"omega = {omega}"] = (*taken, round(exact, 7))
    exact = amplitude(2, "0.3750")
    first = 0.5 * (8 / 3 + 8 / 3) * (2**7 * math.sqrt(6) / 3**5) * (-3 * math.sqrt(3))
    for table_id in TABLES_II_III:
        grid = json_table(DIRECTORY, table_id)
        last = next(row for row in grid if row[0] == "20")
        for cell in last[1:3]:
            if abs(float(_plain(cell)) - exact) > SUM_TOLERANCE:
                failures.append(f"{table_id} N = 20: {cell} against exact {exact:.7f}")
        term = next(row for row in grid if row[0] == "1")[1]
        if abs(float(_plain(term)) - first) > 1e-5:
            failures.append(f"{table_id} N = 1: {term} against the 2p term {first:.6f}")
        out[f"{table_id} N = 20, first term"] = (last[1], last[2], term)
    out["2p term"] = round(first, 6)
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return the N0 = 15 amplitudes and the Ref. 1 column, after the checks."""
    checks()
    files: dict[str, list[Point]] = {}
    for omega, lines in _groups():
        for name, column, line in TAKEN:
            printed = " | ".join([omega, *lines[line][1:]])
            point = Point(omega, omega, _plain(lines[line][column]), "", printed)
            files.setdefault(f"H/{SOURCE}__{name}.csv", []).append(point)
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
