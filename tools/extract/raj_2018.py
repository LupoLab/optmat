# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Raj, Hamaguchi & Witek (2018) and thermally average their H2 and D2 values.

Library tables (checked against the table images):

- Supplement Tables T1 (H2, ``tbl_b0cae9bc8235``) and T5 (D2, ``tbl_3f331311ba2d``):
  the mean polarizability matrix elements at 52 laser wavelengths (nm, the program's
  vacuum wavelengths) and the static limit (printed "Static"), four decimals, a.u.
  Only the diagonal columns <v=0,J=0|alpha|v=0,J=0> and <1,0|alpha|1,0> are extracted
  (the Raman transition elements are not refractive-index quantities).
- Main-paper Table II (``tbl_aa58a97e916a``, JSON only): the authors' own thermal
  rotational average of H2 at 293 K ("Our", averaged up to v = 0, J = 7) at 15
  wavelengths. The Bishop & Cheung (Ref. 75) and Victor & Dalgarno (Ref. 82) columns
  are reproductions of values curated in ``bishop_1980``; they are compared there, not
  extracted.
- Supplement Tables T26 (H2 S0 and O1, ``tbl_8f44df72c173``; Q1 and S1 in the
  uncaptioned continuation ``tbl_db455a894d4d``) and T29 (D2, ``tbl_7c62aa367012``),
  JSON only: the authors' own computed transition energies (cm^-1, "Our" columns),
  from the same rovibrational wave functions the program integrates over. The pure
  rotational S0(J) = E(0, J + 2) - E(0, J) for J = 0-8 give the v = 0 rotational
  levels within each nuclear-spin species, J = 0-10 (the even levels from J = 0, the
  odd levels from J = 1). They are inputs of the average, not data files; each S0 is
  checked against S1(J) - Q1(J + 2) and Q1(J) - O1(J + 2) from the same tables.

Program values (``--program``): the per-state diagonal elements <0,J|alpha|0,J>,
J = 0-10, of H2 and D2 at T1's wavelengths and the static limit, computed with the
authors' own Python module ``rovibME`` (their GitHub package H2-PolarizabilityMatrix
Elements, Ref. 39, the copy at ``$RAJ_2018_PROGRAM``) through
``_raj_2018_program.py``, which runs it unmodified in its pinned environment. The
module's SHA-256 is checked against ``ROVIBME_SHA256``. Written to
``raj_2018__program_J<J>.csv`` with the module's printed line kept in ``printed``.
These files are not re-made by ``extract()`` (the program takes minutes and runs
outside the project environment); ``extract()`` reads them, checks every J = 0 value
against T1 and T5 (to their four printed decimals), and averages them.

Thermal averages (``extract()``): normal H2 and normal D2, i.e. the nuclear-spin
species in their high-temperature proportions g_even : g_odd (H2 1 : 3, para : ortho;
D2 6 : 3, ortho : para), each species Boltzmann-distributed over its own rotational
levels:

    alpha(T) = sum_s f_s sum_{J in s} (2J + 1) exp(-E_J / kT) alpha_J
                         / sum_{J in s} (2J + 1) exp(-E_J / kT),
    f_s = g_s / (g_even + g_odd),

with E_J from S0 as above and kT in cm^-1 from CODATA (``optmat.constants``).
Only v = 0 is included: the v = 1 population (Q1(0) = 4160.9 cm^-1 for H2, 2993.5 for
D2, T26 and T29) is below 5e-7 at 296.15 K. Spin weights: H2 1 (even J) and 3 (odd J)
as in Bishop & Cheung (1980), Eq. (23) text, and HITRAN's H2 lower-state weights
g'' = (2J+1) g_i; D2 even : odd = 2 : 1, stated for D2 (nuclear spin I = 1) in
Unverified_09b870d2f279, p. 15, and the rule 6 : 3 for two spin-1 nuclei in Borysow
et al. (2019) for 14N2 (source notes). Temperatures 273.15, 293.15 and 296.15 K.

``--compare`` runs the program at Bishop & Cheung's, Rychlewski's and Table II's
wavelengths and prints the 293 K averages against theirs (source notes).
"""

import json
import math
import os
import subprocess
import sys
from pathlib import Path

from _json_table import json_table
from _library import Point, number, path, table, write

from optmat.constants import HBAR, K_B, C

DIRECTORY = "Raj_2018_PolarizabilityTensorInvariants"
SUPPLEMENT = f"{DIRECTORY}/supplement_01"
MEAN_TABLES = {"H2": "tbl_b0cae9bc8235", "D2": "tbl_3f331311ba2d"}  # T1, T5
TABLE_II = "tbl_aa58a97e916a"
# Transition energies: (S0/O1 table, Q1/S1 table, first data row of each).
ENERGY_TABLES = {
    "H2": ("tbl_8f44df72c173", 3, "tbl_db455a894d4d", 2),
    "D2": ("tbl_7c62aa367012", 3, "tbl_7c62aa367012", 15),
}
# Raj et al.'s published program (their supplementary repository), found at
# $RAJ_2018_PROGRAM; needed only with --program, to regenerate the program files.
PROGRAM = Path(os.environ.get("RAJ_2018_PROGRAM", "RAJ_2018_PROGRAM-is-not-set"))
ROVIBME_SHA256 = "75823f2ee89d6d649d980e84deaf82480a7424e8ae2b228ead70a78c68197570"
RUNNER = Path(__file__).resolve().with_name("_raj_2018_program.py")
J_MAX = 10  # the highest J whose energy the S0 columns (J = 0-8) give
# Nuclear-spin weights of the even-J and odd-J levels (source notes).
SPIN_WEIGHTS = {"H2": (1, 3), "D2": (6, 3)}
TEMPERATURES = {"273K": "273.15", "293K": "293.15", "296K": "296.15"}
# k / (h c) in cm^-1 per kelvin.
K_OVER_HC = K_B / (2 * math.pi * HBAR * C * 100)


def _wavelength(cell: str) -> str:
    """Return a T1/T5 wavelength cell as printed, "inf" for "Static"."""
    cell = cell.strip()
    return "inf" if cell == "Static" else number(cell)


def _mean_columns(gas: str) -> dict[str, list[Point]]:
    """Return T1's or T5's <0,0|alpha|0,0> and <1,0|alpha|1,0> columns."""
    rows = table(SUPPLEMENT, MEAN_TABLES[gas])[1:]
    files = {}
    for column, stem in ((1, "alpha_v0_J0"), (2, "alpha_v1_J0")):
        files[f"{gas}/raj_2018__{stem}.csv"] = [
            Point(
                _wavelength(c[0]),
                _wavelength(c[0]),
                number(c[column]),
                "",
                " | ".join(c),
            )
            for c in rows
        ]
    return files


def _table_ii() -> list[Point]:
    """Return Table II's thermal average of H2 at 293 K ("Our", rows 21-35)."""
    grid = json_table(DIRECTORY, TABLE_II)
    start = next(i for i, c in enumerate(grid) if "bar{\\alpha} \\rangle av" in c[1])
    points = []
    for cells in grid[start:]:
        if not cells[2]:
            break  # the header of the anisotropy block
        points.append(
            Point(
                number(cells[2]),
                number(cells[2]),
                number(cells[3]),
                "",
                " | ".join(cells),
            )
        )
    return points


def _column(table_id: str, first: int, column: int) -> dict[int, float]:
    """Return a transition-energy column by initial J (the "Our" values, cm^-1)."""
    grid = json_table(SUPPLEMENT, table_id)
    values = {}
    for cells in grid[first : first + 9]:
        if cells[column]:
            values[int(cells[0])] = float(number(cells[column]))
    return values


def rotational_energies(gas: str) -> dict[int, float]:
    """Return E(v=0, J) in cm^-1, J = 0-10, within each nuclear-spin species.

    Built from the S0 column; each S0(J) is checked against S1(J) - Q1(J + 2) and
    Q1(J) - O1(J + 2) from the same tables. The tables agree with themselves to
    7e-4 cm^-1 (D2: S1 - Q1 at J = 2 differs by 3e-4 with four printed decimals, O1
    has three), so the check (1e-3 cm^-1) catches misread digits, not rounding.
    """
    s0_table, s0_first, q1_table, q1_first = ENERGY_TABLES[gas]
    s0 = _column(s0_table, s0_first, 1)
    o1 = _column(s0_table, s0_first, 5)
    q1 = _column(q1_table, q1_first, 1)
    s1 = _column(q1_table, q1_first, 5)
    for j, value in s0.items():
        if j in s1 and j + 2 in q1:
            assert abs(s1[j] - q1[j + 2] - value) <= 1e-3, (gas, j, "S1 - Q1")
        if j in q1 and j + 2 in o1:
            assert abs(q1[j] - o1[j + 2] - value) <= 1e-3, (gas, j, "Q1 - O1")
    energies = {0: 0.0, 1: 0.0}
    for j in range(J_MAX - 1):
        energies[j + 2] = energies[j] + s0[j]
    return energies


def thermal_average(
    gas: str, alphas: dict[int, float], energies: dict[int, float], temperature: float
) -> float:
    """Return the normal-gas thermal average of per-J polarizabilities."""
    kt = K_OVER_HC * temperature
    weights = SPIN_WEIGHTS[gas]
    total = 0.0
    for parity, weight in enumerate(weights):
        levels = [j for j in alphas if j % 2 == parity]
        boltzmann = {j: (2 * j + 1) * math.exp(-energies[j] / kt) for j in levels}
        species = sum(boltzmann[j] * alphas[j] for j in levels) / sum(
            boltzmann.values()
        )
        total += weight / sum(weights) * species
    return total


def _program_file(gas: str, j: int) -> Path:
    return path(f"{gas}/raj_2018__program_J{j}.csv")


def _program_values(gas: str) -> dict[int, dict[str, str]]:
    """Return the committed program values: J -> wavelength row -> printed value."""
    import csv  # noqa: PLC0415

    values = {}
    for j in range(J_MAX + 1):
        with _program_file(gas, j).open(newline="", encoding="utf-8") as handle:
            values[j] = {r["row"]: r["value"] for r in csv.DictReader(handle)}
    return values


def _averages(gas: str, printed_j0: list[Point]) -> dict[str, list[Point]]:
    """Check the program's J = 0 values against T1/T5 and average over J."""
    program = _program_values(gas)
    rows = [p.row for p in printed_j0]
    for point in printed_j0:
        computed = float(program[0][point.row])
        # T1/T5 print four decimals: the program's value must round to the printed one.
        assert abs(computed - float(point.value)) <= 0.5e-4 + 1e-9, (gas, point.row)
    energies = rotational_energies(gas)
    files = {}
    for label, temperature in TEMPERATURES.items():
        points = []
        for row in rows:
            alphas = {j: float(program[j][row]) for j in program}
            average = thermal_average(gas, alphas, energies, float(temperature))
            printed = f"normal {gas}, {temperature} K, v = 0, J = 0-{J_MAX}"
            points.append(Point(row, row, f"{average:.7f}", "", printed))
        files[f"{gas}/raj_2018__normal_{label}.csv"] = points
    return files


def extract() -> dict[str, list[Point]]:
    """Return the tabulated columns and the thermal averages of the program values."""
    files: dict[str, list[Point]] = {}
    for gas in MEAN_TABLES:
        columns = _mean_columns(gas)
        files |= columns
        files |= _averages(gas, columns[f"{gas}/raj_2018__alpha_v0_J0.csv"])
    files["H2/raj_2018__thermal_293K_table_ii.csv"] = _table_ii()
    return files


def run_program(jobs: list[list[object]]) -> list[tuple[list[object], str, str]]:
    """Run the authors' module on ``jobs`` and return (job, value, printed line)."""
    request = json.dumps({"module": str(PROGRAM / "python-module"), "jobs": jobs})
    result = subprocess.run(
        ["uv", "run", "--quiet", "--script", str(RUNNER)],
        input=request,
        capture_output=True,
        text=True,
        check=True,
    )
    header, *lines = result.stdout.splitlines()
    versions = json.loads(header)
    print(f"rovibME run with {versions}")
    if versions["rovibME_sha256"] != ROVIBME_SHA256:
        msg = f"rovibME.py differs from the version recorded: {versions}"
        raise SystemExit(msg)
    return [(job, value, line) for job, value, _, line in map(json.loads, lines)]


def _program_files() -> dict[str, list[Point]]:
    """Compute <0,J|alpha|0,J> for every J and T1's wavelengths with the program."""
    rows = [_wavelength(c[0]) for c in table(SUPPLEMENT, MEAN_TABLES["H2"])[1:]]
    jobs: list[list[object]] = [
        [gas, 0, j, 0, j, "static" if row == "inf" else float(row), "nm"]
        for gas in MEAN_TABLES
        for j in range(J_MAX + 1)
        for row in rows
    ]
    results = iter(run_program(jobs))
    files = {}
    for gas in MEAN_TABLES:
        for j in range(J_MAX + 1):
            points = []
            for row in rows:
                job, value, line = next(results)
                assert job[0] == gas and job[2] == j, job
                points.append(Point(row, row, value, "", line))
            files[f"{gas}/raj_2018__program_J{j}.csv"] = points
    return files


# The 293 K averages compared: (data file, the program's unit for its wavelengths).
COMPARISONS = {
    "H2": (
        ("H2/bishop_1980__alpha_293K.csv", "A"),
        ("H2/rychlewski_1983__alpha_293K.csv", "A"),
        ("H2/raj_2018__thermal_293K_table_ii.csv", "nm"),
    ),
    "D2": (("D2/rychlewski_1983__alpha_293K.csv", "A"),),
}
# The program's range, 0.25-0.0345 hartree, in its own units (rovibME's range check).
PROGRAM_RANGE = {"A": (1822.5341, 13206.7688), "nm": (182.25341, 1320.67688)}


def compare(temperature: float = 293.0) -> None:
    """Print the normal-gas average against Bishop & Cheung, Rychlewski and Table II.

    Runs the program at the other sources' wavelengths within its range (and the
    static limit) for J = 0-10, averages as ``extract()`` does, and prints each value
    with the difference (ours minus theirs).
    """
    import csv  # noqa: PLC0415

    for gas, comparisons in COMPARISONS.items():
        energies = rotational_energies(gas)
        for relative, unit in comparisons:
            with path(relative).open(encoding="utf-8") as handle:
                theirs = {r["row"]: float(r["value"]) for r in csv.DictReader(handle)}
            low, high = PROGRAM_RANGE[unit]
            rows = [r for r in theirs if r == "inf" or low <= float(r) <= high]
            jobs: list[list[object]] = [
                [gas, 0, j, 0, j, "static" if r == "inf" else float(r), unit]
                for r in rows
                for j in range(J_MAX + 1)
            ]
            results = iter(run_program(jobs))
            print(f"{gas} normal {temperature} K against {relative} ({unit}):")
            differences = []
            for row in rows:
                alphas = {}
                for _ in range(J_MAX + 1):
                    job, value, _line = next(results)
                    alphas[job[2]] = float(value)
                ours = thermal_average(gas, alphas, energies, temperature)
                difference = ours - theirs[row]
                differences.append(difference / theirs[row])
                print(
                    f"  {row:>8}  {theirs[row]:.4f}  {ours:.5f}  {difference:+.5f}"
                    f"  {1e2 * difference / theirs[row]:+.3f} %"
                )
            skipped = [r for r in theirs if r not in rows]
            print(
                f"  mean {1e2 * sum(differences) / len(differences):+.3f} %, "
                f"max |{1e2 * max(map(abs, differences)):.3f}| %; "
                f"outside the program's range: {', '.join(skipped) or 'none'}"
            )


def main() -> None:
    """Write the data files; with ``--program``, first run the authors' module.

    With ``--compare``, also print the comparison with the other 293 K averages.
    """
    if "--program" in sys.argv:
        write(_program_files())
    write(extract())
    if "--compare" in sys.argv:
        compare()


if __name__ == "__main__":
    main()
