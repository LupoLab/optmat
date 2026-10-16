# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Yoshino and Freeman (1985): the Rydberg levels of Xe from its VUV absorption.

J. Opt. Soc. Am. B 2, 1268-1274 (1985); journal p. 1268 is PDF page 1. The five
odd-parity J = 1 Rydberg series seen in absorption from the 5p6 1S0 ground level, read
on photographic plates in the second order of a 6.65-m normal-incidence spectrograph
(922-1296 A, vacuum wavelengths), one table each:

* **Table 1** (``tbl_fc5223dd7155``, p. 1270): 5p5 nd[3/2]1, n = 5-66;
* **Table 2** (``tbl_8f0a90e1d919``, p. 1270): 5p5 nd[1/2]1, n = 5-34;
* **Table 3** (``tbl_0c56a9f6ed72``, p. 1271): 5p5 ns[3/2]1, n = 7-47;
* **Table 4** (``tbl_60d55d81bd56``, p. 1271): 5p5 nd'[3/2]1, n = 5-65, autoionising
  above the first limit (n >= 6);
* **Table 5** (``tbl_ca7c8f0fa582``, p. 1272): 5p5 ns'[1/2]1, n = 6-39, autoionising
  for n >= 8.

Each row prints n, the wavelength (A), the wave number (cm^-1) and the effective
quantum number n*. The wave number of a line from the ground level is the energy of
its upper level: one ``level_energy`` file per table, the row and the ``x`` label the
level (n, and the NIST-style configuration), the value the printed wave number with
its digit-group spaces removed. Rows printed in parentheses are the authors'
estimates from the quantum defect, not measurements (Table 3 n = 22, under Lyman beta;
Table 4 n = 6, 7 and Table 5 n = 8, 9, too diffuse): kept in the files and excluded by
the record.

:func:`checks` tests every row against the paper itself, twice: the wave number
against 1e8 over the vacuum wavelength, and n* against the Rydberg formula with the
paper's series limits (97 833.75 and 108 370.8 cm^-1, p. 1271) and the Rydberg
constant of xenon; both within the printed rounding. Six wave numbers fail both and
are misprints (``MISPRINTS``), corrected by the record's errata; the corrected values
pass both. Four of the six are confirmed by NIST (``data/external/nist_asd``): 11s'
and 11d' are NIST's levels with the bibliographic key L6449 (this paper), 10s and
10d[1/2] agree with NIST's levels from other sources within 0.01 cm^-1; the other two
(n = 23 of ns[3/2], n = 31 of nd[1/2]) lie beyond NIST's list. Two wavelengths (Table
1 n = 14, '1030.453' for 1030.435; Table 5 n = 12, '3937.650' for 937.650) and one n*
(Table 4 n = 23, '20.74' for 20.72) are misprinted too; they are only in ``printed``.
It also compares the levels with NIST's (n <= 20): the levels NIST cites by L6449 are
this paper's numbers to every printed digit (28 rows); the others differ by -0.26 to
+0.23 cm^-1.

The tables are plain grids; every cell of the rows that matter was checked against
the table images (2026-10-09).
"""

import csv
import math
from pathlib import Path
from typing import Final

from _library import DATA, Point, table, write
from scipy.constants import Rydberg, atomic_mass

from optmat.constants import M_E
from optmat.gas.eos import equation_of_state

DIRECTORY: Final = "Yoshino_1985_AbsorptionSpectrumXenon"
FIRST_LIMIT: Final = 97833.75
"""The 2P3/2 series limit, '97 833.75 +- 0.05 cm-1' (p. 1271), from Table 1."""
SECOND_LIMIT: Final = 108370.8
"""The 2P1/2 series limit, '108 370.8 +- 0.2 cm-1' (p. 1271), from Table 4."""

SERIES: Final = {
    "nd32": ("tbl_fc5223dd7155", "5p5.(2P*<3/2>).{n}d", "2[3/2]*", FIRST_LIMIT),
    "nd12": ("tbl_8f0a90e1d919", "5p5.(2P*<3/2>).{n}d", "2[1/2]*", FIRST_LIMIT),
    "ns32": ("tbl_0c56a9f6ed72", "5p5.(2P*<3/2>).{n}s", "2[3/2]*", FIRST_LIMIT),
    "ndp32": ("tbl_60d55d81bd56", "5p5.(2P*<1/2>).{n}d", "2[3/2]*", SECOND_LIMIT),
    "nsp12": ("tbl_ca7c8f0fa582", "5p5.(2P*<1/2>).{n}s", "2[1/2]*", SECOND_LIMIT),
}
"""Each table: its library id, the NIST-style configuration, term and series limit."""

MISPRINTS: Final = {
    ("nd12", "10"): ("95912.98", "95912.87"),
    ("nd12", "31"): ("97798.4", "97698.4"),
    ("ns32", "10"): ("94487.09", "94787.09"),
    ("ns32", "23"): ("97350.0", "97530.0"),
    ("ndp32", "11"): ("107209.", "106929."),
    ("nsp12", "11"): ("1061208", "106120.8"),
}
"""The wave numbers the record corrects by errata: the file's string, the correction."""
WAVELENGTH_MISPRINTS: Final = {
    ("nd32", "14"): ("1030.453", "1030.435"),
    ("nsp12", "12"): ("3937.650", "937.650"),
}
"""Misprinted wavelengths, kept only in ``printed``: the cell, the correction."""
NSTAR_MISPRINTS: Final = {("ndp32", "23"): ("20.74", "20.72")}
"""A misprinted n*, kept only in ``printed``: the cell, the correction."""

NIST_KEY: Final = "L6449"
"""The bibliographic key of the NIST levels that are this paper's (n <= 20)."""


def _decimals(text: str) -> int:
    """Return the number of printed decimals of a number string."""
    return len(text.partition(".")[2])


def _wave_number(cell: str) -> str:
    """Return a printed wave number without its digit-group spaces or parentheses."""
    return cell.strip().strip("()").replace(" ", "")


def rows(series: str) -> list[list[str]]:
    """Return a table's data rows, as printed, without the header."""
    grid = table(DIRECTORY, SERIES[series][0])
    if grid[0][0] != "n" or "Wave Number" not in grid[0][2]:
        msg = f"{series}: the header moved: {grid[0]}"
        raise ValueError(msg)
    return grid[1:]


def rydberg_xenon() -> float:
    """Return the Rydberg constant of xenon, cm^-1, with the atom's mean mass."""
    mass = equation_of_state("Xe").molar_mass / 1e-3 * atomic_mass
    return Rydberg / 100 / (1 + M_E / mass)


def _nist() -> dict[tuple[str, str], tuple[float, str]]:
    """Return NIST's odd J = 1 levels of Xe I: (configuration, term) -> (cm^-1, key)."""
    path = DATA / "external" / "nist_asd" / "Xe_I_levels.tsv"
    levels = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["J"].strip('"') == "1":
                key = (row["Configuration"].strip('"'), row["Term"].strip('"'))
                levels[key] = (
                    float(row["Level (cm-1)"].strip('"')),
                    row["Reference"].strip('"'),
                )
    return levels


def _corrected(series: str, cells: list[str], failures: list[str]) -> tuple[str, ...]:
    """Return a row's n, wavelength, printed and corrected wave number, and n*.

    The known misprints of the wavelength and of n* are replaced by their corrections
    (the corrected wave number is empty where the printed one is right); a known
    misprint that the table no longer holds is reported in ``failures``.
    """
    n, wavelength, wave, nstar = (c.strip() for c in cells)
    printed = _wave_number(wave)
    known = MISPRINTS.get((series, n))
    if known and known[0] != printed:
        failures.append(f"{series} n = {n}: the file holds {printed}")
    fixed = WAVELENGTH_MISPRINTS.get((series, n))
    if fixed and fixed[0] != wavelength:
        failures.append(f"{series} n = {n}: the wavelength is {wavelength}")
    wrong = NSTAR_MISPRINTS.get((series, n))
    if wrong and wrong[0] != nstar:
        failures.append(f"{series} n = {n}: n* is {nstar}")
    return (
        n,
        fixed[1] if fixed else wavelength,
        printed,
        known[1] if known else "",
        wrong[1] if wrong else nstar,
    )


def _misses(text: str, lam: str, nstar: str, limit: float) -> tuple[float, float]:
    """Return a wave number's misses of 1e8 / lambda and of n*, over their rounding.

    Each is the miss divided by what the printed rounding allows: a row passes when
    both are at most 1.
    """
    rydberg = rydberg_xenon()
    value = float(text)
    # The rounding of the printed wavelength and wave number: half a unit of the last
    # digit of each, the first through d(1e8/lambda) = 1e8 dlambda / lambda^2.
    half_lam = 0.5 * 10.0 ** -_decimals(lam)
    half_wave = 0.5 * 10.0 ** -_decimals(text)
    allowed = 1e8 * half_lam / float(lam) ** 2 + half_wave + 1e-9
    # n* = sqrt(R / (limit - T)); its rounding, and the wave number's rounding
    # through dn*/dT = n*^3 / (2 R). The four-decimal n* of the lowest members
    # reproduce only to 1-2 units of their last digit (2.2 for nd[1/2] n = 6; the
    # authors' R and rounding are not stated): 2.5 units allowed there, half a unit
    # elsewhere.
    computed = math.sqrt(rydberg / (limit - value)) if value < limit else 0
    places = _decimals(nstar)
    nstar_allowed = (
        (2.5 if places >= 4 else 0.5) * 10.0**-places
        + computed**3 / (2 * rydberg) * half_wave
        + 1e-9
    )
    return (
        abs(value - 1e8 / float(lam)) / allowed,
        abs(computed - float(nstar)) / nstar_allowed,
    )


def _series_checks(series: str, failures: list[str]) -> dict[str, object]:
    """Return one table's checks against its wavelengths, n* and NIST."""
    _, configuration, term, limit = SERIES[series]
    nist = _nist()
    worst_wave, worst_nstar = 0.0, 0.0
    matches, differences = 0, []
    for cells in rows(series):
        if cells[1].strip().startswith("("):
            continue
        n, lam, printed, corrected, nstar = _corrected(series, cells, failures)
        wave_miss, nstar_miss = _misses(printed, lam, nstar, limit)
        passes = wave_miss <= 1 and nstar_miss <= 1
        if corrected:
            if passes:
                failures.append(f"{series} n = {n}: misprint {printed} passes")
            if max(_misses(corrected, lam, nstar, limit)) > 1:
                failures.append(f"{series} n = {n}: correction {corrected} fails")
        else:
            if not passes:
                failures.append(
                    f"{series} n = {n}: {printed} against 1e8/{lam} or n* {nstar} "
                    f"({wave_miss:.2f}, {nstar_miss:.2f} of the rounding)"
                )
            worst_wave = max(worst_wave, wave_miss)
            worst_nstar = max(worst_nstar, nstar_miss)
        level = nist.get((configuration.format(n=n), term))
        if level is not None:
            value = float(corrected or printed)
            if level[1] != NIST_KEY:
                differences.append(value - level[0])
            elif value != level[0]:
                failures.append(f"{series} n = {n}: NIST {level[0]} ({NIST_KEY})")
            else:
                matches += 1
    out: dict[str, object] = {
        f"{series}: worst wave-number miss / rounding": round(worst_wave, 3),
        f"{series}: worst n* miss / rounding": round(worst_nstar, 3),
        f"{series}: rows equal to NIST {NIST_KEY}": matches,
    }
    if differences:
        out[f"{series}: other NIST rows, this - NIST (cm-1)"] = (
            len(differences),
            round(min(differences), 3),
            round(max(differences), 3),
        )
    return out


def checks() -> dict[str, object]:
    """Return the checks of every row against its wavelength, n* and NIST.

    Raises if a row fails a check other than the known misprints, a known misprint
    passes or its correction fails, or a level NIST cites to this paper differs.
    """
    out: dict[str, object] = {}
    failures: list[str] = []
    for series in SERIES:
        out.update(_series_checks(series, failures))
    out["R(Xe) used for n* (cm-1)"] = round(rydberg_xenon(), 3)
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return one file per table: the level energies (cm^-1), rows by n."""
    checks()
    files: dict[str, list[Point]] = {}
    for series, (_, configuration, _, _) in SERIES.items():
        points = []
        for cells in rows(series):
            n = cells[0].strip()
            points.append(
                Point(
                    n,
                    configuration.format(n=n),
                    _wave_number(cells[2]),
                    "",
                    " | ".join(c.strip() for c in cells),
                )
            )
        files[f"Xe/yoshino_1985__levels_{series}.csv"] = points
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
