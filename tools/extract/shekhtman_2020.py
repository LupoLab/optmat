# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shekhtman, Mustafa and Parziale (2020): two-photon cross sections of Kr.

The paper's calculated frequency-integrated two-photon excitation cross sections
sigma_o^(2) (cm^4) from the Kr ground state to eight 5p, 5p' and 6p levels, one per
excitation line from 192.749 to 216.667 nm, are printed in Tables 5, 6 and 7 (basis sets
1, 2 and 3; PDF pp. 10-11, journal pp. 10835-10836), each for two ways of evaluating the
dipole matrix elements: quantum-defect radial functions alone ("Quantum-Defect") and
with the NIST oscillator strengths where available ("Quantum-Defect with Oscillator
Strengths"), and in Table 9 (p. 11/10836) in the single-path approximation. Tables 5-7
have merged header cells and no CSV rendering, so they are read from their JSON cells
(``_json_table``); Tables 3 and 9 from their CSV files. Every cell was checked against
the table images.

Each table row is one excitation line; the level it reaches is taken from Table 3
(p. 9/10834), whose lambda_L column gives the laser wavelength of each level's
two-photon excitation, "as measured in vacuum" (its note c). The row id is the level in
Racah notation (``5p[3/2]2``; a prime for the 2P_1/2 core), the abscissa the printed
lambda_L. Values are written as ``<mantissa>e<exponent>`` (the exponent varies from cell
to cell), the printed cells kept in ``printed``.

Each table also prints the profile-dependent rate coefficient sigma^(2) = sigma_o^(2)
g(2 omega_L) (cm^4 s; Eq. 7) and its ratio to the largest. They are not extracted
(the vocabulary's quantity is sigma_o^(2), not a peak value for one line shape); they
serve as checks (:func:`checks`):

* sigma^(2) / sigma_o^(2) is the Gaussian peak of Eq. (9), g = 2 sqrt(ln 2 / pi) /
  sqrt(2 Delta omega_L^2 + Delta omega_T^2), with the laser width Delta omega_L = 2 pi x
  1350 MHz and the Doppler width of Eq. (10) at the two-photon frequency 2 omega_L for
  krypton at T = 298 K (the temperature the paper states for its Eq. 13 example), both
  as angular frequencies: the printed ratio is 1.004 to 1.013 times it in every cell but
  the three listed in :data:`ANOMALIES` (with g normalised in ordinary frequency instead
  it would be 2 pi times smaller). So the line profile is normalised in the two-photon
  angular frequency.
* the normalised column is sigma^(2) over the table's largest sigma^(2), within the
  rounding of three printed figures, but for :data:`ANOMALIES`.
* the lambda_L of Table 3 is 2 h c / (E - E_G) from the table's own level energies, to
  0.002 nm: vacuum wavelengths (standard air would be 0.07 nm shorter).

Five printed cells disagree with their own row (:data:`ANOMALIES`): Table 5's
quantum-defect sigma^(2) at 214.769 and 216.667 nm (printed 1.98e-46 and 2.24e-46, ten
times too small for their sigma_o^(2) and their normalised ratios 0.404 and 0.455: not
extracted, a finding), Table 6's quantum-defect sigma_o^(2) at 193.494 nm (printed
9.85e-35, which repeats the 204.196 nm cell two rows below; its sigma^(2) and ratio
give 2.26e-35: an erratum of the record), Table 7's quantum-defect normalised ratio at
193.494 nm (printed 0.198; its sigma^(2) gives 0.169: a finding), and Table 9's
sigma_o^(2) at 214.769 nm (printed 8.54e-35; its sigma^(2) and ratio give 8.5e-36: an
erratum).
"""

import math
import re

from _json_table import json_table
from _library import Point, table, write

from optmat.constants import ATOMIC_MASS, K_B, C
from optmat.gas.eos import equation_of_state

DIRECTORY = "Shekhtman_2020_TwoPhotonCross"
TABLE_3 = "tbl_5f23a92b0f93"
BASIS_TABLES = {
    "basis1": "tbl_795599a25fb2",
    "basis2": "tbl_c5a835601f3b",
    "basis3": "tbl_9e3a7f52b1d8",
}
TABLE_9 = "tbl_305221c17bc2"
METHODS = {"qd": 1, "qdos": 4}
"""The column of sigma_o^(2) of each method in Tables 5-7; sigma^(2) and the normalised
ratio follow it."""
FIRST_ROW = 3
"""Tables 5-7: two header rows of method and basis, one of column names."""

LASER_LINEWIDTH_HZ = 1350e6
"""Delta omega_L, '1350 MHz in this work' (p. 3), read as 2 pi x 1350 MHz."""
TEMPERATURE_K = 298.0
"""The room temperature of the paper's Eq. (13) example (p. 3)."""
RATIO_RANGE = (0.995, 1.02)
"""sigma^(2) / (sigma_o^(2) g) in every regular cell: the rounding of two three-figure
numbers (0.5 % each) about the observed mean, 1.008 (the authors' temperature or mass
may differ slightly from ours: 1.008 would be T = 292 K)."""
NORMALISED_TOLERANCE = (0.012, 0.0006)
"""|computed - printed| <= 0.012 printed + 0.0006: two three-figure numbers and the
printed ratio's last digit."""
VACUUM_TOLERANCE_NM = 0.002
"""Table 3's lambda_L, printed to 0.001 nm, against 2 h c / (E - E_G)."""
HC_EV_NM = 1239.841984
"""h c in eV nm (CODATA 2018), for Table 3's wavelengths."""
TEXT_SIGMA_214 = ("basis3", "qdos", "214.769", "4.18e-35")
"""Section 5 (p. 9/10834): 'In Table 7, the calculated cross-section for 214.769 nm
excitation is 4.18 x 10^-35 cm4'."""

ANOMALIES = {
    ("basis1", "qd", "214.769", "sigma"): 10.0,
    ("basis1", "qd", "216.667", "sigma"): 10.0,
    ("basis2", "qd", "193.494", "sigma_o"): None,
    ("basis3", "qd", "193.494", "normalised"): None,
    ("single_path", "sp", "214.769", "sigma_o"): 10.0,
}
"""The printed cells that disagree with their row: (table, method, line, cell) -> the
factor by which the row's sigma^(2) / (sigma_o^(2) g) is too small (a misprinted power
of ten), or None where the digits are wrong. Table 6's sigma_o^(2) at 193.494 nm
(9.85e-35, the 204.196 nm cell repeated) is corrected by an erratum to sigma^(2) / g,
:data:`CORRECTED_BASIS2_193`."""
CORRECTED_BASIS2_193 = 2.26e-35
"""The erratum's sigma_o^(2) of Table 6, quantum defect, 193.494 nm: 7.42e-46 cm^4 s
over g(193.494 nm) times the mean ratio of the regular cells."""

_CELL = re.compile(r"^(\d\.\d+) \\times 10\^\{(-\d+)\}$")
_STATE = re.compile(r"^\(2P(\d)o/2\)(\d)([spdf])2\[(\d/\d)\](\d)$")


def _value(cell: str) -> str:
    r"""Return a printed ``m \times 10^{k}`` cell as ``me<k>``."""
    match = _CELL.match(cell.strip())
    if not match:
        msg = f"not a printed power of ten: {cell!r}"
        raise ValueError(msg)
    return f"{match.group(1)}e{match.group(2)}"


def _float(cell: str) -> float:
    """Return a printed cell's number."""
    return float(_value(cell))


def levels() -> dict[str, tuple[str, float]]:
    """Return Table 3's level of each excitation line: lambda_L -> (label, E in eV).

    The label is the Racah notation, ``5p[3/2]2``, primed for the 2P_1/2 core.
    """
    out = {}
    for cells in table(DIRECTORY, TABLE_3)[1:]:
        if cells[4] == "-":
            continue
        state = cells[1].split()[0]
        match = _STATE.match(state)
        if not match:
            msg = f"Table 3: unexpected state {state!r}"
            raise ValueError(msg)
        core, n, ell, k, j = match.groups()
        prime = "'" if core == "1" else ""
        out[cells[4]] = (f"{n}{ell}{prime}[{k}]{j}", float(cells[2]))
    return out


def ground_energy() -> float:
    """Return Table 3's ground-state energy (eV, below the ionisation limit)."""
    return float(table(DIRECTORY, TABLE_3)[1][2])


def rows() -> dict[tuple[str, str], list[list[str]]]:
    """Return the printed rows of each (table, method): line, sigma_o, sigma, ratio."""
    out: dict[tuple[str, str], list[list[str]]] = {}
    for basis, table_id in BASIS_TABLES.items():
        grid = json_table(DIRECTORY, table_id)
        for method, column in METHODS.items():
            out[basis, method] = [
                [cells[0], *cells[column : column + 3]] for cells in grid[FIRST_ROW:]
            ]
    out["single_path", "sp"] = [
        [cells[0], *cells[3:6]] for cells in table(DIRECTORY, TABLE_9)[1:]
    ]
    return out


def gaussian_peak(wavelength_nm: float) -> float:
    """Return Eq. (9)'s peak line shape g(2 omega_L = omega_T), s (per rad/s)."""
    omega = 2 * math.pi * C / (wavelength_nm * 1e-9)
    mass = equation_of_state("Kr").molar_mass / 1e-3 * ATOMIC_MASS
    doppler = (
        2 * omega * math.sqrt(8 * math.log(2) * K_B * TEMPERATURE_K / (mass * C**2))
    )
    laser = 2 * math.pi * LASER_LINEWIDTH_HZ
    return 2 * math.sqrt(math.log(2) / math.pi) / math.sqrt(2 * laser**2 + doppler**2)


def _files() -> dict[tuple[str, str], str]:
    """Return the data file of each (table, method)."""
    out = {
        (basis, method): f"Kr/shekhtman_2020__{method}_{basis}.csv"
        for basis in BASIS_TABLES
        for method in METHODS
    }
    out["single_path", "sp"] = "Kr/shekhtman_2020__single_path.csv"
    return out


_Cell = tuple[str, str, str, str]
"""(table, method, line, cell), as :data:`ANOMALIES` names a printed cell."""


def _check_table_3(out: dict[str, object], failures: list[str]) -> None:
    """Check Table 3's lambda_L: the vacuum wavelength of half the level's energy."""
    e_g = ground_energy()
    worst = 0.0
    for line, (label, energy) in levels().items():
        computed = 2 * HC_EV_NM / (energy - e_g)
        worst = max(worst, abs(computed - float(line)))
        if abs(computed - float(line)) > VACUUM_TOLERANCE_NM:
            failures.append(f"Table 3 {label}: lambda_L {line}, computed {computed}")
    out["Table 3 lambda_L - 2hc/(E - E_G), worst (nm)"] = worst


def _check_row(
    key: tuple[str, str, str], cells: list[str], largest: float
) -> tuple[float, dict[_Cell, float], list[str]]:
    """Check one printed row: sigma / (sigma_o g), and the normalised ratio.

    Returns the ratio, the anomalies found (cell -> the computed number) and failures.
    """
    _, sigma_o, sigma, normalised = cells
    found: dict[_Cell, float] = {}
    failures = []
    ratio = _float(sigma) / (_float(sigma_o) * gaussian_peak(float(key[2])))
    if not RATIO_RANGE[0] <= ratio <= RATIO_RANGE[1]:
        cell = next((c for c in ("sigma", "sigma_o") if (*key, c) in ANOMALIES), "")
        if cell:
            found[*key, cell] = ratio
        else:
            failures.append(f"{key}: ratio {ratio:.4f}")
    computed = _float(sigma) / largest * (ANOMALIES.get((*key, "sigma")) or 1.0)
    tolerance = NORMALISED_TOLERANCE[0] * float(normalised) + NORMALISED_TOLERANCE[1]
    if abs(computed - float(normalised)) > tolerance:
        if (*key, "normalised") in ANOMALIES:
            found[*key, "normalised"] = computed
        else:
            failures.append(f"{key}: normalised {normalised}, computed {computed:.4f}")
    return ratio, found, failures


def checks() -> dict[str, object]:
    """Return the checks of every table against its own columns, Eq. (9) and Table 3.

    Raises if a cell outside :data:`ANOMALIES` fails, or if an anomaly is not one.
    """
    out: dict[str, object] = {}
    failures: list[str] = []
    _check_table_3(out, failures)
    found: dict[_Cell, float] = {}
    ratios = []
    for (name, method), printed_rows in rows().items():
        largest = max(_float(cells[2]) for cells in printed_rows)
        for cells in printed_rows:
            ratio, anomalies, failed = _check_row(
                (name, method, cells[0]), cells, largest
            )
            found |= anomalies
            failures += failed
            if not any(cell[3] != "normalised" for cell in anomalies):
                ratios.append(ratio)
    mean = sum(ratios) / len(ratios)
    out["sigma/(sigma_o g) of the regular cells (min, mean, max, n)"] = (
        min(ratios),
        mean,
        max(ratios),
        len(ratios),
    )
    for key, factor in ANOMALIES.items():
        if key not in found:
            failures.append(f"{key} is listed as an anomaly but agrees with its row")
            continue
        out[f"anomaly {key}"] = found[key]
        if factor is not None and abs(found[key] * factor / mean - 1) > 0.02:
            failures.append(f"{key}: a factor {factor} does not restore it")
    # Table 6's sigma_o^(2) at 193.494 nm from its sigma^(2).
    restored = _float(rows()["basis2", "qd"][1][2]) / (gaussian_peak(193.494) * mean)
    out["Table 6 QD 193.494 nm sigma_o from sigma^(2) (cm4)"] = restored
    if abs(restored / CORRECTED_BASIS2_193 - 1) > 0.005:
        failures.append(f"Table 6 QD 193.494: sigma^(2)/g gives {restored:.3e}")
    # The value the text quotes.
    basis, method, line, printed = TEXT_SIGMA_214
    cells = next(c for c in rows()[basis, method] if c[0] == line)
    out["text: Table 7 QD+OS sigma_o at 214.769 nm"] = (_value(cells[1]), printed)
    if _value(cells[1]) != printed:
        failures.append(f"text quotes {printed}, Table 7 prints {cells[1]}")
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return sigma_o^(2) of each table and method, after the checks."""
    checks()
    named = levels()
    files: dict[str, list[Point]] = {}
    for key, printed_rows in rows().items():
        points = []
        for cells in printed_rows:
            label = named[cells[0]][0]
            points.append(
                Point(label, cells[0], _value(cells[1]), "", " | ".join(cells))
            )
        files[_files()[key]] = points
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
