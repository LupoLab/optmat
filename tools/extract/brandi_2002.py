# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Brandi, Hogervorst and Ubachs (2002): absolute energies of three Kr levels.

J. Phys. B 35, 1071-1084 (2002); journal p. 1071 is PDF page 1. Resonance-enhanced
photoionisation of krypton in a skimmed atomic beam (1VUV + 1UV and 2UV + 1UV),
isotopes separated by time of flight, the pulsed-dye-amplifier seed calibrated against
I2 hyperfine lines. Table 4 (``tbl_04c6e9828857``, p. 1079) prints, for 86Kr, the
measured energy of three levels above the 4p6 1S0 ground level (the transition wave
number from the ground level, cm^-1, with its uncertainty in parentheses) and the
difference from Kaufman and Humphreys (1969): one ``level_energy`` file per level,
the row and ``x`` the paper's jl-coupling label, the value and ``sigma`` as printed.
The library merges the energy and difference columns into one cell; :func:`_cell`
splits it.

:func:`checks` tests the paper against itself and NIST:

* the text's mean of the differences from Kaufman and Humphreys, '-0.0558 cm-1 with a
  standard deviation of 0.0005 cm-1' (p. 1079), from the three printed differences;
* the energies, moved from 86Kr to 84Kr by the paper's own isotope shifts (Table 2,
  ``tbl_7b89f13ee16c``, p. 1076: Delta^{84,86} = nu^84 - nu^86, MHz), against NIST's
  Kr I levels (``data/external/nist_asd/Kr_I_levels.tsv``, which gives no isotope),
  within 0.001 cm^-1. Unshifted, the 86Kr energies lie 0.004-0.007 cm^-1 above NIST's.

Every cell taken was checked against the table images (2026-10-09).
"""

import csv
import re
import statistics
from typing import Final

from _library import DATA, Point, number, table, write

DIRECTORY: Final = "Brandi_2002_HighResolutionVacuum"
ENERGIES: Final = "tbl_04c6e9828857"
"""Table 4: the 86Kr level energies and the differences from Kaufman and Humphreys."""
SHIFTS: Final = "tbl_7b89f13ee16c"
"""Table 2: the isotope shifts Delta^{84,A} = nu^84 - nu^A (MHz), A the column."""

LEVELS: Final = {
    "5d[1/2]1": ("5d12_1", "4s2.4p5.(2P*<3/2>).5d", "2[1/2]*", "1"),
    "6s[3/2]1": ("6s32_1", "4s2.4p5.(2P*<3/2>).6s", "2[3/2]*", "1"),
    "5p[1/2]0": ("5p12_0", "4s2.4p5.(2P*<3/2>).5p", "2[1/2]", "0"),
}
"""The paper's label: the dataset suffix and NIST's configuration, term and J."""

MEAN_DIFFERENCE: Final = -0.0558
"""'The average value of (nu - nu_KH) is -0.0558 cm-1' (p. 1079)."""
SPREAD: Final = 0.0005
"""'with a standard deviation of 0.0005 cm-1 (15 MHz)' (p. 1079)."""
MHZ_PER_CM: Final = 29979.2458
_CELL: Final = re.compile(r"^\s*([\d.]+)\((\d+)\)\s+(\S+)\s*$")


def _cell(text: str) -> tuple[str, str, str]:
    """Split Table 4's merged cell: the energy, its uncertainty and the difference."""
    match = _CELL.match(text)
    if match is None:
        msg = f"Table 4: cannot read {text!r}"
        raise ValueError(msg)
    energy, digits, difference = match.groups()
    decimals = len(energy.partition(".")[2])
    sigma = f"{int(digits) / 10**decimals:.{decimals}f}"
    return energy, sigma, number(difference)


def _rows() -> dict[str, tuple[list[str], str, str, str]]:
    """Return Table 4 by the paper's label: the cells, energy, sigma and difference."""
    rows = {}
    for cells in table(DIRECTORY, ENERGIES)[1:]:
        energy, sigma, difference = _cell(cells[2])
        rows[cells[0]] = (cells, energy, sigma, difference)
    if set(rows) != set(LEVELS):
        msg = f"Table 4 rows {sorted(rows)}, expected {sorted(LEVELS)}"
        raise ValueError(msg)
    return rows


def extract() -> dict[str, list[Point]]:
    """Return the three 86Kr level energies, one file each."""
    files = {}
    for label, (cells, energy, sigma, _) in _rows().items():
        suffix = LEVELS[label][0]
        point = Point(label, label, energy, sigma, " | ".join(cells))
        files[f"Kr/brandi_2002__level_{suffix}.csv"] = [point]
    return files


def _nist() -> dict[tuple[str, str, str], float]:
    """Return NIST's Kr I level energies (cm^-1) by configuration, term and J."""
    path = DATA / "external" / "nist_asd" / "Kr_I_levels.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {
        (row["Configuration"], row["Term"], row["J"]): float(row["Level (cm-1)"])
        for row in rows
        if row["Term"] != "Limit"
    }


def _shift_86() -> dict[str, float]:
    """Return Table 2's Delta^{84,86} (MHz) by the paper's label, this work's rows."""
    shifts = {}
    for cells in table(DIRECTORY, SHIFTS)[1:]:
        if cells[0] in LEVELS:
            shifts[cells[0]] = float(number(cells[1].partition("(")[0]))
    return shifts


def checks() -> None:
    """Check the printed mean difference and the energies against NIST; raise if not."""
    rows = _rows()
    differences = [float(row[3]) for row in rows.values()]
    mean = statistics.mean(differences)
    spread = statistics.stdev(differences)
    if abs(mean - MEAN_DIFFERENCE) > 0.00005 or abs(spread - SPREAD) > 0.00005:
        msg = f"Table 4 differences: mean {mean:.5f}, sd {spread:.5f}"
        raise ValueError(msg)
    nist = _nist()
    shifts = _shift_86()
    for label, (_, energy, _, _) in rows.items():
        _, configuration, term, j = LEVELS[label]
        # Delta^{84,86} = nu^84 - nu^86, so nu^84 = nu^86 + Delta^{84,86}.
        energy_84 = float(energy) + shifts[label] / MHZ_PER_CM
        reference = nist[configuration, term, j]
        if abs(energy_84 - reference) > 0.001:
            msg = f"{label}: 84Kr {energy_84:.4f} against NIST {reference:.4f} cm^-1"
            raise ValueError(msg)


def main() -> None:
    """Check the paper, then write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
