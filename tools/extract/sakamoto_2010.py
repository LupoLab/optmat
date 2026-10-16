# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Sakamoto et al. (2010), NIFS-DATA-109: the discrete oscillator strengths.

The report tabulates, per species, the discrete oscillator strengths below (and some
resonances above) the first ionisation potential, then the continuum df/dE. The
continuum tables are kept as the NIFS database's text files, ``<gas>_table.txt.gz``
(retrieved on 2026-09-25 and attached byte for byte to the report in the library, with
their URLs); they are not extracted here. The same files also
hold the discrete tables, and those are what this script reads, for He, Ne, Ar, H2,
N2, O2, CO and N2O
(CO2, O3, CH4 and SF6 have none).

Why the NIFS text and not the library's cells: the library's tables of this report are
OCR of pages printed sideways. Most continuum and several discrete tables come out as
reversed strings (``20+E0056 .7`` for ``7.6500E+02``), others merge two tables (the
argon and carbon monoxide pages) or misread digits (helium Table I prints 2.4452E+01,
6.6800E-04, 6.4816E+01 and 6.5108E+01; the cells read 2.4425E+01, 6.8800E-04,
6.416E+01 and a repeated 6.5181E+01). The PDF itself (``original/source_01``) is
typeset text: its text layer (``pdftotext -layout``, report pages 14-258) was compared
with every NIFS file on 2026-09-25. Every number of every file (discrete and
continuum, 5 900 rows) is on the PDF page, each row as one run of cells on one printed
line, and the multisets of numbers agree exactly, species by species. The page images
of every discrete table, and a continuum page per species, were checked by eye.

What the NIFS text loses, restored here (structural, from the page images):

- the footnote marks: helium Table I's last row (``2.0760E-03^a)``: the sum over
  n >= 21) and argon Table I's six (``a)`` ... ``f)``: series sums and the continuum
  between the two thresholds). They are appended to the ``printed`` column as
  ``[mark a)]``;
- the horizontal rules that split a table into groups: neon Table II (the 2s -> np
  resonances at 45-49 eV above the K-edge resonances at 867-870 eV), argon Table II
  (the negative window resonances at 26.6-29.1 eV above the K-edge line at 3203.3 eV)
  and oxygen Table I (the sub-ionisation intervals above the K-edge line at 530.9 eV).
  Each group is its own file.

Rows are numbered ``<table>-<n>``, n counting the table's rows in the NIFS order, which
is the page's reading order (down the left column pair, then the right). Energies (eV)
and strengths are the printed strings; the wavelength column (Angstrom, computed from
the energy by the authors) goes to ``printed`` only. Interval tables (O2 Table I, N2O
Table I) give the integrated strength over an energy range: the row is the range as
printed and the abscissa its lower edge, as for ``kumar_2010_supplement``.
"""

import gzip
import re
from typing import Final

from _library import Point, attachment, write

PAPER: Final = "Sakamoto_2010_OscillatorStrengthSpectra"
_ROMAN: Final = re.compile(r"^Table (?P<table>[IVX]+)\.")
_NUMBER: Final = r"-?\d*\.\d+E[-+]\d+"
_LINE: Final = re.compile(rf"^\s*({_NUMBER})\s+({_NUMBER})\s+({_NUMBER})\s*$")
_INTERVAL: Final = re.compile(
    r"^\s*(?P<low>\d+\.\d+)\s*-\s*(?P<high>\d+\.\d+)(?P<ip>\(IP\))?\s+"
    r"(?P<f>\d*\.\d+)\s+(?P<lam>\S+\s*-\s*\S+)\s*$"
)

# Footnote marks the page images print and the NIFS text drops: (gas, row) -> mark.
MARKS: Final = {
    ("He", "I-20"): "a)",
    ("Ar", "I-5"): "a)",
    ("Ar", "I-9"): "b)",
    ("Ar", "I-13"): "c)",
    ("Ar", "I-14"): "d)",
    ("Ar", "I-18"): "e)",
    ("Ar", "I-21"): "f)",
}

# (gas, table) -> [(dataset id, predicate on the energy in eV)]: one file per group
# of rows that the page separates by a rule.
GROUPS: Final = {
    ("He", "I"): [("discrete_He", None)],
    ("He", "II"): [("excess_He", None)],
    ("Ne", "I"): [("discrete_Ne", None)],
    ("Ne", "II"): [("resonance_2s_Ne", (0.0, 100.0)), ("k_edge_Ne", (100.0, 1e4))],
    ("Ar", "I"): [("discrete_Ar", None)],
    ("Ar", "II"): [("window_Ar", (0.0, 100.0)), ("k_edge_Ar", (100.0, 1e4))],
    ("H2", "I"): [("lyman_H2", None)],
    ("H2", "II"): [("werner_H2", None)],
    ("H2", "III"): [("b_prime_H2", None)],
    ("H2", "IV"): [("d_H2", None)],
    ("H2", "V"): [("rydberg_H2", None)],
    ("H2", "VI"): [("predissociation_H2", None)],
    ("H2", "VII"): [("autoionization_H2", None)],
    ("N2", "I"): [("discrete_N2", None)],
    ("N2", "II"): [("k_edge_N2", None)],
    ("O2", "I"): [("intervals_O2", (0.0, 100.0)), ("k_edge_O2", (100.0, 1e4))],
    ("CO", "I"): [("discrete_CO", None)],
    ("CO", "II"): [("k_edge_CO", None)],
    ("N2O", "I"): [("intervals_N2O", None)],
}
# Expected row counts per dataset: a guard against a changed or misparsed file.
COUNTS: Final = {
    "discrete_He": 20,
    "excess_He": 18,
    "discrete_Ne": 19,
    "resonance_2s_Ne": 6,
    "k_edge_Ne": 4,
    "discrete_Ar": 21,
    "window_Ar": 7,
    "k_edge_Ar": 1,
    "lyman_H2": 37,
    "werner_H2": 14,
    "b_prime_H2": 9,
    "d_H2": 16,
    "rydberg_H2": 20,
    "predissociation_H2": 17,
    "autoionization_H2": 146,
    "discrete_N2": 47,
    "k_edge_N2": 3,
    "intervals_O2": 11,
    "k_edge_O2": 1,
    "discrete_CO": 19,
    "k_edge_CO": 2,
    "intervals_N2O": 5,
}


def _lines(gas: str) -> list[str]:
    """Return the lines of a species' NIFS text table."""
    path = attachment(PAPER, f"{gas}_table.txt.gz")
    with gzip.open(path, "rt", encoding="ascii") as handle:
        return handle.read().splitlines()


def _discrete(gas: str) -> dict[str, list[Point]]:
    """Return the discrete tables of one species, one file per dataset."""
    files: dict[str, list[Point]] = {}
    table = None
    counter: dict[str, int] = {}
    for line in _lines(gas):
        heading = _ROMAN.match(line.strip())
        if heading is not None:
            table = heading["table"]
            if (gas, table) not in GROUPS:
                break  # the continuum table: the discrete tables are done
            counter[table] = 0
            continue
        if table is None:
            continue
        line_match = _LINE.match(line)
        interval = _INTERVAL.match(line)
        if line_match is None and interval is None:
            continue
        counter[table] += 1
        row = f"{table}-{counter[table]}"
        if interval is not None:
            ip = interval["ip"] or ""
            label = f"{interval['low']}-{interval['high']}{ip}"
            energy, value = interval["low"], interval["f"]
            cells = [f"{interval['low']} - {interval['high']}{ip}", interval["f"]]
            cells.append(" ".join(interval["lam"].split()))
            row = f"{table}-{counter[table]} ({label})"
        else:
            assert line_match is not None
            energy, value, wavelength = line_match.groups()
            cells = [energy, value, wavelength]
        printed = f"Table {table} | " + " | ".join(cells)
        mark = MARKS.get((gas, row))
        if mark is not None:
            printed += f" [mark {mark}]"
        for dataset, window in GROUPS[gas, table]:
            if window is None or window[0] <= float(energy) < window[1]:
                relative = f"{gas}/sakamoto_2010__{dataset}.csv"
                files.setdefault(relative, []).append(
                    Point(row, energy, value, "", printed)
                )
                break
    return files


def extract() -> dict[str, list[Point]]:
    """Return every discrete-table file of the eight species that have them."""
    files: dict[str, list[Point]] = {}
    for gas in ("He", "Ne", "Ar", "H2", "N2", "O2", "CO", "N2O"):
        files |= _discrete(gas)
    for relative, points in files.items():
        dataset = relative.split("__")[1].removesuffix(".csv")
        if len(points) != COUNTS[dataset]:
            msg = f"{dataset}: {len(points)} rows, expected {COUNTS[dataset]}"
            raise ValueError(msg)
    if len(files) != len(COUNTS):
        msg = f"{len(files)} files, expected {len(COUNTS)}"
        raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
