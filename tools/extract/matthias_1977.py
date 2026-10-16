# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Matthias et al. (1977): resonance-line lifetimes and f values of Kr and Xe.

Chem. Phys. Lett. 52, 239-244 (1977) = PDF pages 1-6 (journal page = PDF page + 238; a
scanned, raster PDF). Tables 1 (Kr) and 2 (Xe) (PDF p. 4, journal p. 242; library tables
``tbl_a11562a4d474`` and ``tbl_b97e7277fc92``) compile the lifetimes and absorption
oscillator strengths of the two resonance lines, ``^3P_1 -> ^1S_0`` (columns ``^3tau``,
``^3f``) and ``^1P_1 -> ^1S_0`` (``^1tau``, ``^1f``), one row per method and reference.
Taken: the row "this work" (resonance fluorescence after synchrotron excitation,
extrapolated to zero pressure), its two lifetimes and two f values per gas, one file
per gas and kind, the row named ``3P1`` or ``1P1`` and ``x`` the column heading as the
library extracted it (Table 1's "1263 Angstrom" is corrected by an erratum in the
record). The other rows are other papers' values (reproduced, not taken). The Xe 5d
J = 1 lifetime is printed only in the text (an inline dataset of the record). Every
cell was checked against the table images (2026-10-09).

**Checks** (:func:`checks`; the script fails if one is missed): every row of both
tables against the Ladenburg relation ``f = eps0 m_e c lambda^2 (g_u / g_l) / (2 pi
e^2 tau)``, ``g_u / g_l = 3``, at the NIST vacuum wavelengths of the lines
(``src/optmat/data/external/nist_asd/<gas>_I_lines.tsv``): every other row's pair
agrees to within 1 % (2.3 % for the 0.142 of Table 1's Ref. [15], printed to two
significant figures of its lifetime), and the row "this work" does not, except Kr
``^1P_1``: its printed f is 0.963 (Kr ``^3P_1``), 0.937 (Xe ``^3P_1``) and 1.043 (Xe
``^1P_1``) times the f of its own printed lifetime.
"""

import csv
import math
import re
from pathlib import Path

from _library import Point, number, table, write

from optmat.constants import E_CHARGE, EPS0, M_E, C

DIRECTORY = "Matthias_1977_TimeResolvedVUV"
SOURCE = "matthias_1977"
TABLES = {"Kr": "tbl_a11562a4d474", "Xe": "tbl_b97e7277fc92"}
NIST = Path(__file__).resolve().parents[2] / "src/optmat/data/external/nist_asd"
# The NIST level of each line's upper state: configuration's core and outer shell.
UPPER = {
    "Kr": {"3P1": ("(2P*<3/2>).5s", "1"), "1P1": ("(2P*<1/2>).5s", "1")},
    "Xe": {"3P1": ("(2P*<3/2>).6s", "1"), "1P1": ("(2P*<1/2>).6s", "1")},
}
# Columns of each line: lifetime, f.
COLUMNS = {"3P1": (0, 1), "1P1": (2, 3)}
_VALUE = re.compile(r"^\s*([\d.]+)\s*(?:[±+÷-]\s*([\d.]+))?\s*$")
# The row 'this work': printed f over the f of its printed lifetime, to 0.001.
THIS_WORK = {("Kr", "3P1"): 0.963, ("Kr", "1P1"): 1.004, ("Xe", "3P1"): 0.937}
THIS_WORK |= {("Xe", "1P1"): 1.043}


def _split(cell: str) -> tuple[str, str]:
    """Return the value and the uncertainty printed in a cell ("" when none)."""
    match = _VALUE.match(cell)
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    return number(match.group(1)), match.group(2) or ""


def extract() -> dict[str, list[Point]]:
    """Return the lifetimes and f values of the row 'this work', per gas."""
    files = {}
    for gas, table_id in TABLES.items():
        rows = table(DIRECTORY, table_id)
        header, body = rows[0], rows[1:]
        mine = [cells for cells in body if cells[5].strip() == "this work"]
        if len(mine) != 1 or header[5] != "Ref.":
            msg = f"{gas}: {len(mine)} rows 'this work' under {header!r}"
            raise ValueError(msg)
        cells = mine[0]
        printed = " | ".join(cells)
        for kind, index in (("tau", 0), ("f", 1)):
            points = []
            for row, columns in COLUMNS.items():
                value, sigma = _split(cells[columns[index]])
                x = header[columns[index]]
                points.append(Point(row, x, value, sigma, printed))
            files[f"{gas}/{SOURCE}__{kind}_{gas}.csv"] = points
    return files


def _wavelength_nm(gas: str, line: str) -> float:
    """Return the NIST vacuum wavelength of a resonance line, from its levels."""
    core, j = UPPER[gas][line]
    path = NIST / f"{gas}_I_lines.tsv"
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for entry in reader:
            cells = {key: (text or "").strip('"') for key, text in entry.items()}
            lower = cells["Ei(cm-1)"]
            ground = bool(lower) and float(lower) == 0
            if ground and cells["conf_k"].endswith(core) and cells["J_k"] == j:
                return 1e7 / float(cells["Ek(cm-1)"])
    msg = f"{gas} {line}: no NIST line"
    raise LookupError(msg)


def checks() -> None:
    """Check every row of both tables against the Ladenburg relation (module notes)."""
    factor = EPS0 * M_E * C / (2 * math.pi * E_CHARGE**2) * 3
    for gas, table_id in TABLES.items():
        for cells in table(DIRECTORY, table_id)[1:]:
            for line, (left, right) in COLUMNS.items():
                if not cells[left].strip() or not cells[right].strip():
                    continue
                tau = float(_split(cells[left])[0]) * 1e-9
                f = float(_split(cells[right])[0])
                wavelength = _wavelength_nm(gas, line) * 1e-9
                ratio = f / (factor * wavelength**2 / tau)
                if cells[5].strip() == "this work":
                    if round(ratio, 3) != THIS_WORK[gas, line]:
                        msg = f"{gas} {line} this work: ratio {ratio:.4f}"
                        raise ValueError(msg)
                elif abs(ratio - 1) > (0.025 if cells[5] == "[15]" else 0.01):
                    msg = f"{gas} {line} {cells[5]}: f / f(tau) = {ratio:.4f}"
                    raise ValueError(msg)


def main() -> None:
    """Check the tables and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
