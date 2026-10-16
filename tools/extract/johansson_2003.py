# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Johansson et al. (2003): photoionisation of He 1s2p and 1s3p 1P, parallel.

Eur. Phys. J. D 22, 3-11 (2003); journal page = PDF page + 2. A born-digital PDF.

* **Table 1** (p. 8, ``tbl_c7aba2f99989``, JSON only: the state cells span rows): the
  absolute photoionisation cross sections sigma_par (Mb) of the excited 1s2p 1P and
  1s3p 1P levels, pump and probe polarisations parallel, against the probe wavelength
  (nm), with the excess energy E_exc (eV). One file per level; the row is the probe
  wavelength as printed, and ``printed`` keeps the table's row (the state carried down
  from its spanning cell). Checked against the table image.

The lifetimes (p. 7) are inline in the record (``tau``).

**Checks** (:func:`checks`; the script fails if one is missed): the counts; each E_exc
recomputed from the vacuum wavelength, NIST's level energy (the held
``He_I_levels.tsv``) and NIST's He ionisation energy (``ionization_energies.csv``),
within 0.005 eV; Eq. (7) from Eqs. (5)-(6) (the 6j symbols by Racah's formula) and its
sigma_par / sigma_perp at d/s = 13 against the printed 1.6; the probe rate sigma I / h
nu at 355 nm and 5 x 10^8 W/cm2 (14.8 ns^-1, against the printed 13: the record's
doubt) and 13 ns^-1 x 70 ps against the printed 0.9; the widths 3.5, 5.7 and 0.3 pm at
the 1s2p line against the printed 1.3, 2.0 and 0.11 meV; the inline lifetimes against
NIST's 1/A to the ground state (the held ``He_I_lines.tsv``), upper bounds.
"""

import csv
import math
import tomllib
from math import factorial

from _json_table import json_table
from _library import DATA, Point, number, write

from optmat.constants import E_CHARGE

DIRECTORY = "Johansson_2003_TwoColorPump"
SOURCE = "johansson_2003"
TABLE_1 = "tbl_c7aba2f99989"
HC_EV_NM = 1239.841984
"""hc in eV nm (CODATA 2018)."""
CM_EV = HC_EV_NM * 1e-7
"""One cm^-1 in eV."""
NIST = DATA / "external" / "nist_asd"
STATES = {"1s2p ^{1}\\text{P}": "1s2p", "1s3p ^{1}\\text{P}": "1s3p"}
"""Table 1's state cells and the file stems."""
CONFIGURATIONS = {"1s2p": "1s.2p", "1s3p": "1s.3p"}
"""NIST's configuration of each level (term 1P*, J = 1)."""


def _file(level: str) -> str:
    return f"He/{SOURCE}__sigma_{level}.csv"


def _rows() -> list[tuple[str, list[str]]]:
    """Return Table 1's body rows with the state carried down its spanning cell."""
    rows = json_table(DIRECTORY, TABLE_1)
    if rows[0] != [
        "Excited state",
        "\\lambda_{\\text{probe}}(nm)",
        "E_{\\text{exc}}(eV)",
        "\\sigma_{\\parallel}(Mb)",
    ]:
        msg = f"Table 1 header {rows[0]}"
        raise ValueError(msg)
    body, state = [], ""
    for cells in rows[1:]:
        state = cells[0] or state
        body.append((STATES[state], [state, *cells[1:]]))
    return body


def extract() -> dict[str, list[Point]]:
    """Return the two files of Table 1."""
    files: dict[str, list[Point]] = {_file(level): [] for level in STATES.values()}
    for level, cells in _rows():
        wavelength = number(cells[1])
        point = Point(wavelength, wavelength, number(cells[3]), "", " | ".join(cells))
        files[_file(level)].append(point)
    return files


def _level_ev(level: str) -> float:
    """Return NIST's energy (eV) of the level's 1P*1 term."""
    with (NIST / "He_I_levels.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (
                row["Configuration"].strip('"') == CONFIGURATIONS[level]
                and row["Term"].strip('"') == "1P*"
                and row["J"].strip('"') == "1"
            ):
                return float(row["Level (cm-1)"].strip('"')) * CM_EV
    msg = f"no NIST level for {level}"
    raise KeyError(msg)


def _ionisation_ev() -> float:
    """Return NIST's He I ionisation energy (eV)."""
    with (NIST / "ionization_energies.csv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["Sp. Name"].strip('="') == "He I":
                return float(row["Ionization Energy (eV)"].strip('="'))
    msg = "no He I ionisation energy"
    raise KeyError(msg)


def _a_to_ground(level: str) -> float:
    """Return NIST's A (s^-1) of the level's line to 1s2 1S0."""
    with (NIST / "He_I_lines.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (
                row["conf_i"].strip('"') == "1s2"
                and row["conf_k"].strip('"') == CONFIGURATIONS[level]
                and row["term_k"].strip('"') == "1P*"
            ):
                return float(row["Aki(s^-1)"].strip('"'))
    msg = f"no NIST line for {level}"
    raise KeyError(msg)


def _six_j(a: int, b: int, c: int, d: int, e: int, f: int) -> float:
    """Return the 6j symbol {a b c; d e f} by Racah's formula (integer arguments)."""

    def delta(x: int, y: int, z: int) -> float:
        return (
            factorial(x + y - z)
            * factorial(x - y + z)
            * factorial(-x + y + z)
            / factorial(x + y + z + 1)
        )

    pre = math.sqrt(delta(a, b, c) * delta(a, e, f) * delta(d, b, f) * delta(d, e, c))
    total = 0.0
    for t in range(a + b + c + d + e + f + 2):
        low = (t - a - b - c, t - a - e - f, t - d - b - f, t - d - e - c)
        high = (a + b + d + e - t, a + c + d + f - t, b + c + e + f - t)
        if min(low + high) < 0:
            continue
        denominator = math.prod(factorial(k) for k in low + high)
        total += (-1) ** t * factorial(t + 1) / denominator
    return pre * total


def _ratio(d_over_s: float) -> float:
    """Return sigma_par / sigma_perp from Eqs. (5)-(6) for sigma_d / sigma_s."""
    sigma = {0: 1.0, 2: d_over_s}
    iso = sum(sigma.values())
    # Eq. (5) with A20 = -sqrt(2): 1 + sqrt(10/3) sqrt(2) beta202 P2, beta202 of Eq. (6)
    weighted = sum(_six_j(1, 1, 2, 1, 1, L) * s for L, s in sigma.items())
    beta = 3 / iso * math.sqrt(3 / 5) * weighted
    factor = math.sqrt(10 / 3) * math.sqrt(2) * beta
    # Eq. (7): the same with 2 (sigma_s + sigma_d / 10) as the P2 coefficient
    if abs(factor * iso - 2 * (sigma[0] + sigma[2] / 10)) > 1e-12:
        msg = f"Eq. (7) against Eqs. (5)-(6): {factor * iso}"
        raise ValueError(msg)
    return (1 + factor) / (1 - factor / 2)


def _inline_tau() -> dict[str, float]:
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        record = tomllib.load(handle)
    (tau,) = (d for d in record["datasets"] if d["id"] == "tau")
    return {row: float(value) for row, _, value, _ in tau["points"]}


def checks() -> None:
    """Check the counts, E_exc, Eq. (7), the probe rate, the widths and lifetimes."""
    files = extract()
    counts = [len(points) for points in files.values()]
    if counts != [2, 3]:
        msg = f"counts {counts}"
        raise ValueError(msg)
    ionisation = _ionisation_ev()
    for level, cells in _rows():
        excess = HC_EV_NM / float(cells[1]) - (ionisation - _level_ev(level))
        if abs(excess - float(cells[2])) > 0.005:
            msg = f"E_exc {level} {cells[1]} nm: {excess:.4f} against {cells[2]}"
            raise ValueError(msg)
    if abs(_ratio(13.0) - 1.6) > 0.05 or abs(_ratio(20.0) - 1.5) > 1e-9:
        msg = f"sigma_par / sigma_perp: {_ratio(13.0)}, {_ratio(20.0)}"
        raise ValueError(msg)
    photon = HC_EV_NM / 355 * E_CHARGE
    rate = 16.6e-18 * 5e8 / photon / 1e9  # ns^-1
    if abs(rate - 14.83) > 0.01 or abs(13 * 0.070 - 0.9) > 0.02:
        msg = f"probe rate {rate:.3f} ns^-1"
        raise ValueError(msg)
    line_nm = HC_EV_NM / _level_ev("1s2p")
    for pm, mev in ((3.5, 1.3), (5.7, 2.0), (0.3, 0.11)):
        width = HC_EV_NM * pm * 1e-3 / line_nm**2 * 1e3
        if abs(width - mev) > 0.08 * mev:
            msg = f"{pm} pm: {width:.3f} meV against {mev}"
            raise ValueError(msg)
    tau = _inline_tau()
    for level, row in (("1s2p", "1s2p 1P1"), ("1s3p", "1s3p 1P1")):
        bound = 1e9 / _a_to_ground(level)
        if not tau[row] * 0.9 < bound < tau[row] * 1.05:
            msg = f"lifetime {row}: {tau[row]} ns against NIST {bound:.4f}"
            raise ValueError(msg)


def main() -> None:
    """Check the table and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
