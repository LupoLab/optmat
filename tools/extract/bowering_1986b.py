# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Böwering, Bruce and Keto (1986), part II: lifetimes and quenching of Xe 6p.

J. Chem. Phys. 84, 715-726 (1986); journal p. 715 is PDF page 2 (PDF page 1 is the AIP
download cover), so the journal page is the PDF page plus 713. A scanned (raster) PDF:
every cell of Tables I-V was checked against the table images (2026-10-09); the cells
come from the library's ``.json``. Rows are the levels in the paper's Paschen notation
(``2p5``) for 6p, and ``6s'[1/2]1``, ``6s'[1/2]0`` for 1s2 and 1s3 (Table III's names).
Taken, the column 'This work' only:

* **Table I** (PDF p. 6, journal p. 719, ``tbl_f63037673f36``), radiative lifetimes (ns)
  of 2p5, 2p6 and 2p9. The 2p9 cell holds two stacked lines, '30 ± 4' and '36 ± 4^a'
  (footnote a: 'Extracted for data from 1 to 10 Torr'); the library's cell runs them
  together, '30 ± 436 ± 4^{a}', and the script splits it as the table image shows; the
  second goes to a file of its own.
* **Table II** (PDF p. 6, journal p. 719, ``tbl_b0cdf9f5ecef``), total quenching by
  xenon: binary (10^-12 cm^3/s) for the six 6p levels, tertiary (10^-31 cm^6/s) for 2p5.
* **Table III** (PDF p. 9, journal p. 722, ``tbl_29410d84c1ff``), the coupled rates of
  6s'[1/2]1 and 6s'[1/2]0 (footnote a, 'Coupled with Xe 2p10'); each row has its own
  power of ten, so each is a file with its own scale. The 2p10 row repeats Table II.
* **Table IV** (PDF p. 10, journal p. 723, ``tbl_c338b72fc6b5``), the coupled-state
  model's state-to-state rates k12, k13, k10 (from 2p10), k21 (from 6s'[1/2]1) and k31
  (from 6s'[1/2]0), one file each (the channel and the power of ten differ). k12 and k21
  are printed with a kappa, the letter of the equilibrium coefficient (record notes).
  k1tot (Table II's) and k2tot, k3tot (Ref. 4's) are not taken.

**Checks** (:func:`checks`; the script fails if one is missed): the counts; Table IV
against the model of Sec. V with the equilibrium coefficients kappa = (g_l / g_u)
exp(Delta E / kT) at the stated 294 K from NIST's level energies (1.509 and 63.4 against
the printed 1.51 and 63.5): k12 / k21 and k13 / k31, k10 = k1tot - k12 - k13 (0.93 x
10^-10), k21 = k2tot; Eq. (10), k12 / (1 - k2c / k21) = 2.74 x 10^-10 against 'on the
order of' 2.8 x 10^-10; Eq. (8), tau1 = tau2 / (1 + kappa) = 43.8 ns from tau2 = 110 ns
against the record's inline 43.8 ± 8.0; and Table V's cross sections against Table II's
rates, sigma = k / v, with one mean relative speed v (2.97 x 10^4 cm/s, within 1.5 %
for every level), which is the speed at 273 K, not at the measurement's 296 K.
"""

import math
import re
import tomllib

from _json_table import json_table
from _library import DATA, Point, number, write

from optmat.constants import ATOMIC_MASS, HBAR, K_B, C

DIRECTORY = "Bowering_1986_CollisionalDeactivationTwo_5bde1b"
SOURCE = "bowering_1986b"
TABLE_I = "tbl_f63037673f36"
TABLE_II = "tbl_b0cdf9f5ecef"
TABLE_III = "tbl_29410d84c1ff"
TABLE_IV = "tbl_c338b72fc6b5"
TABLE_V = "tbl_4c4f96325404"
# Table I's 2p9 cell: two stacked lines on the table image.
STACKED = {"30 ± 436 ± 4^{a}": ("30 ± 4", "36 ± 4^{a}")}
_PASCHEN = re.compile(r"^(?:\\text\{Xe \})?2p_\{(\d+)\}$")
_VALUE = re.compile(r"^([\d.]+) (?:\\pm|±) ([\d.]+)(?:\^\{?a\}?)?$")
_POWER = re.compile(r"^\(([\d.]+) \\pm ([\d.]+)\) \\times 10\^\{(-\d+)a?\}$")
_MODEL = re.compile(r"^\\?(?:kappa|k)_\{(\d+)\} = ([\d.]+) \\times 10\^\{(-\d+)\}")
# Table IV's index -> (file stem, row: the reactant).
MODEL = {
    "12": ("k12", "2p10"),
    "13": ("k13", "2p10"),
    "10": ("k10", "2p10"),
    "21": ("k21", "6s'[1/2]1"),
    "31": ("k31", "6s'[1/2]0"),
}
# NIST level energies (cm^-1): 2p10, 6s'[1/2]1 (1s2), 6s'[1/2]0 (1s3).
E_2P10, E_1S2, E_1S3 = 77269.145, 77185.041, 76196.767
XENON_U = 131.293


def _file(stem: str) -> str:
    return f"Xe/{SOURCE}__{stem}_Xe.csv"


def _value(cell: str) -> tuple[str, str]:
    match = _VALUE.match(cell.strip())
    if match is None:
        msg = f"unexpected value {cell!r}"
        raise ValueError(msg)
    return number(match.group(1)), match.group(2)


def _lifetimes() -> dict[str, list[Point]]:
    """Return Table I's 'This work': the main lifetimes and 2p9 from 1-10 Torr."""
    grid = json_table(DIRECTORY, TABLE_I)
    row = next(cells for cells in grid if cells[0] == "This work")
    main, low = [], []
    for header, cell in zip(grid[0][2:], row[2:], strict=True):
        level = "2p" + _PASCHEN.match(header).group(1)  # type: ignore[union-attr]
        lines = STACKED.get(cell, (cell,))
        if "".join(lines) != cell:
            msg = f"split of {cell!r}"
            raise ValueError(msg)
        for k, line in enumerate(lines):
            point = Point(level, level, *_value(line), f"{header}: {line}")
            (low if k else main).append(point)
    return {_file("tau"): main, _file("tau_low"): low}


def _quenching() -> dict[str, list[Point]]:
    """Return Table II's 'This work': binary and tertiary rates."""
    grid = json_table(DIRECTORY, TABLE_II)
    section = ""
    two, three = [], []
    for cells in grid[2:]:
        if not cells[0]:
            section = cells[1]
            continue
        level = "2p" + _PASCHEN.match(cells[0]).group(1)  # type: ignore[union-attr]
        printed = " | ".join(cells)
        (three if section else two).append(
            Point(level, level, *_value(cells[1]), printed)
        )
    return {_file("kq_self"): two, _file("kq3_self"): three}


def _coupled() -> dict[str, list[Point]]:
    """Return Table III's coupled rates of the two 6s' levels, one file each."""
    files = {}
    for cells in json_table(DIRECTORY, TABLE_III)[2:]:
        level = cells[0].removeprefix("\\text{Xe } ").replace("_", "")
        match = _POWER.match(cells[1])
        if match is None:
            msg = f"unexpected cell {cells[1]!r}"
            raise ValueError(msg)
        stem = {"6s'[1/2]1": "kc_1s2", "6s'[1/2]0": "kc_1s3"}[level]
        point = Point(level, level, match.group(1), match.group(2), " | ".join(cells))
        files[_file(stem)] = [point]
    return files


def _scales() -> dict[str, int]:
    """Return the power of ten of each coupled and model file, as printed."""
    powers = {}
    for cells in json_table(DIRECTORY, TABLE_III)[2:]:
        match = _POWER.match(cells[1])
        stem = "kc_1s2" if "6s'[1/2]_1" in cells[0] else "kc_1s3"
        powers[stem] = int(match.group(3))  # type: ignore[union-attr]
    for (cell,) in json_table(DIRECTORY, TABLE_IV):
        match = _MODEL.match(cell)
        if match is not None and match.group(1) in MODEL:
            powers[MODEL[match.group(1)][0]] = int(match.group(3))
    return powers


def _model() -> dict[str, list[Point]]:
    """Return Table IV's state-to-state rates, one file each."""
    files = {}
    for (cell,) in json_table(DIRECTORY, TABLE_IV):
        if "tot" in cell:  # k_1tot (Table II's), k_2tot and k_3tot (Ref. 4's)
            continue
        match = _MODEL.match(cell)
        if match is None or match.group(1) not in MODEL:
            msg = f"unexpected cell {cell!r}"
            raise ValueError(msg)
        stem, row = MODEL[match.group(1)]
        files[_file(stem)] = [Point(row, row, number(match.group(2)), "", cell)]
    return files


def extract() -> dict[str, list[Point]]:
    """Return every file of the record."""
    return _lifetimes() | _quenching() | _coupled() | _model()


def _record() -> dict:
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        return tomllib.load(handle)


def checks() -> None:
    """Check the counts, Table IV's model, Eqs. (8) and (10), and Table V."""
    files = extract()
    counts = [len(points) for points in files.values()]
    if counts != [3, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1]:
        msg = f"counts {counts}"
        raise ValueError(msg)
    # The record's scales must be the printed powers of ten.
    scales = {d["id"]: d.get("scale", 1.0) for d in _record()["datasets"]}
    rate: dict[str, float] = {}
    for stem, power in _scales().items():
        if not math.isclose(scales[f"{stem}_Xe"], 10.0**power):
            msg = f"scale of {stem}: {scales[f'{stem}_Xe']} against 1e{power}"
            raise ValueError(msg)
        rate[stem] = float(files[_file(stem)][0].value) * 10.0**power
    total = {p.row: float(p.value) * 1e-12 for p in files[_file("kq_self")]}
    kt = K_B * 294.0 / (2 * math.pi * HBAR * C * 100)  # cm^-1
    kappa12 = math.exp((E_2P10 - E_1S2) / kt)
    kappa13 = math.exp((E_2P10 - E_1S3) / kt) / 3
    model = {
        "kappa12 (1.51)": (kappa12, 1.51, 0.005),
        "kappa13 (63.5)": (kappa13, 63.5, 0.003),
        "k12 / k21": (rate["k12"] / rate["k21"], kappa12, 0.005),
        "k13 / k31": (rate["k13"] / rate["k31"], kappa13, 0.02),
        "k10": (total["2p10"] - rate["k12"] - rate["k13"], rate["k10"], 0.01),
        "Eq. (10)": (rate["k12"] / (1 - rate["kc_1s2"] / rate["k21"]), 2.8e-10, 0.03),
    }
    for name, (computed, printed, tolerance) in model.items():
        if abs(computed / printed - 1) > tolerance:
            msg = f"{name}: {computed:.4g} against {printed}"
            raise ValueError(msg)
    inline = {
        d["id"]: float(d["points"][0][2])
        for d in _record()["datasets"]
        if d.get("table") == "text"
    }
    if abs(110.0 / (1 + kappa12) - inline["tau_2p10_coupled_Xe"]) > 0.05:
        msg = "Eq. (8)"
        raise ValueError(msg)
    # Table V: sigma (angstrom^2) = k / v with one mean relative speed v.
    racah = {"6p[1/2]_0": "2p5", "6p[3/2]_2": "2p6", "6p[3/2]_1": "2p7"}
    racah |= {"6p[5/2]_3": "2p8", "6p[5/2]_2": "2p9", "6p[1/2]_1": "2p10"}
    speeds = [
        total[racah[cells[0]]] / (float(cells[1]) * 1e-16)
        for cells in json_table(DIRECTORY, TABLE_V)[1:]
    ]
    mean = sum(speeds[1:]) / len(speeds[1:])
    if any(abs(v / mean - 1) > 0.015 for v in speeds):
        msg = f"Table V speeds {speeds}"
        raise ValueError(msg)
    mu = XENON_U / 2 * ATOMIC_MASS
    temperature = math.pi * mu * (mean / 100) ** 2 / (8 * K_B)
    if not 268 < temperature < 278:
        msg = f"Table V implies {temperature:.1f} K"
        raise ValueError(msg)


def main() -> None:
    """Check the tables and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
