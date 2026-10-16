# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bishop, Pipin & Cybulski (1991): averaged gamma of H2 and D2.

All images checked against the cells; no differences.

* Tables VI-XI (CSV): the averaged (295 K) electronic gamma_par and gamma_perp, the pure
  vibrational contributions and their totals, for the dc Kerr effect, ESHG and THG
  (THG: gamma_par only; gamma_perp = gamma_par / 3), H2 (VI-VIII) and D2 (IX-XI), at
  hbar omega / E_h = 0.00-0.05 and seven laser frequencies. One file per column:
  ``<process>_<component>_<part>_<gas>`` with part elec, vib or total. Row id and
  wavelength: omega as printed with the thin space dropped ("0.034 54" -> "0.03454").
* Table XIV (cells file): the ratios S = gamma_par^ESHG(total, H2) / (total, D2) and T =
  gamma_par^ESHG(total, H2) / gamma_par^ESHG(He), theory, experiment and "fitted
  experiment", against lambda in angstrom. The experimental uncertainty printed in
  parentheses on the last digits ("0.987(3)") is written as a sigma in the value's unit
  (0.003).
* Table XVI (CSV): gamma^K = (3/2)(gamma_par^Kerr(total) -
  gamma_perp^Kerr(total)) for H2 and D2 against lambda in angstrom.

Table IV (the molecular-frame components at R = 1.4 a0), Tables III, V and XII (alpha, B
and C) and Table XIII (the ratio R = gamma_par / gamma_perp) are not extracted (see the
record); Table XV's fits are formulas in the record.
"""

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Bishop_1991_TheoreticalInvestigationNonlinear"
PROCESS_TABLES = {  # table id -> (gas, process, has gamma_perp)
    "tbl_7138422de0e9": ("H2", "dc_kerr", True),
    "tbl_a26e8934041a": ("H2", "eshg", True),
    "tbl_4827aae9526c": ("H2", "thg", False),
    "tbl_25f1365d9669": ("D2", "dc_kerr", True),
    "tbl_cdf929a9c5ce": ("D2", "eshg", True),
    "tbl_b3b2cb2c1406": ("D2", "thg", False),
}
PARTS = ("elec", "vib", "total")
TABLE_XIV = "tbl_fad13f7ef088"
TABLE_XVI = "tbl_60de52729ff9"


def _prefix(gas: str) -> str:
    return f"{gas}/bishop_1991__"


def _processes() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for table_id, (gas, process, perp) in PROCESS_TABLES.items():
        grid = table(DIRECTORY, table_id)
        if not grid[0][1].startswith("\\overline{\\gamma}_{\\parallel}(\\text{elec})"):
            msg = f"unexpected header in {table_id}: {grid[0]!r}"
            raise ValueError(msg)
        columns = [("par", part) for part in PARTS]
        if perp:
            columns += [("perp", part) for part in PARTS]
        for cells in grid[1:]:
            omega = cells[0].replace(" ", "")
            float(omega)
            for (component, part), value in zip(columns, cells[1:], strict=True):
                float(value)
                name = f"{process}_{component}_{part}_{gas}"
                files.setdefault(f"{_prefix(gas)}{name}.csv", []).append(
                    Point(omega, omega, value, "", " | ".join(cells))
                )
    return files


def _parenthesised(cell: str) -> tuple[str, str]:
    """Split "0.987(3)" into ("0.987", "0.003"): sigma on the last digits."""
    value, _, rest = cell.partition("(")
    digits = rest.removesuffix(")")
    decimals = len(value.partition(".")[2])
    sigma = f"{int(digits) / 10**decimals:.{decimals}f}"
    float(value)
    return value, sigma


def _ratios() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_XIV)
    if grid[1] != ["", "Theory", "Expt.", "Theory", "Expt.", "Fitted Expt.^a"]:
        msg = f"unexpected Table XIV header: {grid[1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        lam = cells[0].replace(" ", "")
        printed = " | ".join(cells)
        for name, cell in (
            ("s_theory", cells[1]),
            ("s_expt", cells[2]),
            ("t_theory", cells[3]),
            ("t_expt", cells[4]),
            ("t_fitted", cells[5]),
        ):
            value, sigma = _parenthesised(cell) if "(" in cell else (cell, "")
            float(value)
            files.setdefault(f"{_prefix('H2')}ratio_{name}_H2.csv", []).append(
                Point(lam, lam, value, sigma, printed)
            )
    return files


def _kerr() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_XVI)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        lam = cells[0].replace(" ", "")
        for gas, value in (("H2", cells[1]), ("D2", cells[2])):
            float(value)
            files.setdefault(f"{_prefix(gas)}gamma_k_total_{gas}.csv", []).append(
                Point(lam, lam, value, "", " | ".join(cells))
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables VI-XI, XIV and XVI."""
    return _processes() | _ratios() | _kerr()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
