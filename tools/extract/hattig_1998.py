# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hattig & Jorgensen (1998): CC gamma of Ne, Ar, N2 and CH4 by process.

Every table used was compared with its image.

* Tables VIII (Ne), IX (Ar), XIII (N2), XV (CH4): CCSD gamma for the dc Kerr effect
  (gamma_K), DFWM (the intensity-dependent index), ESHG and THG, rows lambda/nm with
  omega/a.u. Files ``ccsd_<process>_<gas>``; the row id is lambda as printed and the
  wavelength column is omega (hartree). Table XV prints lambda "590.0" beside omega
  0.077357 a.u. (589.0 nm; Tables VII and XVI print 589.0 for the same row): the row id
  keeps the printed "590.0", the wavelength is omega.
* Tables IV (Ne), V (Ar), VI (N2), VII (CH4): THG by SCF, CCS, CC2 and CCSD. The CCSD
  columns of IV, V and VII repeat the THG columns of VIII, IX and XV (same numbers) and
  are not written twice; N2's (Table VI, at the equilibrium geometry) differs from
  Table XIII (vibrationally averaged) and is written as ``ccsd_re_thg_N2``.
* Table III: static gamma_0 by method. The "effect" rows (differences) are not taken;
  the rows quoting other papers (MCSCF, CCSD(T)^j, CCSD(T)-R12, QCISD(T)) go to
  ``static_others_<gas>``. N2 and CH4 have a second "CCSD unrelaxed" row (footnote c,
  d-aug-cc-pVTZ basis, on an unlabelled line): row id "CCSD unrelaxed (d-aug-cc-pVTZ)".
* Table XII: only the vibrationally averaged static <0|gamma_0|0> of N2 (915.0).
* Tables X (Ne), XI (Ar), XIV (N2), XVI (CH4): the scaled CCSD columns (the paper's
  best estimates, ``best_<process>_<gas>``) and the experimental columns
  (``expt_<process>_<gas>``, reproduced). The library's cells for Tables X and XVI are
  garbled: the plus-minus sign is read as "6", footnote letters run into the numbers,
  and in Table X the dc Kerr CCSD and Expt cells, and in Table XVI the THG CCSD and Expt
  cells, are merged. Those cells are mapped below to what the image shows; each
  mapping asserts the garbled text, so a change in the library stops the script. Table
  XIV's image prints "1327.5 +- 14..3" (a doubled decimal point): read as 14.3.
"""

import re

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Hattig_1998_CoupledClusterCalculations"
PROCESSES = ("dc_kerr", "kerr", "eshg", "thg")  # Kerr, DFWM, ESHG, THG columns
PROCESS_TABLES = {  # gas -> (table id, CSV available)
    "Ne": ("tbl_26296ab9e171", True),
    "Ar": ("tbl_c806d363f4a8", True),
    "N2": ("tbl_633f4bbb50b3", True),
    "CH4": ("tbl_5798b1c436d3", True),
}
THG_TABLES = {  # gas -> (table id, write the CCSD column)
    "Ne": ("tbl_05c0048e936a", False),
    "Ar": ("tbl_1d1b437d62c7", False),
    "N2": ("tbl_1b5718881686", True),
    "CH4": ("tbl_faa012f2f892", False),
}
TABLE_III = "tbl_ae487f4bcb5f"
TABLE_XII = "tbl_bee0d90000c1"
COMPARISON = {
    "Ne": "tbl_1b31198501c7",
    "Ar": "tbl_0741c20d0afb",
    "N2": "tbl_592e2af72196",
    "CH4": "tbl_c2b8213c9cf9",
}
# Garbled cells of Tables X and XVI -> (value, sigma, footnote) as the image shows.
GARBLED: dict[str, tuple[str, str, str]] = {
    # Table X (Ne)
    "111.060.8b": ("111.0", "0.8", "b"),
    "109.960.5b": ("109.9", "0.5", "b"),
    "126.464.1c": ("126.4", "4.1", "c"),
    "122.260.4b": ("122.2", "0.4", "b"),
    "123.560.5b": ("123.5", "0.5", "b"),
    "90645d": ("90", "45", "d"),
    "79.768e": ("79.7", "8", "e"),
    "108618f": ("108", "18", "f"),
    # Table XVI (CH4); dc Kerr column footnote b, ESHG column e (headers)
    "2855621": ("2855", "21", ""),
    "2967616": ("2967", "16", ""),
    "3070611c": ("3070", "11", "c"),
    "3132627": ("3132", "27", ""),
    "3236619": ("3236", "19", ""),
    "3286624": ("3286", "24", ""),
    "3360614d": ("3360", "14", "d"),
    "2603611": ("2603", "11", ""),
    "2833610": ("2833", "10", ""),
    "3311610": ("3311", "10", ""),
    "361369": ("3613", "9", ""),
    "4060623": ("4060", "23", ""),
    "4414624": ("4414", "24", ""),
}
# Merged cells: CCSD value plus an experimental value (Table X dc Kerr, Table XVI THG).
MERGED: dict[str, tuple[str, tuple[str, str, str]]] = {
    "110.2 10168a": ("110.2", ("101", "8", "a")),
    "4162 a 39976120": ("4162", ("3997", "120", "")),
}
_CLEAN = re.compile(r"^([\d.]+)±([\d.]+)(?:\^\{(?:\\text\{)?([a-h])\}?\}?)?$")
_CCSD = re.compile(
    r"^([\d.]+)(?:\s*a|\^\{(?:\\text\{)?a\}?\}|\^\{\\text\{a\}\}|\^\{a\})?$"
)


def _grid(table_id: str, has_csv: bool) -> list[list[str]]:
    return table(DIRECTORY, table_id) if has_csv else json_table(DIRECTORY, table_id)


def _processes() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for gas, (table_id, has_csv) in PROCESS_TABLES.items():
        grid = _grid(table_id, has_csv)
        for cells in grid[1:]:
            lam, omega = cells[0], cells[1]
            float(lam)
            float(omega)
            for process, value in zip(PROCESSES, cells[2:6], strict=True):
                float(value)
                files.setdefault(
                    f"{gas}/hattig_1998__ccsd_{process}_{gas}.csv", []
                ).append(Point(lam, omega, value, "", " | ".join(cells)))
    return files


def _thg_methods() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for gas, (table_id, with_ccsd) in THG_TABLES.items():
        grid = table(DIRECTORY, table_id)
        if grid[0][2:] != ["SCF", "CCS", "CC2", "CCSD"]:
            msg = f"unexpected header in {table_id}: {grid[0]!r}"
            raise ValueError(msg)
        methods = ("scf", "ccs", "cc2", "ccsd") if with_ccsd else ("scf", "ccs", "cc2")
        for cells in grid[1:]:
            for method, value in zip(methods, cells[2:], strict=False):
                float(value)
                # N2's Table VI CCSD is at r_e: not the averaged ccsd_thg_N2 of XIII.
                name = "ccsd_re" if method == "ccsd" else method
                files.setdefault(f"{gas}/hattig_1998__{name}_thg_{gas}.csv", []).append(
                    Point(cells[0], cells[1], value, "", " | ".join(cells))
                )
    return files


def _static() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_III)
    gases = ("Ne", "Ar", "N2", "CH4")
    others = ("MCSCF", "CCSD(T)^j", "CCSD(T)-R12", "QCISD(T)")
    files: dict[str, list[Point]] = {}
    label = ""
    for cells in grid[1:]:
        label = cells[0] or f"{label} (d-aug-cc-pVTZ)"
        if "effect" in label:
            continue
        name = "static_others" if label.startswith(others) else "static"
        row = re.sub(r"\^[a-l]$", "", label)
        for gas, cell in zip(gases, cells[1:], strict=True):
            if not cell:
                continue
            value = re.sub(r"\^[a-c]$", "", cell)
            float(value)
            files.setdefault(f"{gas}/hattig_1998__{name}_{gas}.csv", []).append(
                Point(row, "0", value, "", " | ".join(cells))
            )
    return files


def _vibrational_average() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_XII)
    (cells,) = [r for r in grid if r[0] == "^0ugu0&"]  # the <0|gamma|0> row
    float(cells[1])
    return {
        "N2/hattig_1998__vib_avg_static_N2.csv": [
            Point("<0|gamma|0>", "0", cells[1], "", " | ".join(cells))
        ]
    }


def _expt(cell: str) -> tuple[str, str, str]:
    if cell in GARBLED:
        return GARBLED[cell]
    cell = cell.replace("14..3", "14.3")  # the image's doubled decimal point
    match = _CLEAN.match(cell)
    if match is None:
        msg = f"unparsed experimental cell {cell!r}"
        raise ValueError(msg)
    return match.group(1), match.group(2), match.group(3) or ""


def _ccsd(cell: str) -> str:
    match = _CCSD.match(cell)
    if match is None:
        msg = f"unparsed CCSD cell {cell!r}"
        raise ValueError(msg)
    return match.group(1)


def _comparison() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for gas, table_id in COMPARISON.items():
        grid = json_table(DIRECTORY, table_id)
        for cells in grid[2:]:
            lam, omega = cells[0], cells[1]
            float(omega)
            printed = " | ".join(cells)
            pairs: list[tuple[str, str, str]] = []  # process, CCSD cell, expt cell
            if gas == "Ne":  # dc Kerr CCSD and Expt share column 2
                pairs = [("dc_kerr", cells[2], ""), ("eshg", cells[3], cells[4])]
                pairs.append(("thg", cells[5], cells[6]))
            elif gas == "CH4":  # THG CCSD and Expt share column 6
                pairs = [("dc_kerr", cells[2], cells[3]), ("eshg", cells[4], cells[5])]
                pairs.append(("thg", cells[6], ""))
            else:
                pairs = [
                    ("dc_kerr", cells[2], cells[3]),
                    ("eshg", cells[4], cells[5]),
                    ("thg", cells[6], cells[7]),
                ]
            for process, ccsd_cell, expt in pairs:
                ccsd = ccsd_cell
                if ccsd_cell in MERGED:
                    ccsd, merged = MERGED[ccsd_cell]
                    expt_values: tuple[str, str, str] | None = merged
                else:
                    expt_values = _expt(expt) if expt else None
                if ccsd:
                    files.setdefault(
                        f"{gas}/hattig_1998__best_{process}_{gas}.csv", []
                    ).append(Point(lam, omega, _ccsd(ccsd), "", printed))
                if expt_values is not None:
                    value, sigma, mark = expt_values
                    row = f"{lam}^{mark}" if mark else lam
                    files.setdefault(
                        f"{gas}/hattig_1998__expt_{process}_{gas}.csv", []
                    ).append(Point(row, omega, value, sigma, printed))
    return files


def extract() -> dict[str, list[Point]]:
    """Return every extracted file."""
    return (
        _processes()
        | _thg_methods()
        | _static()
        | _vibrational_average()
        | _comparison()
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
