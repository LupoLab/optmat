# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Sekino & Bartlett (1993): gamma of CO2, N2, H2 and CO by method and process.

Tables II (CO2, ``tbl_b4418052253e``), III (N2, ``tbl_bb8355cdb8fd``), IV (H2,
``tbl_91e0d5d3a5bb``, cells file only) and VII (CO, ``tbl_e455f4ab1704``); images
checked. Only the gamma rows are taken (alpha is linear data): the static components
gamma_xxxx, gamma_zzzz, gamma_xxzz and the static gamma_par, then THG, dc-SHG (ESHG),
IDRI (intensity-dependent refractive index) and EOKE (dc Kerr) gamma_par at omega =
0.0656 a.u. One file per row (process and component); the row ids are the methods
(columns SCF/TDHF, MBPT(2), CCSD, CCSD(T)). The EXP column (experiment, reproduced) and
H2's 'Exact' column (Bishop et al.'s values, reproduced) give one file per row.

Structural notes from the images: the component printed "gamma_xxxz" in the library's
cells is gamma_xxzz in all four tables. Table IV prints two rows per quantity (basis
set 1 [3s2p] above, basis set 2 [3s2p1d] below, footnote a): row ids "<method> basis1"
and "<method> basis2". Footnote markers (^f, ^{\mathrm{g}}) are stripped from the
values and kept in ``printed``.
"""

import re

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Sekino_1993_MolecularHyperpolarizabilities"
TABLES = {  # gas -> (table id, CSV available)
    "CO2": ("tbl_b4418052253e", True),
    "N2": ("tbl_bb8355cdb8fd", True),
    "H2": ("tbl_91e0d5d3a5bb", False),
    "CO": ("tbl_e455f4ab1704", True),
}
METHODS = ("scf_tdhf", "mbpt2", "ccsd", "ccsd_t")
ROWS = (  # printed label prefix -> (file name, wavelength in hartree)
    ("\\gamma_{xxxx}", "static_xxxx", "0"),
    ("\N{GREEK SMALL LETTER GAMMA}_{xxxx}", "static_xxxx", "0"),
    ("\\gamma_{zzzz}", "static_zzzz", "0"),
    ("\N{GREEK SMALL LETTER GAMMA}_{zzzz}", "static_zzzz", "0"),
    ("\\gamma_{xxxz}", "static_xxzz", "0"),
    ("\N{GREEK SMALL LETTER GAMMA}_{xxxz}", "static_xxzz", "0"),
    ("\\gamma_{||}", "static_par", "0"),
    ("\N{GREEK SMALL LETTER GAMMA}_{||}", "static_par", "0"),
    ("THG", "thg", "0.0656"),
    ("dc-SHG", "eshg", "0.0656"),
    ("IDRI", "kerr", "0.0656"),
    ("EOKE", "dc_kerr", "0.0656"),
)
_MARK = re.compile(r"\^\{?(?:\\mathrm\{)?[a-m]\}?\}?$")
_EXPT = re.compile(r"^([\d.]+)(?:±([\d.]+))?\^\{?(?:\\mathrm\{)?([a-m])\}?\}?$")


def _row(label: str) -> tuple[str, str] | None:
    for prefix, name, omega in ROWS:
        if label.startswith(prefix):
            return name, omega
    return None


def extract() -> dict[str, list[Point]]:
    """Return the gamma rows of Tables II, III, IV and VII."""
    files: dict[str, list[Point]] = {}
    for gas, (table_id, has_csv) in TABLES.items():
        grid = (
            table(DIRECTORY, table_id) if has_csv else json_table(DIRECTORY, table_id)
        )
        prefix = f"{gas}/sekino_1993__"
        label = ""
        basis = ""
        for cells in grid[1:]:
            if gas == "H2":  # two rows per quantity: basis 1, then basis 2
                basis = "basis2" if not cells[0] else "basis1"
                label = cells[0] or label
            else:
                label = cells[0]
            found = _row(label)
            if found is None:  # an alpha row
                continue
            name, omega = found
            printed = " | ".join(cells)
            methods = METHODS[:3] if gas == "H2" else METHODS
            for method, cell in zip(methods, cells[1:], strict=False):
                value = _MARK.sub("", cell)
                float(value)
                row = f"{method} {basis}" if basis else method
                files.setdefault(f"{prefix}{name}_{gas}.csv", []).append(
                    Point(row, omega, value, "", printed)
                )
            if gas == "H2" and cells[4]:  # the 'Exact' column (Bishop et al.)
                value = _MARK.sub("", cells[4])
                float(value)
                files.setdefault(f"{prefix}exact_{name}_{gas}.csv", []).append(
                    Point(name, omega, value, "", printed)
                )
            expt = cells[5] if len(cells) > 5 else ""
            if expt:
                match = _EXPT.match(expt)
                if match is None:
                    msg = f"unparsed experimental cell {expt!r}"
                    raise ValueError(msg)
                files.setdefault(f"{prefix}expt_{name}_{gas}.csv", []).append(
                    Point(
                        f"expt^{match.group(3)}",
                        omega,
                        match.group(1),
                        match.group(2) or "",
                        printed,
                    )
                )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
