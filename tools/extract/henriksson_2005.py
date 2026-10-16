# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Henriksson et al. (2005): TDHF and TDDHF one- and two-photon strengths.

J. Chem. Phys. 122, 114106 (2005), pp. 114106-1 to -9 = PDF pages 2-10 (PDF page 1 is
the AIP download cover).

* **Tables I-V** (He, Ne, Ar, Kr, Xe; PDF pp. 5-7, pp. 114106-4 to -6; library tables
  ``tbl_930ae06840e1``, ``tbl_1fcee34b02cc``, ``tbl_37bc1a24b4ae``,
  ``tbl_42d892cee858``, ``tbl_e64a945b01a8``): per configuration and LS term, the
  experimental excitation energy (Moore's tables; not taken), the non-relativistic
  TDHF excitation energy and strength ``delta`` (a.u.), then in separate columns J and
  the relativistic TDDHF energy and strength. ``delta`` is ``delta_OPA`` (Eq. 2, the
  oscillator strength) or ``delta_TPA`` (Eq. 6) 'clear from the atomic selection
  rules' (footnote b): a level of the ground state's (even) parity is reached by two
  photons, an odd one by one. The parity is the configuration's (``sum l`` over the
  electrons), checked against the term's ``o`` mark.
* Files: per gas, TDHF and TDDHF, two-photon (``tdhf_tpa_<gas>``, ``tddhf_tpa_<gas>``)
  and one-photon (``tdhf_opa_<gas>``, ``tddhf_opa_<gas>``). TDHF rows are named by the
  configuration and term as printed (``2p5 3p 1D_2``; ``*`` marks an odd term). The
  TDDHF columns are sorted by the computed energy within each configuration's block of
  rows and carry their own J: they do not line up with the LS terms beside them, so a
  TDDHF row is named by what it prints, ``J = 2, 20.090 eV``. Selection-rule zeros
  (spin-forbidden in TDHF, Delta J = 1 or 3 in TDDHF) are kept as printed.
* ``delta`` cells print digit groups with a thin space (``0.164 47``); the value drops
  it, ``printed`` keeps the row's cells.

Table VI (Rn) is not taken: optmat has no radon.

**The neon table** (Table II, ``tbl_1fcee34b02cc``): the library's extraction merged
each configuration's rows into one cell per column (``'000011.060000012.850'``). The
split, :data:`NEON`, was read from the table image (``tables/tbl_1fcee34b02cc.jpg``)
cell by cell (2026-10-08); the script checks that each column's cells, joined, are
exactly the merged cell, and :func:`checks` finds every row's six numbers, in order,
in the PDF's text layer of the page.

**Checks** (:func:`checks`; the script fails if one is missed): every row of every
table in the PDF's text layer; parity against the ``o`` marks; the TDDHF J of each
block against the LS terms' J (equal for He to Kr; Xe's differ, reported); the text's
numbers (pp. 114106-4, -8): the Ne and Ar 1D2 strengths from the quoted tensor
elements, ``delta = (4/3) S_xy^2``, the He 2 1S0 basis comparison, Ne's quoted
oscillator strengths, argon's 9 % leakage; the sum of the TDDHF strengths of each J
against TDHF ('virtually the same', Ne and Ar). Reported: the text's quotation of
Hättig et al.'s values (hattig_1998b).
"""

import re
from dataclasses import dataclass

import pypdfium2
from _json_table import json_table
from _library import LIBRARY, Point, number, write

DIRECTORY = "Henriksson_2005_TwoPhotonAbsorption"
SOURCE = "henriksson_2005"
PDF = LIBRARY / DIRECTORY / "original" / "source_01" / "bound-rel.pdf"
# Table, library id, gas, PDF page (0-based) for the text-layer check.
TABLES = (
    ("I", "tbl_930ae06840e1", "He", 4),
    ("II", "tbl_1fcee34b02cc", "Ne", 4),
    ("III", "tbl_37bc1a24b4ae", "Ar", 5),
    ("IV", "tbl_42d892cee858", "Kr", 5),
    ("V", "tbl_e64a945b01a8", "Xe", 6),
)
DESIGNATION = "Desig."  # codespell:ignore desig - the tables' heading, as printed
HEADER = [
    ["Configuration", DESIGNATION, "Expt.^a", "TDHF", "", "J", "TDDHF", ""],
    ["", "", "", "ΔE", "δ^b", "", "ΔE", "δ^b"],
]
NEON_HEADER = [
    ["", "", "", "", "TDHF", "", "", "TDDHF"],
    ["Configuration", DESIGNATION, "Expt.a", "DE", "db", "J", "DE", "db"],
]
NEON: dict[str, tuple[str, tuple[tuple[str, ...], ...]]] = {
    # merged configuration cell: (configuration, columns designation, Expt., TDHF dE,
    # TDHF delta, J, TDDHF dE, TDDHF delta), each a tuple of the rows' cells as the
    # image prints them (the designation in the library's flattened form).
    "2p53s": (
        "2p5 3s",
        (
            ("3Po2", "3Po1", "3Po0", "1Po1"),
            ("16.619", "16.671", "16.716", "16.848"),
            ("17.990", "17.990", "17.990", "18.363"),
            ("0", "0", "0", "0.164 47"),
            ("2", "1", "0", "1"),
            ("17.927", "17.997", "18.045", "18.349"),
            ("0", "0.004 329", "0", "0.160 81"),
        ),
    ),
    "2p53p": (
        "2p5 3p",
        (
            ("3S1", "3D3", "3D2", "3D1", "1D2", "1P1", "3P2", "3P0", "3P1", "1S0"),
            (
                "18.382",
                "18.555",
                "18.576",
                "18.613",
                "18.637",
                "18.694",
                "18.704",
                "18.712",
                "18.727",
                "18.966",
            ),
            (
                "19.532",
                "19.967",
                "19.967",
                "19.967",
                "20.129",
                "20.159",
                "20.159",
                "20.159",
                "20.159",
                "20.564",
            ),
            ("0", "0", "0", "0", "11.060", "0", "0", "0", "0", "12.850"),
            ("1", "3", "2", "1", "2", "1", "2", "0", "1", "0"),
            (
                "19.510",
                "19.910",
                "19.950",
                "19.993",
                "20.090",
                "20.137",
                "20.164",
                "20.173",
                "20.186",
                "20.554",
            ),
            (
                "0",
                "0",
                "0.819 4",
                "0",
                "5.608 6",
                "0",
                "4.712 6",
                "0.3097",
                "0",
                "12.591",
            ),
        ),
    ),
    "2p54s": (
        "2p5 4s",
        (
            ("3Po2", "3Po1", "3Po0", "1Po1"),
            ("19.664", "19.688", "16.761", "19.780"),
            ("21.215", "21.215", "21.215", "21.295"),
            ("0", "0", "0", "0.028 633"),
            ("2", "1", "0", "1"),
            ("21.156", "21.197", "21.276", "21.316"),
            ("0", "0.009 680", "0", "0.018974"),
        ),
    ),
}
"""Table II as the image prints it (2p5 4s 3Po0's 'Expt.' 16.761 is printed so)."""
L_VALUES = {"s": 0, "p": 1, "d": 2, "f": 3}


@dataclass(frozen=True)
class Row:
    """One printed row: configuration, term, and the cells of the columns."""

    configuration: str
    term: str
    odd_term: bool
    cells: tuple[str, ...]
    """Expt., TDHF dE, TDHF delta, J, TDDHF dE, TDDHF delta, as printed."""

    @property
    def even(self) -> bool:
        """Whether the configuration has the ground state's (even) parity."""
        # An occupation is written only before a space or the end ('3p5 4s', '1s2s').
        orbitals = re.findall(r"\d([spdf])(\d+(?=\s|$))?", self.configuration)
        return sum(L_VALUES[ell] * int(n or 1) for ell, n in orbitals) % 2 == 0

    @property
    def printed(self) -> str:
        """The row's cells joined, for the ``printed`` column."""
        return " | ".join((self.configuration, self.term, *self.cells))


def _configuration(cell: str) -> str:
    """``3p^{5}4s`` -> ``3p5 4s``; ``1s2s`` stays."""
    return re.sub(r"\^\{(\d+)\}", r"\1 ", cell)


def _term(designation: str) -> tuple[str, bool]:
    """``^3P_2^o`` or the flattened ``3Po2`` -> (``3P*_2``, True)."""
    match = re.fullmatch(r"\^?(\d)([SPDF])(?:_(\d)(\^o)?|(o?)(\d))", designation)
    if match is None:
        msg = f"unexpected designation {designation!r}"
        raise ValueError(msg)
    multiplicity, letter, j1, odd1, odd2, j2 = match.groups()
    odd = bool(odd1 or odd2)
    return f"{multiplicity}{letter}{'*' if odd else ''}_{j1 or j2}", odd


def _neon(grid: list[list[str]]) -> list[Row]:
    rows = []
    for cells in grid[2:]:
        configuration, columns = NEON[cells[0]]
        for merged, split in zip(cells[1:], columns, strict=True):
            if "".join(split) != merged:
                msg = f"Table II {cells[0]}: {split!r} does not join to {merged!r}"
                raise ValueError(msg)
        for designation, *rest in zip(*columns, strict=True):
            term, odd = _term(designation)
            rows.append(Row(configuration, term, odd, tuple(rest)))
    return rows


def rows(table_id: str) -> list[Row]:
    """Return a table's printed rows."""
    grid = json_table(DIRECTORY, table_id)
    if table_id == TABLES[1][1]:
        if grid[:2] != NEON_HEADER:
            msg = f"unexpected header of Table II: {grid[:2]!r}"
            raise ValueError(msg)
        return _neon(grid)
    if grid[:2] != HEADER:
        msg = f"unexpected header of {table_id}: {grid[:2]!r}"
        raise ValueError(msg)
    out = []
    configuration = ""
    for cells in grid[2:]:
        configuration = _configuration(cells[0]) or configuration
        term, odd = _term(cells[1])
        out.append(Row(configuration, term, odd, tuple(cells[2:])))
    return out


def _value(cell: str) -> str:
    return number(cell.replace(" ", ""))


def extract() -> dict[str, list[Point]]:
    """Return the TDHF and TDDHF one- and two-photon strengths, after the checks."""
    checks()
    files: dict[str, list[Point]] = {}
    for _, table_id, gas, _ in TABLES:
        found: dict[str, list[Point]] = {}
        for row in rows(table_id):
            _, _, tdhf, j, energy, tddhf = row.cells
            kind = "tpa" if row.even else "opa"
            label = f"{row.configuration} {row.term}"
            found.setdefault(f"tdhf_{kind}", []).append(
                Point(label, label, _value(tdhf), "", row.printed)
            )
            label = f"J = {j}, {energy} eV"
            found.setdefault(f"tddhf_{kind}", []).append(
                Point(label, label, _value(tddhf), "", row.printed)
            )
        for stem, points in found.items():
            labels = [p.row for p in points]
            if len(set(labels)) != len(labels):
                msg = f"{gas} {stem}: repeated row labels"
                raise ValueError(msg)
            files[f"{gas}/{SOURCE}__{stem}_{gas}.csv"] = points
    return files


# Hättig et al. (hattig_1998b) as Tables III-V there print them: (level, method) ->
# (CCSD excitation energy in eV from its Table I, S^0f in a.u.).
HATTIG = {
    ("Ne 1De", "CCSD"): (18.532, 12.024),
    ("Ne 1Se", "CCSD"): (18.855, 19.915),
    ("Ne 1Se", "CCS"): (20.586, 12.559),
    ("Ar 1De", "CCSD"): (13.218, 104.59),
    ("Ar 1Se", "CCSD"): (13.487, 197.71),
    ("He 1Se", "CCS"): (21.137, 10.824),
    ("He 1Se", "CCSD"): (20.615, 12.306),
}
# What this paper's text quotes from them (p. 114106-8): level -> (dE, delta).
QUOTED_HATTIG = {
    "Ne 1D2": (18.532, 19.915),
    "Ne 1S0": (18.855, 12.559),
    "Ar 1D2": (13.218, 104.59),
    "Ar 1S0": (13.487, 197.71),
    "He 1S0": (20.615, 12.306),
}


def _find(gas: str, column: int, label: str) -> float:
    """Return a printed number: TDHF rows by 'configuration term', TDDHF by 'J, dE'."""
    table_id = next(t for _, t, g, _ in TABLES if g == gas)
    for row in rows(table_id):
        names = (
            f"{row.configuration} {row.term}",
            f"J = {row.cells[3]}, {row.cells[4]} eV",
        )
        if label in names:
            return float(row.cells[column].replace(" ", ""))
    msg = f"{gas}: no row {label!r}"
    raise KeyError(msg)


def checks() -> dict[str, object]:  # noqa: PLR0912 - one list of checks
    """Return the checks of the tables against the PDF and the text; raise on a miss."""
    out: dict[str, object] = {}
    failures = []
    document = pypdfium2.PdfDocument(PDF)
    pages = {
        page: re.sub(r"\s+", "", document[page].get_textpage().get_text_range())
        for page in {p for *_, p in TABLES}
    }
    for name, table_id, gas, page in TABLES:
        printed = rows(table_id)
        missing = [
            row.printed
            for row in printed
            if "".join(c.replace(" ", "") for c in row.cells) not in pages[page]
        ]
        out[f"Table {name} ({gas}): rows found in the PDF text"] = (
            len(printed) - len(missing),
            len(printed),
        )
        if missing:
            failures.append(f"Table {name}: rows not in the PDF text: {missing}")
        if any(row.odd_term == row.even for row in printed):
            failures.append(
                f"Table {name}: a term's parity mark contradicts its configuration"
            )
        # The TDDHF rows of each configuration's block against its LS terms' J.
        for configuration in dict.fromkeys(row.configuration for row in printed):
            block = [row for row in printed if row.configuration == configuration]
            ls = sorted(row.term[-1] for row in block)
            dhf = sorted(row.cells[3] for row in block)
            if ls != dhf:
                out[
                    f"Table {name} {configuration}: TDDHF J differ from the LS terms'"
                ] = (
                    "".join(ls),
                    "".join(dhf),
                )
                if gas != "Xe":
                    failures.append(f"Table {name} {configuration}: TDDHF J != LS J")
    # p. 114106-8: the tensor elements S_xy of the 1D2 levels, delta = (4/3) S_xy^2
    # (summed over M_J, linear parallel polarisation).
    for gas, column, label, s_xy in (
        ("Ne", 2, "2p5 3p 1D_2", 2.8801),
        ("Ne", 5, "J = 2, 20.090 eV", 2.0510),
        ("Ar", 2, "3p5 4p 1D_2", 9.222),
    ):
        table = _find(gas, column, label)
        out[f"{gas} {label}: (4/3) S_xy^2 vs table"] = (4 / 3 * s_xy**2, table)
        # S_xy to 4-5 figures: 2 dS/S relative, at most 1.1e-4 (Ar's 9.222).
        if abs(4 / 3 * s_xy**2 / table - 1) > 0.0002:
            failures.append(f"{gas} {label}: delta is not (4/3) S_xy^2")
    # p. 114106-4: He 'X 1S0 -> 2 1S0 ... from delta = 1.6002 a.u. (see Table I)'.
    he = _find("He", 5, "J = 0, 23.405 eV")
    out["He 2 1S0 TDDHF (text 1.6002)"] = he
    if he != 1.6002:
        failures.append("He 2 1S0 TDDHF is not the text's 1.6002")
    # p. 114106-4: Ne 'X 1S0 -> 1 3P1 ... 0.0043', '1 1P1 ... delta_OPA = 0.1608',
    # 'nonrelativistic value of delta_OPA = 0.1645'.
    for column, label, text in (
        (5, "J = 1, 17.997 eV", 0.0043),
        (5, "J = 1, 18.349 eV", 0.1608),
        (2, "2p5 3s 1P*_1", 0.1645),
    ):
        table = _find("Ne", column, label)
        out[f"Ne {label} (text {text})"] = table
        if round(table, 4) != text:
            failures.append(f"Ne {label}: {table} is not the text's {text}")
    # p. 114106-5: argon's leakage 'amounts to 9%': 3P1 over 1P1 (TDDHF).
    leak = _find("Ar", 5, "J = 1, 11.864 eV") / _find("Ar", 5, "J = 1, 12.206 eV")
    out["Ar leakage 3P1/1P1 (text 9 %)"] = leak
    if round(leak, 2) != 0.09:
        failures.append(f"Ar leakage {leak:.3f} is not 9 %")
    # pp. 114106-7/8: the integrated TDDHF strengths of each J 'virtually the same' as
    # TDHF's (Ne, Ar; Figs. 2, 3).
    for gas in ("Ne", "Ar"):
        table_id = next(t for _, t, g, _ in TABLES if g == gas)
        block = [r for r in rows(table_id) if r.even]
        for j in ("2", "0"):
            tdhf = sum(float(r.cells[2]) for r in block if r.term.endswith(j))
            tddhf = sum(
                float(r.cells[5].replace(" ", "")) for r in block if r.cells[3] == j
            )
            out[f"{gas} J = {j}: sum TDDHF / sum TDHF"] = tddhf / tdhf
            if abs(tddhf / tdhf - 1) > 0.02:
                failures.append(f"{gas} J = {j}: integrated TDDHF departs from TDHF")
    # Reported: the text's quotation of hattig_1998b.
    for level, (energy, delta) in QUOTED_HATTIG.items():
        matches = [k for k, (e, s) in HATTIG.items() if s == delta]
        out[f"text quotes Hattig {level}: {energy} eV, {delta} a.u."] = matches
    if failures:
        raise ValueError("; ".join(failures))
    return out


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
