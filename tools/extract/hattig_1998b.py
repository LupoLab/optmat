# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hättig et al. (1998): coupled-cluster two-photon strengths of He, Ne, Ar.

J. Chem. Phys. 108, 8355-8359 (1998); journal p. 8355 is PDF page 2 (PDF page 1 is the
AIP download cover), so the journal page is the PDF page plus 8353.

* **Tables III, IV, V** (PDF pp. 5-6, journal pp. 8358-8359; library tables
  ``tbl_1bfdadd98b52``, ``tbl_821a4b276e0b``, ``tbl_75eb5ba0dfd3``): the two-photon
  transition strength ``S^{0f}`` (a.u.) and rate constant ``K_{0-f}`` (cm^4 s) from the
  ground state, CCS, CC2 and CCSD (and, for the 1De level, SCF and MCSCF results of
  Refs. 6 and 12), for ``omega_1 = omega_2 = omega_f0 / 2`` and for ``omega_2 = 10.49
  eV``. Taken: the CCS, CC2 and CCSD strengths of the 1Se levels (He 1s2s, Ne 2p5 3p,
  Ar 3p5 4p), one file per gas, method and pair of frequencies, the row named by the
  method as printed and ``x`` by the level (Table III prints none: its caption's "the
  1Se state in He" is Table I's 1s1 2s1 (1Se)). Table V prints the argon levels as
  "2p5 3p1"; the files keep that, and the record corrects it by errata.
* **Table II** (PDF p. 4, journal p. 8357, ``tbl_ee7b52d10a52``): the CCS, CC2 and
  CCSD oscillator strengths of the ground state to the lowest 1P level, one file per
  gas and method, ``x`` the column heading as printed ("X 1Se - 1Pe"; corrected to
  1Po by errata). The experimental rows (Chan et al.) and Sundholm et al.'s SCF and
  CAS rows are not taken.

Not taken: the 1De rows of Tables IV and V, whose ``S^{0f}`` is not the rotationally
averaged strength (below), the K columns (a derived quantity, reproduced by
:func:`checks`), the experimental K, and Table I (computed excitation energies; used
by :func:`checks`).

The tables have two header rows with merged cells and are read from the library's
cells (``.json``); every cell was checked against the table images (2026-10-08). Cells
keep their LaTeX in ``printed``.

**Checks** (:func:`checks`; the script fails if one is missed):

* every K the tables print (CC and the reproduced SCF/MCSCF rows, 1Se and 1De, both
  pairs of frequencies) from its ``S^{0f}`` by Eq. (15) with the method's own excitation
  energy (Table I; for Ne's SCF and MCSCF footnote b's ``omega_1 = omega_2 = 0.3425``
  a.u.) and ``omega_2 = 0.3855`` a.u., within the rounding of the printed K. Eq. (15) as
  printed, ``K = a0^4 t0 (2 pi)^2 alpha^2 omega_1 omega_2 delta^TP``, gives every K a
  factor 2 pi too small: the tables follow ``(2 pi)^3``. For a 1Se level the tables'
  ``delta^TP`` is ``S^{0f}``; for a 1De level it is 4/15 ``S^{0f}``, the rotational
  average (Eqs. 16-19, F = G = H = 2) of one component of the level whose only tensor
  elements are ``S_xy = S_yx``, with ``S^{0f} = S_xy^2``. The cell printed
  "1.19\times 20^{-51}" (Ar, SCF, 10.49 eV) is read as 1.19e-51 (a misprint).
* the SCF and MCSCF ``S^{0f}`` of the 1De level against the squares of the tensor
  elements ``S_xy`` that henriksson_2005 quotes from Sundholm et al. (1994) and Hettema
  et al. (1992).
* the text: ``omega_2 = 10.49 eV = 0.3855`` a.u. and the CCSD and experimental K of
  argon's 1De level quoted on p. 8358.
"""

import math
import re

from _json_table import json_table
from _library import Point, number, table, write

from optmat.constants import (
    AU_ANGULAR_FREQUENCY,
    AU_ENERGY,
    BOHR_RADIUS,
    E_CHARGE,
    EPS0,
    HBAR,
    C,
)

DIRECTORY = "Hattig_1998_CoupledClusterResponse"
SOURCE = "hattig_1998b"
TABLE_I = "tbl_3100671a92e3"
TABLE_II = "tbl_ee7b52d10a52"
# Table, gas.
STRENGTHS = (
    ("tbl_1bfdadd98b52", "He"),
    ("tbl_821a4b276e0b", "Ne"),
    ("tbl_75eb5ba0dfd3", "Ar"),
)
_S = r"S^{0f}/\text{a.u.}"
_K = r"K_{0-f}/\text{cm}^{4} \text{ s}"
HEADER = [
    [
        "Method",
        r"\omega_{1}=\omega_{2}=\frac{1}{2}\omega_{f0}",
        "",
        r"\omega_{2}=10.49 \text{ eV}",
        "",
    ],
    ["", _S, _K, _S, _K],
]
TAKEN = ("CCS", "CC2", "CCSD")
# The two pairs of frequencies: the table's columns of S and K, and the dataset suffix.
PAIRS = (("degenerate", 1, 2), ("w2_10p49eV", 3, 4))
HE_LEVEL = "1s1 2s1 (1Se)"
"""Table III has no level row: its caption's 'the 1Se state in He', the abstract's 'the
first excited 1Se state in helium', Table I's 1s1 2s1 (1Se)."""
HEADING_II = "fX 1Se -1Pe"
TRANSITION_II = "X 1Se - 1Pe"
GASES_II = {"Helium": "He", "Neon": "Ne", "Argon": "Ar"}

# Atomic units (CODATA, via optmat.constants): a0 in cm, t0 = hbar / E_h in s.
A0_CM = BOHR_RADIUS * 100.0
T0 = 1.0 / AU_ANGULAR_FREQUENCY
ALPHA = E_CHARGE**2 / (4 * math.pi * EPS0 * HBAR * C)
HARTREE_EV = AU_ENERGY / E_CHARGE
EQ15 = A0_CM**4 * T0 * (2 * math.pi) ** 2 * ALPHA**2
"""Eq. (15)'s factor as printed, cm^4 s per a.u. of omega_1 omega_2 delta^TP."""
OMEGA2 = 0.3855
"""omega_2 in a.u., as the text gives it (p. 8358: '10.49 eV = 0.3855 a.u.')."""
NE_REF6_OMEGA = 0.3425
"""Table IV, footnote b: K of Ref. 6's SCF and MCSCF 'assuming omega_1 = omega_2 =
0.3425 a.u.'."""
D_FACTOR = 4 / 15
"""delta^TP / S^{0f} of a 1De row: one component, S_xy = S_yx only, F = G = H = 2."""
MISPRINT = r"1.19\times 20^{-51}"
"""Table V, Ar 1De SCF, omega_2 = 10.49 eV: '20' for '10' (the PDF's text and image)."""
TABLE_I_LEVEL = {
    ("He", "1Se"): "1s12s1 (1Se)",
    ("Ne", "1Se"): "2p53p1 (1Se)",
    ("Ne", "1De"): "2p53p1 (1De)",
    ("Ar", "1Se"): "3p54p1 (1Se)",
    ("Ar", "1De"): "3p54p1 (1De)",
}
# Tensor elements S_xy (a.u.) of the 1De level that henriksson_2005 quotes (J. Chem.
# Phys. 122, 114106, pp. 114106-7 and -8) from Hettema et al. (1992; Ref. 6 here) and
# Sundholm et al. (1994; Ref. 12 here), against this paper's S^{0f} of the same rows.
QUOTED_SXY = {
    ("Ne", "SCF"): 2.9663,
    ("Ar", "SCF"): 9.499,
    ("Ar", "MCSCF"): 10.424,
}
QUOTED_SXY_UNMATCHED = {("Ne", "MCSCF"): 3.3450}
"""Henriksson's 'SD RAS' of Hettema et al. (19.170 eV); Table I's MCSCF of Ref. 6 is at
18.874 eV, a different calculation: reported, not checked."""


def _label(latex: str) -> str:
    r"""Return a level label from a table's LaTeX.

    ``2p^{5}3p^{1} (^{1}S^{e})`` becomes ``2p5 3p1 (1Se)``.
    """
    text = re.sub(r"\^\{(\w+)\}(?=\d)", r"\1 ", latex)
    return re.sub(r"[\^_{}]", "", text)


def _cell_number(cell: str) -> str:
    """Return the number of a table cell, without its footnote mark."""
    return number(re.sub(r"\^\{\\text\{[ab]\}\}$", "", cell))


def _k_value(cell: str) -> tuple[float, float]:
    r"""Return a K cell's value in cm^4 s and its relative rounding.

    The cell is ``0.4091 \times 10^{-51}``, with or without a footnote mark; the
    rounding is half a unit of its last printed digit.
    """
    match = re.fullmatch(
        r"([\d.]+)\s*\\times\s*(?:10|(20))\^\{-51(?:\\text\{[ab]\})?\}", cell
    )
    if match is None or (match.group(2) is not None and cell != MISPRINT):
        msg = f"unexpected K cell {cell!r}"
        raise ValueError(msg)
    mantissa = match.group(1)
    decimals = len(mantissa.split(".")[1])
    return float(mantissa) * 1e-51, 0.5 * 10**-decimals / float(mantissa)


def _sections(table_id: str, gas: str) -> list[tuple[str, str, list[str]]]:
    """Return a strength table's method rows as (term, level label, cells)."""
    grid = json_table(DIRECTORY, table_id)
    # Table V spaces its unit slashes ("S^{0f} / a.u."); the others do not.
    header = [[re.sub(r"\s*/\s*", "/", cell) for cell in row] for row in grid[:2]]
    if header != HEADER:
        msg = f"unexpected header of {table_id}: {grid[:2]!r}"
        raise ValueError(msg)
    rows = []
    term, level = ("1Se", HE_LEVEL) if gas == "He" else ("", "")
    for cells in grid[2:]:
        if cells[0].startswith("2p^{5}3p^{1}"):
            level = _label(cells[0])
            term = level[-4:-1]
            continue
        rows.append((term, level, cells))
    return rows


def _excitation_energies() -> dict[tuple[str, str], dict[str, float]]:
    """Return Table I's excitation energies (eV) by (gas, term) and method."""
    rows = table(DIRECTORY, TABLE_I)
    methods = rows[0][1:]
    if methods != ["SCF", "CCS", "CC2", "CCSD", "MCSCF", "Exp."]:
        msg = f"unexpected header of Table I: {rows[0]!r}"
        raise ValueError(msg)
    by_label = {cells[0]: cells[1:] for cells in rows[1:]}
    energies = {}
    for key, label in TABLE_I_LEVEL.items():
        energies[key] = {
            method: float(cell.split()[0])
            for method, cell in zip(methods, by_label[label], strict=True)
            if cell
        }
    return energies


def _strengths() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for table_id, gas in STRENGTHS:
        for term, level, cells in _sections(table_id, gas):
            method = cells[0]
            if term != "1Se" or method not in TAKEN:
                continue
            printed = " | ".join(cells)
            for suffix, s_column, _ in PAIRS:
                value = _cell_number(cells[s_column])
                name = f"{gas}/{SOURCE}__S_{gas}_{method}_{suffix}.csv"
                files[name] = [Point(method, level, value, "", printed)]
    return files


def _oscillator_strengths() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_II)
    if rows[0] != ["Method", HEADING_II]:
        msg = f"unexpected header of Table II: {rows[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in rows[1:]:
        if cells[0] in GASES_II:
            gas = GASES_II[cells[0]]
            continue
        if cells[0] not in TAKEN:
            continue  # Exp. (Chan et al.), SCF and CAS (Ref. 12)
        printed = " | ".join(cells)
        name = f"{gas}/{SOURCE}__f_{gas}_{cells[0]}.csv"
        files[name] = [Point(cells[0], TRANSITION_II, number(cells[1]), "", printed)]
    return files


def _k_rows() -> list[tuple[str, str, str, str, float, float, tuple[float, float]]]:
    """Return every printed K: (gas, term, method, pair, S, omega_1 omega_2, K).

    K is the value in cm^4 s and its relative rounding.
    """
    energies = _excitation_energies()
    out = []
    for table_id, gas in STRENGTHS:
        for term, _, cells in _sections(table_id, gas):
            method = cells[0]
            if method == "Exp.":
                continue
            for suffix, s_column, k_column in PAIRS:
                if not cells[k_column]:
                    continue
                if gas == "Ne" and method in {"SCF", "MCSCF"}:
                    omega_f = 2 * NE_REF6_OMEGA
                else:
                    omega_f = energies[gas, term][method] / HARTREE_EV
                if suffix == "degenerate":
                    product = (omega_f / 2) ** 2
                else:
                    product = (omega_f - OMEGA2) * OMEGA2
                strength = float(_cell_number(cells[s_column]))
                k = _k_value(cells[k_column])
                out.append((gas, term, method, suffix, strength, product, k))
    return out


def _check_k(out: dict[str, object], failures: list[str]) -> None:
    """Reproduce every printed K from its S^{0f} (2 pi times Eq. 15)."""
    ratios = []
    for gas, term, method, suffix, strength, product, (k, rounding) in _k_rows():
        factor = 1.0 if term == "1Se" else D_FACTOR
        reproduced = 2 * math.pi * EQ15 * product * factor * strength
        deviation = reproduced / k - 1
        # The printed K's last digit, and S's and Table I's (0.0005 eV) rounding.
        margin = rounding + 0.001
        out[f"K {gas} {term} {method} {suffix}"] = (reproduced, k, deviation)
        ratios.append(k / (EQ15 * product * factor * strength))
        if abs(deviation) > margin:
            failures.append(f"K {gas} {term} {method} {suffix}: {deviation:+.4f}")
    out["K over Eq. (15) as printed (min, max), 2 pi = 6.2832"] = (
        min(ratios),
        max(ratios),
    )
    # 0.3 %: the largest rounding of a printed K (three figures) with Table I's.
    if not all(abs(r / (2 * math.pi) - 1) < 0.003 for r in ratios):
        failures.append("a K departs from 2 pi times Eq. (15)")


def _check_sxy(out: dict[str, object], failures: list[str]) -> None:
    """Check the 1De S^{0f} of the reproduced rows against the quoted S_xy squared."""
    for table_id, gas in STRENGTHS[1:]:
        for term, _, cells in _sections(table_id, gas):
            key = (gas, cells[0])
            if term != "1De" or key not in QUOTED_SXY | QUOTED_SXY_UNMATCHED:
                continue
            printed = float(_cell_number(cells[1]))
            if key in QUOTED_SXY:
                square = QUOTED_SXY[key] ** 2
                out[f"S {gas} 1De {cells[0]} = S_xy^2"] = (square, printed)
                # S_xy to 4-5 figures: 2 dS/S is at most 1.1e-4 (9.499).
                if abs(square / printed - 1) > 0.0005:
                    failures.append(f"S {gas} 1De {cells[0]} is not S_xy^2")
            else:
                square = QUOTED_SXY_UNMATCHED[key] ** 2
                out[f"S {gas} 1De {cells[0]} vs SD RAS S_xy^2 (unmatched)"] = (
                    square,
                    printed,
                )


def _check_text(out: dict[str, object], failures: list[str]) -> None:
    """Check the numbers the text prints (p. 8358)."""
    # '10.49 eV = 0.3855 a.u.'; 'experimental value of 2.05 x 10^-51 cm4 s agrees well
    # with the CCSD result of 1.83 x 10^-51 cm4 s'.
    out["10.49 eV in a.u."] = 10.49 / HARTREE_EV
    if abs(10.49 / HARTREE_EV - OMEGA2) > 0.00005:
        failures.append("10.49 eV is not 0.3855 a.u.")
    argon = {
        cells[0]: cells[4]
        for term, _, cells in _sections("tbl_75eb5ba0dfd3", "Ar")
        if term == "1De" or cells[0] == "Exp."
    }
    out["Ar 1De K CCSD and Exp. (text)"] = (argon["CCSD"], argon["Exp."])
    if argon["CCSD"] != r"1.83\times 10^{-51}" or not argon["Exp."].startswith(
        r"2.05\times 10^{-51"
    ):
        failures.append("Table V's argon K differ from the text's 1.83 and 2.05e-51")
    # 'For all three atoms we obtain at the CCSD level oscillator strengths ... within
    # the experimental uncertainities' (Table II's Exp. rows, whose '+-' the library's
    # cells render as '6': '0.28060.007 a' is 0.280 +- 0.007).
    gas = ""
    ccsd: dict[str, float] = {}
    for cells in table(DIRECTORY, TABLE_II)[1:]:
        if cells[0] in GASES_II:
            gas = GASES_II[cells[0]]
        elif cells[0] == "CCSD":
            ccsd[gas] = float(cells[1])
        elif cells[0] == "Exp.":
            match = re.fullmatch(r"(\d\.\d{3})6(\d\.\d{3}) [a-d]", cells[1])
            if match is None:
                failures.append(f"unexpected Exp. cell {cells[1]!r} of Table II")
                continue
            f_exp, sigma = float(match.group(1)), float(match.group(2))
            out[f"f {gas} CCSD vs Exp."] = (ccsd[gas], f_exp, sigma)
            if abs(ccsd[gas] - f_exp) > sigma:
                failures.append(f"f {gas}: CCSD outside the experimental uncertainty")


def checks() -> dict[str, object]:
    """Return the checks of the K, the 1De convention and the text; raise on a miss."""
    out: dict[str, object] = {}
    failures: list[str] = []
    _check_k(out, failures)
    _check_sxy(out, failures)
    _check_text(out, failures)
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return the 1Se strengths and the oscillator strengths, after the checks."""
    checks()
    return _strengths() | _oscillator_strengths()


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
