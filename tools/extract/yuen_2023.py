# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Yuen, Modak, Song, Zhao & Lin (2023): O2 in the DM-SDI model (TC-956).

Read from the library's CSV renderings; every cell was checked against the table
images (2026-10-09).

* Table I (p. 2, ``tbl_8972394235c5``): the experimental vertical I_p of O2 to the
  O2+ X 2Pi_g, a 4Pi_u and b 4Sigma_g- states (Baltzer et al. 1992), and the MO-ADK
  C_lm of the O2 orbital each one removes, as sets "{l, ...}" and "{C, ...}". The
  CSV loses the Greek letters of the state labels ("X2\x05g", "a4\x05u", "b4-g-") and
  writes the configurations in Unicode ("1πg-1"); the image shows X^2Pi_g
  (1pi_g^-1), a^4Pi_u (1pi_u^-1) and b^4Sigma_g^- (3sigma_g^-1). The script names the
  rows from the image.
* Table II (p. 2, ``tbl_4676ff61cd02``): the calculated vertical I_p of O2 to eight
  O2^2+ states from the zero-point level of O2, with each state's limit, KER and
  parents. The KER, one file per limit (O+(4S) + O+(2D): states 1, 2, 4, 5; O+(4S) +
  O+(4S): states 3, 6, 7); state 0 is metastable (KER "-"). The script checks that
  I_p - KER is one number per limit, as the caption says the KER was obtained: 35.68 eV
  (4S + 2D) and 32.35 eV (4S + 4S).
* Table III (p. 3, ``tbl_4b5f5b8c6a36``): the effective binding energies of the
  spin-up orbitals of O2+ X (a model input, ``ionization_model_parameter``) and their
  C_lm.

* Table IV (p. 4, ``tbl_c2f42888899a``): the transition dipoles (au) between O2^2+
  states, every one printed along x-hat. The CSV loses the Greek letters of the state
  labels (Delta as a space, Pi as the control character U+0005, Sigma dropped); the
  image shows W 3Delta_u - B 3Pi_g, B 3Sigma_u- - B 3Pi_g, A 3Sigma_u+ - B 3Pi_g and
  1 1Sigma_u- - 1 1Pi_g, from which the script names the rows. One file, material
  O2++ (registered in TC-954).
"""

from _library import Point, number, table, write

DIRECTORY = "Yuen_2023_ModelingSequentialDissociative"
PI = "\N{GREEK SMALL LETTER PI}"
SIGMA = "\N{GREEK SMALL LETTER SIGMA}"
TABLE_1 = "tbl_8972394235c5"
TABLE_2 = "tbl_4676ff61cd02"
TABLE_3 = "tbl_4b5f5b8c6a36"
TABLE_4 = "tbl_c2f42888899a"
# Table IV: the printed pair of state cells -> (row, the pair as imaged).
PAIRS_4 = {
    ("W3 u", "B3\x05g"): ("W-B", "W 3Delta_u - B 3Pi_g"),
    ("B3-u-", "B3\x05g"): ("B'-B", "B 3Sigma_u- - B 3Pi_g"),
    ("A3-u+", "B3\x05g"): ("A-B", "A 3Sigma_u+ - B 3Pi_g"),
    ("11-u-", "11\x05g"): ("1-1", "1 1Sigma_u- - 1 1Pi_g"),
}
X_HAT = "x\N{MODIFIER LETTER CIRCUMFLEX ACCENT}"
# Table I: the printed state and configuration cells -> (row, state as imaged, name).
STATES_1 = {
    # The lost capital Pi is a control character, U+0005, in the CSV.
    ("X2\x05g", f"1{PI}g-1"): ("X", "X 2Pi_g (1pi_g^-1)", "x"),
    ("a4\x05u", f"1{PI}u-1"): ("a", "a 4Pi_u (1pi_u^-1)", "a"),
    ("b4-g-", f"3{SIGMA}g-1"): ("b", "b 4Sigma_g- (3sigma_g^-1)", "b"),
}
# Table III: the printed orbital -> file name.
ORBITALS_3 = {f"1{PI}g": "1pig", f"1{PI}u": "1piu", f"3{SIGMA}g": "3sigmag"}
LIMITS = {"O^+(^4S)+O^+(^2D)": "4s2d", "O^+(^4S)+O^+(^4S)": "4s4s"}
# I_p - KER, eV: the limit above the zero-point level of O2 implied by the table.
IMPLIED = {"4s2d": 35.68, "4s4s": 32.35}
ROUNDING = 0.005  # half the last printed digit of I_p and KER


def _sets(waves: str, values: str, printed: str) -> list[Point]:
    """Return one point per partial wave from the cells "{l, ...}" and "{C, ...}"."""
    ls = [cell.strip() for cell in waves.strip("{}").split(",")]
    cs = [cell.strip() for cell in values.strip("{}").split(",")]
    if len(ls) != len(cs):
        msg = f"unequal sets of l and C_lm: {printed!r}"
        raise ValueError(msg)
    return [
        Point(str(int(l)), str(int(l)), number(c), "", printed)
        for l, c in zip(ls, cs, strict=True)  # noqa: E741 - l is the partial wave
    ]


def _table_1() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_1)
    if rows[0] != ["States", "Config.", "Ip (eV)", "l", "Clm"]:
        msg = f"unexpected header of {TABLE_1}: {rows[0]!r}"
        raise ValueError(msg)
    energies: list[Point] = []
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        row, state, name = STATES_1[cells[0], cells[1]]
        printed = " | ".join(cells)
        energies.append(Point(row, state, number(cells[2]), "", printed))
        files[f"O2/yuen_2023__table1_cl_{name}_O2.csv"] = _sets(
            cells[3], cells[4], printed
        )
    if [p.row for p in energies] != ["X", "a", "b"]:
        msg = f"Table I rows {[p.row for p in energies]!r}"
        raise ValueError(msg)
    return {"O2/yuen_2023__table1_ip_O2.csv": energies} | files


def _table_2() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_2)
    header = ["Index", "State", "Config.", "I_p (eV)", "Limit", "KER (eV)", "Parents"]
    if rows[0] != header:
        msg = f"unexpected header of {TABLE_2}: {rows[0]!r}"
        raise ValueError(msg)
    energies: list[Point] = []
    ker: dict[str, list[Point]] = {"4s2d": [], "4s4s": []}
    for cells in rows[1:]:
        index, state, _, energy, limit, release, _ = cells
        printed = " | ".join(cells)
        energies.append(Point(index, state, number(energy), "", printed))
        if limit == "Metastable":
            if release != "-":
                msg = f"a metastable state with a KER: {cells!r}"
                raise ValueError(msg)
            continue
        name = LIMITS[limit]
        if abs(float(energy) - float(release) - IMPLIED[name]) > 2 * ROUNDING:
            msg = f"I_p - KER is not the limit's {IMPLIED[name]} eV: {cells!r}"
            raise ValueError(msg)
        ker[name].append(Point(index, state, number(release), "", printed))
    if [p.row for p in energies] != [str(i) for i in range(8)]:
        msg = f"Table II rows {[p.row for p in energies]!r}"
        raise ValueError(msg)
    return {
        "O2/yuen_2023__table2_ip_O2.csv": energies,
        "O2/yuen_2023__table2_ker_4s2d_O2.csv": ker["4s2d"],
        "O2/yuen_2023__table2_ker_4s4s_O2.csv": ker["4s4s"],
    }


def _table_3() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_3)
    if rows[0] != ["Orbitals", "Ip (eV)", "l", "Clm"]:
        msg = f"unexpected header of {TABLE_3}: {rows[0]!r}"
        raise ValueError(msg)
    energies: list[Point] = []
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        name = ORBITALS_3[cells[0]]
        printed = " | ".join(cells)
        energies.append(Point(name, name, number(cells[1]), "", printed))
        files[f"O2+/yuen_2023__table3_cl_{name}_O2+.csv"] = _sets(
            cells[2], cells[3], printed
        )
    if [p.row for p in energies] != ["1pig", "1piu", "3sigmag"]:
        msg = f"Table III rows {[p.row for p in energies]!r}"
        raise ValueError(msg)
    return {"O2+/yuen_2023__table3_eb_O2+.csv": energies} | files


def _table_4() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_4)
    if rows[0] != ["State 1", "State 2", "Dipole moment (a.u.)"]:
        msg = f"unexpected header of {TABLE_4}: {rows[0]!r}"
        raise ValueError(msg)
    points = []
    for cells in rows[1:]:
        row, pair = PAIRS_4[cells[0], cells[1]]
        if not cells[2].endswith(X_HAT):
            msg = f"a dipole not along x: {cells!r}"
            raise ValueError(msg)
        value = number(cells[2].removesuffix(X_HAT))
        points.append(Point(row, pair, value, "", " | ".join(cells)))
    if len(points) != len(PAIRS_4):
        msg = f"Table IV gave {len(points)} rows"
        raise ValueError(msg)
    return {"O2++/yuen_2023__table4_dipole_O2++.csv": points}


def extract() -> dict[str, list[Point]]:
    """Return the datasets of Tables I-IV."""
    return _table_1() | _table_2() | _table_3() | _table_4()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
