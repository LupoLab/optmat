# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract van Duijnen & Swart (1998): the "experimental" polarizabilities it uses.

Two tables (JSON cells only; images checked, cells agree):

- Table 1 (``tbl_3b0fee402144``): the experimental mean polarizabilities (au) of
  Thole's original learning set (16 molecules) and control set (6): of the curated
  species carbon monoxide, hydrogen, nitrogen, oxygen and methane.
- Table 11 (``tbl_737b997b82b0``) and its continuation (``tbl_f45b9a639abe``): the 52 +
  18 molecules of the extended fit, "exp" column: of the curated species also CO2, SF6
  and N2O. Its CO, H2, N2 and O2 equal Table 1's (asserted; not written twice); its
  methane is printed 17.86 against Table 1's 17.68 and goes to its own file.

Only the "exp" column is taken; the model and CPHF columns are the authors' fits and
calculations of Thole-model parameters, not polarizability data.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "VanDuijnen_1998_MolecularAtomicPolarizabilities"
TABLE_1 = "tbl_3b0fee402144"
TABLE_11 = ("tbl_737b997b82b0", "tbl_f45b9a639abe")
WAVELENGTH = "5893"  # angstrom, "usually at 5893 A, sodium D line" (p. 4)
TABLE_1_GASES = {
    "carbon monoxide": "CO",
    "hydrogen": "H2",
    "nitrogen": "N2",
    "oxygen": "O2",
    "methane": "CH4",
}
TABLE_11_GASES = {
    "carbon monoxide": "CO",
    "hydrogen": "H2",
    "nitrogen": "N2",
    "oxygen": "O2",
    "methane": "CH4",
    "carbon dioxide": "CO2",
    "sulfur hexafluoride": "SF6",
    "nitrous oxide": "N2O",
}


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and table."""
    files: dict[str, list[Point]] = {}
    first: dict[str, str] = {}
    rows = json_table(DIRECTORY, TABLE_1)
    assert rows[0][:2] == ["molecule", "exp"], rows[0]
    for cells in rows:
        gas = TABLE_1_GASES.get(cells[0].strip())
        if gas is None:
            continue
        first[gas] = number(cells[1])
        files[f"{gas}/van_duijnen_1998__table1_{gas}.csv"] = [
            Point(WAVELENGTH, WAVELENGTH, first[gas], "", " | ".join(cells))
        ]
    for table_id in TABLE_11:
        rows = json_table(DIRECTORY, table_id)
        assert rows[0][:3] == ["model", "", "exp"], rows[0]
        for cells in rows:
            gas = TABLE_11_GASES.get(cells[0].strip())
            if gas is None:
                continue
            value = number(cells[2])
            if gas in first and value == first[gas]:
                continue  # the same value as Table 1
            files[f"{gas}/van_duijnen_1998__table11_{gas}.csv"] = [
                Point(WAVELENGTH, WAVELENGTH, value, "", " | ".join(cells))
            ]
    assert len(first) == len(TABLE_1_GASES), first
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
