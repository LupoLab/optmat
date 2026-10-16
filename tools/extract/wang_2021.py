# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Wang et al. (2021): polarizabilities of He, Ne, Ar and their cations; I_p.

Table 2 (``tbl_d0fc08f74765``; image and page checked) lists, per atom, the paper's
HF and CCSD(T) values, numerical HF (NHF), other calculations ("Calc.") and
experimental values ("Expt."), in atomic units; the Calc. and Expt. cells hold several
values, each with a footnote letter naming its reference (the library joins them with
``\``). Only the curated species (He, Ne, Ar) are extracted: the paper's CCSD(T) value
(column 3), its HF value (column 1), the NHF value (column 2) and each quoted value as
its own dataset, chosen by its footnote letter. The ECC column (CCSD(T) - HF) is not
extracted.

Table 1 (``tbl_c22e704aed88``, CSV; image checked, p. 4): the first ionisation
potentials (eV), CCSD(T)/aug-PC-4 ("Calc.") and experimental ("Exp.^a", Moore 1970),
in three blocks of seven atoms. The library lost the header rows of the second and third
blocks into the experimental cells of the block above ("24.587F": He's 24.587 and F,
the second block's second atom); the script reads each block's atoms from those
suffixes. He, Ne and Ar are taken.

Table 3 (three library tables, image checked): the polarizabilities of the free ions,
rows HF, "Ours." (CCSD(T), per the caption), ECC, CHF, UHF, Calc. and Expt. He+ is in
the first part (``tbl_82ace8ce2dca``, p. 5), Ne+ in the continuation
(``tbl_eede931cbc2d``, p. 6), Ar+ in the unlabelled last part (``tbl_0bf0dd86ab78``,
p. 7, cells only). The first part's cells break numbers with spaces ("0. 281 1" for the
printed 0.2811) and run footnote letters on without a caret ("0. 281 19e"); the spaces
are removed and the letter read as the footnote. In the last part the library reads the
footnote of Ar+'s Calc. value as "ijj"; the image prints "jjj" (Ref. [96]).
"""

import re

from _json_table import json_table
from _library import LIBRARY, Point, number, table, write

DIRECTORY = "Wang_2021_StaticDipolePolarizabilities"
TABLE = "tbl_d0fc08f74765"
GASES = ("He", "Ne", "Ar")
# (gas, footnote letter) -> dataset, for the Calc. (column 4) and Expt. (column 6)
# values.
QUOTED = {
    ("He", "d"): "yan_pachucki_He",
    ("He", "e"): "lach_2004_He",
    ("He", "f"): "expt_newell_grohmann_He",
    ("Ne", "s"): "calc_82_84_Ne",
    ("Ne", "t"): "klopper_2004_Ne",
    ("Ne", "u"): "dalgarno_kingston_1960_Ne",
    ("Ne", "v"): "gaiser_fellmuth_2010_Ne",
    ("Ar", "ee"): "reinsch_meyer_1976_Ar",
    ("Ar", "jj"): "lupinetti_hohm_Ar",
    ("Ar", "kk"): "newell_baird_1965_Ar",
    ("Ar", "ll"): "johnston_1966_Ar",
}
_ENTRY = re.compile(r"^([\d.]+)(?:\s*±\s*([\d.]+))?\^\{?([a-z]+)\}?$")


def extract() -> dict[str, list[Point]]:
    """Return Table 2 for He, Ne and Ar, Table 3 for their cations, Table 1's I_p."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        gas = cells[0].strip()
        if gas not in GASES:
            continue
        printed = " | ".join(cells)
        files[f"{gas}/wang_2021__ccsdt_{gas}.csv"] = [
            Point("1", "inf", number(cells[3]), "", printed)
        ]
        for column in (4, 6):
            for entry in filter(None, (e.strip() for e in cells[column].split("\\"))):
                match = _ENTRY.match(entry)
                if match is None:
                    msg = f"cannot read {entry!r} in the {gas} row"
                    raise ValueError(msg)
                value, sigma, mark = match.groups()
                dataset = QUOTED[gas, mark]
                files[f"{gas}/wang_2021__{dataset}.csv"] = [
                    Point("1", "inf", number(value), sigma or "", printed)
                ]
    return files | _table_2_hf() | _table_1() | _table_3()


# Table 2 columns 1 (HF) and 2 (NHF, Stiehler and Hinze, Ref. [64]).
TABLE_2_COLUMNS = {1: "hf", 2: "nhf"}
TABLE_1 = "tbl_c22e704aed88"
_EXPT_1 = re.compile(r"^([\d.]+)([A-Z][a-z]?)?$")
# Table 3: the part holding each cation, its column, and the quoted values by
# footnote letter.
TABLE_3 = {
    "He+": ("tbl_82ace8ce2dca", "He+"),
    "Ne+": ("tbl_eede931cbc2d", "Ne^+"),
    "Ar+": ("tbl_0bf0dd86ab78", "Ar^+"),
}
TABLE_3_QUOTED = {
    ("He+", "e"): "calc_tang_2012_He+",
    ("Ne+", "dd"): "calc_medved_2000_Ne+",
    ("Ne+", "gg"): "calc_thakkar_2001_Ne+",
    ("Ar+", "rrr"): "chf_stewart_1975_Ar+",
    ("Ar+", "jjj"): "calc_medved_2000_Ar+",
    ("Ar+", "qqq"): "expt_bristow_1972_Ar+",
}
# The Ar+ Calc. footnote: the cells file reads "ijj", the image "jjj".
FOOTNOTE_FIXES = {("Ar+", "ijj"): "jjj"}
_QUOTED_3 = re.compile(r"([\d.]+)(?:\((\d+)\))?(?:\^\{([a-z]+)\}|\^([a-z]+)|([a-z]+))")


def _table_2_hf() -> dict[str, list[Point]]:
    """Return the HF and NHF values of He, Ne and Ar (Table 2)."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        gas = cells[0].strip()
        if gas not in GASES:
            continue
        printed = " | ".join(cells)
        for column, name in TABLE_2_COLUMNS.items():
            files[f"{gas}/wang_2021__{name}_{gas}.csv"] = [
                Point("1", "inf", number(cells[column]), "", printed)
            ]
    return files


def _table_1() -> dict[str, list[Point]]:
    """Return the CCSD(T) and experimental first ionisation potentials."""
    files: dict[str, list[Point]] = {}
    grid = table(DIRECTORY, TABLE_1)
    atoms = [cell.strip() for cell in grid[0][1:]]
    for calc, expt in zip(grid[1::2], grid[2::2], strict=True):
        if calc[0] != "Calc." or not expt[0].startswith("Exp"):
            msg = f"unexpected Table 1 rows {calc!r}, {expt!r}"
            raise ValueError(msg)
        following = []
        for column, atom in enumerate(atoms, start=1):
            cell = expt[column].strip()
            match = _EXPT_1.match(cell) if cell else None
            following.append(match.group(2) if match and match.group(2) else "")
            if atom not in GASES:
                continue
            if match is None:
                msg = f"unreadable Table 1 cell {expt[column]!r}"
                raise ValueError(msg)
            printed = f"{atom} | Calc. {calc[column]} | Exp. {expt[column]}"
            files[f"{atom}/wang_2021__ip_ccsdt_{atom}.csv"] = [
                Point("I_p", "", number(calc[column]), "", printed)
            ]
            files[f"{atom}/wang_2021__ip_expt_{atom}.csv"] = [
                Point("I_p", "", number(match.group(1)), "", printed)
            ]
        atoms = following
    if len(files) != 2 * len(GASES):
        msg = f"Table 1 lost an atom: {sorted(files)}"
        raise ValueError(msg)
    return files


def _table_3_part(table_id: str) -> list[list[str]]:
    """Return a part of Table 3, from its CSV if it has one, else its cells."""
    if (LIBRARY / DIRECTORY / "tables" / f"{table_id}.csv").exists():
        return table(DIRECTORY, table_id)
    return json_table(DIRECTORY, table_id)


def _table_3() -> dict[str, list[Point]]:
    """Return the HF, CCSD(T) and quoted polarizabilities of He+, Ne+ and Ar+."""
    files: dict[str, list[Point]] = {}
    for ion, (table_id, heading) in TABLE_3.items():
        grid = _table_3_part(table_id)
        block = next(
            i for i, cells in enumerate(grid) if heading in (c.strip() for c in cells)
        )
        column = [c.strip() for c in grid[block]].index(heading)
        rows: dict[str, str] = {}
        for cells in grid[block + 1 :]:
            name = cells[0].strip().rstrip(".").split("^")[0].rstrip("d").strip()
            if not name or name in rows:
                break
            rows[name] = cells[column].replace(" ", "")
        printed = f"{heading} | " + " | ".join(f"{k}: {v}" for k, v in rows.items())
        files[f"{ion}/wang_2021__hf_{ion}.csv"] = [
            Point("1", "inf", number(rows["HF"]), "", printed)
        ]
        files[f"{ion}/wang_2021__ccsdt_{ion}.csv"] = [
            Point("1", "inf", number(rows["Ours"]), "", printed)
        ]
        for row in ("CHF", "Calc", "Expt"):
            for match in _QUOTED_3.finditer(rows.get(row, "")):
                value, sigma = match.group(1), match.group(2) or ""
                mark = match.group(3) or match.group(4) or match.group(5)
                mark = FOOTNOTE_FIXES.get((ion, mark), mark)
                dataset = TABLE_3_QUOTED[ion, mark]
                if sigma:  # uncertainty in the last digits, in parentheses
                    decimals = len(value.split(".")[1]) if "." in value else 0
                    sigma = f"{int(sigma) / 10**decimals:.{decimals}f}"
                files[f"{ion}/wang_2021__{dataset}.csv"] = [
                    Point("1", "inf", number(value), sigma, printed)
                ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
