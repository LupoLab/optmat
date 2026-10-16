# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Adair, Chase & Payne (1989): n2, linear index, hyperpolarizabilities.

Adair et al. measured n2 of some fifty crystals by nearly degenerate three-wave mixing
at 1.064/1.071 um. Their n2 are printed in 10^-13 esu, n = n0 + n2 <E^2> (Eq. 1), the
quantity ``nonlinear_index_field`` (nothing is converted; the record's scale is 1e-13).

- Table IV (two parts, ``tbl_55822e89640f`` and ``tbl_2e1f875bba6b``, both JSON-only):
  per row, the linear index n at 1.06 um (reproduced from the handbooks of footnotes
  a-k), the measured n2 and the n2 calculated with the BGO formula (Eq. 6, K = 68).
  An index cell carries its footnote letter, as ``^a`` in the first part and as a bare
  trailing letter in the continuation; CeF3's index is printed ``~1.6`` (approximate;
  the table image shows the tilde). Blank calculated cells (no Abbe number) are skipped.
- Tables I (``tbl_c3b6ccc6cfdc``, JSON-only) and II (``tbl_0f4b0f5d826e``): SF-58's n2
  against the frequency offset Delta and against the laser power, both at 1.064 um
  (the offset and the power are per-row columns). CS2 (Table I) is a liquid: skipped.
- Table III (``tbl_fbb684fc9121``, JSON-only): the first column ("TWM^b", this work)
  repeats Table IV, which the script checks, except the unexplained "(av)" of Al2O3
  (1.25, not Table IV's 1.23 or 1.30), which is extracted on its own; the "(av)" of
  LaF3, CeF3 and YAG equals Table IV's single value. The other columns are literature
  values, one file per material and column; a TWM column's second wavelength is the
  per-row column ``omega2``, from the header.
- Table VIII (``tbl_8652b8e99339``): n2 of the cubic crystals, polarization along [100]
  (footnote a). It is extracted as printed, since three cells differ from Table IV
  (KBr 3.07 vs 2.93, CaO 5.07 vs 5.20, MgO 1.6 vs 1.61; checked on the images).
- Table X (``tbl_5192778e0ba3``): the n2 calculated from the simple oxides' pair
  hyperpolarizabilities; its measured column repeats Table IV (checked by the script).
  Table X prints "Gd3Sc2Al2O12" (image) for Table IV's Gd3Sc2Al3O12 (GSAG): a misprint
  of the formula, mapped to GSAG here; its measured n2 equals Table IV's GSAG row.
- Table V (``tbl_62267fd0462a``, JSON-only): the measured and the (adjusted) theoretical
  hyperpolarizability of the seven alkali halides measured here.
- Tables VII (``tbl_00fcb404f207``) and IX (``tbl_152efc508e07``): the measured
  hyperpolarizability per neutral formula unit of the fluorides, oxides and sulfides,
  with the anion-cation separation as an extra column. Table VII repeats Table V's
  LiF, NaF and KF values; those rows are checked equal and taken from Table V only.
  Table IX's BeO and Sc2/3O rows are pair values inferred by additivity from BeAl2O4
  and GSAG, not the gamma of any crystal: they are not extracted (see the record).

All cells were compared with the table images (2026-09-26): no misread cell.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Adair_1989_NonlinearRefractiveIndex"
TABLE_IV = ("tbl_55822e89640f", "tbl_2e1f875bba6b")
TABLE_V = "tbl_62267fd0462a"
TABLE_VII = "tbl_00fcb404f207"
TABLE_IX = "tbl_152efc508e07"
TABLE_I = "tbl_c3b6ccc6cfdc"
TABLE_II = "tbl_0f4b0f5d826e"
TABLE_III = "tbl_fbb684fc9121"
TABLE_VIII = "tbl_8652b8e99339"
TABLE_X = "tbl_5192778e0ba3"
PUMP_UM = "1.064"

# Table IV row label (as extracted) -> (material, dataset id).
TABLE_IV_ROWS = {
    "LiF[100]": ("LiF", "n_LiF"),
    "NaF[100]": ("NaF", "n_NaF"),
    "KF[100]": ("KF", "n_KF"),
    "NaCl[100]": ("NaCl", "n_NaCl"),
    "KCl[100]": ("KCl", "n_KCl"),
    "NaBr[100]": ("NaBr", "n_NaBr"),
    "KBr[100]": ("KBr", "n_KBr"),
    "MgF_2 (o)": ("MgF2", "n_MgF2_o"),
    "CaF_2[100]": ("CaF2", "n_CaF2"),
    "SrF_2[100]": ("SrF2", "n_SrF2"),
    "CdF_2[100]": ("CdF2", "n_CdF2"),
    "BaF_2[100]": ("BaF2", "n_BaF2"),
    "LaF_3 (o)": ("LaF3", "n_LaF3_o"),
    "CeF_3 (o)": ("CeF3", "n_CeF3_o"),
    "AgCl (polycryst.)": ("AgCl", "n_AgCl"),
    "MgO[100]": ("MgO", "n_MgO"),
    "CaO[100]": ("CaO", "n_CaO"),
    "SrO[110]": ("SrO", "n_SrO"),
    "ZnO (e)": ("ZnO", "n_ZnO_e"),
    "ZnO (o)": ("ZnO", "n_ZnO_o"),
    "Al_2O_3 (o)": ("Al2O3", "n_Al2O3_o"),
    "Al_2O_3 (e)": ("Al2O3", "n_Al2O3_e"),
    "Ga_2O_3": ("Ga2O3", "n_Ga2O3"),
    "Y_2O_3": ("Y2O3", "n_Y2O3"),
    "Er_2O_3": ("Er2O3", "n_Er2O3"),
    "SiO_2 (fused)": ("SiO2", "n_SiO2"),
    "SiO_2 (quartz) (o)": ("SiO2_quartz", "n_SiO2_quartz_o"),
    "SiO_2 (quartz) (e)": ("SiO2_quartz", "n_SiO2_quartz_e"),
    "TiO_2": ("TiO2", "n_TiO2"),
    "ZrO_2": ("ZrO2", "n_ZrO2"),
    "BeAl_2O_4": ("BeAl2O4", "n_BeAl2O4"),
    "MgAl_2O_4": ("MgAl2O4", "n_MgAl2O4"),
    "CaMgSi_2O_6": ("CaMgSi2O6", "n_CaMgSi2O6"),
    "YAlO_3 (\N{GREEK SMALL LETTER GAMMA})": ("YAlO3", "n_YAlO3_gamma"),
    "Y_3Al_5O_{12} (YAG)": ("YAG", "n_YAG"),
    "Gd_3Sc_2Al_3O_{12} (GSAG)": ("GSAG", "n_GSAG"),
    "Gd_3Sc_2Ga_3O_{12} (GSGG)": ("GSGG", "n_GSGG"),
    "Gd_3Ga_5O_{12} (GGG)": ("GGG", "n_GGG"),
    "Y_3Ga_5O_{12} (YGG)": ("YGG", "n_YGG"),
    "La_3Lu_2Ga_3O_{12} (LLGG)": ("LLGG", "n_LLGG"),
    "SrTiO_3": ("SrTiO3", "n_SrTiO3"),
    "CaCO_3 (o)": ("CaCO3", "n_CaCO3_o"),
    "CaCO_3 (e)": ("CaCO3", "n_CaCO3_e"),
    "KH2PO4(o)": ("KDP", "n_KDP_o"),
    "KH2PO4(e)": ("KDP", "n_KDP_e"),
    "K(TiO)PO4(x+y)": ("KTP", "n_KTP_xy"),
    "L-arginine phosphate (y)": ("LAP", "n_LAP_y"),
    "L-arginine phosphate (x+z)": ("LAP", "n_LAP_xz"),
    "KTaO3": ("KTaO3", "n_KTaO3"),
    "CaWO4(o)": ("CaWO4", "n_CaWO4_o"),
    "CaWO4(e)": ("CaWO4", "n_CaWO4_e"),
    "ZnS (e)": ("ZnS", "n_ZnS_e"),
    "ZnS (o)": ("ZnS", "n_ZnS_o"),
    "CdS (e)": ("CdS", "n_CdS_e"),
    "CdS (o)": ("CdS", "n_CdS_o"),
}

# The alkali halides of Table V measured here (the others are theory only).
TABLE_V_MEASURED = ("LiF", "NaF", "KF", "NaCl", "KCl", "NaBr", "KBr")

# Table VII formula unit -> material (LiF, NaF and KF come from Table V).
TABLE_VII_ROWS = {
    "Mg_{1/2}F": "MgF2",
    "Ca_{1/2}F": "CaF2",
    "Sr_{1/2}F": "SrF2",
    "Cd_{1/2}F": "CdF2",
    "Ba_{1/2}F": "BaF2",
    "La_{1/3}F": "LaF3",
    "Ce_{1/3}F": "CeF3",
}

# Table IX formula unit -> material; BeO and Sc_{2/3}O are pair values inferred by
# additivity (footnote b; text p. 3347), not a crystal's gamma, and are left out.
TABLE_IX_ROWS = {
    "MgO": "MgO",
    "CaO": "CaO",
    "ZnO": "ZnO",
    "SrO": "SrO",
    "Al_{2/3}O": "Al2O3",
    "Ga_{2/3}O": "Ga2O3",
    "Y_{2/3}O": "Y2O3",
    "Er_{2/3}O": "Er2O3",
    "fused Si_{1/2}O": "SiO2",
    "Ti_{1/2}O": "TiO2",
    "Zr_{1/2}O": "ZrO2",
    "ZnS": "ZnS",
    "CdS": "CdS",
}
TABLE_IX_SKIPPED = ("BeO", "Sc_{2/3}O")

# The index cell: an optional "~" (CeF3), the number, and the footnote letter, printed
# "^a" in the first part of Table IV and as a bare trailing letter in the continuation.
_INDEX_CELL = re.compile(r"^(~?)([\d.]+)\^?([a-k])$")
# A footnote letter trailing a number (Table VII "1.99^b", Table IX "2.26c", "1040.0d").
_MARKED = re.compile(r"^([\d.]+)\^?([a-d])?$")


def _unmarked(cell: str) -> str:
    """Return the number of a cell without its trailing footnote letter."""
    match = _MARKED.match(cell.strip())
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    return number(match.group(1))


def _index() -> dict[str, list[Point]]:
    files = {}
    seen = set()
    for table_id in TABLE_IV:
        for cells in json_table(DIRECTORY, table_id)[2:]:
            material, dataset = TABLE_IV_ROWS[cells[0]]
            seen.add(cells[0])
            match = _INDEX_CELL.match(cells[1])
            if match is None:
                msg = f"unexpected Table IV index cell {cells[1]!r}"
                raise ValueError(msg)
            point = Point("1", "1.06", number(match.group(2)), "", " | ".join(cells))
            files[f"{material}/adair_1989__{dataset}.csv"] = [point]
    if seen != set(TABLE_IV_ROWS):
        msg = f"Table IV rows missing: {set(TABLE_IV_ROWS) - seen}"
        raise ValueError(msg)
    return files


def _gammas() -> dict[str, list[Point]]:
    files = {}
    expt = {}
    for cells in json_table(DIRECTORY, TABLE_V)[2:]:
        if cells[0] not in TABLE_V_MEASURED:
            continue
        printed = " | ".join(cells)
        expt[cells[0]] = number(cells[6])
        files[f"{cells[0]}/adair_1989__gamma_expt_{cells[0]}.csv"] = [
            Point("1", "1.064", number(cells[6]), "", printed)
        ]
        files[f"{cells[0]}/adair_1989__gamma_theory_{cells[0]}.csv"] = [
            Point("1", "1.06", number(cells[5]), "", printed)
        ]
    for cells in table(DIRECTORY, TABLE_VII)[1:]:
        printed = " | ".join(cells)
        if cells[0] in expt:
            # Table VII repeats Table V's measured gamma for LiF, NaF and KF.
            if float(_unmarked(cells[2])) != float(expt[cells[0]]):
                msg = f"Table VII {cells[0]} differs from Table V"
                raise ValueError(msg)
            continue
        material = TABLE_VII_ROWS[cells[0]]
        files[f"{material}/adair_1989__gamma_expt_{material}.csv"] = [
            Point(
                "1",
                "1.064",
                _unmarked(cells[2]),
                "",
                printed,
                {"separation": _unmarked(cells[1])},
            )
        ]
    for cells in table(DIRECTORY, TABLE_IX)[1:]:
        if cells[0] in TABLE_IX_SKIPPED:
            continue
        material = TABLE_IX_ROWS[cells[0]]
        files[f"{material}/adair_1989__gamma_expt_{material}.csv"] = [
            Point(
                "1",
                "1.064",
                _unmarked(cells[2]),
                "",
                " | ".join(cells),
                {"separation": _unmarked(cells[1])},
            )
        ]
    return files


def _n2_table_iv() -> tuple[dict[str, list[Point]], dict[str, str]]:
    """Return Table IV's measured and calculated n2, and the measured value by row."""
    files = {}
    measured = {}
    for table_id in TABLE_IV:
        for cells in json_table(DIRECTORY, table_id)[2:]:
            material, dataset = TABLE_IV_ROWS[cells[0]]
            suffix = dataset.removeprefix("n_")
            printed = " | ".join(cells)
            measured[cells[0]] = number(cells[3])
            files[f"{material}/adair_1989__n2_{suffix}.csv"] = [
                Point("1", PUMP_UM, number(cells[3]), "", printed)
            ]
            if cells[4]:
                files[f"{material}/adair_1989__n2_bgo_{suffix}.csv"] = [
                    Point("1", PUMP_UM, number(cells[4]), "", printed)
                ]
    return files, measured


def _n2_sf58() -> dict[str, list[Point]]:
    """Return SF-58's n2 against the offset Delta (Table I) and the power (Table II)."""
    offset = [
        Point(str(i), PUMP_UM, number(c[1]), "", " | ".join(c), {"delta": number(c[0])})
        for i, c in enumerate(json_table(DIRECTORY, TABLE_I)[2:], start=1)
    ]
    power = [
        Point(str(i), PUMP_UM, number(c[1]), "", " | ".join(c), {"power": number(c[0])})
        for i, c in enumerate(table(DIRECTORY, TABLE_II)[1:], start=1)
    ]
    return {
        "SF-58/adair_1989__n2_offset_SF-58.csv": offset,
        "SF-58/adair_1989__n2_power_SF-58.csv": power,
    }


# Table III row label -> (Table IV row label, material); Al2O3 has no single row.
TABLE_III_ROWS = {
    "LiF": ("LiF[100]", "LiF"),
    "NaF": ("NaF[100]", "NaF"),
    "NaCl": ("NaCl[100]", "NaCl"),
    "NaBr": ("NaBr[100]", "NaBr"),
    "KCl": ("KCl[100]", "KCl"),
    "KBr": ("KBr[100]", "KBr"),
    "CaF_2": ("CaF_2[100]", "CaF2"),
    "SrF_2": ("SrF_2[100]", "SrF2"),
    "CdF_2": ("CdF_2[100]", "CdF2"),
    "BaF_2": ("BaF_2[100]", "BaF2"),
    "LaF_3": ("LaF_3 (o)", "LaF3"),
    "CeF_3": ("CeF_3 (o)", "CeF3"),
    "Al_2O_3": ("", "Al2O3"),
    "fused SiO_2": ("SiO_2 (fused)", "SiO2"),
    "Y_3Al_5O_{12}": ("Y_3Al_5O_{12} (YAG)", "YAG"),
    "MgO": ("MgO[100]", "MgO"),
}
# Table III literature column -> (dataset prefix, wavelength um, second wavelength um
# of a three-wave-mixing column or None), from the column headers (image checked).
TABLE_III_COLUMNS = {
    2: ("tri", "1.06", None),
    3: ("maker_terhune", "0.693", "0.746"),
    4: ("levenson_bloembergen", "0.575", "0.610"),
    5: ("lynch", "0.575", "0.592"),
    6: ("levenson", "0.560", "0.590"),
    7: ("smith_pdf", "1.06", None),
}
_AVERAGE = " (av)"


def _n2_table_iii(measured: dict[str, str]) -> dict[str, list[Point]]:
    """Return Table III's literature columns and its one value not in Table IV."""
    files = {}
    for cells in json_table(DIRECTORY, TABLE_III)[2:]:
        row, material = TABLE_III_ROWS[cells[0]]
        printed = " | ".join(cells)
        this_work = number(cells[1].removesuffix(_AVERAGE))
        if row:
            # "This work" repeats Table IV's measured value (also where "(av)").
            if float(this_work) != float(measured[row]):
                msg = f"Table III {cells[0]} differs from Table IV"
                raise ValueError(msg)
        else:
            # Al2O3 "1.25 (av)": in no other table (Table IV has o 1.23, e 1.30).
            files[f"{material}/adair_1989__n2_av_{material}.csv"] = [
                Point("1", "1.06", this_work, "", printed, {"omega2": "1.07"})
            ]
        for column, (prefix, wavelength, second) in TABLE_III_COLUMNS.items():
            if not cells[column]:
                continue
            extra = {} if second is None else {"omega2": second}
            point = Point("1", wavelength, number(cells[column]), "", printed, extra)
            files[f"{material}/adair_1989__n2_{prefix}_{material}.csv"] = [point]
    return files


# Table VIII row label -> material (all cubic; polarization along [100]).
TABLE_VIII_ROWS = {
    "LiF": "LiF",
    "NaF": "NaF",
    "KF": "KF",
    "NaCl": "NaCl",
    "KCl": "KCl",
    "NaBr": "NaBr",
    "KBr": "KBr",
    "CaF_2": "CaF2",
    "SrF_2": "SrF2",
    "CdF_2": "CdF2",
    "BaF_2": "BaF2",
    "MgO": "MgO",
    "CaO": "CaO",
}


def _n2_table_viii() -> dict[str, list[Point]]:
    """Return Table VIII's n2 column, one file per crystal."""
    files = {}
    for cells in table(DIRECTORY, TABLE_VIII)[1:]:
        material = TABLE_VIII_ROWS[cells[0]]
        files[f"{material}/adair_1989__n2_table8_{material}.csv"] = [
            Point("1", PUMP_UM, number(cells[1]), "", " | ".join(cells))
        ]
    return files


# Table X row label -> (Table IV row label, material). Table X prints
# "Gd3Sc2Al2O12" (image checked) for Table IV's Gd3Sc2Al3O12, GSAG: a misprint.
TABLE_X_ROWS = {
    r"\text{BeAl}_{2}\text{O}_{4}^{\text{a}}": ("BeAl_2O_4", "BeAl2O4"),
    r"\text{MgAl}_{2}\text{O}_{4}": ("MgAl_2O_4", "MgAl2O4"),
    r"\text{CaMgSi}_{2}\text{O}_{6}": ("CaMgSi_2O_6", "CaMgSi2O6"),
    r"\text{YAlO}_{3}": ("YAlO_3 (\N{GREEK SMALL LETTER GAMMA})", "YAlO3"),
    r"\text{Y}_{3}\text{Al}_{5}\text{O}_{12}": ("Y_3Al_5O_{12} (YAG)", "YAG"),
    r"\text{Gd}_{3}\text{Sc}_{2}\text{Al}_{2}\text{O}_{12}^{\text{a}}": (
        "Gd_3Sc_2Al_3O_{12} (GSAG)",
        "GSAG",
    ),
    r"\text{Gd}_{3}\text{Sc}_{2}\text{Ga}_{3}\text{O}_{12}": (
        "Gd_3Sc_2Ga_3O_{12} (GSGG)",
        "GSGG",
    ),
    r"\text{Gd}_{3}\text{Ga}_{5}\text{O}_{12}": ("Gd_3Ga_5O_{12} (GGG)", "GGG"),
    r"\text{Y}_{3}\text{Ga}_{5}\text{O}_{12}": ("Y_3Ga_5O_{12} (YGG)", "YGG"),
    r"\text{La}_{3}\text{Lu}_{2}\text{Ga}_{3}\text{O}_{12}": (
        "La_3Lu_2Ga_3O_{12} (LLGG)",
        "LLGG",
    ),
    r"\text{SrTiO}_{3}": ("SrTiO_3", "SrTiO3"),
}


def _n2_table_x(measured: dict[str, str]) -> dict[str, list[Point]]:
    """Return Table X's n2 calculated by additivity (its measured n2 is Table IV's)."""
    files = {}
    for cells in table(DIRECTORY, TABLE_X)[1:]:
        row, material = TABLE_X_ROWS[cells[0]]
        if float(number(cells[2])) != float(measured[row]):
            msg = f"Table X {material} differs from Table IV"
            raise ValueError(msg)
        if cells[3]:
            files[f"{material}/adair_1989__n2_additivity_{material}.csv"] = [
                Point("1", PUMP_UM, number(cells[3]), "", " | ".join(cells))
            ]
    return files


def _n2() -> dict[str, list[Point]]:
    """Return every n2 file (Tables I-IV, VIII and X)."""
    files, measured = _n2_table_iv()
    return {
        **files,
        **_n2_sf58(),
        **_n2_table_iii(measured),
        **_n2_table_viii(),
        **_n2_table_x(measured),
    }


def extract() -> dict[str, list[Point]]:
    """Return the n2, Table IV's linear index and the hyperpolarizabilities."""
    return {**_index(), **_gammas(), **_n2()}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
