# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Leonard (1974): the inert-gas compilation (Tables I-VI).

All tables are read from the library's ``tables/<id>.json`` cell files, because the
library writes no CSV for tables with merged cells (every table here but three).

Table I gives, per gas and per original author, the absolute ``(n - 1) x 10^6`` at
``lambda_vac = 5462.258`` A recalculated by Leonard to his STP (273.16 K, 1 atm, second
virial coefficients), with its probable error. Table II gives the dispersion of each
author's measurements, rescaled to the recommended absolute value of Table I, at
wavelengths **in standard air (15 degrees C)**. Tables III-VI hold values measured by
other techniques, the Verdet constants and the static polarizabilities.

Structural fixes, all read from the table images:

- Table II: wherever the original prints one wavelength for several values (a brace),
  the library puts the wavelength on the wrong row, and for krypton it merges three
  stacked cells into one (``2446.9052379.9902345.55``) and moves the braced values into
  the wavelength cell. The *values* and (except for krypton) the *reference numbers*
  are in the library in the printed reading order; only the wavelengths are misplaced.
  The script therefore takes the values in reading order from the library and assigns
  each its wavelength and reference from an explicit list read from the image
  (``TABLE_II``), asserting that every wavelength exists in the library's column, that
  the value counts agree and (except for krypton) that the library's reference column
  equals the image's reference sequence.
- Table II's rows marked ``*`` (5460.740 A in air) repeat the recommended absolute
  values of Table I; they are not extracted twice (the ``recommended`` datasets are
  taken from Table I, and the script asserts the two agree).
- Table I: the brace joining He Refs. 7 and 18 (one recalculated value for two quoted
  ones) and the brace joining the two Ar rows of Ref. 29 are restored; the underlined
  (recommended) values are listed in ``RECOMMENDED``.
- Table III's header is garbled in the library (the part labels (a)-(d) and the gas
  columns are merged); the column order He, Ne, Ar, Kr, Xe is read from the image.
- Table VI's brace joining the two He values of Ref. 82 is restored.

A value whose reference number is printed in parentheses is "discrepant or low
accuracy" (Table II legend); the datasets of those authors are ``excluded``.
"""

import json
import re
from dataclasses import replace

from _library import LIBRARY, Point, number, write

DIRECTORY = "Leonard_1974_RefractiveIndicesVerdet"
SOURCE = "leonard_1974"

TABLE_I = "tbl_c30192225841"
TABLE_II = {
    "He": "tbl_9d80fa2bdd50",
    "Ne": "tbl_b21803ae40ad",
    "Ar": "tbl_b21803ae40ad",
    "Kr": "tbl_bd18889f4119",
    "Xe": "tbl_688a696b8c54",
}
TABLE_III_A = "tbl_f7fa4d9454e9"
TABLE_III_B = "tbl_10c2bb7d578a"
TABLE_III_C = "tbl_088071093ce3"
TABLE_III_D = "tbl_3f0c43ea8e18"
TABLE_IV = "tbl_ef7f5e3531dc"
TABLE_V = "tbl_1bf9e7617fb4"
TABLE_VI = "tbl_034a45f4d3ef"

# Leonard's reference number -> dataset name of that author's Table II values.
AUTHORS = {
    "8": "scheel_schmidt_1908",
    "9": "tausz_gorlacher_1931",
    "10": "damkohler_1934",
    "11": "larsen_1934",
    "12": "kronjager_1936",
    "13": "koch_1949",
    "14": "peck_fisher_1964",
    "15": "mansfield_peck_1969",
    "18": "cuthbertson_1911",
    "19": "cuthbertson_1932",
    "25": "herrmann_1908",
    "26": "burton_1908",
    "27": "koch_1913",
    "28": "ahrberg_1909",
    "29": "rusch_1923",
    "30": "quarder_1924",
}

# Table II, read from the table images: for each gas and each half of the printed
# table, one line per printed wavelength (in air at 15 degrees C, as printed) followed
# by the reference number of each value printed against it, top to bottom, with the
# parentheses of a discrepant value. A line with several references is a brace; ``3x10``
# is a brace of three values that share one printed reference (Damkohler's krypton).
TABLE_II_KEYS = {
    ("He", 0): """
        20581.277 15
        16940.570 15
        15295.963 15
        15289.355 15
        14752.461 15
        13718.567 15
        10139.787 15
        9657.779 15
        9224.493 15
        9122.962 15
        8264.517 15
        7425.113 15
        6678.149 (8)
        6652.783 26
        6438.4696 (25) 18
        5895.923 26
        5875.618 (8) 26
        5790.654 (8) 18
        5769.59 (25) (8) 18
        5460.740 *
        5210.067 18
        5085.824 18 15
        5015.675 (8) 26
        4916.036 27
        4861.327 26
        4799.918 18 15
    """,
    ("He", 3): """
        4358.35 (25) (8) 19
        4108.066 27
        4077.811 19
        4046.561 19
        3983.977 27
        3906.410 19
        3861.08 27
        3680.008 27
        3663.276 19
        3543.079 27
        3341.478 27 19
        3131.833 19
        3125.663 19
        3021.499 19
        2967.278 27 19
        2925.406 27 19
        2893.595 27 19
        2856.969 27
        2759.712 27
        2752.775 27 19
        2674.99 27
        2576.295 27
        2536.519 27
        2464.068 27
        2446.905 27
        2379.990 27
        2302.09 27
    """,
    ("Ne", 0): """
        6438.4696 18
        5790.654 18
        5769.59 18
        5460.740 *
        5210.067 18
        5085.824 18
        4916.036 19
        4799.918 18
        4358.35 19
    """,
    ("Ne", 3): """
        4077.811 19
        4046.561 19
        3906.410 19
        3663.276 19
        3341.478 19
        3131.833 19
        3021.499 19
        2967.278 19
        2893.595 19
    """,
    ("Ar", 0): """
        67628.5 29
        20581.277 14
        16940.570 14
        15295.963 14
        15289.355 14
        14752.461 14
        13718.567 14
        10139.787 14
        9657.779 14
        9224.493 14
        9122.962 14
        8264.517 14
        7425.113 14
        7032.409 14
        6562.783 26 29 (9)
        6438.4696 (28) 18
        5895.923 26
        5875.618 26 (9)
        5790.654 (28) 18 29
        5782.158 (30)
        5769.59 (28) 18 29
        5700.240 (30)
        5677.17 11
        5460.740 *
        5218.202 (30)
        5210.067 18
        5153.235 (30)
    """,
    ("Ar", 3): """
        5105.541 (30)
        5085.824 18 14
        5015.675 26
        4916.036 11
        4861.327 26
        4799.918 18 14
        4678.155 14
        4651.134 (30)
        4358.35 (28) (9)
        4275.131 (30)
        4108.066 11
        3801.658 11
        3349.292 (30)
        3341.478 11
        2961.165 (30)
        2925.406 11
        2893.565 11
        2824.369 (30)
        2766.371 (30)
        2759.712 11
        2752.775 11
        2674.99 11
        2576.295 11
        2492.146 (30)
        2464.068 11
        2446.905 11
        2441.637 (30)
        2379.990 11
        2345.55 11
        2302.09 11
    """,
    ("Kr", 0): """
        6707.844 18
        6562.783 10
        6462.5 10
        6438.4696 18
        6431.2 10
        6287.5 10
        6234.37 13
        6131.2 10
        6123.27 13
        6072.62 13
        5987.5 10
        5890.6 10
        5850.0 10
        5790.654 18
        5769.59 18
        5712.5 10
        5677.17 13
        5662.5 10
        5575.0 10
        5460.740 *
        5337.5 10
        5268.7 3x10
        5225.0 10
        5210.067 18
        5125.0 10
        5087.5 10
        5085.824 18 12
        4984.4 10
        4918.7 2x10
        4916.036 13
        4799.918 18 12
        4756.2 10
        4731.2 10
        4712.5 10
        4681.3 10
        4650.0 10
        4618.8 10
        4587.5 10
    """,
    ("Kr", 3): """
        4531.2 10
        4512.5 10
        4493.7 10
        4413.042 12
        4406.2 10
        4108.066 13
        3801.658 13
        3403.653 12
        3341.478 13
        2980.628 12
        2925.406 13
        2893.595 13
        2856.969 13
        2836.907 12
        2759.712 13
        2752.775 13
        2748.58 12
        2674.99 13
        2576.295 13
        2573.09 12
        2464.068 13
        2446.905 13
        2379.990 13
        2345.55 13
        2329.282 12
        2302.09 13
        2265.017 12
        2194.63 12
        2144.382 12
        2025.51 12
    """,
    ("Xe", 0): """
        6707.844 18
        6562.783 10
        6473.5 10
        6441.2 10
        6438.4696 18
        6408.8 10
        6376.5 10
        6353.0 10
        6308.8 10
        6234.37 13
        6201.5 10
        6144.7 10
        6123.27 13
        6072.62 13
        6032.4 10
        5932.4 10
        5898.5 10
        5808.4 10
        5790.654 18
        5779.4 10
        5769.59 18
        5682.4 10
        5677.17 13
        5655.9 10
        5566.2 10
        5535.3 10
        5477.9 10
        5460.740 *
        5450.0 10
        5426.5 10
        5347.1 10
        5314.7 10
        5241.2 10
        5214.7 10
        5210.067 18
        5141.2 10
        5717.7 10
        5085.824 18 (12)
        5041.2 10
        5019.1 10
        4953.0 10
        4932.4 10
        4916.036 13
        4864.7 10
    """,
    ("Xe", 3): """
        4850.0 10
        4799.918 18 (12)
        4782.4 10
        4766.2 10
        4729.4 10
        4711.8 10
        4696.5 10
        4666.2 10
        4632.4 10
        4614.7 10
        4588.3 10
        4573.5 10
        4541.2 10
        4529.4 10
        4525.0 10
        4508.9 10
        4413.042 (12)
        4108.066 13
        3801.658 13
        3403.653 (12)
        3341.478 13
        2980.628 (12)
        2925.406 13
        2893.595 13
        2856.969 13
        2836.907 (12)
        2759.712 13
        2752.775 13
        2748.58 (12)
        2674.99 13
        2576.295 13
        2573.09 (12)
        2464.068 13
        2446.905 13
        2379.990 13
        2345.55 13
        2329.282 (12)
        2302.09 13
        2265.017 (12)
        2194.63 (12)
        2144.382 (12)
    """,
}

# Printed wavelengths the library lost entirely (the image shows them).
NOT_IN_LIBRARY = {("Ar", "4861.327")}  # 283.9 (Ref. 26), right half of argon

# Table II rows of each half-table (the Ne and Ar parts share one library table).
TABLE_II_ROWS = {
    "He": slice(2, None),
    "Ne": slice(1, 10),
    "Ar": slice(11, None),
    "Kr": slice(1, None),
    "Xe": slice(1, None),
}

# Single Table II values split into their own dataset (see the source record):
# Rusch's absolute infrared value, which Leonard says "cannot be properly corrected to
# STP" (p. 4), and Larsen's 4108.066 A value, which falls 1.3 % below its neighbours.
SPLIT = {
    ("Ar", "29", "67628.5"): "rusch_1923_Ar_ir",
    ("Ar", "11", "4108.066"): "larsen_1934_Ar_4108",
}

# The underlined (recommended) value of each gas in Table I, by row id.
RECOMMENDED = {"He": "ref23", "Ne": "ref24", "Ar": "ref14", "Kr": "av", "Xe": "av"}


def grid(table_id: str) -> list[list[str]]:
    """Return a library table's cell grid from its JSON file (merged cells once)."""
    path = LIBRARY / DIRECTORY / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))
    rows = [[""] * cells["columns"] for _ in range(cells["rows"])]
    for cell in cells["cells"]:
        rows[cell["row"]][cell["column"]] = cell["text"].strip()
    return rows


def _plain(cell: str) -> str:
    """Return a printed value without the parentheses of a discrepant value."""
    return number(cell.replace("(", "").replace(")", "").replace("*", ""))


def _keys(block: str) -> list[tuple[str, str]]:
    """Expand an image-key block into one (wavelength, reference) per value."""
    keys = []
    for line in block.strip().splitlines():
        wavelength, *refs = line.split()
        for ref in refs:
            count, _, bare = ref.rpartition("x")
            keys.extend((wavelength, bare) for _ in range(int(count or 1)))
    return keys


def _printed_refs(block: str) -> str:
    """Return the reference numbers as printed (a shared brace once), concatenated."""
    return "".join(
        ref.rpartition("x")[2]
        for line in block.strip().splitlines()
        for ref in line.split()[1:]
    )


# The library's krypton reference column (read cell by cell and concatenated) is the
# printed sequence with the last references lost where it merged stacked rows: two
# ``10`` at the end of the left half, one ``12`` at the end of the right half.
KR_LOST_REFS = {0: "1010", 3: "12"}


def _values_refs(gas: str, half: int) -> tuple[list[str], list[str], str]:
    """Return the value cells and reference cells of one half-table, reading order.

    Also returns the text of the half's wavelength column, for membership checks.
    """
    rows = grid(TABLE_II[gas])[TABLE_II_ROWS[gas]]
    values, refs, wavelengths = [], [], []
    for cells in rows:
        wavelength, value, ref = cells[half : half + 3]
        wavelengths.append(wavelength)
        if gas == "Kr":
            # The library moves braced krypton values into the wavelength cell
            # (duplicating or displacing the value cell) and merges stacked cells;
            # every krypton value is printed with two decimals.
            if r"\begin{array}" in wavelength:
                values.extend(
                    re.findall(r"\d+\.\d+", wavelength.split(r"\begin{array}{c}")[1])
                )
            elif value:
                values.extend(re.findall(r"\d{3}\.\d{2}", value))
            refs.append(ref)
            continue
        if value:
            values.append(value)
            refs.append(ref)
    return values, refs, " ".join(wavelengths)


def _rows(entries: list[tuple[str, str, str]]) -> list[Point]:
    """Turn (wavelength, cell, printed) entries into points with unique row ids.

    A wavelength with several values of one author (Damkohler's krypton readings at
    5268.7 and 4918.7 A) gets rows ``5268.7a``, ``5268.7b``, ... in printed order.
    """
    counts: dict[str, int] = {}
    for wavelength, _, _ in entries:
        counts[wavelength] = counts.get(wavelength, 0) + 1
    seen: dict[str, int] = {}
    points = []
    for wavelength, cell, printed in entries:
        row = wavelength
        if counts[wavelength] > 1:
            row += "abcdefgh"[seen.get(wavelength, 0)]
            seen[wavelength] = seen.get(wavelength, 0) + 1
        points.append(Point(row, wavelength, _plain(cell), "", printed))
    return points


def _table_ii() -> tuple[dict[str, list[Point]], dict[str, str]]:
    """Return one data file per gas and original author of Table II.

    Also returns the value of each gas's starred (recommended) row.
    """
    entries: dict[str, list[tuple[str, str, str]]] = {}
    starred = {}
    for (gas, half), block in TABLE_II_KEYS.items():
        values, library_refs, column = _values_refs(gas, half)
        keys = _keys(block)
        if len(keys) != len(values):
            msg = f"{gas}/{half}: {len(keys)} keys for {len(values)} values"
            raise ValueError(msg)
        if gas == "Kr":
            same = "".join(library_refs) + KR_LOST_REFS[half] == _printed_refs(block)
        else:
            same = [ref for _, ref in keys] == library_refs
        if not same:
            msg = f"{gas}/{half}: reference sequence differs from the library's"
            raise ValueError(msg)
        for (wavelength, ref), cell in zip(keys, values, strict=True):
            if wavelength not in column and (gas, wavelength) not in NOT_IN_LIBRARY:
                msg = f"{gas}: wavelength {wavelength} not in the library column"
                raise ValueError(msg)
            if ref == "*":
                starred[gas] = _plain(cell)
                continue
            bare = ref.strip("()")
            dataset = SPLIT.get((gas, bare, wavelength), f"{AUTHORS[bare]}_{gas}")
            relative = f"{gas}/{SOURCE}__{dataset}.csv"
            entries.setdefault(relative, []).append(
                (wavelength, cell, f"{wavelength} | {cell} | {ref}")
            )
    return {k: _rows(v) for k, v in entries.items()}, starred


def _table_i() -> tuple[dict[str, list[Point]], dict[str, Point]]:
    """Return Table I's recalculated values per gas, and the recommended rows."""
    rows = grid(TABLE_I)[1:]
    by_gas: dict[str, list[Point]] = {}
    gas = ""
    for cells in rows:
        gas = cells[0] or gas
        _quoted, recalculated, accuracy, ref = cells[1:5]
        points = by_gas.setdefault(gas, [])
        printed = " | ".join(cells)
        if not recalculated:
            # He, Ref. 18: braced with Ref. 7 in the image; one recalculated value.
            points[-1] = replace(points[-1], row=f"{points[-1].row}_{ref}")
            continue
        if not ref:
            # Ar, second row of Ref. 29 (a brace in the image).
            row = f"{points[-1].row}b"
            points[-1] = replace(points[-1], row=f"{points[-1].row}a")
        elif ref == "Av.":
            row = "av"
        else:
            row = "ref" + ref.replace(", ", "_")
        sigma = number(accuracy) if accuracy else ""
        points.append(
            Point(row, "5462.258", number(recalculated.rstrip("*")), sigma, printed)
        )
    recommended = {
        g: next(p for p in by_gas[g] if p.row == RECOMMENDED[g]) for g in by_gas
    }
    files = {
        f"{g}/{SOURCE}__absolute_{g}.csv": [p for p in points if p.row != "av"]
        for g, points in by_gas.items()
    }
    return files, recommended


def _value_sigma(cell: str) -> tuple[str, str]:
    r"""Split ``266.0 ± 8`` or ``0.0742 \pm .002`` into value and sigma."""
    parts = re.split(r"\s*(?:±|\\pm)\s*", cell.replace("^{\\dagger}", ""))
    value = number(parts[0])
    sigma = number(parts[1]) if len(parts) > 1 else ""
    return value, sigma


def _table_iii() -> dict[str, list[Point]]:
    """Return Table III (a)-(d), one data file per gas and reference."""
    files: dict[str, list[Point]] = {}
    gases = ("He", "Ne", "Ar", "Kr", "Xe")  # columns 1-5, read from the image

    def add(part: str, gas: str, cells: list[str], wavelength: str, cell: str) -> None:
        value, sigma = _value_sigma(cell)
        relative = f"{gas}/{SOURCE}__{part}_{gas}_ref{cells[6]}.csv"
        files.setdefault(relative, []).append(
            Point(wavelength, wavelength, value, sigma, " | ".join(cells))
        )

    for cells in grid(TABLE_III_A)[2:]:
        wavelength = cells[0].split("±")[0].strip()
        for gas, cell in zip(gases, cells[1:6], strict=True):
            if cell not in {"", "-"}:
                add("rayleigh", gas, cells, wavelength, cell)
    for cells in grid(TABLE_III_B):
        # Ratios to argon (the argon cells are the unit); ref 38's argon value is
        # relative to nitrogen, whose ratio the dagger note gives (added below).
        for gas, cell in zip(gases, cells[1:6], strict=True):
            if gas != "Ar" and cell not in {"", "-"}:
                add("ratio", gas, cells, cells[0], cell)
    for cells in grid(TABLE_III_C):
        wavelength = cells[0].split("±")[0].strip()
        add("cerenkov", "Ar", cells, wavelength, cells[3])
    for cells in grid(TABLE_III_D):
        add("hook", "Ar", cells, cells[0], cells[3])
    return files


def _table_iv() -> dict[str, list[Point]]:
    """Return Table IV (kappa for argon), one data file per reference."""
    files: dict[str, list[Point]] = {}
    for cells in grid(TABLE_IV)[2:]:
        value, sigma = _value_sigma(cells[1])
        relative = f"Ar/{SOURCE}__kappa_Ar_ref{cells[2]}.csv"
        files.setdefault(relative, []).append(
            Point(cells[0], cells[0], value, sigma, " | ".join(cells))
        )
    return files


VERDET_GROUPS = {
    "56, 57": "ingersoll_liebenberg",
    "52, 53, 54": "mallemann_1951",
    "48": "mallemann_early",
    "48, 49, 51": "mallemann_early",
}


def _table_v() -> dict[str, list[Point]]:
    """Return Table V (Verdet constants), one data file per gas and laboratory."""
    files: dict[str, list[Point]] = {}
    gases = ("He", "Ne", "Ar", "Kr", "Xe")
    wavelength = ""
    for cells in grid(TABLE_V)[2:]:
        wavelength = cells[0] or wavelength
        group = VERDET_GROUPS[cells[6]]
        for gas, cell in zip(gases, cells[1:6], strict=True):
            if cell in {"", "-"}:
                continue
            relative = f"{gas}/{SOURCE}__verdet_{group}_{gas}.csv"
            files.setdefault(relative, []).append(
                Point(wavelength, wavelength, number(cell), "", " | ".join(cells))
            )
    return files


def _table_vi() -> dict[str, list[Point]]:
    """Return Table VI (static polarizabilities), per gas and kind of value."""
    files: dict[str, list[Point]] = {}
    rows = grid(TABLE_VI)
    sections = [(rows[2:25], ("He", "Ne", "Ar")), (rows[26:], ("Kr", "Xe"))]
    for block, gases in sections:
        for column, gas in enumerate(gases):
            kind = "experiment"
            last_ref = ""
            for cells in block:
                value, ref = cells[2 * column], cells[2 * column + 1]
                if not value:
                    continue
                if value.startswith("("):
                    kind = "ri_extrapolated" if "extrapolated" in value else "theory"
                    continue
                if not ref:  # He, second value of Ref. 82 (a brace in the image)
                    ref = last_ref
                last_ref = ref
                relative = f"{gas}/{SOURCE}__polarizability_{kind}_{gas}.csv"
                points = files.setdefault(relative, [])
                row = "ref" + ref.replace(",", "_")
                if any(p.row == row for p in points):
                    first = next(p for p in points if p.row == row)
                    points[points.index(first)] = replace(first, row=row + "a")
                    row += "b"
                points.append(Point(row, "inf", number(value), "", f"{value} | {ref}"))
    return files


def extract() -> dict[str, list[Point]]:
    """Return every data file of the source."""
    files, starred = _table_ii()
    absolute, recommended = _table_i()
    for gas, point in recommended.items():
        if starred[gas] != point.value:
            msg = f"{gas}: Table II's starred value {starred[gas]} is not {point.value}"
            raise ValueError(msg)
    files.update(absolute)
    for gas, point in recommended.items():
        files[f"{gas}/{SOURCE}__recommended_{gas}.csv"] = [point]
    files.update(_table_iii())
    files.update(_table_iv())
    files.update(_table_v())
    files.update(_table_vi())
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
