# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Jansonas et al. (2022): n2 of 15 infrared materials at 1.03, 2.2, 3.2 um.

Table 2 (1.03 um) and Table 3 (2.2 and 3.2 um) give the effective n2 (n2_eff, which in a
noncentrosymmetric crystal includes the cascaded second-order contribution) of this
work, and a column of literature values; Table 4 gives n2_Kerr, the effective n2 with
the computed cascaded part removed, for GaSe, AgGaS2, ZnGeP2, KTA and LGS, and
literature n2_Kerr the authors computed from published n2_eff. All in 10^-14 cm^2/W.

Structural fixes, read from the table images (``tables/tbl_8fc389324d09.jpg``,
``tbl_df33b8fa6ec8.jpg``, ``tbl_59933d97ff63.jpg``):

- The sub-header row "(this work) ... (literature)" is merged into the first data row
  ("(this work)0.185 ± 0.023"); the labels are stripped.
- A literature cell holds several values, each printed on its own line as
  "value [± sigma]_(wavelength) [refs]"; the extracted cell runs them together
  ("0.81 ± 0.23(2µm) [29]1 ± 0.15(2.4µm) [29]"), so the script splits it with a pattern
  and checks that nothing is left over.
- The superscript b (Table 2: "same value for both polarization modes"; Table 4:
  "calculated numerically by simulating the propagation") is extracted as a leading
  "b"; it is dropped from the value and kept in ``printed``. Likewise the "≈" of one
  literature value.
- Table 2's "<0.64" (LGS (o), Ref. [48]) is an upper bound, not a value: not extracted.

The pulse durations, "in vicinity of 230 fs, 170 fs and 190 fs for 1.03 um, 2.2 um and
3.2 um pulses" (p. 30515), are written as the per-row column ``pulse_fs``.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Jansonas_2022_InterferometricMeasurementsNonlinear"
TABLE_2 = "tbl_8fc389324d09"
TABLE_3 = "tbl_df33b8fa6ec8"
TABLE_4 = "tbl_59933d97ff63"

# The row label as printed -> (material id, file key).
ROWS = {
    "GaSe": ("GaSe", "GaSe"),
    "AGS (o)": ("AgGaS2", "AgGaS2_o"),  # codespell:ignore ags
    "AGS (e)": ("AgGaS2", "AgGaS2_e"),  # codespell:ignore ags
    "ZGP (o)": ("ZnGeP2", "ZnGeP2_o"),
    "ZGP (e)": ("ZnGeP2", "ZnGeP2_e"),
    "SBN-61 (o)": ("SBN-61", "SBN-61_o"),
    "KTA (o)": ("KTA", "KTA_o"),
    "KTA (e)": ("KTA", "KTA_e"),
    "LGS (o)": ("LGS", "LGS_o"),
    "LGS (e)": ("LGS", "LGS_e"),
    "GaAs": ("GaAs", "GaAs"),
    "BGSe": ("BGSe", "BGSe"),
    "Si": ("Si", "Si"),
    "ZnS": ("ZnS", "ZnS"),
    "ZnTe": ("ZnTe", "ZnTe"),
    "ZnSe": ("ZnSe", "ZnSe"),
    "KRS-5": ("KRS-5", "KRS-5"),
    "KRS-6": ("KRS-6", "KRS-6"),
    "UVFS": ("SiO2", "SiO2"),
}

# Reference numbers as printed -> the file-name suffix of the reproduced dataset.
REFERENCES = {
    "28": "patwardhan_2021",
    "29": "suminas_2020",
    "30": "pigeon_2020",
    "32": "bristow_2007",
    "33": "li_2017",
    "36": "ensley_2019",
    "45,46": "li_2000_2001",
    "47": "mero_2019",
    "48": "namboodiri_2021",
    "49": "sheik_bahae_1991",
    "50": "milam_1998",
    "52": "marcinkeviciute_2018",
}

# Pulse duration by wavelength (um, as the tables print it), p. 30515.
PULSE_FS = {"1.03": "230", "2.2": "170", "3.2": "190"}

MARKS = ("(this work)", "(literature)")
DASH = "-"
_WAVELENGTH = re.compile(r"([\d.]+)\s*µm")
_ITEM = re.compile(
    r"\s*(?P<mark>b)?\s*(?P<approx>≈)?\s*(?P<bound><)?(?P<value>[\d.]+)"
    r"(?:\s*±\s*(?P<sigma>[\d.]+))?\((?P<um>[\d.]+)µm\)\s*\[(?P<refs>[\d,]+)\]\s*"
)


def _clean(cell: str) -> str:
    """Strip the merged sub-header labels (module docstring)."""
    text = cell.strip()
    for mark in MARKS:
        text = text.removeprefix(mark)
    return text.strip()


def _measured(cell: str) -> tuple[str, str] | None:
    """Return (value, sigma) of a this-work cell, or None for an empty or dashed one."""
    text = _clean(cell).removeprefix("b").strip()
    if text in {"", DASH}:
        return None
    value, sigma = text.split("±")
    return number(value), number(sigma)


def _literature(cell: str) -> list[re.Match[str]]:
    """Split a literature cell into its values, checking that nothing is left over."""
    text = _clean(cell)
    items = list(_ITEM.finditer(text))
    if "".join(item.group(0) for item in items).strip() != text:
        msg = f"unparsed literature cell {cell!r}"
        raise ValueError(msg)
    return items


def _table(
    files: dict[str, list[Point]], table_id: str, kind: str, lit_kind: str
) -> None:
    """Add one table's this-work columns and its literature column to ``files``."""
    rows = table(DIRECTORY, table_id)
    header = rows[0]
    wavelengths = {}
    for index in range(1, len(header) - 1):
        match = _WAVELENGTH.search(header[index])
        wavelengths[index] = match.group(1) if match else "1.03"  # Table 2: caption
    for number_, cells in enumerate(rows[1:], start=1):
        material, key = ROWS[cells[0].strip()]
        printed = " | ".join(cells)
        tag = table_id[-4:]
        for index, um in wavelengths.items():
            measured = _measured(cells[index])
            if measured is None:
                continue
            files.setdefault(f"{material}/jansonas_2022__{kind}_{key}.csv", []).append(
                Point(
                    f"{tag}-{number_}-{um}",
                    um,
                    *measured,
                    printed,
                    {"pulse_fs": PULSE_FS[um]},
                )
            )
        for k, item in enumerate(_literature(cells[-1]), start=1):
            if item.group("bound"):
                continue  # an upper bound (module docstring)
            name = REFERENCES[item.group("refs")]
            relative = f"{material}/jansonas_2022__{lit_kind}_{key}_{name}.csv"
            files.setdefault(relative, []).append(
                Point(
                    f"{tag}-{number_}-{k}",
                    item.group("um"),
                    number(item.group("value")),
                    number(item.group("sigma")) if item.group("sigma") else "",
                    printed,
                )
            )


def extract() -> dict[str, list[Point]]:
    """Return n2_eff, n2_Kerr and the literature values, per material and wave."""
    files: dict[str, list[Point]] = {}
    _table(files, TABLE_2, "eff", "lit_eff")
    _table(files, TABLE_3, "eff", "lit_eff")
    _table(files, TABLE_4, "kerr", "lit_kerr")
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
