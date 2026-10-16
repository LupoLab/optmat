# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Liu et al. (1978): two-photon absorption at 355 and 266 nm (Tables I, II).

Table I (``tbl_16c90fa7a1af``, p. 4625, 355 nm) and Table II (``tbl_c0af8ebd30c9``,
p. 4626, 266 nm) give per sample: the band gap E_g, the surface normal, the direction
of the laser polarization, the self-focusing critical power P_cr, the sample length l,
beta in cm/MW with its relative uncertainty ("7.29 x 10^{-3} +- 20%" in Table I,
"3.5 \times 10^{-3} \pm 25\%" in Table II), and the index n; Table II also the
effective mass m*/m of the Keldysh fit. Some materials have several rows (another
polarization, another sample length, a coloured sample); each row is its own beta
file. beta cells printed "<" are upper bounds: kept (the number in ``value``, the raw
cell in ``printed``), and the record marks them excluded. beta and m*/m are stored in
e-notation (``7.29e-3``); the relative uncertainty of a beta cell stays in ``printed``,
and the record states it as the dataset's relative uncertainty.

The band gap, P_cr and n are extracted once per material and wavelength, from the
first row that prints them; a gap printed in both tables is checked equal. Index cells
of the uniaxial crystals print both rays, "1.548 (n_0) 1.499 (n_e)": split into the
o and e files. The two rows of deuterated CDA (Table I) print different n_o (1.60 and
1.59): both are kept, as rows "perp_z" and "z" of the o-ray file.

Structural notes (images checked, every cell agrees with the extraction):
- Table I prints sapphire as "A_2O_3" (the image shows "A2O3"): a misprint of the
  label Al2O3, mapped here; no value is affected.
- Footnote letters e-h of Table I are extracted onto the names (``ADA^e``); stripped.
- Table II's KDP surface-normal cell "1.06 \mu \rightarrow 0.53 \mucut" is the image's
  "1.06 µ → 0.53 µ cut" (a doubling-cut crystal); kept in ``printed``.
"""

import re

from _library import MINUS, Point, number, table, write

DIRECTORY = "Liu_1978_AbsoluteTwoPhoton"
TABLE_I = "tbl_16c90fa7a1af"
TABLE_II = "tbl_c0af8ebd30c9"

# (table, row index in the body) -> (beta dataset id, material).
BETA_ROWS = {
    (TABLE_I, 0): ("beta_355_KI", "KI"),
    (TABLE_I, 1): ("beta_355_RbBr", "RbBr"),
    (TABLE_I, 2): ("beta_355_RbI", "RbI"),
    (TABLE_I, 3): ("beta_355_RbI_011", "RbI"),
    (TABLE_I, 4): ("beta_355_RbI_circular", "RbI"),
    (TABLE_I, 5): ("beta_355_ADA", "ADA"),
    (TABLE_I, 6): ("beta_355_ADP", "ADP"),
    (TABLE_I, 7): ("beta_355_DCDA_perp_z", "DCDA"),
    (TABLE_I, 8): ("beta_355_DCDA", "DCDA"),
    (TABLE_I, 9): ("beta_355_CDA", "CDA"),
    (TABLE_I, 10): ("beta_355_KDA", "KDA"),
    (TABLE_I, 11): ("beta_355_DKDA", "DKDA"),
    (TABLE_I, 12): ("beta_355_KDP", "KDP"),
    (TABLE_I, 13): ("beta_355_DKDP", "DKDP"),
    (TABLE_I, 14): ("beta_355_RDA", "RDA"),
    (TABLE_I, 15): ("beta_355_RDP", "RDP"),
    (TABLE_I, 16): ("beta_355_Al2O3", "Al2O3"),
    (TABLE_I, 17): ("beta_355_SiO2", "SiO2"),
    (TABLE_II, 0): ("beta_266_NaCl", "NaCl"),
    (TABLE_II, 1): ("beta_266_NaCl_011", "NaCl"),
    (TABLE_II, 2): ("beta_266_NaBr", "NaBr"),
    (TABLE_II, 3): ("beta_266_KCl", "KCl"),
    (TABLE_II, 4): ("beta_266_KCl_F_centers", "KCl"),
    (TABLE_II, 5): ("beta_266_KCl_l210", "KCl"),
    (TABLE_II, 6): ("beta_266_KBr", "KBr"),
    (TABLE_II, 7): ("beta_266_KI", "KI"),
    (TABLE_II, 8): ("beta_266_RbCl", "RbCl"),
    (TABLE_II, 9): ("beta_266_RbCl_l191", "RbCl"),
    (TABLE_II, 10): ("beta_266_RbBr", "RbBr"),
    (TABLE_II, 11): ("beta_266_RbI", "RbI"),
    (TABLE_II, 12): ("beta_266_CaCO3", "CaCO3"),
    (TABLE_II, 13): ("beta_266_KDP", "KDP"),
    (TABLE_II, 14): ("beta_266_ADP", "ADP"),
    (TABLE_II, 15): ("beta_266_Al2O3", "Al2O3"),
    (TABLE_II, 16): ("beta_266_Al2O3_perp_c", "Al2O3"),
    (TABLE_II, 17): ("beta_266_SiO2", "SiO2"),
    (TABLE_II, 18): ("beta_266_SiO2_quartz", "SiO2_quartz"),
    (TABLE_II, 19): ("beta_266_LiF", "LiF"),
    (TABLE_II, 20): ("beta_266_CaF2", "CaF2"),
}
WAVELENGTH = {TABLE_I: "355", TABLE_II: "266"}
UNIAXIAL = {"ADP", "DCDA", "KDA", "KDP", "DKDP", "RDA", "RDP", "CaCO3", "Al2O3"}

TIMES = "\N{MULTIPLICATION SIGN}"
PLUS_MINUS = "\N{PLUS-MINUS SIGN}"
_BETA = re.compile(
    r"^(?P<bound><)?(?P<mantissa>[\d.]+)\s*(?:" + TIMES + r"|\\times)\s*10\^\{"
    r"(?P<power>[-\d]+)\}(?:\s*(?:" + PLUS_MINUS + r"|\\pm)\s*(?P<percent>\d+)\\%)?$"
)
_POWER = re.compile(r"^([\d.]+)\s*\\times\s*10\^\{([-\d]+)\}$")
_RAYS = re.compile(r"^([\d.]+) \(n_[0o]\)(?: ([\d.]+) \(n_e\))?$")


def _beta(cell: str) -> str:
    """Return a beta cell ``[<] m x 10^{k} [+- p%]`` as ``mek`` (rest printed)."""
    match = _BETA.match(cell.strip().replace(MINUS, "-"))
    if match is None:
        msg = f"unexpected beta cell {cell!r}"
        raise ValueError(msg)
    value = f"{match['mantissa']}e{match['power']}"
    float(value)
    return value


def _mass(cell: str) -> str:
    r"""Return an m*/m cell, ``5.3 \times 10^{-3}`` as ``5.3e-3``, or a plain number."""
    match = _POWER.match(cell.strip())
    if match is None:
        return number(cell)
    return f"{match.group(1)}e{match.group(2)}"


def _index(cell: str) -> dict[str, str]:
    """Return the index of each ray printed in an n cell (o and e, or one value)."""
    match = _RAYS.match(cell.strip())
    if match is None:
        return {"": number(cell)}
    rays = {"o": number(match.group(1))}
    if match.group(2):
        rays["e"] = number(match.group(2))
    return rays


def _indices(
    files: dict[str, list[Point]], material: str, dataset: str, cell: str, row: Point
) -> None:
    """Add the index (per ray) printed in a row's n cell to the material's n files."""
    for ray, n in _index(cell).items():
        if material in UNIAXIAL and not ray:
            msg = f"{material}: an index without its ray"
            raise ValueError(msg)
        suffix = f"_{ray}" if ray else ""
        label = row.x
        if dataset.startswith("beta_355_DCDA"):
            # The two deuterated-CDA rows print different indices at one wavelength.
            label = "perp_z" if dataset.endswith("perp_z") else "z"
        files.setdefault(f"{material}/liu_1978__n_{material}{suffix}.csv", []).append(
            Point(label, row.x, n, "", row.printed)
        )


def _critical_power(
    files: dict[str, list[Point]], material: str, cell: str, row: Point
) -> None:
    """Add a row's P_cr to the material's file, once per wavelength."""
    name = f"{material}/liu_1978__critical_power_{material}.csv"
    if all(p.row != row.x for p in files.get(name, [])):
        files.setdefault(name, []).append(
            Point(row.x, row.x, number(cell), "", row.printed)
        )


def extract() -> dict[str, list[Point]]:
    """Return beta per row, and E_g, P_cr, n and m*/m per material and wavelength."""
    files: dict[str, list[Point]] = {}
    gaps: dict[str, tuple[str, str]] = {}
    for table_id in (TABLE_I, TABLE_II):
        rows = table(DIRECTORY, table_id)
        header, body = rows[0], rows[1:]
        if not header[6].startswith(r"\beta(cm/MW)"):
            msg = f"unexpected header {header!r}"
            raise ValueError(msg)
        wavelength = WAVELENGTH[table_id]
        for position, cells in enumerate(body):
            dataset, material = BETA_ROWS[table_id, position]
            row = Point("1", wavelength, _beta(cells[6]), "", " | ".join(cells))
            files[f"{material}/liu_1978__{dataset}.csv"] = [row]
            if (
                cells[1]
                and gaps.setdefault(material, (cells[1], row.printed))[0] != (cells[1])
            ):
                msg = f"{material}: band gaps differ between the tables"
                raise ValueError(msg)
            if cells[4]:
                _critical_power(files, material, cells[4], row)
            if cells[-1]:
                _indices(files, material, dataset, cells[-1], row)
            if table_id == TABLE_II and cells[7]:
                files[f"{material}/liu_1978__effective_mass_{material}.csv"] = [
                    Point("1", "", _mass(cells[7]), "", row.printed)
                ]
    for material, (gap, printed) in gaps.items():
        files[f"{material}/liu_1978__band_gap_{material}.csv"] = [
            Point("1", "", number(gap), "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
