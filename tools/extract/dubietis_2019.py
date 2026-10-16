# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Dubietis & Couairon (2019): the material tables of the supercontinuum book.

Four tables, all checked cell by cell against their images (they agree):

- Table 4.1 (p. 50, reproduced from Dubietis et al., Lith. J. Phys. 57, 113 (2017)):
  band gap U_g, transmission range, n2 and n0 at 800 nm, and the zero-GVD wavelength of
  twelve media. Water (a liquid) is skipped. The transmission range and the zero-GVD
  wavelength have no quantity and are not extracted.
- Table 4.2 (p. 52, from Gallais & Commandre, Appl. Opt. 53, A186 (2014)): single-shot
  damage fluences at 1030 nm, 500 fs; the materials are the header row.
- Table 5.1 (p. 67, glasses) and Table 5.2 (p. 82, semiconductors): U_g and n2 (the
  wavelength of n2 is not stated).

Each value is its own one-point dataset (one quantity, one material). A band gap has no
abscissa; n0 and n2 of Table 4.1 are at the caption's 800 nm; n2 of Tables 5.1 and 5.2
at a wavelength the book does not state.
"""

from _library import Point, number, table, write

DIRECTORY = "Dubietis_2019_UltrafastSupercontinuumGeneration"
TABLE_4_1 = "tbl_03e88cdc5b30"
TABLE_4_2 = "tbl_cbdf49c2bfe4"
TABLE_5_1 = "tbl_04e87619e528"
TABLE_5_2 = "tbl_87900098e048"

# The printed material names and the registry's ids; None is skipped (water, a liquid).
MATERIALS = {
    "LiF": "LiF",
    "CaF2": "CaF2",
    "Al2O3": "Al2O3",
    "BaF2": "BaF2",
    "SiO2 (FS)": "SiO2",
    "SiO2": "SiO2",
    "KDP": "KDP",
    "H2O": None,
    "YAG": "YAG",
    "β-BBO": "BBO",
    "BK7": "BK7",
    "KGW": "KGW",
    "YVO4": "YVO4",
    "Fused silica": "SiO2",
    "ZBLAN": "ZBLAN",
    "La glass": "La_glass",
    "Te glass": "Te_glass",
    "As2S3": "As2S3",
    "ZnS": "ZnS",
    "ZnSe": "ZnSe",
    "GaAs": "GaAs",
    "Si": "Si",
}
NOT_STATED = "not stated"


def _file(material: str, dataset: str) -> str:
    return f"{material}/dubietis_2019__{dataset}.csv"


def _gap_and_n2(
    table_id: str, prefix: str, n2_wavelength: str, n0: bool
) -> dict[str, list[Point]]:
    """Return one band-gap, one n2 and (Table 4.1) one n0 dataset per material."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, table_id)[1:]:
        material = MATERIALS[cells[0]]
        if material is None:
            continue
        printed = " | ".join(cells)
        files[_file(material, f"{prefix}_gap_{material}")] = [
            Point("1", "", number(cells[1]), "", printed)
        ]
        files[_file(material, f"{prefix}_n2_{material}")] = [
            Point("1", n2_wavelength, number(cells[3]), "", printed)
        ]
        if n0:
            files[_file(material, f"{prefix}_n0_{material}")] = [
                Point("1", "800", number(cells[4]), "", printed)
            ]
    return files


def _damage() -> dict[str, list[Point]]:
    """Return Table 4.2's fluences: materials in the header row, values below."""
    header, values = table(DIRECTORY, TABLE_4_2)
    files = {}
    for name, value in zip(header[1:], values[1:], strict=True):
        material = MATERIALS[name]
        if material is None:
            msg = f"Table 4.2 names a skipped material {name!r}"
            raise ValueError(msg)
        printed = f"{name} | {values[0]} | {value}"
        files[_file(material, f"t4_2_damage_{material}")] = [
            Point("1", "1030", number(value), "", printed)
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return every dataset's file."""
    return {
        **_gap_and_n2(TABLE_4_1, "t4_1", "800", n0=True),
        **_damage(),
        **_gap_and_n2(TABLE_5_1, "t5_1", NOT_STATED, n0=False),
        **_gap_and_n2(TABLE_5_2, "t5_2", NOT_STATED, n0=False),
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
