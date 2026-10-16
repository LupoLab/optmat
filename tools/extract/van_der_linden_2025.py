# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract van der Linden et al. (2025): n2 of CaF2, SiO2, N-BK7 and KBr (Table 1).

Table 1 pairs each TWINS measurement of this work (at the 804 nm laser, printed as 800
nm) with one literature Z-scan value. The table image (``tables/tbl_a71e46ef0edc.jpg``)
agrees with the extracted cells; the sample name is printed only on the first row of
each pair, so a literature row takes the material of the row above. Each measurement is
its own dataset (one material each); each literature value is a reproduced dataset named
by its reference.
"""

from _library import Point, number, table, write

DIRECTORY = "VanderLinden_2025_NonlinearRefractiveIndex"
TABLE = "tbl_a71e46ef0edc"

# The sample as printed -> the registered material id.
MATERIALS = {"CaF2": "CaF2", "SiO2": "SiO2", "NBK-7": "N-BK7", "KBr": "KBr"}

# The reference cell -> the dataset-id prefix: [27] is Ensley & Bambha (2019), [41]
# DeSalvo et al. (1996) (the paper's reference list).
REFERENCES = {"This work": "twins", "[27]": "ensley_2019", "[41]": "desalvo_1996"}


def _value_sigma(cell: str) -> tuple[str, str]:
    """Split a printed "value ± sigma"."""
    value, sigma = cell.split("±")
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return one file per measurement and per reproduced literature value."""
    files: dict[str, list[Point]] = {}
    sample = ""
    for index, cells in enumerate(table(DIRECTORY, TABLE)[1:], start=1):
        sample = cells[0].strip() or sample
        material = MATERIALS[sample]
        reference = REFERENCES[cells[4].strip()]
        # A reproduced value belongs to the material its original names: Ensley and
        # Bambha's glass is "BK7" (ensley_2019), though it sits on the NBK-7 row here.
        if reference == "ensley_2019" and material == "N-BK7":
            material = "BK7"
        value, sigma = _value_sigma(cells[1])
        name = f"{reference}_{material}"
        relative = f"{material}/van_der_linden_2025__{name}.csv"
        point = Point(str(index), number(cells[2]), value, sigma, " | ".join(cells))
        files.setdefault(relative, []).append(point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
