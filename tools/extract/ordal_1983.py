# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Ordal et al. (1983): the compiled n, k, epsilon of Al, Cu, Au and Ag.

Tables 1 (Al), 2 (Cu), 3 (Au) and 5 (Ag) of Appl. Opt. 22, 1099 (1983), pp. 1103-1109
(PDF pages 5-11), are compilations: each is a run of blocks, one per original source,
each block under its own caption ("TABLE 2. Cu, COPPER (Continued)") and reference
line, with the columns omega (cm-1), lambda (um), -epsilon_1, epsilon_2, n and k to
three significant figures. The library splits the pages into one table per block
(39 in all, 155 files), but its captions are not reliable: a block's caption and
reference line are often attached to the block above or below it, or lost. Which
block is which was therefore read from the page images, where each caption sits
directly above its block (all checked 2026-09-28):

- ``tbl_6418218be243`` is labelled "TABLE 3. Au, GOLD" by the library but is the
  second block of p. 1105, copper after Dold and Mecke (the gold table starts with
  the third block, ``tbl_d345d262b8b8``, which the library leaves unlabelled);
- ``tbl_7fb8ff587269`` carries Motulevich and Shubin's reference, correctly, but
  ``tbl_9c35180248e0`` is Padalka and Shklyarevskii's (the library attaches that
  reference to the block above as a note);
- ``tbl_fd6106702469``, ``tbl_d345d262b8b8``, ``tbl_7bbc87402022``,
  ``tbl_d3ef91129a37``, ``tbl_21d434b5eb92``, ``tbl_ea7ae411dc7b`` and
  ``tbl_3358809f5c2b`` have no caption in the library.

The script checks each block's first wavenumber and row count against the image so a
re-extraction that reorders the blocks fails loudly.

The header cells for epsilon are misread by the library ("-€ 1", "-61", "-<1" for the
image's "-epsilon_1"): a structural fix, the columns are taken by position after
checking the omega, lambda, n and k headers. The fourth column is -epsilon_1: its values
are stored as printed and the dataset's scale of -1 makes them epsilon_1. lambda has no
dataset; it stays in ``printed`` with the row.

Every row of these twenty blocks satisfies epsilon_1 = n^2 - k^2 and epsilon_2 = 2nk,
and lambda = 1e4 / omega, to the printed rounding (checked 2026-09-28); unlike Ordal et
al. (1988), no |epsilon| is misprinted by a factor of ten.
"""

from typing import NamedTuple

from _library import Point, number, table, write

DIRECTORY = "Ordal_1983_OpticalPropertiesMetals"


class Block(NamedTuple):
    """One source's block of a metal's table, as the page image shows it."""

    metal: str
    key: str
    table_id: str
    first: str
    rows: int


BLOCKS = (
    # Table 1, Al (pp. 1103-1104).
    Block("Al", "shiles", "tbl_fd6106702469", "3.23E+02", 26),
    Block("Al", "bennett", "tbl_c2e00320bb13", "3.13E+02", 25),
    Block("Al", "schulz", "tbl_7ec42378518f", "1.05E+04", 12),
    # Table 2, Cu (pp. 1104-1105).
    Block("Cu", "schulz", "tbl_18533a8f9142", "1.05E+04", 11),
    Block("Cu", "lenham", "tbl_97a509391e26", "5.00E+02", 10),
    Block("Cu", "robusto", "tbl_594ca77df023", "1.56E+04", 7),
    Block("Cu", "hagemann", "tbl_8b93c46d68a7", "8.07E+02", 10),
    Block("Cu", "dold", "tbl_6418218be243", "1.00E+03", 11),
    # Table 3, Au (pp. 1105-1107).
    Block("Au", "bennett", "tbl_d345d262b8b8", "3.13E+02", 18),
    Block("Au", "schulz", "tbl_17b239fd0da5", "1.05E+04", 7),
    Block("Au", "motulevich", "tbl_7fb8ff587269", "8.33E+02", 11),
    Block("Au", "padalka", "tbl_9c35180248e0", "9.09E+02", 11),
    Block("Au", "bolotin", "tbl_7bbc87402022", "1.05E+03", 12),
    Block("Au", "brandli", "tbl_d3ef91129a37", "3.14E+01", 14),
    Block("Au", "weaver", "tbl_21d434b5eb92", "8.07E+02", 21),
    # Table 5, Ag (pp. 1108-1109).
    Block("Ag", "bennett", "tbl_ea7ae411dc7b", "3.13E+02", 19),
    Block("Ag", "schulz", "tbl_3358809f5c2b", "1.05E+04", 12),
    Block("Ag", "hagemann", "tbl_b0d84ee51ea1", "8.07E+02", 8),
)
# Dataset prefix -> column (the fourth column is -epsilon_1).
COLUMNS = {"n": 4, "k": 5, "eps1": 2, "eps2": 3}
HEADER = {0: "\\omega (cm^{-1})", 1: "\\lambda (\\mu m)", 4: "n", 5: "k"}


def extract() -> dict[str, list[Point]]:
    """Return n, k, -epsilon_1 and epsilon_2 of each block against the wavenumber."""
    files: dict[str, list[Point]] = {}
    for block in BLOCKS:
        rows = table(DIRECTORY, block.table_id)
        header, body = rows[0], rows[1:]
        if any(header[i] != text for i, text in HEADER.items()) or len(header) != 6:
            msg = f"{block.table_id}: unexpected header {header}"
            raise ValueError(msg)
        if body[0][0] != block.first or len(body) != block.rows:
            msg = f"{block.table_id}: not the {block.metal} {block.key} block"
            raise ValueError(msg)
        stem = f"{block.metal}/ordal_1983__{{}}_{block.metal}_{block.key}.csv"
        for name in COLUMNS:
            files[stem.format(name)] = []
        for cells in body:
            wavenumber = number(cells[0])
            printed = " | ".join(cells)
            for name, column in COLUMNS.items():
                files[stem.format(name)].append(
                    Point(wavenumber, wavenumber, number(cells[column]), "", printed)
                )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
