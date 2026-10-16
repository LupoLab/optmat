# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Potts & Williams (1974): He II photoelectron bands of N2, CO, CO2 and N2O.

Table 1 (``tbl_5eb7440f6f19``, p. 7) and Table 3 (``tbl_278510f78d78``, p. 14) have no
CSV rendering; their cells are read from the JSON (images checked). Each molecule is a
block of four rows: "Ionization potential in eV" (band maxima or vertical ionisation
energies, per the caption note), "Relative intensity of band", "Ionized orbital or
ionization process" and "XPS IP in eV" (other authors' X-ray values). One column per
band. Only the curated species are taken: N2 and CO (Table 1) and CO2 and N2O
(Table 3); HCN, COS and CS2 are not curated (see the record).

The caption note (p. 7) says the ionisation energies "in the 584 region have been taken
from Turner et al. [17]": the bands below the He I line (21.22 eV) are Turner, Baker,
Baker & Brundle's values, reprinted; those above are this paper's He II measurements.
Each molecule therefore gets two files, ``<gas>/potts_1974__hei_<gas>.csv`` (below
21.22 eV) and ``<gas>/potts_1974__heii_<gas>.csv`` (above).

Structural fixes:
- A footnote mark on a value ("25.3^{a}") is dropped from the value and kept in
  ``printed``.
- In Table 3 the first orbital of each block is merged into the row's label cell
  ("Ionized orbital or ionization process 1\pi_g", followed by an empty cell); the
  image shows it under the first band. It is moved back to the first column.

Rows are labelled by the band's assignment as printed ("2sigma_g", "1pi_u"); a band
assigned "CI" (a two-electron, configuration-interaction satellite) is labelled
"CI_<energy>" to keep the labels unique. ``printed`` holds the band's four cells.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Potts_1974_ObservationForbiddenTransitions"
TABLE_1 = "tbl_5eb7440f6f19"
TABLE_3 = "tbl_278510f78d78"
HE_I_EV = 21.22  # the He I 584 A line: bands below it are Turner et al.'s (p. 7)
BLOCKS = {
    (TABLE_1, "a N_{2}"): "N2",
    (TABLE_1, "b CO"): "CO",
    (TABLE_3, "a CO_2"): "CO2",
    (TABLE_3, "d N_2O"): "N2O",
}
SKIPPED = {(TABLE_1, "c HCN"), (TABLE_3, "b COS"), (TABLE_3, "c CS_2")}
_FOOTNOTE = re.compile(r"\^\{?[a-f]\}?$")
_MERGED_ORBITAL = re.compile(r"^(Ionized orbital or ionization process)\s+(\S+)$")


def _orbital(cell: str) -> str:
    r"""Return a printed orbital label as plain text: 2\sigma_{g} -> 2sigma_g."""
    text = cell.replace("\\", "").replace("{", "").replace("}", "")
    return text.replace(" ", "")


def _blocks(table_id: str) -> dict[str, dict[str, list[str]]]:
    """Return each molecule's block: the row name (first words) -> its band cells."""
    blocks: dict[str, dict[str, list[str]]] = {}
    current = ""
    for printed_cells in json_table(DIRECTORY, table_id):
        cells = list(printed_cells)
        label = cells[0].strip()
        if all(not cell.strip() for cell in cells[1:]):
            if label != "Ionization":  # Table 3 splits this label over two rows
                current = label
                blocks[current] = {}
            continue
        if match := _MERGED_ORBITAL.match(label):
            label = match.group(1)
            if cells[1].strip():
                msg = f"{table_id} {current}: orbital cell not empty: {cells!r}"
                raise ValueError(msg)
            cells = [label, match.group(2), *cells[2:]]
        key = label.split()[0].lower()  # ionization/potential, relative, ionized, xps
        key = "ionization" if key == "potential" else key
        blocks[current][key] = [cell.strip() for cell in cells[1:]]
    return blocks


def extract() -> dict[str, list[Point]]:
    """Return the band positions of the curated molecules, split at the He I line."""
    files: dict[str, list[Point]] = {}
    for table_id in (TABLE_1, TABLE_3):
        blocks = _blocks(table_id)
        if set(blocks) != {name for t, name in set(BLOCKS) | SKIPPED if t == table_id}:
            msg = f"unexpected blocks in {table_id}: {sorted(blocks)}"
            raise ValueError(msg)
        for name, rows in blocks.items():
            if (table_id, name) in SKIPPED:
                continue
            gas = BLOCKS[table_id, name]
            hei: list[Point] = []
            heii: list[Point] = []
            for column, cell in enumerate(rows["ionization"]):
                if not cell:
                    continue
                value = number(_FOOTNOTE.sub("", cell))
                orbital = _orbital(rows["ionized"][column])
                label = f"CI_{value}" if orbital == "CI" else orbital
                printed = " | ".join(
                    rows[key][column]
                    for key in ("ionization", "relative", "ionized", "xps")
                )
                point = Point(label, orbital, value, "", printed)
                (hei if float(value) < HE_I_EV else heii).append(point)
            files[f"{gas}/potts_1974__hei_{gas}.csv"] = hei
            files[f"{gas}/potts_1974__heii_{gas}.csv"] = heii
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
