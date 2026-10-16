# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Matsui et al. (2003): O2 Schumann-Runge (v', 0) band oscillator strengths.

Table 5 (``tbl_0c3259b06159``, CSV; image checked) gives, for the bands (12,0)-(17,0),
the total and extended integrated cross-sections (10^-17 cm^2 cm^-1; not extracted, see
the record) and the band oscillator strengths (10^-5) of this work and of five earlier
works: LBC (Lewis et al. 1986), YFEP (Yoshino et al. 1983), GGMBL (Gies et al. 1981),
HCMK (Huebner et al. 1975) and Bethke (1959). Structural fixes, from the image:

- the second header line (Total, Ext., Present, LBC, YFEP, GGMBL, HCMK, Bethke, with
  their note marks) is merged into the (12,0) cells ("Totala2.04", "Present2.38"): the
  number is taken from the end of each cell;
- the LBC and YFEP columns are merged into one cell column: in the (12,0) row the cell
  "LBCc YFEPd2.44 2.43" holds LBC 2.44 and YFEP 2.43; the other rows have no YFEP
  value (the image's YFEP column is empty below (12,0)), so their single number is LBC.

The rows are the bands "(v',0)", the abscissa the index v'. The mantissas are as
printed; the 10^-5 factor is the dataset scale.
"""

import re

from _library import Point, table, write

DIRECTORY = "Matsui_2003_HighResolutionAbsorption"
TABLE_5 = "tbl_0c3259b06159"
PREFIX = "O2/matsui_2003__"
COLUMNS = ("present", "lbc", "yfep", "ggmbl", "hcmk", "bethke")

_NUMBER = re.compile(r"\d+\.\d+")


def extract() -> dict[str, list[Point]]:
    """Return this work's and the five reproduced band oscillator strengths."""
    files: dict[str, list[Point]] = {c: [] for c in COLUMNS}
    for cells in table(DIRECTORY, TABLE_5)[1:]:
        band, _, _, present, merged, ggmbl, hcmk, bethke = (c.strip() for c in cells)
        match = re.fullmatch(r"\((\d+),0\)", band)
        if match is None:
            msg = f"unexpected band {band!r}"
            raise ValueError(msg)
        upper = match.group(1)
        lbc_yfep = _NUMBER.findall(merged)
        lbc = lbc_yfep[0] if lbc_yfep else ""
        yfep = lbc_yfep[1] if len(lbc_yfep) > 1 else ""
        values = {
            "present": present,
            "lbc": lbc,
            "yfep": yfep,
            "ggmbl": ggmbl,
            "hcmk": hcmk,
            "bethke": bethke,
        }
        printed = " | ".join(c.strip() for c in cells)
        for column, cell in values.items():
            numbers = _NUMBER.findall(cell)
            if not numbers:
                continue
            files[column].append(Point(band, upper, numbers[-1], "", printed))
    return {f"{PREFIX}{column}.csv": points for column, points in files.items()}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
