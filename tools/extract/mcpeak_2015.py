# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract McPeak et al. (2015): epsilon of the template-stripped Ag, Au, Cu, Al films.

The paper and the library's Supporting Information PDF print the optical constants in
figures only; the "digital files containing the actual data points" of the Supporting
Information are attached to the paper in the library, with their origin:

- ``ph5004237_si_002.txt``, "Optical Data for template stripped Silver, Gold and
  Copper": wavelength (nm), then e1 and e2 of Ag, Au and Cu, 300-1700 nm in 10 nm
  steps (141 rows);
- ``ph5004237_si_003.txt``, "Optical Data for template stripped Aluminum": wavelength
  (nm), e1 and e2, 200-1700 nm in 5 nm steps (287 rows: 1360-1425 nm are absent, the
  file jumps from 1355 to 1430 nm).

Each file opens with six ``#`` comment lines (the laboratory, the date, the paper's
title, the contents and the column names, checked here) and has tab-separated rows;
the first file ends with a line of empty tab-separated cells, which is skipped. The
numbers are kept as the files print them; the row id is the wavelength.
"""

from _library import Point, attachment, number, write

PAPER = "McPeak_2015_PlasmonicFilmsCan"

# File -> (the column-name comment, the metals in column order, the row count).
CONTENTS = {
    "ph5004237_si_002.txt": (
        "# Wavelength(nm) Ag(e1) Ag(e2) Au(e1) Au(e2) Cu(e1) Cu(e2)",
        ("Ag", "Au", "Cu"),
        141,
    ),
    "ph5004237_si_003.txt": ("# Wavelength (nm) e1 e2", ("Al",), 287),
}


def _rows(name: str, header: str, count: int) -> list[list[str]]:
    """Return a file's data rows as lists of cells, after checking its header."""
    lines = attachment(PAPER, name).read_text(encoding="utf-8").splitlines()
    comments = [line for line in lines if line.startswith("#")]
    if len(comments) != 6 or comments[-1].strip() != header:
        msg = f"{name}: unexpected header {comments}"
        raise ValueError(msg)
    rows = [
        line.split("\t")
        for line in lines
        if not line.startswith("#") and line.replace("\t", "").strip()
    ]
    if len(rows) != count:
        msg = f"{name}: {len(rows)} rows, not {count}"
        raise ValueError(msg)
    return rows


def extract() -> dict[str, list[Point]]:
    """Return epsilon_1 and epsilon_2 of each film against the wavelength."""
    files: dict[str, list[Point]] = {}
    for name, (header, metals, count) in CONTENTS.items():
        rows = _rows(name, header, count)
        for index, metal in enumerate(metals):
            for part, column in (("eps1", 1 + 2 * index), ("eps2", 2 + 2 * index)):
                files[f"{metal}/mcpeak_2015__{part}_{metal}.csv"] = [
                    Point(
                        number(cells[0]),
                        number(cells[0]),
                        number(cells[column]),
                        "",
                        " | ".join(cell.strip() for cell in cells),
                    )
                    for cells in rows
                ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
