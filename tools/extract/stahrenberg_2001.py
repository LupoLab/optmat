# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Stahrenberg et al. (2001): n and k of Ag and Cu, as the database holds them.

The paper prints its dielectric functions in figures only. The numbers were given by N.
Esser to the refractiveindex.info database ("Numerical data kindly provided by Prof. Dr.
Norbert Esser"), whose files ``main/Ag/nk/Stahrenberg.yml`` and
``main/Cu/nk/Stahrenberg.yml`` are kept in ``data/external/refractiveindex_info/`` (see
its ``PROVENANCE.md``); they were chosen as this source's numbers (2026-09-28). Each
file has one ``tabulated nk`` block: rows of wavelength (um), n and k, separated by
spaces, Ag 0.12782-0.49594 um (361 rows) and Cu 0.13776-0.49594 um (326 rows). The
script checks the files' REFERENCES (the paper and the note) and block type, and keeps
the numbers as the files print them; the row id is the wavelength.
"""

from _library import DATA, Point, number, write

DATABASE = DATA / "external" / "refractiveindex_info" / "main"
ROWS = {"Ag": 361, "Cu": 326}
PROVIDED = "(Numerical data kindly provided by Prof. Dr. Norbert Esser)"
DOI = "https://doi.org/10.1103/PhysRevB.64.115111"


def _rows(metal: str) -> list[list[str]]:
    """Return the file's rows of wavelength, n and k, after checking its preamble."""
    name = f"{metal}/nk/Stahrenberg.yml"
    lines = (DATABASE / name).read_text(encoding="utf-8").splitlines()
    text = "\n".join(lines)
    if PROVIDED not in text or DOI not in text or text.count("type:") != 1:
        msg = f"{name}: not the expected single Stahrenberg nk block"
        raise ValueError(msg)
    start = lines.index("  - type: tabulated nk") + 1
    if lines[start].strip() != "data: |":
        msg = f"{name}: unexpected block {lines[start - 1 : start + 1]}"
        raise ValueError(msg)
    rows = [line.split() for line in lines[start + 1 :] if line.strip()]
    if len(rows) != ROWS[metal] or any(len(cells) != 3 for cells in rows):
        msg = f"{name}: {len(rows)} rows, not {ROWS[metal]} of three numbers"
        raise ValueError(msg)
    return rows


def extract() -> dict[str, list[Point]]:
    """Return n and k of each metal against the wavelength."""
    files: dict[str, list[Point]] = {}
    for metal in ROWS:
        rows = _rows(metal)
        for part, column in (("n", 1), ("k", 2)):
            files[f"{metal}/stahrenberg_2001__{part}_{metal}.csv"] = [
                Point(
                    number(cells[0]),
                    number(cells[0]),
                    number(cells[column]),
                    "",
                    " ".join(cells),
                )
                for cells in rows
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
