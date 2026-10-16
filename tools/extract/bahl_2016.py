# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Bahl, Wright & Kolesik (2016): static SAE ionisation rates of Ne, Ar, Kr, Xe.

The paper's Supplemental Material (ref. [22], p. 5: "the Supplemental Material [22]
includes these data sets") is attached to the paper in the library byte for byte, with
its origin: per atom an ``<ATOM>_imE_nlp.dat`` and an ``<ATOM>_pa_corr.dat``. Each
``_imE_nlp`` file has three whitespace-separated columns and no header: the field F
(au), the imaginary part of the energy of the metastable ground state, Im E_0(F) (au,
negative), and the nonlinear dipole moment D_NL(F) (Eq. (31), au). The order is the
paper's list of tabulated quantities ("D_NL(F), E_0(F), and M(F,0)", p. 5) read with
the file's name; the reading is checked: Im E is never positive and |Im E| rises with
F, and D_NL/F^3 is constant at low field (Ne 11.3, Ar 159, Kr 355, Xe 1275 au at 0.001
to 0.003 au), as a third-order response must be.

One file per gas, ``rate_<gas>``: the row id and ``x`` are F as printed; the value is
Im E_0 as printed (the dataset's scale, -2, makes it the rate Gamma = -2 Im E_0); the
``printed`` cell keeps the file's whole line, the nonlinear dipole included.

He is not extracted. Its field column is printed to three decimals while its grid is
(k - 1/2) x 0.001 au (D_NL/F^3 is constant only with that grid: 4.67 to 4.75 au over
rows 2 to 30, against 0.6 to 4.5 with the printed fields), so 20 printed fields occur
twice and every printed field is off by 0.0005 au: the rates cannot be placed on the
field axis from the file (the source record says more).

The ``_pa_corr`` files (the post-adiabatic correction M(F,0) of Eq. (26)) are not
extracted.
"""

from _library import Point, attachment, number, write

PAPER = "Bahl_2016_NonlinearOpticalResponse"
# Gas -> the attached file and its number of rows (Ne lacks 0.061 and 0.091 au, Xe
# 0.017 and 0.027 au).
FILES = {
    "Ne": ("NE_imE_nlp.dat", 149),
    "Ar": ("AR_imE_nlp.dat", 151),
    "Kr": ("KR_imE_nlp.dat", 131),
    "Xe": ("XE_imE_nlp.dat", 99),
}
COLUMNS = 3


def _rows(name: str, count: int) -> list[list[str]]:
    """Return a file's rows as cells, after checking their shape and order."""
    lines = attachment(PAPER, name).read_text(encoding="utf-8").splitlines()
    rows = [line.split() for line in lines if line.strip()]
    if len(rows) != count or {len(cells) for cells in rows} != {COLUMNS}:
        msg = f"{name}: {len(rows)} rows of {sorted({len(c) for c in rows})} cells"
        raise ValueError(msg)
    fields = [float(cells[0]) for cells in rows]
    if fields != sorted(set(fields)):
        msg = f"{name}: the fields are not strictly increasing"
        raise ValueError(msg)
    if any(float(cells[1]) > 0 for cells in rows):
        msg = f"{name}: a positive Im E"
        raise ValueError(msg)
    return rows


def extract() -> dict[str, list[Point]]:
    """Return Im E_0 against F for each gas."""
    files: dict[str, list[Point]] = {}
    for gas, (name, count) in FILES.items():
        files[f"{gas}/bahl_2016__rate_{gas}.csv"] = [
            Point(
                number(cells[0]),
                number(cells[0]),
                number(cells[1]),
                "",
                " | ".join(cells),
            )
            for cells in _rows(name, count)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
