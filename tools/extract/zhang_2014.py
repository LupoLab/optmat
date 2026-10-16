# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Zhang, Lan & Lu (2014): the modified barrier-suppression formula, its data.

* Table XI (``tbl_5fd51895d9aa``): a_1, a_2, a_3 of Eq. (8) per target (columns H, He,
  Ne, Ar, Kr, Xe, Ar+), one file per curated target (``tableXI_<target>``), rows
  labelled by the coefficient. The a_1 row is a clean cell per target. STRUCTURAL FIX:
  the a_2 and a_3 rows are merged into one row labelled "aa23" whose cells interleave
  the characters of the two printed values (H: "--00..9006903334" for -0.909 33 and
  -0.060 34). The two values are transcribed from the table image (p. 4) and the script
  checks that each merged cell is exactly an interleaving of the two transcribed
  strings, each in order (and that the a_1 cells are the transcription); a mismatch
  stops the extraction.
* Tables IX and X (``tbl_9167dfe0cbd8``, ``tbl_b692e7df546d``): the half static
  ionisation rates Gamma/2 (au) of Ar and Ar+ against the field amplitude (au), with
  the intensity (W/cm^2) printed beside it, for the GSZ, Muller and Tong-Lin potentials,
  one file per potential (``tableIX_<potential>_Ar``, ``tableX_<potential>_Ar+``).
  STRUCTURAL FIX: both tables are scrambled in the library's text layer (each column's
  fifteen rows are one cell of interleaved characters). They are transcribed below from
  the table images (p. 3), and the script checks each transcribed column two ways:
  it is exactly the multiset of characters of the scrambled text-layer cell, and it is
  exactly the library's OCR reading of the table (``alternatives`` in the JSON cells),
  which the library did not use because it disagreed with the scrambled text layer. A
  mismatch of either stops the extraction.
* Tables I-III (``tbl_8a51c93568de``, ``tbl_312150b826a7``, ``tbl_b5fc44745413``): the
  parameters of the three single-active-electron potentials, Tong-Lin (Eq. (5)), GSZ
  (Eq. (6)) and Muller (Eq. (7)), one file per curated target (``tableI_<target>``,
  ``tableII_<target>``, ``tableIII_<target>``), rows labelled by the parameter. Clean
  CSVs, checked against the images.

H is not curated and is skipped. "Ne^+", "Ar^+", "Ne+" and "Ar+" are the ions Ne+ and
Ar+.
"""

import json
from collections import Counter

from _library import LIBRARY, Point, number, table, write

DIRECTORY = "Zhang_2014_EmpiricalFormulaOver"
TABLE_I = "tbl_8a51c93568de"
TABLE_II = "tbl_312150b826a7"
TABLE_III = "tbl_b5fc44745413"
TABLE_IX = "tbl_9167dfe0cbd8"
TABLE_X = "tbl_b692e7df546d"
TABLE_XI = "tbl_5fd51895d9aa"
TIMES = "\N{MULTIPLICATION SIGN}"
TARGETS = {
    "He": "He",
    "Ne": "Ne",
    "Ar": "Ar",
    "Kr": "Kr",
    "Xe": "Xe",
    "Ne+": "Ne+",
    "Ar+": "Ar+",
    "Ne^+": "Ne+",
    "Ar^+": "Ar+",
}
SKIPPED = {"H"}

# Table XI's a_2 and a_3 rows as printed (the image, p. 4), per column; the a_1 row
# as printed, to check the clean cells.
TABLE_XI_IMAGE: dict[str, tuple[str, str, str]] = {
    "H": ("0.11714", "-0.90933", "-0.06034"),
    "He": ("0.13550", "-0.86210", "0.021562"),
    "Ne": ("0.10061", "-1.04832", "-0.07542"),
    "Ar": ("0.16178", "-1.50441", "0.32127"),
    "Kr": ("0.14640", "-1.36533", "0.02055"),
    "Xe": ("0.21080", "-1.88482", "0.574281"),
    "Ar+": ("0.30441", "-2.70461", "1.40821"),
}

# A rate or intensity as printed: (mantissa, exponent) for "mantissa x 10^exponent", or
# (value, None) for a plain decimal.
Printed = tuple[str, str | None]
# A row of Tables IX-X: the field (au), the intensity, the three potentials' rates.
RateRow = tuple[str, Printed, tuple[Printed, ...]]
# Tables IX and X as printed (the images, p. 3), one row per line: field amplitude
# (au), intensity (W/cm^2), Gamma/2 for the GSZ, Muller and Tong-Lin potentials. A
# number printed "a x 10^k" is written "aek" here; digit groups printed with a thin
# space ("0.053 29") are written without it.
TABLE_IX_TEXT = """
    0.05329  1.000e14  8.7871e-6  8.2616e-6  1.0464e-5
    0.09399  3.111e14  2.6916e-3  2.6573e-3  3.3095e-3
    0.12177  5.222e14  0.01110  0.013229  0.01110
    0.14430  7.333e14  0.02166  0.021782  0.02503
    0.16376  9.444e14  0.03208  0.032306  0.03606
    0.18114  1.156e15  0.04143  0.041682  0.04543
    0.19699  1.367e15  0.04941  0.049612  0.05305
    0.21166  1.578e15  0.05607  0.056152  0.05917
    0.22537  1.789e15  0.06158  0.061511  0.06411
    0.23830  2.000e15  0.06617  0.065925  0.06816
    0.24993  2.200e15  0.06984  0.069433  0.07137
    0.26776  2.525e15  0.07476  0.074113  0.07568
    0.28447  2.850e15  0.07880  0.077925  0.07924
    0.30025  3.175e15  0.08220  0.081178  0.08230
    0.31524  3.500e15  0.08523  0.084068  0.08505
"""
TABLE_X_TEXT = """
    0.07536  2.000e14  6.0344e-8  7.2791e-8  6.6251e-8
    0.12306  5.333e14  2.8655e-4  3.8221e-4  3.7933e-4
    0.15687  8.667e14  3.5719e-3  4.4861e-3  4.5811e-3
    0.18459  1.200e15  0.01136  0.01377  0.01413
    0.20865  1.533e15  0.02180  0.02584  0.02646
    0.23022  1.867e15  0.03291  0.03828  0.03905
    0.24993  2.200e15  0.04338  0.04956  0.05034
    0.26820  2.533e15  0.05254  0.05896  0.05969
    0.28530  2.867e15  0.06022  0.06645  0.06710
    0.30143  3.200e15  0.06651  0.07231  0.07290
    0.31524  3.500e15  0.07117  0.07650  0.07705
    0.33701  4.000e15  0.07733  0.08192  0.08242
    0.35745  4.500e15  0.08209  0.08607  0.08654
    0.37678  5.000e15  0.08594  0.08950  0.08991
    0.39517  5.500e15  0.08923  0.09253  0.09287
"""


def _printed(text: str) -> Printed:
    """Return "aek" as ("a", "k"), a plain decimal as (it, None)."""
    mantissa, _, exponent = text.partition("e")
    return (mantissa, exponent or None)


def _image(text: str) -> list[RateRow]:
    """Return the rows of a transcribed rate table."""
    rows: list[RateRow] = []
    for line in text.strip().splitlines():
        field, intensity, *rates = line.split()
        rows.append((field, _printed(intensity), tuple(_printed(r) for r in rates)))
    return rows


TABLE_IX_IMAGE = _image(TABLE_IX_TEXT)
TABLE_X_IMAGE = _image(TABLE_X_TEXT)
POTENTIALS = ("gsz", "muller", "tonglin")
RATE_HEADERS = ("GSZ", "Muller", "Tong-Lin")


def _interleaves(merged: str, first: str, second: str) -> bool:
    """Return whether ``merged`` interleaves ``first`` and ``second``, each in order."""
    if len(merged) != len(first) + len(second):
        return False
    # reachable[j]: merged[: i + j] interleaves first[:i] and second[:j]
    reachable = [True] * (len(second) + 1)
    for j in range(1, len(second) + 1):
        reachable[j] = reachable[j - 1] and second[j - 1] == merged[j - 1]
    for i in range(1, len(first) + 1):
        reachable[0] = reachable[0] and first[i - 1] == merged[i - 1]
        for j in range(1, len(second) + 1):
            k = i + j - 1
            reachable[j] = (reachable[j] and first[i - 1] == merged[k]) or (
                reachable[j - 1] and second[j - 1] == merged[k]
            )
    return reachable[len(second)]


def _table_xi() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_XI)
    header = [cell.strip() for cell in grid[0]]
    if header != ["Parameter", *TABLE_XI_IMAGE] or [r[0] for r in grid[1:]] != [
        "a1",
        "aa23",
    ]:
        msg = f"unexpected Table XI layout {grid!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for column, heading in enumerate(header[1:], start=1):
        a_1, a_2, a_3 = TABLE_XI_IMAGE[heading]
        if number(grid[1][column]) != a_1:
            msg = f"Table XI {heading} a_1: {grid[1][column]!r} is not {a_1!r}"
            raise ValueError(msg)
        merged = grid[2][column].strip()
        if not _interleaves(merged, a_2, a_3):
            msg = f"Table XI {heading}: {merged!r} does not interleave {a_2!r}, {a_3!r}"
            raise ValueError(msg)
        if heading in SKIPPED:
            continue
        target = TARGETS[heading]
        printed = f"{heading} | a1 {grid[1][column]} | a2, a3 merged {merged}"
        files[f"{target}/zhang_2014__tableXI_{target}.csv"] = [
            Point("a_1", "a_1", a_1, "", printed),
            Point("a_2", "a_2", a_2, "", f"{printed} (from the image)"),
            Point("a_3", "a_3", a_3, "", f"{printed} (from the image)"),
        ]
    return files


def _as_text(item: Printed) -> str:
    """Return a printed number as the characters the table shows ("a", times, "10k")."""
    mantissa, exponent = item
    return mantissa if exponent is None else f"{mantissa}{TIMES}10{exponent}"


def _as_latex(item: Printed) -> str:
    r"""Return a printed number as the library's OCR writes it ("a \times 10^{k}")."""
    mantissa, exponent = item
    return mantissa if exponent is None else rf"{mantissa} \times 10^{{{exponent}}}"


def _as_number(item: Printed) -> str:
    """Return a printed number as a plain number string ("ae-k")."""
    mantissa, exponent = item
    return mantissa if exponent is None else f"{mantissa}e{exponent}"


def _ocr_rows(table_id: str) -> list[list[str]]:
    """Return the rows of the library's OCR alternative of a table, cells normalised."""
    path = LIBRARY / DIRECTORY / "tables" / f"{table_id}.json"
    record = json.loads(path.read_text(encoding="utf-8"))
    (alternative,) = [a for a in record["alternatives"] if a.get("source") == "ocr"]
    rows = []
    for row in alternative["html"].split("<tr>")[1:]:
        cells = [cell.split("</td>")[0] for cell in row.split("<td>")[1:]]
        rows.append(
            [
                cell.replace("<eq>", "").replace("</eq>", "").replace(r"\,", "")
                for cell in cells
            ]
        )
    return rows


def _rates(
    table_id: str,
    image: list[RateRow],
    material: str,
    stem: str,
) -> dict[str, list[Point]]:
    grid = table(DIRECTORY, table_id)
    if len(grid) != 2 or grid[0][0] != "Field amplitude (intensity)":
        msg = f"unexpected {stem} layout {grid[0]!r}"
        raise ValueError(msg)
    if [cell.split(": ")[-1] for cell in grid[0][1:]] != list(RATE_HEADERS):
        msg = f"unexpected {stem} header {grid[0]!r}"
        raise ValueError(msg)
    # Check 1: each scrambled text-layer column is the multiset of the transcription's
    # characters (spaces dropped).
    transcribed = [
        "".join(f"{row[0]}({_as_text(row[1])})" for row in image),
        *("".join(_as_text(row[2][k]) for row in image) for k in range(3)),
    ]
    for column, (cell, text) in enumerate(zip(grid[1], transcribed, strict=True)):
        if Counter(cell.replace(" ", "")) != Counter(text):
            msg = f"{stem} column {column}: {cell!r} does not scramble {text!r}"
            raise ValueError(msg)
    # Check 2: the transcription is the library's OCR reading, cell for cell.
    ocr = _ocr_rows(table_id)
    expected = [
        [
            f"{row[0]} ({_as_latex(row[1])})",
            *(_as_latex(rate) for rate in row[2]),
        ]
        for row in image
    ]
    if ocr[1:] != expected:
        msg = f"{stem}: the transcription differs from the OCR reading {ocr[1:]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for k, potential in enumerate(POTENTIALS):
        points = []
        for row in image:
            printed = " | ".join(
                [
                    f"{row[0]} ({_as_text(row[1])})",
                    *(_as_text(rate) for rate in row[2]),
                ]
            )
            points.append(
                Point(
                    row[0],
                    row[0],
                    _as_number(row[2][k]),
                    "",
                    f"{printed} (from the image)",
                    {"intensity": _as_number(row[1])},
                )
            )
        files[f"{material}/zhang_2014__{stem}_{potential}_{material}.csv"] = points
    return files


def _parameters(
    table_id: str, stem: str, names: tuple[str, ...]
) -> dict[str, list[Point]]:
    """Return one file per curated target of a potential-parameter table."""
    grid = table(DIRECTORY, table_id)
    if tuple(cell.strip() for cell in grid[0][1:]) != names:
        msg = f"unexpected {stem} header {grid[0]!r}"
        raise ValueError(msg)
    labels = tuple(name.replace("Zc", "Z_c") for name in names)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        heading = cells[0].strip()
        if heading in SKIPPED:
            continue
        target = TARGETS[heading]
        printed = " | ".join(cells)
        files[f"{target}/zhang_2014__{stem}_{target}.csv"] = [
            Point(label, label, number(cell), "", printed)
            for label, cell in zip(labels, cells[1:], strict=True)
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the files of Tables I-III and IX-XI."""
    return (
        _table_xi()
        | _rates(TABLE_IX, TABLE_IX_IMAGE, "Ar", "tableIX")
        | _rates(TABLE_X, TABLE_X_IMAGE, "Ar+", "tableX")
        | _parameters(
            TABLE_I, "tableI", ("Z_c", "a_1", "a_2", "a_3", "a_4", "a_5", "a_6")
        )
        | _parameters(TABLE_II, "tableII", ("Z", "Zc", "H", "d"))
        | _parameters(TABLE_III, "tableIII", ("Z", "Zc", "A", "B", "C"))
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
