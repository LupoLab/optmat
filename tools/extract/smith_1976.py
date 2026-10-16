# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Smith, Huber & Parkinson (1976): n - 1 of He, H2, O2, CO and Kr, 167-298 nm.

All tables give ``(n - 1) x 10^6`` at vacuum wavelengths (nm), reduced to 273.15 K and
760 Torr. The table images were compared with the extracted cells.

- Table IV, helium (``tbl_e2f6e80b5c2b``, cells file only): the "Identification"
  header spans two cells (element and multiplet), so "This experiment" is column 3,
  Victor's theory column 4 and "Other measurements" column 6, whose footnote letter
  names the source (d Abjean et al. 1970, e Huber and Tondello 1974, f Bideau-Mehu et
  al., unpublished); each goes to its own dataset. The theory column covers every row,
  including the other authors' wavelengths.
- Table V, hydrogen (``tbl_5b5868c91642``): this experiment, the paper's two-term
  Sellmeier interpolation (a check), Ford's theory, Koch (1913) and Kirn (1921). The
  Koch cells marked ``h``, "Extrapolated values", are not extracted.
- Table VI, oxygen (``tbl_1bad6e7f544c``): the ``b`` marks on 191.85 and 184.30 nm
  ("these wavelengths ... do not correspond to spectral lines", +-0.1 nm) are stripped
  from the wavelengths; those two rows carry a larger uncertainty in the record.
- Table VIII, krypton (``tbl_8156b916136f``): this experiment, the interpolation, the
  Bideau-Mehu et al. values (d) and the Kronjager/Koch values (f) as listed by Leonard.
- Table VII, carbon monoxide (``tbl_a6e86bb4461b``): the "This expt." column holds
  both measured values and, marked ``d``, values "computed from three-term Lagrange fit
  to our data" at Koch's wavelengths (image checked: the ``d`` marks are on exactly
  the ten rows with a "Prev. expt." value). The two go to separate datasets; column 2
  is Koch (1914) as listed in the International Critical Tables. The Delta (%) column
  is derived and not extracted.
"""

import json
import re

from _library import LIBRARY, Point, number, table, write

DIRECTORY = "Smith_1976_RefractivitiesH2He"
HELIUM = "tbl_e2f6e80b5c2b"
HYDROGEN = "tbl_5b5868c91642"
OXYGEN = "tbl_1bad6e7f544c"
KRYPTON = "tbl_8156b916136f"
CARBON_MONOXIDE = "tbl_a6e86bb4461b"

_MARK = re.compile(r"\^\{?([a-h])\}?$")
# Table VII's cells print their marks as ``^{\text{d}}``.
_TEXT_MARK = re.compile(r"\^\{\\text\{d\}\}$")

# Helium "Other measurements" footnote letter -> dataset.
HELIUM_OTHERS = {
    "d": "abjean_1970_He",
    "e": "huber_tondello_1974_He",
    "f": "bideau_mehu_unpublished_He",
}
# Column index -> dataset, for the hydrogen and krypton tables.
HYDROGEN_COLUMNS = {
    1: "experiment_H2",
    2: "interpolated_H2",
    4: "ford_theory_H2",
    6: "koch_1913_H2",
    8: "kirn_1921_H2",
}
KRYPTON_COLUMNS = {
    1: "experiment_Kr",
    2: "interpolated_Kr",
    4: "bideau_mehu_1981_Kr",
    6: "kronjager_koch_Kr",
}


def _cells_grid(directory: str, table_id: str) -> list[list[str]]:
    """Return a table's rows from its cells file, for tables without a CSV."""
    path = LIBRARY / directory / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))["cells"]
    n_rows = 1 + max(c["row"] for c in cells)
    n_columns = max(c["column"] + c["column_span"] for c in cells)
    grid = [[""] * n_columns for _ in range(n_rows)]
    for cell in cells:
        grid[cell["row"]][cell["column"]] = cell["text"]
    return grid


def _split_mark(cell: str) -> tuple[str, str]:
    """Return a cell's text without its footnote letter, and the letter (or "")."""
    text = cell.strip()
    match = _MARK.search(text)
    if match is None:
        return text, ""
    return text[: match.start()].strip(), match.group(1)


def _path(gas: str, dataset: str) -> str:
    return f"{gas}/smith_1976__{dataset}.csv"


def _helium() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in _cells_grid(DIRECTORY, HELIUM)[1:]:
        wavelength, printed = cells[0].strip(), " | ".join(cells)
        if cells[3].strip():
            files.setdefault(_path("He", "experiment_He"), []).append(
                Point(wavelength, wavelength, number(cells[3]), "", printed)
            )
        files.setdefault(_path("He", "victor_theory_He"), []).append(
            Point(wavelength, wavelength, number(cells[4]), "", printed)
        )
        if cells[6].strip():
            value, mark = _split_mark(cells[6])
            files.setdefault(_path("He", HELIUM_OTHERS[mark]), []).append(
                Point(wavelength, wavelength, number(value), "", printed)
            )
    return files


def _by_column(
    gas: str, table_id: str, columns: dict[int, str], skip_marks: dict[int, str]
) -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, table_id)[1:]:
        wavelength, printed = cells[0].strip(), " | ".join(cells)
        for column, dataset in columns.items():
            value, mark = _split_mark(cells[column])
            if not value or (mark and mark == skip_marks.get(column)):
                continue
            files.setdefault(_path(gas, dataset), []).append(
                Point(wavelength, wavelength, number(value), "", printed)
            )
    return files


def _oxygen() -> list[Point]:
    points = []
    for cells in table(DIRECTORY, OXYGEN)[1:]:
        wavelength, _ = _split_mark(cells[0])
        points.append(
            Point(wavelength, wavelength, number(cells[1]), "", " | ".join(cells))
        )
    return points


def _carbon_monoxide() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, CARBON_MONOXIDE)[1:]:
        wavelength, printed = cells[0].strip(), " | ".join(cells)
        value = cells[1].strip()
        # Cells marked d are the paper's Lagrange interpolation, not measurements.
        marked = _TEXT_MARK.search(value)
        if marked:
            value = value[: marked.start()]
        dataset = "interpolated_CO" if marked else "experiment_CO"
        files.setdefault(_path("CO", dataset), []).append(
            Point(wavelength, wavelength, number(value), "", printed)
        )
        if cells[2].strip():
            files.setdefault(_path("CO", "koch_1914_CO"), []).append(
                Point(wavelength, wavelength, number(cells[2]), "", printed)
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the helium, hydrogen, oxygen, carbon-monoxide and krypton datasets."""
    files = _helium()
    # Koch's hydrogen cells marked h are extrapolations.
    files |= _by_column("H2", HYDROGEN, HYDROGEN_COLUMNS, {6: "h"})
    files[_path("O2", "experiment_O2")] = _oxygen()
    files |= _by_column("Kr", KRYPTON, KRYPTON_COLUMNS, {})
    files |= _carbon_monoxide()
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
