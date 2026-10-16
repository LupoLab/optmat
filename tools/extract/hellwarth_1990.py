# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hellwarth, Pennington & Henesian (1990): dc Kerr, sigma and c coefficients.

All values are in 1e-19 esu at 20 C and 760 Torr.

- Table I (p. 2770, JSON cells): (n_par - n_perp) / (2 pi E0^2) per gas and wavelength,
  in three columns: "Normalized to He" (a ratio), "Absolute" and "Average". A letter
  after each value names its origin (footnotes a-g); it decides the data file. A
  parenthesis holds the "limits of error of last figures given by source": 251(36) is
  251 +- 36 and 1.99(13) is 1.99 +- 0.13, converted here to the printed unit.
- Table II (p. 2771, CSV): the electronic parameter sigma, errors in parentheses as
  absolute numbers (36.5(1.4) is 36.5 +- 1.4; He 1.76(0.05) is the ~3 % of Table I's
  note g).
- Tables V and VII (pp. 2772, 2774, JSON): sigma + B0 extrapolated to the
  long-wavelength limit from each wavelength's Kerr value, with formula A and B; the
  parentheses there are relative weights of the average (the footnotes say so), not
  errors: they are kept only in the ``printed`` column.
- Tables VI and VIII (pp. 2773, 2775, JSON): the degenerate four-wave-mixing c
  coefficients c1122, c1111 (formula A) and c1122, c1212, c1111 (formula B) at 261-1053
  nm and the long-wavelength limit (LWL, written ``inf``). Air rows are skipped.

The cells of Tables I, VI and VIII match the images. Tables V and VII need structural
fixes, commented where they are made.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Hellwarth_1990_IndicesGoverningOptical"
KERR = "tbl_223eea4ff99e"
SIGMA = "tbl_0c2959460533"
LWL_A = "tbl_2ee7689e92b6"
DFWM_A = "tbl_03c27c55bc84"
LWL_B = "tbl_3b8085ea2a30"
DFWM_B = "tbl_757adefa584f"

GASES = {"N_{2}": "N2", "O_{2}": "O2", "Ar": "Ar", "H_{2}": "H2", "He": "He"}
GASES_PLAIN = {"N_2": "N2", "O_2": "O2", "Ar": "Ar", "H_2": "H2", "He": "He"}
CELL = re.compile(r"^([\d.]+)(?:\(([\d.]+)\))?(?:\^?\{?([a-g])\}?)?$")
RATIO = {"a": "kerr_ratio_carusotto", "e": "kerr_ratio_buckingham_dunmur"}
ABSOLUTE = {
    "b": "kerr_he_scaled",
    "c": "kerr_carusotto",
    "d": "kerr_buckingham_1970",
    "f": "kerr_buckingham_orr",
    "g": "kerr_bishop_lam",
}


def _last_figures(value: str, error: str) -> str:
    """Return an error given in units of the value's last figure, in its unit."""
    decimals = len(value.split(".")[1]) if "." in value else 0
    return f"{int(error) / 10**decimals:.{decimals}f}"


def _kerr() -> dict[str, list[Point]]:
    """Table I, split by column and footnote letter."""
    files: dict[str, list[Point]] = {}
    gas = wavelength = ""
    for cells in json_table(DIRECTORY, KERR)[2:]:
        gas = GASES.get(cells[0], gas)
        wavelength = {"": wavelength, "∞": "inf"}.get(cells[1], cells[1])
        printed = " | ".join(cells)
        for column, names in ((2, RATIO), (3, ABSOLUTE), (4, None)):
            if not cells[column] or (gas == "He" and column == 2):
                continue  # He normalized to itself (1) is not a datum
            match = CELL.match(cells[column])
            if match is None:
                raise ValueError(cells[column])
            value, error, letter = match.groups()
            sigma = _last_figures(value, error) if error else ""
            name = "kerr_average" if names is None else names[letter]
            path = f"{gas}/hellwarth_1990__{name}_{gas}.csv"
            files.setdefault(path, []).append(
                Point(wavelength, wavelength, number(value), sigma, printed)
            )
    return files


def _sigma() -> dict[str, list[Point]]:
    """Table II: sigma per gas, in the long-wavelength limit."""
    files = {}
    for cells in table(DIRECTORY, SIGMA)[1:]:
        gas = GASES[cells[0]]
        match = re.match(r"^([\d.]+)\(([\d.]+)\)", cells[1])
        if match is None:
            raise ValueError(cells[1])
        value, error = match.groups()
        path = f"{gas}/hellwarth_1990__sigma_{gas}.csv"
        files[path] = [Point("LWL", "inf", value, error, " | ".join(cells))]
    return files


def _weighted(cell: str) -> str:
    """Return the value of a Table V/VII cell, dropping the weight in parentheses."""
    return number(re.sub(r"\(\d+\)$", "", cell))


def _lwl_a() -> dict[str, list[Point]]:
    """Table V: sigma + B0 from formula A, per data wavelength, and the averages.

    Structural fix: the table image prints each gas's average on its own line below the
    gas's rows (N2 242, O2 517, Ar 44.4 and "(47.4)^f", H2 55.8, He 1.76^g); the library
    cells put them in the "average" column of the wrong rows (242 beside O2, 44.4 beside
    H2, ...). They are reassigned in their printed order. The He rows (Carusotto's and
    Buckingham and Dunmur's He, divided by F^A, and the theory value) are not extracted.
    """
    rows = json_table(DIRECTORY, LWL_A)[1:]
    averages = [cells[4] for cells in rows if cells[4]]
    if averages != ["242", "517", "44.4(47.4)f", "55.8", "1.76^g"]:
        raise ValueError(averages)
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in rows:
        gas = GASES_PLAIN.get(cells[0], gas)
        if gas == "He" or not cells[3]:
            continue
        path = f"{gas}/hellwarth_1990__kerr_lwl_A_{gas}.csv"
        files.setdefault(path, []).append(
            Point(cells[1], "inf", _weighted(cells[3]), "", " | ".join(cells))
        )
    for gas, cell in zip(("N2", "O2", "Ar", "H2"), averages, strict=False):
        path = f"{gas}/hellwarth_1990__kerr_lwl_A_{gas}.csv"
        printed = f"{gas} average (printed below the gas's rows) | {cell}"
        if gas == "Ar":  # "44.4 (47.4)^f": the average, and the value used (from sigma)
            files[path].append(Point("average", "inf", "44.4", "", printed))
            files[path].append(Point("used", "inf", "47.4", "", printed))
        else:
            files[path].append(Point("average", "inf", number(cell), "", printed))
    return files


def _lwl_b() -> dict[str, list[Point]]:
    """Table VII: sigma + B0 from formula B, per data wavelength, and the averages.

    Structural fix: the image prints H2's second value, 49.7(3), on a line with no
    wavelength and no F factor; it is the 633-nm Kerr value of Table I (52.5, footnote
    f) divided by F_1221^B(633 nm) = 1.057 (omega_0 = 100e3 cm^-1), so its row is 633.
    The He rows are dispersed He Kerr values (1.76 x F, footnote c), not limits,
    and are not extracted.
    """
    files: dict[str, list[Point]] = {}
    gas = wavelength = ""
    for cells in json_table(DIRECTORY, LWL_B)[2:]:
        gas = GASES_PLAIN.get(cells[0], gas)
        if gas == "He":
            continue
        wavelength = cells[1] or ("633" if gas == "H2" else wavelength)
        path = f"{gas}/hellwarth_1990__kerr_lwl_B_{gas}.csv"
        printed = " | ".join(cells)
        if cells[3]:
            files.setdefault(path, []).append(
                Point(wavelength, "inf", _weighted(cells[3]), "", printed)
            )
        if cells[4] == "(47.4)b":  # Ar: the value used (from sigma), footnote b
            files[path].append(Point("used", "inf", "47.4", "", printed))
        elif cells[4]:
            files[path].append(Point("average", "inf", number(cells[4]), "", printed))
    return files


def _dfwm(table_id: str, columns: dict[int, str], formula: str) -> dict:
    """Tables VI and VIII: c coefficients per gas and wavelength (air skipped)."""
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in json_table(DIRECTORY, table_id)[1:]:
        if cells[0]:
            gas = GASES_PLAIN.get(cells[0], "air")
            continue
        if gas == "air":
            continue  # air is not a curated species
        row = cells[1]
        wavelength = "inf" if row == "LWL" else row
        for column, component in columns.items():
            path = f"{gas}/hellwarth_1990__dfwm_{formula}_{component}_{gas}.csv"
            files.setdefault(path, []).append(
                Point(row, wavelength, number(cells[column]), "", " | ".join(cells))
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return every data file of the source."""
    return {
        **_kerr(),
        **_sigma(),
        **_lwl_a(),
        **_lwl_b(),
        **_dfwm(DFWM_A, {3: "c1122", 4: "c1111"}, "A"),
        **_dfwm(DFWM_B, {2: "c1122", 3: "c1212", 4: "c1111"}, "B"),
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
