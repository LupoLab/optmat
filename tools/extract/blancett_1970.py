# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Blancett, Hall & Canfield (1970): Burnett isotherms of He-Ar mixtures.

Table I (pp. 82-87, six library tables, one per page, cells files only) prints, for
each composition, two Burnett runs at 323.15, 273.15 and 223.15 K: the observed
pressure P (atm), Z = Pv/RT and sigma(Z) x 10^5, one column triple per temperature.
Table II (p. 88, ``tbl_c09ea288c262``) prints the virial coefficients B_k
((cc/g mol)^(k-1)) of each isotherm with their standard deviations.

Extracted here:

- the four mixtures' Z up to 150 atm, one file per composition and temperature
  (``Ar/blancett_1970__Z_xAr_<x>_<T>K.csv``), the row id the printed pressure, and the
  run in the printed cell; sigma(Z) is written in e-notation of the printed digits
  (``14`` becomes ``14e-5``), the raw cell kept in ``printed``. The pure-helium and
  pure-argon isotherms are not extracted (optmat takes the pure gases from their
  reference equations of state);
- the second virial coefficient B_2 (k = 2) of every composition, the mixtures'
  (``Ar/blancett_1970__B_xAr_<x>.csv``) and the pure gases' (``He/..._B_He.csv``,
  ``Ar/..._B_Ar.csv``), against temperature.

The pair is filed under argon (``Ar`` sorts before ``He``) with helium as the partner,
and the mole fraction is argon's, x_Ar = 1 - x_He: the paper gives mole percent helium
(21.99, 41.05, 59.35, 80.00), so x_Ar is 0.7801, 0.5895, 0.4065 and 0.2000.

What the page images show (checked 2026-10-01 against the table crops and 300 dpi
renderings of the pdf pages): the library's scan clips the long edges of the landscape
pages. Pages 83, 85 and 87 have lost their composition labels (and, at least on p. 85,
a row cut through at the top); p. 82 shows a further row cut through below 12.3080 atm.
The compositions of the unlabelled pages are therefore set here, not read: p. 83 is
80.00 % helium, p. 85 41.05 % helium, p. 87 100 % argon. Each identification is
established by recomputing every printed Z up to 150 atm from that composition's
Table II virial series at the printed pressure: the assigned composition reproduces
all of them within 1e-5 (one unit of the last printed digit), any other is off by at
least 9e-3. Table II's sigma cells of the form ``0.0925.`` carry the page's trailing
period (as on the image), dropped from the number.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Blancett_1970_IsothermsHeAr"
TABLE_II = "tbl_c09ea288c262"
TEMPERATURES = ("323.15", "273.15", "223.15")
MAX_PRESSURE_ATM = 150.0
"""The cut: optmat's domain is below about 100 bar."""

# Table I, page by page: (library table, the composition label it prints or, for a
# page whose label the scan has lost, None, and the composition it holds).
TABLE_I = (
    ("tbl_53dbe22a4283", "100% helium", "100% helium"),
    ("tbl_63fcad0741b0", None, "80.00% helium"),
    ("tbl_44bf2fa88c63", "59.35% helium", "59.35% helium"),
    ("tbl_4f20e52886e0", None, "41.05% helium"),
    ("tbl_25b7d6538540", "21.99% helium", "21.99% helium"),
    ("tbl_bd9a3ca61717", None, "100% argon"),
)
# Argon's mole fraction for each mixture, from the printed mole percent of helium.
X_AR = {
    "80.00% helium": "0.2000",
    "59.35% helium": "0.4065",
    "41.05% helium": "0.5895",
    "21.99% helium": "0.7801",
}
PURE = {
    "100% helium": "He/blancett_1970__B_He.csv",
    "100% argon": "Ar/blancett_1970__B_Ar.csv",
}


def _is_data(cells: list[str]) -> bool:
    """Return whether a row of Table I holds values (not a header or a label)."""
    return bool(cells[1]) and cells[0][:1].isdigit()


def compressibility() -> dict[str, list[Point]]:
    """Return the mixtures' Z up to the cut, one file per composition and isotherm."""
    files: dict[str, list[Point]] = {}
    for table_id, label, composition in TABLE_I:
        grid = json_table(DIRECTORY, table_id)
        labels = [cells[0] for cells in grid if cells[0] and not _is_data(cells)]
        if label is None and any("%" in text for text in labels):
            msg = f"{table_id}: a composition label appeared: {labels!r}"
            raise ValueError(msg)
        if label is not None and label not in labels:
            msg = f"{table_id}: no label {label!r} in {labels!r}"
            raise ValueError(msg)
        if composition not in X_AR:
            continue
        rows = [cells for cells in grid if _is_data(cells)]
        for column, temperature in enumerate(TEMPERATURES):
            name = f"Ar/blancett_1970__Z_xAr_{X_AR[composition]}_{temperature}K.csv"
            points: list[Point] = []
            run, previous = 1, None
            for cells in rows:
                pressure, z, sigma = (
                    number(c) for c in cells[3 * column : 3 * column + 3]
                )
                if previous is not None and float(pressure) > previous:
                    run += 1  # a run starts again at its highest pressure
                previous = float(pressure)
                if float(pressure) > MAX_PRESSURE_ATM:
                    continue
                printed = f"run {run} | {pressure} | {z} | {sigma}"
                points.append(Point(pressure, pressure, z, f"{sigma}e-5", printed))
            if run != 2 or len({p.row for p in points}) != len(points):
                msg = f"{name}: {run} runs or repeated pressures"
                raise ValueError(msg)
            files[name] = points
    return files


def second_virial() -> dict[str, list[Point]]:
    """Return Table II's B_2 of each composition against temperature."""
    grid = json_table(DIRECTORY, TABLE_II)
    if grid[0][0] != "Composition" or grid[1][2] != "B_k":
        msg = f"unexpected Table II header {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        if not cells[0]:
            continue  # k > 2: the row of a composition's label is its k = 2 row
        composition = cells[0]
        name = PURE.get(composition)
        if name is None:
            name = f"Ar/blancett_1970__B_xAr_{X_AR[composition]}.csv"
        points = []
        for column, temperature in enumerate(TEMPERATURES):
            k, value, sigma = cells[1 + 3 * column : 4 + 3 * column]
            if k != "2":
                msg = f"Table II {composition}: k = {k!r}, not 2"
                raise ValueError(msg)
            printed = f"{composition} | {temperature} K | k = 2 | {value} | {sigma}"
            points.append(
                Point(
                    temperature,
                    temperature,
                    number(value),
                    number(sigma.removesuffix(".")),
                    printed,
                )
            )
        files[name] = points
    if len(files) != len(X_AR) + len(PURE):
        msg = f"Table II: {sorted(files)}"
        raise ValueError(msg)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the mixtures' Z (Table I) and every composition's B_2 (Table II)."""
    return compressibility() | second_virial()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
