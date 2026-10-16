# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Franta et al. (2016): silica n and k from the newAD2 tables; Table 3's E_g.

The paper prints its optical constants only in figures (Figs. 20-23). The authors
publish the models' tables on the newAD2 site; the files used here are attached to the
paper in the library, with their URLs: the Lithosil Q2 plate
(``LithosilQ2-SPIE9890.Enk``) and the films X2546 and X2551 of this paper (the SPIE 9890
versions, ``SiO2-X2546-SPIE9890.Enk`` and ``SiO2-X2551-SPIE9890.Enk``). Each has 4001
rows of photon energy (eV), n and k, 0.01-100 eV on a logarithmic grid. The rows are
kept exactly as the files print them; the row id is the energy.

Which rows (the fitted ranges as the site states them):

- Lithosil Q2 was fitted from 80 cm-1 (0.0099 eV, below the table's first row) to 8.5
  eV. Above 8.5 eV the model was fitted to the shape of Palik's (Philipp's) tables, so
  those rows are separate datasets, reproduced from ``philipp_1985``.
- The films were fitted from 80 cm-1 to 45 eV; the rows above 45 eV are the model's
  extrapolation and are not extracted.

Table 3 (``tbl_a2a84171965e``, CSV; image checked) gives the band-gap parameter E_g of
the three materials, extracted as one band_gap dataset per sample. The image prints E_h
and E_x once, centred across the three columns (shared values), which the CSV puts under
X2551 only; they are not extracted (the record's formula lists them).
"""

from _library import Point, attachment, number, table, write

DIRECTORY = "Franta_2016_OpticalCharacterizationSiO"
TABLE_3 = "tbl_a2a84171965e"
PREFIX = "SiO2/franta_2016__"

LITHOSIL_FITTED_EV = 8.5
FILM_FITTED_EV = 45.0

# Table 3 column -> sample.
SAMPLES = {2: "X2546", 3: "X2551", 4: "LithosilQ2"}


def _enk(name: str) -> list[tuple[str, str, str, str]]:
    """Return a table's rows: energy, n and k as printed, and the raw line."""
    rows = []
    for line in attachment(DIRECTORY, name).read_text(encoding="utf-8").splitlines():
        cells = line.split()
        if not cells:
            continue
        energy, n, k = (number(cell) for cell in cells)
        rows.append((energy, n, k, " ".join(cells)))
    if len(rows) != 4001:
        msg = f"{name}: {len(rows)} rows, not 4001"
        raise ValueError(msg)
    return rows


def _nk(
    rows: list[tuple[str, str, str, str]], n_file: str, k_file: str
) -> dict[str, list[Point]]:
    return {
        f"{PREFIX}{n_file}.csv": [Point(e, e, n, "", raw) for e, n, _, raw in rows],
        f"{PREFIX}{k_file}.csv": [Point(e, e, k, "", raw) for e, _, k, raw in rows],
    }


def _band_gaps() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_3)
    if rows[0][2:5] != list(SAMPLES.values()):
        msg = f"unexpected Table 3 header {rows[0]!r}"
        raise ValueError(msg)
    (gap,) = [cells for cells in rows if cells[0] == "Eg"]
    printed = " | ".join(gap)
    return {
        f"{PREFIX}band_gap_{sample}.csv": [
            Point(sample, sample, number(gap[column]), "", printed)
        ]
        for column, sample in SAMPLES.items()
    }


def extract() -> dict[str, list[Point]]:
    """Return the n and k tables of the three materials and their band gaps."""
    files = _band_gaps()
    lithosil = _enk("LithosilQ2-SPIE9890.Enk")
    fitted = [r for r in lithosil if float(r[0]) <= LITHOSIL_FITTED_EV]
    palik = [r for r in lithosil if float(r[0]) > LITHOSIL_FITTED_EV]
    files |= _nk(fitted, "lithosil_q2_n", "lithosil_q2_k")
    files |= _nk(palik, "lithosil_q2_n_palik", "lithosil_q2_k_palik")
    for film in ("X2546", "X2551"):
        rows = _enk(f"SiO2-{film}-SPIE9890.Enk")
        rows = [r for r in rows if float(r[0]) <= FILM_FITTED_EV]
        files |= _nk(rows, f"film_{film}_n", f"film_{film}_k")
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
