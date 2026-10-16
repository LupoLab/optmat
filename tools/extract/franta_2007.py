# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Franta, Nečas & Zajíčková (2007): the PDOS model parameters and band gaps.

The three tables are the fitted parameters of the parameterised-density-of-states
(PDOS) model:

- Table 1 (``tbl_363cc8a41358``): six DLC films (as deposited and annealed at 280-510
  C), Q, E_g and E_h of the pi and sigma bands, then quantities derived from them and
  from other measurements (the thickness ratio, alpha, beta, X_H, the sp3/sp2 ratio);
- Table 2 (``tbl_54c954601fc2``): a typical a-Si:H film, Eqs. (19)-(23);
- Table 3 (``tbl_e67d945ff5b3``): an SiO2-like film, Eqs. (9), (10), (24), (25).

The fitted parameters are the record's ``[[formulas]]`` (one per film), transcribed from
the images; this script checks each recorded coefficient against the table cells, so
that the record and the library cannot drift apart. The only vocabulary quantity is the
band gap (``band_gap``, a model parameter): E_g_pi of each DLC film, E_g of the a-Si:H
and of the SiO2-like film, three data files.

All cells checked against the table images (2026-09-27): they agree. Table 1's "-" are
the rows that have no value for the as-deposited film.
"""

import tomllib

from _library import DATA, Point, number, table, write

DIRECTORY = "Franta_2007_ModelsDielectricResponse"
TABLE_1 = "tbl_363cc8a41358"
TABLE_2 = "tbl_54c954601fc2"
TABLE_3 = "tbl_e67d945ff5b3"
SIGMA = "\N{GREEK SMALL LETTER SIGMA}"
FILMS = ("as deposited", "280 ◦C", "325 ◦C", "406 ◦C", "460 ◦C", "510 ◦C")
FILM_IDS = ("as_deposited", "280C", "325C", "406C", "460C", "510C")
DLC_ROWS = {
    "Qπ [eV3/2]": "Q_pi",
    "Egπ [eV]": "E_gpi",
    "Ehπ [eV]": "E_hpi",
    f"Q{SIGMA} [eV3/2]": "Q_sigma",
    f"Eg{SIGMA} [eV]": "E_gsigma",
    f"Eh{SIGMA} [eV]": "E_hsigma",
}
A_SI_COLUMNS = (
    "A_xi", "E_g", "E_h", "A_1", "B_1", "E_1", "A_2", "B_2", "E_2", "A_lambda",
    "E_lambda",
)  # fmt: skip
SIO2_COLUMNS = (
    "Q_xi", "E_g", "E_h", "A_delta1", "B_delta1", "E_delta1", "A_delta2", "B_delta2",
    "E_delta2", "A_p", "B_p", "E_p",
)  # fmt: skip


def _record() -> dict[str, dict[str, object]]:
    """Return the record's formulas by id."""
    path = DATA / "sources" / "franta_2007.toml"
    with path.open("rb") as handle:
        return {entry["id"]: entry for entry in tomllib.load(handle)["formulas"]}


def _check(formula: dict[str, object], expected: dict[str, str]) -> None:
    """Check a recorded formula's coefficients against the cells."""
    if formula["coefficients"] != expected:
        msg = f"{formula['id']}: the record differs from the table, {expected}"
        raise ValueError(msg)


def _dlc(formulas: dict[str, dict[str, object]]) -> list[Point]:
    rows = table(DIRECTORY, TABLE_1)
    if tuple(rows[0][1:]) != FILMS:
        msg = f"unexpected Table 1 header {rows[0]!r}"
        raise ValueError(msg)
    by_label = {cells[0]: cells for cells in rows[1:]}
    points = []
    for column, (film, film_id) in enumerate(zip(FILMS, FILM_IDS, strict=True), 1):
        expected = {
            name: number(by_label[label][column]) for label, name in DLC_ROWS.items()
        }
        _check(formulas[f"pdos_dlc_{film_id}"], expected)
        cells = by_label["Egπ [eV]"]
        printed = f"{cells[0]} | {film} | {cells[column]}"
        label = film.replace(" ◦C", " C")
        points.append(Point(film_id, label, expected["E_gpi"], "", printed))
    return points


def _single(table_id: str, names: tuple[str, ...], formula: dict[str, object]) -> Point:
    rows = table(DIRECTORY, table_id)
    if len(rows) != 3 or len(rows[2]) != len(names):
        msg = f"unexpected layout of {table_id}: {rows!r}"
        raise ValueError(msg)
    expected = {name: number(cell) for name, cell in zip(names, rows[2], strict=True)}
    _check(formula, expected)
    printed = " | ".join(f"{h} ({u}) {v}" for h, u, v in zip(*rows, strict=True))
    return Point("1", "", expected["E_g"], "", printed)


def extract() -> dict[str, list[Point]]:
    """Return the band gaps, after checking the recorded formulas."""
    formulas = _record()
    return {
        "a-C_H/franta_2007__dlc_band_gap_pi.csv": _dlc(formulas),
        "a-Si_H/franta_2007__a_si_h_band_gap.csv": [
            _single(TABLE_2, A_SI_COLUMNS, formulas["pdos_a_si_h"])
        ],
        "SiO2-like_film_HMDSO/franta_2007__sio2_like_band_gap.csv": [
            _single(TABLE_3, SIO2_COLUMNS, formulas["pdos_sio2_like"])
        ],
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
