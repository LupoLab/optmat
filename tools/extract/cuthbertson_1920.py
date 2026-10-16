# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract C. and M. Cuthbertson (1920): CO2, CO and CH4 dispersion at 480-671 nm.

Unlike the Cuthbertsons' 1910 noble-gas tables, these refractivities are **not**
doubled: Table I's 45,050 x 10^-8 at 5461 A is the absolute value 4505 x 10^-7 of the
text (p. 1), mu - 1 of carbon dioxide itself. Wavelengths (``lambda x 10^8``, A) are
the air wavelengths of the Li, Cd, Hg and Ag lines rounded to 1 A (5770 for the Hg
5769.6 line). Values print a thousands comma (``44,714``), dropped here.

- Table I (``tbl_475eaffb74ab``, cells file only; image checked), carbon dioxide:
  observed and calculated ``(mu - 1) 10^8``, the difference, and the observed values
  reduced to the hydrogen number density (column 5) with two calculated values
  (column 6). Only the observed and calculated columns are extracted: columns 5 and 6
  are columns 2 and 3 times a constant factor (see the source record).
- Table III (``tbl_ae9680fa9a81``, cells file only; image checked), methane: observed
  and calculated ``(mu - 1) 10^8``.
- Table II (``tbl_5e46be9ef072``, cells file only; image checked), carbon monoxide:
  observed and calculated ``(mu - 1) 10^8``, laid out as Tables I and III (two header
  rows, the "(mu - 1) 10^8" heading spanning both value columns).
- The unlabelled table of earlier carbon-dioxide values at 5461 A
  (``tbl_b6b36b76ac3b``; image checked), ``(mu - 1) 10^7``: one reproduced dataset
  per observer; G. W. Walker's ``4510 +- 5`` keeps its +- as sigma.

- The unlabelled table of earlier carbon-monoxide values at 5461 A
  (``tbl_fe180e2332d2``; image checked), ``(mu - 1) 10^7``: the observer in column 0,
  a leader of dots in column 1 and the value in column 2; one reproduced dataset per
  observer.

The earlier methane values (``tbl_860051328440``) state no wavelength and are not
extracted; nor is the table of lambda_0 of earlier carbon-monoxide work
(``tbl_ceb2362c6e63``), derived numbers.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Cuthbertson_1920_RefractionDispersionCarbon"
TABLES = {
    "CO2": "tbl_475eaffb74ab",
    "CH4": "tbl_ae9680fa9a81",
    "CO": "tbl_5e46be9ef072",
}
EARLIER = "tbl_b6b36b76ac3b"
EARLIER_CO = "tbl_fe180e2332d2"
# Earlier observer as printed -> dataset id.
OBSERVERS = {
    "Mascart": "mascart_CO2",
    "Perreau": "perreau_CO2",
    "G. W. Walker": "walker_CO2",
    "Stuckert": "stuckert_CO2",
    "Grüschke": "grueschke_CO2",
    "Koch": "koch_CO2",
}
_PLUS_MINUS = re.compile(r"^\s*([\d.]+)\s*\\pm\s*([\d.]+)\s*$")


def _value(cell: str) -> str:
    """Return a printed value without its thousands comma."""
    return number(cell.replace(",", ""))


def _dispersion(gas: str) -> dict[str, list[Point]]:
    observed, calculated = [], []
    for cells in json_table(DIRECTORY, TABLES[gas])[2:]:  # two header rows
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        observed.append(Point(wavelength, wavelength, _value(cells[1]), "", printed))
        calculated.append(Point(wavelength, wavelength, _value(cells[2]), "", printed))
    prefix = f"{gas}/cuthbertson_1920__"
    return {
        f"{prefix}experiment_{gas}.csv": observed,
        f"{prefix}fit_{gas}.csv": calculated,
    }


def _earlier() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, EARLIER)[1:]:
        observer = cells[0].strip()
        match = _PLUS_MINUS.match(cells[1])
        value, sigma = (match.groups()) if match else (number(cells[1]), "")
        files[f"CO2/cuthbertson_1920__{OBSERVERS[observer]}.csv"] = [
            Point("5461", "5461", value, sigma, " | ".join(cells))
        ]
    return files


def _earlier_co() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, EARLIER_CO)[1:]:
        dataset = OBSERVERS[cells[0].strip()].replace("_CO2", "_CO")
        files[f"CO/cuthbertson_1920__{dataset}.csv"] = [
            Point("5461", "5461", number(cells[2]), "", " | ".join(cells))
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables I to III and the earlier CO2 and CO values as data files."""
    return (
        _dispersion("CO2")
        | _dispersion("CH4")
        | _earlier()
        | _dispersion("CO")
        | _earlier_co()
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
