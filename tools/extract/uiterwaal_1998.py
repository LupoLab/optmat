# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Check Uiterwaal et al. (1998)'s record against the paper: Table I and Sec. VI.

The source's values are printed in Table I (``tbl_02103ed7cb91``, the generalized
cross sections, image checked) and in the text of Sec. VI (the saturation intensities,
p. 399); they are recorded inline in the source record and this script writes no data
file. It checks:

* each inline cross section and its printed errors against Table I's Experiment row,
  whose cells read ``(0.9^{-0.4}_{+9}) \times 10^{-152}``: value, superscript and
  subscript, power (the record reads the larger error as the upper one);
* each inline saturation intensity and its error against the sentence of Sec. VI;
* Eq. (9), sigma_(N) (I_sat / hbar omega)^N t_eff = 1, with t_eff = int F^N dt for the
  authors' sech^2 pulse of 500 fs FWHM, int sech^(2N)(t/T) dt = T sqrt(pi) Gamma(N) /
  Gamma(N + 1/2), T = FWHM / (2 arccosh sqrt 2): the saturation intensities give the
  printed cross sections within ``AGREEMENT``.
"""

import math
import re
import tomllib
from typing import Any

from _library import DATA, LIBRARY, Point, table

from optmat.constants import E_CHARGE

DIRECTORY = "Uiterwaal_1998_GeneralizedMultiphotonIonization"
SOURCE = "uiterwaal_1998"
TABLE_I = "tbl_02103ed7cb91"
# Table I's column per gas, and the order N.
COLUMNS = {"He": (1, 5), "Ar": (3, 4), "Kr": (4, 3), "Xe": (5, 3)}
_CELL = re.compile(
    r"^\((?P<value>[\d.]+)\^\{(?P<sup>[-+][\d.]+)\}_\{(?P<sub>[-+][\d.]+)\}\)"
    r" \\times 10\^\{(?P<power>-\d+)\}$"
)
_SENTENCE = re.compile(
    r"I_\{\\mathrm\{sat,(?P<gas>He|Ar|Kr)\}\}=\((?P<value>[\d.]+)\\pm"
    r"(?P<error>[\d.]+)\)\\times10\^\{(?P<power>\d+)\}"
)
PHOTON_ENERGY_EV = 4.99  # p. 397: 'our photon energy of hbar omega = 4.99 eV'
FWHM_S = 500e-15
# Eq. (9) against Table I: the cross sections are printed to one or two digits and the
# saturation intensities to two, an N-th power apart.
AGREEMENT = 0.25


def _record() -> dict[str, dict[str, Any]]:
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        return {d["id"]: d for d in tomllib.load(handle)["datasets"]}


def _text() -> str:
    r"""Return the paper's text with whitespace and math delimiters removed.

    The extraction spaces out some formulae and splits others ('I _ { \mathrm { s a t ,
    A r } } = ( 1 . 0 \pm 0 . 1 5 )$ $\times 1 0 ^ { 1 4 }'); without spaces and '$'
    every sentence reads the same way.
    """
    text = (LIBRARY / DIRECTORY / "paper.md").read_text(encoding="utf-8")
    return re.sub(r"[\s$]", "", text)


def effective_duration(photons: int) -> float:
    """Return t_eff = int F^N dt for a sech^2 pulse of 500 fs FWHM, seconds."""
    width = FWHM_S / (2 * math.acosh(math.sqrt(2)))
    return width * math.sqrt(math.pi) * math.gamma(photons) / math.gamma(photons + 0.5)


def check_table() -> None:
    """Raise unless the record's cross sections and errors are Table I's."""
    grid = table(DIRECTORY, TABLE_I)
    (row,) = [cells for cells in grid if cells[0] == "Experiment"]
    record = _record()
    for gas, (column, photons) in COLUMNS.items():
        match = _CELL.match(row[column].strip())
        if match is None:
            msg = f"unexpected Table I cell {row[column]!r}"
            raise ValueError(msg)
        dataset = record[f"sigma{photons}_{gas}"]
        (point,) = dataset["points"]
        errors = sorted((match["sup"], match["sub"]), key=float)
        parameters = dataset["parameters"]
        printed = (match["value"], f"1e{match['power']}", errors[1][1:], errors[0][1:])
        recorded = (
            point[2],
            str(dataset["scale"]),
            str(parameters["error_plus"]),
            str(parameters["error_minus"]),
        )
        if [float(x) for x in printed] != [float(x) for x in recorded]:
            msg = f"{gas}: Table I prints {printed}, the record has {recorded}"
            raise ValueError(msg)


def saturation_intensities() -> dict[str, tuple[float, float]]:
    """Return Sec. VI's saturation intensities and errors, W/cm2, as printed."""
    found = {
        m["gas"]: (
            float(f"{m['value']}e{m['power']}"),
            float(f"{m['error']}e{m['power']}"),
        )
        for m in _SENTENCE.finditer(_text())
    }
    if set(found) != {"He", "Ar", "Kr"}:
        msg = f"Sec. VI's saturation intensities not found: {sorted(found)}"
        raise ValueError(msg)
    return found


def check_text() -> None:
    """Raise unless the record's saturation intensities are Sec. VI's."""
    record = _record()
    for gas, (value, error) in saturation_intensities().items():
        dataset = record[f"isat_{gas}"]
        (point,) = dataset["points"]
        scale = float(dataset["scale"])
        recorded = (float(point[2]) * scale, float(point[3]) * scale)
        if not all(
            math.isclose(a, b) for a, b in zip(recorded, (value, error), strict=True)
        ):
            msg = f"{gas}: Sec. VI prints {(value, error)}, the record has {recorded}"
            raise ValueError(msg)


def reproduced() -> dict[str, float]:
    """Return sigma_(N) from each saturation intensity by Eq. (9), cm^2N s^(N-1)."""
    photon = PHOTON_ENERGY_EV * E_CHARGE
    out = {}
    for gas, (intensity, _) in saturation_intensities().items():
        photons = COLUMNS[gas][1]
        flux = intensity / photon  # cm^-2 s^-1
        out[gas] = 1 / (flux**photons * effective_duration(photons))
    return out


def extract() -> dict[str, list[Point]]:
    """Check the record against the paper; there are no data files to write."""
    check_table()
    check_text()
    record = _record()
    for gas, sigma in reproduced().items():
        dataset = record[f"sigma{COLUMNS[gas][1]}_{gas}"]
        (point,) = dataset["points"]
        printed = float(point[2]) * float(dataset["scale"])
        if abs(sigma / printed - 1) > AGREEMENT:
            msg = f"{gas}: Eq. (9) gives {sigma:.3g}, Table I prints {printed:.3g}"
            raise ValueError(msg)
    return {}


def main() -> None:
    """Run the checks and print the reproduced cross sections."""
    extract()
    for gas, sigma in reproduced().items():
        print(f"{gas}: sigma_({COLUMNS[gas][1]}) = {sigma:.3g} by Eq. (9)")


if __name__ == "__main__":
    main()
