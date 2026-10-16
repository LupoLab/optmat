# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Proulx and Shakeshaft (1993): two- and excess-photon ionisation of He.

The PDF's page 1 is the IOPscience download banner; journal pp. L7-L14 are PDF pages
2-9, scanned. Two tables hold the two-photon ionisation of ground-state He (velocity
gauge, linear polarisation, the He+ ion left in its ground state):

* **Table 1** (p. L11, ``tbl_677dfb58bbf0``): Gamma / I^2 in atomic units, below the
  one-photon threshold, at the photoelectron energies epsilon = 0.15 to 0.75 au;
  Gamma_a from the amplitude of Eq. (8), Gamma_b from the induced width of the ground
  level. Gamma_a is taken; Gamma_b checks it. The abscissa is converted from the
  photoelectron energy to the photon energy, hbar omega = (epsilon + I_p) / 2 with I_p
  the calculation's own ionisation energy (:data:`IONISATION`, from the binding
  energies the text prints); ``printed`` keeps the row as printed.
* **Table 2** (p. L11, ``tbl_8d4370f4f253``): the two-photon cross section, 10^-53
  cm^4 s, at omega = 1.0 to 1.3 au, above the one-photon threshold: 'Present' is
  taken; 'Bachau et al' is bachau_1990's Table I, checked equal and not taken.

The cells of both tables were checked against the table images (all agree).

**Checks** (:func:`checks`): Gamma_b against Gamma_a; Table 2's Bachau column against
bachau_1990's file; the atomic unit of intensity from the text's 'about 10^11 s^-1'
for an 'average' rate of 1 x I^2 au at 10^13 W/cm2 (p. L14), which fixes Gamma / I^2's
unit; and, independently, Table 1 as cross sections against saenz_1999's partial
curves summed (TC-873), within :data:`SAENZ`. Reported: Table 1 against
vanderhart_2005's dashed line (Saenz's numbers again); Table 2 against the circles of
nikolopoulos_2001's Fig. 4 (TC-875), read with that record's calibrated vector axes.
"""

import csv
import functools

import nikolopoulos_2001 as nikolopoulos
import numpy as np
from _library import DATA, Point, table, write

DIRECTORY = "Proulx_1993_TwoExcessPhoton"
SOURCE = "proulx_1993"
TABLE_1 = "tbl_677dfb58bbf0"
TABLE_2 = "tbl_8d4370f4f253"
HEADER_1 = ["\\epsilon", "{\\Gamma }_{a}/{I}^{2}", "{\\Gamma }_{b}/{I}^{2}"]
HEADER_2 = ["ω", "Bachau et al", "Present"]
FILES = {
    "table1": "He/proulx_1993__tbl1_rate2_He.csv",
    "table2": "He/proulx_1993__tbl2_sigma2_He.csv",
}
BINDING = (2.90283, 2.90316)
"""'binding energies of between 2.90283 and 2.90316 au for the ground state of He'
(p. L9), depending on the basis's kappa."""
IONISATION = (sum(BINDING) / 2) - 2.0
"""The calculation's ionisation energy, au: the binding energy less He+'s 2 au (the
midpoint; the range is +-0.00017 au, +-0.00008 au in the photon energy)."""

# CODATA 2018.
HARTREE_EV = 27.211386245988
HARTREE_J = 4.3597447222071e-18
TIME_AU = 2.4188843265857e-17
BOHR_CM = 5.29177210903e-9
INTENSITY_AU = HARTREE_J / TIME_AU / BOHR_CM**2
"""E_h / (hbar a_0^2) in W/cm2, 6.436e15: the atomic unit of intensity that the text's
rates imply (:func:`checks`)."""


@functools.cache
def _rows(table_id: str, header: tuple[str, ...]) -> tuple[tuple[str, ...], ...]:
    """Return a table's rows, after checking its header and its cells are numbers."""
    grid = table(DIRECTORY, table_id)
    if grid[0] != list(header):
        msg = f"unexpected header of {table_id}: {grid[0]!r}"
        raise ValueError(msg)
    for cells in grid[1:]:
        for cell in cells:
            float(cell)
    return tuple(tuple(cells) for cells in grid[1:])


def photon_energy(epsilon: float) -> float:
    """Return the photon energy (au) of two-photon ionisation to epsilon (au)."""
    return (epsilon + IONISATION) / 2


def _table_1() -> list[Point]:
    """Return Table 1's Gamma_a / I^2 against the photon energy (au)."""
    return [
        Point(
            cells[0],
            f"{photon_energy(float(cells[0])):.5f}",
            cells[1],
            "",
            " | ".join(cells),
        )
        for cells in _rows(TABLE_1, tuple(HEADER_1))
    ]


def _table_2() -> list[Point]:
    """Return Table 2's 'Present' column against the photon energy (au)."""
    return [
        Point(cells[0], cells[0], cells[2], "", " | ".join(cells))
        for cells in _rows(TABLE_2, tuple(HEADER_2))
    ]


def cross_section(rate: float, omega: float) -> float:
    """Return sigma_2 (cm^4 s) from Gamma / I^2 (au) at a photon energy (au).

    sigma_2 = Gamma / Phi^2 = (Gamma / I^2) omega^2 in atomic units, Phi = I / omega,
    I in units of :data:`INTENSITY_AU`; then a_0^4 t_au.
    """
    return rate * omega**2 * BOHR_CM**4 * TIME_AU


WIDTH = 0.004
"""Gamma_b against Gamma_a (the table: 0.04 to 0.20 % apart)."""
TEXT_RATE = 1e11
"""p. L14: 'the resonant two-photon rate is about 10^11 s^-1 if we take an "average"
resonant two-photon rate ... of 1 x I^2 au', at 10^13 W/cm2 (one figure)."""


SAENZ = 0.20
"""Table 1 as sigma_2 against saenz_1999 (0.88 to 1.12; two calculations, saenz_1999's
axis shifted by its authors and read from a figure): a factor of 1.2. The other
atomic unit of intensity would put them 30 times apart."""


def checks() -> dict[str, object]:
    """Return the checks of Tables 1 and 2 and of the intensity unit; raise on a miss.

    The margins are :data:`WIDTH`, :data:`SAENZ` and the text's one figure.
    """
    out: dict[str, object] = {}
    failures = []
    widths = {
        cells[0]: float(cells[2]) / float(cells[1]) - 1
        for cells in _rows(TABLE_1, tuple(HEADER_1))
    }
    out["Gamma_b / Gamma_a - 1"] = widths
    if max(abs(v) for v in widths.values()) > WIDTH:
        failures.append(f"Gamma_b differs from Gamma_a: {widths}")
    path = DATA / "values" / "ionisation" / "He" / "bachau_1990__tbl1_sigma2_L2_He.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        bachau = {r["row"]: r["value"] for r in csv.DictReader(handle)}
    for cells in _rows(TABLE_2, tuple(HEADER_2)):
        if bachau.get(cells[0]) != cells[1]:
            failures.append(f"Table 2's Bachau et al at {cells[0]}: {cells[1]}")
    out["Table 2's Bachau column equals bachau_1990's Table I"] = not failures
    # 1 x I^2 au at 1e13 W/cm2 as a rate, s^-1: with I in E_h / (hbar a_0^2) it is
    # 1.0e11 (the text's 'about 10^11'); with the other common unit, the intensity of
    # a field of 1 au amplitude (3.51e16 W/cm2), it would be 3.4e9.
    rate = (1e13 / INTENSITY_AU) ** 2 / TIME_AU
    out["1 x I^2 au at 1e13 W/cm2, s^-1"] = rate
    if round(np.log10(rate)) != round(np.log10(TEXT_RATE)):
        failures.append(f"1 x I^2 au at 1e13 W/cm2 is {rate:.3g} s^-1, not ~1e11")
    # The unit again, independently: Table 1 as sigma_2 against Saenz and
    # Lambropoulos's perturbative cross sections (saenz_1999, 1S + 1D).
    agreement = {e: sigma / saenz for e, (sigma, saenz, _) in comparisons().items()}
    out["Table 1 as sigma_2 over saenz_1999 summed"] = agreement
    if max(abs(np.log(v)) for v in agreement.values()) > np.log(1 + SAENZ):
        failures.append(f"Table 1 as sigma_2 against saenz_1999: {agreement}")
    if failures:
        raise ValueError("; ".join(failures))
    return out


def _saenz(name: str) -> tuple[np.ndarray, np.ndarray]:
    """Return a committed curve: its abscissae and values."""
    path = DATA / "values" / "ionisation" / "He" / f"{name}_He.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return (
        np.array([float(r["x"]) for r in rows]),
        np.array([float(r["value"]) for r in rows]),
    )


def comparisons() -> dict[str, tuple[float, float, float]]:
    """Return Table 1 as sigma_2 against saenz_1999 (summed) and vanderhart_2005.

    Per row (epsilon): (sigma_2, cm^4 s; saenz_1999's 1S + 1D at omega x 27.2114 eV on
    its own axis; vanderhart_2005's dashed line, the same calculation, at omega au).
    """
    s_x, s_1s = _saenz("saenz_1999__fig1a_sigma2_1S")
    d_x, d_1d = _saenz("saenz_1999__fig1b_sigma2_1D")
    v_x, v_pt = _saenz("vanderhart_2005__fig2_sigma2_saenz")
    out = {}
    for cells in _rows(TABLE_1, tuple(HEADER_1)):
        omega = photon_energy(float(cells[0]))
        sigma = cross_section(float(cells[1]), omega)
        energy = omega * HARTREE_EV
        summed = float(np.interp(energy, s_x, s_1s) + np.interp(energy, d_x, d_1d))
        out[cells[0]] = (sigma, summed, float(np.interp(omega, v_x, v_pt)))
    return out


GREEN = (85, 185, 83)
"""The fill of nikolopoulos_2001's Fig. 4 circles ('He+, Proulx and Shakesaft,1993')."""


def nikolopoulos_2001() -> dict[str, tuple[float, float, float]]:
    """Return nikolopoulos_2001's Fig. 4 circles against Table 2's 'Present' column.

    Per circle, by Table 2's omega (au): (its photon energy, eV; its value, cm^4 s;
    Table 2's value). The circles are drawn filled and outlined; the legend's is left
    out.
    """
    figure = nikolopoulos.FIGURES["fig4"]
    left, bottom, right, top = nikolopoulos.frame("fig4")
    legend = figure.legend
    if legend is None:
        msg = "nikolopoulos_2001's Fig. 4 has no legend box"
        raise ValueError(msg)
    centres = set()
    for drawn in nikolopoulos._drawn(figure.page):
        if drawn.fill != GREEN:
            continue
        west, south, east, north = drawn.box
        x, y = (west + east) / 2, (south + north) / 2
        inside = left <= x <= right and bottom <= y <= top
        if inside and not (legend[0] <= x <= legend[2] and legend[1] <= y <= legend[3]):
            centres.add((round(x, 3), round(y, 3)))
    x_axis, y_axis = nikolopoulos.axes("fig4")
    present = {
        float(cells[0]): float(cells[2]) * 1e-53
        for cells in _rows(TABLE_2, tuple(HEADER_2))
    }
    out = {}
    for x, y in sorted(centres):
        energy = float(x_axis(x))
        omega = min(present, key=lambda w: abs(w * HARTREE_EV - energy))
        out[f"{omega} au"] = (energy, float(y_axis(y)), present[omega])
    if len(out) != len(centres):
        msg = "two of nikolopoulos_2001's circles match one row of Table 2"
        raise ValueError(msg)
    return out


def extract() -> dict[str, list[Point]]:
    """Return Table 1's Gamma_a / I^2 and Table 2's cross sections, after the checks."""
    checks()
    return {FILES["table1"]: _table_1(), FILES["table2"]: _table_2()}


def main() -> None:
    """Write the data files, and print the checks and the comparisons."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")
    for epsilon, (sigma, saenz, dashed) in comparisons().items():
        print(
            f"Table 1 epsilon {epsilon} (omega {photon_energy(float(epsilon)):.4f} au):"
            f" sigma_2 {sigma:.4e}; saenz_1999 1S+1D {saenz:.4e} (ratio "
            f"{sigma / saenz:.3f}); vanderhart_2005 dashed {dashed:.4e} (ratio "
            f"{sigma / dashed:.3f})"
        )
    for label, (energy, drawn, printed) in nikolopoulos_2001().items():
        print(
            f"nikolopoulos_2001 Fig. 4 circle at {energy:.3f} eV ({label}): "
            f"{drawn:.4e}, Table 2 {printed:.4e}, ratio {drawn / printed:.4f}"
        )


if __name__ == "__main__":
    main()
