# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ritchie (1977): two-photon ionisation cross sections of He, Tables I and II.

The PDF is a scan (one 4938 x 6686 bilevel image per page, 612 dpi) with a text layer
(Times-Roman, not embedded). Every cell of both tables was checked against the table
images and the page image at its native pixels (all agree; the specks over some
decimal points of Table II's E column are dust on the scan, not marks of the numbers).

* **Table I** (p. 2083, ``tbl_2491aff218c3``): 'Values of sigma^(2) x I_0/2I in a_0^2
  for two-photon ionization of He in the 1 1S state, approximated by Hartree-Fock
  theory', at 22 photon energies E_p (au, the wavelength in angstrom beside it) from
  0.456 to 0.868, below the one-photon threshold and across the 1s2p and 1s3p
  resonances. Two columns, NE (nonexchange) and E (exchange); each cell holds a value
  and, in parentheses, the value 'obtained ignoring final-state, perturbed-initial-state
  nonorthogonality'. Dipole-velocity form.
* **Table II** (p. 2083, ``tbl_7b95e928b8ad``): the same with the one-parameter
  variational ground state (zeta = 27/16), at ten photon energies from 0.456 to 0.799.

One file per table, column and value (the plain value or the parenthesised one): eight
files. The row ids are the photon energies as printed; ``printed`` keeps the row.

**Checks** (:func:`checks`; the script fails if one is missed): each printed wavelength
against the photon energy beside it (hc / E_h = 455.634 angstrom au, within the two
roundings); the text's numbers (pp. 2081-2085): the first unshifted resonance at E_p =
0.793 au for the orbital energy -0.918 au and at 0.771 for -0.896 (nu = 2 in Eq. A1:
E_p = |epsilon^(0)| - 1/8), the orbital energy -0.896 of zeta = 27/16 (zeta^2 / 2 - 2
zeta + 5 zeta / 8), the accurate ionisation potential 0.904 au and first-resonance
position 0.780 au (nist_asd_2026's He ionisation energy and 1s2p 1P1 level, rounded),
and the 'about 20 %' dynamical error at threshold (1 - E / 2 NE at 0.456 au).

**Reported** (:func:`comparisons`): Table I's E column converted to the generalised
cross section sigma_2 (Gamma = sigma_2 Phi^2, Phi = I / hbar omega), sigma_2 = 2 v a_0^2
E_p E_h / I_0 for a printed v, against the curated He calculations at the
off-resonance photon energies. The conversion is the caption's definition with Eq. 1
(sigma^(2) the intensity-dependent cross section, rate sigma^(2) Phi); lhuillier_1986's
Fig. 13 draws Table I's E column as W / I^2 and agrees with it (that record's checks).
"""

import functools
import math
import re

import numpy as np
from _library import Point, table, write

from optmat.constants import AU_ENERGY, BOHR_RADIUS, E_CHARGE, HBAR, C
from optmat.datasets import DATA_DIRECTORY, Dataset, load_source

DIRECTORY = "Ritchie_1977_CalculationCrossSections"
SOURCE = "ritchie_1977"
TABLES = {"tbl1": "tbl_2491aff218c3", "tbl2": "tbl_7b95e928b8ad"}
HEADER = ["E_p (a.u.)", "NE", "E"]
ROWS = {"tbl1": 22, "tbl2": 10}
COLUMNS = {"NE": 1, "E": 2}
FILES = {
    (name, column, which): f"He/ritchie_1977__{name}_{column}{suffix}_He.csv"
    for name in TABLES
    for column in COLUMNS
    for which, suffix in (("plain", ""), ("nonorth", "_nonorth"))
}
"""The data files, by (table, column, plain or parenthesised value)."""

_CELL = re.compile(r"^(\d+\.\d+)\((\d+\.\d+)\)$")
_ENERGY = re.compile(r"^(\d\.\d{3})(?: \((\d+) Å\))?$")


@functools.cache
def rows(name: str) -> tuple[tuple[str, str, str, str, str, str], ...]:
    """Return a table's rows, checked: (energy, wavelength, NE, NE (), E, E ()).

    The wavelength is '' in Table II, which prints none.
    """
    grid = table(DIRECTORY, TABLES[name])
    if grid[0] != HEADER:
        msg = f"{name}: unexpected header {grid[0]!r}"
        raise ValueError(msg)
    out = []
    for cells in grid[1:]:
        energy = _ENERGY.match(cells[0])
        ne = _CELL.match(cells[1])
        e = _CELL.match(cells[2])
        if energy is None or ne is None or e is None:
            msg = f"{name}: unexpected row {cells!r}"
            raise ValueError(msg)
        if (energy.group(2) is None) != (name == "tbl2"):
            msg = f"{name}: a wavelength where none is expected, or none: {cells!r}"
            raise ValueError(msg)
        out.append(
            (
                energy.group(1),
                energy.group(2) or "",
                ne.group(1),
                ne.group(2),
                e.group(1),
                e.group(2),
            )
        )
    if len(out) != ROWS[name]:
        msg = f"{name}: {len(out)} rows, not {ROWS[name]}"
        raise ValueError(msg)
    return tuple(out)


def _printed(row: tuple[str, ...]) -> str:
    energy, wavelength, ne, ne_p, e, e_p = row
    first = f"{energy} ({wavelength} Å)" if wavelength else energy
    return f"{first} | {ne}({ne_p}) | {e}({e_p})"


def points(name: str, column: str, which: str) -> list[Point]:
    """Return one column's plain or parenthesised values as data points."""
    index = {"NE": 2, "E": 4}[column] + (which == "nonorth")
    return [Point(row[0], row[0], row[index], "", _printed(row)) for row in rows(name)]


HC_ANGSTROM_AU = 2 * math.pi * HBAR * C / AU_ENERGY * 1e10
"""hc / E_h in angstrom (455.634): the vacuum wavelength of a photon of 1 hartree."""
HARTREE_EV = AU_ENERGY / E_CHARGE
TEXT = {
    "unshifted resonance, Table I": (0.918, 0.793),
    "unshifted resonance, Table II": (0.896, 0.771),
}
"""pp. 2084-2085: the orbital energy |epsilon^(0)| (au) and the first resonance's E_p
(au) 'unshifted by the atomic field' (the first pole of Eq. A1, nu = 2)."""
ZETA = 27 / 16
ORBITAL_II = -0.896
"""p. 2084: epsilon^(0) = zeta^2 / 2 - 2 zeta + 5 zeta / 8 'or -0.896 a.u.'."""
ACCURATE_IP = 0.904
"""p. 2081: 'the accurate ionization potential V_I = 0.904 a.u.'."""
ACCURATE_RESONANCE = 0.780
"""p. 2084: the first resonance from 'the accurate positions of the 1 1S and 2 1P
eigenstates of He occurs at 0.780 a.u.'."""
DYNAMICAL = 0.20
"""p. 2085: the dynamical error of the nonexchange model 'is found to be about 20 % at
threshold' (the deviation of E / NE from 2)."""


def _nist(key: str) -> float:
    """Return a nist_asd_2026 He value in hartree (the IE in eV, the level in cm-1)."""
    record = load_source(DATA_DIRECTORY / "sources" / "nist_asd_2026.toml")
    found = next(d for d in record.datasets if d.key == f"nist_asd_2026/{key}")
    value = float(found.scaled_values()[0])
    if found.unit == "eV":
        return value / HARTREE_EV
    return value * 100 * HC_ANGSTROM_AU * 1e-10


def checks() -> dict[str, object]:
    """Return the checks of Tables I and II against the text; raise if one is missed."""
    out: dict[str, object] = {}
    # The wavelengths: E_p and lambda are each rounded (to 0.001 au and to 1 angstrom).
    worst = 0.0
    for row in rows("tbl1"):
        energy, wavelength = float(row[0]), float(row[1])
        allowed = 0.0005 + HC_ANGSTROM_AU * 0.5 / wavelength**2
        miss = abs(HC_ANGSTROM_AU / wavelength - energy)
        if miss > allowed:
            msg = f"Table I: {row[0]} au is not {row[1]} angstrom"
            raise ValueError(msg)
        worst = max(worst, miss / allowed)
    out["Table I wavelengths, worst miss / allowed"] = worst
    for label, (orbital, resonance) in TEXT.items():
        unshifted = orbital - 1 / (2 * 2**2)
        if round(unshifted, 3) != resonance:
            msg = f"{label}: |epsilon0| - 1/8 = {unshifted:.4f}, not {resonance}"
            raise ValueError(msg)
        out[label] = unshifted
    orbital = ZETA**2 / 2 - 2 * ZETA + 5 * ZETA / 8
    if round(orbital, 3) != ORBITAL_II:
        msg = f"zeta = 27/16 gives epsilon0 = {orbital:.5f}, not -0.896"
        raise ValueError(msg)
    out["epsilon0 for zeta = 27/16"] = orbital
    ip, level = _nist("ie_He_I"), _nist("level_He_1s2p_1P1")
    if round(ip, 3) != ACCURATE_IP or round(level, 3) != ACCURATE_RESONANCE:
        msg = f"NIST: IP {ip:.5f}, 2 1P {level:.5f} au; the text 0.904, 0.780"
        raise ValueError(msg)
    out["NIST IP and 2 1P, au"] = (ip, level)
    for name in TABLES:
        first = rows(name)[0]
        error = 1 - float(first[4]) / (2 * float(first[2]))
        if round(error, 2) != DYNAMICAL:
            msg = f"{name}: 1 - E / 2 NE at {first[0]} au is {error:.3f}, not 0.20"
            raise ValueError(msg)
        out[f"{name}: 1 - E / 2 NE at {first[0]} au"] = error
    return out


I0_W_CM2 = 14.038e16
"""Eq. 1: 'I_0 = 14.038 x 10^16 W/cm^2'."""


def sigma2_cm4_s(value: float, energy_au: float) -> float:
    """Return a printed v = sigma^(2) I_0 / 2I as the generalised sigma_2, cm^4 s.

    sigma^(2) = v a_0^2 (2I / I_0) is a cross section (cm^2) for the rate sigma^(2)
    Phi, so sigma_2 = sigma^(2) / Phi = 2 v a_0^2 hbar omega / I_0 (I in W/cm2).
    """
    a0_cm = BOHR_RADIUS * 100
    return 2 * value * a0_cm**2 * energy_au * AU_ENERGY / I0_W_CM2


def w_over_i2(value: float, energy_au: float) -> float:
    """Return a printed v as W / I^2 = sigma_2 / (hbar omega)^2, s^-1 W^-2 cm^4."""
    return sigma2_cm4_s(value, energy_au) / (energy_au * AU_ENERGY) ** 2


OFF_RESONANCE = ("0.456", "0.506", "0.570", "0.651", "0.701")
"""Table I's rows below the 1s2p resonance region (E_p < 0.75 au)."""
OTHERS = {
    "mihelic_2021 (L = 0 + 2)": (
        "mihelic_2021/fig3a_sigma2_L0",
        "mihelic_2021/fig3a_sigma2_L2",
    ),
    "saenz_1999 (1S + 1D, velocity)": (
        "saenz_1999/fig1a_sigma2_1S",
        "saenz_1999/fig1b_sigma2_1D",
    ),
    "nikolopoulos_2001 (S + D)": (
        "nikolopoulos_2001/fig4_sigma2_S",
        "nikolopoulos_2001/fig4_sigma2_D",
    ),
}
"""The curated He curves compared (their partials summed), interpolated in log sigma
within their own range."""
BACHAU = "bachau_1990/tbl1_sigma2_L2"
"""Compared at its 0.5 au row only (its others lie above the one-photon threshold)."""


def _dataset(key: str) -> Dataset:
    """Return another record's dataset by its key, ``<source>/<id>``."""
    source, _ = key.split("/")
    record = load_source(DATA_DIRECTORY / "sources" / f"{source}.toml")
    return next(d for d in record.datasets if d.key == key)


def _at(key: str, energy_ev: float) -> float | None:
    """Return a cm4 s dataset at a photon energy (eV), None outside its range."""
    found = _dataset(key)
    x = 2 * math.pi * HBAR * C / E_CHARGE / found.vacuum_wavelength()
    y = found.scaled_values()
    order = np.argsort(x)
    x, y = x[order], y[order]
    if not x[0] <= energy_ev <= x[-1]:
        return None
    return float(np.exp(np.interp(energy_ev, x, np.log(y))))


def model_cm4_s(energy_au: float) -> float:
    """Return OEP 6's He model (mihelic_2021, TC-943) at a photon energy, cm^4 s."""
    from optmat.gas.two_photon import two_photon_ionisation_model  # noqa: PLC0415

    omega = energy_au * AU_ENERGY / HBAR
    return float(two_photon_ionisation_model("He").cross_section(omega)) * 1e8


def comparisons() -> dict[str, dict[str, float | None]]:
    """Return Table I's E column (as sigma_2) over each other calculation, by E_p."""
    table = {row[0]: float(row[4]) for row in rows("tbl1")}
    out: dict[str, dict[str, float | None]] = {}
    for energy in OFF_RESONANCE:
        ours = sigma2_cm4_s(table[energy], float(energy))
        ev = float(energy) * HARTREE_EV
        ratios: dict[str, float | None] = {"sigma_2 (cm4 s)": ours, "eV": ev}
        for label, keys in OTHERS.items():
            parts = [_at(key, ev) for key in keys]
            known = [p for p in parts if p is not None]
            ratios[label] = ours / sum(known) if len(known) == len(parts) else None
        ratios["OEP 6 He model"] = ours / model_cm4_s(float(energy))
        out[energy] = ratios
    # bachau_1990's 0.5 au row against Table I interpolated (log sigma) to 0.5 au.
    bachau = _dataset(BACHAU)
    pairs = zip(bachau.abscissa_values(), bachau.scaled_values(), strict=True)
    half = {float(x): float(v) for x, v in pairs}[0.5]
    low, high = (sigma2_cm4_s(table[e], float(e)) for e in ("0.456", "0.506"))
    ours = float(np.exp(np.interp(0.5, [0.456, 0.506], np.log([low, high]))))
    out["0.500"] = {"sigma_2 (cm4 s)": ours, "bachau_1990 (Table I, L^2)": ours / half}
    return out


def extract() -> dict[str, list[Point]]:
    """Return the eight columns of Tables I and II, after the checks."""
    checks()
    return {
        relative: points(name, column, which)
        for (name, column, which), relative in FILES.items()
    }


def main() -> None:
    """Write the data files, and print the checks and the comparison."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")
    for energy, ratios in comparisons().items():
        print(f"Table I E at {energy} au: {ratios}")


if __name__ == "__main__":
    main()
