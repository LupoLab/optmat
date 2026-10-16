# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Fink and Zoller (1989): two-photon ionisation of atomic hydrogen 1s.

The paper is a theory of the quantum-defect (MQDT) parametrisation of perturbative
two-photon ionisation amplitudes; its only numbers are for atomic hydrogen. The PDF is
a scan (one 1-bit image per page, 5169 x 7129 pixels, about 600 dpi; the library's
tables are read, every cell checked against the table images, all agree).

* **Table I** (p. 2939, ``tbl_790c596397dd``): H 1s, circular polarisation, at nine
  vacuum wavelengths from 924.06 to 300 angstrom: the smooth amplitude parameters
  t^s, t^o, t_3 (a.u.), the total generalised cross section sigma^(2) (cm^4 s,
  printed as ``5.530[-52]``, the bracket a power of ten) and Ref. 13's sigma^(2) for
  comparison. Only this paper's sigma^(2) is written (``value`` the number with an
  ``e`` exponent, ``printed`` the cell); the row ids are the wavelengths as printed.
* Tables II and III (pp. 2939-2940, ``tbl_b632fc4fdf31``, ``tbl_af1a4e1d9ee8``): the
  same for H 2s and 3s, excited initial states; read for the checks only.

Checks (:func:`checks`; the script fails if one is missed): every row of Tables I-III
reproduces its printed sigma^(2) from its printed t^s, t^o, t_3: sigma^(2) = 2 pi (2 pi
alpha)^2 omega^2 |T|^2 a0^4 t0 (Eq. 2.1) with T = Omega (t^(2) + t_3) / omega^4 (Eq.
3.2), Omega^2 = 2/15 (p. 2937, s states, circular light), t^(2) = t^s + i t^o above
the one-photon threshold and t^s - cot(pi nu) t^o below it (Eq. 3.9; nu the effective
quantum number of the intermediate energy E_a + omega, E_a = -1/(2 n^2)), omega the
photon energy in hartree from the wavelength (infinite nuclear mass: 911.25 angstrom
is omega = 0.50001); within :data:`ABOVE` above threshold, :data:`BELOW` for the first
row of each table, just below it (where cot(pi nu) is steep); and this paper's
sigma^(2) agrees with Ref. 13's within :data:`REFERENCE` wherever both are printed.
"""

import math

from _library import Point, table, write

from optmat.constants import AU_ENERGY, BOHR_RADIUS, E_CHARGE, EPS0, HBAR, C

DIRECTORY = "Fink_1989_QuantumDefectParametrization"
TABLES = {1: "tbl_790c596397dd", 2: "tbl_b632fc4fdf31", 3: "tbl_af1a4e1d9ee8"}
"""The library tables of Tables I, II and III, by the initial state's n (1s, 2s, 3s)."""
FILE = "H/fink_1989__tab1_sigma2_1s.csv"
ALPHA = E_CHARGE**2 / (4 * math.pi * EPS0 * HBAR * C)
"""The fine-structure constant (CODATA, through scipy)."""
ABOVE = 0.003
"""Above threshold, the reproduction of sigma^(2) is within this, relative: the
parameters are printed to 1e-4 a.u., a few 1e-3 of the amplitude in Table III (t^s +
t_3 about -0.012); found within 2e-3."""
BELOW = 0.02
"""Just below threshold, within this: the first row of each table is at nu = 8.50,
where cot(pi nu) is small and steep (dcot/dnu = -pi), so the rounding of the
wavelength and of t^o gives a few 1e-3 of sigma; found within 0.6 %."""
REFERENCE = 0.01
"""Against Ref. 13's three-digit sigma^(2), within this; found within 0.9 %."""


def _value(cell: str) -> float:
    """Return a cell printed as ``5.530[-52]`` as a number."""
    mantissa, exponent = cell.strip().rstrip("]").split("[")
    return float(f"{mantissa}e{exponent}")


def _text(cell: str) -> str:
    """Return a ``5.530[-52]`` cell with an ``e`` exponent, its digits as printed."""
    mantissa, exponent = cell.strip().rstrip("]").split("[")
    float(mantissa)
    return f"{mantissa}e{int(exponent)}"


def _sigma(n: int, wavelength: float, ts: float, to: float, t3: float) -> float:
    """Return sigma^(2) (cm^4 s) by Eqs. 2.1, 3.2 and 3.9: H ns, circular light."""
    omega = 2 * math.pi * HBAR * C / (wavelength * 1e-10) / AU_ENERGY
    energy = -1 / (2 * n**2) + omega
    if energy > 0:
        amplitude = abs(complex(ts + t3, to)) ** 2
    else:
        nu = 1 / math.sqrt(-2 * energy)
        amplitude = (ts + t3 - to / math.tan(math.pi * nu)) ** 2
    t_squared = (2 / 15) * amplitude / omega**8
    unit = (BOHR_RADIUS * 100) ** 4 * HBAR / AU_ENERGY
    return 2 * math.pi * (2 * math.pi * ALPHA) ** 2 * omega**2 * t_squared * unit


def checks() -> dict[str, object]:
    """Return the checks of Tables I-III; raise if one fails its margin."""
    out: dict[str, object] = {}
    failures: list[str] = []
    for n, table_id in TABLES.items():
        name = f"Table {'I' * n}"
        rows = table(DIRECTORY, table_id)[1:]
        worst_above, worst_reference = 0.0, 0.0
        for k, (wavelength, ts, to, t3, sigma, reference) in enumerate(rows):
            computed = _sigma(n, float(wavelength), float(ts), float(to), float(t3))
            ratio = computed / _value(sigma)
            if k == 0:
                out[f"{name} {wavelength} A (below threshold)"] = round(ratio, 4)
                if abs(ratio - 1) > BELOW:
                    failures.append(f"{name} {wavelength}: {ratio:.4f}")
            else:
                worst_above = max(worst_above, abs(ratio - 1))
            if reference.strip():
                ref = _value(sigma) / _value(reference)
                worst_reference = max(worst_reference, abs(ref - 1))
        out[f"{name} worst reproduction above threshold"] = round(worst_above, 5)
        out[f"{name} worst against Ref. 13"] = round(worst_reference, 4)
        if worst_above > ABOVE:
            failures.append(f"Table {'I' * n}: reproduced within {worst_above:.4f}")
        if worst_reference > REFERENCE:
            failures.append(f"Table {'I' * n}: Ref. 13 within {worst_reference:.4f}")
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return Table I's sigma^(2) of H 1s, after the checks."""
    checks()
    rows = table(DIRECTORY, TABLES[1])[1:]
    return {
        FILE: [
            Point(row[0].strip(), row[0].strip(), _text(row[4]), "", row[4].strip())
            for row in rows
        ]
    }


def main() -> None:
    """Write the data file, and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"{label}: {result}")


if __name__ == "__main__":
    main()
