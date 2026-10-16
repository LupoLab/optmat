# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Record solprop's ionisation rates of the molecules as a test fixture.

    cd ../solprop && uv run python ../optmat/tools/solprop_molecular_rates.py

runs in *solprop's* environment (this package never imports solprop) and writes
``tests/data/solprop_molecular_rates.json``, with the solprop commit it came from: for
every molecule, the rates of each of solprop's three molecular forms on a grid of
fields, with the ionisation energy and the settings each was evaluated with:

- ``moadk``: ``IsotropicMOADK`` with each tabulation of the structure coefficients
  (H2, D2, N2, O2), instantaneous;
- ``ppt``: ``PPT`` with Talebpour et al.'s effective charge (N2, O2), at 800 and 400
  nm, instantaneous and cycle averaged, with solprop's defaults (occupancy 2; no
  barrier-suppression, Stark or dipole factor exists for a molecule);
- ``adk``: ``ADK``, which for a molecule without quantum numbers is the form built on
  the ionisation energy alone (H2, D2, CH4, SF6, N2O) and for N2 and O2 takes the
  effective charge; instantaneous and cycle averaged, occupancy 2.

These are the parity targets of the molecular rates' literature forms in
``optmat.gas.molecular_ionisation_rates``.
"""

import json
import subprocess
from pathlib import Path

import numpy as np

# solprop is not a dependency of this package: this tool runs in solprop's environment.
from solprop.exceptions import MaterialError  # pyright: ignore[reportMissingImports]
from solprop.nonlinear.ionisation import (  # pyright: ignore[reportMissingImports]
    ADK,
    PPT,
    IsotropicMOADK,
)
from solprop.physdata import Gas  # pyright: ignore[reportMissingImports]
from solprop.physdata.ionisation_data import (  # pyright: ignore[reportMissingImports]
    ionisation_potential,
)

TESTS = Path(__file__).resolve().parents[1] / "tests" / "data"
MOLECULES = ("H2", "D2", "N2", "O2", "CH4", "SF6", "N2O")
FIELDS = np.geomspace(8e9, 6e10, 10)
WAVELENGTHS = (800e-9, 400e-9)


def _commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path.cwd(),
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()


def _rates(rate: object) -> list[float]:
    return [float(value) for value in rate(FIELDS)]  # pyright: ignore[reportCallIssue]


def _molecule(name: str) -> dict[str, object]:
    gas = Gas(name)
    entry: dict[str, object] = {"potential_j": ionisation_potential(gas)}
    moadk = {}
    for source in ("zhao2010", "tong2002"):
        try:
            rate = IsotropicMOADK(gas, coefficients=source)
        except MaterialError:
            continue
        moadk[source] = {"occupancy": rate.occupancy, "instantaneous": _rates(rate)}
    if moadk:
        entry["moadk"] = moadk
    adk = ADK(gas)
    entry["adk"] = {
        "occupancy": adk.occupancy,
        "potential_only": adk.potential_only,
        "instantaneous": _rates(adk),
        "cycle_averaged": _rates(ADK(gas, cycle_averaged=True)),
    }
    ppt = {}
    for wavelength in WAVELENGTHS:
        try:
            rate = PPT(gas, wavelength)
        except MaterialError:
            continue
        ppt[f"{wavelength * 1e9:.0f}"] = {
            "occupancy": rate.occupancy,
            "settings": [list(pair) for pair in rate.parameters],
            "instantaneous": _rates(rate),
            "cycle_averaged": _rates(PPT(gas, wavelength, cycle_averaged=True)),
        }
    if ppt:
        entry["ppt"] = ppt
    return entry


def main() -> None:
    """Evaluate solprop's molecular rates and write the fixture."""
    TESTS.mkdir(exist_ok=True)
    fixture = {
        "solprop_commit": _commit(),
        "fields_v_per_m": [float(field) for field in FIELDS],
        "molecules": {name: _molecule(name) for name in MOLECULES},
    }
    (TESTS / "solprop_molecular_rates.json").write_text(
        json.dumps(fixture, indent=1) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
