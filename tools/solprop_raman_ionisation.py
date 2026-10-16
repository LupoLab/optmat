# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Record solprop's current Raman and ionisation data as test fixtures.

    cd ../solprop && uv run python ../optmat/tools/solprop_raman_ionisation.py

runs in *solprop's* environment (this package never imports solprop) and writes, with
the solprop commit they came from:

- ``tests/data/solprop_raman.json``: for every gas with Raman data in
  ``solprop.physdata.raman_data``, the rotational constants (B and D in m^-1, the
  nuclear-spin weights, the anisotropy of the volume polarizability in m^3) and the
  vibrational ones (angular frequency in rad/s, the derivative of the mean volume
  polarizability in m^2, the reduced mass in kg), each with its linewidth law (FWHM in
  Hz as ``constant + inverse_density / n + density * n``, n in amagat) and citation;
- ``tests/data/solprop_ionisation.json``: every gas's first ionisation potential in J,
  the quantum numbers of the ADK and PPT rates, the static polarizabilities of the
  neutral and the cation in C m^2/V, the MO-ADK coefficients of both tabulations, the
  barrier-suppression coefficients and the dry-air fractions.

These are the parity targets of the Raman and ionisation work (TC-439): every value
here is to have a curated counterpart in optmat, or a stated reason why not.
"""

import json
import subprocess
from dataclasses import asdict
from pathlib import Path

# solprop is not a dependency of this package: this tool runs in solprop's environment.
from solprop.exceptions import MaterialError  # pyright: ignore[reportMissingImports]
from solprop.nonlinear import ionisation  # pyright: ignore[reportMissingImports]
from solprop.physdata import (  # pyright: ignore[reportMissingImports]
    Gas,  # pyright: ignore[reportMissingImports]
    ionisation_data,
    raman_data,
)
from solprop.physdata.mixtures import (  # pyright: ignore[reportMissingImports]
    AIR,
    AIR_NITROGEN_OXYGEN,
)

TESTS = Path(__file__).resolve().parents[1] / "tests" / "data"


def _commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path.cwd(),
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()


def _raman() -> dict[str, object]:
    entries: dict[str, object] = {}
    for gas in Gas:
        status = raman_data.raman_data_status(gas)
        if status != "available":
            entries[gas.name] = {"status": status}
            continue
        model = raman_data.raman_parameters(gas)
        entries[gas.name] = {
            "status": status,
            "rotational": (
                None if model.rotational is None else asdict(model.rotational)
            ),
            "vibrational": (
                None if model.vibrational is None else asdict(model.vibrational)
            ),
        }
    return entries


def _ionisation() -> dict[str, object]:
    gases: dict[str, object] = {}
    for gas in Gas:
        entry: dict[str, object] = {}
        try:
            entry["potential_j"] = ionisation_data.ionisation_potential(gas)
        except MaterialError as error:
            entry["potential_j"] = None
            entry["potential_missing"] = str(error)
        if ionisation_data.has_quantum_numbers(gas):
            entry["quantum_numbers"] = asdict(ionisation_data.quantum_numbers(gas))
        try:
            entry["polarizability_si"] = ionisation_data.polarisability(gas)
            entry["ion_polarizability_si"] = ionisation_data.polarisability(
                gas, ion=True
            )
        except MaterialError:
            pass
        for source in ("tong2002", "zhao2010"):
            try:
                coefficients = ionisation_data.moadk_coefficients(gas, source=source)
            except MaterialError:
                continue
            entry[f"moadk_{source}"] = asdict(coefficients)
        if gas in ionisation.BSI_TONG_LIN_ALPHA:
            entry["bsi_tong_lin_alpha"] = ionisation.BSI_TONG_LIN_ALPHA[gas]
        if gas in ionisation.BSI_ZHANG_COEFFICIENTS:
            entry["bsi_zhang"] = list(ionisation.BSI_ZHANG_COEFFICIENTS[gas])
        gases[gas.name] = entry
    return {
        "gases": gases,
        "moadk_default_gases": sorted(g.name for g in ionisation.MOADK_DEFAULT_GASES),
        "default_occupancy": ionisation._DEFAULT_OCCUPANCY,
        "air": {str(k): v for k, v in AIR.items()},
        "air_nitrogen_oxygen": {str(k): v for k, v in AIR_NITROGEN_OXYGEN.items()},
    }


def main() -> None:
    """Evaluate solprop's Raman and ionisation data and write the two fixtures."""
    commit = _commit()
    TESTS.mkdir(exist_ok=True)
    (TESTS / "solprop_raman.json").write_text(
        json.dumps({"solprop_commit": commit, "raman": _raman()}, indent=1) + "\n",
        encoding="utf-8",
    )
    (TESTS / "solprop_ionisation.json").write_text(
        json.dumps({"solprop_commit": commit, **_ionisation()}, indent=1) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
