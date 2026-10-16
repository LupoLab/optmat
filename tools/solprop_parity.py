# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Record solprop's current gas Sellmeier values as a test fixture.

    cd ../solprop && uv run python ../optmat/tools/solprop_parity.py

runs in *solprop's* environment (this package never imports solprop) and writes
``tests/data/solprop_sellmeier.json``: for every (gas, source) of
``solprop.physdata.sellmeier``, its reference state, validity range and the
susceptibility ``chi`` at that state at eight wavelengths across the range, with the
solprop commit it came from. ``tests/test_solprop_parity.py`` compares the formulas of
this package with them.
"""

import json
import subprocess
from pathlib import Path

import numpy as np

# solprop is not a dependency of this package: this tool runs in solprop's environment.
from solprop.physdata import Gas  # pyright: ignore[reportMissingImports]
from solprop.physdata.sellmeier import (  # pyright: ignore[reportMissingImports]
    gamma,
    sources,
)

C = 299792458.0
OUT = Path(__file__).resolve().parents[1] / "tests" / "data" / "solprop_sellmeier.json"


def main() -> None:
    """Evaluate every solprop gas model at its reference state and write the fixture."""
    solprop = Path.cwd()
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=solprop,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    entries = []
    for gas in Gas:
        for source in sources(gas):
            model = gamma(gas, source=source)
            low, high = model.wavelength_range
            wavelength = np.geomspace(low, high, 8)
            omega = 2 * np.pi * C / wavelength
            # gamma * rho_ref is the formula's chi at its reference state.
            chi = model.susceptibility(omega) * model.reference_density
            entries.append(
                {
                    "gas": gas.name,
                    "source": source,
                    "pressure_bar": model.reference_state.pressure_bar,
                    "temperature": model.reference_state.temperature,
                    "wavelength_m": wavelength.tolist(),
                    "chi": np.asarray(chi).tolist(),
                }
            )
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps({"solprop_commit": commit, "models": entries}, indent=1))
    print(f"wrote {len(entries)} models from solprop {commit[:8]} to {OUT}")


if __name__ == "__main__":
    main()
