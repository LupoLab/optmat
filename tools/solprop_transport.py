# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Record solprop's current transport properties as a test fixture (OEP 1, Phase 3).

    cd ../solprop && PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \\
        ../optmat/tools/solprop_transport.py

runs in *solprop's* environment (this package never imports solprop), reading it only,
and writes ``tests/data/solprop_transport.json``, with the solprop commit it came from:
for every registry gas, its transport temperature range, its molar mass (kg/mol), and
at temperatures across the range its dilute viscosity (Pa s) and dilute heat-capacity
ratio. They are the baseline of solprop SEP 3's changed-results report, compared with
``optmat.gas.transport`` on the page *Transport of the dilute gases*.
"""

import json
import subprocess
from pathlib import Path

import numpy as np

# solprop is not a dependency of this package: this tool runs in solprop's environment.
from solprop.physdata import Gas  # pyright: ignore[reportMissingImports]
from solprop.physdata.transport import (  # pyright: ignore[reportMissingImports]
    heat_capacity_ratio,
    molar_mass,
    temperature_range,
    viscosity,
)

OUT = Path(__file__).resolve().parents[1] / "tests" / "data" / "solprop_transport.json"


def main() -> None:
    """Evaluate every gas's transport properties and write the fixture."""
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path.cwd(),
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    gases = {}
    for gas in Gas:
        if gas.name == "Air":
            continue
        low, high = temperature_range(gas)
        temperatures = [float(t) for t in np.linspace(low, high, 9)]
        gases[gas.name] = {
            "temperature_range_k": [low, high],
            "molar_mass_kg_per_mol": molar_mass(gas),
            "temperatures_k": temperatures,
            "viscosity_pa_s": [viscosity(gas, t) for t in temperatures],
            "heat_capacity_ratio": [heat_capacity_ratio(gas, t) for t in temperatures],
        }
    OUT.write_text(
        json.dumps({"solprop_commit": commit, "gases": gases}, indent=1) + "\n"
    )


if __name__ == "__main__":
    main()
