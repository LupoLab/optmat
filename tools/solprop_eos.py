# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Record solprop's current equation-of-state densities as a test fixture (OEP 1, 3).

    cd ../solprop && PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \\
        ../optmat/tools/solprop_eos.py

runs in *solprop's* environment (this package never imports solprop), reading it only,
and writes ``tests/data/solprop_eos.json``, with the solprop commit it came from: for
every registry gas, the number density (m^-3) at a set of pressures (bar) at 294 K and
273.15 K, and the pressure its inverse returns at that density. They are the baseline
of solprop SEP 3's changed-results report, compared with ``optmat.gas.eos``.
"""

import json
import subprocess
from pathlib import Path

# solprop is not a dependency of this package: this tool runs in solprop's environment.
from solprop.exceptions import MaterialError  # pyright: ignore[reportMissingImports]
from solprop.physdata import Gas  # pyright: ignore[reportMissingImports]
from solprop.physdata.equation_of_state import (  # pyright: ignore[reportMissingImports]
    density,
    pressure,
)

OUT = Path(__file__).resolve().parents[1] / "tests" / "data" / "solprop_eos.json"
PRESSURES_BAR = (0.1, 1.0, 10.0, 30.0)
TEMPERATURES_K = (273.15, 294.0)


def main() -> None:
    """Evaluate every gas's densities and their round trip, and write the fixture."""
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path.cwd(),
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    gases = {}
    for gas in Gas:
        rows = []
        for temperature in TEMPERATURES_K:
            for p in PRESSURES_BAR:
                try:
                    n = density(gas, p, temperature)
                except MaterialError as exc:  # the refusal is recorded
                    rows.append(
                        {
                            "pressure_bar": p,
                            "temperature_k": temperature,
                            "refused": str(exc),
                        }
                    )
                    continue
                rows.append(
                    {
                        "pressure_bar": p,
                        "temperature_k": temperature,
                        "number_density_m3": n,
                        "round_trip_bar": pressure(gas, n, temperature),
                    }
                )
        gases[gas.name] = rows
    OUT.write_text(
        json.dumps({"solprop_commit": commit, "gases": gases}, indent=1) + "\n"
    )


if __name__ == "__main__":
    main()
