# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Record solprop's Raman response functions as a test fixture (TC-690).

    cd ../solprop && PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \\
        ../optmat/tools/solprop_raman_response.py

runs in *solprop's* environment (this package never imports solprop), reading it only,
and writes ``tests/data/solprop_raman_response.json`` with the solprop commit it came
from: for each gas with a Raman response in solprop, at 295 K, each term's impulse
response ``h(t)`` (rotation, vibration; their sum is the response) in m^5 V^-2 s^-1
at a set of delays and densities; and fused silica's normalised Hollenbeck & Cantrell
shape and its Lin & Agrawal shapes h_a and h_b.

These are the comparison of optmat's responses with the ones solprop uses today, on
the manual page *The Raman response*; the differences are the register's (TC-690).
"""

import json
import subprocess
from pathlib import Path

import numpy as np

# solprop is not a dependency of this package: this tool runs in solprop's environment.
from solprop.nonlinear.raman import (  # pyright: ignore[reportMissingImports]
    MolecularOptions,
    SilicaTensorResponse,
    response_for,
)
from solprop.physdata import Gas  # pyright: ignore[reportMissingImports]

TESTS = Path(__file__).resolve().parents[1] / "tests" / "data"

AMAGAT = 101325.0 / (1.380649e-23 * 273.15)
"""One amagat, m^-3, as optmat defines it (optmat.gas.linewidths.AMAGAT)."""

GASES = ("N2", "O2", "H2", "D2", "CH4", "SF6")
DENSITIES_AMAGAT = (1.0, 10.0)
TEMPERATURE = 295.0

# Delays: the first half picosecond at 1 fs (an eighth of H2's vibrational period),
# then to 60 ps at 50 fs (the rotational revivals, 8.4 ps for N2, and the decay).
DELAYS = np.concatenate(
    [np.arange(0.0, 0.5e-12, 1e-15), np.arange(0.5e-12, 60e-12, 50e-15)]
)
SILICA_DELAYS = np.arange(0.0, 1.5e-12, 1e-15)


def _rounded(values: np.ndarray) -> list[float]:
    """Return values to 8 significant figures, keeping the fixture small."""
    return [float(f"{v:.8g}") for v in np.asarray(values, dtype=float)]


def _commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path.cwd(),
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()


def _gas(name: str) -> dict[str, object]:
    entry: dict[str, object] = {}
    parts: dict[str, MolecularOptions] = {}
    for part in ("rotation", "vibration"):
        other = "vibration" if part == "rotation" else "rotation"
        try:
            response_for(
                Gas[name], MolecularOptions(temperature=TEMPERATURE, **{other: False})
            )
        except Exception:
            continue
        parts[part] = MolecularOptions(temperature=TEMPERATURE, **{other: False})
    for part, part_options in parts.items():
        response = response_for(Gas[name], part_options)
        entry[part] = {
            str(n): _rounded(response.impulse(DELAYS, n * AMAGAT))
            for n in DENSITIES_AMAGAT
        }
    return entry


def main() -> None:
    """Write the fixture."""
    silica = response_for("SiO2")
    tensor = SilicaTensorResponse()
    shapes = tensor.shapes(SILICA_DELAYS)
    data = {
        "solprop_commit": _commit(),
        "temperature_k": TEMPERATURE,
        "amagat_m3": AMAGAT,
        "delays_s": _rounded(DELAYS),
        "gases": {name: _gas(name) for name in GASES},
        "silica": {
            "delays_s": _rounded(SILICA_DELAYS),
            "hollenbeck_impulse": _rounded(silica.impulse(SILICA_DELAYS, 1.0)),
            "lin_agrawal_shapes": [_rounded(row) for row in shapes],
            "lin_agrawal_fraction": tensor.fraction,
        },
    }
    path = TESTS / "solprop_raman_response.json"
    path.write_text(json.dumps(data) + "\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
