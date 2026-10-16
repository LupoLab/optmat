# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Record solprop's current solid indices as a test fixture.

    cd ../solprop && uv run python ../optmat/tools/solprop_solids.py

runs in *solprop's* environment (this package never imports solprop) and writes
``tests/data/solprop_solids.json``: for every window glass and source of
``solprop.physdata.glass`` (CaF2, MgF2, BaF2, sapphire, N-BK7), fused silica's sources,
and the stitched silver and aluminium tables of ``solprop.physdata.metals``, the model's
name, source, citation, validity range and complex index n + ik at 24 wavelengths across
the range, with the solprop commit it came from. The parity tests compare optmat's
recommended solid models with them (TC-544).
"""

import json
import subprocess
from pathlib import Path

import numpy as np

# solprop is not a dependency of this package: this tool runs in solprop's environment.
from solprop.physdata import (  # pyright: ignore[reportMissingImports]
    aluminium,
    fused_silica,
    fused_silica_sources,
    glass,
    silver,
)
from solprop.physdata.glasses import (  # pyright: ignore[reportMissingImports]
    _WINDOW_GLASSES,  # pyright: ignore[reportPrivateUsage]
)

C = 299792458.0
OUT = Path(__file__).resolve().parents[1] / "tests" / "data" / "solprop_solids.json"
MATERIALS = {
    "caf2": "CaF2",
    "mgf2": "MgF2",
    "bk7": "N-BK7",
    "sapphire": "Al2O3",
    "baf2": "BaF2",
}
"""solprop's window-glass keys and the optmat material each is."""


def _entry(model: object, material: str, source: str) -> dict[str, object]:
    """Evaluate one solprop model across its range.

    Parameters
    ----------
    model : solprop.physdata.glasses.GlassIndex
        The model.
    material : str
        The optmat material id.
    source : str
        solprop's source name.

    Returns
    -------
    dict
        The fixture entry.
    """
    low, high = model.wavelength_range  # pyright: ignore[reportAttributeAccessIssue]
    wavelength = np.geomspace(low * (1 + 1e-9), high * (1 - 1e-9), 24)
    index = np.asarray(model.index(2 * np.pi * C / wavelength))  # pyright: ignore[reportAttributeAccessIssue]
    return {
        "material": material,
        "source": source,
        "citation": getattr(model, "citation", ""),
        "wavelength_range_m": [low, high],
        "wavelength_m": wavelength.tolist(),
        "n": np.real(index).tolist(),
        "k": np.imag(index).tolist(),
    }


def main() -> None:
    """Evaluate every solprop solid model and write the fixture."""
    solprop = Path.cwd()
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=solprop,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    entries = []
    for key, material in MATERIALS.items():
        for model in _WINDOW_GLASSES[key]:
            entry = _entry(glass(key, source=model.source), material, model.source)
            entry["sellmeier_b_lambda_um"] = [list(p) for p in model.poles]
            entries.append(entry)
    entries.extend(
        _entry(fused_silica(source=source), "SiO2", source)
        for source in fused_silica_sources()
    )
    entries.append(_entry(silver(), "Ag", "stitched"))
    entries.append(_entry(aluminium(), "Al", "stitched"))
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(
        json.dumps({"solprop_commit": commit, "models": entries}, indent=1),
        encoding="utf-8",
    )
    print(f"wrote {len(entries)} models from solprop {commit[:8]} to {OUT}")


if __name__ == "__main__":
    main()
