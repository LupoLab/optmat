# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Record solprop's current third-order coefficients as a test fixture.

    cd ../solprop && uv run python ../optmat/tools/solprop_gamma3.py

runs in *solprop's* environment (this package never imports solprop) and writes
``tests/data/solprop_gamma3.json``: for every gas with a coefficient in
``solprop.physdata.nonlinear_data``, its source and ``gamma3`` in m^5/V^2 (such that
``chi3 = gamma3 * rho``), with the solprop commit it came from.
``tests/test_nonlinear_accuracy_page.py`` compares them with the fits here.
"""

import json
import subprocess
from pathlib import Path

# solprop is not a dependency of this package: this tool runs in solprop's environment.
from solprop.physdata import Gas  # pyright: ignore[reportMissingImports]
from solprop.physdata.nonlinear_data import (  # pyright: ignore[reportMissingImports]
    gamma3,
    gamma3_sources,
)

OUT = Path(__file__).resolve().parents[1] / "tests" / "data" / "solprop_gamma3.json"


def main() -> None:
    """Evaluate every solprop gamma3 and write the fixture."""
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path.cwd(),
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    entries = {}
    for gas in Gas:
        for source in gamma3_sources(gas):
            entries[gas.name] = {"source": source, "gamma3": gamma3(gas, source=source)}
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(
        json.dumps({"solprop_commit": commit, "gamma3": entries}, indent=1) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
