# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Derive and store the two-photon ionisation models (OEP 6, Phase 1).

    uv run python tools/fit/two_photon_record.py [<gas> ...]

reads ``data/gas/two_photon.toml`` and the curated datasets
(``optmat.gas._two_photon_derivation``), prints a summary and writes
``data/gas/models/two_photon/<gas>_<polarization>.json`` for every polarisation the
selection approves; with no gas, every gas. Re-run it after a change to the selection
or to a selected dataset: a test re-derives every stored model and compares.
"""

import math
import sys

from optmat.constants import C
from optmat.datasets import DATA_DIRECTORY
from optmat.gas._two_photon_derivation import derive
from optmat.gas.two_photon import _DIRECTORY, POLARIZATIONS, _dumps, _selection


def _nm(omega: float) -> str:
    """Return an angular frequency as a vacuum wavelength, nm."""
    return f"{2 * math.pi * C / omega * 1e9:.2f}"


def main(gas: str) -> None:
    """Derive, report and write one gas's models."""
    entry = _selection()[gas]
    if "missing" in entry:
        print(f"{gas}: no model: {' '.join(entry['missing'].split())[:100]}...")
        return
    for polarization in (p for p in POLARIZATIONS if p in entry):
        model = derive(gas, polarization)
        print(f"{gas} {polarization}: {len(model.omega)} vertices, {model.status}")
        for r in model.regions:
            upper = "inf" if math.isinf(r.upper) else _nm(r.upper)
            lower = "inf" if r.lower == 0 else _nm(r.lower)
            print(f"  {lower:>8} to {upper:>8} nm  {r.label}")
        for s in model.uncertainty.spreads:
            print(
                f"  {s.label}: {s.lower_ratio:.3f} to {s.upper_ratio:.3f} "
                f"({s.points} points, {s.role})"
            )
        out = _DIRECTORY / f"{gas}_{polarization}.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(_dumps(model), encoding="utf-8")
        print(f"  wrote {out.relative_to(DATA_DIRECTORY)} ({model.digest[:12]})")


if __name__ == "__main__":
    for name in sys.argv[1:] or list(_selection()):
        main(name)
