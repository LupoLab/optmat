# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Derive and store the bound two-photon level records (OEP 6, Phase 2a).

    uv run python tools/fit/two_photon_bound_record.py [<gas> ...]

reads ``data/gas/two_photon_bound.toml``, the NIST level lists held as returned and the
curated datasets (``optmat.gas._two_photon_bound_derivation``), prints a summary and
writes ``data/gas/models/two_photon_bound/<gas>.json``; with no gas, every noble gas.
Re-run it after a change to the selection or to a selected dataset: a test re-derives
every stored record and compares.
"""

import sys

from optmat.datasets import DATA_DIRECTORY
from optmat.gas._two_photon_bound_derivation import derive_bound
from optmat.gas.two_photon_bound import _DIRECTORY, _GASES, _dumps


def main(gas: str) -> None:
    """Derive, report and write one gas's record."""
    record = derive_bound(gas)
    linear = sum(level.reach == "linear" for level in record.levels)
    print(
        f"{gas}: {len(record.levels)} levels ({linear} J = 0, 2), "
        f"{len(record.strengths)} strengths, {len(record.photoionisation)} "
        f"photoionisation sets; {record.assessment}"
    )
    for s in record.strengths:
        where = s.level.split("|", 1)[-1] if s.level else "unassigned"
        delta = "-" if s.delta is None else f"{s.delta:.4g}"
        sigma0 = "-" if s.sigma0 is None else f"{s.sigma0:.3e}"
        print(f"  {s.dataset:38s} {where:14s} {s.kind:10s}", end=" ")
        print(f"delta {delta:>9} sigma0 {sigma0}")
    out = _DIRECTORY / f"{gas}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(_dumps(record), encoding="utf-8")
    print(f"  wrote {out.relative_to(DATA_DIRECTORY)} ({record.digest[:12]})")


if __name__ == "__main__":
    for name in sys.argv[1:] or _GASES:
        main(name)
