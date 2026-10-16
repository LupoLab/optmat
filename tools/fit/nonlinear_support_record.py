# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Derive and store the hyperpolarizability support records (OEP 5).

    uv run python tools/fit/nonlinear_support_record.py [<gas> ...]

reads each gas's stored fit, ``data/gas/nonlinear_support.toml`` and the curated
datasets (``optmat.gas._nonlinear_support``), prints a summary and writes
``data/gas/models/nonlinear_support/<gas>.json``; with no gas, every gas with a stored
fit. Re-run it after a refit, or a change to the selection or a selected dataset: a test
re-derives every stored record and compares.
"""

import sys

from optmat.datasets import DATA_DIRECTORY
from optmat.gas._nonlinear_support import derive_support
from optmat.gas.nonlinear import (
    _SUPPORT_DIRECTORY,
    SupportLimit,
    _support_dumps,
    hyperpolarizability_models,
)


def _nm(limit: SupportLimit | None, factor: float = 1.0) -> str:
    """Return a limit as the wavelength of one leg (times a factor), nm, or a dash."""
    if limit is None:
        return "-"
    return f"{factor * 1e9 / limit.wavenumber:.1f} ({limit.kind})"


def main(gas: str) -> None:
    """Derive, report and write one gas's record."""
    record = derive_support(gas)
    edges = [record.band_edge(p) for p in ("linear", "any")]
    print(
        f"{gas}: {record.fitted_points} points; largest leg "
        f"{_nm(record.leg_max)} from {', '.join(record.leg_max.datasets)}; "
        f"pair {_nm(record.pair_max, 2)}"
    )
    print(
        f"  E1 {_nm(record.e1)}; E2 {_nm(record.e2, 2)}; E2 any "
        f"{_nm(record.e2_any, 2)}; band edge "
        + " / ".join("-" if e is None else f"{1e9 / e:.1f} nm" for e in edges)
    )
    for name, reason in record.missing_reasons.items():
        print(f"  {name} missing: {reason}")
    out = _SUPPORT_DIRECTORY / f"{gas}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(_support_dumps(record), encoding="utf-8")
    print(f"  wrote {out.relative_to(DATA_DIRECTORY)} ({record.digest[:12]})")


if __name__ == "__main__":
    names = sys.argv[1:] or sorted(g.value for g in hyperpolarizability_models())
    for name in names:
        main(name)
