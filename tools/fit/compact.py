# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Derive a compact model from its stored parent fit and store it (OEP 2).

    uv run python tools/fit/compact.py <gas> [--fewer]

reads ``data/gas/compact/<gas>.toml`` and the stored ``fit/<gas>``, derives the compact
vibrational poles, prints the residuals and the reproduction of the parent's data and
writes ``data/gas/models/compact/<gas>.json``, which the package loads as
``compact/<gas>`` (``optmat.gas.models.compact``). ``--fewer`` also reports the
reduction with one pole fewer in each searched region, which must fail. A test checks
the stored model against its parent and configuration; the ``fits`` lane re-derives it.
Run it after every refit of the parent.
"""

import dataclasses
import json
import sys
import warnings

from optmat.datasets import DATA_DIRECTORY
from optmat.gas.models.compact import (
    COMPACT_MODEL_DIRECTORY,
    compact_record,
    derive,
    load_compact_config,
    reproduction,
)
from optmat.gas.models.fitted import fitted_models

HARTREE_CM = 219474.6313632


def main(gas: str, *, fewer: bool = False) -> None:
    """Derive, report and write one compact model."""
    config = load_compact_config(gas)
    parent = fitted_models()[config.parent]
    reduction = derive(parent, config)
    worst = ", ".join(f"{k} {v:.3g}" for k, v in reduction.worst.items())
    count = len(reduction.energies)
    print(f"compact/{gas}: {count} poles, worst |residual|/sigma {worst}")
    tolerance = config.tolerance_sigma
    print(
        f"  passed {reduction.passed} (tolerance {tolerance}); start {reduction.start}"
    )
    m0, m2 = reduction.moment_errors
    print(f"  moments M0, M-2 relative errors {m0:.2g}, {m2:.2g}")
    for f, e, region, kept in zip(
        reduction.strengths,
        reduction.energies,
        reduction.region,
        reduction.kept,
        strict=True,
    ):
        tag = " (kept)" if kept else ""
        wavenumber = e * HARTREE_CM
        print(f"  region {region}: F = {f:.10e}  E = {e:.10e} hartree", end="")
        print(f" ({wavenumber:.2f} cm-1){tag}")
    if fewer:
        for index, count in enumerate(config.poles):
            if count == 0:
                continue
            poles = list(config.poles)
            poles[index] -= 1
            try:
                trial = derive(parent, dataclasses.replace(config, poles=tuple(poles)))
            except Exception as error:
                print(f"  one fewer in region {index}: {error}")
                continue
            worst = ", ".join(f"{k} {v:.3g}" for k, v in trial.worst.items())
            passed = trial.passed
            print(f"  one fewer in region {index} ({poles}): {worst}; passed {passed}")
    out = COMPACT_MODEL_DIRECTORY / f"{gas}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(compact_record(parent, config, reduction), indent=1) + "\n"
    )
    print(f"wrote {out.relative_to(DATA_DIRECTORY)}")
    # The reproduction needs the stored model; load it afresh.
    from optmat.gas.models import compact  # noqa: PLC0415

    compact.compact_models.cache_clear()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = reproduction(compact.compact_models()[f"compact/{gas}"])
    print("  reproduction (change in each datum's sigma):")
    for key, entry in report.items():
        if "static_relative" in entry:
            print(f"    {key:50s} static, relative {entry['static_relative']:.2g}")
            continue
        beyond = entry.get("beyond_max_sigma")
        extra = f"  beyond {beyond:.3g}" if beyond is not None else ""
        print(
            f"    {key:50s} n = {int(entry['n']):3d}  max {entry['max_sigma']:.3g}  "
            f"chi2 change {entry['chi2_change']:+.3g}{extra}"
        )
    if not reduction.passed:
        sys.exit(f"compact/{gas} misses the tolerance: give a region more poles")


if __name__ == "__main__":
    main(sys.argv[1], fewer="--fewer" in sys.argv[2:])
