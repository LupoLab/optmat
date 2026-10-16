# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Leave one dataset (or calculation) out of a hyperpolarizability fit.

    uv run python tools/fit/nonlinear_robustness.py <gas> [<gas> ...]

refits ``data/gas/nonlinear/<gas>.toml`` once without each fitted dataset, or without
each scale group (a calculation's processes leave together), and writes
``data/gas/models/nonlinear_robustness/<gas>.json``: per dropped item, the relative
change of the Kerr gamma at 400, 800 and 1800 nm and of the THG gamma at 800 nm, with
the refit's degrees of freedom. A gas with a single fitted item has no entries. A refit
with fewer values than parameters (negative degrees of freedom) is decided by the priors
alone and is not a measure of the uncertainty. The page *Third-order accuracy* tabulates
them.
"""

import dataclasses
import json
import sys
import warnings

import numpy as np

from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.nonlinear import (
    HyperpolarizabilityModel,
    fit_hyperpolarizability,
    load_nonlinear_config,
)

OUT = GAS_DATA_DIRECTORY / "models" / "nonlinear_robustness"
POINTS = (("kerr", 400.0), ("kerr", 800.0), ("kerr", 1800.0), ("thg", 800.0))


def _values(parameters: np.ndarray) -> np.ndarray:
    """Evaluate the probe points of a parameter vector."""
    model = HyperpolarizabilityModel(
        gas=load_nonlinear_config("He").gas,
        parameters=tuple(float(p) for p in parameters),  # pyright: ignore[reportArgumentType]
        covariance=np.zeros((4, 4)),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.array([float(model.gamma_au(p, w * 1e-9)) for p, w in POINTS])


def main(gas: str) -> None:
    """Refit one gas without each item and write the shifts."""
    config = load_nonlinear_config(gas)
    full = fit_hyperpolarizability(config)
    base = _values(full.parameters)
    items: dict[str, list] = {}
    for options in config.fitted:
        items.setdefault(options.scale_group or options.key, []).append(options)
    dropped = {}
    if len(items) > 1:
        for name, members in items.items():
            kept = tuple(o for o in config.fitted if o not in members)
            refit = fit_hyperpolarizability(dataclasses.replace(config, fitted=kept))
            shift = _values(refit.parameters) / base - 1
            dropped[name] = {
                "shift_relative": shift.tolist(),
                "reduced_chi2": refit.chi2 / max(refit.dof, 1),
                "dof": refit.dof,
            }
    record = {
        "gas": gas,
        "config_sha256": config.digest,
        "points": [[p, w] for p, w in POINTS],
        "dropped": dropped,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{gas}.json").write_text(json.dumps(record, indent=1) + "\n")
    worst = max(
        (max(abs(x) for x in e["shift_relative"]) for e in dropped.values()),
        default=0.0,
    )
    print(f"{gas}: {len(dropped)} items, largest shift {worst * 100:.3g} %")


if __name__ == "__main__":
    for name in sys.argv[1:]:
        main(name)
