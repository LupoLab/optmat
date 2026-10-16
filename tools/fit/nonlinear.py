# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Fit a gas's second hyperpolarizability and store the model.

    uv run python tools/fit/nonlinear.py <gas> [<gas> ...]

reads ``data/gas/nonlinear/<gas>.toml``, fits the dispersion model
(``optmat.gas.nonlinear``), prints a summary and writes
``data/gas/models/nonlinear/<gas>.json``. A test re-runs every stored fit and compares.
"""

import json
import sys

import numpy as np

from optmat.datasets import DATA_DIRECTORY
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.nonlinear import dumps, fit_hyperpolarizability, load_nonlinear_config


def main(gas: str) -> None:
    """Fit, report and write one gas."""
    config = load_nonlinear_config(gas)
    result = fit_hyperpolarizability(config)
    text = dumps(result)
    record = json.loads(text)
    sigma = np.sqrt(np.diag(result.covariance))
    names = ("gamma0 (au)", "A (1e-10 cm2)", "B (1e-20 cm4)", "B' (1e-20 cm4)")
    print(f"{gas}: chi2 = {result.chi2:.2f} for {result.dof} dof")
    for name, value, error in zip(names, result.parameters, sigma, strict=True):
        print(f"  {name:16s} {value:12.5g} +- {error:.3g}")
    for key, (n, chi2) in result.dataset_chi2.items():
        scale = result.scales[key]
        print(f"  {key:50s} n = {n:3d} chi2 = {chi2:9.2f} scale = {scale:.4f}")
    print("  validators:")
    for key, entry in record["validation"].items():
        if "note" in entry:
            print(f"    {key:50s} {entry['note']}")
        else:
            print(
                f"    {key:50s} n = {entry['n']:3d} data/model = "
                f"{entry['mean_ratio']:.4f} chi2 = {entry['chi2']:.2f}"
            )
    out = GAS_DATA_DIRECTORY / "models" / "nonlinear" / f"{gas}.json"
    out.write_text(text, encoding="utf-8")
    print(f"wrote {out.relative_to(DATA_DIRECTORY)}")


if __name__ == "__main__":
    for name in sys.argv[1:]:
        main(name)
