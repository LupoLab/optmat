# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Run a fit configuration and store the result.

    uv run python tools/fit/run.py <name>

reads ``data/gas/fits/<name>.toml``, fits, prints a summary and writes
``data/gas/models/fitted/<name>.json``, which the package loads as a model
(``optmat.gas.models.fitted``). A test re-runs every stored fit and compares.
"""

import sys

import numpy as np

from optmat.datasets import DATA_DIRECTORY
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.fitting import dumps, fit, load_config, validation


def main(name: str) -> None:
    """Fit, report and write one configuration."""
    config = load_config(name)
    result = fit(config)
    print(
        f"{config.name}: chi2 = {result.chi2:.2f} for {result.dof} dof "
        f"(reduced {result.reduced_chi2:.3f})"
    )
    for key, (n, chi2) in result.dataset_chi2.items():
        scale = result.scales[key]
        print(f"  {key:45s} n = {n:3d}  chi2 = {chi2:9.2f}  scale = {scale:.6f}")
    sigma = np.sqrt(np.diag(result.covariance))
    model = result.model
    for label, f, e in zip(model.labels, model.strengths, model.energies, strict=True):
        print(f"  {label:20s} f = {f:.6g}  E = {e:.6g} hartree ({45.56335 / e:.3f} nm)")
    print("  parameter sigmas:", " ".join(f"{s:.2g}" for s in sigma))
    print("  validators:")
    for key, entry in validation(result).items():
        if "note" in entry:
            print(f"    {key:45s} {entry['note']}")
        else:
            print(
                f"    {key:45s} n = {entry['n']:3d}  chi2 = {entry['chi2']:9.2f}  "
                "mean deviation = "
                f"{float(entry['mean_relative_deviation']) * 1e6:9.1f} ppm"
            )
    out = GAS_DATA_DIRECTORY / "models" / "fitted" / f"{config.name}.json"
    out.write_text(dumps(result))
    print(f"wrote {out.relative_to(DATA_DIRECTORY)}")


if __name__ == "__main__":
    main(sys.argv[1])
