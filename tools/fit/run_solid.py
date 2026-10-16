# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Run a solid fit configuration and store the result.

    uv run python tools/fit/run_solid.py <name>

reads ``data/solid/fits/<name>.toml``, fits the universal dispersion model's free
parameters to the configured datasets, prints a summary and writes
``data/solid/models/fitted/<name>.json``, which the package loads as the model
``fit/<name>`` (``optmat.solid.models.fitted``).
"""

import sys

import numpy as np

from optmat.datasets import DATA_DIRECTORY
from optmat.solid import SOLID_DATA_DIRECTORY
from optmat.solid.fitting import dumps, fit, load_config, validation


def main(name: str) -> None:
    """Fit, report and write one configuration."""
    config = load_config(name)
    result = fit(config)
    print(
        f"{config.name}: chi2 = {result.chi2:.2f} for {result.dof} dof "
        f"(reduced {result.chi2 / result.dof:.3f})"
    )
    for key, (n, chi2) in result.dataset_chi2.items():
        print(f"  {key:50s} n = {n:4d}  chi2 = {chi2:10.2f}")
    for key, offset in result.offsets.items():
        print(f"  offset {key:43s} {offset:+.3e}")
    sigma = np.sqrt(np.diag(result.covariance))
    for name_, s in zip(config.free, sigma, strict=False):
        print(f"  {name_:10s} = {result.values[name_]:.10g} +- {s:.3g}")
    print("  validators:")
    for key, entry in validation(result).items():
        print(
            f"    {key:50s} n = {int(entry['n']):4d}  chi2 = {entry['chi2']:10.2f}  "
            f"mean deviation = {entry['mean_deviation']:+.2e}"
        )
    out = SOLID_DATA_DIRECTORY / "models" / "fitted" / f"{config.name}.json"
    out.write_text(dumps(result))
    print(f"wrote {out.relative_to(DATA_DIRECTORY)}")


if __name__ == "__main__":
    main(sys.argv[1])
