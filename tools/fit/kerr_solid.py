# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Run a solid's Kerr fit and store the result.

    uv run python tools/fit/kerr_solid.py <name>

reads ``data/solid/kerr/<name>.toml``, fits n2(lambda), the Raman fraction and the
calibration scales, prints a summary and writes ``data/solid/kerr/<name>.json``
(``optmat.solid.kerr``).
"""

import sys

import numpy as np

from optmat.solid.kerr import KERR_DIRECTORY, dumps, fit, load_config


def main(name: str) -> None:
    """Fit, report and write one configuration."""
    config = load_config(name)
    result = fit(config)
    model = result.model
    print(f"{config.name}: chi2 = {result.chi2:.2f} for {result.dof} dof")
    for key, (n, chi2) in result.dataset_chi2.items():
        print(f"  {key:45s} n = {n:2d}  chi2 = {chi2:7.2f}")
    for group, scale in model.scales.items():
        print(f"  scale {group:25s} {scale:.3f}")
    sigma = np.sqrt(np.diag(model.covariance))
    for label, value, s in zip(
        ("N", "a", "b", "f_R"), model.parameters, sigma, strict=True
    ):
        print(f"  {label:4s} = {value:.6g} +- {s:.3g}")
    for wavelength in (400e-9, 800e-9, 1030e-9, 1550e-9, 2000e-9):
        print(
            f"  n2({wavelength * 1e9:.0f} nm) = {model.n2(wavelength):.4g} "
            f"+- {model.n2_sigma(wavelength):.2g} m2/W"
        )
    print("  validators (value / model):")
    for key, ratios in result.validation.items():
        print(f"    {key:45s} {' '.join(f'{r:.3f}' for r in ratios)}")
    out = KERR_DIRECTORY / f"{config.name}.json"
    out.write_text(dumps(result))
    print(f"wrote {out}")


if __name__ == "__main__":
    main(sys.argv[1])
