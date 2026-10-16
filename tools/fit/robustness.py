# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Leave one dataset out of a fit and record how far the model moves.

    uv run python tools/fit/robustness.py <name> [<name> ...]

refits ``data/gas/fits/<name>.toml`` once without each of its fitted datasets and writes
``data/gas/models/robustness/<name>.json``: for each dropped dataset, the relative
change of the polarizability at a few wavelengths across the model's range, in units of
the full fit's one-sigma uncertainty, the relative change of the group-velocity
dispersion (``d^2 (omega alpha) / d omega^2``, which ``beta_2`` is proportional to), and
the refit's reduced chi-squared. A dataset whose removal moves the model by more than
its uncertainty carries the fit there alone (or pulls it against the others). The
gallery page *Robustness* plots the files.
"""

import dataclasses
import json
import sys

import numpy as np

from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.fitting import HARTREE_NM, fit, load_config
from optmat.gas.models.fitted import fitted_models

OUT = GAS_DATA_DIRECTORY / "models" / "robustness"


def main(name: str) -> None:
    """Refit one configuration without each dataset and write the shifts."""
    config = load_config(name)
    full = fitted_models()[f"fit/{name}"]
    if full.record["config_sha256"] != config.digest:
        msg = f"{name}: the stored fit is stale; run tools/fit/run.py {name} first"
        raise SystemExit(msg)
    # Each refit starts from the full solution and follows it to the nearest minimum
    # (starts=0), so that a shift measures the dataset's influence and not a jump
    # between the near-degenerate minima the multi-start search can land in (TC-101).
    full_result = fit(config)
    shortest = config.range_nm[0]
    longest = config.range_nm[1]
    # The accuracy page's wavelengths (343, 800, 1800 nm) and a few across the range;
    # those outside the model's declared range are left out.
    wanted = np.array(
        [1.05 * shortest, 2 * shortest, 343.0, 400.0, 800.0, 1600.0, 1800.0]
    )
    wavelength_nm = np.unique(
        wanted[(wanted >= 1.05 * shortest) & (wanted <= 0.999 * longest)]
    )
    omega = HARTREE_NM / wavelength_nm
    alpha = full.polarizability_au(omega)
    sigma = full.polarizability_sigma_au(omega)
    gvd = full.dispersion_au(omega, 2)
    dropped = {}
    for options in config.fitted:
        reduced = dataclasses.replace(
            config, fitted=tuple(o for o in config.fitted if o is not options)
        )
        result = fit(reduced, starts=0, start=full_result)
        shift = result.model.polarizability_au(omega) - alpha
        dropped[options.key] = {
            "shift_relative": (shift / alpha).tolist(),
            "shift_sigma": (shift / sigma).tolist(),
            "gvd_shift_relative": (
                result.model.dispersion_au(omega, 2) / gvd - 1
            ).tolist(),
            "reduced_chi2": result.reduced_chi2,
        }
        worst = float(np.max(np.abs(shift / sigma)))
        print(f"  {options.key:45s} max |shift| = {worst:6.2f} sigma")
    OUT.mkdir(exist_ok=True)
    record = {
        "name": name,
        "config_sha256": config.digest,
        "wavelength_nm": wavelength_nm.tolist(),
        "sigma_relative": (sigma / alpha).tolist(),
        "gvd_sigma_relative": full.dispersion_sigma_relative(omega, 2).tolist(),
        "dropped": dropped,
    }
    (OUT / f"{name}.json").write_text(json.dumps(record, indent=1))
    print(f"wrote models/robustness/{name}.json")


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(arg)
        main(arg)
