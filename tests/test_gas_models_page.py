# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The table of ``docs/gas_models.md`` is what the stored fits give today."""

import re
from pathlib import Path
from typing import cast

import numpy as np
import pytest

from optmat.constants import AU_POLARIZABILITY
from optmat.gas.fitting import HARTREE_NM
from optmat.gas.models.fitted import fitted_models

PAGE = (Path(__file__).parents[1] / "docs" / "gas_models.md").read_text()
TABLE = PAGE.split("<!-- gas-models:start -->")[1].split("<!-- gas-models:end -->")[0]
ROWS = {
    cells[0]: cells
    for line in TABLE.strip().splitlines()[2:]
    for cells in [[c.strip() for c in line.strip("|").split("|")]]
}


def test_the_table_covers_every_fit():
    assert set(ROWS) == {name.split("/")[1] for name in fitted_models()}


@pytest.mark.parametrize("gas", sorted(ROWS))
def test_the_row_is_current(gas):
    model = fitted_models()[f"fit/{gas}"]
    cells = ROWS[gas]
    (low, high), *rest = model.intervals
    upper = "static" if np.isinf(high) else f"{high * 1e9:.0f}"
    # A fit with OEP 3's windows: its first interval and how many follow, to where.
    more = (
        f" + {len(rest)} intervals to {model.intervals[-1][1] * 1e9:.0f}"
        if rest
        else ""
    )
    assert cells[2] == f"{low * 1e9:.0f}-{upper}{more}"
    alpha = float(np.real(model.polarizability(np.array([800e-9]))[0]))
    sigma = (
        float(model.polarizability_sigma_au(np.array([HARTREE_NM / 800.0]))[0])
        * AU_POLARIZABILITY
        / alpha
    )
    number = float(re.sub(r"[^0-9.e+-]", "", cells[3]))
    printed = number * (1e-6 if "ppm" in cells[3] else 1e-2)
    # Printed to two significant figures.
    assert printed == pytest.approx(float(f"{sigma:.2g}"), rel=1e-9)
    assert float(cells[4]) == pytest.approx(model.reduced_chi2, abs=0.0051)


def test_the_molecules_intervals_are_the_derived_ones():
    from optmat.gas.infrared import infrared_windows, qualified_windows  # noqa: PLC0415

    table = PAGE.split("<!-- molecule-intervals:start -->")[1].split(
        "<!-- molecule-intervals:end -->"
    )[0]
    printed = [
        [c.strip() for c in line.strip("|").split("|")]
        for line in table.strip().splitlines()[2:]
    ]
    expected = []
    for gas in ("CO", "CO2", "CH4", "N2O"):
        intervals, reasons = qualified_windows(gas)
        model = fitted_models()[f"fit/{gas}"]
        # The page's intervals are the stored model's.
        assert [
            list(i) for i in cast("list[list[float]]", model.record["intervals_nm"])
        ] == [list(i) for i in intervals]
        windows = [r for r in infrared_windows(gas) if r.kind == "window"]
        for k, (low, high) in enumerate(intervals):
            low_cm, high_cm = 1e7 / high, 1e7 / low
            over = [
                r
                for r in windows
                if r.high_cm > low_cm + 1e-6 and r.low_cm < high_cm - 1e-6
            ]
            peak = max(float(r.details["max_lorentz_peak_m"]) for r in over)
            expected.append(
                [
                    gas,
                    f"{low / 1e3:.4g}-{high / 1e3:.4g}",
                    f"{max(r.max_absorption_m for r in over):.2g}",
                    f"{peak:.2g}",
                    reasons[k] if k < len(reasons) else "\N{EM DASH}",
                ]
            )
    assert printed == expected
