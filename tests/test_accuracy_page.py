# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The accuracy page's table is what the stored fits and robustness records give today.

Regenerate the numbers with the gallery page *Dispersion uncertainty* (which prints
them) after refitting, and the leave-one-out records with ``tools/fit/robustness.py``.
"""

import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from optmat.conditions import State
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.fitting import HARTREE_NM
from optmat.gas.models.fitted import fitted_models

ROOT = Path(__file__).parents[1]
PAGE = (ROOT / "docs" / "accuracy.md").read_text(encoding="utf-8")
TABLE = PAGE.split("<!-- accuracy-table:start -->")[1].split(
    "<!-- accuracy-table:end -->"
)[0]
ROWS = {
    cells[0]: cells
    for line in TABLE.strip().splitlines()[2:]
    for cells in [[c.strip() for c in line.strip("|").split("|")]]
}


def _percent(cell):
    return float(re.sub(r"[^0-9.+-]", "", cell)) / 100


def _record(gas):
    path = GAS_DATA_DIRECTORY / "models" / "robustness" / f"{gas}.json"
    return json.loads(path.read_text(encoding="utf-8"))


def test_the_table_covers_every_fit():
    assert set(ROWS) == {name.split("/")[1] for name in fitted_models()}


@pytest.mark.parametrize("gas", sorted(ROWS))
def test_statistical_uncertainty_is_current(gas):
    fit = fitted_models()[f"fit/{gas}"]
    sigma = fit.dispersion_sigma_relative(np.array([HARTREE_NM / 800.0]), 2)[0]
    # Printed to two significant figures.
    assert _percent(ROWS[gas][1]) == pytest.approx(sigma, rel=0.05, abs=0.0)


@pytest.mark.parametrize("gas", sorted(ROWS))
def test_leave_one_out_shift_is_current(gas):
    record = _record(gas)
    column = record["wavelength_nm"].index(800.0)
    shifts = {
        key: abs(entry["gvd_shift_relative"][column])
        for key, entry in record["dropped"].items()
    }
    largest = max(shifts, key=shifts.__getitem__)
    assert ROWS[gas][3] == largest
    # Printed to two significant figures.
    assert _percent(ROWS[gas][2]) == pytest.approx(shifts[largest], rel=0.05, abs=0.0)


@pytest.mark.parametrize("gas", sorted(ROWS))
def test_solprop_formula_deviation_is_current(gas):
    historical = (ROOT / "tests/data/historical_documentation_tables.md").read_text(
        encoding="utf-8"
    )
    block = historical.split("<!-- accuracy-table:start -->")[1].split(
        "<!-- accuracy-table:end -->"
    )[0]
    rows = {
        cells[0]: cells
        for line in block.strip().splitlines()[2:]
        for cells in [[c.strip() for c in line.strip("|").split("|")]]
    }
    sys.path.insert(0, str(ROOT / "docs" / "gallery"))
    try:
        from _helpers import (  # noqa: PLC0415  # pyright: ignore[reportMissingImports]
            group_velocity_dispersion,
            solprop_defaults,
        )
    finally:
        sys.path.pop(0)
    state = State(1e5, 293.15)
    wavelength = np.array([800e-9])
    fit = fitted_models()[f"fit/{gas}"]
    if gas not in solprop_defaults():  # CO2: solprop has no formula for it
        assert rows[gas][5] == "n/a"
        return
    formula = solprop_defaults()[gas]
    ratio = (
        group_velocity_dispersion(formula, wavelength, state)
        / group_velocity_dispersion(fit, wavelength, state)
    )[0]
    # Printed to 0.01 %.
    assert _percent(rows[gas][5]) == pytest.approx(ratio - 1, abs=6e-5)


BANDS = PAGE.split("<!-- accuracy-bands:start -->")[1].split(
    "<!-- accuracy-bands:end -->"
)[0]
BAND_ROWS = {
    cells[0]: cells
    for line in BANDS.strip().splitlines()[2:]
    for cells in [[c.strip() for c in line.strip("|").split("|")]]
}


def test_the_band_table_covers_every_fit():
    assert set(BAND_ROWS) == {name.split("/")[1] for name in fitted_models()}


@pytest.mark.parametrize("gas", sorted(BAND_ROWS))
@pytest.mark.parametrize(("wavelength", "column"), [(343.0, 1), (1800.0, 4)])
def test_band_uncertainties_are_current(gas, wavelength, column):
    record = _record(gas)
    cells = BAND_ROWS[gas]
    if wavelength not in record["wavelength_nm"]:
        # Outside the model's declared range: the page says so.
        assert cells[column : column + 3] == ["n/a"] * 3
        return
    fit = fitted_models()[f"fit/{gas}"]
    sigma = fit.dispersion_sigma_relative(np.array([HARTREE_NM / wavelength]), 2)[0]
    index = record["wavelength_nm"].index(wavelength)
    largest = max(
        abs(entry["gvd_shift_relative"][index]) for entry in record["dropped"].values()
    )
    # Printed to two significant figures.
    assert _percent(cells[column]) == pytest.approx(sigma, rel=0.05, abs=0.0)
    assert _percent(cells[column + 1]) == pytest.approx(largest, rel=0.05, abs=0.0)
