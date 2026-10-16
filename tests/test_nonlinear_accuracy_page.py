# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The third-order pages' tables are what the stored fits and records give today.

``docs/nonlinear_accuracy.md`` tabulates, per gas, the statistical and leave-one-out
uncertainties of the Kerr and THG gamma, the working estimate the page's rule gives, the
validators and solprop's coefficients against the fits; ``docs/nonlinear.md`` tabulates
the fitted parameters. After a refit, regenerate the leave-one-out records with
``tools/fit/nonlinear_robustness.py`` and the solprop fixture with
``tools/solprop_gamma3.py`` (in solprop's environment), and update the tables.
"""

import ast
import json
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pytest

from optmat.constants import EPS0
from optmat.exceptions import ModelRangeWarning
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.nonlinear import (
    AU_HYPERPOLARIZABILITY,
    hyperpolarizability_models,
    load_nonlinear_config,
)

ROOT = Path(__file__).parents[1]
PAGE = (ROOT / "docs" / "nonlinear_accuracy.md").read_text()
MANUAL = (ROOT / "docs/user-guide/nonlinear-models.md").read_text()
MODELS = {gas.value: model for gas, model in hyperpolarizability_models().items()}
# The working estimate is the largest measure rounded up to the next of these (%).
STEPS = (0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50)
POINTS = [["kerr", 400.0], ["kerr", 800.0], ["kerr", 1800.0], ["thg", 800.0]]


def _table(text, name):
    block = text.split(f"<!-- {name}:start -->")[1].split(f"<!-- {name}:end -->")[0]
    return {
        cells[0]: cells
        for line in block.strip().splitlines()[2:]
        for cells in [[c.strip() for c in line.strip("|").split("|")]]
    }


def _percent(cell):
    return float(re.sub(r"[^0-9.eE+-]", "", cell))


def _record(gas):
    path = GAS_DATA_DIRECTORY / "models" / "nonlinear_robustness" / f"{gas}.json"
    return json.loads(path.read_text())


def _stored(gas):
    path = GAS_DATA_DIRECTORY / "models" / "nonlinear" / f"{gas}.json"
    return json.loads(path.read_text())


def _measures(gas, index, floor):
    """Return the statistical, leave-one-out (item, value) and estimate, in %."""
    model = MODELS[gas]
    process, nm = POINTS[index]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        statistical = float(
            model.gamma_sigma_au(process, nm * 1e-9)[0]
            / model.gamma_au(process, nm * 1e-9)
        )
    shifts = {
        key: abs(entry["shift_relative"][index]) * 100
        for key, entry in _record(gas)["dropped"].items()
        if entry["dof"] >= 0
    }
    largest = max(shifts, key=shifts.__getitem__) if shifts else None
    worst = max(statistical * 100, shifts.get(largest, 0.0), floor)
    estimate = next(step for step in STEPS if step >= float(f"{worst:.2g}"))
    return statistical * 100, largest, shifts.get(largest), estimate


def _check(cell, value):
    if value is None:
        assert cell == "none"
    else:
        # Printed to two significant figures.
        assert _percent(cell) == pytest.approx(float(f"{value:.2g}"), rel=1e-9)


AT_800 = _table(PAGE, "nonlinear-accuracy-800")
BANDS = _table(PAGE, "nonlinear-accuracy-bands")


def test_the_tables_cover_every_fit():
    assert set(AT_800) == set(MODELS)
    assert set(BANDS) == set(MODELS)


@pytest.mark.parametrize("gas", sorted(MODELS))
def test_the_robustness_record_is_of_the_current_configuration(gas):
    record = _record(gas)
    assert record["config_sha256"] == load_nonlinear_config(gas).digest
    assert record["points"] == POINTS


@pytest.mark.parametrize("gas", sorted(MODELS))
def test_the_800_nm_row_is_current(gas):
    row = AT_800[gas]
    floor = _percent(row[6])
    statistical, largest, shift, estimate = _measures(gas, 1, floor)
    _check(row[1], statistical)
    _check(row[2], shift)
    assert row[3] == (f"`{largest}`" if largest else "none")
    assert _percent(row[7]) == estimate
    if gas == "O3":  # the harmonic of 800 nm is in the Hartley band
        assert row[4:6] == ["n/a", "n/a"]
        assert row[8] == "n/a (resonant)"
        return
    statistical, _, shift, estimate = _measures(gas, 3, floor)
    _check(row[4], statistical)
    _check(row[5], shift)
    assert _percent(row[8]) == estimate


@pytest.mark.parametrize("gas", sorted(MODELS))
@pytest.mark.parametrize(("index", "column"), [(0, 1), (2, 4)])
def test_the_band_rows_are_current(gas, index, column):
    cells = BANDS[gas][column : column + 3]
    low = MODELS[gas].wavelength_range[0] * 1e9
    if POINTS[index][1] < low * (1 - 1e-9):
        assert cells == ["n/a"] * 3
        return
    floor = _percent(AT_800[gas][6])
    statistical, _, shift, estimate = _measures(gas, index, floor)
    _check(cells[0], statistical)
    _check(cells[1], shift)
    assert _percent(cells[2]) == estimate


def test_the_gallery_estimates_are_the_pages():
    source = (
        ROOT / "docs" / "gallery" / "nonlinear" / "plot_nonlinear_accuracy.py"
    ).read_text()
    node = next(
        n
        for n in ast.parse(source).body
        if isinstance(n, ast.Assign)
        and getattr(n.targets[0], "id", None) == "ESTIMATES"
    )
    estimates = eval(  # a literal dict with np.nan
        compile(ast.Expression(node.value), "gallery", "eval"), {"np": np}
    )
    for gas, (at_400, at_800, at_1800) in estimates.items():
        if math.isnan(at_400):
            assert BANDS[gas][3] == "n/a"
        else:
            assert _percent(BANDS[gas][3]) == at_400
        assert _percent(AT_800[gas][7]) == at_800
        assert _percent(BANDS[gas][6]) == at_1800


def test_the_validator_ratios_are_current():
    rows = _table(PAGE, "nonlinear-validators")
    printed = {}
    for cells in rows.values():
        source = cells[1].strip("`")
        for entry in cells[2].split("; "):
            match = re.fullmatch(r"(\w+)(?: \((\w+)\))? ([0-9.]+)", entry)
            if match:
                printed[(source, *match.groups()[:2])] = float(match.group(3))
    compared = 0
    for gas in MODELS:
        validation = _stored(gas)["validation"]
        for key, entry in validation.items():
            if "note" in entry:
                continue
            source, dataset = key.split("/")
            same = sum(1 for k in validation if k.startswith(source + "/"))
            value = printed[(source, gas, dataset if same > 1 else None)]
            assert value == pytest.approx(entry["mean_ratio"], abs=0.0051), key
            compared += 1
    assert compared == len(printed)


def test_solprops_ratios_are_current():
    rows = _table(
        (ROOT / "tests/data/historical_documentation_tables.md").read_text(),
        "nonlinear-solprop",
    )
    fixture = json.loads((ROOT / "tests" / "data" / "solprop_gamma3.json").read_text())
    solprop = fixture["gamma3"]
    assert set(rows) == set(MODELS)
    for gas, cells in rows.items():
        if gas not in solprop:
            assert cells[1] == "none"
            assert cells[2:] == ["n/a"] * 5
            continue
        gamma = 6 * EPS0 * solprop[gas]["gamma3"] / AU_HYPERPOLARIZABILITY
        assert cells[1] == solprop[gas]["source"]
        assert float(cells[2]) == pytest.approx(gamma, rel=5e-4)
        for column, (process, nm) in enumerate(POINTS, start=3):
            with warnings.catch_warnings():
                # N2O's 400 nm and THG values extrapolate its borrowed dispersion,
                # as the page says.
                warnings.simplefilter("ignore", ModelRangeWarning)
                model = float(MODELS[gas].gamma_au(process, nm * 1e-9))
            ratio = gamma / model
            # Printed to 0.001.
            assert float(cells[column]) == pytest.approx(ratio, abs=5.1e-4)


def test_the_manual_fit_table_is_current():
    rows = _table(MANUAL, "nonlinear-fits")
    assert set(rows) == set(MODELS)
    for gas, cells in rows.items():
        record = _stored(gas)
        parameters = list(record["parameters"].values())
        assert float(cells[1]) == pytest.approx(parameters[0], rel=6e-4)
        for column, value in zip(cells[2:5], parameters[1:], strict=True):
            assert float(column) == pytest.approx(value, abs=5.1e-4)
        assert int(cells[6].split("/")[1]) == record["dof"]
        w2 = record["max_w2_fitted_cm-2"]
        assert float(cells[7]) == pytest.approx(w2 / 1e9, rel=6e-3)
