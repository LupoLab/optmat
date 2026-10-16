# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The rates of Luna.jl and solprop as the fixture records them (TC-815)."""

import json
import warnings
from pathlib import Path

import numpy as np
import pytest

from optmat.exceptions import ModelRangeWarning
from optmat.gas.ionisation_rates import ionisation_rate_model

FIXTURE = Path(__file__).parent / "data" / "legacy_ionisation_rates.json"
GASES = ("He", "Ne", "Ar", "Kr", "Xe")


@pytest.fixture(scope="module")
def legacy():
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def _rates(section, gas, name="rate"):
    return np.array(
        [np.nan if r is None else r for r in section["gases"][gas][name]], dtype=float
    )


def test_the_codes_are_recorded(legacy):
    assert legacy["luna"]["commit"].startswith("de9db7c")
    assert legacy["solprop"]["commit"].startswith("05172a4")
    assert legacy["luna"]["clean"]
    assert legacy["solprop"]["clean"]


@pytest.mark.parametrize("gas", GASES)
def test_solprop_is_luna_with_tong_and_lins_factor(legacy, gas):
    luna = _rates(legacy["luna"], gas)
    kept = np.isfinite(luna)
    without = _rates(legacy["solprop"], gas, "rate_bsi_none")
    assert without[kept] == pytest.approx(luna[kept], rel=5e-7, abs=0.0)
    # The factor only lowers the rate.
    assert np.all(_rates(legacy["solprop"], gas)[kept] < luna[kept])


@pytest.mark.parametrize("gas", GASES)
def test_the_recorded_optmat_rates_are_the_current_default(legacy, gas):
    # The fixture's optmat section is regenerated when the default changes.
    field = np.array(legacy["field_v_per_m"])
    model = ionisation_rate_model(gas, legacy["wavelength"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        current = model.instantaneous_rate(field)
    recorded = np.array(legacy["optmat"]["gases"][gas]["default"]["rate"], dtype=float)
    assert recorded == pytest.approx(current, rel=1e-9, abs=0.0)


def test_luna_refuses_strong_fields_and_solprop_does_not(legacy):
    assert legacy["luna"]["gases"]["Xe"]["refused"]
    assert np.all(np.isfinite(_rates(legacy["solprop"], "Xe")))
