# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The solid tables' rule outside their knots (OEP 1, section 6.10)."""

import warnings

import numpy as np
import pytest

from optmat.exceptions import DataError, ModelRangeWarning, StateError
from optmat.materials import as_material
from optmat.solid.models.registry import model
from optmat.solid.models.tabulated import TabulatedModel


def _silver() -> TabulatedModel:
    silver = model("table/Ag")
    assert isinstance(silver, TabulatedModel)
    return silver


def test_the_tables_hold_by_default():
    silver = _silver()
    assert silver.out_of_range == "hold"
    high = silver.knots[-1]
    with pytest.warns(ModelRangeWarning):
        held = silver.refractive_index(2 * high)
    np.testing.assert_allclose(held, silver.refractive_index(high), rtol=1e-12)


def test_a_caller_can_demand_refusal():
    silver = _silver()
    strict = silver.with_out_of_range("raise")
    assert strict.status == silver.status
    assert strict.identity.digest != silver.identity.digest
    inside = silver.knots[len(silver.knots) // 2]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert strict.refractive_index(inside) == silver.refractive_index(inside)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        with pytest.raises(StateError, match="refuses wavelengths outside"):
            strict.refractive_index(silver.knots[-1] * 2)
        with pytest.raises(StateError):
            strict.dpermittivity_domega(silver.knots[0] / 2)
        # A NaN is not refused: it is not a wavelength outside the table.
        assert np.isnan(strict.refractive_index(np.nan))
    with pytest.raises(DataError, match="'hold' or 'raise'"):
        silver.with_out_of_range("clip")


def test_bk7_and_n_bk7_are_distinct():
    assert as_material("BK7") != as_material("N-BK7")
