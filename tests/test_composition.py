# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The composition of dry air (TC-548)."""

import math
from pathlib import Path

import pytest

from optmat.exceptions import OptmatError
from optmat.gas.composition import composition


def test_dry_air_is_cipm_2007():
    air = composition("dry air")
    fractions = {c.material: c.mole_fraction for c in air}
    # Picard et al. 2008, Table 1, as printed.
    assert fractions["N2"] == 0.780848
    assert fractions["O2"] == 0.209390
    assert fractions["Ar"] == 0.009332
    assert fractions["CO2"] == 0.00040
    assert len(air) == 12
    # The printed table sums to 1 - 2.9e-6 (its rounding); summing twelve floats adds
    # at most 1e-15.
    assert math.isclose(sum(fractions.values()), 1 - 2.9e-6, abs_tol=1e-12)
    assert [c.material for c in air][:3] == ["N2", "O2", "Ar"]


def test_the_trace_constituents_are_marked_as_reprints():
    air = {c.material: c for c in composition("dry air (CIPM-2007)")}
    assert air["Ar"].primary
    assert air["Ar"].sigma == 3e-6
    assert not air["Ne"].primary
    assert math.isnan(air["N2"].sigma)


def test_the_page_prints_the_fractions():
    page = (Path(__file__).parents[1] / "docs" / "models.md").read_text(
        encoding="utf-8"
    )
    for printed in ("0.780848", "0.209390", "0.009332 (u = 3e-6)", "0.00040"):
        assert printed in page


def test_an_unknown_mixture_is_refused():
    with pytest.raises(OptmatError, match="no curated composition"):
        composition("martian air")
