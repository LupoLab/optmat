# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The derived band rules, windows and bound of the molecular gases (OEP 3, section 8).

The tables are written by ``tools/hitran/windows.py`` from the HITRAN line lists, which
the tests cannot reach. The tests check what the tables must satisfy by construction,
and that they agree with the band totals ``band_strengths.py`` derived independently
from the same lines.
"""

import csv
import itertools
import math

import numpy as np
import pytest
from scipy import constants

from optmat.datasets import DATA_DIRECTORY
from optmat.exceptions import DataError
from optmat.gas.infrared import (
    band_systems,
    bound_cells,
    derived_provenance,
    infrared_windows,
)

GASES = ("CO", "CO2", "CH4", "N2O")
HARTREE_WAVENUMBER = (
    constants.physical_constants["hartree-inverse meter relationship"][0] / 100
)
PI_RE_CM = math.pi * constants.physical_constants["classical electron radius"][0] * 100
REASONS = {"absorption", "representation", "representation: kept band", "width"}


def _band_totals(gas):
    path = DATA_DIRECTORY / "external" / "hitran" / f"{gas}_bands.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    m0 = sum(float(r["strength_cm_per_molecule"]) for r in rows) / PI_RE_CM
    m_minus2 = (
        sum(float(r["strength_over_nu2_cm3_per_molecule"]) for r in rows)
        / PI_RE_CM
        * HARTREE_WAVENUMBER**2
    )
    return m0, m_minus2


@pytest.mark.parametrize("gas", GASES)
def test_every_rule_is_positive_with_its_nodes_in_its_span(gas):
    for system in band_systems(gas):
        assert all(s > 0 for s in system.strengths), system
        if system.treatment == "folded":
            continue
        low, high = (x / HARTREE_WAVENUMBER for x in system.span_cm)
        # The nodes of a Gaussian rule lie inside its measure's support.
        assert all(
            low * (1 - 1e-12) <= e <= high * (1 + 1e-12) for e in system.energies_au
        )


@pytest.mark.parametrize("gas", GASES)
def test_each_rule_keeps_its_lines_two_moments(gas):
    for system in band_systems(gas):
        f, e = np.array(system.strengths), np.array(system.energies_au)
        # M_0 and M_-2 of the poles are those of the lines, for a rule and a fold; the
        # table's decimal digits allow 1e-12.
        assert f.sum() == pytest.approx(system.m0, rel=1e-12, abs=0.0)
        assert np.sum(f / e**2) == pytest.approx(system.m_minus2_au, rel=1e-12, abs=0.0)


@pytest.mark.parametrize("gas", GASES)
def test_the_rules_sum_to_the_band_tables_totals(gas):
    # The band table (band_strengths.py) and the rules (windows.py) are derived
    # separately from the same lines; their totals agree to the tables' six digits.
    m0, m_minus2 = _band_totals(gas)
    systems = band_systems(gas)
    assert sum(s.m0 for s in systems) == pytest.approx(m0, rel=1e-5, abs=0.0)
    assert sum(s.m_minus2_au for s in systems) == pytest.approx(
        m_minus2, rel=1e-5, abs=0.0
    )


@pytest.mark.parametrize("gas", GASES)
def test_the_windows_and_gaps_tile_the_range_and_every_gap_has_a_reason(gas):
    records = infrared_windows(gas)
    for a, b in itertools.pairwise(records):
        assert a.high_cm == b.low_cm
        assert a.kind != b.kind
    for record in records:
        assert record.low_cm < record.high_cm
        if record.kind == "gap":
            assert record.reason in REASONS
        else:
            assert record.kind == "window"
            assert record.high_cm - record.low_cm >= 50
            assert record.max_absorption_m < 0.1
            # The tolerance qualifies a window beyond the optical range; inside it the
            # ratio is reported only (OEP 3 Phase 1, the user's decision).
            for key in ("max_b0_over_tolerance", "max_b2_over_tolerance"):
                if key in record.details:
                    assert float(record.details[key]) <= 1


@pytest.mark.parametrize("gas", GASES)
def test_no_pole_lies_in_a_window(gas):
    # An analytic check, independent of any grid.
    poles = np.array([e for s in band_systems(gas) for e in s.energies_au])
    poles_cm = poles * HARTREE_WAVENUMBER
    for record in infrared_windows(gas):
        if record.kind == "window":
            inside = (poles_cm > record.low_cm) & (poles_cm < record.high_cm)
            assert not inside.any(), (record, poles_cm[inside])


@pytest.mark.parametrize("gas", GASES)
def test_the_bound_covers_every_window_with_contiguous_cells(gas):
    cells = bound_cells(gas)
    assert np.all(cells["low_cm"] < cells["high_cm"])
    assert np.all(cells["b0"] >= 0) and np.all(cells["b2"] >= 0)
    windows = [r for r in infrared_windows(gas) if r.kind == "window"]
    for record in windows:
        inside = (cells["low_cm"] >= record.low_cm - 1e-3) & (
            cells["high_cm"] <= record.high_cm + 1e-3
        )
        low, high = cells["low_cm"][inside], cells["high_cm"][inside]
        assert low[0] == pytest.approx(record.low_cm, abs=1e-3)
        assert high[-1] == pytest.approx(record.high_cm, abs=1e-3)
        np.testing.assert_allclose(low[1:], high[:-1], atol=1e-9)


@pytest.mark.parametrize("gas", GASES)
def test_the_tables_of_a_gas_share_one_derivation(gas):
    lines = derived_provenance(gas)
    assert any("sha256" in line for line in lines)
    assert any(line.startswith("# settings") for line in lines)
    for table in ("windows", "bound", "absorption"):
        path = DATA_DIRECTORY / "external" / "hitran" / f"{gas}_{table}.csv"
        head = [
            x
            for x in path.read_text(encoding="utf-8").splitlines()
            if x.startswith("#")
        ]
        assert head[: len(lines)] == list(lines)


def test_a_gas_without_derived_tables_is_refused():
    with pytest.raises(DataError, match="no derived"):
        band_systems("N2")
