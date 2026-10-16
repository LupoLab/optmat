# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Several qualified intervals per fit, and the representation bound (OEP 3, Phase 2).

A configuration with ``band_nodes`` takes its infrared terms from the derived rules and
keeps the rotational system as a term class of its own, which electric observations
leave out; ``windows = true`` makes the derived windows the qualified intervals, with
their gaps' reasons in the range warning. The stored record gains these fields only
when a fit uses them, so every present fit is unchanged.
"""

import json
import math
import re
from typing import cast

import numpy as np
import pytest

from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas.fitting import (
    FIT_DIRECTORY,
    DatasetOptions,
    Observation,
    _derived_record,  # pyright: ignore[reportPrivateUsage]
    _model_values,  # pyright: ignore[reportPrivateUsage]
    load_config,
)
from optmat.gas.infrared import (
    band_systems,
    bound_cells,
    derived_digest,
    qualified_windows,
)
from optmat.gas.models.fitted import fitted_models

HARTREE_NM = 45.56335252767


def _with_keys(tmp_path, gas, lines):
    """Load a gas's configuration with its OEP 3 keys replaced by ``lines``."""
    text = (FIT_DIRECTORY / f"{gas}.toml").read_text(encoding="utf-8")
    text = re.sub(
        r"^(infrared_max_wavenumber|band_nodes|windows) = .*\n",
        "",
        text,
        flags=re.MULTILINE,
    )
    path = tmp_path / f"{gas}.toml"
    path.write_text(text.replace("[model]\n", "[model]\n" + lines, 1), encoding="utf-8")
    return load_config(path)


def test_band_nodes_takes_the_rules_and_the_rotational_system(tmp_path):
    config = _with_keys(tmp_path, "CO", "band_nodes = 8\n")
    systems = band_systems("CO")
    rules = [s for s in systems if s.treatment != "rotational"]
    assert len(config.infrared) == sum(len(s.strengths) for s in rules)
    rotational = next(s for s in systems if s.treatment == "rotational")
    assert [t[1] for t in config.rotational] == list(rotational.strengths)
    assert config.fixed_terms[: len(config.rotational)] == config.rotational
    assert config.derived_sha256 == derived_digest("CO")
    assert config.intervals_nm == (config.range_nm,)
    record = _derived_record(config)
    assert record["n_rotational"] == len(config.rotational)
    assert record["band_nodes"] == 8


def test_band_nodes_is_checked(tmp_path):
    with pytest.raises(DataError, match="rules of 8 nodes"):
        _with_keys(tmp_path, "CO", "band_nodes = 4\n")
    with pytest.raises(DataError, match="replaces infrared_max_wavenumber"):
        _with_keys(tmp_path, "CO", "band_nodes = 8\ninfrared_max_wavenumber = 5000\n")


def test_windows_become_the_intervals(tmp_path):
    config = _with_keys(tmp_path, "CO2", "windows = true\n")
    intervals, reasons = qualified_windows("CO2")
    assert config.intervals_nm == intervals
    assert config.gap_reasons == reasons
    assert len(reasons) == len(intervals) - 1
    assert _derived_record(config)["intervals_nm"] == [list(i) for i in intervals]


def test_a_fit_without_the_keys_records_nothing_new(tmp_path):
    config = _with_keys(tmp_path, "CO2", "infrared_max_wavenumber = 7000.0\n")
    assert _derived_record(config) == {}
    assert config.intervals_nm == (config.range_nm,)
    assert not config.rotational


def test_the_rotational_terms_leave_electric_observations(tmp_path):
    config = _with_keys(tmp_path, "CO", "band_nodes = 8\n")
    omega = np.array([0.0, 0.01])
    strengths, energies = np.array([1.0]), np.array([0.5])

    def value(**options):
        obs = Observation(
            key="x",
            omega=omega,
            alpha=np.zeros(2),
            sigma=np.ones(2),
            systematic=None,
            options=DatasetOptions(key="x", **options),
        )
        return _model_values(strengths, energies, obs, config)

    rotational = sum(t[1] / (t[2] ** 2 - omega**2) for t in config.rotational)
    refractive, electric = value(), value(electric=True)
    np.testing.assert_allclose(
        refractive - electric, rotational + config.chi0, rtol=1e-12, atol=1e-18
    )
    assert value(infrared=False)[0] < refractive[0]


@pytest.fixture
def windowed():
    """CO2's stored fit, refitted with the derived windows (OEP 3, Phase 3)."""
    co2 = fitted_models()["fit/CO2"]
    intervals_nm, reasons = qualified_windows("CO2")
    assert [list(i) for i in cast("list[list[float]]", co2.record["intervals_nm"])] == [
        list(i) for i in intervals_nm
    ]
    assert co2.gap_reasons == reasons
    return co2


def test_the_warning_names_the_gaps_reason(windowed):
    (_, low), (high, _) = windowed.intervals[0], windowed.intervals[1]
    reason = windowed.gap_reasons[0]
    with pytest.warns(ModelRangeWarning, match=re.escape(f"is a gap ({reason})")):
        windowed.polarizability(math.sqrt(low * high))
    single = fitted_models()["fit/N2"]
    with pytest.warns(ModelRangeWarning) as caught:
        single.polarizability(single.wavelength_range[0] / 2)
    assert "is a gap" not in str(caught[0].message)
    assert single.gap_reasons == ()
    assert single.n_rotational == 0
    assert fitted_models()["fit/CO"].n_rotational > 0


def test_the_representation_bound(windowed):
    cells = bound_cells("CO2")
    k = len(cells["low_cm"]) // 2
    wavenumber = 0.5 * (cells["low_cm"][k] + cells["high_cm"][k])
    omega = np.array([wavenumber / 219474.6313632])
    for order in (0, 2, 3):
        bound = windowed.representation_bound_au(omega, order)
        assert bound[0] == cells[f"b{order}"][k]
    (_, low), (high, _) = windowed.intervals[0], windowed.intervals[1]
    gap = np.array([HARTREE_NM * 1e-9 / math.sqrt(low * high)])
    with pytest.warns(ModelRangeWarning):
        assert np.isnan(windowed.representation_bound_au(gap, 2)[0])
    with pytest.raises(ValueError, match="order 0, 2 or 3"):
        windowed.representation_bound_au(omega, 1)
    with pytest.raises(DataError, match="no representation bound"):
        fitted_models()["fit/N2"].representation_bound_au(omega, 0)


def test_the_bound_enters_the_identity(windowed):
    record = windowed.identity_record()
    assert windowed.derived_sha256 == derived_digest("CO2")
    assert windowed.derived_sha256 in json.dumps(record, default=str)


def test_beyond_every_interval_is_not_a_gap(windowed):
    with pytest.warns(ModelRangeWarning) as caught:
        windowed.polarizability(windowed.wavelength_range[1] * 2)
    assert "is a gap" not in str(caught[0].message)


def test_a_gas_without_derived_tables():
    with pytest.raises(DataError, match="no derived band_rules table for O3"):
        derived_digest("O3")
