# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The compact infrared terms of the molecular fits (OEP 2)."""

import dataclasses
import json
import math
import warnings
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest

from optmat import _validity
from optmat.conditions import ONE_ATM_20C, State
from optmat.constants import AU_ANGULAR_FREQUENCY, C
from optmat.exceptions import DataError, ModelRangeWarning, OptmatError
from optmat.gas.eos import number_density
from optmat.gas.fitting import HARTREE_NM
from optmat.gas.gases import Gas
from optmat.gas.models import compact
from optmat.gas.models.compact import (
    COMPACT_CONFIG_DIRECTORY,
    COMPACT_MODEL_DIRECTORY,
    CompactConfig,
    absorbing_regions,
    coefficients,
    compact_intervals,
    compact_models,
    compact_record,
    derive,
    format_coefficients,
    load_compact_config,
    parent_digest,
    pole_dispersion_au,
    reproduction,
    residuals_sigma,
    statistical_sigma_au,
    vibrational_terms,
)
from optmat.gas.models.fitted import FittedModel, fitted_models
from optmat.gas.models.registry import default_model, model, models

GASES = sorted(p.stem for p in COMPACT_CONFIG_DIRECTORY.glob("*.toml"))


def _values(stored, key):
    """Return a list from a stored compact record, typed."""
    return cast("list[Any]", stored.record[key])


def _pair(gas):
    """Return a compact model, its parent and its configuration."""
    config = load_compact_config(gas)
    return compact_models()[f"compact/{gas}"], fitted_models()[config.parent], config


def test_every_molecule_with_fixed_infrared_terms_has_a_compact_model():
    # O3's infrared is its measured absorption, not fixed poles (OEP 2, section 2).
    molecules = {
        name.removeprefix("fit/")
        for name, fit in fitted_models().items()
        if fit.n_fixed > fit.n_rotational and fit.absorption is None
    }
    assert set(GASES) == molecules == {"CH4", "CO", "CO2", "N2O", "SF6"}
    assert sorted(p.stem for p in COMPACT_MODEL_DIRECTORY.glob("*.json")) == GASES


@pytest.mark.parametrize("gas", GASES)
def test_stored_compact_model_is_fresh(gas):
    stored, parent, config = _pair(gas)
    assert stored.record["parent_sha256"] == parent_digest(config.parent)
    assert stored.record["parent_config_sha256"] == parent.record["config_sha256"]
    assert stored.record["config_sha256"] == config.digest
    # The residuals recorded are those of the stored poles on the check grid today.
    found = residuals_sigma(
        parent,
        stored.intervals,
        _values(stored, "strengths"),
        _values(stored, "energies_hartree"),
        int(cast("int", stored.record["check_points"])),
    )
    for key, value in stored.residuals_sigma.items():
        # Equal to rounding; a kept model's residuals are rounding alone, near 1e-13.
        assert float(np.max(np.abs(found[key]))) == pytest.approx(
            value, rel=1e-9, abs=1e-9
        )


@pytest.mark.parametrize("gas", GASES)
def test_compact_model_meets_the_tolerance(gas):
    stored, _, config = _pair(gas)
    # 0.1 of the parent's statistical sigma in alpha, D2 and D3 on the check grid (OEP
    # 2, section 7), over the compact intervals.
    assert stored.record["passed"]
    assert max(stored.residuals_sigma.values()) < config.tolerance_sigma == 0.1


@pytest.mark.parametrize("gas", GASES)
def test_compact_model_keeps_the_moments_and_the_electronic_terms(gas):
    stored, parent, _ = _pair(gas)
    f, e = vibrational_terms(parent)
    g, h = (np.array(stored.record[k]) for k in ("strengths", "energies_hartree"))
    assert np.all(g > 0)
    # Exact to the inner solve's 1e-10 (OEP 2, section 8).
    assert g.sum() == pytest.approx(f.sum(), rel=1e-10)
    assert np.sum(g / h**2) == pytest.approx(np.sum(f / e**2), rel=1e-10)
    n = len(parent.strengths) - parent.n_fixed + parent.n_rotational
    assert stored.strengths[:n] == parent.strengths[:n]
    assert stored.energies[:n] == parent.energies[:n]
    assert stored.labels[:n] == parent.labels[:n]
    np.testing.assert_array_equal(stored.covariance, parent.covariance)
    assert stored.reduced_chi2 == parent.reduced_chi2
    assert stored.n_lines == parent.n_lines
    assert stored.n_rotational == parent.n_rotational


@pytest.mark.parametrize("gas", GASES)
def test_no_compact_pole_lies_inside_an_interval(gas):
    # An analytic check, independent of the grid: every pole's wavelength is outside
    # every compact interval, and inside its absorbing region's bounds.
    stored, parent, _ = _pair(gas)
    energies = np.array(_values(stored, "energies_hartree"))
    wavelength = HARTREE_NM / energies * 1e-9
    assert np.all(_validity.outside(wavelength, stored.intervals))
    regions = _values(stored, "regions_hartree")
    for energy, index in zip(energies, _values(stored, "region"), strict=True):
        low, high = regions[index]
        assert low <= energy <= high
    # The intervals are the parent's, inside its optical range.
    low, high = _values(parent, "range_nm")
    for a, b in stored.intervals:
        assert low * 1e-9 <= a < b <= high * 1e-9
        assert not np.any(_validity.outside(np.array([a, b]), parent.intervals))


@pytest.mark.parametrize("gas", GASES)
def test_compact_model_reproduces_the_parents_data(gas):
    stored, _, _ = _pair(gas)
    report = reproduction(stored)
    assert report
    for key, entry in report.items():
        # The constrained M-2 (OEP 2, section 7).
        assert entry.get("static_relative", 0.0) < 1e-10, key
        # 0.1 of each datum's own sigma, inside the compact intervals.
        assert entry.get("max_sigma", 0.0) < 0.1, key


@pytest.mark.parametrize("gas", GASES)
def test_the_tolerance_holds_at_ten_bar(gas):
    # The tolerance is a dilute-limit one; at 10 bar the relative error in n - 1 is the
    # relative error in alpha to the order of the susceptibility (OEP 2, section 4).
    stored, parent, _ = _pair(gas)
    (a, b), *_ = stored.intervals
    wavelength = np.geomspace(a, b, 200)
    state = State(10e5, 293.15)
    dn = (
        stored.refractivity(wavelength, state) / parent.refractivity(wavelength, state)
        - 1
    )
    da = stored.polarizability(wavelength) / parent.polarizability(wavelength) - 1
    np.testing.assert_array_less(np.abs(dn - da), 2e-3 * np.abs(da) + 1e-15)


def test_compact_models_join_the_registry_and_are_never_a_default():
    for gas in GASES:
        assert f"compact/{gas}" in models(gas)
        assert model(f"compact/{gas}").status == "comparison"
        assert default_model(gas).name == f"fit/{gas}"
    with pytest.raises(OptmatError, match="no model named"):
        model("compact/N2")


def test_beyond_its_intervals_the_warning_names_the_parent():
    co = compact_models()["compact/CO"]
    with pytest.warns(ModelRangeWarning, match="fit/CO is qualified there"):
        co.polarizability(3e-6)  # in fit/CO's window at 2.56-4.34 um
    with pytest.warns(ModelRangeWarning) as caught:
        co.polarizability(2.4e-6)  # where fit/CO is not qualified either
    assert "qualified there" not in str(caught[0].message)
    assert co.parent == "fit/CO"
    assert co._sources() == fitted_models()["fit/CO"]._sources()


def test_sf6_keeps_its_parents_poles():
    sf6, parent, _ = _pair("SF6")
    # The same poles, in ascending energy.
    assert sorted(zip(sf6.energies, sf6.strengths, strict=True)) == sorted(
        zip(parent.energies, parent.strengths, strict=True)
    )
    assert all(_values(sf6, "kept"))
    assert sf6.labels[-1] == "compact (the parent's) 3"


def test_the_fitting_page_table_is_the_stored_reductions():
    page = (Path(__file__).parents[1] / "docs" / "fitting.md").read_text()
    rows = {
        cells[0].strip("`"): cells[1:]
        for line in page.splitlines()
        if line.startswith("| `compact/")
        for cells in [[c.strip() for c in line.strip("|").split("|")]]
    }
    assert set(rows) == {f"compact/{gas}" for gas in GASES}
    for name, (interval, parent_poles, poles, *worst) in rows.items():
        stored = compact_models()[name]
        parent = fitted_models()[stored.parent]
        low, high = stored.intervals[0]
        assert interval == f"{low * 1e9:.0f}-{high * 1e9:.0f} nm"
        assert int(parent_poles) == len(vibrational_terms(parent)[0])
        assert int(poles) == len(_values(stored, "strengths"))
        for cell, key in zip(worst, ("alpha", "D2", "D3"), strict=True):
            assert float(cell) == pytest.approx(
                stored.residuals_sigma[key], rel=0.06, abs=1e-12
            )


# Synthetic parents cover the reduction's branches cheaply.


def _parent(vibrational, *, intervals=((200e-9, 800e-9), (1000e-9, 3000e-9))):
    """A fit of one effective term, a rotational term and some vibrational poles."""
    f = [1.0, 1e-8, *(p[0] for p in vibrational)]
    e = [0.5, 1e-4, *(p[1] for p in vibrational)]
    labels = [
        "effective 1",
        "rotational",
        *(f"infrared {k}" for k in range(len(vibrational))),
    ]
    return FittedModel(
        gas=Gas.CO,
        name="fit/synthetic",
        citation="",
        wavelength_range=(intervals[0][0], intervals[-1][1]),
        intervals=intervals,
        kind="fitted",
        strengths=tuple(f),
        energies=tuple(e),
        labels=tuple(labels),
        n_fixed=len(f) - 1,
        covariance=np.diag([1e-8, 1e-10]),
        reduced_chi2=2.0,
        record={
            "range_nm": [intervals[0][0] * 1e9, intervals[-1][1] * 1e9],
            "n_rotational": 1,
            "config_sha256": "synthetic",
            "gap_reasons": ["absorption: synthetic"] * (len(intervals) - 1),
        },
    )


# Below 3000 nm (0.0152 hartree), between 800 and 1000 nm, and above 200 nm.
_BELOW = [(2e-6, 0.008), (5e-6, 0.010), (1e-6, 0.012)]
_GAP = [(1e-9, 0.050), (2e-9, 0.051), (1e-9, 0.052)]
_ABOVE = [(1e-10, 0.300)]


def _config(poles, **changes):
    return dataclasses.replace(
        CompactConfig(gas=Gas.CO, description="", parent="fit/CO", poles=poles),
        **changes,
    )


def test_the_regions_lie_below_between_and_above_the_intervals():
    parent = _parent(_BELOW + _GAP + _ABOVE)
    regions = absorbing_regions(parent, parent.intervals)
    assert [(r.low, r.high) for r in regions] == [
        (0.008, 0.012),
        (0.050, 0.052),
        (0.3, 0.3),
    ]
    assert len(regions[1].energies) == 3


def test_a_pole_inside_an_interval_is_refused():
    parent = _parent([(1e-6, 0.008), (1e-9, 0.1)])  # 0.1 hartree is 456 nm
    with pytest.raises(DataError, match="inside a compact interval"):
        absorbing_regions(parent, parent.intervals)


def test_the_intervals_are_the_parents_optical_ones():
    parent = _parent(_BELOW)
    assert compact_intervals(parent) == parent.intervals
    assert compact_intervals(parent, 900.0) == ((200e-9, 800e-9),)
    with pytest.raises(DataError, match="no qualified interval"):
        compact_intervals(parent, 150.0)


def test_a_reduction_searches_some_regions_and_keeps_others():
    parent = _parent(_BELOW + _GAP + _ABOVE)
    reduction = derive(
        parent, _config((2, 2, 1), tolerance_sigma=1e12, random_starts=1)
    )
    assert reduction.passed
    assert len(reduction.energies) <= 5
    # The region above, given one pole for its one, keeps it.
    assert reduction.kept[-1]
    assert max(reduction.energies) == 0.3
    assert max(abs(x) for x in reduction.moment_errors) < 1e-10
    strict = derive(parent, _config((2, 2, 1), tolerance_sigma=0.0))
    assert not strict.passed


def test_a_region_given_every_pole_keeps_them_unchanged():
    parent = _parent(_BELOW + _ABOVE)
    reduction = derive(parent, _config((3, 1)))
    assert reduction.start == "the parent's poles"
    assert reduction.energies == tuple(e for _, e in _BELOW + _ABOVE)
    assert all(reduction.kept)
    assert max(reduction.worst.values()) < 1e-6  # rounding alone


def test_one_searched_pole_sits_at_the_moments_mean():
    parent = _parent(_BELOW)
    f, e = vibrational_terms(parent)
    reduction = derive(parent, _config((1,), random_starts=2))
    # Tang's pole for all the bands together.
    assert reduction.energies == pytest.approx((math.sqrt(f.sum() / np.sum(f / e**2)),))


def test_poles_that_cannot_keep_both_moments_are_refused():
    parent = _parent(_BELOW + _GAP + _ABOVE)
    with pytest.raises(DataError, match="cannot keep both moments"):
        derive(parent, _config((0, 0, 1)))  # one kept pole, above the mean
    with pytest.raises(DataError, match="cannot keep both moments"):
        derive(parent, _config((0, 1, 0)))  # one pole, bounded above the mean
    with pytest.raises(DataError, match="3 regions"):
        derive(parent, _config((2,)))


def test_pole_dispersion_is_the_derivative_of_omega_alpha():
    omega = np.array([0.02, 0.05])
    f, e = np.array([1e-6]), np.array([0.03])
    step = 1e-5
    value = pole_dispersion_au(f, e, omega, 2)[0]

    def omega_alpha(w):
        return w * pole_dispersion_au(f, e, w, 0)[0]

    # The second central difference, to its O(step^2) truncation.
    numeric = (
        omega_alpha(omega + step) - 2 * omega_alpha(omega) + omega_alpha(omega - step)
    ) / step**2
    np.testing.assert_allclose(value, numeric, rtol=1e-5)


def test_a_compact_record_names_the_gap_between_two_intervals():
    parent = _parent(_BELOW + _GAP + _ABOVE)
    config = _config((2, 3, 1), digest="abc")
    record = compact_record(parent, config, derive(parent, config))
    assert record["gap_reasons"] == ["absorption: synthetic"]
    assert record["parent_sha256"] == parent_digest("fit/CO")
    assert record["config_sha256"] == "abc"
    assert json.loads(json.dumps(record)) == record


def test_a_configuration_is_checked(tmp_path):
    path = tmp_path / "X.toml"
    head = 'gas = "CO"\ndescription = ""\nparent = "fit/CO"\n'
    path.write_text(head + "poles = [2]\nlongest_nm = 1500\n")
    config = load_compact_config(path)
    assert config.longest_nm == 1500.0
    assert config.poles == (2,)
    assert load_compact_config("CO").longest_nm is None
    path.write_text('gas = "CO"\ndescription = ""\nparent = "fit/CO"\npoles = [-1]\n')
    with pytest.raises(DataError, match="cannot have -1 poles"):
        load_compact_config(path)
    path.write_text(
        'gas = "CO"\ndescription = ""\nparent = "fit/CO"\npoles = [1]\nK = 2\n'
    )
    with pytest.raises(DataError, match="unknown keys"):
        load_compact_config(path)


@pytest.mark.fits
@pytest.mark.parametrize("gas", GASES)
def test_stored_compact_model_matches_a_fresh_derivation(gas):
    stored, parent, config = _pair(gas)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        fresh = derive(parent, config)
    _assert_rederived(stored.record, fresh, config)


# The derivation has a fixed grid, starts and seed, but the poles of a searched region
# are nearly degenerate: re-derivations of the committed code move CO2's strengths by
# 4.4 % and N2O's by 0.5 % from the stored ones, while their worst residuals on the
# check grid move by at most 2.2e-3 sigma of the parent's band (TC-1090). The check is
# on those residuals (0.01 sigma), the pole structure and the moments.
RESIDUAL_SIGMA = 0.01


def _assert_rederived(record, fresh, config):
    """Fail unless a fresh derivation is the stored compact model's equal."""
    assert len(fresh.strengths) == len(record["strengths"])
    assert list(fresh.region) == list(record["region"])
    assert list(fresh.kept) == list(record["kept"])
    assert fresh.start == record["start"]
    assert fresh.passed and record["passed"]
    assert max(abs(e) for e in fresh.moment_errors) < 1e-12
    for quantity, value in record["residuals_sigma"].items():
        assert fresh.worst[quantity] < config.tolerance_sigma
        assert abs(fresh.worst[quantity] - value) <= RESIDUAL_SIGMA, quantity


def test_a_re_derivation_equal_to_the_record_is_accepted():
    stored, _, config = _pair("CO2")
    record = cast("dict[str, Any]", stored.record)
    same = compact.Reduction(
        strengths=tuple(record["strengths"]),
        energies=tuple(record["energies_hartree"]),
        region=tuple(record["region"]),
        kept=tuple(record["kept"]),
        worst=dict(record["residuals_sigma"]),
        moment_errors=tuple(record["moment_relative_errors"]),
        start=record["start"],
        passed=record["passed"],
    )
    _assert_rederived(record, same, config)
    with pytest.raises(AssertionError):
        worse = {k: v + 0.02 for k, v in record["residuals_sigma"].items()}
        _assert_rederived(record, dataclasses.replace(same, worst=worse), config)
    with pytest.raises(AssertionError):
        _assert_rederived(record, dataclasses.replace(same, start="random 1"), config)


# The printed forms (OEP 2, Phase 2).


def _parse(text):
    """Return the rows of a formatted table: label, kind and two numbers."""
    rows = []
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) == 4 and cells[1] in {
            "line",
            "effective",
            "infrared",
            "rotational",
        }:
            rows.append((cells[0], cells[1], float(cells[2]), float(cells[3])))
    return rows


@pytest.mark.parametrize("gas", GASES)
def test_the_printed_polarizability_form_round_trips(gas):
    stored, parent, _ = _pair(gas)
    table = coefficients(gas)
    rows = _parse(format_coefficients(table))
    assert [r[0] for r in rows] == list(table.labels)
    assert [r[1] for r in rows] == list(table.kinds)
    f = np.array([r[2] for r in rows])
    energies = (
        2 * np.pi * C / (np.array([r[3] for r in rows]) * 1e-6) / AU_ANGULAR_FREQUENCY
    )
    omega = np.concatenate(
        [HARTREE_NM / np.geomspace(a * 1e9, b * 1e9, 1600) for a, b in stored.intervals]
    )
    sigma = statistical_sigma_au(parent, omega)
    for m in (0, 2, 3):
        printed = pole_dispersion_au(f, energies, omega, m).sum(axis=0)
        exact = pole_dispersion_au(
            np.array(stored.strengths), np.array(stored.energies), omega, m
        ).sum(axis=0)
        # The rounding budget: 0.01 of the parent's sigma (OEP 2, section 7).
        assert np.max(np.abs(printed - exact) / sigma[m]) < 0.01
    assert table.chi0 == 0.0
    assert table.model_digest == stored.identity.digest
    assert table.parent_digest == parent.identity.digest


@pytest.mark.parametrize("gas", GASES)
def test_the_printed_sellmeier_form_round_trips(gas):
    state = ONE_ATM_20C
    table = coefficients(
        gas, "sellmeier", pressure=state.pressure, temperature=state.temperature
    )
    assert table.number_density == pytest.approx(number_density(gas, state), rel=1e-12)
    rows = _parse(format_coefficients(table))
    b = np.array([r[2] for r in rows])
    c = np.array([r[3] for r in rows]) * 1e-12
    wavelength = np.concatenate([np.geomspace(a, z, 400) for a, z in table.intervals])
    x2 = wavelength[:, None] ** 2
    printed = np.sum(b * x2 / (x2 - c), axis=1)
    stored = np.sum(np.array(table.B) * x2 / (x2 - np.array(table.C)), axis=1)
    np.testing.assert_allclose(printed, stored, rtol=1e-8)
    # The stated error against the full parent, which includes Lorentz-Lorenz and the
    # virial (about (n - 1) / 6, OEP 2, section 8).
    full = fitted_models()[table.parent].refractivity(wavelength, state, warn=False)
    assert np.max(np.abs(printed / full - 1)) <= table.stated_error["n_minus_1"] * 1.001
    assert 5e-5 < table.stated_error["n_minus_1"] < 3e-4
    assert table.stated_error["beta2"] < 1e-3


@pytest.mark.parametrize("gas", GASES)
def test_the_printed_covariance_gives_the_parents_sigma(gas):
    table = coefficients(gas)
    parent = fitted_models()[table.parent]
    omega = np.geomspace(0.01, 0.15, 7)
    index = {label: k for k, label in enumerate(table.labels)}
    gradient = []
    for name in table.covariance_parameters:
        kind, label = name.split(" ", 1)
        k = index[label]
        f, lam = table.strengths[k], table.wavelengths[k]
        e = 2 * np.pi * C / lam / AU_ANGULAR_FREQUENCY
        if kind == "f":
            gradient.append(1 / (e**2 - omega**2))
        else:  # d alpha / d lambda = (d alpha / d E) (-E / lambda)
            gradient.append(-2 * f * e / (e**2 - omega**2) ** 2 * (-e / lam))
    g = np.array(gradient)
    sigma = np.sqrt(np.einsum("in,ij,jn->n", g, table.covariance, g))
    # Equal to rounding: CO2's effective terms are correlated to 0.9998 (condition
    # number 8e11), so its sigma is a cancellation, rounded differently in the two
    # parameterisations (5e-7).
    np.testing.assert_allclose(sigma, parent.polarizability_sigma_au(omega), rtol=1e-5)
    sellmeier = coefficients(gas, "sellmeier", pressure=1e5, temperature=293.15)
    assert sellmeier.covariance.shape == table.covariance.shape
    assert all(p.split(" ")[0] in {"B", "C"} for p in sellmeier.covariance_parameters)
    assert not sellmeier.covariance.flags.writeable


def test_the_printed_forms_are_refused_where_there_is_none():
    with pytest.raises(OptmatError, match="no compact model"):
        coefficients("N2")
    with pytest.raises(ValueError, match="form is"):
        coefficients("CO", "cauchy")
    with pytest.raises(ValueError, match="needs a pressure"):
        coefficients("CO", "sellmeier", pressure=1e5)


def test_the_gas_models_page_prints_the_coefficients():
    page = (Path(__file__).parents[1] / "docs" / "gas_models.md").read_text()
    block = page.split("<!-- compact-coefficients:start -->\n")[1]
    block = block.split("\n<!-- compact-coefficients:end -->")[0]
    state = ONE_ATM_20C
    expected = []
    for gas in ("CO", "CO2", "CH4", "N2O", "SF6"):
        expected.append(format_coefficients(coefficients(gas)))
        sellmeier = coefficients(
            gas, "sellmeier", pressure=state.pressure, temperature=state.temperature
        )
        expected.append(format_coefficients(sellmeier))
    assert block == "\n\n".join(expected)
