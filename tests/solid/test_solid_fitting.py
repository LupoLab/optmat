# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import dataclasses
import json

import numpy as np
import pytest

from optmat.datasets import datasets
from optmat.solid import SOLID_DATA_DIRECTORY, fitting
from optmat.solid.fitting import (
    DatasetOptions,
    _observations,  # pyright: ignore[reportPrivateUsage]
    dumps,
    fit,
    load_config,
)
from optmat.solid.models.fitted import FittedModel
from optmat.solid.models.formulas import FormulaModel
from optmat.solid.models.registry import default_model, models
from optmat.solid.udm import from_parameters

CONFIG = """
material = "SiO2"
description = "test: recover two perturbed parameters of Franta's model from its own n"
range_nm = [30.0, 125000.0]

[model]
start = "franta_2016/udm_lithosil_q2"
free = ["Nvx", "Ex"]

[model.scales]
OH = ["Nph36", "Nph37", "Nph38", "Nph39", "Nph40", "Nph41"]

[model.add]
values = { Nvx = 700.0, Ex = 17.5 }

[[fitted]]
key = "franta_2016/lithosil_q2_n"
range_nm = [150.0, 5000.0]
sigma = 1e-5
offset_sigma = 1e-3

[[validators]]
key = "franta_2016/lithosil_q2_k"
range_nm = [5000.0, 125000.0]
sigma = 0.1
"""


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    path = tmp_path_factory.mktemp("fits") / "recover.toml"
    path.write_text(CONFIG)
    return fit(load_config(path))


def test_a_fit_recovers_the_parameters_that_made_its_data(result):
    # The data are the model's own n (printed to 7 digits: rounding 5e-8, far below
    # sigma), so the fit returns the file's values to well within their uncertainty.
    sigma = np.sqrt(np.diag(result.covariance))
    assert result.values["Nvx"] == pytest.approx(735.4107760577147, abs=3 * sigma[0])
    assert result.values["Ex"] == pytest.approx(16.98219466591156, abs=3 * sigma[1])
    assert result.values["Nvx"] == pytest.approx(735.4107760577147, rel=1e-3)
    assert result.offsets["franta_2016/lithosil_q2_n"] == pytest.approx(0, abs=1e-6)
    # The OH scale is untouched by these data: it stays at 1 within 1e-3.
    assert result.values["Nph38"] == pytest.approx(4.529073137525764e-06, rel=1e-3)


def test_the_stored_record_rebuilds_the_model_with_its_uncertainty(result):
    record = json.loads(dumps(result))
    assert record["free"] == ["Nvx", "Ex"]
    assert set(record["validation"]) == {"franta_2016/lithosil_q2_k"}
    attributes = tuple(record["attributes"])
    values = record["values"]
    model = FittedModel(
        material="SiO2",
        name="fit/test",
        citation="",
        wavelength_range=(30e-9, 125e-6),
        kind="fitted",
        udm=from_parameters(attributes, values),
        attributes=attributes,
        values=values,
        free=tuple(record["free"]),
        scales={k: tuple(v) for k, v in record["scales"].items()},
        covariance=np.asarray(record["covariance"]),
    )
    sigma_n, sigma_k = model.index_sigma(np.array([200e-9, 800e-9]))
    assert np.all(sigma_n > 0)
    assert np.all(np.isfinite(sigma_k))
    # Two parameters fitted to 2000 points at 1e-5: the propagated sigma of n in the
    # visible is of the order of that sigma (1.7e-5 here), not of Franta's 1e-3.
    assert sigma_n[1] < 1e-4


def test_the_recommended_silica_model_is_the_refit():
    chosen = default_model("fused silica")
    assert chosen.name == "fit/SiO2"
    assert set(models("SiO2")) >= {
        "fit/SiO2",
        "franta_2016/udm_lithosil_q2",
        "leviton_2008/sellmeier",
        "malitson_1965/sellmeier",
        "arosa_2020/sellmeier_fit3",
    }


def test_the_absolute_refractometry_agrees_and_the_refit_follows_it():
    # At 800 nm and 20 degrees C (TC-399): Leviton & Frey's absolute formula and
    # Malitson's, taken from room air to vacuum, agree to 1e-5 in n and 2e-5 in group
    # index; the refit follows them to 2e-5 in n and 3e-5 in group index, and Franta's
    # published model is 1.6e-4 lower in n and 3e-4 in group index. Regression checks
    # of the finding at the precision the sources state (1e-5).
    all_models = models("SiO2")
    formula = all_models["leviton_2008/sellmeier"]
    assert isinstance(formula, FormulaModel)
    leviton = formula.at_temperature(293.15)
    malitson = all_models["malitson_1965/sellmeier"]
    refit = default_model("SiO2")
    franta = all_models["franta_2016/udm_lithosil_q2"]
    n = {
        name: float(np.real(model.refractive_index(800e-9)))
        for name, model in [
            ("leviton", leviton),
            ("malitson", malitson),
            ("refit", refit),
            ("franta", franta),
        ]
    }
    assert n["malitson"] == pytest.approx(n["leviton"], abs=1e-5)
    assert n["refit"] == pytest.approx(n["leviton"], abs=2e-5)
    assert n["franta"] - n["leviton"] == pytest.approx(-1.6e-4, abs=2e-5)
    assert malitson.group_index(800e-9) == pytest.approx(
        leviton.group_index(800e-9), abs=2e-5
    )
    assert refit.group_index(800e-9) == pytest.approx(
        leviton.group_index(800e-9), abs=3e-5
    )
    assert franta.group_index(800e-9) - leviton.group_index(800e-9) == pytest.approx(
        -3.1e-4, abs=2e-5
    )


def test_the_recommended_silicon_model_is_the_refit():
    assert default_model("Si").name == "fit/Si"
    assert set(models("Si")) >= {
        "fit/Si",
        "franta_2017/udm_c_si",
        "frey_2006/sellmeier_Si",
    }


def test_the_silicon_refit_follows_freys_absolute_index():
    # TC-1088: Frey et al.'s 295 K table taken to 20 C with their own dn/dT. The
    # refit is within 1.4e-4 of it (their stated +-1e-4 plus the model's shape); the
    # published newAD2 model is 0.0158 above at 2 um.
    keyed = {d.key: d for d in datasets()}
    table = keyed["frey_2006/n_Si_295K"]
    slope = keyed["frey_2006/dndt_Si_295K"].scaled_values()
    wavelength = table.vacuum_wavelength()
    frey = table.scaled_values() - slope * (295 - 293.15)
    keep = wavelength > 1.15e-6
    refit = np.real(default_model("Si").refractive_index(wavelength[keep]))
    assert np.max(np.abs(refit - frey[keep])) < 1.4e-4
    newad2 = models("Si")["franta_2017/udm_c_si"]
    at_2um = np.argmin(np.abs(wavelength - 2e-6))
    assert float(
        np.real(newad2.refractive_index(wavelength[at_2um])) - frey[at_2um]
    ) == pytest.approx(0.0158, abs=2e-4)


def test_an_absorption_coefficient_is_fitted_as_k():
    # k = alpha lambda / (4 pi) for a natural-log alpha in cm-1; Schinke et al.'s
    # 1.7e4 cm-1 at 500 nm is k = 0.0676, and the relative uncertainty is kept.
    options = DatasetOptions(key="schinke_2015/alpha_295K")
    observed = _observations(options)
    dataset = next(d for d in datasets() if d.key == options.key)
    alpha = dataset.scaled_values() * 100
    wavelength = dataset.vacuum_wavelength()
    k = alpha * wavelength / (4 * np.pi)
    np.testing.assert_allclose(np.exp(observed.value), k, rtol=1e-12)
    relative = dataset.standard_uncertainty() / dataset.scaled_values()
    np.testing.assert_allclose(observed.sigma, relative, rtol=1e-12)
    assert observed.component == "log_k"


def test_a_decadic_absorption_coefficient_per_metre_is_converted(monkeypatch):
    # The same values read as a base-10 alpha in m-1 give k ln(10) / 100 times.
    dataset = next(d for d in datasets() if d.key == "schinke_2015/alpha_295K")
    decadic = dataclasses.replace(
        dataset, unit="m-1", parameters={**dataset.parameters, "base": "10"}
    )
    natural = _observations(DatasetOptions(key=dataset.key))
    monkeypatch.setattr(fitting, "datasets", lambda: [decadic])
    converted = _observations(DatasetOptions(key=dataset.key))
    np.testing.assert_allclose(
        converted.value - natural.value, np.log(np.log(10) / 100), rtol=1e-12
    )


# A fresh refit against the stored record (TC-1089). The minimum is flat along a few
# parameter combinations, and least_squares stops on its 1e-8 tolerances at a point
# that depends on the path: three settings of the unknowns' scaling move silica's
# parameters by up to 3.3e-3 of their standard uncertainties, n by 1.9e-6 relative (in
# the extreme ultraviolet; 2e-9 in the window), k by 2.9e-5 relative and chi2 by 2e-9
# relative, all at one objective. The bounds are 3-5 times those, far below the models'
# uncertainties; the stored records are not rewritten to meet them.
PARAMETER_SIGMAS = 1e-2
INDEX_RELATIVE = 1e-5
EXTINCTION_RELATIVE = 1e-4
CHI2_RELATIVE = 1e-8
FITTED = ("SiO2", "Si")


def _stored(name):
    path = SOLID_DATA_DIRECTORY / "models" / "fitted" / f"{name}.json"
    return json.loads(path.read_text())


def _unknowns(config, record):
    """Return the free parameters, the group scales and the offsets, in fit order."""
    _, start = fitting._start(config)  # pyright: ignore[reportPrivateUsage]
    scales = [record["values"][m[0]] / start[m[0]] for m in config.scales.values()]
    free = [record["values"][p] for p in config.free]
    return np.array(free + scales + list(record["offsets"].values()))


def _assert_refit(config, stored, fresh):
    """Fail unless a fresh fit record agrees with the stored one within round-off."""
    assert fresh["config_sha256"] == stored["config_sha256"]
    assert fresh["dof"] == stored["dof"]
    assert list(fresh["offsets"]) == list(stored["offsets"])
    assert {k: v[0] for k, v in fresh["dataset_chi2"].items()} == {
        k: v[0] for k, v in stored["dataset_chi2"].items()
    }
    sigma = np.sqrt(np.diag(np.asarray(stored["covariance"])))
    shift = np.abs(_unknowns(config, fresh) - _unknowns(config, stored)) / sigma
    assert shift.max() <= PARAMETER_SIGMAS, shift.max()
    low, high = config.range_nm
    wavelength = np.geomspace(low * 1e-9, high * 1e-9, 2000)
    old = fitting._index(  # pyright: ignore[reportPrivateUsage]
        from_parameters(tuple(stored["attributes"]), stored["values"]), wavelength
    )
    new = fitting._index(  # pyright: ignore[reportPrivateUsage]
        from_parameters(tuple(fresh["attributes"]), fresh["values"]), wavelength
    )
    np.testing.assert_allclose(np.real(new), np.real(old), rtol=INDEX_RELATIVE, atol=0)
    absorbing = np.imag(old) > 1e-12  # below, k is the model's tail, not a fitted value
    np.testing.assert_allclose(
        np.imag(new)[absorbing],
        np.imag(old)[absorbing],
        rtol=EXTINCTION_RELATIVE,
        atol=0,
    )
    assert abs(fresh["chi2"] - stored["chi2"]) <= CHI2_RELATIVE * max(
        1.0, stored["chi2"]
    )


@pytest.mark.parametrize("name", FITTED)
def test_the_stored_fit_is_current(name):
    assert _stored(name)["config_sha256"] == load_config(name).digest


@pytest.mark.fits
@pytest.mark.parametrize("name", FITTED)
def test_the_stored_fit_matches_a_fresh_fit(name):
    config = load_config(name)
    _assert_refit(config, _stored(name), json.loads(dumps(fit(config))))


def _moved(record, key, factor):
    moved = json.loads(json.dumps(record))
    moved["values"][key] *= factor
    return moved


def test_round_off_in_a_refit_is_accepted():
    config, stored = load_config("SiO2"), _stored("SiO2")
    # Nvc's sigma is 0.28: 1e-4 is 4e-4 sigma, the size of the path dependence.
    _assert_refit(config, stored, _moved(stored, "Nvc", 1 + 1e-4 / 358.65))


@pytest.mark.parametrize(
    "change",
    [
        lambda r: _moved(r, "Nvc", 1 + 0.05 * 0.28 / 358.65),  # 0.05 sigma
        lambda r: {**r, "chi2": r["chi2"] * (1 + 1e-7)},
        lambda r: {**r, "dof": r["dof"] + 1},
        lambda r: {**r, "dataset_chi2": dict(list(r["dataset_chi2"].items())[1:])},
        lambda r: _moved(r, "Nph12", 1 + 1e-3),  # the 9 um band's group scale
    ],
    ids=["parameter", "objective", "dof", "selection", "scale"],
)
def test_a_real_change_in_a_refit_is_rejected(change):
    config, stored = load_config("SiO2"), _stored("SiO2")
    with pytest.raises(AssertionError):
        _assert_refit(config, stored, change(stored))
