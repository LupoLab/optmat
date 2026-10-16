# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Synthetic solid tables and fit observations with explicit physical units."""

import dataclasses

import numpy as np
import pytest

import optmat.solid.models.tabulated as tables
from optmat.datasets import Uncertainty, load_source
from optmat.exceptions import DataError, OptmatError
from optmat.solid import fitting
from optmat.solid.models import registry


@pytest.fixture
def pair(tmp_path, monkeypatch):
    path = tmp_path / "synthetic.toml"
    path.write_text(
        'citation="synthetic"\nkind="measurement"\n'
        '[[datasets]]\nid="n"\nmaterial="SiO2"\nquantity="refractive_index"\nunit="1"\n'
        'abscissa_unit="nm"\npoints=[["a","400","1.5"],["b","500","1.6"],'
        '["c","600","1.7"]]\n'
        '[[datasets]]\nid="k"\nmaterial="SiO2"\nquantity="extinction_coefficient"\nunit="1"\n'
        'abscissa_unit="nm"\npoints=[["c","600","0.03"],["a","400","0.01"],'
        '["b","500","0.02"]]\n',
        encoding="utf-8",
    )
    source = load_source(path)
    _serve(monkeypatch, *source.datasets)
    monkeypatch.setattr(tables._catalogue, "source", {source.id: source}.__getitem__)
    return source.datasets


def _serve(monkeypatch, *served):
    # The tables read their datasets by key from the catalogue; serve these instead.
    by_key = {d.key: d for d in served}
    monkeypatch.setattr(tables._catalogue, "by_key", lambda: by_key)


def _table(tmp_path, extra=""):
    path = tmp_path / "table.toml"
    path.write_text(
        'material="SiO2"\n[[segments]]\nreal="synthetic/n"\n'
        'imaginary="synthetic/k"\nrange_nm=[300,700]\nreason="synthetic"\n' + extra,
        encoding="utf-8",
    )
    return path


def test_unsorted_complex_components_align_by_wavelength(pair, tmp_path):
    model = tables.load_table(_table(tmp_path))
    # Each imaginary value must align with its own wavelength, not its input row.
    np.testing.assert_allclose(
        model.index, [1.5 + 0.01j, 1.6 + 0.02j, 1.7 + 0.03j], rtol=2e-15, atol=0
    )
    np.testing.assert_allclose(
        model.refractive_index(model.knots), model.index, rtol=2e-15, atol=0
    )


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"abscissae": ("401", "501", "601")}, "share no wavelength"),
        ({"material": "Al"}, "not both SiO2"),
        ({"values": ("-1", "1.6", "1.7")}, "n or k <= 0"),
        ({"values": ("nan", "nan", "nan")}, "no point"),
        ({"abscissae": ("400", "400", "600")}, "repeated wavelengths"),
    ],
)
def test_unusable_tabulated_segments_fail_explicitly(
    pair, tmp_path, monkeypatch, changes, message
):
    first, second = pair
    first = dataclasses.replace(first, **changes)
    _serve(monkeypatch, first, second)
    with pytest.raises(DataError, match=message):
        tables.load_table(_table(tmp_path))


@pytest.mark.parametrize(
    ("body", "message"),
    [
        ("segments=[]", r"needs \[\[segments\]\]"),
        ("segments=[1]", "every segment needs"),
        ('[[segments]]\nreason="  "', "every segment needs"),
        (
            '[[segments]]\nhenke="al"\nrange_nm=[1,2]\nreason="test"',
            "needs an atom_density",
        ),
    ],
)
def test_table_configuration_requires_segments_and_provenance(tmp_path, body, message):
    path = tmp_path / "table.toml"
    path.write_text('material="SiO2"\n' + body, encoding="utf-8")
    with pytest.raises(DataError, match=message):
        tables.load_table(path)


def test_air_and_thermo_optic_corrections_have_the_right_direction(pair, monkeypatch):
    first, _ = pair
    first = dataclasses.replace(
        first,
        parameters={"index_medium": "air", "temperature_k": 300.0},
        uncertainty=Uncertainty(absolute=0.002),
    )
    coefficient = dataclasses.replace(
        first,
        id="dn_dt",
        quantity="thermo_optic_coefficient",
        unit="1/K",
        parameters={},
        values=("0.0001", "0.0002", "0.0003"),
    )
    monkeypatch.setattr(fitting, "datasets", lambda: (first, coefficient))
    options = fitting.DatasetOptions(
        first.key,
        thermo_optic=coefficient.key,
        temperature_k=290,
        extra=0.003,
        offset_sigma=0.01,
        offset_group="instrument",
    )
    observed = fitting._observations(options)
    # Evaluate the documented Edlen 1966 standard-air and density equations directly,
    # without the shared air-conversion helper; 300 K is inside its stated range.
    s2 = (1000 / np.array([400, 500, 600])) ** 2
    standard = (8342.13 + 2406030 / (130 - s2) + 15997 / (38.9 - s2)) * 1e-8
    t = 300 - 273.15
    density = 760 / 720.775 * (1 + 760 * (0.817 - 0.0133 * t) * 1e-6)
    air = 1 + standard * density / (1 + 0.0036610 * t)
    # Sample is 10 K hotter: subtract dn/dT*10 after converting to absolute n.
    expected = np.array([1.5, 1.6, 1.7]) * air - np.array([0.001, 0.002, 0.003])
    np.testing.assert_allclose(observed.value, expected, rtol=2e-15, atol=0)
    np.testing.assert_allclose(
        observed.sigma, np.sqrt(0.002**2 + 0.003**2), rtol=2e-15, atol=0
    )
    assert observed.offset == "instrument"
    assert observed.offset_sigma == 0.01
    with pytest.raises(DataError, match="thermo_optic needs"):
        fitting._observations(dataclasses.replace(options, temperature_k=None))
    with pytest.raises(DataError, match="thermo_optic needs"):
        fitting._observations(dataclasses.replace(options, thermo_optic=first.key))
    monkeypatch.setattr(
        fitting,
        "datasets",
        lambda: (dataclasses.replace(first, parameters={"index_medium": "air"}),),
    )
    with pytest.raises(DataError, match="needs air_temperature_k"):
        fitting._observations(fitting.DatasetOptions(first.key))


def test_log_extinction_filters_nonpositive_values_and_uses_relative_errors(
    pair, monkeypatch
):
    _, second = pair
    second = dataclasses.replace(
        second, values=("-1", "0.01", "0.02"), uncertainty=Uncertainty(relative=0.1)
    )
    monkeypatch.setattr(fitting, "datasets", lambda: (second,))
    obs = fitting._observations(fitting.DatasetOptions(second.key, offset_sigma=1.0))
    np.testing.assert_allclose(obs.value, np.log([0.01, 0.02]), rtol=2e-15, atol=0)
    np.testing.assert_allclose(obs.sigma, [0.1, 0.1], rtol=2e-15, atol=0)
    assert obs.component == "log_k"
    assert obs.offset_sigma is None
    assert obs.offset == ""


def test_unfittable_or_unweighted_points_require_an_explicit_decision(
    pair, monkeypatch
):
    first, _ = pair
    monkeypatch.setattr(fitting, "datasets", lambda: (first,))
    with pytest.raises(DataError, match="points without an uncertainty"):
        fitting._observations(fitting.DatasetOptions(first.key))
    obs = fitting._observations(fitting.DatasetOptions(first.key, sigma=0.01))
    np.testing.assert_array_equal(obs.sigma, [0.01, 0.01, 0.01])
    first = dataclasses.replace(first, quantity="permittivity_real")
    with pytest.raises(DataError, match="cannot be fitted"):
        fitting._observations(fitting.DatasetOptions(first.key, sigma=0.01))


@pytest.mark.parametrize(
    ("transform", "message"),
    [
        (lambda s: "typo=1\n" + s, "unknown top-level"),
        (lambda s: s.replace("[model]", "[model]\ntypo=1"), r"unknown \[model\]"),
        (lambda s: s + '[[fitted]]\nkey="missing"\n', "unknown dataset"),
        (lambda s: s + '[[fitted]]\nkey="missing"\ntypo=1\n', "unknown dataset keys"),
    ],
)
def test_solid_fit_configuration_rejects_unknown_keys(tmp_path, transform, message):
    text = 'material="SiO2"\nrange_nm=[200,1000]\n[model]\nstart="test"\nfree=[]\n'
    path = tmp_path / "config.toml"
    path.write_text(transform(text), encoding="utf-8")
    with pytest.raises(DataError, match=message):
        fitting.load_config(path)


def test_solid_selection_requires_an_available_model_and_appropriate_axis(monkeypatch):
    with pytest.raises(OptmatError, match="no solid model named"):
        registry.model("missing")
    with pytest.raises(OptmatError, match="axis"):
        registry.default_model("sapphire")
    with pytest.raises(OptmatError, match="axis"):
        registry.default_model("sapphire", axis="x")
    monkeypatch.setattr(registry, "_defaults", lambda: {})
    with pytest.raises(OptmatError, match="recommended"):
        registry.default_model("SiO2")
    with pytest.raises(OptmatError, match="recommended"):
        registry.default_reason("SiO2")


def test_index_deviations_add_the_reference_and_unknown_starts_are_rejected(
    pair, monkeypatch
):
    from optmat.solid.models.formulas import FormulaModel  # noqa: PLC0415

    first, _ = pair
    deviation = dataclasses.replace(
        first,
        quantity="refractive_index_deviation",
        parameters={"reference": "synthetic/reference"},
        values=("0.001", "0.002", "0.003"),
    )
    reference = FormulaModel(
        material="SiO2",
        name="synthetic/reference",
        citation="synthetic",
        wavelength_range=(200e-9, 1e-6),
        kind="published",
        form="n2_sellmeier",
        terms={"B": (1.25,), "C": (0.0,)},
    )
    monkeypatch.setattr(
        fitting, "published_models", lambda: {reference.name: reference}
    )
    monkeypatch.setattr(fitting, "datasets", lambda: (deviation,))
    obs = fitting._observations(
        fitting.DatasetOptions(first.key, reference=reference.name, sigma=0.01)
    )
    # The reference has exactly n=sqrt(1+1.25)=1.5 at every wavelength.
    np.testing.assert_allclose(obs.value, [1.501, 1.502, 1.503], rtol=2e-15, atol=0)
    config = fitting.FitConfig(
        "synthetic", "SiO2", "synthetic", (200, 1000), "missing/model", ()
    )
    with pytest.raises(DataError, match="not a curated UDM formula"):
        fitting.fit(config)


def test_henke_segment_checks_density_and_usable_range(pair, tmp_path, monkeypatch):
    first, _ = pair
    segment = tables.Segment(
        "", "", (0.5, 3), "synthetic", henke="synthetic", atom_density=first.key
    )
    with pytest.raises(DataError, match="not one atom density"):
        tables._henke_points(segment, "SiO2", "synthetic")
    density = dataclasses.replace(
        first,
        quantity="atom_density",
        unit="m-3",
        abscissa="none",
        abscissa_unit="",
        rows=("0",),
        abscissae=("",),
        values=("1e28",),
        sigmas=("",),
        parameters={"counts": "atoms"},
    )
    _serve(monkeypatch, density)
    monkeypatch.setattr(tables, "HENKE_DIRECTORY", tmp_path)
    path = tmp_path / "synthetic.nff"
    path.write_text("energy f1 f2\n1000 1 0.1\n2000 1 0.1\n", encoding="utf-8")
    wavelength, index = tables._henke_points(segment, "SiO2", "synthetic")
    assert np.all(np.diff(wavelength) > 0)
    assert np.all(np.real(index) < 1)
    assert np.all(np.imag(index) > 0)
    with pytest.raises(DataError, match="no point of henke"):
        tables._henke_points(
            dataclasses.replace(segment, range_nm=(10, 20)), "SiO2", "synthetic"
        )
    path.write_text("energy f1 f2\n1000 1 -1\n2000 1 -1\n", encoding="utf-8")
    with pytest.raises(DataError, match="n or k <= 0"):
        tables._henke_points(segment, "SiO2", "synthetic")


def test_fitting_absorption_and_index_with_parameter_priors_uses_both_components():
    # Franta's independently supplied, seven-digit computed n and k tables, already
    # checked against the primary supplemental program by the existing UDM tests.
    config = fitting.load_config("SiO2")
    config = dataclasses.replace(
        config,
        free=("Nvx",),
        scales={},
        added_attributes=(),
        added_values={"Nvx": 700.0},
        priors={"Nvx": (735.4107760577147, 0.1)},
        fitted=(
            fitting.DatasetOptions(
                "franta_2016/lithosil_q2_n", range_nm=(600, 1000), sigma=1e-5
            ),
            fitting.DatasetOptions(
                "franta_2016/lithosil_q2_k", range_nm=(6000, 8000), sigma=0.01
            ),
        ),
        validators=(
            fitting.DatasetOptions(
                "franta_2016/lithosil_q2_n", range_nm=(1000, 1500), sigma=1e-5
            ),
        ),
    )
    result = fitting.fit(config)
    # The source's rounding error is much smaller than the assigned uncertainty.
    assert result.values["Nvx"] == pytest.approx(735.4107760577147, rel=1e-6, abs=0.0)
    assert result.chi2 < 1
    assert set(result.dataset_chi2) == {o.key for o in config.fitted}
    report = fitting.validation(result)
    assert report[config.validators[0].key]["n"] > 0
    with pytest.raises(DataError, match="repeat the start model"):
        fitting.fit(dataclasses.replace(config, added_attributes=("he",)))


def test_solid_fit_inflates_covariance_when_scatter_exceeds_stated_errors(
    tmp_path, monkeypatch
):
    path = tmp_path / "synthetic.toml"
    path.write_text(
        'citation="synthetic"\nkind="theory"\n'
        '[[formulas]]\nid="start"\nterms={attributes="ex=0",Eg="8",Eh="20",Nvc="100"}\n'
        '[[datasets]]\nid="repeat"\nmaterial="SiO2"\nquantity="refractive_index"\n'
        'unit="1"\nabscissa_unit="nm"\npoints=[["a","500","1.499"],["b","500","1.501"]]\n'
        "uncertainty={absolute=0.0001}\n",
        encoding="utf-8",
    )
    source = load_source(path)
    dataset = source.datasets[0]
    monkeypatch.setattr(fitting, "sources", lambda: {source.id: source})
    monkeypatch.setattr(fitting, "datasets", lambda: (dataset,))
    config = fitting.FitConfig(
        "synthetic",
        "SiO2",
        "synthetic repeated measurements",
        (400, 600),
        "synthetic/start",
        ("Nvc",),
        fitted=(fitting.DatasetOptions(dataset.key),),
    )
    first = fitting.fit(config)
    dataset = dataclasses.replace(dataset, values=("1.498", "1.502"))
    wider = fitting.fit(config)
    # Equal weights at one wavelength put the fitted index at their arithmetic mean.
    # Double symmetric scatter: same optimum, four times the Birge variance.
    assert wider.values["Nvc"] == pytest.approx(first.values["Nvc"], rel=1e-7, abs=0.0)
    assert first.chi2 == pytest.approx(200, rel=1e-7, abs=0.0)
    assert wider.chi2 == pytest.approx(800, rel=1e-7, abs=0.0)
    assert wider.covariance[0, 0] / first.covariance[0, 0] == pytest.approx(
        4, rel=1e-6, abs=0.0
    )
