# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Ionisation record boundaries and model assembly, using controlled inputs."""

import dataclasses
import math

import numpy as np
import pytest

import optmat.gas.ionisation_rates as rates
import optmat.gas.molecular_ionisation_rates as molecular
import optmat.gas.molecular_ionisation_validation as validation
from optmat.constants import AU_FIELD, E_CHARGE
from optmat.datasets import load_source
from optmat.exceptions import DataError
from optmat.gas import ionisation
from optmat.gas.gases import Gas


@pytest.fixture
def sample(tmp_path):
    path = tmp_path / "synthetic.toml"
    path.write_text(
        'citation="synthetic"\nkind="theory"\n[[datasets]]\nid="levels"\n'
        'material="Ar+"\nquantity="level_energy"\nunit="eV"\nabscissa="label"\n'
        'points=[["ground","0","0"],["upper","1","2"]]\n'
        'parameters={configuration="synthetic",term="synthetic",j="various"}\n',
        encoding="utf-8",
    )
    return load_source(path).datasets[0]


@pytest.mark.parametrize("reader", [ionisation._row, rates._value])
def test_references_require_an_existing_unambiguous_row(sample, reader):
    by_key = {sample.key: sample}
    with pytest.raises(DataError, match="unknown dataset"):
        reader(by_key, "missing")
    with pytest.raises(DataError, match="name one as"):
        reader(by_key, sample.key)
    with pytest.raises(DataError, match="has no row"):
        reader(by_key, sample.key + "#missing")
    assert reader(by_key, sample.key + "#upper") == (sample, 1)


def test_ionic_channel_weights_and_energy_units_follow_the_reference(sample):
    params = ionisation.ionisation_parameters("Ar")
    channels = rates._channels(
        params, {"fine_structure_level": sample.key + "#upper"}, {sample.key: sample}
    )
    assert [c.weight for c in channels] == [2 / 3, 1 / 3]
    assert channels[0].ionization_energy == params.ionization_energy
    assert channels[1].ionization_energy == params.ionization_energy + 2 * E_CHARGE
    assert sample.key + "#upper" in channels[1].source
    with pytest.raises(DataError, match="not a level of He"):
        rates._channels(
            ionisation.ionisation_parameters("He"),
            {"fine_structure_level": sample.key + "#upper"},
            {sample.key: sample},
        )
    wrong = dataclasses.replace(
        sample, quantity="ionization_energy", parameters={"kind": "adiabatic"}
    )
    with pytest.raises(DataError, match="not a level energy"):
        rates._level_energy(wrong, 1)
    with pytest.raises(DataError, match="not an ionization_energy"):
        ionisation._energy(sample, 1)


def test_structure_coefficients_reject_wrong_quantities_labels_and_duplicates(sample):
    with pytest.raises(DataError, match="not structure coefficients"):
        ionisation.structure_coefficients(sample)
    table = dataclasses.replace(
        sample,
        material="Ar",
        quantity="tunnelling_structure_coefficient",
        unit="1",
        parameters={
            "model": "synthetic",
            "orbital": "p",
            "m": 0,
            "normalisation": "synthetic",
        },
    )
    with pytest.raises(DataError, match="partial wave l"):
        ionisation.structure_coefficients(
            dataclasses.replace(table, abscissae=("p", "s"))
        )
    with pytest.raises(DataError, match="twice"):
        ionisation.combined_structure_coefficients([table, table])
    parsed = ionisation.structure_coefficients(table)
    assert parsed.coefficients == {0: 0, 1: 2}


def test_configuration_subtables_and_vertical_energy_are_checked(sample):
    with pytest.raises(DataError, match="must be a table"):
        ionisation._strings({"reasons": []}, "reasons")
    with pytest.raises(DataError, match="structure must be a table"):
        ionisation._structures({"structure": []}, {})
    energy = dataclasses.replace(
        sample,
        material="Ar",
        quantity="ionization_energy",
        parameters={"kind": "adiabatic"},
    )
    reference = energy.key + "#upper"
    with pytest.raises(DataError, match="not a vertical energy"):
        ionisation._parameters(
            Gas.Ar,
            {"ionization_energy": reference, "vertical_ionization_energy": reference},
            {energy.key: energy},
        )
    assert rates._text({"reasons": []}, "reasons") == {}
    assert molecular._text([]) == {}


def test_rate_configuration_requires_atomic_parameters_and_curated_corrections(
    monkeypatch,
):
    model = {
        "gases": {"Ar": {"electrons": 1}},
        "structure": "hartree",
        "static_correction": "none",
        "coulomb": "static",
    }
    with pytest.raises(DataError, match="does not cover"):
        rates._build(Gas.He, "synthetic", model, {})
    with pytest.raises(DataError, match="not an atom"):
        rates._build(Gas.N2, "synthetic", model | {"gases": {"N2": {}}}, {})
    original = ionisation.ionisation_parameters("Ar")
    missing = dataclasses.replace(
        original, barrier_suppression_alpha=None, barrier_suppression_zhang=None
    )
    monkeypatch.setattr(rates, "ionisation_parameters", lambda gas: missing)
    for correction, message in [
        ("tong-lin", "no Tong-Lin constant"),
        ("zhang", "no Zhang constants"),
    ]:
        with pytest.raises(DataError, match=message):
            rates._build(
                Gas.Ar, "synthetic", model | {"static_correction": correction}, {}
            )
    monkeypatch.setattr(rates, "ionisation_parameters", lambda gas: original)
    built = rates._build(
        Gas.Ar, "synthetic", model | {"static_correction": "zhang"}, {}
    )
    assert built.static_constants == original.barrier_suppression_zhang
    assert (
        built.sources["static_correction"] == original.barrier_suppression_zhang_source
    )
    static = built.at_wavelength(math.inf)
    uncorrected = dataclasses.replace(
        static, static_correction="none", static_constants=()
    )
    fields = np.array([0.02, 0.04])
    energy = original.ionization_energy_au
    a1, a2, a3 = built.static_constants
    scaled = fields / (energy**2 / 4)
    # Zhang's clipped exponential correction, as documented in the function's equation.
    expected = np.exp(np.minimum(a1 * scaled**2 + a2 * scaled + a3, 0))
    np.testing.assert_allclose(
        static.instantaneous_rate(fields * AU_FIELD)
        / uncorrected.instantaneous_rate(fields * AU_FIELD),
        expected,
        rtol=1e-12,
        atol=0,
    )


def test_rate_parameters_reject_nonpositive_wavelength_and_incompatible_structure():
    model = rates.ionisation_rate_model("Ar", 800e-9)
    with pytest.raises(ValueError, match="wavelength must be positive"):
        model.at_wavelength(0.0)
    with pytest.raises(ValueError, match="structure_coefficient goes with"):
        dataclasses.replace(model, structure="hartree")
    assert model.wavelength_range[0] > 0
    assert model.wavelength_range[1] == math.inf


@pytest.mark.parametrize(
    "name", ["wavelength", "ionization_energy", "charge", "electrons"]
)
def test_molecular_rates_require_positive_physical_parameters(name):
    model = molecular.molecular_rate_model("N2", 800e-9)
    with pytest.raises(ValueError, match=f"{name} must be positive"):
        dataclasses.replace(model, **{name: 0.0})


def test_molecular_correction_constants_are_validated():
    model = molecular.molecular_rate_model("N2", 800e-9)
    with pytest.raises(ValueError, match="alpha must be positive"):
        dataclasses.replace(model, static_correction="tong-lin", tong_lin_alpha=0.0)
    with pytest.raises(ValueError, match="must not be negative"):
        dataclasses.replace(model, vibrational_constant=-1.0)
    with pytest.raises(DataError, match="is an atom"):
        molecular._build("synthetic", {}, Gas.Ar)
    with pytest.raises(DataError, match="unknown form"):
        molecular._build("synthetic", {"form": "missing"}, Gas.N2)
    with pytest.raises(DataError, match="no missing structure coefficients"):
        molecular._build("synthetic", {"form": "moadk", "structure": "missing"}, Gas.N2)


def test_missing_polarizability_and_dataset_references_are_explicit():
    original = ionisation.ionisation_parameters("Ar")
    assert (
        dataclasses.replace(original, ion_polarizability=None).polarizability_difference
        is None
    )
    with pytest.raises(DataError, match="unknown dataset"):
        ionisation._dataset({}, "missing")
    with pytest.raises(ValueError, match="exceeds l"):
        rates._angular_factor(0, 1)
    zeros = np.zeros((2, 3))
    np.testing.assert_array_equal(
        rates.ionisation_rate_model("Ar", 800e-9).instantaneous_rate(zeros), zeros
    )


def test_channel_sum_exhaustion_requests_the_integral(monkeypatch):
    # Lower only the iteration budget to force the documented convergence failure.
    monkeypatch.setattr(rates, "_SUM_MAXIMUM", 1)
    with pytest.raises(DataError, match="use channel_sum='integral'"):
        rates._channel_sum(
            0, np.array([1e-4]), np.array([1e-4]), np.array([2.5]), 1e-12
        )


def test_molecular_tunnelling_average_has_its_analytic_factor():
    model = dataclasses.replace(
        molecular.molecular_rate_model("N2", 800e-9), averaging="tunnelling"
    )
    field = np.array([0.02, 0.03])
    kappa = math.sqrt(2 * model.ionization_energy_au)
    # Saddle-point average of exp[-2*kappa^3/(3 F |cos theta|)].
    factor = np.sqrt(3 * field / (np.pi * kappa**3))
    np.testing.assert_allclose(
        model.cycle_averaged_rate(field * AU_FIELD)
        / model.instantaneous_rate(field * AU_FIELD),
        factor,
        rtol=1e-12,
        atol=0,
    )
    from optmat.exceptions import ModelRangeWarning  # noqa: PLC0415

    with pytest.warns(ModelRangeWarning, match="not tested below"):
        model.at_wavelength(100e-9).instantaneous_rate(0.02 * AU_FIELD)


def test_molecular_assembly_rejects_missing_structures_and_effective_charge(
    monkeypatch,
):
    params = ionisation.ionisation_parameters("N2")
    with pytest.raises(DataError, match="no missing structure"):
        molecular._build_molecular("synthetic", {"structure": "missing"}, Gas.N2, {})
    valid_structure = next(iter(params.structure))
    entry = {"ionization_energy": params.ionization_energy_source}
    with pytest.raises(DataError, match="unknown correction"):
        molecular._build_molecular(
            "synthetic",
            {"structure": valid_structure, "static_correction": "missing"},
            Gas.N2,
            entry,
        )
    monkeypatch.setattr(
        molecular,
        "ionisation_parameters",
        lambda gas: dataclasses.replace(params, effective_charge=None),
    )
    with pytest.raises(DataError, match="no effective charge"):
        molecular._build("synthetic", {"form": "ppt"}, Gas.N2)


def test_vibrational_validation_joins_only_shared_field_rows(
    sample, monkeypatch, curated
):
    fixed = dataclasses.replace(
        sample,
        source="majety_2015",
        id="table2_parallel_saenz_2000_H2",
        rows=("0.1", "0.2"),
        values=("2", "3"),
    )
    averaged = dataclasses.replace(
        sample,
        source="saenz_2002",
        id="tableI_a_v0",
        rows=("0.2", "0.3"),
        values=("6", "99"),
    )
    curated(fixed, averaged)
    fields, ratios = validation._saenz_ratio()
    np.testing.assert_array_equal(fields, [0.2])
    np.testing.assert_array_equal(ratios, [2.0])
