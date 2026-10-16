# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Synthetic formula records exercise construction without altering curated sources."""

import copy
import math

import numpy as np
import pytest

import optmat.gas.models.formulas as gas_formulas
import optmat.gas.models.registry as gas_registry
import optmat.materials as registry
import optmat.solid.models.formulas as solid_formulas
from optmat.conditions import State
from optmat.constants import AU_ANGULAR_FREQUENCY, AU_POLARIZABILITY, C
from optmat.exceptions import DataError, OptmatError


@pytest.fixture(params=[gas_formulas, solid_formulas], ids=["gas", "solid"])
def formula_record(request):
    module = request.param
    return module, {
        "id": "synthetic",
        "material": "Ar" if module is gas_formulas else "SiO2",
        "form": "n2_sellmeier",
        "range_nm": [200, 2000],
        "state": "1atm_0C",
        "terms": {"B": ["0.50"], "C": ["0.01"]},
    }


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"range_nm": None}, "needs range_nm"),
        ({"range_nm": [200]}, "needs range_nm"),
        ({"terms": []}, "needs a terms table"),
    ],
)
def test_formula_records_require_terms_and_a_wavelength_range(
    formula_record, changes, message
):
    module, entry = formula_record
    with pytest.raises(DataError, match=message):
        module._model("test", "synthetic", entry | changes)


def _correction(**changes):
    return {
        "term": "B",
        "index": 0,
        "printed": "0.50",
        "corrected": "0.75",
        "reason": "invented typo",
        "evidence": "synthetic fixture",
    } | changes


@pytest.mark.parametrize(
    ("correction", "message"),
    [
        ({"term": "B"}, "correction lacks"),
        (_correction(term="A"), "missing coefficient"),
        (_correction(index=-1), "missing coefficient"),
        (_correction(index=1), "missing coefficient"),
        (_correction(printed="0.5"), "not the recorded"),
    ],
)
def test_corrections_check_evidence_index_and_exact_printing(
    formula_record, correction, message
):
    module, entry = formula_record
    with pytest.raises(DataError, match=message):
        module._model("test", "synthetic", entry | {"corrections": [correction]})


def test_a_formula_correction_changes_only_the_constructed_model(formula_record):
    module, entry = formula_record
    entry["corrections"] = [_correction()]
    original = copy.deepcopy(entry)
    model = module._model("test", "synthetic", entry)
    assert model.terms["B"] == (0.75,)
    assert entry == original
    assert model.name == "test/synthetic"
    assert model.citation == "synthetic"


@pytest.mark.parametrize(
    ("changes", "message"),
    [({"form": "imaginary"}, "unknown form"), ({"state": None}, "needs its state")],
)
def test_gas_formulas_reject_unknown_forms_and_missing_reference_states(
    changes, message
):
    entry = {
        "id": "test",
        "material": "Ar",
        "form": "dn_sigma2",
        "terms": {},
        "range_nm": [200, 2000],
        "state": "1atm_0C",
    } | changes
    with pytest.raises(DataError, match=message):
        gas_formulas._model("test", "synthetic", entry)


@pytest.mark.parametrize(
    "energy_unit,energy", [("hartree", 1), ("rydberg", 2), ("eV", 27.211386245988)]
)
def test_oscillator_energy_units_give_the_same_dynamic_and_static_response(
    energy_unit, energy
):
    # One invented oscillator: alpha(omega) = 2/(1-omega^2) au.
    entry = {
        "id": "test",
        "material": "Ar",
        "form": "oscillator",
        "range_nm": [200, "inf"],
        "energy_unit": energy_unit,
        "terms": {"f": "2", "E": str(energy)},
    }
    model = gas_formulas._model("test", "synthetic", entry)
    wavelengths = np.array([500e-9, 1e-6, np.inf])
    omega = 2 * math.pi * C / wavelengths / AU_ANGULAR_FREQUENCY
    # Only floating-point rounding in unit conversion and evaluation.
    np.testing.assert_allclose(
        model.polarizability(wavelengths) / AU_POLARIZABILITY,
        2 / (1 - omega**2),
        rtol=3e-15,
        atol=0,
    )


def test_polarizability_series_uses_even_frequency_powers():
    model = gas_formulas._model(
        "test",
        "synthetic",
        {
            "id": "test",
            "material": "Ar",
            "form": "alpha_series",
            "range_nm": [200, "inf"],
            "terms": {"c": ["2", "3", "5"]},
        },
    )
    wavelength = np.array([400e-9, 800e-9, np.inf])
    w = 2 * math.pi * C / wavelength / AU_ANGULAR_FREQUENCY
    np.testing.assert_allclose(
        model.polarizability(wavelength) / AU_POLARIZABILITY,
        2 + 3 * w**2 + 5 * w**4,
        rtol=3e-15,
        atol=0,
    )  # rounding only


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"terms": {"B": [1], "C": []}}, "equal length"),
        ({"terms": {"B": [1]}}, "equal length"),
        ({"wavelength_medium": "water"}, "unknown wavelength_medium"),
        ({"index_medium": "water"}, "index_medium must be"),
        ({"index_medium": "air"}, "index_medium must be"),
        ({"form": "n2_sellmeier_temperature", "terms": {"S1": [1]}}, "one length"),
        (
            {"form": "n2_sellmeier_temperature", "terms": {"S1": [1], "L1": [0.1]}},
            "needs temperature_range_k",
        ),
    ],
)
def test_solid_formula_metadata_is_validated(changes, message):
    entry = {
        "id": "test",
        "material": "SiO2",
        "form": "n2_sellmeier",
        "range_nm": [200, 2000],
        "terms": {"B": [1], "C": [0.01]},
    } | changes
    with pytest.raises(DataError, match=message):
        solid_formulas._model("test", "synthetic", entry)


def test_solid_temperature_metadata_and_reference_medium():
    entry = {
        "id": "test",
        "material": "SiO2",
        "form": "n2_sellmeier_temperature",
        "range_nm": [200, 2000],
        "terms": {"S1": [1, 0.001], "L1": [0.1, 0]},
        "temperature_range_k": [200, 400],
        "reference_temperature_k": 300,
        "index_medium": "nitrogen",
        "nitrogen_temperature_k": 290,
    }
    model = solid_formulas._model("test", "synthetic", entry)
    assert isinstance(model, solid_formulas.FormulaModel)
    assert model.air == State(101325.0, 290.0)
    assert model.reference_gas == "nitrogen"
    assert model.temperature == 300
    assert model.temperature_range == (200, 400)
    with pytest.raises(DataError, match="must lie in"):
        solid_formulas._model(
            "test", "synthetic", entry | {"reference_temperature_k": 401}
        )
    # With an absolute index, at 300 K the single-pole epsilon is 1 + 1.3*l^2/(l^2-.01).
    absolute = solid_formulas._model(
        "test", "synthetic", entry | {"index_medium": "vacuum"}
    )
    assert absolute.permittivity(500e-9) == pytest.approx(
        1 + 1.3 * 0.25 / 0.24, rel=2e-15
    )


@pytest.fixture
def temporary_registry(tmp_path, monkeypatch):
    path = tmp_path / "materials.toml"
    monkeypatch.setattr(registry, "MATERIALS_FILE", path)
    registry.materials.cache_clear()
    registry._names.cache_clear()
    yield path
    registry.materials.cache_clear()
    registry._names.cache_clear()


def test_registry_rejects_ambiguous_aliases_and_unknown_fields(temporary_registry):
    temporary_registry.write_text(
        '[A]\nname="a"\nkind="gas"\naliases=["shared"]\n'
        '[B]\nname="b"\nkind="gas"\naliases=["shared"]\n',
        encoding="utf-8",
    )
    with pytest.raises(OptmatError, match="claimed by A and B"):
        registry.as_material("shared")
    registry.materials.cache_clear()
    registry._names.cache_clear()
    temporary_registry.write_text(
        '[A]\nname="a"\nkind="gas"\nmisspelled="field"\n', encoding="utf-8"
    )
    with pytest.raises(OptmatError, match=r"unknown keys.*misspelled"):
        registry.materials()


def test_missing_recommendations_do_not_silently_choose_a_model(monkeypatch):
    monkeypatch.setattr(gas_registry, "_defaults", lambda: {})
    with pytest.raises(OptmatError, match="no dispersion model is recommended"):
        gas_registry.default_model("Ar")


def test_oscillator_terms_have_matching_lengths_positive_energies_and_valid_orders():
    from optmat.gas.gases import Gas  # noqa: PLC0415
    from optmat.gas.models.oscillator import OscillatorModel  # noqa: PLC0415

    common = {
        "gas": Gas.He,
        "name": "synthetic",
        "citation": "synthetic",
        "wavelength_range": (200e-9, math.inf),
        "kind": "fitted",
    }
    with pytest.raises(ValueError, match="differ in length"):
        OscillatorModel(**common, strengths=(1.0,), energies=(1.0,))
    with pytest.raises(ValueError, match="every energy must be positive"):
        OscillatorModel(**common, strengths=(1.0,), energies=(0.0,), labels=("pole",))
    valid = OscillatorModel(
        **common, strengths=(1.0,), energies=(1.0,), labels=("pole",)
    )
    with pytest.raises(ValueError, match="order must be non-negative"):
        valid.dispersion_au(np.array([0.1]), -1)


def test_comparison_formula_five_and_transparent_complex_fallback():
    from optmat.gas.comparison import DatabaseModel, _blocks  # noqa: PLC0415
    from optmat.gas.gases import Gas  # noqa: PLC0415

    model = DatabaseModel(
        gas=Gas.He,
        name="synthetic",
        citation="synthetic",
        wavelength_range=(200e-9, 2e-6),
        kind="comparison",
        formula_type="formula 5",
        coefficients=(1.0, 1e-5, -2.0),
        state=State(101325, 300),
    )
    wavelengths = np.array([500e-9, 1e-6])
    # Formula 5 here is n-1 = 1e-5/lambda_um^2 at the stated reference conditions.
    np.testing.assert_allclose(
        model.refractivity(wavelengths, model.state), [4e-5, 1e-5], rtol=1e-11, atol=0
    )
    complex_alpha = model._complex_polarizability(wavelengths)
    np.testing.assert_array_equal(np.imag(complex_alpha), [0, 0])
    np.testing.assert_array_equal(
        np.real(complex_alpha), model.polarizability(wavelengths)
    )
    with pytest.raises(DataError, match="expected one data block"):
        _blocks({"file": "synthetic", "data": []})
