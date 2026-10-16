# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Every model of the surface has an identity and a status, owns read-only storage,
and declares its qualified intervals (OEP 1, sections 4.1, 4.8 and 6.2).
"""

import dataclasses
import pickle
import warnings
from collections.abc import Mapping

import numpy as np
import pytest

from optmat.exceptions import ModelRangeWarning
from optmat.gas.gases import Gas
from optmat.gas.ionisation_rates import ionisation_rate_models
from optmat.gas.models.formulas import FormulaModel as GasFormula
from optmat.gas.models.oscillator import OscillatorModel
from optmat.gas.models.registry import _all_models as gas_models
from optmat.gas.models.registry import model as gas_model
from optmat.gas.molecular_ionisation_rates import molecular_rate_models
from optmat.gas.nonlinear import hyperpolarizability_models
from optmat.gas.raman_response import raman_response
from optmat.identity import STATUSES, FrozenMapping, Identified
from optmat.solid.kerr import kerr_models
from optmat.solid.models._defaults import curated
from optmat.solid.models.formulas import FormulaModel as SolidFormula
from optmat.solid.models.registry import _all_models as solid_models
from optmat.solid.raman_response import silica_raman_responses


def _surface():
    yield from gas_models().values()
    yield from solid_models().values()
    yield from hyperpolarizability_models().values()
    for registry in (ionisation_rate_models(), molecular_rate_models()):
        for by_gas in registry.values():
            yield from by_gas.values()
    for gas in ("N2", "CO2", "dry air"):
        yield raman_response(gas, 295.0, wavelength=800e-9)
    yield from kerr_models().values()
    yield from silica_raman_responses().values()


@pytest.fixture(scope="module")
def surface():
    return list(_surface())


def _subclasses(cls):
    for sub in cls.__subclasses__():
        yield sub
        yield from _subclasses(sub)


def _oscillator(**changes):
    base = {
        "gas": Gas.Ar,
        "name": "synthetic",
        "citation": "synthetic",
        "wavelength_range": (200e-9, 2e-6),
        "kind": "fitted",
        "strengths": (1.0,),
        "energies": (1.0,),
        "labels": ("pole",),
    }
    return OscillatorModel(**(base | changes))


def _sellmeier(**changes):
    base = {
        "material": "SiO2",
        "name": "test/sellmeier",
        "citation": "",
        "wavelength_range": (210e-9, 3710e-9),
        "kind": "published",
        "form": "n2_sellmeier_wavelength",
        "terms": {"B": (0.6961663,), "L": (0.0684043,)},
    }
    return SolidFormula(**(base | changes))


def test_every_class_with_an_identity_declares_its_own_version():
    # The package's own classes; tests define toy subclasses of their own.
    classes = [c for c in _subclasses(Identified) if c.__module__.startswith("optmat.")]
    assert len(classes) >= 17
    for cls in classes:
        assert "EVALUATION_VERSION" in vars(cls), cls
        assert isinstance(cls.IDENTITY_KIND, str), cls


def test_every_model_has_a_distinct_identity_and_a_curated_status(surface):
    digests = set()
    for model in surface:
        identity = model.identity
        assert len(identity.digest) == 64
        assert identity.status in STATUSES
        assert identity.status != "user", identity.name
        assert identity.name
        assert identity.material
        digests.add(identity.digest)
    assert len(digests) == len(surface)


@pytest.mark.parametrize(
    ("get", "status"),
    [
        (lambda: gas_model("fit/N2"), "validated"),
        (lambda: gas_model("abjean_1971/sellmeier_Xe"), "comparison"),
        (lambda: molecular_rate_models()["default"][Gas.SF6], "untested"),
        (lambda: molecular_rate_models()["default"][Gas.N2], "validated"),
        (lambda: molecular_rate_models()["adk_hydrogenic"][Gas.N2], "comparison"),
        (lambda: ionisation_rate_models()["default"][Gas.Kr], "validated"),
        (lambda: ionisation_rate_models()["ppt_ilkov"][Gas.Ar], "comparison"),
        (lambda: raman_response("CO2", 295.0), "provisional"),
        (lambda: raman_response("dry air", 295.0), "validated"),
        (lambda: silica_raman_responses()["hollenbeck_2002"], "validated"),
        (lambda: silica_raman_responses()["lin_2006"], "comparison"),
        (lambda: kerr_models()["kerr/SiO2"], "validated"),
        (lambda: hyperpolarizability_models()[Gas.Ar], "validated"),
    ],
)
def test_statuses(get, status):
    assert get().status == status


def test_a_provisional_status_says_why():
    assert "TC-542" in raman_response("N2O", 295.0).status_reason
    assert "TC-678" in molecular_rate_models()["default"][Gas.SF6].status_reason


def test_a_material_without_a_default_is_compared_with_none():
    (only,) = curated({"x/y": _sellmeier(material="CaF2", name="x/y")}).values()
    assert only.status == "comparison"
    assert only.status_reason.startswith("kept for comparison")


def test_the_storage_of_every_model_is_read_only(surface):
    for model in surface:
        for field in dataclasses.fields(model):
            value = getattr(model, field.name)
            if isinstance(value, np.ndarray):
                assert not value.flags.writeable, (model.identity.name, field.name)
            if isinstance(value, Mapping):
                assert isinstance(value, FrozenMapping), (
                    model.identity.name,
                    field.name,
                )


@pytest.mark.parametrize(
    "get",
    [
        lambda: gas_model("fit/O3"),
        lambda: ionisation_rate_models()["default"][Gas.Ar],
        lambda: raman_response("N2", 295.0),
        lambda: silica_raman_responses()["hollenbeck_2002"],
        lambda: next(iter(solid_models().values())),
    ],
)
def test_pickling_keeps_the_digest_and_read_only_storage(get):
    model = get()
    copy = pickle.loads(pickle.dumps(model))
    assert copy.identity.digest == model.identity.digest
    for field in dataclasses.fields(copy):
        value = getattr(copy, field.name)
        if isinstance(value, np.ndarray):
            assert not value.flags.writeable


def test_changing_a_parameter_changes_the_digest():
    fit = gas_model("fit/Ar")
    assert isinstance(fit, OscillatorModel)
    moved = dataclasses.replace(fit, offset_au=fit.offset_au + 1e-12)
    assert moved.identity.digest != fit.identity.digest
    rate = ionisation_rate_models()["default"][Gas.Ar]
    assert rate.at_wavelength(400e-9).identity.digest != rate.identity.digest
    kerr = kerr_models()["kerr/SiO2"]
    other = dataclasses.replace(kerr, parameters=(*kerr.parameters[:3], 0.2))
    assert other.identity.digest != kerr.identity.digest


def test_the_components_a_model_reads_enter_its_digest():
    silica = silica_raman_responses()["hollenbeck_2002"]
    kerr = kerr_models()["kerr/SiO2"]
    record = silica.identity_record()
    assert record["extras"]["kerr_model"]["digest"] == kerr.identity.digest
    assert (
        record["extras"]["index_model"]["digest"] == silica.index_model.identity.digest
    )
    rebound = dataclasses.replace(silica)
    rebound.__dict__["kerr_model"] = dataclasses.replace(  # pyright: ignore[reportIndexIssue]
        kerr, parameters=(*kerr.parameters[:3], 0.2)
    )
    assert rebound.identity.digest != silica.identity.digest
    nitrogen = [
        m
        for m in solid_models().values()
        if getattr(m, "reference_gas", "") == "nitrogen"
    ]
    assert nitrogen
    for model in nitrogen:
        extras = model.identity_record()["extras"]
        assert (
            extras["reference_model"]["digest"] == gas_model("fit/N2").identity.digest
        )


def test_a_formula_resolves_its_reference_state_once():
    formula = gas_model("abjean_1971/sellmeier_Xe")
    assert isinstance(formula, GasFormula)
    extras = formula.identity_record()["extras"]
    assert extras["reference_density"] == formula.reference_density > 0
    assert extras["refractivity_virial"] == formula.refractivity_virial


def test_the_sources_of_a_model_without_records():
    assert _oscillator().identity.sources == ()
    assert _sellmeier(name="plain").identity.sources == ()
    assert _sellmeier().identity.sources == ("test",)


def test_the_material_of_a_raman_response():
    assert raman_response(Gas.N2, 295.0).material == "N2"
    mixture = raman_response({"N2": 0.8, "O2": 0.2}, 295.0)
    assert mixture.material == "N2+O2"
    assert mixture.identity.name == "raman_response/N2+O2@295K"


class TestIntervals:
    def test_every_curated_model_spans_its_range(self, surface):
        # One interval, the range, except the molecules fitted with OEP 3's windows.
        several = set()
        for model in surface:
            if hasattr(model, "intervals"):
                assert model.intervals[0][0] == model.wavelength_range[0]
                assert model.intervals[-1][1] == model.wavelength_range[1]
                if len(model.intervals) > 1:
                    several.add(model.identity.name)
        assert several == {
            "fit/CO",
            "fit/CO2",
            "fit/CH4",
            "fit/N2O",
        }

    def test_a_gap_between_intervals_warns(self):
        model = _oscillator(intervals=((200e-9, 500e-9), (600e-9, 2e-6)))
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            model.polarizability([300e-9, 800e-9, np.nan])
        with pytest.warns(ModelRangeWarning, match="200-500, 600-2000 nm"):
            model.polarizability(550e-9)

    def test_a_solid_gap_warns(self):
        model = _sellmeier(intervals=((210e-9, 1e-6), (2e-6, 3710e-9)))
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            model.permittivity([500e-9, 3e-6])
        with pytest.warns(ModelRangeWarning, match="210-1000, 2000-3710 nm"):
            model.permittivity(1.5e-6)

    @pytest.mark.parametrize(
        ("intervals", "match"),
        [
            (((500e-9, 200e-9),), "low < high"),
            (((200e-9, 700e-9), (600e-9, 2e-6)), "disjoint and ascending"),
            (((200e-9, 500e-9), (500e-9, 2e-6)), "disjoint and ascending"),
            (((300e-9, 2e-6),), "do not span"),
        ],
    )
    def test_bad_intervals_are_refused(self, intervals, match):
        with pytest.raises(ValueError, match=match):
            _oscillator(intervals=intervals)

    def test_the_intervals_enter_the_digest(self):
        split = _oscillator(intervals=((200e-9, 500e-9), (600e-9, 2e-6)))
        assert split.identity.digest != _oscillator().identity.digest
