# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Explicit parameters and overrides (OEP 1, sections 6.3, 7.4, 7.5 and 8).

``from_parameters`` with a curated model's inputs reproduces it exactly, digest and
numbers; ``with_parameters(x=v)`` equals ``from_parameters`` of the changed inputs,
derived values included; a derived field is refused; an override records what it
changed, the curated digest and the status ``user``, and drops the stated
uncertainty. A custom species is a plain name.
"""

import dataclasses
import math
from typing import Any

import numpy as np
import pytest

from optmat.conditions import State
from optmat.constants import AU_POLARIZABILITY, E_CHARGE, HBAR, M_E
from optmat.exceptions import DataError, OptmatError
from optmat.gas.comparison import DatabaseModel
from optmat.gas.eos import equation_of_state
from optmat.gas.ionisation_rates import (
    Channel,
    IonisationRateModel,
    ionisation_rate_models,
)
from optmat.gas.models.fitted import fitted_models
from optmat.gas.models.formulas import FormulaModel as GasFormula
from optmat.gas.models.oscillator import OscillatorModel
from optmat.gas.models.registry import _all_models as gas_models
from optmat.gas.molecular_ionisation_rates import (
    MolecularRateModel,
    molecular_rate_models,
)
from optmat.gas.raman_response import GasRamanResponse, raman_response
from optmat.identity import describe
from optmat.solid.models.fitted import FittedModel as SolidFit
from optmat.solid.models.formulas import FormulaModel as SolidFormula
from optmat.solid.models.registry import _all_models as solid_models
from optmat.solid.models.tabulated import TabulatedModel
from optmat.solid.raman_response import silica_raman_responses

FIELDS = np.array([1e10, 3e10, 6e10])  # V/m


def _rate_inputs(m):
    return {
        "ionization_energy": min(c.ionization_energy for c in m.channels),
        "orbital_l": m.orbital_l,
        "electrons": m.electrons,
        "structure": m.structure,
        "structure_coefficient": m.structure_coefficient,
        "coulomb": m.coulomb,
        "static_correction": m.static_correction,
        "static_constants": m.static_constants,
        "channels": m.channels,
        "channel_sum": m.channel_sum,
        "averaging": m.averaging,
        "uncertainty_factor": m.uncertainty_factor,
    }


class TestTheMechanism:
    def test_an_override_records_itself(self):
        base = ionisation_rate_models()["default"][
            next(iter(ionisation_rate_models()["default"]))
        ]
        once = base.with_parameters(electrons=1.5)
        assert once.status == "user"
        assert once.identity.base == base.identity.digest
        assert once.identity.modified == (("electrons", "1.5"),)
        assert "overridden" in once.status_reason
        assert once.identity.digest != base.identity.digest
        twice = once.with_parameters(orbital_l=0)
        assert twice.identity.base == base.identity.digest
        assert twice.identity.modified == (("electrons", "1.5"), ("orbital_l", "0"))

    def test_names_are_checked(self):
        model = molecular_rate_models()["default"][
            next(iter(molecular_rate_models()["default"]))
        ]
        with pytest.raises(OptmatError, match="at least one change"):
            model.with_parameters()
        with pytest.raises(OptmatError, match="derived from the inputs"):
            model.with_parameters(structure_factors={0: 1.0})
        with pytest.raises(OptmatError, match="no input gas"):
            model.with_parameters(gas="Ar")

    def test_a_model_without_inputs_refuses_every_override(self):
        with pytest.raises(OptmatError, match="its inputs are none"):
            equation_of_state("Ar").with_parameters(fluid="Neon")

    def test_the_default_rebuild_replaces_the_fields(self):
        from typing import ClassVar  # noqa: PLC0415

        from optmat.identity import Identified  # noqa: PLC0415

        @dataclasses.dataclass(frozen=True)
        class Plain(Identified):
            IDENTITY_KIND: ClassVar[str] = "test"
            EVALUATION_VERSION: ClassVar[int] = 1
            INPUTS: ClassVar[tuple[str, ...]] = ("value",)
            value: float

            def _identity_label(self):
                return "plain", "x", ()

        changed = Plain(1.0).with_parameters(value=2.0)
        assert changed.value == 2.0
        assert changed.status == "user"

    def test_long_values_are_described_by_their_digest(self):
        assert describe(2.4) == "2.4"
        assert describe(np.arange(100.0)).startswith("sha256:")


class TestAtomicRates:
    @pytest.mark.parametrize("name", sorted(ionisation_rate_models()))
    def test_from_parameters_reproduces_every_curated_model(self, name):
        for m in ionisation_rate_models()[name].values():
            rebuilt = IonisationRateModel.from_parameters(
                m.gas,
                m.wavelength,
                model=m.model,
                sources=m.sources,
                reasons=m.reasons,
                **_rate_inputs(m),
            )
            assert rebuilt.identity.digest == m.identity.digest

    @pytest.mark.parametrize(
        ("name", "value"),
        [
            ("electrons", 1.7),
            ("orbital_l", 0),
            ("structure_coefficient", 2.4),
            ("coulomb", "ilkov"),
            ("static_constants", (7.0,)),
            ("channel_sum", "sum"),
            ("averaging", "tunnelling"),
            ("uncertainty_factor", 3.0),
            ("ionization_energy", 16.0 * E_CHARGE),
        ],
    )
    def test_an_override_is_from_parameters_of_the_changed_inputs(self, name, value):
        base = ionisation_rate_models()["default"][
            next(g for g in ionisation_rate_models()["default"] if g.value == "Ar")
        ]
        if name == "structure_coefficient":
            base = base.with_parameters(
                structure="tong_2002", structure_coefficient=2.0
            )
        changed = base.with_parameters(**{name: value})
        inputs = _rate_inputs(base) | {name: value}
        direct = IonisationRateModel.from_parameters(
            base.gas,
            base.wavelength,
            model=base.model,
            sources=base.sources,
            reasons=base.reasons,
            **inputs,
        )
        if name != "uncertainty_factor":
            direct = dataclasses.replace(direct, uncertainty_factor=math.nan)
        assert changed.identity.digest == direct.identity.digest
        np.testing.assert_array_equal(
            changed.cycle_averaged_rate(FIELDS), direct.cycle_averaged_rate(FIELDS)
        )

    def test_the_uncertainty_is_dropped_unless_given(self):
        base = ionisation_rate_models()["default"][
            next(iter(ionisation_rate_models()["default"]))
        ]
        assert math.isfinite(base.uncertainty_factor)
        assert math.isnan(base.with_parameters(electrons=1.1).uncertainty_factor)
        given = base.with_parameters(electrons=1.1, uncertainty_factor=4.0)
        assert given.uncertainty_factor == 4.0

    def test_the_energy_moves_every_channel(self):
        argon = next(
            m for g, m in ionisation_rate_models()["default"].items() if g.value == "Ar"
        )
        assert len(argon.channels) == 2
        split = (
            argon.channels[1].ionization_energy - argon.channels[0].ionization_energy
        )
        moved = argon.with_parameters(ionization_energy=16.0 * E_CHARGE)
        assert moved.channels[0].ionization_energy == 16.0 * E_CHARGE
        assert moved.channels[1].ionization_energy - moved.channels[
            0
        ].ionization_energy == pytest.approx(split, rel=1e-12, abs=0.0)

    def test_a_custom_species(self):
        model = IonisationRateModel.from_parameters(
            "Ar*",
            800e-9,
            ionization_energy=4.0 * E_CHARGE,
            orbital_l=0,
            electrons=1.0,
            static_correction="none",
            citation="an excited atom",
        )
        assert model.identity.material == "Ar*"
        assert model.identity.name == "rate/user/Ar*@800nm"
        assert model.citation.startswith("an excited atom; ")
        assert float(model.cycle_averaged_rate(1e9)) > 0

    @pytest.mark.parametrize(
        ("changes", "message"),
        [
            ({"orbital_l": -1}, "must not be negative"),
            ({"static_constants": (-1.0,)}, "alpha must be positive"),
            ({"channels": (Channel("x", 0.0, 1.0, ""),)}, "positive ionisation energy"),
        ],
    )
    def test_the_added_checks(self, changes, message):
        inputs = {
            "ionization_energy": 10 * E_CHARGE,
            "orbital_l": 1,
            "electrons": 1.0,
            "static_constants": (9.0,),
        } | changes
        if "channels" in changes:
            inputs["ionization_energy"] = 0.0
        with pytest.raises(ValueError, match=message):
            IonisationRateModel.from_parameters("X", 800e-9, **inputs)


class TestMolecularRates:
    @pytest.mark.parametrize("name", sorted(molecular_rate_models()))
    def test_from_parameters_reproduces_every_curated_model(self, name):
        for m in molecular_rate_models()[name].values():
            inputs = {
                k: getattr(m, k) for k in m.INPUTS if k != "structure_coefficients"
            }
            factors: dict[str, Any] = (
                {"structure_factors": m.structure_factors}
                if m.structure_factors
                else {}
            )
            rebuilt = MolecularRateModel.from_parameters(
                m.gas,
                m.wavelength,
                model=m.model,
                sources=m.sources,
                reasons=m.reasons,
                **inputs,
                **factors,
            )
            assert rebuilt.identity.digest == m.identity.digest

    def test_the_derived_values_are_rebuilt(self):
        ppt = next(iter(molecular_rate_models()["ppt_talebpour"].values()))
        higher = ppt.with_parameters(ionization_energy=ppt.ionization_energy * 1.1)
        assert higher.coefficient_squared != ppt.coefficient_squared
        direct = MolecularRateModel.from_parameters(
            ppt.gas,
            ppt.wavelength,
            form=ppt.form,
            charge=ppt.charge,
            ionization_energy=ppt.ionization_energy * 1.1,
            ionization_energy_kind=ppt.ionization_energy_kind,
            coulomb=ppt.coulomb,
            model=ppt.model,
            sources=ppt.sources,
            reasons=ppt.reasons,
        )
        assert higher.coefficient_squared == direct.coefficient_squared
        moadk = next(iter(molecular_rate_models()["moadk_tong_2002"].values()))
        coefficients = {0: 2.0, 2: 1.0}
        changed = moadk.with_parameters(structure_coefficients=coefficients)
        from optmat.gas._tunnelling import isotropic_structure_factors  # noqa: PLC0415

        assert dict(changed.structure_factors) == isotropic_structure_factors(
            coefficients
        )

    def test_with_ionization_energy_is_an_override(self):
        n2 = molecular_rate_models()["default"][
            next(g for g in molecular_rate_models()["default"] if g.value == "N2")
        ]
        other = n2.with_ionization_energy(15.0 * E_CHARGE, "a test")
        assert other.status == "user"
        assert dict(other.identity.modified)["ionization_energy_kind"] == '"a test"'

    @pytest.mark.parametrize(
        ("form", "structure", "message"),
        [
            ("moadk", {}, "takes the structure coefficients"),
            ("ppt", {"structure_coefficients": {0: 1.0}}, "takes no structure"),
            (
                "moadk",
                {"structure_coefficients": {0: 1.0}, "structure_factors": {0: 1.0}},
                "once",
            ),
        ],
    )
    def test_the_structure_is_checked(self, form, structure, message):
        with pytest.raises(ValueError, match=message):
            MolecularRateModel.from_parameters(
                "X2", 800e-9, form=form, ionization_energy=15 * E_CHARGE, **structure
            )


class TestGasDispersion:
    def test_from_parameters_is_si(self):
        f = [0.8, 2.0]
        energies = np.array([12.0, 25.0]) * E_CHARGE
        offset = -1e-43
        model = OscillatorModel.from_parameters(
            "X",
            f,
            energies,
            intervals=[(100e-9, math.inf)],
            offset=offset,
            citation="a test",
        )
        wavelength = 600e-9
        omega = 2 * math.pi * 299792458.0 / wavelength
        poles = energies / HBAR
        expected = (
            E_CHARGE**2
            / M_E
            * sum(fk / (wk**2 - omega**2) for fk, wk in zip(f, poles, strict=True))
            + offset
        )
        # The model goes through atomic units, the expectation through e^2 / m_e in SI:
        # 1.6e-12 relative apart.
        assert float(model.polarizability(wavelength)) == pytest.approx(
            expected, rel=1e-11, abs=0.0
        )
        assert model.offset_au == pytest.approx(offset / AU_POLARIZABILITY)
        assert model.labels == ("pole 1", "pole 2")

    def test_a_custom_species_is_an_ideal_gas_without_a_virial(self):
        model = OscillatorModel.from_parameters(
            "X", [1.0], [15 * E_CHARGE], intervals=[(100e-9, math.inf)]
        )
        assert model.identity.material == "X"
        assert model.refractivity_virial == 0.0
        state = State(1e5, 294.0)
        n = model.refractivity(800e-9, state)
        alpha = float(model.polarizability(800e-9))
        from optmat.gas.eos import ideal_number_density  # noqa: PLC0415
        from optmat.gas.lorentz_lorenz import refractivity  # noqa: PLC0415

        assert float(n) == pytest.approx(
            float(refractivity(alpha, ideal_number_density(state), 0.0))
        )
        np.testing.assert_allclose(
            np.real(model.refractive_index(800e-9, state)) - 1, n, rtol=1e-12
        )

    def test_the_intervals_set_the_range(self):
        argon = gas_models()["fit/Ar"]
        narrowed = argon.with_parameters(intervals=[(200e-9, 1e-6)])
        assert narrowed.wavelength_range == (200e-9, 1e-6)
        with pytest.raises(OptmatError, match="derived"):
            argon.with_parameters(wavelength_range=(1.0, 2.0))

    def test_the_oscillator_inputs(self):
        argon = fitted_models()["fit/Ar"]
        bigger = argon.with_parameters(strengths=tuple(2 * f for f in argon.strengths))
        low = 1e-6
        ratio = bigger.polarizability(low) / argon.polarizability(low)
        # Twice the oscillators; the magnetic offset (chi_0 < 0) is not doubled.
        assert float(ratio) == pytest.approx(2.0, rel=1e-4, abs=0.0)

    def test_the_fits_infrared_inputs(self):
        co = fitted_models()["fit/CO"]
        n = len(co.strengths) - co.n_fixed
        scaled = co.with_parameters(infrared_scale=2.0)
        np.testing.assert_array_equal(
            scaled.strengths[n:], np.array(co.strengths[n:]) * 2
        )
        assert scaled.strengths[:n] == co.strengths[:n]
        assert np.all(np.isnan(scaled.covariance))
        assert math.isnan(scaled.reduced_chi2)
        assert np.all(np.isnan(scaled.polarizability_sigma_au(np.array([0.05]))))
        replaced = co.with_parameters(infrared=([0.01], [0.01]))
        assert replaced.n_fixed == 1
        assert replaced.labels[-1] == "infrared 1"
        with pytest.raises(ValueError, match="number of terms is fixed"):
            co.with_parameters(strengths=co.strengths[:-1])
        with pytest.raises(OptmatError, match="derived"):
            co.with_parameters(n_fixed=0)

    def test_the_formula_and_database_inputs(self):
        formula = next(m for m in gas_models().values() if type(m) is GasFormula)
        name = next(iter(formula.terms))
        changed = formula.with_parameters(
            terms=dict(formula.terms)
            | {name: tuple(1.01 * v for v in formula.terms[name])}
        )
        assert changed.identity.digest != formula.identity.digest
        database = next(m for m in gas_models().values() if type(m) is DatabaseModel)
        other = database.with_parameters(
            coefficients=tuple(1.01 * c for c in database.coefficients)
        )
        assert other.status == "user"


class TestSolids:
    def test_the_sellmeier_constructor(self):
        silica = SolidFormula.sellmeier(
            "glass",
            {
                "A": (0.0,),
                "B": (0.6961663, 0.4079426, 0.8974794),
                "C": (0.0684043**2, 0.1162414**2, 9.896161**2),
            },
            (210e-9, 3.71e-6),
        )
        assert silica.identity.material == "glass"
        assert float(np.real(silica.refractive_index(1e-6))) == pytest.approx(
            1.450, abs=1e-3
        )
        relative = SolidFormula.sellmeier(
            "glass",
            {"B": (1.0,), "C": (0.01,)},
            (300e-9, 2e-6),
            reference_medium=State(101325.0, 293.15),
        )
        absolute = SolidFormula.sellmeier(
            "glass", {"B": (1.0,), "C": (0.01,)}, (300e-9, 2e-6)
        )
        assert float(np.real(relative.refractive_index(1e-6))) > float(
            np.real(absolute.refractive_index(1e-6))
        )
        with pytest.raises(ValueError, match="B and C of equal lengths"):
            SolidFormula.sellmeier("glass", {"B": (1.0,), "C": ()}, (300e-9, 2e-6))
        formula = next(m for m in solid_models().values() if type(m) is SolidFormula)
        assert (
            formula.with_parameters(
                intervals=[
                    (formula.wavelength_range[0], formula.wavelength_range[0] * 1.5)
                ]
            ).status
            == "user"
        )

    def test_the_table_constructor_and_overrides(self):
        table = TabulatedModel.from_arrays(
            "X",
            np.array([1e-6, 2e-6, 3e-6]),
            np.array([1.5, 1.4, 1.3]),
            np.array([1e-3, 2e-3, 3e-3]),
        )
        assert table.wavelength_range == (1e-6, 3e-6)
        with pytest.raises(ValueError, match="increasing wavelengths"):
            TabulatedModel.from_arrays(
                "X", np.array([2e-6, 1e-6]), np.ones(2), np.ones(2)
            )
        gold = solid_models()["table/Au"]
        assert isinstance(gold, TabulatedModel)
        shorter = gold.with_parameters(knots=gold.knots[:10], index=gold.index[:10])
        assert shorter.wavelength_range == (gold.knots[0], gold.knots[9])
        assert shorter.segments == ()
        chosen = shorter.with_parameters(intervals=[(gold.knots[0], gold.knots[5])])
        assert chosen.wavelength_range[1] == gold.knots[5]

    def test_the_fitted_udm(self):
        fit = next(m for m in solid_models().values() if type(m) is SolidFit)
        name = next(iter(fit.values))
        changed = fit.with_parameters(values={name: fit.values[name] * 1.01})
        assert changed.values[name] == pytest.approx(fit.values[name] * 1.01)
        assert changed.udm is not fit.udm
        assert np.all(np.isnan(changed.covariance))
        with pytest.raises(ValueError, match="no parameter"):
            fit.with_parameters(values={"nonsense": 1.0})
        low, high = fit.wavelength_range
        narrowed = fit.with_parameters(intervals=[(low, (low + high) / 2)])
        assert narrowed.udm is fit.udm
        assert narrowed.wavelength_range == (low, (low + high) / 2)


class TestRamanResponses:
    @pytest.mark.parametrize("gas", ["N2", "H2", "dry air", "CO2"])
    def test_from_lines_round_trips_a_curated_response(self, gas):
        response = raman_response(gas, 295.0, wavelength=800e-9)
        rebuilt = dataclasses.replace(
            GasRamanResponse.from_lines(**response._lines()),
            status=response.status,
            status_reason=response.status_reason,
        )
        assert rebuilt.identity.digest == response.identity.digest

    def test_width_scale_keeps_the_groups_and_moves_the_range(self):
        n2 = raman_response("N2", 295.0)
        narrow = n2.with_parameters(width_scale=0.8)
        np.testing.assert_array_equal(narrow.damping_group, n2.damping_group)
        np.testing.assert_allclose(narrow.width_broadening, 0.8 * n2.width_broadening)
        lines = n2._lines()
        for name in (
            "width_constant",
            "width_dicke",
            "width_broadening",
            "doppler_width",
        ):
            lines[name] = np.asarray(lines[name]) * 0.8
        direct = GasRamanResponse.from_lines(**lines)
        assert narrow.density_range == direct.density_range
        assert narrow.identity.digest == direct.identity.digest

    def test_a_width_law_regroups_its_band(self):
        n2 = raman_response("N2", 295.0)
        law = n2.with_parameters(width_law={"rotation": (1e9, 0.0, 2e9, 0.0)})
        rotation = np.array(law.band) == "rotation"
        assert np.all(law.width_broadening[rotation] == 2e9)
        assert len(set(law.damping_group[rotation])) == 1
        lines = n2._lines()
        for name, value in zip(
            ("width_constant", "width_dicke", "width_broadening", "doppler_width"),
            (1e9, 0.0, 2e9, 0.0),
            strict=True,
        ):
            lines[name] = np.where(rotation, value, lines[name])
        assert (
            law.identity.digest == GasRamanResponse.from_lines(**lines).identity.digest
        )
        with pytest.raises(ValueError, match="no band"):
            n2.with_parameters(width_law={"bending": (1.0, 0.0, 1.0, 0.0)})
        air = raman_response("dry air", 295.0)
        o2 = air.with_parameters(width_law={"O2 rotation": (1e9, 0.0, 2e9, 0.0)})
        chosen = (np.array(o2.species) == "O2") & (np.array(o2.band) == "rotation")
        assert np.all(o2.width_broadening[chosen] == 2e9)
        assert not np.any(o2.width_broadening[~chosen] == 2e9)

    def test_the_amplitude_scale(self):
        n2 = raman_response("N2", 295.0)
        np.testing.assert_allclose(
            n2.with_parameters(amplitude_scale=2.0).amplitude, 2 * n2.amplitude
        )

    def test_from_lines_checks_its_lines(self):
        common = ([1e13], [0.0], [1e-40], [1e9], [0.0], [1e8])
        with pytest.raises(ValueError, match="Dicke term needs its Doppler cap"):
            GasRamanResponse.from_lines(
                "X2",
                295.0,
                [1e13],
                [0.0],
                [1e-40],
                [1e9],
                [1e8],
                [1e8],
                [0.0],
                [0.0],
                mass=30.0,
            )
        with pytest.raises(ValueError, match="no molecular mass"):
            GasRamanResponse.from_lines("X2", 295.0, *common, [0.0], [0.0])
        with pytest.raises(ValueError, match="one entry per line"):
            GasRamanResponse.from_lines(
                "X2", 295.0, [1e13, 2e13], *common[1:], [0.0], [0.0], mass=30.0
            )
        mapped = GasRamanResponse.from_lines(
            "X2", 295.0, *common, [0.0], [0.0], mass={"X2": 30.0}, citation="a test"
        )
        assert mapped.identity.material == "X2"
        assert mapped.sources["citation"] == "a test"


class TestSilica:
    def test_an_override_unpairs_the_fraction(self):
        responses = silica_raman_responses()
        paired = next(r for r in responses.values() if r.paired)
        name = next(iter(paired.parameters))
        changed = paired.with_parameters(
            parameters={name: paired.parameters[name] * 1.01}
        )
        assert not changed.paired
        assert paired.with_parameters(raman_fraction=0.2).raman_fraction == 0.2
        with pytest.raises(OptmatError, match="derived"):
            paired.with_parameters(paired=True)
        with pytest.raises(DataError, match="no shape constant"):
            paired.with_parameters(parameters={"nonsense": 1.0})
