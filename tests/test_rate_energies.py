# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The channel rates and the energy per ionisation event (OEP 1, sections 4.6, 6.8, 8).

``channel_rates`` sum to ``instantaneous_rate``; ``energy_loss_rate`` is the channel
sum of energy times rate, 0 at zero field; ``energy_per_event`` is their ratio, the
lowest channel's energy at zero field and where the rates underflow, and for a single
channel the rate's own ionisation energy.
"""

import math

import numpy as np
import pytest

from optmat.constants import E_CHARGE
from optmat.gas.ionisation_rates import IonisationRateModel, ionisation_rate_models
from optmat.gas.molecular_ionisation_rates import molecular_rate_models

FIELDS = np.array([0.0, 5e9, 2e10, 5e10])  # V/m


def _atoms():
    return [
        m.at_wavelength(800e-9) for m in ionisation_rate_models()["default"].values()
    ]


@pytest.mark.parametrize("model", _atoms(), ids=lambda m: str(m.gas))
def test_the_channels_sum_to_the_rate(model):
    fields = FIELDS[model.maximum_field > FIELDS]
    channels = model.channel_rates(fields)
    assert channels.shape == (len(model.channels), len(fields))
    total = model.instantaneous_rate(fields)
    # The same terms, exponentiated one by one and summed: roundoff.
    np.testing.assert_allclose(channels.sum(axis=0), total, rtol=1e-12)
    loss = model.energy_loss_rate(fields)
    energies = np.array([c.ionization_energy for c in model.channels])
    np.testing.assert_allclose(loss, energies @ channels, rtol=1e-14)
    assert loss[0] == 0 and total[0] == 0
    per_event = model.energy_per_event(fields)
    np.testing.assert_allclose(per_event[1:], loss[1:] / total[1:], rtol=1e-12)
    assert per_event[0] == energies.min()


def test_argons_mean_energy_moves_with_the_field():
    argon = next(
        m for g, m in ionisation_rate_models()["default"].items() if g.value == "Ar"
    )
    energies = argon.at_wavelength(800e-9).energy_per_event([5e9, 5e10]) / E_CHARGE
    # OEP 1, section 4.6: 15.797 and 15.814 eV, between channels of 15.760 and 15.937.
    np.testing.assert_allclose(energies, [15.798, 15.814], atol=5e-4)


def test_the_mean_energy_is_defined_where_the_rates_underflow():
    # In a static field the tunnelling exponent underflows at weak fields.
    xenon = next(
        m for g, m in ionisation_rate_models()["default"].items() if g.value == "Xe"
    ).at_wavelength(math.inf)
    weak = np.array([1e7, 1e8])
    assert np.all(xenon.instantaneous_rate(weak) == 0)
    energies = xenon.energy_per_event(weak)
    lowest = min(c.ionization_energy for c in xenon.channels)
    assert np.all(np.isfinite(energies))
    # The higher channel's share is suppressed by its tunnelling exponent: at 1e7 V/m
    # the mean is the lowest channel's to far below 1e-12.
    assert energies[0] == pytest.approx(lowest, rel=1e-12, abs=0.0)


def test_a_single_channel_model_is_constant():
    for models in molecular_rate_models().values():
        for m in models.values():
            fields = np.array([0.0, 0.3, 0.9]) * m.maximum_field
            np.testing.assert_array_equal(
                m.energy_per_event(fields), m.ionization_energy
            )
            np.testing.assert_array_equal(
                m.channel_rates(fields)[0], m.instantaneous_rate(fields)
            )
            np.testing.assert_allclose(
                m.energy_loss_rate(fields),
                m.ionization_energy * m.instantaneous_rate(fields),
            )
    single = IonisationRateModel.from_parameters(
        "X",
        800e-9,
        ionization_energy=12 * E_CHARGE,
        orbital_l=1,
        electrons=1.0,
        static_constants=(9.0,),
    )
    assert np.all(single.energy_per_event([0.0, 1e10]) == 12 * E_CHARGE)


def test_the_energy_kinds():
    defaults = molecular_rate_models()["default"]
    assert {m.energy_kind for m in defaults.values()} == {"vertical"}
    atoms = ionisation_rate_models()["default"]
    assert {m.energy_kind for m in atoms.values()} == {"adiabatic"}
    argon = next(iter(atoms.values()))
    assert (
        argon.with_parameters(ionization_energy=16 * E_CHARGE).energy_kind
        == "not stated"
    )
    custom = IonisationRateModel.from_parameters(
        "X",
        800e-9,
        ionization_energy=12 * E_CHARGE,
        orbital_l=1,
        electrons=1.0,
        static_constants=(9.0,),
    )
    assert custom.energy_kind == "not stated"
    assert not math.isnan(argon.at_wavelength(math.inf).energy_per_event(1e10))


def test_the_channel_rates_warn_beyond_the_fitted_field():
    from optmat.exceptions import ModelRangeWarning  # noqa: PLC0415

    argon = next(iter(ionisation_rate_models()["default"].values()))
    with pytest.warns(ModelRangeWarning):
        argon.channel_rates(argon.maximum_field * 2)


def test_a_field_of_zero_alone():
    argon = next(iter(ionisation_rate_models()["default"].values()))
    assert np.all(argon.channel_rates(0.0) == 0)
    assert argon.channel_rates(0.0).shape == (len(argon.channels),)
    lowest = min(c.ionization_energy for c in argon.channels)
    assert float(argon.energy_per_event(0.0)) == lowest
