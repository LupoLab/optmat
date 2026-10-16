# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The silica nonlinear bundles (OEP 1, sections 4.7, 6.9, 7.10 and 8).

Each bundle's total n2 is the Kerr fit's and its parts add up to it; ``paired`` means
the Kerr fit's own f_R; each bundle's gain check against Stolen and Ippen is recorded,
within their stated uncertainty; the bundle binds its parts' identities.
"""

import dataclasses

import numpy as np
import pytest

from optmat.exceptions import DataError
from optmat.solid.bundle import silica_bundle, silica_bundles
from optmat.solid.kerr import kerr_models

WAVELENGTHS = np.array([526e-9, 800e-9, 1550e-9])


@pytest.mark.parametrize("name", sorted(silica_bundles()))
def test_the_parts_add_up_to_the_kerr_fits_total(name):
    bundle = silica_bundle(name)
    kerr = kerr_models()["kerr/SiO2"]
    np.testing.assert_array_equal(bundle.n2(WAVELENGTHS), kerr.n2(WAVELENGTHS))
    np.testing.assert_allclose(
        bundle.n2_electronic(WAVELENGTHS) + bundle.n2_raman(WAVELENGTHS),
        bundle.n2(WAVELENGTHS),
        rtol=1e-15,
    )
    if bundle.paired:
        assert bundle.raman_fraction == kerr.raman_fraction
    assert bundle.n2_electronic(1030e-9).shape == ()
    assert bundle.n2_raman(1030e-9).shape == ()


@pytest.mark.parametrize("name", sorted(silica_bundles()))
def test_every_bundle_records_its_gain_check(name):
    check = silica_bundle(name).gain_check
    assert check.reference == "stolen_1973/gain_amplifier"
    # Within one stated uncertainty of the measurement (10 %).
    assert abs(check.predicted - check.measured) < check.sigma


def test_the_gains_of_section_4_7():
    # OEP 1, section 4.7: 1.35 and 1.38e-11 cm/W, -9.7 % and -7.9 %.
    default = silica_bundle("default").gain_check
    tensor = silica_bundle("tensor").gain_check
    assert default.predicted * 100 == pytest.approx(1.35e-11, abs=0.005e-11)
    assert tensor.predicted * 100 == pytest.approx(1.38e-11, abs=0.005e-11)
    assert default.deviation == pytest.approx(-0.097, abs=5e-4)
    assert tensor.deviation == pytest.approx(-0.079, abs=5e-4)


def test_the_names_and_statuses():
    bundles = silica_bundles()
    assert set(bundles) == {
        "default",
        "tensor",
        "hollenbeck_2002",
        "lin_2006",
        "blow_1989",
    }
    assert bundles["default"].response.name == "hollenbeck_2002"
    assert bundles["default"].paired
    assert bundles["tensor"].response.name == "lin_2006"
    assert not bundles["tensor"].paired
    assert bundles["tensor"].raman_fraction == pytest.approx(0.245, abs=5e-4)
    assert {k for k, b in bundles.items() if b.status == "validated"} == {
        "default",
        "tensor",
        "hollenbeck_2002",
        "lin_2006",
    }
    with pytest.raises(DataError, match="no silica bundle"):
        silica_bundle("nonsense")


def test_the_response_and_its_parts():
    t = np.arange(0.0, 3e-12, 0.05e-15)
    for name in ("default", "tensor"):
        h = silica_bundle(name).h_r(t)
        assert np.trapezoid(h, t) == pytest.approx(1.0, abs=1e-3)
    silica_bundle("tensor").h_r(t, "b")
    with pytest.raises(DataError):
        silica_bundle("default").h_r(t, "b")


def test_the_identity_binds_the_parts():
    bundle = silica_bundle("default")
    record = bundle.identity_record()
    assert record["fields"]["kerr"]["digest"] == bundle.kerr.identity.digest
    assert record["fields"]["response"]["digest"] == bundle.response.identity.digest
    other = dataclasses.replace(bundle, response=silica_bundle("tensor").response)
    assert other.identity.digest != bundle.identity.digest
    assert bundle.identity.name == "silica_bundle/default"
    assert bundle.identity.kind == "silica_bundle"
