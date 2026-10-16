# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The Raman response's damping groups, line identifiers and tensor parts (OEP 1, 4.5).

Group members share their width law and shift exactly; the groups' responses sum to
the response; line identifiers are unique and stable; each group's decay rate has the
limit of section 4.5 at zero density; and the copolarised and cross-polarised elements
follow from the isotropic and anisotropic shares, as do solprop's (1, 0) and
(1/4, 3/4) weights.
"""

import importlib.util
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas.linewidths import AMAGAT
from optmat.gas.raman_response import raman_response

GASES = ("N2", "O2", "H2", "D2", "CH4", "SF6", "CO2", "N2O", "dry air")
SNAPSHOT = Path(__file__).parent / "data" / "raman_line_ids.json"


@pytest.fixture(scope="module", params=GASES)
def response(request):
    return raman_response(request.param, 295.0)


def test_group_members_share_their_law_exactly(response):
    for group in response.groups():
        k = np.array(group.lines)
        for name in (
            "width_constant",
            "width_dicke",
            "width_broadening",
            "doppler_width",
            "shift",
        ):
            values = getattr(response, name)[k]
            assert np.all(values == values[0]), name
    assert sorted(i for g in response.groups() for i in g.lines) == list(
        range(len(response))
    )


def test_the_groups_sum_to_the_response(response):
    low, high = response.density_range
    densities = np.geomspace(max(low, 1e22), min(high, 1e27), 100)
    t = np.linspace(0.0, 2e-12, 400)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        for density in densities[::11]:
            for component in ("parallel", "a"):
                total = response.response(t, density, component)
                summed = sum(
                    g.response(t, density, component) for g in response.groups()
                )
                scale = np.max(np.abs(total)) or 1.0
                # The same terms summed in another order: roundoff of the sum.
                assert np.max(np.abs(summed - total)) <= 1e-12 * scale * len(response)
    group = response.groups()[0]
    rates = group.decay_rate(densities)
    expected = math.pi * response._fwhm(densities / AMAGAT)[group.lines[0]]
    np.testing.assert_allclose(rates, expected, rtol=1e-14)


def test_the_line_identifiers_are_unique_and_stable(response):
    ids = response.line_id
    assert len(set(ids)) == len(ids)
    recorded = json.loads(SNAPSHOT.read_text())
    gas = next(
        g for g in GASES if raman_response(g, 295.0).material == response.material
    )
    assert list(ids) == recorded[gas]


def test_the_snapshot_is_the_tools():
    path = Path(__file__).parents[1] / "tools" / "raman_line_ids.py"
    spec = importlib.util.spec_from_file_location("raman_line_ids", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["raman_line_ids"] = module
    spec.loader.exec_module(module)
    assert module.snapshot() == json.loads(SNAPSHOT.read_text())


def test_the_zero_density_limit(response):
    for group in response.groups():
        limit = float(group.decay_rate(0.0))
        expected = group.gamma_constant + (
            group.gamma_doppler if group.gamma_dicke > 0 else 0.0
        )
        assert limit == expected
        # The law tends to its limit: at 1e10 m^-3 the Dicke term is at its cap and
        # only the broadening, gamma_beta N, remains above it.
        near = float(group.decay_rate(1e10))
        assert near == pytest.approx(limit + group.gamma_broadening * 1e10, rel=1e-12)


def test_the_limits_section_4_5_states():
    for gas in ("N2", "O2", "dry air", "N2O", "CO2"):
        assert all(
            float(g.decay_rate(0.0)) == 0 for g in raman_response(gas, 295.0).groups()
        )
    d2 = [float(g.decay_rate(0.0)) for g in raman_response("D2", 295.0).groups()]
    assert min(d2) == 0 < max(d2)
    h2 = [
        float(g.decay_rate(0.0)) / math.pi for g in raman_response("H2", 295.0).groups()
    ]
    assert min(h2) == pytest.approx(152.8e6, rel=1e-3)
    for gas in ("CH4", "SF6"):
        assert all(
            float(g.decay_rate(0.0)) > 0 for g in raman_response(gas, 295.0).groups()
        )


def test_only_hydrogens_q_lines_shift():
    for gas in GASES:
        shifted = [g for g in raman_response(gas, 295.0).groups() if g.shifted]
        if gas == "H2":
            assert [g.line_ids[0].split("/")[1:3] for g in shifted] == [["v1", "Q"]]
        else:
            assert shifted == []


def test_a_group_without_a_shift_has_a_fixed_row():
    n2 = raman_response("N2", 295.0)
    group = n2.groups()[3]
    assert not group.shifted
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        frequencies = n2.frequencies(np.array([0.5, 5.0]) * AMAGAT)
    k = list(group.lines)
    np.testing.assert_array_equal(frequencies[k, 0], frequencies[k, 1])


def test_the_tensor_elements(response):
    iso, aniso = response.isotropic, response.anisotropic
    # Copolarised (x, x, x, x): a + b = iso + aniso; cross-polarised (x, x, y, y): a.
    np.testing.assert_array_equal(
        response.component_amplitudes("parallel"), iso + aniso
    )
    np.testing.assert_allclose(response.component_amplitudes("a"), iso - aniso / 2)
    # solprop's transverse trace and traceless weights, (1, 0) and (1/4, 3/4): their
    # sum is the copolarised element and their difference the cross-polarised one.
    trace = iso * 1 + aniso * 0.25
    traceless = iso * 0 + aniso * 0.75
    np.testing.assert_allclose(trace + traceless, iso + aniso)
    np.testing.assert_allclose(trace - traceless, iso - aniso / 2)


def test_densities_are_checked():
    n2 = raman_response("N2", 295.0)
    group = n2.groups()[0]
    with pytest.raises(DataError, match=">= 0"):
        group.decay_rate(-1.0)
    with pytest.raises(DataError, match=">= 0"):
        n2.check_densities([math.nan])
    low, high = n2.density_range
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        n2.check_densities([low * 2, high / 2])
        group.decay_rate(0.0)
    with pytest.warns(ModelRangeWarning, match="density below"):
        n2.check_densities([0.0])
    with pytest.warns(ModelRangeWarning, match="density above"):
        n2.check_densities([high * 2])
