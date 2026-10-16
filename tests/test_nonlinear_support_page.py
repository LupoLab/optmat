# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The support of the hyperpolarizability fits: the tool, its claims and the tables.

``tools/fit/nonlinear_support.py`` computes the legs of the fitted data, the resonance
limits, the monotonicity of the fits and the per-leg caps (TC-864 to TC-867). These
tests check its numerics against the module and the analytic statements of the manual,
and that the tables of *Third-order nonlinearity* and *Third-order accuracy* are its
output (regenerate them with ``tools/fit/nonlinear_support.py --update-docs``).
"""

import functools
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from optmat.gas.nonlinear import frequency_sums, hyperpolarizability_models

ROOT = Path(__file__).parents[1]


def _tool():
    path = ROOT / "tools" / "fit" / "nonlinear_support.py"
    spec = importlib.util.spec_from_file_location("nonlinear_support", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("nonlinear_support", module)
    spec.loader.exec_module(module)
    return module


T = _tool()


@functools.cache
def _result() -> dict:
    result = {gas: T.analyse(gas) for gas in T.GASES}
    result["cases"] = [T.case(*c) for c in T.CASES]
    # The tables read the result after a JSON round trip, as the tool writes it.
    return json.loads(json.dumps(result, default=float))


def _block(page: str, name: str) -> str:
    text = (ROOT / "docs" / page).read_text(encoding="utf-8")
    return text.split(f"<!-- {name}:start -->\n")[1].split(f"\n<!-- {name}:end -->")[0]


@pytest.mark.parametrize(
    ("process", "other"),
    [
        ("static", None),
        ("kerr", None),
        ("thg", None),
        ("eshg", None),
        ("dc_kerr", None),
        ("cars", 9000.0),
        ("pump_probe", 17241.0),
    ],
)
def test_the_signed_legs_sum_to_zero_and_give_the_modules_sums(process, other):
    nu = np.array([12500.0, 21831.0])
    second = None if other is None else np.full_like(nu, other)
    legs = T.signed_legs(process, nu, second)
    w2, w4 = frequency_sums(process, nu, second)
    # Exact in floating point up to rounding of sums of four squares: 1e-12 relative.
    np.testing.assert_allclose(legs.sum(axis=0), 0.0, atol=1e-9)
    np.testing.assert_allclose(np.sum(legs**2, axis=0), w2, rtol=1e-12)
    np.testing.assert_allclose(np.sum(legs**4, axis=0), w4, rtol=1e-12)


def test_a_pair_and_its_complement_have_the_same_sum():
    legs = T.signed_legs("cars", np.array([18797.0]), np.array([15000.0]))
    sums = sorted(
        abs(legs[a, 0] + legs[b, 0]) for a in range(4) for b in range(a + 1, 4)
    )
    np.testing.assert_allclose(sums[0::2], sums[1::2], rtol=1e-12)
    assert T.pair_sums(legs)[0] == pytest.approx(sums[-1], rel=1e-12, abs=0.0)


def test_no_quartet_within_a_w2_has_a_leg_beyond_the_third_harmonics():
    # Cauchy-Schwarz: nu_j^2 <= 3 W2 / 4, with equality for the third harmonic.
    rng = np.random.default_rng(1)
    inputs = rng.normal(size=(3, 5000))
    legs = np.vstack([-inputs.sum(axis=0), inputs])
    assert np.all(
        np.max(legs**2, axis=0) <= 0.75 * np.sum(legs**2, axis=0) * (1 + 1e-12)
    )
    thg = T.signed_legs("thg", np.array([1.0]))
    assert np.max(thg**2) == pytest.approx(0.75 * np.sum(thg**2), rel=1e-12, abs=0.0)


@pytest.mark.parametrize("gas", T.GASES)
def test_every_fit_increases_with_each_squared_leg(gas):
    model = hyperpolarizability_models()[T._gas(gas)]
    _, a, b, b_prime = model.parameters
    # The manual's bound: A > 0, B >= 0 and B + B' > 0 make the slope at least A.
    assert a > 0
    assert b >= 0
    assert b + b_prime > 0
    w2_max = float(model.record["max_w2_fitted_cm-2"])  # pyright: ignore[reportArgumentType]
    assert T.slope_minimum(model, 4 * w2_max) == pytest.approx(a, rel=1e-12, abs=0.0)


@pytest.mark.parametrize("gas", T.GASES)
def test_the_cap_at_the_third_harmonics_leg_leaves_the_data_and_the_fitted_w2_exact(
    gas,
):
    result = _result()[gas]["caps"]["b'"]
    assert result["fitted_changed"] == 0
    assert result["in_w2_worst"] == 0.0


def test_argons_largest_fitted_leg_is_the_harmonic_of_shelton_s_visible_eshg():
    # Not Haettig & Joergensen's THG through the 650 nm window, which keeps 694.3 nm and
    # longer: harmonics at 231.4 nm and longer (TC-864).
    argon = _result()["Ar"]
    source = argon["leg_max_source"]
    assert source["key"] == "shelton_1990/table3_Ar_literature"
    assert source["process"] == "eshg"
    assert argon["leg_max_nm"] == pytest.approx(229.0, abs=0.05)
    thg = [v for v in argon["fitted_values"] if v["key"] == "hattig_1998/ccsd_thg_Ar"]
    assert min(v["wavelength_nm"] for v in thg) == pytest.approx(694.3, abs=0.05)
    assert max(max(abs(x) for x in v["legs"]) for v in thg) == pytest.approx(
        3e7 / 694.3, rel=1e-4, abs=0.0
    )


def test_heliums_limits_are_solprops_sep_12():
    # SEP 12, section 4.4: two-photon 1s2s 1S at a 120.3 nm fundamental, three-photon
    # 1s2p 1P at 175 nm, three-photon ionisation below 151 nm (to 0.05 nm as printed).
    helium = _result()["He"]["resonances"]
    assert 2e7 / helium["E2"] == pytest.approx(120.3, abs=0.05)
    assert 3e7 / helium["E1"] == pytest.approx(175.3, abs=0.05)
    assert 3e7 / helium["IE"] == pytest.approx(151.3, abs=0.05)


def test_the_kerr_diagonal_cap_moves_argons_266_nm_quartet_beyond_the_carrier():
    (case,) = [c for c in _result()["cases"] if c["name"].startswith("Ar 400")]
    exact = case["uncapped_au"]
    assert case["w2_over_max"] < 1
    assert case["capped_au"]["b'"] == exact
    kerr_diagonal = case["capped_au"]["a"] / exact - 1
    carrier = case["carrier_au"]["400.0"] / exact - 1
    # -9.3 % against -4.3 % (the manual, to 0.1 %).
    assert kerr_diagonal == pytest.approx(-0.093, abs=1e-3)
    assert carrier == pytest.approx(-0.043, abs=1e-3)


def test_heliums_spm_at_120_nm_is_past_the_two_photon_resonance():
    (case,) = [c for c in _result()["cases"] if c["name"] == "He SPM at 120 nm"]
    assert case["past_resonance"]
    assert case["pair_over_e2"] > 1


@pytest.mark.parametrize(
    ("page", "name"),
    [
        ("user-guide/nonlinear-support.md", "nonlinear-support"),
        ("user-guide/nonlinear-support.md", "nonlinear-resonances"),
        ("user-guide/nonlinear-support.md", "nonlinear-derived"),
        ("user-guide/nonlinear-extrapolation.md", "nonlinear-caps"),
        ("user-guide/nonlinear-extrapolation.md", "nonlinear-carrier"),
        ("user-guide/nonlinear-extrapolation.md", "nonlinear-beyond"),
        ("user-guide/nonlinear-extrapolation.md", "nonlinear-cases"),
    ],
)
def test_the_pages_tables_are_the_tools(page, name):
    assert _block(page, name) == T.render(_result())[name]


def test_the_derived_band_edge_is_the_first_resonance_of_an_spm_quartet():
    for gas in T.GASES:
        result = _result()[gas]
        edge = result["derived"]["band_edge"]
        limits = result["resonances"]
        if limits["E1"] is None or limits["E2"] is None:
            assert edge is None
            continue
        # A leg at E1 or an SPM pair 2 nu at E2, whichever comes first.
        assert edge == pytest.approx(
            min(limits["E1"], limits["E2"] / 2), rel=1e-12, abs=0.0
        )
        assert edge <= limits["E1"]
        assert 2 * edge <= limits["E2"] * (1 + 1e-12)


def test_every_table_has_its_own_rows():
    # Each row has the header's number of cells, and the quartet table lists the
    # quartets: guards against one table's rows rendered under another's header.
    tables = T.render(_result())
    for name, table in tables.items():
        header, _, *rows = table.splitlines()
        width = header.count("|")
        assert all(row.count("|") == width for row in rows), name
    quartets = [
        row.split(" | ")[0].lstrip("| ")
        for row in tables["nonlinear-cases"].splitlines()[2:]
    ]
    assert quartets == [c[0] for c in T.CASES]


def test_e2_of_any_polarisation_is_at_most_the_scalar_e2():
    # E2_any is the minimum over all eligible two-photon levels, a superset of the
    # scalar field's (TC-869).
    for gas in T.GASES:
        limits = _result()[gas]["resonances"]
        if limits["E2"] is not None:
            assert limits["E2_any"] <= limits["E2"]


def test_absorption_evidence_never_forms_a_derived_band_edge():
    for gas in ("CH4", "SF6", "N2O", "O3"):
        result = _result()[gas]
        assert result["resonances"]["E1_kind"] != "level"
        assert result["derived"]["band_edge"] is None
