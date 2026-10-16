# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Golden values and identity checks (OEP 1, section 4.8 and Q32).

Gas Raman responses have derived floating-point identity fields. Compare their
canonical records numerically, retaining exact versions, classes and discrete
metadata; their public SHA-256 identities remain exact. Other models keep exact
digest comparisons. Evaluated outputs are checked independently in every case.
The samples are regression checks, not a proof of equivalence at every argument.
"""

import copy
import importlib.util
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pytest

from optmat.identity import sha256_digest

ROOT = Path(__file__).parents[1]
RECORDED = json.loads((ROOT / "tests" / "data" / "golden_models.json").read_text())
RAMAN_IDENTITIES = json.loads(
    (ROOT / "tests" / "data" / "golden_raman_identities.json").read_text()
)
# Arithmetic tolerance, not experimental uncertainty. Different libm/BLAS arithmetic
# changes the derived Raman widths: the investigated decay-rate difference was below
# 3e-14 relative. On an AVX-512 Xeon (Linux) the O3 fit's dispersion derivative, taken
# through the FFT of its absorption, differs from the macOS record by 1.2e-12 relative
# (2026-10-16): 1e-11, also for floating identity fields. No absolute floor: a
# formerly zero coefficient must stay zero, however small its unit.
RTOL = 1e-11
# Wider tolerances, each with its reason. O3's real refractivity is the Kramers-Kronig
# transform of its measured absorption, taken by FFT on a fine wavenumber grid, and its
# derivative amplifies that transform's rounding: Linux differs from the macOS record by
# up to 9.9e-9 relative (CI's runners; 7.4e-9 on an AVX-512 Xeon; 2026-10-16), far below
# the measurements' uncertainty. 1e-7 leaves a factor of ten.
LOOSER = {("gas/fit/O3", "dpolarizability_domega"): 1e-7}
# The canonical serializer represents these integer-valued arrays as floats.
EXACT_FIELDS = {"versions", "j", "damping_group"}
# A model's identity digest hashes its numbers, some of them computed (an FFT, a
# least-squares fit, a derived grid), so it is exact only where it was recorded: on
# Linux (OpenBLAS, glibc's libm) the O3 absorption's and the molecular rates' digests
# differ in the last bits of those numbers (CI, 2026-10-16). Elsewhere the values and
# the arguments are still held to RTOL, and the digest to its own record's.
RECORDED_ON = ("Darwin", "arm64")


def _digests_hold_here():
    return (platform.system(), platform.machine()) == RECORDED_ON


def _plain(value):
    return json.loads(json.dumps(value))


def _tool():
    path = ROOT / "tools" / "golden_models.py"
    spec = importlib.util.spec_from_file_location("golden_models", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("golden_models", module)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def fresh():
    return {(r["id"], r["quantity"]): r for r in _tool().records()}


def _decode(values):
    return np.asarray(
        values if not isinstance(values, list) else _floats(values), dtype=float
    )


def _floats(values):
    return [_floats(v) if isinstance(v, list) else float(v) for v in values]


def _assert_identity(actual, expected, path="identity"):
    assert type(actual) is type(expected), path
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), path
        for key in expected:
            if key in EXACT_FIELDS:
                assert actual[key] == expected[key], f"{path}.{key}"
            else:
                _assert_identity(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert len(actual) == len(expected), path
        for i, (now, before) in enumerate(zip(actual, expected, strict=True)):
            _assert_identity(now, before, f"{path}[{i}]")
    elif isinstance(expected, float):
        np.testing.assert_allclose(actual, expected, rtol=RTOL, atol=0, err_msg=path)
    else:
        assert actual == expected, path


def _assert_golden(record, now):
    assert now["class"] == record["class"]
    # Arguments derived on a grid (a geometric spacing) differ in the last bit too.
    _assert_identity(_plain(now["arguments"]), _plain(record["arguments"]), "arguments")
    if record["id"] in RAMAN_IDENTITIES:
        expected = RAMAN_IDENTITIES[record["id"]]
        assert sha256_digest(expected) == record["digest"], "recorded digest"
        assert sha256_digest(now["identity_record"]) == now["digest"], "runtime digest"
        _assert_identity(now["identity_record"], expected)
    elif _digests_hold_here():
        assert record.get("digest") == now.get("digest"), (
            f"{record['id']}: identity changed; review and re-record with "
            "tools/golden_models.py, with a changelog line"
        )
    np.testing.assert_allclose(
        _decode(now["values"]),
        _decode(record["values"]),
        rtol=LOOSER.get((record["id"], record["quantity"]), RTOL),
        atol=0,
        err_msg=(
            f"{record['id']}: output changed beyond arithmetic tolerance; review the "
            "EVALUATION_VERSION of the class or evaluator module that changed"
        ),
    )


def test_every_model_is_recorded(fresh):
    assert {(r["id"], r["quantity"]) for r in RECORDED} == set(fresh)
    assert set(RAMAN_IDENTITIES) == {
        r["id"] for r in fresh.values() if "identity_record" in r
    }


@pytest.mark.parametrize(
    "record", RECORDED, ids=[f"{r['id']}:{r['quantity']}" for r in RECORDED]
)
def test_golden_values(record, fresh):
    _assert_golden(record, fresh[(record["id"], record["quantity"])])


@pytest.fixture
def raman_pair():
    # Start with the reference record itself, so mutation tests do not depend on
    # which platform generated it or on the fitter's arithmetic on this machine.
    before = next(r for r in RECORDED if r["id"] == "raman_response/CH4@295K,800nm")
    now = copy.deepcopy(before)
    now["identity_record"] = copy.deepcopy(RAMAN_IDENTITIES[before["id"]])
    return before, now


def test_raman_roundoff_changes_exact_digest_but_passes_regression(raman_pair):
    before, now = raman_pair
    widths = now["identity_record"]["fields"]["width_broadening"]
    widths[0] = float(np.nextafter(widths[0], np.inf))
    now["digest"] = sha256_digest(now["identity_record"])
    assert now["digest"] != before["digest"]
    _assert_golden(before, now)


@pytest.mark.parametrize(
    ("field", "replacement"),
    [
        ("width_broadening", [4e8]),
        ("isotropic", [4e-37]),  # Tiny SI values still need a relative comparison.
        ("width_dicke", [1e-30]),  # A new nonzero term is not round-off about zero.
        ("band", ["nu2"]),
        ("j", [0.0]),
        ("damping_group", [1e-15]),  # Discrete even though serialized as a float.
        ("species", ["N2"]),
        ("sources", {}),
        ("reasons", {}),
        ("density_range", [0.0, "nan"]),
    ],
)
def test_raman_parameter_and_metadata_regressions_are_rejected(
    raman_pair, field, replacement
):
    before, now = raman_pair
    now["identity_record"]["fields"][field] = replacement
    now["digest"] = sha256_digest(now["identity_record"])
    with pytest.raises(AssertionError, match=field):
        _assert_golden(before, now)


@pytest.mark.parametrize("key", ["class", "versions"])
def test_raman_evaluation_identity_changes_are_rejected(raman_pair, key):
    before, now = raman_pair
    now["identity_record"][key] = "different" if key == "class" else {"evaluator": 2}
    now["digest"] = sha256_digest(now["identity_record"])
    with pytest.raises(AssertionError, match=key):
        _assert_golden(before, now)


def test_raman_output_regression_is_rejected_with_unchanged_identity(raman_pair):
    before, now = raman_pair
    now["values"][0] *= 1.00001
    with pytest.raises(AssertionError, match="output changed"):
        _assert_golden(before, now)


def test_raman_inconsistent_digest_is_rejected(raman_pair):
    before, now = raman_pair
    now["digest"] = "0" * 64
    with pytest.raises(AssertionError, match="runtime digest"):
        _assert_golden(before, now)


def test_other_models_still_require_exact_digests(monkeypatch):
    # Where the records were made; elsewhere the digest is not compared at all.
    monkeypatch.setattr(sys.modules[__name__], "_digests_hold_here", lambda: True)
    before = next(
        r for r in RECORDED if r["id"] not in RAMAN_IDENTITIES and "digest" in r
    )
    with pytest.raises(AssertionError, match="identity changed"):
        _assert_golden(before, {**before, "digest": "0" * 64})


def test_generator_writes_values_and_deduplicated_identities(
    tmp_path, monkeypatch, fresh
):
    tool = _tool()
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    monkeypatch.setattr(tool, "OUTPUT", tmp_path / "values.json")
    monkeypatch.setattr(tool, "RAMAN_OUTPUT", tmp_path / "identities.json")
    monkeypatch.setattr(tool, "records", lambda: copy.deepcopy(list(fresh.values())))
    tool.main()
    values = json.loads(tool.OUTPUT.read_text())
    identities = json.loads(tool.RAMAN_OUTPUT.read_text())
    assert set(identities) == set(RAMAN_IDENTITIES)
    assert len(values) == len(fresh)
    for record in values:
        assert "identity_record" not in record
        if record["id"] in identities:
            assert sha256_digest(identities[record["id"]]) == record["digest"]
        original = fresh[(record["id"], record["quantity"])]
        np.testing.assert_array_equal(
            _decode(record["values"]), _decode(original["values"])
        )
    first_values, first_identities = (
        tool.OUTPUT.read_bytes(),
        tool.RAMAN_OUTPUT.read_bytes(),
    )
    tool.main()
    assert tool.OUTPUT.read_bytes() == first_values
    assert tool.RAMAN_OUTPUT.read_bytes() == first_identities
