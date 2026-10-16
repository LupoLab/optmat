# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The hyperpolarizability support records (OEP 5): derivation, facts and identity.

Each gas's record (``data/gas/models/nonlinear_support/<gas>.json``) is derived by
``tools/fit/nonlinear_support_record.py`` from the stored fit,
``data/gas/nonlinear_support.toml`` and the curated datasets. These tests check that the
stored records are fresh, that their facts reproduce their sources, that their extents
agree with the fitted values the manual's analysis (``tools/fit/nonlinear_support.py``,
which reads its limits from the records) enumerates with its own code, and that a record
has its own identity and leaves the model's untouched.
"""

import dataclasses
import functools
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from optmat.constants import HBAR, C
from optmat.exceptions import DataError
from optmat.gas import _nonlinear_support as derivation
from optmat.gas import nonlinear
from optmat.gas.gases import Gas
from optmat.gas.ionisation import ionisation_parameters
from optmat.gas.nonlinear import (
    AbsorptionEvidence,
    HyperpolarizabilitySupport,
    SupportLimit,
    SupportPoint,
    SupportUncertainty,
    frequency_sums,
    hyperpolarizability_models,
)
from optmat.identity import FrozenMapping

ROOT = Path(__file__).parents[1]
GASES = sorted(g.value for g in hyperpolarizability_models())
STORED = ROOT / "src" / "optmat" / "data" / "gas" / "models" / "nonlinear_support"
# Sums of four products of floats: they agree to rounding, about 1e-16 relative.
RTOL = 1e-12


def _tool():
    path = ROOT / "tools" / "fit" / "nonlinear_support.py"
    spec = importlib.util.spec_from_file_location("nonlinear_support", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("nonlinear_support", module)
    spec.loader.exec_module(module)
    return module


@functools.cache
def _record(gas: str) -> HyperpolarizabilitySupport:
    record = hyperpolarizability_models()[Gas(gas)].support
    assert record is not None
    return record


def _nm(wavenumber: float, photons: float = 1.0) -> float:
    """The vacuum wavelength of one photon of a wavenumber shared by n photons, nm."""
    return photons * 1e9 / wavenumber


def _printed(key: str) -> float:
    source, _, ident = key.partition("/")
    dataset = next(d for d in derivation.datasets() if d.key == key)
    assert dataset.source == source
    assert dataset.id == ident
    return float(dataset.values[0])


@pytest.fixture
def stored_directory(monkeypatch, tmp_path):
    """Point the loader at an empty directory, and forget the cached records."""
    monkeypatch.setattr(nonlinear, "_SUPPORT_DIRECTORY", tmp_path)
    nonlinear._stored_supports.cache_clear()
    yield tmp_path
    nonlinear._stored_supports.cache_clear()


# ------------------------------------------------------------------- freshness


@pytest.mark.parametrize("gas", GASES)
def test_the_stored_record_is_a_fresh_derivation(gas):
    fresh = nonlinear._support_dumps(derivation.derive_support(gas))
    assert (STORED / f"{gas}.json").read_text(encoding="utf-8") == fresh, (
        f"re-derive it: uv run python tools/fit/nonlinear_support_record.py {gas}"
    )


def test_every_fit_has_a_selection_and_a_record():
    assert sorted(derivation._selection()) == GASES
    assert sorted(p.stem for p in STORED.glob("*.json")) == GASES


@pytest.mark.parametrize("gas", GASES)
def test_a_stored_record_reads_back_to_the_same_digest(gas):
    record = _record(gas)
    again = nonlinear._support_loads(nonlinear._support_dumps(record))
    assert again.digest == record.digest
    assert again.identity.status == "validated"


# --------------------------------------------------------------- fitted extents


@pytest.mark.parametrize("gas", GASES)
def test_no_fitted_point_has_a_leg_beyond_the_derived_bound(gas):
    record = _record(gas)
    points = derivation.fitted_points(gas)
    assert len(points) == record.fitted_points
    for point, w2 in points:
        assert sum(math.fsum([x]) for x in point.legs) == pytest.approx(0.0, abs=1e-6)
        assert w2 == pytest.approx(sum(x**2 for x in point.legs), rel=RTOL, abs=0.0)
        assert max(abs(x) for x in point.legs) <= record.leg_max_in_w2 * (1 + RTOL)
    assert record.leg_max.wavenumber <= record.leg_max_in_w2 * (1 + RTOL)


def test_the_third_harmonic_attains_the_derived_bound():
    # Cauchy-Schwarz: nu_k^2 <= 3 W2 / 4, with equality for the third harmonic's output.
    legs = derivation.signed_legs("thg", np.array([1.0]))[:, 0]
    assert max(legs**2) == pytest.approx(0.75 * np.sum(legs**2), rel=RTOL, abs=0.0)
    rng = np.random.default_rng(5)
    inputs = rng.normal(size=(3, 4000))
    quartets = np.vstack([-inputs.sum(axis=0), inputs])
    assert np.all(np.max(quartets**2, axis=0) <= 0.75 * np.sum(quartets**2, axis=0))


@pytest.mark.parametrize("gas", GASES)
def test_the_extents_agree_with_the_analysis_of_the_manual(gas):
    # The analysis tool enumerates the fitted values with its own code.
    support = _tool().support(gas)
    record = _record(gas)
    assert len(support.fitted) == record.fitted_points
    assert record.leg_max.wavenumber == pytest.approx(
        100 * support.leg_max, rel=RTOL, abs=0.0
    )
    assert record.pair_max.wavenumber == pytest.approx(
        100 * support.pair_max, rel=RTOL, abs=0.0
    )
    assert record.max_w2 == pytest.approx(1e4 * support.max_w2, rel=RTOL, abs=0.0)


@pytest.mark.parametrize("gas", GASES)
def test_the_witnesses_attain_their_extents(gas):
    record = _record(gas)
    model = hyperpolarizability_models()[Gas(gas)]
    # W2 in SI: m^-2 is 1e4 cm^-2, exactly the stored fit's largest W2.
    assert record.max_w2 == 1e4 * model.record["max_w2_fitted_cm-2"]  # pyright: ignore[reportOperatorIssue]
    for point in record.max_w2_witnesses:
        assert sum(x**2 for x in point.legs) == pytest.approx(
            record.max_w2, rel=RTOL, abs=0.0
        )
    for point in record.leg_max.witnesses:
        assert max(map(abs, point.legs)) == pytest.approx(
            record.leg_max.wavenumber, rel=RTOL, abs=0.0
        )
    for point in record.pair_max.witnesses:
        assert derivation._pair_max(point.legs) == pytest.approx(
            record.pair_max.wavenumber, rel=RTOL, abs=0.0
        )
    assert record.leg_max.datasets == tuple(
        dict.fromkeys(p.dataset for p in record.leg_max.witnesses)
    )


def test_argons_largest_leg_is_the_harmonic_of_shelton_s_visible_eshg():
    record = _record("Ar")
    (witness,) = record.leg_max.witnesses
    assert witness.process == "eshg"
    assert witness.dataset == "shelton_1990/table3_Ar_literature"
    # The ESHG harmonic of 458.1 nm, as printed (0.05 nm, half the last digit).
    assert _nm(record.leg_max.wavenumber) == pytest.approx(229.0, abs=0.05)
    assert witness.legs == pytest.approx(
        (-2 / witness.wavelength, 1 / witness.wavelength, 1 / witness.wavelength, 0.0),
        rel=RTOL,
        abs=0.0,
    )


def test_the_windows_leave_out_the_calculated_harmonics_of_short_fundamentals():
    # Haettig & Joergensen's argon THG enters at 650 nm and longer (TC-316).
    points = [p for p, _ in derivation.fitted_points("Ar")]
    thg = [p for p in points if p.dataset == "hattig_1998/ccsd_thg_Ar"]
    assert thg
    # 694.3 nm as printed: half its last digit.
    assert min(p.wavelength for p in thg) == pytest.approx(694.3e-9, abs=0.05e-9)


def test_a_pump_probe_point_keeps_both_its_wavelengths():
    # No fitted value is static or CARS today; N2O's one value is a pump-probe index.
    (point,) = [p for p, _ in derivation.fitted_points("N2O")]
    assert point.process == "pump_probe"
    assert point.second_wavelength is not None
    nu, other = 1 / point.wavelength, 1 / point.second_wavelength
    legs = derivation.signed_legs("pump_probe", np.array([nu]), np.array([other]))
    assert point.legs == pytest.approx(tuple(legs[:, 0]), rel=RTOL, abs=0.0)
    assert point.second_wavelength == pytest.approx(580e-9, rel=RTOL, abs=0.0)


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
    legs = derivation.signed_legs(process, nu, second)
    w2, w4 = frequency_sums(process, nu, second)
    np.testing.assert_allclose(legs.sum(axis=0), 0.0, atol=1e-9)
    np.testing.assert_allclose(np.sum(legs**2, axis=0), w2, rtol=RTOL)
    np.testing.assert_allclose(np.sum(legs**4, axis=0), w4, rtol=RTOL)


@pytest.mark.parametrize(
    ("process", "other"), [("cars", None), ("pump_probe", None), ("sfg", None)]
)
def test_signed_legs_refuse_an_unknown_or_incomplete_process(process, other):
    with pytest.raises(DataError):
        derivation.signed_legs(process, np.array([1.0]), other)


def test_the_six_pair_sums_of_a_quartet_are_enumerated():
    legs = (-5.0, 1.0, 3.0, 1.0)
    pairs = [abs(legs[a] + legs[b]) for a in range(4) for b in range(a + 1, 4)]
    assert len(pairs) == 6
    assert derivation._pair_max(legs) == max(pairs) == 4.0


def test_ties_are_all_witnesses_in_order():
    points = [
        SupportPoint("a/b", str(i), "kerr", 1e-6, None, (-1.0, 1.0, 1.0, -1.0))
        for i in range(3)
    ]
    top, witnesses = derivation._extent(points, [2.0, 1.0, 2.0 * (1 - 1e-15)])
    assert top == 2.0
    assert [w.row for w in witnesses] == ["0", "2"]


# ------------------------------------------------------------------- the limits


def test_heliums_limits_are_solprops_sep_12():
    # SEP 12, section 4.4: two-photon 1s2s 1S at a 120.3 nm fundamental, three-photon
    # 1s2p 1P at 175 nm, three-photon ionisation below 151 nm (0.05 nm, as printed).
    record = _record("He")
    assert record.e2 is not None
    assert record.e1 is not None
    assert _nm(record.e2.wavenumber, 2) == pytest.approx(120.3, abs=0.05)
    assert _nm(record.e1.wavenumber, 3) == pytest.approx(175.3, abs=0.05)
    assert _nm(record.ionisation.wavenumber, 3) == pytest.approx(151.3, abs=0.05)


def test_xenons_two_photon_level_is_lehmbergs_resonance():
    # Lehmberg et al.: the two-photon resonance with 6p[1/2]0 at 249.6 nm (TC-854).
    record = _record("Xe")
    assert record.e2 is not None
    (j0,) = [c for c in record.e2.candidates if c.datasets[0].endswith("_even_J0")]
    assert _nm(j0.wavenumber, 2) == pytest.approx(249.6, abs=0.05)
    assert record.e2.wavenumber < j0.wavenumber


@pytest.mark.parametrize("gas", GASES)
def test_every_level_is_its_printed_value_in_si(gas):
    record = _record(gas)
    for limit in (record.e1, record.e2, record.e2_any):
        if limit is None or limit.kind != "level":
            continue
        (key,) = limit.datasets
        # cm^-1 to m^-1: times 100, exactly as the loader scales the printed string.
        assert limit.wavenumber == 100 * _printed(key)


@pytest.mark.parametrize(
    ("gas", "name", "chain"),
    [
        ("H2", "e2", ("H2_E", "H2_B")),
        ("N2", "e1", ("N2_b", "N2_a", "N2_ap")),
        ("N2", "e2", ("N2_a", "N2_ap")),
        ("O2", "e2", ("O2_C",)),
        ("CO", "e2_any", ("CO_I",)),
    ],
)
def test_a_band_origin_is_the_sum_of_its_chain_of_0_0_bands(gas, name, chain):
    limit = getattr(_record(gas), name)
    keys = tuple(f"nist_webbook_2026/state_{s}" for s in chain)
    assert limit.kind == "band_origin"
    assert limit.datasets == keys
    datasets = {d.key: d for d in derivation.datasets()}
    total = sum(float(datasets[k].parameters["nu00_cm1"]) for k in keys)  # pyright: ignore[reportArgumentType]
    assert limit.wavenumber == pytest.approx(100 * total, rel=RTOL, abs=0.0)
    assert not limit.uncertainty.stated
    assert any("hot-band" in q for q in limit.qualifications)


@pytest.mark.parametrize("gas", [g for g in GASES if g != "O3"])
def test_the_ionisation_energy_is_the_ionisation_modules(gas):
    parameters = ionisation_parameters(gas)
    record = _record(gas)
    assert record.ionisation.datasets == (parameters.ionization_energy_source,)
    hc = 2 * math.pi * HBAR * C
    assert record.ionisation.wavenumber == pytest.approx(
        parameters.ionization_energy / hc, rel=RTOL, abs=0.0
    )


def test_ozones_ionisation_is_the_lower_measured_threshold():
    record = _record("O3")
    # 12.66 eV over hc, to the printed 0.01 eV.
    assert record.ionisation.wavenumber * 1.239841984e-6 == pytest.approx(
        12.66, abs=0.005
    )


@pytest.mark.parametrize("gas", GASES)
def test_e2_of_any_polarisation_is_at_most_e2(gas):
    record = _record(gas)
    if record.e2 is not None:
        assert record.e2_any is not None
        assert record.e2_any.wavenumber <= record.e2.wavenumber


@pytest.mark.parametrize(
    ("gas", "e2", "e2_any"),
    [
        ("Ar", "level_Ar_even_J2", "level_Ar_even_J1"),
        ("Ne", "level_Ne_even_J2", "level_Ne_even_J1"),
        ("CO", "state_CO_A", "state_CO_I"),
    ],
)
def test_the_smallest_candidate_wins(gas, e2, e2_any):
    record = _record(gas)
    assert record.e2 is not None
    assert record.e2_any is not None
    assert record.e2.datasets[0].endswith(e2)
    assert record.e2_any.datasets[0].endswith(e2_any)
    assert record.e2_any.wavenumber == min(
        c.wavenumber for c in record.e2_any.candidates
    )


def test_uncertainties_are_stated_or_said_to_be_absent():
    neon, argon = _record("Ne"), _record("Ar")
    assert neon.e1 is not None
    assert argon.e1 is not None
    # NIST ASD's 0.0010 cm^-1 for neon's level, in m^-1.
    uncertainty = neon.e1.uncertainty
    assert isinstance(uncertainty, SupportUncertainty)
    assert uncertainty.stated
    assert uncertainty.value == pytest.approx(0.1, rel=RTOL, abs=0.0)
    assert uncertainty.level == "not stated"
    unstated = argon.e1.uncertainty
    assert isinstance(unstated, SupportUncertainty)
    assert not unstated.stated
    assert unstated.value is None
    assert "states none" in unstated.provenance
    assert neon.leg_max.uncertainty is None


@pytest.mark.parametrize(
    ("gas", "missing"),
    [
        ("CO2", {"e1", "e2", "e2_any"}),
        ("CH4", {"e2", "e2_any"}),
        ("SF6", {"e2", "e2_any"}),
        ("N2O", {"e2", "e2_any"}),
        ("O3", {"e2", "e2_any"}),
    ],
)
def test_the_polyatomic_molecules_state_what_is_missing(gas, missing):
    record = _record(gas)
    assert set(record.missing_reasons) == missing
    assert all("curated" in r for r in record.missing_reasons.values())
    if record.e1 is not None:
        assert record.e1.kind == "absorption_evidence"
        assert record.e1.uncertainty is None
        assert record.e1.evidence is not None
    assert record.band_edge() is None
    assert record.band_edge("any") is None


@pytest.mark.parametrize(
    ("gas", "detected", "nm"),
    [
        ("CH4", False, 197.7),
        ("SF6", False, 197.7),
        ("N2O", True, 217.5),
        ("O3", True, 1046.8),
    ],
)
def test_absorption_evidence_records_its_window(gas, detected, nm):
    limit = _record(gas).e1
    assert limit is not None
    assert limit.evidence is not None
    assert limit.evidence.detected is detected
    low, high = limit.evidence.sampled
    assert low <= limit.wavenumber <= high
    assert _nm(limit.wavenumber) == pytest.approx(nm, abs=0.05)


# ------------------------------------------------------------- derived quantities


@pytest.mark.parametrize("gas", GASES)
def test_the_band_edge_is_the_smaller_of_e1_and_half_e2(gas):
    record = _record(gas)
    for polarisation, e2 in (("linear", record.e2), ("any", record.e2_any)):
        edge = record.band_edge(polarisation)
        if record.e1 is None or e2 is None or record.e1.kind == "absorption_evidence":
            assert edge is None
        else:
            assert edge == min(record.e1.wavenumber, e2.wavenumber / 2)
    assert record.leg_max_in_w2 == math.sqrt(0.75 * record.max_w2)


@pytest.mark.parametrize(
    ("gas", "linear", "any_"),
    [
        ("He", 120.3, 120.3),
        ("Ne", 133.5, 134.9),
        ("Ar", 189.4, 192.1),
        ("Kr", 216.7, 219.4),
        ("Xe", 256.0, 258.8),
        ("H2", 201.7, 201.7),
        ("D2", 201.1, 201.1),
        ("N2", 290.1, 290.1),
        ("O2", 304.6, 304.6),
        ("CO", 308.9, 309.9),
    ],
)
def test_the_band_edges_are_the_manuals(gas, linear, any_):
    # The manual's table nonlinear-derived, to its printed 0.1 nm.
    record = _record(gas)
    assert _nm(record.band_edge()) == pytest.approx(linear, abs=0.05)  # pyright: ignore[reportArgumentType]
    assert _nm(record.band_edge("any")) == pytest.approx(any_, abs=0.05)  # pyright: ignore[reportArgumentType]


def test_an_unknown_polarisation_is_refused():
    with pytest.raises(DataError, match="polarisation"):
        _record("Ar").band_edge("circular")


def test_a_bound_on_every_leg_bounds_every_pair_including_the_outputs():
    rng = np.random.default_rng(7)
    inputs = rng.uniform(-1, 1, size=(3, 4000))
    quartets = np.vstack([-inputs.sum(axis=0), inputs])
    bound = np.max(np.abs(quartets), axis=0)  # the outputs included
    pairs = np.array(
        [quartets[a] + quartets[b] for a in range(4) for b in range(a + 1, 4)]
    )
    assert np.all(np.max(np.abs(pairs), axis=0) <= 2 * bound * (1 + RTOL))


def test_two_photons_of_equal_frequency_do_not_reach_an_even_j_1_level():
    # DeMille et al. 1999, eqs. 1-4: the antisymmetric part of the two-photon amplitude
    # is proportional to e_a x e_b and to 1/(W - w_a) - 1/(W - w_b), W the intermediate
    # level; it vanishes for equal frequencies, so the SPM quartet at the "any" edge of
    # argon, set by its even J = 1 level, is not a resonance of that level.
    def antisymmetric(level: float, w_a: float, w_b: float) -> float:
        return 1 / (level - w_a) - 1 / (level - w_b)

    record = _record("Ar")
    assert record.e1 is not None
    assert record.e2_any is not None
    assert record.e2_any.datasets[0].endswith("_even_J1")
    level, w = record.e1.wavenumber, record.e2_any.wavenumber / 2
    assert antisymmetric(level, w, w) == 0.0
    assert antisymmetric(level, 0.9 * w, 1.1 * w) != 0.0


# ---------------------------------------------------------------------- identity


def test_the_support_is_not_part_of_the_models_identity():
    model = hyperpolarizability_models()[Gas.Ar]
    assert "support" not in model.identity_record()["fields"]
    assert _record("Ar").model_digest == model.identity.digest
    assert _record("Ar").digest != model.identity.digest
    assert _record("Ar").identity.kind == "gas_hyperpolarizability_support"
    assert "nist_asd_2026" in _record("Ar").identity.sources


@pytest.mark.parametrize(
    "change",
    [
        lambda r: {"max_w2_witnesses": r.max_w2_witnesses[:0]},
        lambda r: {"e1": dataclasses.replace(r.e1, kind="lower_bound")},
        lambda r: {
            "e1": dataclasses.replace(
                r.e1, uncertainty=dataclasses.replace(r.e1.uncertainty, stated=True)
            )
        },
        lambda r: {"input_digests": {**r.input_digests, "selection": "0" * 64}},
        lambda r: {"e1": dataclasses.replace(r.e1, datasets=("other/dataset",))},
        lambda r: {"e2": None, "missing_reasons": {"e2": "not curated"}},
    ],
)
def test_any_change_to_a_record_changes_its_digest(change):
    record = _record("Ar")
    changed = dataclasses.replace(record, **change(record))
    assert changed.digest != record.digest


def test_a_changed_missing_reason_changes_the_digest():
    record = _record("CO2")
    changed = dataclasses.replace(
        record, missing_reasons={**record.missing_reasons, "e1": "another reason"}
    )
    assert changed.digest != record.digest


def test_a_record_and_its_parts_are_immutable():
    record = _record("Ar")
    with pytest.raises(dataclasses.FrozenInstanceError):
        record.max_w2 = 0.0  # pyright: ignore[reportAttributeAccessIssue]
    with pytest.raises(dataclasses.FrozenInstanceError):
        record.leg_max.wavenumber = 0.0  # pyright: ignore[reportAttributeAccessIssue]
    assert isinstance(record.leg_max.witnesses, tuple)
    assert isinstance(record.missing_reasons, FrozenMapping)


def test_a_users_model_has_no_support_and_a_curated_one_must():
    model = hyperpolarizability_models()[Gas.Ar]
    a, b, c, d = model.parameters
    user = dataclasses.replace(model, parameters=(a, b * 1.01, c, d), status="user")
    assert user.support is None
    same = dataclasses.replace(model, status="user")
    # The same polynomial is described by the same record, whatever its status.
    assert same.support is not None
    assert same.support.digest == _record("Ar").digest
    curated = dataclasses.replace(user, status="validated")
    with pytest.raises(DataError, match="re-derive"):
        _ = curated.support


def test_a_missing_shipped_record_is_an_error(stored_directory):
    with pytest.raises(DataError, match="re-derive"):
        _ = hyperpolarizability_models()[Gas.Ar].support


def _write(directory: Path, name: str, text: str) -> None:
    (directory / f"{name}.json").write_text(text, encoding="utf-8")


@pytest.mark.parametrize(
    ("gas", "edit", "message"),
    [
        ("CH4", lambda d: d.pop("max_w2"), "exactly the fields"),
        ("CH4", lambda d: d.update(derivation_version=0), "derivation version"),
        ("CH4", lambda d: d.update(leg_max=None), "always given"),
        ("Ar", lambda d: d["e1"].update(kind="guess"), "unknown support limit kind"),
        ("CH4", lambda d: d["e1"].update(evidence=None), "absorption evidence"),
        ("CH4", lambda d: d["missing_reasons"].update(e3="?"), "no limit"),
        ("CH4", lambda d: d["missing_reasons"].pop("e2"), "missing reason"),
        ("Ar", lambda d: d["e1"]["uncertainty"].update(extra=1), "exactly the fields"),
        ("CH4", lambda d: d["e1"]["evidence"].pop("detected"), "exactly the fields"),
    ],
)
def test_a_corrupted_shipped_record_is_an_error(stored_directory, gas, edit, message):
    data = json.loads((STORED / f"{gas}.json").read_text(encoding="utf-8"))
    edit(data)
    _write(stored_directory, gas, json.dumps(data))
    with pytest.raises(DataError, match=message):
        nonlinear._stored_supports()


def test_a_record_under_another_gas_name_is_an_error(stored_directory):
    _write(stored_directory, "He", (STORED / "Ar.json").read_text(encoding="utf-8"))
    with pytest.raises(DataError, match="holds the support record of Ar"):
        nonlinear._stored_supports()


def test_a_static_point_stores_its_infinite_wavelength_as_null():
    point = SupportPoint("a/b", "1", "static", math.inf, None, (0.0, 0.0, 0.0, 0.0))
    stored = nonlinear._plain(point)
    assert stored["wavelength"] is None
    assert nonlinear._point(stored, "test") == point


def test_evidence_and_only_evidence_carries_a_measurement():
    evidence = AbsorptionEvidence(True, (1.0, 2.0), "detected")
    with pytest.raises(DataError):
        SupportLimit(1.0, "level", (), "", "", (), evidence=evidence)
    with pytest.raises(DataError):
        SupportLimit(1.0, "absorption_evidence", (), "", "", ())


# --------------------------------------------------- the derivation refuses errors


def _models_with(gas: str, **record_changes):
    model = hyperpolarizability_models()[Gas(gas)]
    record = {**model.record, **record_changes}
    return {
        **hyperpolarizability_models(),
        Gas(gas): dataclasses.replace(model, record=record),
    }


def _zero_counts() -> dict[str, list[float]]:
    fitted = hyperpolarizability_models()[Gas.Ar].record["dataset_chi2"]
    assert isinstance(fitted, dict | FrozenMapping)
    return {str(k): [0, 0.0] for k in fitted}


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        (lambda: {"config_sha256": "0" * 64}, "not of the current configuration"),
        (lambda: {"max_w2_fitted_cm-2": 1.0}, "largest W2"),
        (lambda: {"dataset_chi2": _zero_counts()}, "disagree with the stored fit"),
    ],
)
def test_a_disagreement_with_the_stored_fit_is_an_error(monkeypatch, changes, message):
    models = _models_with("Ar", **changes())
    monkeypatch.setattr(derivation, "hyperpolarizability_models", lambda: models)
    with pytest.raises(DataError, match=message):
        derivation.fitted_points("Ar")


def test_an_unreadable_fitted_dataset_is_an_error_not_a_smaller_support(monkeypatch):
    def refuse(dataset, options):
        msg = f"{dataset.key}: cannot be converted"
        raise DataError(msg)

    monkeypatch.setattr(derivation, "gamma_observation", refuse)
    with pytest.raises(DataError, match="cannot be converted"):
        derivation.derive_support("Ar")


def test_legs_that_do_not_give_the_fits_w2_are_an_error(monkeypatch):
    real = derivation.signed_legs
    monkeypatch.setattr(derivation, "signed_legs", lambda *a: 2 * real(*a))
    with pytest.raises(DataError, match="do not give the fit's W2"):
        derivation.fitted_points("Ar")


def _selection_with(monkeypatch, gas: str, entry):
    selection = {**derivation._selection(), gas: entry}
    monkeypatch.setattr(derivation, "_selection", lambda: selection)


@pytest.mark.parametrize(
    ("edit", "message"),
    [
        (lambda e: e.update(extra=1), "no valid entry"),
        (lambda e: e["e1"].update(kind="guess"), "no known kind"),
        (lambda e: e["e1"].pop("candidates"), "needs candidates"),
        (
            lambda e: e["e1"].update(candidates=["nist_asd_2026/ie_Ar_I"]),
            "not a level_energy",
        ),
        (
            lambda e: e["e2"].update(kind="band_origin"),
            "not a electronic_state_energy",
        ),
        (
            lambda e: e.update(e1={"kind": "absorption_evidence", "dataset": "x/y"}),
            "needs exactly",
        ),
        (
            lambda e: e.update(
                e1={
                    "kind": "absorption_evidence",
                    "dataset": "no/such",
                    "at": "lowest",
                    "detected": True,
                    "registers": [],
                    "qualifications": [],
                }
            ),
            "no dataset",
        ),
    ],
)
def test_a_malformed_selection_is_an_error(monkeypatch, edit, message):
    entry = json.loads(json.dumps(derivation._selection()["Ar"]))
    edit(entry)
    _selection_with(monkeypatch, "Ar", entry)
    with pytest.raises(DataError, match=message):
        derivation.derive_support("Ar")


def test_a_gas_without_a_selection_is_an_error(monkeypatch):
    selection = {k: v for k, v in derivation._selection().items() if k != "Ar"}
    monkeypatch.setattr(derivation, "_selection", lambda: selection)
    with pytest.raises(DataError, match="no valid entry"):
        derivation.derive_support("Ar")


def test_a_level_must_be_one_value_in_cm1_and_a_band_needs_its_origin():
    datasets = {d.key: d for d in derivation.datasets()}
    key = "nist_asd_2026/level_Ar_odd_J1"
    wrong = {key: dataclasses.replace(datasets[key], unit="eV")}
    with pytest.raises(DataError, match="one level in cm-1"):
        derivation._level(wrong, key)
    band = "nist_webbook_2026/state_H2_B"
    parameters = {k: v for k, v in datasets[band].parameters.items() if k != "nu00_cm1"}
    bare = {band: dataclasses.replace(datasets[band], parameters=parameters)}
    with pytest.raises(DataError, match="no 0-0 band"):
        derivation._origin(bare, band)
