# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The bound two-photon level records (OEP 6, Phase 2a): derivation, convention, checks.

Each record (``data/gas/models/two_photon_bound/<gas>.json``) is derived by
``tools/fit/two_photon_bound_record.py`` from ``data/gas/two_photon_bound.toml``, the
NIST level lists held as returned and the curated datasets. These tests check that the
records are fresh, the convention of TC-952 (sigma_0 = 2 pi^3 alpha^2 a0^4 omega_au^2
delta, per level), the agreements the assessments rest on, and the refusals.
"""

import dataclasses
import json
import math

import pytest
from scipy import constants

from optmat.constants import AU_ANGULAR_FREQUENCY
from optmat.exceptions import DataError, UnavailableModelError
from optmat.gas import _two_photon_bound_derivation as derivation
from optmat.gas import two_photon_bound
from optmat.gas.two_photon_bound import (
    SIGMA0_PER_DELTA,
    BoundLevel,
    sigma0_from_delta,
    two_photon_bound_levels,
)

GASES = ("He", "Ne", "Ar", "Kr", "Xe")
XE_6P = "5p5.(2P*<3/2>).6p"
XE_6P1 = "5p5.(2P*<1/2>).6p"


def _sigma0(record, dataset: str, level: BoundLevel) -> float:
    (strength,) = [s for s in record.strengths_of(level) if s.dataset == dataset]
    assert strength.sigma0 is not None
    return strength.sigma0


@pytest.mark.parametrize("gas", GASES)
def test_the_stored_records_are_fresh(gas):
    stored = (two_photon_bound._DIRECTORY / f"{gas}.json").read_text(encoding="utf-8")
    assert two_photon_bound._dumps(derivation.derive_bound(gas)) == stored


def test_the_convention_constant():
    # 2 pi^3 alpha^2 a0^4 (TC-952): 2.58949e-44 m^4 per delta at omega = 1 au; from
    # scipy's tabulated CODATA alpha and a0 independently, to their rounding raised to
    # the sixth power (9e-12 relative apart; 1e-10).
    alpha = constants.fine_structure
    a0 = constants.physical_constants["Bohr radius"][0]
    assert (
        pytest.approx(2 * math.pi**3 * alpha**2 * a0**4, rel=1e-10, abs=0.0)
        == SIGMA0_PER_DELTA
    )
    assert sigma0_from_delta(1.0, AU_ANGULAR_FREQUENCY) == pytest.approx(
        2.58949e-44, rel=1e-5, abs=0.0
    )


@pytest.mark.parametrize("gas", GASES)
def test_the_levels_are_nists_even_levels_below_threshold(gas):
    record = two_photon_bound_levels(gas)
    assert record.levels
    for level in record.levels:
        assert "*" not in level.term
        assert level.j in {0, 1, 2}
        assert 0 < level.energy < record.ionisation
        assert level.reach == ("other" if level.j == 1 else "linear")
        if gas == "He":
            assert level.term.startswith("1")
    energies = [level.energy for level in record.levels]
    assert energies == sorted(energies)


def test_xenons_6p_level_and_its_wavelength():
    # 80118.962 cm^-1 (NIST): two photons of 249.63 nm.
    level = two_photon_bound_levels("Xe").level(XE_6P, "2[1/2]", 0)
    assert level.printed == "80118.962"
    assert level.wavelength * 1e9 == pytest.approx(249.63, abs=0.005)


def test_xenon_calculations_against_measurements():
    # TC-952: Henriksson et al.'s TDDHF per level over Kroll and Bischel's G2-corrected
    # sigma_tot and Drag et al.'s; read to the derivation's two printed digits.
    xe = two_photon_bound_levels("Xe")
    ratios = []
    for config, term, j, measured in (
        (XE_6P, "2[1/2]", 0, "kroll_1990/sigma_tot"),
        (XE_6P, "2[3/2]", 2, "kroll_1990/sigma_tot"),
        (XE_6P, "2[5/2]", 2, "kroll_1990/sigma_tot"),
        (XE_6P1, "2[3/2]", 2, "drag_2021/sigma_6p1_32_2"),
        (XE_6P1, "2[1/2]", 0, "drag_2021/sigma_6p1_12_0"),
    ):
        level = xe.level(config, term, j)
        calc = _sigma0(xe, "henriksson_2005/tddhf_tpa_Xe", level)
        ratios.append(calc / _sigma0(xe, measured, level))
    assert [round(r, 2) for r in ratios] == [1.30, 0.83, 0.98, 1.17, 1.19]


def test_the_per_level_convention_matches_one_experiments_ratio():
    # Drag et al.'s 6p' J = 0 over J = 2, 1.38 (+0.21 -0.18), against Henriksson et
    # al.'s per level, 1.41; one component would give 7.0 (TC-952).
    xe = two_photon_bound_levels("Xe")
    j0, j2 = xe.level(XE_6P1, "2[1/2]", 0), xe.level(XE_6P1, "2[3/2]", 2)
    measured = _sigma0(xe, "drag_2021/sigma_6p1_12_0", j0) / _sigma0(
        xe, "drag_2021/sigma_6p1_32_2", j2
    )
    per_level = _sigma0(xe, "henriksson_2005/tddhf_tpa_Xe", j0) / _sigma0(
        xe, "henriksson_2005/tddhf_tpa_Xe", j2
    )
    assert measured == pytest.approx(1.38, abs=0.01)
    assert per_level == pytest.approx(1.41, abs=0.01)
    assert 1.38 - 0.18 <= per_level <= 1.38 + 0.21


def test_heliums_1s2s_strengths():
    he = two_photon_bound_levels("He")
    level = he.level("1s.2s", "1S", 0)
    deltas = {s.dataset: s.delta for s in he.strengths_of(level) if s.delta}
    assert deltas["hattig_1998b/S_He_CCSD_degenerate"] == pytest.approx(12.306)
    spread = max(deltas.values()) / min(deltas.values())
    assert spread < 1.19  # CCS 10.824 to CCSD 12.306, TDDHF 10.424
    assert he.assessment == "checked"


def test_kryptons_measurements_disagree():
    # Shekhtman et al.'s calculation is 94 times Khambatta et al.'s measurement at
    # 6p[3/2]2 (TC-952): the record holds both and calls krypton unresolved.
    kr = two_photon_bound_levels("Kr")
    level = kr.level("4s2.4p5.(2P*<3/2>).6p", "2[3/2]", 2)
    ratio = _sigma0(kr, "shekhtman_2020/qdos_basis3", level) / _sigma0(
        kr, "khambatta_1988/sigma0_6p32_2", level
    )
    assert ratio == pytest.approx(94, rel=0.01, abs=0.0)
    assert kr.assessment == "unresolved"
    kinds = {s.kind for s in kr.strengths_of(level)}
    assert kinds == {"calculated", "measured"}


def test_unassigned_values_keep_their_reasons():
    for gas in GASES:
        record = two_photon_bound_levels(gas)
        for item in (*record.strengths, *record.photoionisation):
            assert bool(item.level) != bool(item.unassigned), (gas, item.dataset)
        for s in record.strengths:
            assert (s.sigma0 is None) == (not s.level)


def test_photoionisation_in_si():
    xe = two_photon_bound_levels("Xe")
    level = xe.level(XE_6P, "2[1/2]", 0)
    (kroll,) = [
        p for p in xe.photoionisation_of(level) if p.dataset.startswith("kroll")
    ]
    # Kroll and Bischel's 4.3e-18 cm^2 at 248 nm, 1.9e-18 at 193 nm.
    assert kroll.sigma == pytest.approx((4.3e-22, 1.9e-22))
    omega = kroll.omega[0]
    assert omega is not None
    assert 2 * math.pi * constants.c / omega == pytest.approx(248e-9)
    he = two_photon_bound_levels("He")
    (d3,) = he.photoionisation_of(he.level("1s.3d", "1D", 2))
    assert d3.omega == (None,)  # Hummer and Storey state no 1D energy (TC-951)


def test_refusals_and_lookups():
    with pytest.raises(UnavailableModelError, match="noble gases"):
        two_photon_bound_levels("N2")
    xe = two_photon_bound_levels("Xe")
    with pytest.raises(DataError, match="no listed level"):
        xe.level(XE_6P, "2[9/2]", 4)
    assert xe.digest == xe.identity.digest
    assert xe.identity.name == "two_photon_bound/Xe"
    assert {"kroll_1990", "nist_asd_2026"} <= set(xe.identity.sources)


def test_the_record_checks_its_levels():
    xe = two_photon_bound_levels("Xe")
    with pytest.raises(DataError, match="unknown levels"):
        dataclasses.replace(xe, levels=xe.levels[:1])
    with pytest.raises(DataError, match="assessment"):
        dataclasses.replace(xe, assessment="good")


@pytest.fixture
def stored(monkeypatch, tmp_path):
    """Point the loader at an empty directory, and forget the cached records."""
    monkeypatch.setattr(two_photon_bound, "_DIRECTORY", tmp_path)
    two_photon_bound._stored.cache_clear()
    yield tmp_path
    two_photon_bound._stored.cache_clear()


def test_stored_record_errors(stored):
    with pytest.raises(DataError, match="missing"):
        two_photon_bound_levels("Ar")
    two_photon_bound._stored.cache_clear()
    (stored / "Ar.json").write_text(
        two_photon_bound._dumps(derivation.derive_bound("Ne")), encoding="utf-8"
    )
    with pytest.raises(DataError, match="holds the record of Ne"):
        two_photon_bound_levels("Ar")
    two_photon_bound._stored.cache_clear()
    raw = json.loads(two_photon_bound._dumps(derivation.derive_bound("Ar")))
    (stored / "Ar.json").write_text(
        json.dumps(raw | {"derivation_version": 0}), encoding="utf-8"
    )
    with pytest.raises(DataError, match="re-derive"):
        two_photon_bound_levels("Ar")
    two_photon_bound._stored.cache_clear()
    del raw["levels"]
    (stored / "Ar.json").write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(DataError, match="needs exactly the fields"):
        two_photon_bound_levels("Ar")


def _with(monkeypatch, gas: str, **entry):
    selection = derivation._selection()
    changed = selection | {gas: dict(selection[gas]) | entry}
    monkeypatch.setattr(derivation, "_selection", lambda: changed)


def test_one_component_strengths_are_scaled_to_the_level(monkeypatch):
    # A J = 2 level's single tensor component times 4/3 is the per-level strength
    # (TC-952); here applied to a printed value to test the arithmetic only.
    entry = {
        "dataset": "hattig_1998b/S_He_CCSD_degenerate",
        "row": "CCSD",
        "level": "1s.2s|1S|0",
        "convention": "one_component",
        "reason": "test",
    }
    _with(monkeypatch, "He", strengths=[entry], photoionisation=[])
    (strength,) = derivation.derive_bound("He").strengths
    assert strength.delta == pytest.approx(12.306 * 4 / 3)
    assert strength.derived


@pytest.mark.parametrize(
    ("kind", "entry", "match"),
    [
        (
            "strengths",
            {"dataset": "nobody_2000/x", "row": "1", "convention": "per_level"},
            "unknown dataset",
        ),
        (
            "strengths",
            {
                "dataset": "kroll_1990/sigma_tot",
                "row": "2p5",
                "level": "nowhere|x|0",
                "convention": "sigma0",
            },
            "no listed level",
        ),
        (
            "strengths",
            {"dataset": "kroll_1990/sigma_tot", "row": "2p7", "convention": "sigma0"},
            "no row",
        ),
        (
            "strengths",
            {
                "dataset": "kroll_1990/sigma_tot",
                "row": "2p5",
                "convention": "per_level",
            },
            "convention",
        ),
        (
            "photoionisation",
            {"dataset": "kroll_1990/sigma_tot", "level": "nowhere|x|0"},
            "no listed level",
        ),
        (
            "photoionisation",
            {"dataset": "kroll_1990/sigma_tot"},
            "not an excited level",
        ),
    ],
)
def test_the_derivation_refuses_what_it_cannot_use(monkeypatch, kind, entry, match):
    other = "photoionisation" if kind == "strengths" else "strengths"
    _with(monkeypatch, "Xe", **{kind: [entry], other: []})
    with pytest.raises(DataError, match=match):
        derivation.derive_bound("Xe")
