# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import pytest

from optmat.datasets import datasets
from optmat.gas.gases import Gas
from optmat.gas.virial import virial_ratio, virial_sources


def test_every_choice_points_at_primary_data_of_its_gas():
    by_key = {d.key: d for d in datasets()}
    for gas, (b_key, a_key) in virial_sources().items():
        assert by_key[b_key].material == gas
        assert by_key[a_key].material == gas
        assert by_key[b_key].quantity == "refractivity_virial_coefficient"
        assert by_key[a_key].quantity == "molar_refractivity"
        assert by_key[b_key].use == "primary"


def test_ratio_is_b_over_a_in_si():
    # Argon: B_R = 1.76 cm^6/mol^2, A_R = 4.196 cm^3/mol (Achtermann et al. 1991).
    assert virial_ratio("Ar") == pytest.approx(1.76e-12 / 4.196e-6, rel=1e-12, abs=0.0)


def test_gases_without_a_curated_value_have_no_pair_term():
    # O2's b is a published formula, not a pair of datasets (TC-794).
    missing = set(Gas) - set(virial_sources()) - {Gas.O2}
    assert Gas.D2 in missing
    for gas in missing:
        assert virial_ratio(gas) == 0.0


def test_helium_pair_term_is_negative_and_small():
    # Helium's B_R is negative; at 1 atm (4.46e-5 mol/cm^3) the correction is about
    # -6e-6 of n - 1.
    ratio = virial_ratio("He")
    assert -2e-7 < ratio < -1e-7


def test_the_like_terms_keep_their_uncertainty():
    from optmat.gas.virial import like_term  # noqa: PLC0415

    argon = like_term("Ar")
    # Achtermann et al. 1991: B_R = 1.76 +- 0.05 (3 sigma), A_R = 4.196 cm^3/mol.
    assert argon.sigma == pytest.approx(0.05 / 3 * 1e-12 / 4.196e-6, rel=1e-3, abs=0.0)
    assert (argon.status, argon.frequency, argon.temperature) == (
        "curated",
        "632.99",
        323.0,
    )
    krypton = like_term("Kr")
    assert (krypton.frequency, krypton.contribution) == ("632.99", "induced")
    deuterium = like_term("D2")
    assert (deuterium.status, deuterium.b, deuterium.contribution) == (
        "unknown",
        0.0,
        "",
    )
    # O2: Egan et al.'s b_R at 632.99 nm and 303 K, 0.194 + 0.00641 sigma^2 = 0.2100
    # cm^3/mol, with 50 % (TC-794).
    oxygen = like_term("O2")
    assert oxygen.status == "provisional"
    assert oxygen.b == pytest.approx(
        (0.194 + 0.00641 / 0.63299**2) * 1e-6, rel=1e-12, abs=0.0
    )
    assert oxygen.sigma == pytest.approx(0.5 * oxygen.b, rel=1e-12, abs=0.0)


def test_a_static_value_must_state_its_induced_part(monkeypatch):
    from optmat.exceptions import DataError  # noqa: PLC0415
    from optmat.gas import virial  # noqa: PLC0415

    # Krypton pointed at Orcutt and Cole's static B_eps, without its contribution.
    record = {k: dict(v) for k, v in virial._record().items()}
    record["Kr"] |= {
        "B_R": "orcutt_1967/B_this_work_Kr",
        "A_R": "orcutt_1967/A_this_work_Kr",
    }
    monkeypatch.setattr(virial, "_record", lambda: record)
    monkeypatch.setattr(
        virial,
        "virial_sources",
        lambda: {Gas.Kr: (record["Kr"]["B_R"], record["Kr"]["A_R"])},
    )
    virial.like_term.cache_clear()
    try:
        with pytest.raises(DataError, match="induced part is stated"):
            virial.like_term("Kr")
    finally:
        virial.like_term.cache_clear()


def test_the_third_terms():
    from optmat.gas.virial import pair_rules, third_term  # noqa: PLC0415

    # Achtermann et al. 1991, Table II: argon's C_R = -86.6(20) cm^9/mol^3 at 323 K.
    argon = third_term("Ar")
    assert argon.ratio == pytest.approx(86.6e-18 / 4.196e-6, rel=1e-3, abs=0.0)
    # Helium: the measured bound |C_R| <= 1 cm^9/mol^3.
    assert third_term("He").ratio == pytest.approx(1e-18 / 0.5213e-6, rel=1e-3, abs=0.0)
    # Neon, without a measured C_R: the declared scale on its B_R.
    neon = third_term("Ne")
    # Achtermann et al. 1993: B_R = -0.11(2), A_R = 1.0012 at 303 K.
    b_r, a_r = -0.11e-12, 1.0012e-6
    scale = float(pair_rules()["third_ratio"])
    assert neon.ratio == pytest.approx(scale * b_r**2 / a_r**2, rel=1e-2, abs=0.0)
    assert third_term("CO").status == "none"


def test_a_formula_must_be_curated(monkeypatch):
    from optmat.exceptions import DataError  # noqa: PLC0415
    from optmat.gas import virial  # noqa: PLC0415

    record = {k: dict(v) for k, v in virial._record().items()}
    record["O2"]["b_formula"] = "egan_2026/no_such_formula"
    monkeypatch.setattr(virial, "_record", lambda: record)
    virial._formula_b.cache_clear()
    try:
        with pytest.raises(DataError, match="no formula"):
            virial._formula_b("O2")
    finally:
        virial._formula_b.cache_clear()
