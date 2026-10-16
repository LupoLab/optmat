# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The optics of gas mixtures at finite density (OEP 4, Phase 1; its section 8).

One species is the pure law exactly, value and derivative; the unlike pairs are the mean
of the like terms with the declared uncertainty; the optical domain is the lower of the
third-term and pair-term bounds; the interface refuses what section 6 refuses.
"""

import functools
import itertools
import math
import warnings

import numpy as np
import pytest

from optmat.conditions import State
from optmat.constants import EPS0, N_A, C
from optmat.exceptions import ModelAbsorptionWarning, ModelRangeWarning, OptmatError
from optmat.gas.gases import Gas
from optmat.gas.mixture import MixtureOptics, mixture_optics
from optmat.gas.models.registry import default_model
from optmat.gas.virial import like_term, pair_rules, third_term

STATE = State(pressure=20e5, temperature=294.0)
"""20 bar: the pair terms are 1e-4 to 1e-3 of y, far above rounding."""


def _inside(model):
    """Return three wavelengths inside a model's qualified intervals, m."""
    low, high = model.intervals[0]
    return np.geomspace(low * 1.5, min(high, 3e-6) / 1.5, 3)


def _density(gas):
    """Return a pure gas's number density at STATE, as its model takes it."""
    return default_model(gas)._density(STATE, ideal=False)


@pytest.mark.parametrize("gas", list(Gas))
def test_one_species_is_the_pure_law(gas):
    model = default_model(gas)
    optics = mixture_optics([gas])
    wl = _inside(model)
    densities = {gas: _density(gas)}
    with warnings.catch_warnings():
        # A pure unknown-B_R gas is beyond the mixture domain at 20 bar (400 mol/m^3);
        # O3 absorbs in the ultraviolet, and both its laws say so.
        warnings.simplefilter("ignore", ModelRangeWarning)
        warnings.simplefilter("ignore", ModelAbsorptionWarning)
        got = optics.refractivity(wl, densities)
        index = np.sqrt(1 + optics.susceptibility(wl, densities))
        derivative = optics.dsusceptibility_domega(wl, densities)
        pure = model.refractivity(wl, STATE)
    # The same terms in a different order: agreement to a few ulps of n - 1.
    np.testing.assert_allclose(got, pure, rtol=1e-13)
    np.testing.assert_allclose(index, model.refractive_index(wl, STATE), rtol=1e-15)
    # d chi/d omega = (d chi/d y)(d y/d alpha)(d alpha/d omega), from the pure law.
    n, b = densities[gas], model.refractivity_virial
    alpha = model.complex_polarizability(wl)
    y = n * alpha / (3 * EPS0) * (1 + b * n / N_A)
    chain = 3 / (1 - y) ** 2 * n / (3 * EPS0) * (1 + b * n / N_A)
    np.testing.assert_allclose(
        derivative, chain * model.dpolarizability_domega(wl), rtol=1e-13
    )


def test_an_absorbing_constituent_keeps_its_imaginary_part():
    # O3 at 255 nm (the Hartley band): the complex law gives absorption, the real none.
    model = default_model("O3")
    optics = mixture_optics(["O3"])
    state = State(pressure=1e3, temperature=294.0)
    density = model._density(state, ideal=False)
    chi = optics.susceptibility(255e-9, {"O3": density})
    assert np.imag(chi) > 0
    np.testing.assert_allclose(
        np.sqrt(1 + chi), model.refractive_index(255e-9, state), rtol=1e-15
    )
    # The real convention leaves the absorption out, and warns as the pure model does.
    with pytest.warns(ModelAbsorptionWarning):
        real = optics.refractivity(255e-9, {"O3": density})
    with pytest.warns(ModelAbsorptionWarning):
        # The same refractivity by two routes: equal to the last bit on macOS, one
        # ulp apart on an AVX-512 Xeon (2026-10-16).
        assert real == pytest.approx(
            model.refractivity(255e-9, state), rel=1e-14, abs=0.0
        )
    # The conventions differ at second order in the imaginary part. With y = y' + i y'',
    # Re chi = chi(y') - 3 y''^2 and Re sqrt(1 + chi) - 1 = Re chi / 2 - Re(chi^2) / 8,
    # so the complex law falls below the real one by 3 y''^2 / 8 = (Im chi)^2 / 24, to a
    # relative 1e-5 (the next order, chi).
    difference = float(np.real(np.sqrt(1 + chi))) - 1 - real
    assert difference == pytest.approx(-(float(np.imag(chi)) ** 2) / 24, rel=1e-4)


def test_the_review_figures_for_ozone():
    # OEP 4 review R8: O3 at 255 nm, 1 bar and 294 K, printed to 8 and 5 figures.
    model = default_model("O3")
    density = model._density(State(pressure=1e5, temperature=294.0), ideal=False)
    optics = mixture_optics(["O3"])
    index = np.sqrt(1 + optics.susceptibility(255e-9, {"O3": density}))
    with pytest.warns(ModelAbsorptionWarning):
        real = float(optics.refractivity(255e-9, {"O3": density}))
    assert f"{real:.7e}" == "4.5393908e-04"
    assert f"{float(np.real(index)) - 1:.7e}" == "4.5388635e-04"
    assert f"{float(np.real(index)) - 1 - real:.4e}" == "-5.2724e-08"
    assert f"{float(np.imag(index)):.7e}" == "5.6189256e-04"


def test_the_pairs_follow_the_rule():
    optics = mixture_optics(["N2", "CO", "Ar", "He"])
    assert len(optics.pairs) == 10
    like = {p.species: p for p in optics.pairs if p.species == p.partner}
    assert {s: p.status for s, p in like.items()} == {
        "Ar": "curated",
        "He": "curated",
        "N2": "curated",
        "CO": "unknown",
    }
    assert like["CO"].b == 0.0
    assert like["CO"].sigma == float(pair_rules()["unknown_b"]) * 1e-6
    for p in optics.pairs:
        if p.species != p.partner:
            unknown = "CO" in (p.species, p.partner)
            assert p.status == ("unknown" if unknown else "estimate")
            assert p.source == "rule: mean"
            assert math.isnan(p.b)
        assert p.species <= p.partner
    assert like["Ar"].b == like_term("Ar").b
    assert optics.status == "provisional"
    assert "CO" in optics.status_reason
    assert mixture_optics(["He", "Ar"]).status == "validated"
    # O2's like term is Egan et al.'s published estimate: its pairs are estimates,
    # and its mixtures provisional (TC-794).
    air = mixture_optics(["N2", "O2", "Ar"])
    statuses = {(p.species, p.partner): p.status for p in air.pairs}
    assert statuses[("O2", "O2")] == "provisional"
    assert statuses[("N2", "O2")] == "estimate"
    assert air.status == "provisional"
    assert "O2 a provisional estimate" in air.status_reason


def test_an_unlike_term_is_the_mean_of_the_like_terms():
    # For 50/50 He/Ar, y minus the pure-gas parts is 2 B_st x_He x_Ar rho^2, with
    # B_st = (b_He A_He + b_Ar A_Ar) / 2.
    optics = mixture_optics(["He", "Ar"])
    n = 2.5e26
    wl = 800e-9
    chi = optics.susceptibility(wl, {"He": n, "Ar": n})
    y = chi / (chi + 3)
    rho = n / N_A
    a = {
        g: N_A * default_model(g).complex_polarizability(wl) / (3 * EPS0)
        for g in ("He", "Ar")
    }
    b = {g: like_term(g).b for g in ("He", "Ar")}
    pure = sum(a[g] * rho * (1 + b[g] * rho) for g in a)
    cross = 2 * 0.5 * (b["He"] * a["He"] + b["Ar"] * a["Ar"]) * rho**2
    # y is rebuilt from chi: rounding of y, about 1e-16 of y, against 5e-3 of y here.
    np.testing.assert_allclose(y - pure, cross, rtol=1e-9)


def test_the_uncertainty_is_the_larger_like_term():
    optics = mixture_optics(["He", "Ar"])
    n = 2.5e26
    wl = 800e-9
    got = optics.refractivity_uncertainty(wl, {"He": n, "Ar": n})
    rho = n / N_A
    a = {
        g: N_A * default_model(g).polarizability(wl) / (3 * EPS0) for g in ("He", "Ar")
    }
    sigma = max(abs(like_term(g).b) * a[g] for g in a)
    # d(n - 1)/d y at this y, by a difference of the law itself.
    dn = optics.refractivity(wl, {"He": n, "Ar": n})
    y = ((1 + dn) ** 2 - 1) / ((1 + dn) ** 2 + 2)
    # A second-order difference: truncation (step/y)^2, 1e-8; rounding 1e-16 / step,
    # 1e-9 of the slope.
    step = 1e-4 * y

    def index(value):
        return math.sqrt((1 + 2 * value) / (1 - value)) - 1

    slope = (index(y + step) - index(y - step)) / (2 * step)
    assert got == pytest.approx(slope * 2 * sigma * rho**2, rel=1e-6)
    # A curated like term alone carries no sensitivity: the pure gas has none.
    assert mixture_optics(["Ar"]).refractivity_uncertainty(wl, {"Ar": n}) == 0.0


def _omega(wavelength):
    return 2 * math.pi * C / wavelength


@pytest.mark.parametrize("species", [("He", "Ar"), ("N2", "O2", "Ar"), ("CO2", "Xe")])
def test_the_frequency_derivative_against_differences(species):
    optics = mixture_optics(species)
    densities = dict.fromkeys(species, 3e26)
    wl = 900e-9
    omega = _omega(wl)
    h = 1e-3 * omega
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        f = [
            optics.susceptibility(2 * math.pi * C / (omega + k * h), densities)
            for k in (-2, -1, 1, 2)
        ]
        got = optics.dsusceptibility_domega(wl, densities)
    # Fourth-order central difference: truncation h^4 chi^(5)/30, about 1e-12 of the
    # derivative this far from resonance; rounding 1e-16 chi / h, below 1e-11.
    fd = (f[0] - 8 * f[1] + 8 * f[2] - f[3]) / (12 * h)
    np.testing.assert_allclose(got, fd, rtol=1e-8)


@pytest.mark.parametrize("species", [("He", "Ar"), ("N2", "O2", "Ar"), ("SF6", "He")])
def test_the_density_derivative_against_differences(species):
    optics = mixture_optics(species)
    densities = {s: (i + 1) * 2e26 for i, s in enumerate(species)}
    wl = 800e-9
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        got = optics.dsusceptibility_ddensity(wl, densities)
        for s in species:
            h = 1e-3 * densities[s]

            def chi(k, s=s, h=h):
                return optics.susceptibility(wl, densities | {s: densities[s] + k * h})

            # chi is a rational function of the densities with its pole at y = 1:
            # the fourth-order difference is exact to rounding, 1e-16 chi / h.
            fd = (chi(-2) - 8 * chi(-1) + 8 * chi(1) - chi(2)) / (12 * h)
            np.testing.assert_allclose(got[s], fd, rtol=1e-9)
    assert set(got) == set(species)


def test_the_densities_broadcast_with_the_wavelength():
    optics = mixture_optics(["He", "Ar"])
    wl = np.array([600e-9, 800e-9])[:, None]
    n = np.array([0.0, 1e25, 2e25])
    got = optics.refractivity(wl, {"He": n, "Ar": 1e25})
    assert got.shape == (2, 3)
    for i, j in itertools.product(range(2), range(3)):
        one = optics.refractivity(wl[i, 0], {"He": n[j], "Ar": 1e25})
        assert got[i, j] == one
    derivatives = optics.dsusceptibility_ddensity(wl, {"He": n, "Ar": 1e25})
    assert derivatives["He"].shape == (2, 3)


def test_zero_densities():
    optics = mixture_optics(["He", "Ar"])
    with warnings.catch_warnings():
        warnings.simplefilter("error", ModelRangeWarning)
        assert optics.susceptibility(800e-9, {"He": 0.0, "Ar": 0.0}) == 0
        alone = optics.refractivity(800e-9, {"He": 0.0, "Ar": 1e25})
    pure = mixture_optics(["Ar"]).refractivity(800e-9, {"Ar": 1e25})
    assert alone == pure
    # The intervals are the species', not the densities': zero still warns outside.
    with pytest.warns(ModelRangeWarning, match="common range"):
        optics.refractivity(50e-9, {"He": 0.0, "Ar": 0.0})


@pytest.mark.parametrize(
    ("densities", "match"),
    [
        ({"He": 1e25}, "missing \\['Ar'\\]"),
        ({"He": 1e25, "Ar": 1e25, "N2": 1e25}, "extra \\['N2'\\]"),
        ({"He": -1.0, "Ar": 1e25}, "must not be negative"),
    ],
)
def test_the_keys_are_checked(densities, match):
    with pytest.raises(OptmatError, match=match):
        mixture_optics(["He", "Ar"]).refractivity(800e-9, densities)


def test_the_species_are_checked():
    with pytest.raises(OptmatError, match="distinct species"):
        mixture_optics([])
    with pytest.raises(OptmatError, match="distinct species"):
        mixture_optics(["Ar", Gas.Ar])
    with pytest.raises(OptmatError, match="unknown gas"):
        mixture_optics(["Ar", "argon"])


def test_order_and_symmetry():
    one = mixture_optics(["N2", "Ar", "He"])
    assert one is mixture_optics([Gas.He, "N2", Gas.Ar])
    densities = {"N2": 1e26, "Ar": 2e26, "He": 3e26}
    reordered = {"He": 3e26, "N2": 1e26, "Ar": 2e26}
    assert one.refractivity(700e-9, densities) == one.refractivity(700e-9, reordered)
    assert one.identity.digest != mixture_optics(["N2", "Ar"]).identity.digest


def test_the_dilute_limit():
    # At 1e20 m^-3 the pair and local-field terms are below 1e-8 of chi.
    optics = mixture_optics(["N2", "O2", "Ar"])
    densities = {"N2": 7.8e19, "O2": 2.1e19, "Ar": 1e18}
    got = optics.susceptibility(633e-9, densities)
    dilute = sum(
        n * default_model(g).complex_polarizability(633e-9) / EPS0
        for g, n in densities.items()
    )
    np.testing.assert_allclose(got, dilute, rtol=1e-8)


def test_the_third_term_bound():
    fraction = float(pair_rules()["density_fraction"])
    argon = mixture_optics(["Ar"])
    assert third_term("Ar").status == "measured"
    expected = math.sqrt(fraction / third_term("Ar").ratio) * N_A
    assert argon.density_bound({"Ar": 1.0}) == pytest.approx(expected, rel=1e-15)
    # OEP 4 section 4.3: about 7000 and 6700 mol/m^3.
    assert argon.density_bound({"Ar": 1.0}) / N_A == pytest.approx(6961, abs=1)
    assert mixture_optics(["N2"]).density_bound({"N2": 1.0}) / N_A == pytest.approx(
        6675, abs=1
    )
    assert third_term("He").status == "bound"
    assert third_term("Ne").status == "scale"
    assert third_term("CO").status == "none"
    assert math.isnan(third_term("CO").ratio)
    # O2's b, an estimate: the declared scale on it.
    assert third_term("O2").status == "scale"


def test_the_pair_bound():
    fraction = float(pair_rules()["density_fraction"])
    # An unknown like term alone: b rho reaches 1e-3 at 1e-3 / 2 cm^3/mol.
    deuterium = mixture_optics(["D2"])
    scale = float(pair_rules()["unknown_b"]) * 1e-6
    assert deuterium.density_bound({"D2": 1.0}) / N_A == pytest.approx(
        fraction / scale, rel=1e-12
    )
    # 50/50 He/Ar: OEP 4 section 7.2 puts it near 2700 mol/m^3 at one wavelength;
    # the smallest over the qualified interval is lower.
    he_ar = mixture_optics(["He", "Ar"])
    bound = he_ar.density_bound({"He": 0.5, "Ar": 0.5}) / N_A
    assert 2300 < bound < 2750
    # An absent species does not bound it.
    assert he_ar.density_bound({"Ar": 1.0}) == mixture_optics(["Ar"]).density_bound(
        {"Ar": 1.0}
    )


def test_the_domain_warns_beyond_its_bound():
    optics = mixture_optics(["He", "Ar"])
    bound = optics.density_bound({"He": 0.5, "Ar": 0.5})
    with warnings.catch_warnings():
        warnings.simplefilter("error", ModelRangeWarning)
        optics.refractivity(800e-9, {"He": 0.49 * bound, "Ar": 0.49 * bound})
    with pytest.warns(ModelRangeWarning, match="optical domain"):
        optics.refractivity(800e-9, {"He": 0.51 * bound, "Ar": 0.51 * bound})


@pytest.mark.parametrize(
    ("composition", "match"),
    [
        ({"He": 0.5, "N2": 0.5}, "not a composition"),
        ({"He": 1.5, "Ar": -0.5}, "not a composition"),
        ({"He": 0.5, "Ar": 0.4}, "not one"),
    ],
)
def test_the_composition_is_checked(composition, match):
    with pytest.raises(OptmatError, match=match):
        mixture_optics(["He", "Ar"]).density_bound(composition)


def test_the_intervals_are_the_intersection():
    optics = mixture_optics(["He", "CO2"])
    assert optics.intervals == default_model("CO2").intervals
    assert isinstance(optics, MixtureOptics)
    assert optics.identity.kind == "gas_refractive_mixture"
    assert optics.identity.name == "mixture/CO2+He"


def _beta2(function, wavelength):
    """Return the GVD of n - 1 = function(wavelength), fs^2/m, by differences."""
    omega = 2 * math.pi * C / wavelength
    h = 2e-3 * omega
    n = [float(function(2 * math.pi * C / (omega + k * h))) for k in range(-2, 3)]
    first = (n[0] - 8 * n[1] + 8 * n[3] - n[4]) / (12 * h)
    second = (-n[0] + 16 * n[1] - 30 * n[2] + 16 * n[3] - n[4]) / (12 * h * h)
    return (2 * first + omega * second) / C * 1e30


def test_the_pages_sensitivity_table():
    from pathlib import Path  # noqa: PLC0415

    from optmat.constants import K_B  # noqa: PLC0415
    from optmat.gas.composition import composition  # noqa: PLC0415

    page = (Path(__file__).parents[1] / "docs" / "mixtures.md").read_text(
        encoding="utf-8"
    )
    mixtures = {
        "He/Ar": {"He": 0.5, "Ar": 0.5},
        "Ne/Ar": {"Ne": 0.5, "Ar": 0.5},
        "He/Ne": {"He": 0.5, "Ne": 0.5},
        "Air": {
            c.material: c.mole_fraction for c in composition("dry air (N2, O2, Ar)")
        },
    }
    rows = [
        line.split("|")[1:-1] for line in page.splitlines() if line.startswith("| ")
    ]
    rows = [[c.strip() for c in r] for r in rows if r[0].strip() in mixtures]
    assert len(rows) == 8
    for name, pressure, dn, spread, relative, gvd, gvd_spread in rows:
        x = mixtures[name]
        optics = mixture_optics(list(x))
        total = float(pressure.split()[0]) * 1e5 / (K_B * 294.0)
        densities = {s: total * v for s, v in x.items()}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ModelRangeWarning)
            got = float(optics.refractivity(800e-9, densities))
            u = float(optics.refractivity_uncertainty(800e-9, densities))
            b2 = _beta2(
                functools.partial(optics.refractivity, densities=densities), 800e-9
            )
            ub2 = _beta2(
                functools.partial(optics.refractivity_uncertainty, densities=densities),
                800e-9,
            )
        # Each to the figures printed.
        assert got == pytest.approx(float(dn), abs=5e-4 * float(dn))
        assert u == pytest.approx(float(spread), rel=0.05)
        assert u / got == pytest.approx(float(relative), rel=0.05)
        assert b2 == pytest.approx(float(gvd), abs=0.05)
        assert ub2 == pytest.approx(float(gvd_spread), rel=0.05)


def _variant(kind, wavelength, densities, optics):
    """Return n - 1 by the full law, without pair terms, or the dilute law."""
    if kind == "full":
        return float(optics.refractivity(wavelength, densities))
    alpha = {s: default_model(s).polarizability(wavelength) for s in densities}
    if kind == "dilute":
        chi = sum(densities[s] * alpha[s] / EPS0 for s in densities)
    else:
        y = sum(densities[s] * alpha[s] / (3 * EPS0) for s in densities)
        chi = 3 * y / (1 - y)
    return float(chi / (1 + np.sqrt(1 + chi)))


def _variant_gvd(kind, densities, optics):
    """Return the GVD (fs^2/m) at 800 nm of a variant, by a centred difference."""
    omega = 2 * math.pi * C / 800e-9
    step = omega * 1e-3
    beta = [
        (omega + j * step)
        / C
        * _variant(kind, 2 * math.pi * C / (omega + j * step), densities, optics)
        for j in (-1, 0, 1)
    ]
    return (beta[0] - 2 * beta[1] + beta[2]) / step**2 * 1e30


def test_the_pages_he_ar_sensitivity_table():
    from pathlib import Path  # noqa: PLC0415

    from optmat.gas.eos import equation_of_state  # noqa: PLC0415
    from optmat.gas.mixture_eos import mixture_equation_of_state  # noqa: PLC0415

    page = (Path(__file__).parents[1] / "docs" / "mixtures.md").read_text(
        encoding="utf-8"
    )
    rows = [
        [c.strip() for c in line.split("|")[1:-1]]
        for line in page.splitlines()
        if line.startswith("| ") and " bar | " in line and "e-3 /" in line
    ]
    assert len(rows) == 3
    x = {"He": 0.5, "Ar": 0.5}
    optics = mixture_optics(list(x))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        default = mixture_equation_of_state(x)
        multi = mixture_equation_of_state(x, method="multi-fluid")
        for pressure, dn, dilute, pairs, pp, mf, estimate in rows:
            p = float(pressure.split()[0]) * 1e5
            state = State(pressure=p, temperature=294.0)
            full_d = default.partial_densities(state)
            full = _variant("full", 800e-9, full_d, optics)
            pp_d = {
                s: equation_of_state(s).number_density(
                    State(pressure=v * p, temperature=294.0)
                )
                for s, v in x.items()
            }
            got = {
                "dilute": _variant("dilute", 800e-9, full_d, optics) / full - 1,
                "pairs": _variant("pairs", 800e-9, full_d, optics) / full - 1,
                "pp": _variant("full", 800e-9, pp_d, optics) / full - 1,
                "mf": _variant("full", 800e-9, multi.partial_densities(state), optics)
                / full
                - 1,
            }
            assert full == pytest.approx(float(dn), rel=5e-4)
            full_gvd = _variant_gvd("full", full_d, optics)
            gvd = {
                "dilute": _variant_gvd("dilute", full_d, optics) / full_gvd - 1,
                "pairs": _variant_gvd("pairs", full_d, optics) / full_gvd - 1,
                "pp": _variant_gvd("full", pp_d, optics) / full_gvd - 1,
                "mf": _variant_gvd("full", multi.partial_densities(state), optics)
                / full_gvd
                - 1,
            }
            # Each cell's n - 1 and GVD, to the two figures printed.
            for key, cell in (
                ("dilute", dilute),
                ("pairs", pairs),
                ("pp", pp),
                ("mf", mf),
            ):
                n_cell, gvd_cell = (
                    float(c.replace("\N{MINUS SIGN}", "-")) for c in cell.split("/")
                )
                assert got[key] == pytest.approx(n_cell, rel=0.05), (pressure, key)
                assert gvd[key] == pytest.approx(gvd_cell, rel=0.05), (pressure, key)
            spread = float(optics.refractivity_uncertainty(800e-9, full_d)) / full
            assert spread == pytest.approx(float(estimate.strip("±")), rel=0.05)
