# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The dilute gases' transport properties and Wilke's rule (OEP 1, 4.3, 6.5 and 8).

Each viscosity must be its curated paper's: CoolProp's correlations reproduce the
printed dilute-gas forms, Huber's form reproduces his check values, and D2's
correlation its table, smoothly. The noble gases' heat capacity is (5/2) R. Where
solprop takes its values from the same sources (CoolProp at one bar, the WebBook tables
at 0.1 MPa), ours differ from them by the density term alone. Wilke's rule reduces to
the pure gas and does not depend on the order of the gases; temperatures outside a
range are refused.
"""

import csv
import dataclasses
import math

import numpy as np
import pytest

from optmat.constants import K_B, N_A
from optmat.datasets import DATA_DIRECTORY, sources
from optmat.exceptions import DataError, OptmatError, StateError
from optmat.gas import transport
from optmat.gas.air import DEFAULT_AIR
from optmat.gas.composition import composition
from optmat.gas.eos import _abstract_state, _coolprop
from optmat.gas.transport import (
    TRANSPORT_GASES,
    heat_capacity,
    mixture_viscosity,
    molar_mass,
    transport_model,
    viscosity,
)

R = N_A * K_B
WEBBOOK = DATA_DIRECTORY / "external" / "nist_webbook"


def _formula(source, formula_id):
    return transport._formula(source, formula_id)


def _lennard_jones(c, t):
    # Lemmon and Jacobsen's Eq. (2), Muzny et al.'s Eqs. (3)-(4): uPa s.
    ln_t = np.log(t / c["epsilon_k"])
    b = [c[k] for k in sorted(c) if len(k) == 2 and k[0] in "ab" and k[1].isdigit()]
    omega = np.exp(sum(bi * ln_t**i for i, bi in enumerate(b)))
    return c["prefactor"] * np.sqrt(c["molar_mass"] * t) / (c["sigma_nm"] ** 2 * omega)


def _arp(c, t):
    micro_g_per_cm_s = (
        c["a"] * t ** c["b"] * np.exp(c["c"] / t - c["d"] / t**2 - c["e"])
    )
    return micro_g_per_cm_s * 0.1


def _laesecke(c, t):
    cube = t ** (1 / 3)
    denominator = (
        c["a0"]
        + c["a1"] * t ** (1 / 6)
        + c["a2"] * np.exp(c["a3"] * cube)
        + (c["a4"] + c["a5"] * cube) / np.exp(cube)
        + c["a6"] * np.sqrt(t)
    )
    return c["scale"] * np.sqrt(t) / denominator * 1e3


def _quinones(c, t):
    tr = t / c["critical_temperature"]
    terms = [c[f"d{i}"] * tr ** (i / 4) for i in range(5) if f"d{i}" in c]
    return sum(terms) * 1e3


PRINTED = {
    "He": ("arp_1998", "eq8", _arp),
    "N2": ("lemmon_2004", "eq2_N2", _lennard_jones),
    "Ar": ("lemmon_2004", "eq2_Ar", _lennard_jones),
    "O2": ("lemmon_2004", "eq2_O2", _lennard_jones),
    "H2": ("muzny_2013", "eq3_eq4", _lennard_jones),
    "CO2": ("laesecke_2017", "eq4", _laesecke),
    "CH4": ("quinones_cisneros_2006", "eq24_CH4", _quinones),
    "SF6": ("quinones_cisneros_2012", "eq6", _quinones),
}


@pytest.mark.parametrize("gas", sorted(PRINTED))
def test_coolprop_reproduces_the_printed_dilute_form(gas):
    source, formula_id, form = PRINTED[gas]
    model = transport_model(gas)
    assert model.viscosity_source == source
    low, high = model.temperature_range("viscosity")
    t = np.geomspace(max(low, 120.0), min(high, 900.0), 9)
    printed = form(_formula(source, formula_id), t) * 1e-6
    # The same formula in double precision: agreement to rounding, below 1e-9.
    np.testing.assert_allclose(model.viscosity(t), printed, rtol=1e-9)


@pytest.mark.parametrize(
    ("gas", "temperature", "check"),
    # Huber, Table 7 (pdf p. 274): the dilute viscosities eta*, uPa s.
    [
        ("Ne", 40.0, 6.006075),
        ("Kr", 188.5, 16.32999),
        ("Xe", 260.8, 20.11099),
        ("N2O", 278.6, 13.86694),
    ],
)
def test_huber_form_reproduces_his_check_values(gas, temperature, check):
    parameters = transport_model(gas).viscosity_parameters
    value = transport._chapman_enskog(parameters, np.array(temperature))
    # Table 7 prints seven digits: 3e-7 covers its rounding.
    assert float(value) == pytest.approx(check * 1e-6, rel=3e-7)


def _webbook(gas):
    rows = []
    for name in (f"{gas}_0.1MPa_200-400K.tsv", f"{gas}_0.1MPa_293.15-298.15K.tsv"):
        with (WEBBOOK / name).open(newline="", encoding="utf-8") as handle:
            rows += list(csv.DictReader(handle, delimiter="\t"))
    t = np.array([float(r["Temperature (K)"]) for r in rows])
    eta = np.array([float(r["Viscosity (uPa*s)"]) for r in rows]) * 1e-6
    order = np.argsort(t)
    return t[order], eta[order]


@pytest.mark.parametrize(
    ("gas", "offset"),
    [("Ne", 1.3e-4), ("Kr", 1.34e-3), ("Xe", 2.11e-3), ("N2O", 2.05e-3)],
)
def test_the_webbook_tables_lie_above_by_the_density_term(gas, offset):
    t, table = _webbook(gas)
    excess = table / viscosity(gas, t) - 1
    # The tables are the same models at 0.1 MPa (solprop's source): positive and
    # falling with temperature (the density term), at most 0.54 % (Xe at 200 K), below
    # every stated uncertainty (0.7 % for Xe). Each 5 K step may rise by the tables'
    # rounding (five digits: 3e-5); at 298.15 K, the offsets of the record.
    stated = {"Ne": 0.02, "Kr": 0.03, "Xe": 0.007, "N2O": 0.05}[gas]
    assert np.all((excess > 0) & (excess < min(6e-3, stated)))
    assert np.all(np.diff(excess[t % 5 == 0]) < 3e-5)
    assert excess[t == 298.15][0] == pytest.approx(offset, abs=4e-5)


@pytest.mark.parametrize("gas", ["He", "Ar", "N2", "H2", "O2", "CH4", "SF6"])
def test_solprops_one_bar_values_lie_above_by_the_density_term(gas):
    # solprop takes these gases' viscosities from CoolProp at one bar over 200-600 K:
    # the same correlation, with its density terms. They move it by up to 1.7e-3, at
    # 200 K (N2; TC-762), below every stated uncertainty (0.5 % for N2 and Ar); N2 at
    # 294 K is 1.7612e-5 against 1.7599e-5 Pa s (OEP 1, section 4.3). SF6, from
    # solprop's 225 K, near its triple point and saturation at one bar, lies 4.7e-3
    # below, within its stated 1.5 % there.
    model = transport_model(gas)
    low, _ = model.temperature_range("viscosity")
    state = _abstract_state(model.viscosity_parameters["fluid"])
    bound = 5e-3 if gas == "SF6" else 2e-3
    for t in np.linspace(max(low, 225.0 if gas == "SF6" else 200.0), 600.0, 9):
        state.update(_coolprop().PT_INPUTS, 1e5, float(t))
        assert abs(state.viscosity() / float(model.viscosity(t)) - 1) < bound
    if gas == "N2":
        assert float(model.viscosity(294.0)) == pytest.approx(1.7599e-5, abs=5e-10)


def _deuterium_table():
    (dataset,) = [d for d in sources()["assael_1987"].datasets if d.id == "table2"]
    printed = [str(v) for v in dataset.values]
    t = np.asarray(dataset.abscissa_values(), dtype=float)
    table = np.asarray([float(v) for v in printed]) * dataset.scale
    unit = np.array([10.0 ** -len(v.split(".")[1]) for v in printed]) * dataset.scale
    return t, table, unit


def test_the_deuterium_correlation_reproduces_its_table():
    # TC-765: assael_1986's Eqs. (1) and (3) with the D2 mass, times assael_1987's
    # f(T*) at every temperature, give Table 2 within 1.5 units of its last digit;
    # the 90 K cell (5.410, degraded in the scan) is 6.5 units off. D2 is served
    # from the correlation (TC-839).
    model = transport_model("D2")
    assert model.viscosity_form == "assael"
    t, table, unit = _deuterium_table()
    units = np.abs(model.viscosity(t) - table) / unit
    assert np.all(units[t != 90] < 1.5)
    assert units[t == 90][0] == pytest.approx(6.5, abs=0.1)
    assert model.status == "validated"
    assert model.identity.sources == ("assael_1986", "assael_1987", "richardson_2014")


def test_the_deuterium_viscosity_is_smooth():
    # solprop SEP 3 (Q33) represents the viscosity by Chebyshev series over 60-600 K
    # and wider. The cubic spline through Table 2 stalled at 4e-4 to 8e-4 whatever the
    # degree (its knots' four-digit rounding, 2e-4, and the 90 K cell); the correlation
    # is analytic, and degree 20 in ln T reaches 1e-11 of ln eta, the rounding of the
    # fit (TC-839).
    model = transport_model("D2")
    for low, high in ((60.0, 600.0), (30.0, 2200.0)):
        x = np.linspace(-1.0, 1.0, 2001)
        t = np.clip(np.exp(np.log(low) + (x + 1) / 2 * np.log(high / low)), low, high)
        y = np.log(model.viscosity(t))
        series = np.polynomial.chebyshev.Chebyshev.fit(x, y, 20)
        assert np.max(np.abs(series(x) - y)) < 1e-9


def test_the_noble_gases_heat_capacity_is_five_halves_r():
    for gas in ("He", "Ne", "Ar", "Kr", "Xe"):
        model = transport_model(gas)
        low, high = model.temperature_range("heat_capacity")
        t = np.geomspace(low, high, 7)
        # The equation's c_p0/R is 5/2 to rounding; optmat's R is CODATA's.
        np.testing.assert_allclose(model.heat_capacity(t), 2.5 * R, rtol=1e-14)


def test_sulfur_hexafluorides_heat_capacity_is_guder_and_wagners():
    c = _formula("guder_2009", "eq15")
    t = np.array([230.0, 294.0, 400.0, 625.0])
    tau = c["critical_temperature"] / t
    reduced = (
        1
        + c["n3"]
        + sum(
            c[f"n{i}"]
            * (c[f"theta{i}"] * tau) ** 2
            * np.exp(c[f"theta{i}"] * tau)
            / np.expm1(c[f"theta{i}"] * tau) ** 2
            for i in (4, 5, 6)
        )
    )
    np.testing.assert_allclose(heat_capacity("SF6", t), reduced * R, rtol=1e-14)


def test_oxygens_heat_capacity_stops_where_it_leaves_its_source():
    # CoolProp's ideal part is a refit of Stewart et al.'s Eq. (11): within 3.3e-6 to
    # 1500 K; from 1800 K beyond the equation's 0.003 J/(mol K), so the range stops at
    # 1700 K (TC-765).
    c = _formula("stewart_1991", "eq11")
    t = np.array([60.0, 100.0, 294.0, 1000.0, 1500.0])
    u, eta = c["N6"] / t, c["N7"] / t
    reduced = (
        c["N1"] * t**-1.5
        + c["N2"]
        + c["N3"] * t**2
        + c["N4"] * u**2 * np.exp(u) / np.expm1(u) ** 2
        + c["N5"] * (2 / 3) * eta**2 * np.exp(-eta) / (1 + (2 / 3) * np.exp(-eta)) ** 2
    )
    np.testing.assert_allclose(heat_capacity("O2", t), reduced * R, rtol=3.5e-6)
    assert transport_model("O2").temperature_range("heat_capacity") == (54.361, 1700.0)
    with pytest.raises(StateError, match="heat capacity of O2"):
        heat_capacity("O2", 1800.0)


def test_nitrous_oxides_heat_capacity_is_lemmon_and_spans():
    c = _formula("lemmon_2006", "eq10_N2O")
    t = np.array([200.0, 294.0, 500.0])
    reduced = c["c0"] + sum(
        c[f"v{k}"]
        * (c[f"u{k}"] / t) ** 2
        * np.exp(c[f"u{k}"] / t)
        / np.expm1(c[f"u{k}"] / t) ** 2
        for k in (1, 2, 3)
    )
    np.testing.assert_allclose(heat_capacity("N2O", t), reduced * R, rtol=1e-14)


def test_the_molar_mass_is_the_equation_of_states():
    assert molar_mass("Ar") == 0.039948
    assert molar_mass("N2O") == pytest.approx(44.0128e-3, rel=1e-15)
    parameters = transport_model("Kr").viscosity_parameters
    assert parameters["molar_mass"] * 1e-3 == pytest.approx(molar_mass("Kr"))


class TestWilke:
    def test_it_reduces_to_the_pure_gas(self):
        for gas in ("He", "Xe", "N2"):
            assert float(mixture_viscosity({gas: 1.0}, 300.0)) == float(
                viscosity(gas, 300.0)
            )
            together = mixture_viscosity({gas: 0.5, "Ar": 0.5}, 300.0)
            assert float(mixture_viscosity({gas: 1.0, "Ar": 0.0}, 300.0)) == float(
                viscosity(gas, 300.0)
            )
            assert float(together) > 0

    def test_it_is_symmetric_under_relabelling(self):
        t = np.array([250.0, 300.0, 400.0])
        a = mixture_viscosity({"He": 0.3, "Ar": 0.5, "Xe": 0.2}, t)
        b = mixture_viscosity({"Xe": 0.2, "He": 0.3, "Ar": 0.5}, t)
        # The same terms summed in another order.
        np.testing.assert_allclose(a, b, rtol=1e-15)

    def test_the_binary_form(self):
        # Wilke's Eq. (11), the rule for two gases, written out.
        eta = {g: float(viscosity(g, 293.15)) for g in ("He", "Ar")}
        m = {g: molar_mass(g) for g in ("He", "Ar")}

        def phi(i, j):
            top = (1 + math.sqrt(eta[i] / eta[j]) * (m[j] / m[i]) ** 0.25) ** 2
            return top / math.sqrt(8 * (1 + m[i] / m[j]))

        # His identity phi_21 = phi_12 (eta_2/eta_1)(M_1/M_2) (wilke_1950's notes).
        assert phi("Ar", "He") == pytest.approx(
            phi("He", "Ar") * eta["Ar"] / eta["He"] * m["He"] / m["Ar"], rel=1e-14
        )
        x = 0.4
        expected = x * eta["He"] / (x + (1 - x) * phi("He", "Ar")) + (1 - x) * eta[
            "Ar"
        ] / (x * phi("Ar", "He") + 1 - x)
        assert float(mixture_viscosity({"He": x, "Ar": 1 - x}, 293.15)) == (
            pytest.approx(expected, rel=1e-14)
        )

    def test_the_named_compositions(self):
        air = {c.material: c.mole_fraction for c in composition(DEFAULT_AIR)}
        assert float(mixture_viscosity("dry air", 294.0)) == float(
            mixture_viscosity(air, 294.0)
        )
        # Lemmon and Jacobsen's pseudo-pure air (1 % in the dilute gas): Wilke's rule
        # over N2, O2 and Ar lies 9e-4 below it at 294 K, well inside both accuracies.
        state = _abstract_state("Air")
        state.update(_coolprop().DmolarT_INPUTS, 1e-6, 294.0)
        assert float(mixture_viscosity("dry air", 294.0)) == pytest.approx(
            state.viscosity(), rel=2e-3
        )
        with pytest.raises(OptmatError, match="no transport model for CO"):
            mixture_viscosity("dry air (CIPM-2007)", 294.0)

    def test_bad_fractions_are_refused(self):
        with pytest.raises(OptmatError, match="sum to one"):
            mixture_viscosity({"He": 0.5, "Ar": 0.4}, 300.0)
        with pytest.raises(OptmatError, match="not negative"):
            mixture_viscosity({"He": 1.5, "Ar": -0.5}, 300.0)
        with pytest.raises(OptmatError, match="not negative"):
            mixture_viscosity({"He": math.nan}, 300.0)


class TestRanges:
    def test_outside_a_range_is_refused(self):
        with pytest.raises(
            StateError, match="viscosity of He is defined for 110-1500 K"
        ):
            viscosity("He", 100.0)
        with pytest.raises(StateError, match="heat capacity of CO2"):
            heat_capacity("CO2", 210.0)
        with pytest.raises(StateError, match="not nan K"):
            viscosity("Ar", [300.0, math.nan])
        with pytest.raises(StateError):
            mixture_viscosity({"He": 0.5, "Ar": 0.5}, 3000.0)

    def test_the_ranges_follow_the_sources(self):
        assert transport_model("CO2").temperature_range("viscosity") == (200.0, 2000.0)
        assert transport_model("D2").temperature_range("viscosity") == (30.0, 2200.0)
        xenon = transport_model("Xe")
        assert xenon.temperature_range("viscosity") == (
            xenon.equation_of_state.temperature_range
        )
        for gas in TRANSPORT_GASES:
            model = transport_model(gas)
            if gas != "O2":
                assert model.temperature_range("heat_capacity") == (
                    model.equation_of_state.temperature_range
                )
        with pytest.raises(OptmatError, match="unknown quantity"):
            xenon.temperature_range("conductivity")

    def test_arrays_keep_their_shape(self):
        t = np.array([[250.0, 300.0], [350.0, 400.0]])
        for gas in ("Ne", "D2", "Ar"):
            assert viscosity(gas, t).shape == (2, 2)
            assert np.all(np.diff(viscosity(gas, t).ravel()) > 0)
        assert heat_capacity("N2", t).shape == (2, 2)
        assert mixture_viscosity("dry air", t).shape == (2, 2)


def test_every_entry_names_curated_records():
    records = sources()
    for gas in TRANSPORT_GASES:
        model = transport_model(gas)
        for source in model.identity.sources:
            assert source in records
        assert (model.status == "validated") == (gas != "CH4")
        assert model.viscosity_uncertainty


def test_gases_without_a_model_are_refused():
    with pytest.raises(OptmatError, match="no transport model for O3"):
        transport_model("O3")


def test_the_identity_binds_the_equation_of_state():
    model = transport_model("Ar")
    identity = model.identity
    assert identity.kind == "transport"
    assert identity.name == "transport/Ar"
    assert identity.sources == ("lemmon_2004", "tegeler_1999")
    record = model.identity_record()
    assert record["fields"]["equation_of_state"]["digest"] == (
        model.equation_of_state.identity.digest
    )
    other = dataclasses.replace(model, viscosity_range=(100.0, 1000.0))
    assert other.identity.digest != identity.digest
    assert transport_model("SF6").identity.sources == (
        "guder_2009",
        "quinones_cisneros_2012",
    )
    assert transport_model("Ne").identity.sources == ("huber_2018",)


def test_the_model_checks_its_form_and_range():
    model = transport_model("Ar")
    with pytest.raises(DataError, match="unknown viscosity form"):
        dataclasses.replace(model, viscosity_form="guess")
    with pytest.raises(DataError, match="empty viscosity range"):
        dataclasses.replace(model, viscosity_range=(300.0, 200.0))


def test_missing_formulas_are_errors():
    with pytest.raises(DataError, match="no formula"):
        transport._formula("huber_2018", "eq6_eq7_He")
