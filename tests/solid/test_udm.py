# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import itertools
import math
import os
from pathlib import Path

import numpy as np
import pytest
from scipy import integrate

from optmat.datasets import sources
from optmat.exceptions import DataError
from optmat.solid import udm

# The published models whose newAD2 parameter files the source records hold as terms.
CURATED = {
    "LithosilQ2-SPIE9890": ("franta_2016", "udm_lithosil_q2"),
    "SiO2-X2551-SPIE9628": ("franta_2015b", "udm_SiO2_X2551"),
    "MgF2-X2935-SPIE9628": ("franta_2015b", "udm_MgF2_X2935"),
    "Al2O3-X2890-SPIE9628": ("franta_2015b", "udm_Al2O3_X2890"),
    "c-Si": ("franta_2017", "udm_c_si"),
}
SILICA = ["LithosilQ2-SPIE9890"]
OTHERS = ["SiO2-X2551-SPIE9628", "MgF2-X2935-SPIE9628", "Al2O3-X2890-SPIE9628", "c-Si"]
# The authors' files, attached to the papers in the curators' library (the library
# lane): the tables, and the models of the SPIE 9890 films, which the source records
# do not hold.
LIBRARY_FILES = {
    "LithosilQ2-SPIE9890": "Franta_2016_OpticalCharacterizationSiO",
    "SiO2-X2546-SPIE9890": "Franta_2016_OpticalCharacterizationSiO",
    "SiO2-X2551-SPIE9890": "Franta_2016_OpticalCharacterizationSiO",
    "SiO2-X2551-SPIE9628": "Franta_2015_DispersionModelOptical",
    "MgF2-X2935-SPIE9628": "Franta_2015_DispersionModelOptical",
    "Al2O3-X2890-SPIE9628": "Franta_2015_DispersionModelOptical",
    "c-Si": "Franta_2017_TemperatureDependentDispersion",
}
LIBRARY_SILICA = ["LithosilQ2-SPIE9890", "SiO2-X2546-SPIE9890", "SiO2-X2551-SPIE9890"]

# The authors' tables are computed on the grid E = 0.01 exp(i ln(1e4) / 4000) eV and
# printed to seven significant digits; E itself is printed rounded, which on a steep
# phonon peak moves n by many units of the last digit, so the exact grid is used.
GRID = 0.01 * np.exp(np.arange(4001) * math.log(1e4) / 4000)


def _terms(name) -> dict[str, str]:
    source, formula = CURATED[name]
    entry = next(f for f in sources()[source].formulas if f["id"] == formula)
    terms = entry["terms"]
    assert isinstance(terms, dict)
    return {str(k): str(v) for k, v in terms.items()}


def _model(name):
    terms = _terms(name)
    attributes = udm.parse_attributes(terms["attributes"], name)
    values = {k: float(v) for k, v in terms.items() if k != "attributes"}
    return udm.from_parameters(attributes, values)


def _attached(name, suffix):
    """Return an attached newAD2 file from the library, or skip without it."""
    if "PAPER_LIBRARY" not in os.environ:
        pytest.skip("set PAPER_LIBRARY to the paper library to read the newAD2 files")
    paper = Path(os.environ["PAPER_LIBRARY"]) / LIBRARY_FILES[name]
    found = sorted(paper.glob(f"data/*/{name}{suffix}"))
    assert len(found) == 1, (paper, name)
    return found[0]


def _library_model(name):
    return udm.from_file(_attached(name, ".par"))


def _nk(eps):
    # k from eps_i / (2 n): numpy's complex square root loses a k below ~1e-16 n.
    n = np.sqrt((np.abs(eps) + eps.real) / 2)
    return n, eps.imag / (2 * n)


def _last_digit_units(ours, printed):
    unit = 1e-6 * 10 ** np.floor(np.log10(np.maximum(printed, 1e-300)))
    return np.abs(ours - printed) / unit


@pytest.mark.library
@pytest.mark.parametrize("name", sorted(LIBRARY_FILES))
def test_the_authors_tables_are_reproduced_above_1_ev(name):
    table = np.loadtxt(_attached(name, ".Enk"))
    np.testing.assert_allclose(table[:, 0], GRID, rtol=6e-7)
    n, k = _nk(_library_model(name).permittivity(GRID))
    above = GRID > 1
    # Half a unit of the seventh printed digit (plus float slack): exact reproduction.
    assert _last_digit_units(n, table[:, 1])[above].max() <= 0.5001
    assert _last_digit_units(k, table[:, 2])[above].max() <= 0.5001


@pytest.mark.library
@pytest.mark.parametrize("name", LIBRARY_SILICA)
def test_silica_tables_below_1_ev_differ_only_by_the_programs_arithmetic(name):
    # Below 1 eV the authors' program replaces the higher-energy term's real part by
    # its E -> 0 limit under 5e-3 E_xi (0.1 eV), and its closed forms lose digits to
    # cancellation; both move the printed n and k by up to 2 units of the last digit.
    # This implementation is exact there (checked against 50-digit arithmetic, TC-321).
    table = np.loadtxt(_attached(name, ".Enk"))
    n, k = _nk(_library_model(name).permittivity(GRID))
    below = GRID <= 1
    assert _last_digit_units(n, table[:, 1])[below].max() < 2.0
    assert _last_digit_units(k, table[:, 2])[below].max() < 2.1


@pytest.mark.library
def test_the_urbach_tail_follows_the_program_not_eq_33():
    # Reading Eq. (33)'s E_m = (E_g + E_h)/2 as the tail's upper energy moves n by
    # percents; the program's E_h reproduces the table (TC-321).
    import dataclasses  # noqa: PLC0415

    model = _model("LithosilQ2-SPIE9890")
    terms = tuple(
        dataclasses.replace(
            t, shape=(t.shape[0], t.shape[1], (t.shape[1] + t.shape[2]) / 2)
        )
        if t.kind == "urbach"
        else t
        for t in model.terms
    )
    table = np.loadtxt(_attached("LithosilQ2-SPIE9890", ".Enk"))
    n, _ = _nk(udm.UniversalDispersionModel(terms).permittivity(GRID))
    assert np.max(np.abs(n / table[:, 1] - 1)) > 0.01


def _kramers_kronig(model, energy, edges):
    # eps_r(E) - 1 = (2/pi) P int eps_i(x) x / (x^2 - E^2) dx, the pole handled by
    # quad's Cauchy weight on 1/(x - E) with x / (x + E) in the integrand.
    def imag(x):
        return float(np.imag(model.permittivity(np.array([x])))[0])

    total = 0.0
    pieces = [*edges, 2000.0]
    for low, high in itertools.pairwise(pieces):
        if low < energy < high:
            total += integrate.quad(
                lambda x: imag(x) * x / (x + energy),
                low,
                high,
                weight="cauchy",
                wvar=energy,
                limit=400,
            )[0]
        else:
            total += integrate.quad(
                lambda x: imag(x) * x / (x * x - energy**2), low, high, limit=400
            )[0]
    return 2 / math.pi * total


@pytest.mark.parametrize("energy", [0.5, 2.0, 6.0, 10.0])
def test_the_real_part_is_the_kramers_kronig_transform_of_the_imaginary(energy):
    # The electronic terms of Lithosil Q2 only (the phonon peaks are narrow and many).
    # Beyond 2000 eV only the higher-energy term absorbs, eps_i = 3 E_xi N (x -
    # E_xi)^2 / x^5; its contribution (about 1e-6) is added analytically to order
    # (E/x)^2. Quadrature: 1e-8.
    model = _model("LithosilQ2-SPIE9890")
    electronic = udm.UniversalDispersionModel(
        tuple(t for t in model.terms if t.kind != "gaussian")
    )
    eg, eh = 8.36270918663193, 39.6024678126438
    exi = next(t.shape[0] for t in electronic.terms if t.kind == "higher")
    edges = [eg / 2, eg, 9.5, 12.0, 14.0, 20.0, exi, eh, eh + eg / 2, 100.0]
    expected = _kramers_kronig(electronic, energy, sorted(edges))
    strength = next(t.strength for t in electronic.terms if t.kind == "higher")
    top = 2000.0
    tail = (
        3
        * exi
        * strength
        * (1 / (3 * top**3) - exi / (2 * top**4) + (exi**2 + energy**2) / (5 * top**5))
    )
    expected += 2 / math.pi * tail
    actual = np.real(electronic.permittivity(np.array([energy])))[0] - 1
    assert actual == pytest.approx(expected, rel=1e-8)


@pytest.mark.parametrize("name", SILICA)
def test_every_term_has_unit_transition_strength(name):
    # Eq. (3): int eps_i^0(E) E dE = 1 for each normalised term; quadrature 1e-9.
    for term in _model(name).terms:
        shape = udm._SHAPES[term.kind]

        def strength(x, shape=shape, term=term):
            return float(shape(np.array([x]), *term.shape)[1][0]) * x

        if term.kind == "gaussian":
            centre, width = term.shape
            low, high = max(1e-9, centre - 40 * width), centre + 40 * width
            points = [centre]
        else:
            edges = udm._edges(term)
            low, high, points = edges[0] / 2, 1e5, edges
        total = integrate.quad(strength, low, high, points=points, limit=500)[0]
        if term.kind == "higher":
            # The tail beyond `high`, exactly: int 3 E_xi (x - E_xi)^2 / x^4 dx.
            exi = term.shape[0]
            total += 3 * exi * (1 / high - exi / high**2 + exi**2 / (3 * high**3))
        assert total == pytest.approx(1.0, rel=1e-9), term


@pytest.mark.parametrize("energy", [0.03, 0.137, 1.0, 5.0, 9.0, 11.3, 30.0])
def test_the_derivative_is_exact(energy):
    # Complex-step against a fourth-order central difference with h = 1e-4 E, whose
    # truncation error is below 1e-9 relative where the function is smooth.
    model = _model("LithosilQ2-SPIE9890")
    h = 1e-4 * energy
    e = energy + h * np.array([-2.0, -1.0, 1.0, 2.0])
    f = model.permittivity(e)
    numeric = (f[0] - 8 * f[1] + 8 * f[2] - f[3]) / (12 * h)
    exact = model.dpermittivity_denergy(np.array([energy]))[0]
    assert abs(exact - numeric) <= 1e-7 * abs(exact) + 1e-12


def test_the_model_is_passive():
    for name in SILICA + OTHERS:
        assert np.all(np.imag(_model(name).permittivity(GRID)) >= 0), name


def test_the_low_energy_series_joins_the_closed_form():
    # At the switch the even series and the closed form agree to the latter's
    # rounding (the series' truncation is below 1e-14 there).
    for name in SILICA + OTHERS:
        for term in _model(name).terms:
            if term.kind == "gaussian":
                continue
            energy = np.array([udm._SMALL[term.kind] * term.lowest])
            closed = udm._SHAPES[term.kind](energy, *term.shape)[0][0]
            series = sum(
                m * energy[0] ** (2 * k) for k, m in enumerate(udm._moments(term))
            )
            assert series == pytest.approx(closed, rel=1e-10), (name, term.kind)


@pytest.mark.library
@pytest.mark.parametrize("name", sorted(CURATED))
def test_the_source_records_hold_the_parameter_files_values(name):
    # The terms are the file's strings, exactly: its attributes and every value.
    lines = _attached(name, ".par").read_text(encoding="utf-8").splitlines()
    attributes, values = udm.parse_parameters("\n".join(lines), name)
    terms = _terms(name)
    assert tuple(terms["attributes"].split(":")) == attributes
    assert {k: float(v) for k, v in terms.items() if k != "attributes"} == values
    printed = {
        line.split("=")[0].strip(): line.split("=")[1].strip()
        for line in lines[1:]
        if "=" in line and line.split("=")[0].strip() in values
    }
    assert {k: v for k, v in terms.items() if k != "attributes"} == printed


def test_an_unsupported_attribute_is_refused(tmp_path):
    par = tmp_path / "x.par"
    par.write_text(
        "dEg:ut\nNvc = 1\nEg = 5\nEh = 20\nEp = 0.1\nNut = 1\nEu = 0.1\n",
        encoding="utf-8",
    )
    with pytest.raises(DataError, match="not supported"):
        udm.read_parameter_file(par)


@pytest.mark.parametrize("name", ["LithosilQ2-SPIE9890", "MgF2-X2935-SPIE9628"])
def test_the_moments_quadrature_agrees_with_adaptive_quadrature(name):
    # The composite Gauss-Legendre moments against scipy's adaptive quadrature at
    # 1e-13; they agree to 1e-14 (checked for all seven files, 2026-09-26).
    for term in _model(name).terms:
        if term.kind == "gaussian":
            continue
        shape = udm._SHAPES[term.kind]
        start, edges = term.lowest / 2, udm._edges(term)
        top = 50 * max(edges)
        points = [e for e in edges if start < e < top]
        if term.kind == "exciton":
            points.append(term.shape[2])
        for k in (0, 3, 11):
            power = 2 * k + 1

            def integrand(x, shape=shape, term=term, power=power):
                return float(shape(np.array([x]), *term.shape)[1][0]) / x**power

            value = integrate.quad(
                integrand, start, top, points=points, limit=1000, epsabs=0, epsrel=1e-13
            )[0]
            if term.kind == "higher":
                exi, q = term.shape[0], power + 2
                value += (
                    3
                    * exi
                    * (
                        1 / (q * top**q)
                        - 2 * exi / ((q + 1) * top ** (q + 1))
                        + exi**2 / ((q + 2) * top ** (q + 2))
                    )
                )
            assert udm._moments(term)[k] == pytest.approx(
                2 / math.pi * value, rel=1e-12
            )
