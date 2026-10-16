# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import dataclasses
import math

import numpy as np
import pytest

from optmat.gas.fitting import HARTREE_NM
from optmat.gas.models.fitted import fitted_models


def _beta_k(model, omega, order):
    """d^k (omega alpha) / d omega^k by a polynomial fit on a nine-point stencil."""
    # 2 % of omega: a narrower stencil makes the third derivative round-off limited
    # (5e-3 gave 3e-3 scatter against the parameter steps below; 2e-2 gives 3e-5).
    step = omega * 2e-2
    offsets = np.arange(-4, 5)
    values = (omega + offsets * step) * model.polarizability_au(omega + offsets * step)
    coefficients = np.polyfit(offsets, values, 8)[::-1]
    return float(coefficients[order] * math.factorial(order) / step**order)


@pytest.mark.parametrize("order", [2, 3])
def test_dispersion_uncertainty_matches_a_numerical_propagation(order):
    model = fitted_models()["fit/Ar"]
    omega = HARTREE_NM / 800.0
    strengths = np.array(model.strengths)
    energies = np.array(model.energies)
    parameters = np.concatenate([strengths, energies[model.n_lines :]])
    gradient = []
    for i in range(len(parameters)):
        shifted = []
        for sign in (1, -1):
            p = parameters.copy()
            p[i] *= 1 + sign * 1e-4
            varied = dataclasses.replace(
                model,
                strengths=tuple(p[: len(strengths)]),
                energies=(*energies[: model.n_lines], *p[len(strengths) :]),
            )
            shifted.append(_beta_k(varied, omega, order))
        gradient.append((shifted[0] - shifted[1]) / (2e-4 * parameters[i]))
    gradient = np.array(gradient)
    covariance = np.asarray(model.covariance) * max(1.0, model.reduced_chi2)
    expected = np.sqrt(gradient @ covariance @ gradient) / abs(
        _beta_k(model, omega, order)
    )
    analytic = model.dispersion_sigma_relative(np.array([omega]), order)
    # The numerical route (stencil derivatives, central parameter differences) is good
    # to a few 1e-4 of itself; the analytic one is exact.
    assert analytic[0] == pytest.approx(expected, rel=2e-3, abs=0.0)


def test_dispersion_uncertainty_needs_order_two_or_more():
    with pytest.raises(ValueError, match="at least 2"):
        fitted_models()["fit/He"].dispersion_sigma_relative(np.array([0.05]), 1)


def test_heliums_dispersion_coefficient_meets_yang_et_als_measurement():
    """fit/He's refractive A_2 against Yang, Stone & Egan's two-colour value (TC-705).

    The measurement (yang_2025, formula ``he_A2``) sees the retardation term: without it
    the fit was 4.8 of the measurement's sigma away, with it 0.6.
    """
    import re  # noqa: PLC0415

    from scipy import constants  # noqa: PLC0415

    from optmat.datasets import sources  # noqa: PLC0415

    formula = next(f for f in sources()["yang_2025"].formulas if f["id"] == "he_A2")
    coefficients = formula["coefficients"]
    assert isinstance(coefficients, dict)
    measured = float(coefficients["A_2"])
    digits = re.search(r"A_2 1\.332468\((\d)\)", str(formula["printed"]))
    assert digits is not None
    sigma = int(digits[1]) * 1e-10  # the last printed digit, 10^-4 x 10^-6
    model = fitted_models()["fit/He"]
    # The refractive polarizability's omega^2 coefficient, converted to A_2 of Yang et
    # al.'s Eq. (4): A_2 = (4 pi / 3) N_A a0^3 S (2 pi 10^14 Hz hbar / E_h)^2, cm^3/mol.
    s4 = model.cauchy_moment(-4) + model.retardation_au
    a0_cm = constants.physical_constants["Bohr radius"][0] * 100
    molar = 4 * math.pi / 3 * constants.Avogadro * a0_cm**3
    frequency = (
        2
        * math.pi
        * 1e14
        * constants.hbar
        / constants.physical_constants["Hartree energy"][0]
    )
    a2 = s4 * molar * frequency**2
    # Two of the measurement's sigma: the fit's own uncertainty of S(-4) is 4 ppm, about
    # 0.7 of that sigma (6e-10 of 1.33e-4 is 4.5 ppm).
    assert abs(a2 - measured) < 2 * sigma
