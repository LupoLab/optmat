# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import numpy as np
import pytest

from optmat.exceptions import DataError
from optmat.gas.linewidths import (
    AMAGAT,
    WAVENUMBER_TO_HZ,
    Broadening,
    ExponentialGapLaw,
    LineTable,
    LineTableLaw,
    LinewidthLaw,
)

J = np.array([0, 1, 2])


def test_a_linewidth_law_adds_dicke_and_collisional_terms():
    law = LinewidthLaw(
        dicke=300e6,
        broadening={
            "self": Broadening(50e6, 298.0, exponent=0.5),
            "Ar": Broadening(100e6, 298.0),
        },
        substitutes={"He": "Ar"},
    )
    density = 2 * AMAGAT
    np.testing.assert_allclose(law.fwhm(J, 298.0, density), 150e6 + 100e6)
    # The temperature law scales the collisional term only.
    np.testing.assert_allclose(law.fwhm(J, 4 * 298.0, density), 150e6 + 50e6)
    mixed = law.fwhm(J, 298.0, density, {"self": 0.5, "He": 0.5})
    np.testing.assert_allclose(mixed, 150e6 + 0.5 * 100e6 + 0.5 * 200e6)
    with pytest.raises(DataError, match="perturber"):
        law.fwhm(J, 298.0, density, {"N2": 1.0})
    with pytest.raises(DataError, match="positive density"):
        law.fwhm(J, 298.0, 0.0)


def test_bischel_and_dyers_polynomial_law():
    # Their Eq. (2) at 474 K: B = 51.8 + 0.152 x 176 + 4.85e-4 x 176^2 = 93.6 MHz/amagat
    # (their Table I prints 94.0 +- 1.3) and A = 309 (474/298)^0.92 = 473.6 MHz amagat.
    law = LinewidthLaw(
        dicke=309e6,
        dicke_temperature=298.0,
        dicke_exponent=0.92,
        broadening={
            "self": Broadening(51.8e6, 298.0, linear=0.152e6, quadratic=4.85e-4 * 1e6)
        },
    )
    width = law.fwhm(np.array([1]), 474.0, AMAGAT)[0]
    assert width == pytest.approx((473.6 + 93.6) * 1e6, rel=2e-3)


def test_a_line_table_interpolates_and_holds_its_ends():
    table = LineTable((2, 4), (4e9, 3e9), 295.0, (0.3, 0.5))
    np.testing.assert_allclose(table.at(np.array([0, 3, 9]), 295.0), [4e9, 3.5e9, 3e9])
    # At 2 x T0 the exponent of J = 3 is 0.4: (1/2)^0.4.
    assert table.at(np.array([3]), 590.0)[0] == pytest.approx(3.5e9 * 0.5**0.4)
    law = LineTableLaw({"self": table})
    assert law.fwhm(np.array([2]), 295.0, 3 * AMAGAT)[0] == pytest.approx(12e9)
    with pytest.raises(DataError, match="increase"):
        LineTable((4, 2), (1.0, 1.0), 295.0, (0.0, 0.0))


def test_the_exponential_gap_law_reproduces_its_evaluation():
    # Koszykowski et al.'s modified law with their constants, evaluated independently
    # (koszykowski_1987's formula record, eq44, B = 1.990 cm-1): FWHM 126.3, 99.0, 91.6,
    # 88.4 x 1e-3 cm-1/atm for J = 0, 2, 4, 6 at 295 K, printed to 0.1.
    law = ExponentialGapLaw(
        alpha=0.023,
        beta=1.67,
        delta=1.26,
        a=1.5,
        exponent=1.346,
        reference_temperature=295.0,
        b=1.990,
    )
    one_atm = AMAGAT * 273.15 / 295.0
    width = law.fwhm(np.array([0, 2, 4, 6]), 295.0, one_atm) / WAVENUMBER_TO_HZ * 1e3
    np.testing.assert_allclose(width, [126.3, 99.0, 91.6, 88.4], atol=0.06)
    with pytest.raises(DataError, match="self-broadening"):
        law.fwhm(np.array([0]), 295.0, one_atm, {"Ar": 1.0})
