# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The mixture models against the measured data, under TC-778's rule (OEP 4, section 8).

The pair records' qualification is what the comparison gives; the fallback the defaults
use agrees with the measured mixture densities and with CIPM-2007's air; the effective
cross coefficient of a multi-fluid model depends on composition.
"""

import warnings

import numpy as np
import pytest

from optmat.conditions import State
from optmat.datasets import datasets
from optmat.exceptions import ModelRangeWarning
from optmat.gas import mixture_eos
from optmat.gas.mixture_eos import mixture_equation_of_state
from optmat.gas.mixture_validation import (
    comparisons,
    effective_cross_virial,
    pair_results,
)

REQUIRED = (
    ("Ar", "He"),
    ("Ar", "Ne"),
    ("He", "Ne"),
    ("N2", "O2"),
    ("Ar", "N2"),
    ("Ar", "O2"),
)
"""The pairs of OEP 4's required mixtures: He/Ar, Ne/Ar, He/Ne and dry air."""


def test_the_records_are_the_comparisons_result():
    results = pair_results()
    records = mixture_eos._pair_records()
    for pair, record in records.items():
        assert record["qualified"] == results[pair].qualified, pair
    # No required pair's multi-fluid model qualifies (TC-779).
    assert not any(results[p].qualified for p in REQUIRED)


def test_the_rule_counts_only_stated_primary_points():
    for point in comparisons():
        if point.counted:
            assert not np.isnan(point.sigma)
        assert point.passed <= point.counted
        assert 250.0 <= point.temperature <= 350.0
    counted = {p.dataset for p in comparisons() if p.counted}
    # Brewer and Vaughn's B_12 rests on 1929 pure values and states no uncertainty:
    # their E is counted instead.
    assert "brewer_1969/E_Ar_He" in counted
    assert "brewer_1969/B12_Ar_He" not in counted


def test_the_effective_cross_coefficient_depends_on_composition():
    # Gernert's limitation (TC-743): for He-Ar at 298 K it runs from +16 to -4 cm^3/mol.
    low = effective_cross_virial("Ar", "He", 298.15, 0.1)
    high = effective_cross_virial("Ar", "He", 298.15, 0.9)
    assert low * 1e6 == pytest.approx(16.0, abs=0.1)
    assert high * 1e6 == pytest.approx(-4.2, abs=0.1)


@pytest.mark.parametrize("dataset", ["vogl_1972/Z_x0494_40C", "vogl_1972/Z_x0494_0C"])
def test_the_default_fallback_against_vogl_and_hall(dataset):
    # The He-Ne default (virial, B_12 from Brewer and Vaughn's E) against Vogl and
    # Hall's Z of the 49.4 % mixture, within twice its stated sigma.
    data = {d.key: d for d in datasets()}[dataset]
    temperature = float(str(data.parameters["temperature_k"]))
    eos = mixture_equation_of_state({"He": 0.494, "Ne": 0.506})
    assert (eos.method, eos.status) == ("virial", "validated")
    passed = [
        abs(eos.compressibility(State(pressure=float(p), temperature=temperature)) - z)
        <= 2 * s
        for p, z, s in zip(
            data.abscissa_values(),
            data.scaled_values(),
            data.standard_uncertainty(),
            strict=True,
        )
    ]
    # 18 of 18 at 313 K, 17 of 18 at 273 K (TC-779).
    assert sum(passed) >= len(passed) - 1


def _cipm_worst(eos):
    """Return the largest |Z - Z_CIPM| over 600-1100 hPa and 17-27 C, dry air."""
    a0, a1, a2, d = 1.58123e-6, -2.9331e-8, 1.1043e-10, 1.83e-11
    worst = 0.0
    for pressure in (6e4, 8e4, 101325.0, 1.1e5):
        for celsius in (17.0, 20.0, 27.0):
            temperature = celsius + 273.15
            cipm = (
                1
                - pressure / temperature * (a0 + a1 * celsius + a2 * celsius**2)
                + (pressure / temperature) ** 2 * d
            )
            got = eos.compressibility(State(pressure=pressure, temperature=temperature))
            worst = max(worst, abs(got - cipm))
    return worst


def test_the_default_air_against_cipm_2007():
    # The author's default for air (TC-796): the multi-fluid model at every
    # temperature, provisional. Against CIPM-2007, u(Z) = 1.5e-5: within 8.6e-6.
    with pytest.warns(ModelRangeWarning, match="provisional"):
        eos = mixture_equation_of_state("dry air (N2, O2, Ar)")
    assert (eos.method, eos.status) == ("multi-fluid", "provisional")
    assert eos.temperature_range[0] < 250.0 and eos.temperature_range[1] > 350.0
    assert _cipm_worst(eos) < 1e-5
    # The fallback, at Martin et al.'s 290-320 K only, within 1.8e-6 (TC-779).
    with pytest.warns(ModelRangeWarning, match="provisional"):
        fallback = mixture_equation_of_state("dry air (N2, O2, Ar)", method="virial")
    assert fallback.temperature_range == (290.0, 320.0)
    assert _cipm_worst(fallback) < 2e-6


def test_the_fallback_against_the_multi_fluid_model():
    # Section 8: measured and reported. At 10 bar and 294 K the fallback's density is
    # 2.1e-3 below the multi-fluid model's for He/Ar (the multi-fluid cross-virial
    # error) and 1.8e-4 above it for air (TC-779).
    state = State(pressure=1e6, temperature=294.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        for composition, size in (
            ({"He": 0.5, "Ar": 0.5}, 2.5e-3),
            ("dry air (N2, O2, Ar)", 2.5e-4),
        ):
            fallback = mixture_equation_of_state(composition)
            multi = mixture_equation_of_state(composition, method="multi-fluid")
            ratio = fallback.number_density(state) / multi.number_density(state) - 1
            assert abs(ratio) < size


def test_the_fallback_records_are_the_comparisons_result():
    from optmat.gas.mixture_validation import fallback_results  # noqa: PLC0415

    results = fallback_results()
    records = mixture_eos._pair_records()
    for pair in REQUIRED:
        result = results.get(pair)
        expected = result.qualified if result else False
        assert records[pair]["fallback_qualified"] == expected, pair
    assert results[("He", "Ne")].qualified
    assert not results[("Ar", "He")].qualified


def test_the_defaults_status():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        statuses = {
            name: mixture_equation_of_state(c).status
            for name, c in (
                ("He/Ne", {"He": 0.5, "Ne": 0.5}),
                ("He/Ar", {"He": 0.5, "Ar": 0.5}),
                ("Ne/Ar", {"Ne": 0.5, "Ar": 0.5}),
                ("air", "dry air (N2, O2, Ar)"),
            )
        }
    assert statuses == {
        "He/Ne": "validated",
        "He/Ar": "provisional",
        "Ne/Ar": "provisional",
        "air": "provisional",
    }


def test_a_density_abscissa_and_a_mass_density(monkeypatch):
    import dataclasses  # noqa: PLC0415

    from optmat.gas import mixture_validation  # noqa: PLC0415

    # Vogl and Hall's 40 C isotherm, restated against its model densities.
    data = {d.key: d for d in datasets()}["vogl_1972/Z_x0494_40C"]
    eos = mixture_equation_of_state({"He": 0.494, "Ne": 0.506})
    densities = [
        eos.number_density(State(pressure=float(p), temperature=313.15))
        for p in data.abscissa_values()[:3]
    ]
    restated = dataclasses.replace(
        data,
        abscissa="density",
        abscissa_unit="m-3",
        abscissae=tuple(f"{n:.12e}" for n in densities),
        values=data.values[:3],
        sigmas=data.sigmas[:3],
        rows=data.rows[:3],
        excluded_rows={},
        columns={},
    )
    points = mixture_validation._density_points(restated, "virial")
    assert [p.abscissa for p in points] == pytest.approx(densities, rel=1e-11, abs=0.0)
    # 1000 mol/m^3 of helium is 4.002602 kg/m^3, 4.002602e-3 g/cm^3.
    molar = 1000.0 * 6.02214076e23
    mass = mixture_validation._in_unit(molar, "kg/m3", {"He": 1.0})
    assert mass == pytest.approx(4.002602, rel=1e-12, abs=0.0)
    assert mixture_validation._in_unit(molar, "g/cm3", {"He": 1.0}) == pytest.approx(
        4.002602e-3, rel=1e-12, abs=0.0
    )


def test_a_pair_without_a_cross_virial_has_no_fallback_comparison(monkeypatch):
    from optmat.gas import mixture_validation  # noqa: PLC0415

    records = {
        k: v for k, v in mixture_eos._pair_records().items() if k != ("He", "Ne")
    }
    monkeypatch.setattr(mixture_validation, "_pair_records", lambda: records)
    points = mixture_validation.fallback_comparisons.__wrapped__()
    assert ("He", "Ne") not in {p.pair for p in points}


def test_a_fallback_outside_its_measured_temperatures_is_not_compared(monkeypatch):
    from optmat.gas import mixture_validation  # noqa: PLC0415

    # With Martin et al.'s E (290-320 K) as He-Ar's cross coefficient, Blancett et
    # al.'s 273.15 and 323.15 K isotherms lie outside the fallback's range.
    records = dict(mixture_eos._pair_records())
    records["Ar", "He"] = dict(records["Ar", "He"]) | {
        "cross_virial": "martin_1982/E_Ar_He"
    }
    monkeypatch.setattr(mixture_eos, "_pair_records", lambda: records)
    data = {d.key: d for d in datasets()}["blancett_1970/Z_xAr_0.4065_323.15K"]
    assert mixture_validation._density_points(data, "virial") == []
