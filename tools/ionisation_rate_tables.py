# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The tables of the noble-gas ionisation pages.

The tables live in ``docs/ionisation_rates.md`` and the
``docs/user-guide/ionisation-*.md`` pages it leads to; the historical comparisons
(``ionisation-rate-legacy``) live in ``tests/data/historical_documentation_tables.md``.

    uv run python tools/ionisation_rate_tables.py

prints the tables; ``tests/test_ionisation_rates_page.py`` checks that the pages hold
exactly what this script prints, so the pages are regenerated whenever a model, a
constant or a reference changes (and its text is updated where a finding changes):

- ``ionisation-rate-models``: what each named model is made of;
- ``ionisation-rate-electrons``: the constants of the default model per gas;
- ``ionisation-rate-atoms``: published many-electron calculations over one electron in
  their pulses;
- ``ionisation-rate-weak-field``: the atoms' weak-field rate over the model's one
  electron, from Tolstikhina et al.'s coefficients;
- ``ionisation-rate-static``: the static rates against their references;
- ``ionisation-rate-helium``: helium at 780 and 390 nm against Parker et al.'s rates;
- ``ionisation-rate-helium-two-electron``: helium at 390 nm with two electrons, haCC
  and one electron, and Parker et al.'s rate over the same pulse;
- ``ionisation-rate-helium-matched``: helium's two-electron, haCC and one-electron
  calculations against each other on box-converged measures, and against Sato et al.;
- ``ionisation-rate-xenon-rmt``: RMT's xenon against one electron in its field, and
  the default against it;
- ``ionisation-rate-xenon-measured``: the default against Hartmann et al.'s measured
  neutral abundance;
- ``ionisation-rate-rmatrix``: He, Ne and Ar at 390 nm against van der Hart's
  all-electron R-matrix rates;
- ``ionisation-rate-tdse``: the one-electron yields, with each Coulomb factor;
- ``ionisation-rate-tdse42``: argon's one-electron yield of a 42 fs pulse;
- ``ionisation-rate-saturation``: Uiterwaal et al.'s saturation intensities at 248.6
  nm;
- ``ionisation-rate-wahlstrand``, ``ionisation-rate-wallace``,
  ``ionisation-rate-wallace-points``: the measured yields;
- ``ionisation-rate-accuracy``: the stated uncertainty of each gas;
- ``ionisation-rate-legacy``: the default against Luna.jl's and solprop's rates;
- ``ionisation-rate-phase``: the phase shifts of Nie et al. and Nambu et al.;
- ``ionisation-rate-luna``, ``ionisation-rate-luna-tests``: the default against
  Luna.jl's rates and its tests' yields.
"""

import dataclasses
import importlib.util
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np

from optmat.constants import AU_ANGULAR_FREQUENCY, AU_FIELD, E_CHARGE, EPS0, C
from optmat.datasets import datasets
from optmat.exceptions import ModelRangeWarning
from optmat.gas import ionisation_validation as validation
from optmat.gas.ionisation_rates import ionisation_rate_model, ionisation_rate_models

ATOMS = ("He", "Ne", "Ar", "Kr", "Xe")
COULOMB_MODELS = (
    ("one_electron", "`one_electron` (the default's factor)"),
    ("one_electron_static", "static"),
    ("ppt_ilkov", "`ppt_ilkov`"),
    ("ppt_popruzhenko", "`ppt_popruzhenko`"),
)


def _table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def _wcm2(intensity: float) -> str:
    return f"{intensity / 1e4:.2e}"


def models_table() -> str:
    """Return what each named model is made of.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for name, models in ionisation_rate_models().items():
        model = models[next(iter(models))]
        electrons = sorted({m.electrons for m in models.values()})
        rows.append(
            [
                f"`{name}`",
                model.coulomb,
                model.static_correction,
                model.structure,
                "yes" if any(len(m.channels) > 1 for m in models.values()) else "no",
                ", ".join(f"{e:g}" for e in electrons),
            ]
        )
    return _table(
        [
            "Model",
            "Coulomb factor",
            "Static correction",
            "Coefficient",
            "Fine structure",
            "N_eff",
        ],
        rows,
    )


def electrons_table() -> str:
    """Return the constants of the default model per gas.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for gas in ATOMS:
        model = ionisation_rate_model(gas, 800e-9)
        upper = model.channels[-1].ionization_energy / E_CHARGE
        barrier = model.barrier_suppression_field
        rows.append(
            [
                gas,
                f"{model.electrons:g}",
                f"{model.structure_coefficient:g}",
                f"{model.static_constants[0]:g}",
                f"{model.channels[0].ionization_energy / E_CHARGE:.4f}",
                "" if len(model.channels) == 1 else f"{upper:.4f}",
                f"{barrier / AU_FIELD:.4f}",
                _wcm2(0.5 * EPS0 * C * barrier**2),
            ]
        )
    return _table(
        [
            "Gas",
            "N_eff",
            "C_l",
            "alpha",
            "I_p (eV)",
            "I_p to 2P1/2 (eV)",
            "F_b (au)",
            "I at F_b (W/cm2)",
        ],
        rows,
    )


def static_table() -> str:
    """Return the static rates of the noble gases over their references.

    Returns
    -------
    str
        A Markdown table: the model over the reference at the fields the page shows.
    """
    fields = (0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.15, 0.2, 0.25, 0.3)
    rows = []
    for gas in ATOMS:
        for comparison in validation.static_rates(gas):
            field = comparison.x / AU_FIELD
            cells = []
            for wanted in fields:
                near = np.flatnonzero(np.abs(field / wanted - 1) < 0.07)
                cells.append(
                    "" if near.size == 0 else f"{comparison.ratio[near[0]]:.2f}"
                )
            rows.append([gas, comparison.label, f"`{comparison.model}`", *cells])
    return _table(
        ["Gas", "Reference", "Model", *(f"{f:g} au" for f in fields)],
        rows,
    )


def helium_table() -> str:
    """Return helium's rate at 780 and 390 nm over Parker et al.'s.

    Returns
    -------
    str
        A Markdown table, with each Coulomb factor.
    """
    intensities = np.array([1, 2, 3, 4, 6, 8, 10, 14]) * 1e18
    rows = []
    for wavelength in (780, 390):
        rate = ionisation_rate_model("He", wavelength * 1e-9)
        gamma = rate.keldysh(validation.amplitude(intensities))
        expected = validation.parker_rate_au(intensities, wavelength)
        columns = []
        for coulomb in ("frequency", "static", "ilkov", "popruzhenko"):
            varied = dataclasses.replace(rate, coulomb=coulomb)
            predicted = varied.cycle_averaged_rate(validation.amplitude(intensities))
            columns.append(predicted / AU_ANGULAR_FREQUENCY / expected)
        for i, intensity in enumerate(intensities):
            rows.append(
                [
                    f"{wavelength} nm",
                    _wcm2(float(intensity)),
                    f"{gamma[i]:.2f}",
                    f"{expected[i]:.2e}",
                    *(f"{column[i]:.2f}" for column in columns),
                ]
            )
    return _table(
        [
            "Wavelength",
            "I (W/cm2)",
            "gamma",
            "Parker (au)",
            "`default`",
            "with the static factor",
            "with Ilkov's factor",
            "with Popruzhenko's",
        ],
        rows,
    )


RMATRIX_ROWS = {
    "vanderhart_2006/fig6_floquet_He": (1.1, 1.5, 2.0, 2.4),
    "vanderhart_2006/fig6_floquet_Ne": (1.1, 1.5, 2.0),
    "vanderhart_2006/fig6_floquet_Ar": (0.65, 1.0, 1.5, 1.9),
    "vanderhart_2007/table1_one_state_Ar": (0.48, 0.75, 1.5),
    "vanderhart_2007/fig4_td_rmatrix_Ar": (0.25, 1.25, 1.75, 2.5, 3.0),
}


def rmatrix_table() -> str:
    """Return the models at 390 nm over van der Hart's all-electron rates.

    Returns
    -------
    str
        A Markdown table: selected intensities of each reference, with helium's
        two-electron rate of Parker et al. over the same reference.
    """
    rows = []
    for gas in ("He", "Ne", "Ar"):
        one = {
            c.reference: c for c in validation.rmatrix_rates(gas, "one_electron_static")
        }
        candidate = {
            c.reference: c for c in validation.rmatrix_rates(gas, "one_electron")
        }
        for comparison in validation.rmatrix_rates(gas):
            chosen = RMATRIX_ROWS[comparison.reference]
            gamma = ionisation_rate_model(gas, 390e-9).keldysh(
                validation.amplitude(comparison.x)
            )
            for wanted in chosen:
                k = int(np.argmin(np.abs(comparison.x / 1e18 - wanted)))
                parker = (
                    validation.parker_rate_au(comparison.x[k : k + 1], 390)[0]
                    / comparison.expected[k]
                    if gas == "He"
                    else None
                )
                rows.append(
                    [
                        gas,
                        f"`{comparison.reference}`",
                        _wcm2(float(comparison.x[k])),
                        f"{gamma[k]:.2f}",
                        f"{comparison.expected[k] * AU_ANGULAR_FREQUENCY:.2e}",
                        f"{comparison.ratio[k]:.2f}",
                        f"{one[comparison.reference].ratio[k]:.2f}",
                        f"{candidate[comparison.reference].ratio[k]:.2f}",
                        "" if parker is None else f"{parker:.2f}",
                    ]
                )
    return _table(
        [
            "Gas",
            "Reference",
            "I (W/cm2)",
            "gamma",
            "Reference (1/s)",
            "`default`",
            "`one_electron_static`",
            "`one_electron` (the default's factor)",
            "Parker et al.",
        ],
        rows,
    )


def long_pulse_table() -> str:
    """Return the models over argon's one-electron yield of a 42 fs pulse.

    Returns
    -------
    str
        A Markdown table: the calculation's ionised fraction, the models' yields over
        it, and it over Wahlstrand et al.'s measured yield of the atom.
    """
    one = validation.long_pulse_yield("one_electron_static")
    candidate = validation.long_pulse_yield("one_electron")
    default = validation.long_pulse_yield("default")
    measured = validation.measured_yields("Ar", intensities=list(one.x))
    gamma = ionisation_rate_model("Ar", 800e-9).keldysh(validation.amplitude(one.x))
    rows = [
        [
            _wcm2(float(one.x[k])),
            f"{gamma[k]:.2f}",
            f"{one.expected[k]:.3e}",
            f"{one.ratio[k]:.2f}",
            f"{candidate.ratio[k]:.2f}",
            f"{default.ratio[k]:.2f}",
            f"{one.expected[k] / measured.expected[k]:.2f}",
        ]
        for k in range(len(one.x))
    ]
    return _table(
        [
            "I (W/cm2)",
            "gamma",
            "TDSE: ionised fraction",
            "`one_electron_static`",
            "`one_electron`",
            "`default`",
            "TDSE / measured (the atom)",
        ],
        rows,
    )


def tdse_table() -> str:
    """Return the one-electron yields over the calculations', per Coulomb factor.

    Returns
    -------
    str
        A Markdown table: model over the calculation's ionised fraction, for the
        cycle-averaged rate on the envelope and, for the default's factor, the
        instantaneous rate on the carrier; every calculation of
        :data:`~optmat.gas.ionisation_validation.TDSE_YIELDS`.
    """
    rows = []
    for gas in validation.TDSE_YIELDS:
        envelope = {
            name: validation.tdse_yields(gas, name, carrier=False)
            for name, _ in COULOMB_MODELS
        }
        carrier = validation.tdse_yields(gas, "one_electron", carrier=True)
        for k, ionised in enumerate(envelope["one_electron"]):
            dataset = next(d for d in datasets() if d.key == ionised.reference)
            relative = dataset.uncertainty.relative_rows
            sigma = dict(
                zip(
                    dataset.abscissa_values(),
                    (relative[row] for row in dataset.rows),
                    strict=True,
                )
            )
            rate = ionisation_rate_model(gas, ionised.wavelength, "one_electron_static")
            for i, intensity in enumerate(ionised.x):
                rows.append(
                    [
                        gas,
                        f"{ionised.wavelength * 1e9:.0f}",
                        _wcm2(float(intensity)),
                        f"{float(rate.keldysh(validation.amplitude(intensity))):.2f}",
                        f"{ionised.expected[i]:.3e}",
                        f"{100 * sigma[intensity]:.1f}",
                        *(
                            f"{envelope[name][k].ratio[i]:.2f}"
                            for name, _ in COULOMB_MODELS
                        ),
                        f"{carrier[k].ratio[i]:.2f}",
                    ]
                )
    return _table(
        [
            "Gas",
            "nm",
            "I (W/cm2)",
            "gamma",
            "TDSE: ionised fraction",
            "± (%)",
            *(label for _, label in COULOMB_MODELS),
            "the default's, on the carrier",
        ],
        rows,
    )


def wahlstrand_table() -> str:
    """Return the yields of a 42 fs pulse over Wahlstrand et al.'s.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for gas in ("Ar", "Kr", "Xe"):
        peaks = list(validation.measured_yields(gas).x[[0, 3, 6, 9]])
        default = validation.measured_yields(gas, intensities=peaks)
        others = [
            validation.measured_yields(gas, name, intensities=peaks)
            for name in ("default_static", "ppt_ilkov", "ppt_popruzhenko")
        ]
        rate = ionisation_rate_model(gas, default.wavelength)
        for i, peak in enumerate(peaks):
            rows.append(
                [
                    gas,
                    f"{peak / 1e16:.0f}",
                    f"{float(rate.keldysh(validation.amplitude(peak))):.2f}",
                    f"{default.expected[i]:.2e}",
                    f"{default.ratio[i]:.2f}",
                    *(f"{other.ratio[i]:.2f}" for other in others),
                ]
            )
    return _table(
        [
            "Gas",
            "I (TW/cm2)",
            "gamma",
            "Measured yield",
            "`default`",
            "`default_static`",
            "`ppt_ilkov`",
            "`ppt_popruzhenko`",
        ],
        rows,
    )


def wallace_table() -> str:
    """Return the intensity scales at which the models match Wallace et al.'s yields.

    Returns
    -------
    str
        A Markdown table.
    """
    names = (
        "default",
        "default_static",
        "one_electron_static",
        "ppt_ilkov",
        "ppt_popruzhenko",
    )
    rows = []
    for gas in ("Ar", "Kr", "Xe"):
        tests = [validation.intensity_scale_test(gas, name) for name in names]
        rows.append(
            [
                gas,
                f"{100 * tests[0].calibration:.1f} %",
                *(f"{test.scale:.2f}" for test in tests),
            ]
        )
    return _table(["Gas", "Calibration", *(f"`{name}`" for name in names)], rows)


def wallace_points_table() -> str:
    """Return the models against the measured points of krypton's yield.

    Returns
    -------
    str
        A Markdown table.
    """
    names = (
        "default",
        "default_static",
        "one_electron_static",
        "ppt_ilkov",
        "ppt_popruzhenko",
    )
    rows = []
    for name in names:
        test = validation.measured_points_test("Kr", name)
        rows.append(
            [
                f"`{name}`",
                f"{test.scale:.2f}",
                f"{test.chi_square:.0f}",
                f"{test.chi_square_unscaled:.0f}",
                f"{100 * test.spread:.1f} %",
            ]
        )
    return _table(
        ["Model", "Scale", "chi2 at it", "chi2 at a scale of 1", "Spread"], rows
    )


def accuracy_table() -> str:
    """Return the stated uncertainty of the default model per gas.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for gas in ATOMS:
        model = ionisation_rate_model(gas, 800e-9)
        assert math.isfinite(model.uncertainty_factor)
        rows.append(
            [
                gas,
                f"{model.uncertainty_factor:g}",
                model.reasons["uncertainty_factor"],
            ]
        )
    return _table(["Gas", "Factor", "Basis"], rows)


def phase_shift_table() -> str:
    """Return the models against the phase shifts of Nie et al. and Nambu et al.

    Returns
    -------
    str
        A Markdown table, one row per dataset whose beams are stated: the scale on the
        measurement's intensities at which each model fits, and the default's phase
        over the measured one at the measurement's intensities (geometric mean and
        range).
    """
    rows = []
    references = sorted(
        (d for d in datasets() if d.quantity == "ionization_phase_shift"),
        key=lambda d: (d.source, float(str(d.parameters["wavelength_nm"])), d.id),
    )
    for dataset in references:
        if "m2" not in dataset.parameters:
            continue
        names = ("default", "default_static", "ppt_ilkov", "ppt_popruzhenko")
        if dataset.material not in ATOMS:
            names = ("default",)
        tests = {name: validation.phase_shift_test(dataset.key, name) for name in names}
        default = tests["default"]
        rows.append(
            [
                f"`{dataset.key}`",
                f"{float(str(dataset.parameters['wavelength_nm'])):.0f}",
                str(dataset.parameters["polarization"]),
                f"{_wcm2(float(np.min(default.intensity)))} to "
                f"{_wcm2(float(np.max(default.intensity)))}",
                *(
                    f"{tests[name].scale:.2f}" if name in tests else ""
                    for name in (
                        "default",
                        "default_static",
                        "ppt_ilkov",
                        "ppt_popruzhenko",
                    )
                ),
                f"{float(np.exp(np.mean(np.log(default.ratio)))):.2f} "
                f"({float(np.min(default.ratio)):.2f} to "
                f"{float(np.max(default.ratio)):.2f})",
            ]
        )
    return _table(
        [
            "Dataset",
            "nm",
            "Polarisation",
            "I (W/cm2)",
            "scale: `default`",
            "`default_static`",
            "`ppt_ilkov`",
            "`ppt_popruzhenko`",
            "`default` / measured",
        ],
        rows,
    )


LEGACY = (
    Path(__file__).resolve().parents[1]
    / "tests"
    / "data"
    / "legacy_ionisation_rates.json"
)
LEGACY_INTENSITIES = (5e13, 1e14, 2e14, 5e14)  # W/cm^2


def _legacy_index(values: list[float], intensity: float) -> int:
    (index,) = np.flatnonzero(np.isclose(values, intensity, rtol=1e-9))
    return int(index)


def legacy_table() -> str:
    """Return optmat's default against the rates of Luna.jl and solprop at 800 nm.

    Returns
    -------
    str
        A Markdown table: the ratios of the instantaneous rates at a field of the
        intensity's amplitude, and the yields of a 30 fs Gaussian pulse of that peak
        intensity. Luna's and solprop's numbers are the fixture's
        (``tools/legacy_ionisation_rates.jl`` and ``.py``); optmat's are computed.
    """
    legacy = json.loads(LEGACY.read_text(encoding="utf-8"))
    luna, solprop = legacy["luna"]["gases"], legacy["solprop"]["gases"]
    peaks = legacy["pulses"]["peak_intensity_w_cm2"]
    rows = []
    for gas in ATOMS:
        model = ionisation_rate_model(gas, legacy["wavelength"])
        for intensity in LEGACY_INTENSITIES:
            i = _legacy_index(legacy["intensity_w_cm2"], intensity)
            field = legacy["field_v_per_m"][i]
            ours = float(model.instantaneous_rate(field))
            theirs = luna[gas]["rate"][i]
            mine = solprop[gas]["rate"][i]
            j = _legacy_index(peaks, intensity)
            yields = [
                luna[gas]["pulses"]["30fs"]["yield"][j],
                solprop[gas]["pulses"]["30fs"]["yield"][j],
                validation.pulse_yield(model, intensity * 1e4, 30e-15, carrier=True),
            ]
            refused = theirs is None
            rows.append(
                [
                    gas,
                    _wcm2(intensity * 1e4),
                    f"{float(model.keldysh(field)):.2f}",
                    "" if refused else f"{mine / theirs:.2f}",
                    "" if refused else f"{ours / theirs:.2f}",
                    f"{ours / mine:.2f}",
                    *("" if y is None else f"{y:.2g}" for y in yields),
                ]
            )
    return _table(
        [
            "Gas",
            "I (W/cm2)",
            "gamma",
            "rate: solprop / Luna",
            "optmat / Luna",
            "optmat / solprop",
            "30 fs yield: Luna",
            "solprop",
            "optmat",
        ],
        rows,
    )


HYDROGEN_FORMS = ("static", "frequency", "popruzhenko")
HYDROGEN_BINS = ((0.0, 1.0), (1.0, 2.0), (2.0, math.inf))


def _binned(ratio: np.ndarray, gamma: np.ndarray) -> list[str]:
    cells = []
    for low, high in HYDROGEN_BINS:
        inside = (gamma >= low) & (gamma < high)
        cells.append(
            f"{float(np.exp(np.mean(np.log(ratio[inside])))):.2f}"
            if inside.any()
            else ""
        )
    return cells


def hydrogen_table() -> str:
    """Return the hydrogen model against the one-electron calculations of hydrogen.

    Returns
    -------
    str
        A Markdown table: for each calculation and Coulomb factor, the geometric mean
        of model over calculation where the yield is below 0.75, in three ranges of
        the Keldysh parameter (:func:`~optmat.gas.ionisation_validation.hydrogen_yields`
        and ``hydrogen_spectra``).
    """
    by_form = {
        form: validation.hydrogen_yields(form) + validation.hydrogen_spectra(form)
        for form in HYDROGEN_FORMS
    }
    rows = []
    for k, first in enumerate(by_form["static"]):
        if math.isnan(first.wavelength):
            field = float(
                datasets_by_key()[first.reference].parameters["intensity_w_cm2"]
            )
            amplitude = validation.amplitude(np.full(first.x.shape, field * 1e4))
            wavelengths = 2 * math.pi * C / (first.x * AU_ANGULAR_FREQUENCY)
            gamma = np.array(
                [
                    float(validation.hydrogen_model(w).keldysh(a))
                    for w, a in zip(wavelengths, amplitude, strict=True)
                ]
            )
            where = f"{first.x.min():.2f} to {first.x.max():.2f} au"
        else:
            gamma = validation.hydrogen_model(first.wavelength).keldysh(
                validation.amplitude(first.x)
            )
            where = f"{first.wavelength * 1e9:.0f} nm"
        below = first.expected < 0.75
        cells = []
        for form in HYDROGEN_FORMS:
            cells += _binned(by_form[form][k].ratio[below], gamma[below])
        rows.append([f"`{first.reference}`", where, *cells])
    headings = [
        f"{form}: gamma {low:g} to {high:g}".replace(" to inf", " and above")
        for form in HYDROGEN_FORMS
        for low, high in HYDROGEN_BINS
    ]
    return _table(["Calculation", "Wavelength or frequency", *headings], rows)


def datasets_by_key() -> dict:
    """Return the curated datasets by key."""
    return {d.key: d for d in datasets()}


def _measured_tool():  # noqa: ANN202 - a module loaded from its path
    """Return ``tools/legacy_ionisation_measured.py`` as a module."""
    name = "legacy_ionisation_measured"
    if name not in sys.modules:
        path = Path(__file__).with_name(f"{name}.py")
        spec = importlib.util.spec_from_file_location(name, path)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


def legacy_tests_table() -> str:
    """Return Luna's, solprop's and optmat's rates through the tests at 800 nm.

    Returns
    -------
    str
        A Markdown table: for each absolute measurement at 800 nm, the factor on its
        intensities that fits each code best (1 is a match; below 1, the code's rate
        is high), or the code's rate over the measured one; and the rate over Parker
        et al.'s two-electron helium. The codes' rates are the fixture's, interpolated
        (``tools/legacy_ionisation_measured.py``, TC-815, TC-829).
    """
    tool = _measured_tool()
    legacy = json.loads(LEGACY.read_text(encoding="utf-8"))
    codes = (*tool.CODES, "optmat")
    rows = []
    for gas in tool.GASES:
        cells = []
        for code in codes:
            with tool.substituted(legacy, code, gas):
                cells.append(validation.intensity_scale_test(gas))
        label = "Wallace et al. 2016: intensity factor (calibration "
        label += f"{100 * cells[0].calibration:.1f} %)"
        rows.append([label, gas, *(f"{c.scale:.2f}" for c in cells)])
    for key in tool.phase_shift_references():
        dataset = next(d for d in datasets() if d.key == key)
        cells = []
        for code in codes:
            with tool.substituted(legacy, code, dataset.material):
                cells.append(validation.phase_shift_test(key).scale)
        rows.append(
            [
                f"`{key}`: intensity factor",
                dataset.material,
                *(f"{c:.2f}" for c in cells),
            ]
        )
    for gas in tool.GASES:
        cells = []
        for code in codes:
            with tool.substituted(legacy, code, gas):
                ratio = validation.measured_yields(gas).ratio
            cells.append(f"{ratio.min():.2f} to {ratio.max():.2f}")
        rows.append(["Wahlstrand et al. 2018: rate over measured", gas, *cells])
    cells = []
    for code in (*tool.CODES, "optmat_800nm", "optmat"):
        with tool.substituted(legacy, code, "He"):
            ratio = tool.parker_ratios()
        cells.append(f"{ratio.min():.2f} to {ratio.max():.2f}")
    rows.append(
        [
            "Parker et al. 2007 (780 nm, two electrons, to 1e15 W/cm2): rate over "
            "theirs, the codes at 800 nm",
            "He",
            cells[0],
            cells[1],
            f"{cells[2]} (at 780 nm: {cells[3]})",
        ]
    )
    return _table(["Test", "Gas", "Luna", "solprop", "optmat"], rows)


def helium_two_electron_table() -> str:
    """Return helium at 390 nm: two electrons, haCC, one electron and Parker et al.

    Returns
    -------
    str
        A Markdown table: per calculation and intensity, the norm-loss 1000 au after the
        pulse, over one electron at the same ionisation energy, and Parker et al.'s
        rate integrated over the pulse.
    """
    by_key = {d.key: d for d in datasets()}
    rows = []
    for ratio in validation.helium_two_electron_ratios():
        for i, intensity in enumerate(ratio.intensity):
            parker = validation.parker_pulse_yield(float(intensity), 6e-15, 390)
            rows.append(
                [
                    f"`{ratio.reference.split('/')[1]}`",
                    _wcm2(float(intensity)),
                    f"{ratio.atom[i]:.4g}",
                    f"{ratio.one_electron[i]:.4g}",
                    f"{ratio.ratio[i]:.2f}",
                    f"{parker:.4g}",
                    f"{ratio.atom[i] / parker:.2f}",
                ]
            )
    box = by_key["travers_2026_trecx/twoelectron_normloss_gauss6_390nm_He_box60"]
    intensity, value = box.abscissa_values()[0], box.scaled_values()[0]
    parker = validation.parker_pulse_yield(float(intensity), 6e-15, 390)
    rows.append(
        [
            "`twoelectron_normloss_gauss6_390nm_He_box60`",
            _wcm2(float(intensity)),
            f"{value:.4g}",
            "",
            "",
            f"{parker:.4g}",
            f"{value / parker:.2f}",
        ]
    )
    return _table(
        [
            "Calculation",
            "I (W/cm2)",
            "Norm-loss",
            "One electron, same Ip",
            "Over one electron",
            "Parker et al. over the pulse",
            "Over Parker et al.",
        ],
        rows,
    )


def helium_matched_table() -> str:
    """Return helium's calculations against each other on box-converged measures.

    Returns
    -------
    str
        A Markdown table: per pair of :data:`~optmat.gas.ionisation_validation.
        HELIUM_MATCHED` and intensity, the two yields, their ratio and its rate form;
        and Parker et al.'s rate integrated over the 390 nm pulse against the two
        electrons.
    """
    rows = []
    ratios = validation.helium_matched_ratios()
    for ratio in ratios:
        for i, intensity in enumerate(ratio.intensity):
            rows.append(
                [
                    f"`{ratio.reference.split('/')[1]}`",
                    f"`{ratio.against}`",
                    _wcm2(float(intensity)),
                    f"{ratio.value[i]:.4g}",
                    f"{ratio.against_value[i]:.4g}",
                    f"{ratio.ratio[i]:.3f}",
                    f"{ratio.rate_ratio[i]:.3f}",
                ]
            )
    two = next(r for r in ratios if "twoelectron_above4" in r.reference)
    for i, intensity in enumerate(two.intensity):
        parker = validation.parker_pulse_yield(float(intensity), 6e-15, 390)
        rows.append(
            [
                "Parker et al.'s rate over the pulse",
                f"`{two.reference}`",
                _wcm2(float(intensity)),
                f"{parker:.4g}",
                f"{two.value[i]:.4g}",
                f"{parker / two.value[i]:.3f}",
                f"{math.log1p(-parker) / math.log1p(-two.value[i]):.3f}",
            ]
        )
    return _table(
        [
            "Calculation",
            "Against",
            "I (W/cm2)",
            "Yield",
            "Against's",
            "Ratio",
            "ln(1 - Y) ratio",
        ],
        rows,
    )


def xenon_rmt_table() -> str:
    """Return RMT's xenon against one electron in its field, and the default against it.

    Returns
    -------
    str
        A Markdown table: per intensity label, RMT's ground-state survival, one
        electron's with m = 0 and m = 1, the atom's depletion in m = 0 electrons and
        over the shell's six independent electrons, and the default's over RMT's.
    """
    rmt = validation.xenon_survival()[0]
    default = {float(i): ratio for i, ratio in zip(rmt.x, rmt.ratio, strict=True)}
    rows = []
    for ratio in validation.xenon_rmt_ratios():
        unscaled = ratio.reference.endswith("@unscaled")
        for i, intensity in enumerate(ratio.intensity):
            rows.append(
                [
                    f"{float(intensity) / 1e18:g}",
                    "-0.4448" if unscaled else "-0.444",
                    f"{ratio.atom[i]:.4f}",
                    f"{ratio.m0[i]:.4f}",
                    f"{ratio.m1[i]:.4f}",
                    f"{ratio.per_m0[i]:.2f}",
                    f"{ratio.per_shell[i]:.3f}",
                    "" if unscaled else f"{default[float(intensity)]:.2f}",
                ]
            )
    return _table(
        [
            "Label (1e14 W/cm2)",
            "5p energy (au)",
            "RMT, P",
            "One electron, m = 0",
            "m = 1",
            "ln P / ln P_0",
            "ln P / (2 ln P_0 + 4 ln P_1)",
            "Default / RMT (-ln P)",
        ],
        rows,
    )


def xenon_measured_table() -> str:
    """Return the default against Hartmann et al.'s measured neutral abundance.

    Returns
    -------
    str
        A Markdown table: per peak intensity, the measured abundance, its logarithm and
        the default's at the peak intensity in the measured 4.5 fs pulse.
    """
    comparison = validation.xenon_survival()[1]
    rows = [
        [
            _wcm2(float(intensity)),
            f"{math.exp(-expected):.3f}",
            f"{expected:.3f}",
            f"{predicted:.3f}",
            f"{predicted / expected:.2f}",
        ]
        for intensity, expected, predicted in zip(
            comparison.x, comparison.expected, comparison.predicted, strict=True
        )
    ]
    return _table(
        [
            "I (W/cm2)",
            "Measured, neutral",
            "-ln P",
            "Default, -ln P",
            "Default / measured",
        ],
        rows,
    )


def published_atoms_table() -> str:
    """Return published many-electron calculations over one electron in their pulses.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for ratio in validation.published_atom_ratios():
        dataset = next(d for d in datasets() if d.key == ratio.reference)
        for i, intensity in enumerate(ratio.intensity):
            rows.append(
                [
                    f"`{ratio.reference}`",
                    dataset.material,
                    f"{float(str(dataset.parameters['wavelength_nm'])):.0f}",
                    _wcm2(float(intensity)),
                    f"{ratio.atom[i]:.4g}",
                    f"{ratio.one_electron[i]:.4g}",
                    f"{ratio.ratio[i]:.2f}",
                ]
            )
    return _table(
        [
            "Published",
            "Gas",
            "nm",
            "I (W/cm2)",
            "Published, outside 20 au",
            "One electron (tRecX)",
            "Ratio",
        ],
        rows,
    )


def weak_field_table() -> str:
    """Return the atoms' weak-field coefficients against the model's one electron.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for factor in validation.weak_field_factors():
        model = ionisation_rate_model(factor.gas, 800e-9)
        rows.append(
            [
                factor.gas.value,
                f"{factor.kappa:.4f}",
                f"{factor.atom:.3f}",
                f"{factor.one_electron:.3f}",
                f"{factor.factor:.2f}",
                f"{model.electrons:g}",
            ]
        )
    return _table(
        [
            "Gas",
            "kappa (HF, au)",
            "g, the atom (HF)",
            "g, the model's one electron",
            "atom / one electron at weak field",
            "N_eff (default)",
        ],
        rows,
    )


def saturation_table() -> str:
    """Return the models against Uiterwaal et al.'s saturation intensities.

    Returns
    -------
    str
        A Markdown table: per gas, the Keldysh parameter and the measured saturation
        intensity, and per model the integral of its rate over the pulse there (1 if
        it matches) and its own saturation intensity over the measured one.
    """
    names = ("default", "default_static", "one_electron", "ppt_popruzhenko")
    rows = []
    for gas in validation.SATURATION_INTENSITIES:
        tests = [validation.saturation_test(gas, name) for name in names]
        first = tests[0]
        rows.append(
            [
                gas,
                f"{first.keldysh:.1f}",
                f"{_wcm2(first.measured)} ({100 * first.sigma / first.measured:.0f} %)",
                *(f"{t.exponent:.2g} ({t.scale:.2f})" for t in tests),
            ]
        )
    return _table(
        [
            "Gas",
            "gamma",
            "I_sat measured (W/cm2)",
            *(f"`{name}`: rate x time (I_sat scale)" for name in names),
        ],
        rows,
    )


def _columns(table: str, indices: tuple[int, ...]) -> str:
    """Select columns without changing the snapshot comparison's arithmetic."""
    return "\n".join(
        "| " + " | ".join(line.strip("|").split("|")[i].strip() for i in indices) + " |"
        for line in table.splitlines()
    )


def luna_table() -> str:
    """Return the recorded Luna snapshot comparison for the public manual."""
    return _columns(legacy_table(), (0, 1, 2, 4, 6, 8))


def luna_tests_table() -> str:
    """Return Luna and optmat through the same recorded validation conditions."""
    return _columns(legacy_tests_table(), (0, 1, 2, 4))


TABLES = {
    "ionisation-rate-models": models_table,
    "ionisation-rate-electrons": electrons_table,
    "ionisation-rate-weak-field": weak_field_table,
    "ionisation-rate-atoms": published_atoms_table,
    "ionisation-rate-static": static_table,
    "ionisation-rate-helium": helium_table,
    "ionisation-rate-helium-two-electron": helium_two_electron_table,
    "ionisation-rate-helium-matched": helium_matched_table,
    "ionisation-rate-xenon-rmt": xenon_rmt_table,
    "ionisation-rate-xenon-measured": xenon_measured_table,
    "ionisation-rate-rmatrix": rmatrix_table,
    "ionisation-rate-tdse": tdse_table,
    "ionisation-rate-hydrogen": hydrogen_table,
    "ionisation-rate-tdse42": long_pulse_table,
    "ionisation-rate-saturation": saturation_table,
    "ionisation-rate-wahlstrand": wahlstrand_table,
    "ionisation-rate-wallace": wallace_table,
    "ionisation-rate-wallace-points": wallace_points_table,
    "ionisation-rate-accuracy": accuracy_table,
    "ionisation-rate-legacy": legacy_table,
    "ionisation-rate-legacy-tests": legacy_tests_table,
    "ionisation-rate-phase": phase_shift_table,
    "ionisation-rate-luna": luna_table,
    "ionisation-rate-luna-tests": luna_tests_table,
}


def main() -> None:
    """Print the tables."""
    with warnings.catch_warnings():
        # The yields are taken into saturation, beyond the fitted range of the static
        # correction, where the ionised fraction is 1 whatever the rate; and the
        # literature's forms have no static correction, so their range ends at the
        # barrier-suppression field.
        warnings.simplefilter("ignore", ModelRangeWarning)
        for name, function in TABLES.items():
            print(f"<!-- {name}:start -->\n{function()}\n<!-- {name}:end -->\n")


if __name__ == "__main__":
    main()
