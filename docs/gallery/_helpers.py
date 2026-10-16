# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Shared plotting for the gallery pages (not a gallery page itself)."""

import dataclasses
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from optmat.conditions import ONE_ATM_0C, State
from optmat.constants import AU_ANGULAR_FREQUENCY, AU_POLARIZABILITY, EPS0, C
from optmat.datasets import Dataset, datasets, sources
from optmat.exceptions import ModelRangeWarning
from optmat.gas import reduction
from optmat.gas.comparison import database_models
from optmat.gas.gases import Gas
from optmat.gas.ionisation_rates import ionisation_rate_model
from optmat.gas.models.base import Model
from optmat.gas.models.formulas import published_models
from optmat.gas.models.registry import default_model
from optmat.gas.molecular_ionisation_rates import molecular_rate_model
from optmat.gas.nonlinear import AU_HYPERPOLARIZABILITY, hyperpolarizability_models
from optmat.vocabulary import POLARIZABILITY_QUANTITIES

# A categorical palette, distinct in hue and lightness; datasets past the eighth take
# the next marker shape rather than a generated hue.
COLOURS = [
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#4a3aa7",
    "#e34948",
]
MARKERS = ["o", "s", "D", "^", "v", "P", "X", "<", ">", "h"]
LINE_COLOURS = ["#0b0b0b", "#52514e", "#86b6ef", "#184f95", "#c3c2b7"]


def plottable(gas: str) -> list[Dataset]:
    """Return a gas's datasets that determine the polarizability at a wavelength."""
    return [
        d
        for d in datasets(gas)
        if d.quantity in POLARIZABILITY_QUANTITIES
        and d.use not in {"check", "excluded"}
        and d.abscissa == "wavelength"
    ]


def _quiet(model: Model, wavelength: np.ndarray) -> np.ndarray:
    """Evaluate a model's refractivity at 0 C, 1 atm without range warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        return model.refractivity(wavelength, ONE_ATM_0C)


KIND_MARKERS = {"measurement": "o", "theory": "^", "compilation": "s", "fit": "D"}


def _data(
    top: Axes, bottom: Axes | None, gas: str, reference: Model | None, edge: float
) -> list[float]:
    """Draw every plottable dataset and its deviation; return the deviations (%).

    A source's own values take the source's colour (in the fixed categorical order)
    and a marker for the source's kind; every reproduced value is a small grey open
    marker, labelled once.
    """
    deviations: list[float] = []
    colours: dict[str, str] = {}
    reproduced_labelled = False
    for dataset in plottable(gas):
        wavelength = dataset.vacuum_wavelength()
        shown = np.where(np.isinf(wavelength), edge, wavelength) * 1e9
        dn = reduction.refractivity_at(dataset, ONE_ATM_0C)
        sigma = np.nan_to_num(
            reduction.polarizability_uncertainty(dataset)
            / reduction.polarizability(dataset)
            * dn
        )
        kind = sources()[dataset.source].kind
        if dataset.use == "reproduced":
            style = {
                "marker": ".",
                "color": "#9a9994",
                "markersize": 5,
                "zorder": 1,
                "markerfacecolor": "none",
                "elinewidth": 0.5,
            }
            label = None if reproduced_labelled else "reproduced from elsewhere"
            reproduced_labelled = True
        else:
            if dataset.source not in colours:
                colours[dataset.source] = COLOURS[len(colours) % len(COLOURS)]
                label = dataset.source
            else:
                label = None
            style = {
                "marker": KIND_MARKERS.get(kind, "o"),
                "color": colours[dataset.source],
                "markersize": 5,
                "zorder": 3,
                "elinewidth": 0.8,
            }
            if len(colours) > len(COLOURS):
                style["markerfacecolor"] = "white"
        style["linestyle"] = "none"
        top.errorbar(shown, dn * 1e6, yerr=sigma * 1e6, label=label, **style)
        if bottom is not None and reference is not None:
            ref = _quiet(reference, wavelength)
            deviation = (dn / ref - 1) * 100
            deviations.extend(deviation[np.isfinite(deviation)])
            bottom.errorbar(shown, deviation, yerr=sigma / ref * 100, **style)
    return deviations


def _models(
    top: Axes, bottom: Axes | None, gas: str, reference: Model | None, xlim: tuple
):
    """Draw every published, legacy and comparison model over its stated range."""
    grid = np.geomspace(*xlim, 800)
    models = [m for m in published_models().values() if m.gas == gas]
    models += [m for m in database_models().values() if m.gas == gas]
    for number, model in enumerate(models):
        low, high = model.wavelength_range
        inside = grid[(grid >= low) & (grid <= high)]
        if inside.size == 0:
            continue
        dn = model.refractivity(inside, warn=False)
        style = {
            "linewidth": 1.1,
            "linestyle": {"comparison": ":", "legacy": "--"}.get(model.kind, "-"),
            "color": LINE_COLOURS[number % len(LINE_COLOURS)],
        }
        top.plot(inside * 1e9, dn * 1e6, **style)
        if bottom is not None and reference is not None:
            ref = _quiet(reference, inside)
            bottom.plot(inside * 1e9, (dn / ref - 1) * 100, label=model.name, **style)


def data_overview(
    gas: str, reference: str | None, *, xlim: tuple[float, float]
) -> Figure:
    """Plot every direct dataset of a gas and the published models.

    Top: n - 1 of the real gas at 0 degrees C and 1 atm implied by each value, against
    vacuum wavelength, with the published formulas (solid), Luna's legacy formulas
    (dashed) and refractiveindex.info's (dotted) over their stated ranges. Bottom: the
    deviation of each value and model from a reference model, in per cent. Filled
    markers are a source's own values, open markers values reproduced from elsewhere;
    static values (lambda = inf) are drawn at the right edge.
    """
    reference_model = None if reference is None else published_models()[reference]
    if reference_model is None:
        fig, top = plt.subplots(figsize=(11, 5), layout="constrained")
        bottom = None
    else:
        fig, (top, bottom) = plt.subplots(
            2,
            1,
            sharex=True,
            figsize=(11, 8.5),
            layout="constrained",
            height_ratios=(3, 2),
        )
    edge = xlim[1] * 0.9
    deviations = _data(top, bottom, gas, reference_model, edge)
    _models(top, bottom, gas, reference_model, xlim)
    top.set_xscale("log")
    top.set_xlim(*(np.array(xlim) * 1e9))
    top.set_ylabel(r"$(n - 1) \times 10^6$ at 0 °C, 1 atm (real gas)")
    top.set_title(f"{gas}: curated data and published models")
    top.legend(
        fontsize=6.5,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        title="data (o measured, ^ calculated, s compiled)",
        title_fontsize=7,
    )
    axes = [top]
    if bottom is not None:
        bottom.legend(
            fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1.0), title="models"
        )
        bottom.set_ylabel(f"deviation from\n{reference} (%)")
        if deviations:
            # Robust limits: most of the data, not the few points near a resonance or
            # from an outlying source, which fall outside the panel.
            low, high = np.percentile(deviations, [3, 97])
            pad = max(0.25 * (high - low), 0.2)
            bottom.set_ylim(low - pad, high + pad)
        axes.append(bottom)
    axes[-1].set_xlabel("vacuum wavelength (nm); static values at the right edge")
    for axis in axes:
        axis.grid(which="both", alpha=0.25)
    return fig


def _fit_data(top: Axes, bottom: Axes, model, config, xlim: tuple) -> None:
    """Draw direct data with the fitter's conventions, windows and uncertainties.

    Configured data are put on the refractive convention by adding the difference
    between the full model and its prediction of that observation (TC-729). Thus the
    percentage deviations and the covariance band use the same reference. Bottom
    residuals have the opposite sign to the optimiser: positive means data above fit.
    """
    from optmat.gas.fitting import (  # noqa: PLC0415 - gallery only
        _model_values,
        _observation,
    )

    fitted = {o.key: o for o in config.fitted}
    validators = {o.key: o for o in config.validators}
    scales = model.record["scales"]
    n_fitted = len(model.strengths) - model.n_fixed
    strengths = np.asarray(model.strengths)[:n_fitted]
    energies = np.asarray(model.energies)[:n_fitted]
    edge = xlim[1] * 0.9
    colours: dict[str, str] = {}
    for dataset in plottable(config.gas.value):
        wavelength = dataset.vacuum_wavelength()
        if dataset.key in fitted:
            role = "fitted"
            options = fitted[dataset.key]
            scale = float(scales[dataset.key])
        elif dataset.key in validators:
            role = "validator"
            # Like fitting.validation, permit validators without stated uncertainty.
            # The tiny allowance keeps _observation valid without a visible error bar.
            options = dataclasses.replace(
                validators[dataset.key],
                extra_relative=max(validators[dataset.key].extra_relative, 1e-300),
            )
            scale = 1.0
        else:
            role = "other"
            options = None
            scale = 1.0
        if options is not None:
            # The fitter's reference temperature brings datasets with a temperature
            # option to it, as in the fit (TC-730).
            obs = _observation(dataset, options, config.reference_temperature)
            wavelength = np.divide(
                2 * np.pi * C / AU_ANGULAR_FREQUENCY,
                obs.omega,
                out=np.full_like(obs.omega, np.inf),
                where=obs.omega != 0,
            )
            reference = model.polarizability_au(obs.omega)
            prediction = _model_values(strengths, energies, obs, config)
            difference = obs.alpha / scale - prediction
            sigma = obs.sigma / scale
        else:
            reference = model.polarizability(wavelength, warn=False) / AU_POLARIZABILITY
            difference = (
                reduction.polarizability(dataset) / AU_POLARIZABILITY - reference
            )
            sigma = (
                np.nan_to_num(reduction.polarizability_uncertainty(dataset))
                / AU_POLARIZABILITY
            )
        shown = np.where(np.isinf(wavelength), edge, wavelength) * 1e9
        style: dict = {"marker": ".", "color": "#9a9994", "markersize": 4, "zorder": 1}
        label = None
        if role != "other":
            colour = colours.setdefault(dataset.key, COLOURS[len(colours) % 8])
            style = {
                "marker": "o" if role == "fitted" else "s",
                "color": colour,
                "markersize": 5,
                "zorder": 3,
                "markerfacecolor": colour if role == "fitted" else "none",
            }
            label = f"{dataset.key} ({role})"
        top.errorbar(
            shown,
            difference / reference * 100,
            yerr=sigma / np.abs(reference) * 100,
            linestyle="none",
            elinewidth=0.7,
            label=label,
            **style,
        )
        if role == "fitted":
            bottom.plot(
                shown, difference / sigma, linestyle="none", label=dataset.key, **style
            )


def _fit_models(axis: Axes, gas: str, grid: np.ndarray, alpha: np.ndarray) -> None:
    """Draw every published, legacy and comparison model's deviation from the fit."""
    models = [m for m in published_models().values() if m.gas == gas]
    models += [m for m in database_models().values() if m.gas == gas]
    for number, other in enumerate(models):
        low, high = other.wavelength_range
        inside = (grid >= low) & (grid <= high)
        if not np.any(inside):
            continue
        ratio = other.polarizability(grid[inside], warn=False) / (
            alpha[inside] * AU_POLARIZABILITY
        )
        axis.plot(
            grid[inside] * 1e9,
            (ratio - 1) * 100,
            label=other.name,
            linestyle={"comparison": ":", "legacy": "--"}.get(other.kind, "-"),
            color=LINE_COLOURS[number % len(LINE_COLOURS)],
            linewidth=1.1,
        )


def fit_overview(name: str, *, xlim: tuple[float, float]) -> Figure:
    """Plot a stored fit against its data, validators and the published models.

    Top: direct data's percentage deviations on the refractive convention, with the
    fit's one-sigma band (shaded). Configured datasets use the fitter's wavelength
    windows, corrections and uncertainties, with fitted scales removed; other data
    use their stated uncertainties. Middle: published, legacy and comparison models.
    Bottom: the fitted direct observations' residuals in the uncertainty the fit used,
    positive for data above the prediction. Moments, ratios and priors are not shown.
    """
    from optmat.gas.fitting import load_config  # noqa: PLC0415 - gallery only
    from optmat.gas.models.fitted import fitted_models  # noqa: PLC0415
    from optmat.gas.models.oscillator import angular_frequency_au  # noqa: PLC0415

    model = fitted_models()[f"fit/{name}"]
    config = load_config(name)
    fig, (top, middle, bottom) = plt.subplots(
        3,
        1,
        sharex=True,
        figsize=(11, 11),
        layout="constrained",
        height_ratios=(3, 2, 2),
    )
    grid = np.geomspace(*xlim, 800)
    omega = angular_frequency_au(grid)
    alpha = model.polarizability_au(omega)
    band = model.polarizability_sigma_au(omega) / alpha * 100
    for axis in (top, middle):
        axis.fill_between(
            grid * 1e9, -band, band, color="#cde2fb", label="fit, 1 sigma"
        )
        axis.axhline(0.0, color="#2a78d6", linewidth=1.0)
    _fit_data(top, bottom, model, config, xlim)
    _fit_models(middle, config.gas.value, grid, alpha)
    top.set_xscale("log")
    top.set_xlim(*(np.array(xlim) * 1e9))
    top.set_ylabel("data / fit - 1 (%)")
    middle.set_ylabel("model / fit - 1 (%)")
    bottom.set_ylabel("fitted residuals (sigma)")
    bottom.set_xlabel("vacuum wavelength (nm); static values at the right edge")
    chi2 = model.reduced_chi2
    top.set_title(f"{config.gas.value}: fit {name} (reduced chi2 {chi2:.2f})")
    top.legend(fontsize=6.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    middle.legend(fontsize=6.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    for axis, limit in ((top, 1.0), (middle, 2.0), (bottom, 6.0)):
        axis.set_ylim(-limit, limit)
        axis.grid(which="both", alpha=0.25)
    return fig


def solprop_defaults() -> dict[str, Model]:
    """Return solprop's current default formula for each gas fitted here.

    The mapping is the one ``tests/test_solprop_parity.py`` checks (TC-100);
    solprop reads Ermolov's helium coefficients at 1 atm, so that model is re-referred.

    Returns
    -------
    dict of str to Model
        By gas name.
    """
    formulas = published_models()
    ermolov = dataclasses.replace(
        formulas["luna_physdata/jct_He"],
        name="Ermolov at 1 atm (solprop)",
        state=ONE_ATM_0C,
    )
    return {
        "He": ermolov,
        "Ne": formulas["borzsonyi_2008/sellmeier_Ne"],
        "Ar": formulas["luna_physdata/kingston_fit_Ar"],
        "Kr": formulas["luna_physdata/kingston_fit_Kr_1atm"],
        "Xe": formulas["borzsonyi_2008/sellmeier_Xe"],
        "N2": formulas["borzsonyi_2008/sellmeier_N2"],
        "H2": formulas["peck_1977/eq7"],
        # solprop applies the hydrogen formula to deuterium (Orr 1936).
        "D2": formulas["peck_1977/eq7"],
        "O2": formulas["kren_2011/kren_O2"],
        "CH4": formulas["he_2021/CH4"],
        "SF6": formulas["he_2021/SF6"],
        "N2O": formulas["he_2021/N2O"],
    }


def solprop_gamma() -> dict[str, float]:
    """Return solprop's current third-order coefficient of each gas as gamma_T (au).

    solprop uses one frozen ``gamma3 = gamma_T / (6 eps0)`` per gas for every process
    and wavelength; the values are the fixture ``tests/data/solprop_gamma3.json``
    (``tools/solprop_gamma3.py``, run in solprop's environment).

    Returns
    -------
    dict of str to float
        By gas name, atomic units.
    """
    path = Path(__file__).parents[2] / "tests" / "data" / "solprop_gamma3.json"
    record = json.loads(path.read_text())["gamma3"]
    return {
        gas: 6 * EPS0 * entry["gamma3"] / AU_HYPERPOLARIZABILITY
        for gas, entry in record.items()
    }


def group_velocity_dispersion(
    model: Model, wavelength: np.ndarray, state: State
) -> np.ndarray:
    """Return beta_2 in fs^2/m by a centred second difference in angular frequency.

    Parameters
    ----------
    model : Model
        The dispersion model.
    wavelength : numpy.ndarray
        Vacuum wavelengths, m.
    state : State
        The gas state.

    Returns
    -------
    numpy.ndarray
        beta_2, fs^2/m.
    """
    omega = 2 * np.pi * C / wavelength
    step = omega * 1e-3
    # The vacuum part of beta = omega n / c is linear in omega and drops out of the
    # second difference; leaving it out avoids cancelling it in floating point.
    beta = [
        (omega + j * step) / C * _state_refractivity(model, omega + j * step, state)
        for j in (-1, 0, 1)
    ]
    return (beta[0] - 2 * beta[1] + beta[2]) / step**2 * 1e30


def _state_refractivity(model: Model, omega: np.ndarray, state: State) -> np.ndarray:
    """Return n - 1 at a state for angular frequencies, without range warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        return model.refractivity(2 * np.pi * C / omega, state)


# Luna.jl's material parameters (the page *Coming from Luna*): the fixture
# tests/data/luna_master.json, written by tools/luna_reference.jl in Luna's environment.
LUNA_STATE = State(1e5, 293.15)
LUNA_MOLECULES = ("H2", "D2")


def luna() -> dict:
    """Return the fixture of Luna's values (``tools/luna_reference.jl``)."""
    path = Path(__file__).parents[2] / "tests" / "data" / "luna_master.json"
    return json.loads(path.read_text())


def _luna_at(entry: dict, key: str, wavelength: np.ndarray) -> np.ndarray:
    grid = np.array(entry["wavelength"])
    index = [int(np.argmin(np.abs(grid - w))) for w in np.atleast_1d(wavelength)]
    return np.array(entry[key])[index]


def luna_refractivity_ratio(gas: str, wavelength: float, pressure_bar: float) -> float:
    """Return Luna's n - 1 over optmat's at 20 C (800 nm unless at 1 bar)."""
    entry = luna()["gases"][gas]
    state = State(pressure_bar * 1e5, LUNA_STATE.temperature)
    if pressure_bar == 1.0:
        theirs = float(_luna_at(entry, "refractivity", wavelength)[0])
    else:
        theirs = entry["refractivity_800nm_by_pressure"][f"{pressure_bar}"]
    return theirs / float(default_model(gas).refractivity(wavelength, state))


def luna_gvd_ratio(gas: str, wavelength: np.ndarray) -> np.ndarray:
    """Return Luna's group-velocity dispersion over optmat's at 1 bar and 20 C."""
    theirs = _luna_at(luna()["gases"][gas], "beta2", wavelength) * 1e30
    ours = group_velocity_dispersion(default_model(gas), wavelength, LUNA_STATE)
    return theirs / ours


def luna_kerr_ratio(gas: str, wavelength: np.ndarray, process: str = "kerr"):
    """Return Luna's frozen gamma3 over optmat's coefficient for a process."""
    model = hyperpolarizability_models()[Gas(gas)]
    return luna()["gases"][gas]["gamma3"] / model.gamma3(process, wavelength)


def luna_rates(gas: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return intensity (W/cm^2), Luna's and optmat's cycle-averaged rates at 800 nm.

    Luna's is its instantaneous rate averaged over a cycle, as a propagation that
    resolves the carrier applies it; the intensities stop at the largest field of
    optmat's default model.
    """
    record = luna()["gases"][gas]["rate"]
    intensity = np.array(record["intensity_w_cm2"])
    field = np.sqrt(2 * intensity * 1e4 / (C * EPS0))
    if gas in LUNA_MOLECULES:
        model = molecular_rate_model(gas, record["wavelength"])
    else:
        model = ionisation_rate_model(gas, record["wavelength"])
    keep = field <= model.maximum_field
    ours = model.cycle_averaged_rate(field[keep])
    return intensity[keep], np.array(record["cycle_averaged"])[keep], ours
