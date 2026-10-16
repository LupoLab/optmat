# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Test the ionisation-rate models against the curated references (TC-660).

Nothing in :mod:`optmat.gas.ionisation_rates` is fitted to a rate or a yield, except
helium's effective number of electrons, which is set on Parker et al.'s static rates
(TC-658). Every reference is therefore a validator, and this module predicts each one
from a model and returns the comparison. The references, in the order of their rank
(TC-439):

- **static rates** from all-electron calculations (helium) and from exact solutions
  for one electron in a model potential (argon);
- **rates in a laser field** from the two-electron calculations of helium, as the
  fits their authors publish (Parker et al. 2007, eq. 5);
- **measured yields**: Wahlstrand et al.'s (2018) absolute yields of argon, krypton
  and xenon for a 42 fs pulse, as their published fits, and Wallace et al.'s (2016)
  yields for a 5.5 fs pulse on a free scale, whose intensities are fixed by atomic
  hydrogen;
- **yields from time-dependent Schrödinger calculations** for one electron in Tong
  and Lin's model potentials (``travers_2026_trecx``), to compare with the
  ``one_electron`` model.

A yield is the ionised fraction after a pulse, :math:`1 - \exp(-\int W\,dt)`
(:func:`ionised_fraction`), with the rate taken along the pulse in one of two ways
(:func:`pulse_yield`): the instantaneous rate on the carrier-resolved field, or the
cycle-averaged rate on the envelope.

SI throughout: intensities in W/m\ :sup:`2` (cycle averaged, :math:`I = \varepsilon_0
c F^2/2` for an amplitude :math:`F`), durations in s.
"""

import math
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from jaxtyping import Float
from numpy.typing import ArrayLike

from optmat.constants import (
    AU_ANGULAR_FREQUENCY,
    AU_ENERGY,
    AU_FIELD,
    E_CHARGE,
    EPS0,
    M_E,
    C,
)
from optmat.datasets import Dataset, datasets, sources
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas.gases import Gas, as_gas
from optmat.gas.ionisation_rates import IonisationRateModel, ionisation_rate_model
from optmat.gas.molecular_ionisation_rates import (
    MolecularRateModel,
    molecular_rate_model,
)

__all__ = [
    "HACC_REFERENCES",
    "HACC_YIELDS",
    "HELIUM_MATCHED",
    "HELIUM_TWO_ELECTRON",
    "HYDROGEN_ALPHA",
    "HYDROGEN_SPECTRA",
    "HYDROGEN_YIELDS",
    "IP_MATCHED_YIELDS",
    "LONG_PULSE_YIELD",
    "PUBLISHED_ATOMS",
    "RMATRIX_RATES",
    "SATURATION_INTENSITIES",
    "STATIC_REFERENCES",
    "TDSE_YIELDS",
    "WALLACE_TDSE",
    "WEAK_FIELD_ATOMS",
    "XENON_RMT",
    "XENON_SURVIVAL",
    "AtomRatio",
    "PhaseShiftTest",
    "PointsTest",
    "RateComparison",
    "SaturationTest",
    "ScaleTest",
    "SurvivalRatio",
    "WeakFieldFactor",
    "YieldRatio",
    "amplitude",
    "circular_yield",
    "cycle_pulse",
    "focal_average",
    "gaussian_pulse",
    "hacc_ratios",
    "helium_laser_rates",
    "helium_matched_ratios",
    "helium_two_electron_ratios",
    "hydrogen_model",
    "hydrogen_spectra",
    "hydrogen_yields",
    "intensity_scale_test",
    "ionised_fraction",
    "long_pulse_yield",
    "measured_points_test",
    "measured_yields",
    "parker_pulse_yield",
    "parker_rate_au",
    "phase_shift_test",
    "plasma_phase_shift",
    "published_atom_ratios",
    "pulse_yield",
    "rmatrix_rates",
    "saturation_test",
    "static_rates",
    "tdse_intensity_scale_test",
    "tdse_yields",
    "wallace_yield",
    "weak_field_factors",
    "xenon_rmt_ratios",
    "xenon_survival",
]

_SAMPLES_PER_CYCLE = 96
"""Samples of the carrier per optical cycle in a pulse's yield.

The instantaneous rate is sharply peaked at the crests of the field: with a
tunnelling exponent of 15 its width is a tenth of a cycle. 96 samples put ten across
it, and the trapezoid rule on a smooth periodic integrand then converges
geometrically; the yield changes by less than 1e-4 when the number is doubled.
"""

_FIELD_LIMIT = 2.0
"""Largest field compared at, in barrier-suppression fields.

Tong and Lin's correction holds to twice the barrier-suppression field. The models
without a correction hold to less, and are compared over the same fields so that the
comparisons can be set side by side; beyond their range a pulse's yield is saturated.
"""

_PULSE_HALF_WIDTH = 3.0
"""Half of the time window of a Gaussian pulse, in FWHM of its intensity.

At three FWHM the intensity has fallen by exp(-100): nothing ionises beyond it.
"""

STATIC_REFERENCES: Mapping[str, tuple[tuple[str, str, str], ...]] = {
    "He": (
        ("parker_2009/table1_present", "default", "Parker et al. 2009, two electrons"),
        ("scrinzi_1999/table1_He", "default", "Scrinzi et al. 1999, two electrons"),
        ("majety_2015/table1_hacc5_He", "default", "Majety & Scrinzi 2015, haCC(5)"),
        ("jagau_2016/tableS3_fci_He", "default", "Jagau 2016, FCI"),
    ),
    "Ne": (
        ("jagau_2016/tableS4_ccsd_Ne", "default", "Jagau 2016, CCSD, all electrons"),
        ("bahl_2016/rate_Ne", "one_electron", "Bahl et al. 2016, one electron"),
    ),
    "Ar": (
        ("jagau_2016/tableS5_ccsd_Ar", "default", "Jagau 2016, CCSD, all electrons"),
        (
            "zhang_2014/tableIX_tonglin_Ar",
            "one_electron",
            "Zhang et al. 2014, one electron, Tong-Lin potential",
        ),
        (
            "zhang_2014/tableIX_gsz_Ar",
            "one_electron",
            "Zhang et al. 2014, one electron, GSZ potential",
        ),
        (
            "zhang_2014/tableIX_muller_Ar",
            "one_electron",
            "Zhang et al. 2014, one electron, Muller potential",
        ),
        (
            "travers_2026_trecx/static_Ar_m0",
            "one_electron",
            "tRecX, one electron, Tong-Lin potential",
        ),
        ("bahl_2016/rate_Ar", "one_electron", "Bahl et al. 2016, one electron"),
    ),
    "Kr": (
        ("bahl_2016/rate_Kr", "one_electron", "Bahl et al. 2016, one electron"),
        (
            "travers_2026_trecx/static_Kr_m0",
            "one_electron",
            "tRecX, one electron, Zhang et al.'s GSZ potential",
        ),
        (
            "santra_2006/total_Kr",
            "default",
            "Santra et al. 2006, spin-orbit, all 4p electrons",
        ),
    ),
    "Xe": (
        ("bahl_2016/rate_Xe", "one_electron", "Bahl et al. 2016, one electron"),
        (
            "travers_2026_trecx/static_Xe_m0",
            "one_electron",
            "tRecX, one electron, Zhang et al.'s Tong-Lin-form potential",
        ),
    ),
}
"""The static references of each gas: dataset, the model to compare, and a label.

An all-electron rate tests the default; the rate of one electron in a model
potential tests the ``one_electron`` model, which leaves out the effective number of
electrons and the fine structure of the ion.
"""


@dataclass(frozen=True)
class RateComparison:
    """A model's prediction of one reference.

    Attributes
    ----------
    gas : Gas or str
        The gas, or ``"H"`` for atomic hydrogen (:func:`hydrogen_yields`).
    reference : str
        The dataset or formula compared with, ``source/id``.
    label : str
        What the reference is.
    model : str
        The name of the rate model.
    quantity : str
        ``static rate`` (atomic units, against the field in V/m), ``rate``
        (atomic units, cycle averaged, against the intensity in W/m^2), ``yield``
        (against the peak intensity in W/m^2) or ``depletion`` (-ln P, P
        the ground-state survival, against the peak intensity in W/m^2).
    wavelength : float
        Wavelength of the field, m; ``math.inf`` for a static field.
    x : numpy.ndarray
        The field or the intensity.
    expected : numpy.ndarray
        The reference's values.
    predicted : numpy.ndarray
        The model's.
    """

    gas: Gas | str
    reference: str
    label: str
    model: str
    quantity: str
    wavelength: float
    x: Float[np.ndarray, " n"]
    expected: Float[np.ndarray, " n"]
    predicted: Float[np.ndarray, " n"]

    @property
    def ratio(self) -> Float[np.ndarray, " n"]:
        """The model's value over the reference's.

        Returns
        -------
        numpy.ndarray
            ``predicted / expected``.
        """
        return self.predicted / self.expected


def _evaluate(
    model: IonisationRateModel | MolecularRateModel,
    field_strength: Float[np.ndarray, " *shape"],
    *,
    cycle_averaged: bool,
) -> Float[np.ndarray, " *shape"]:
    """Return a model's rate without its warning beyond the fitted range.

    Parameters
    ----------
    model : IonisationRateModel or MolecularRateModel
        The rate.
    field_strength : numpy.ndarray
        Field, V/m.
    cycle_averaged : bool
        The rate against the amplitude, or the instantaneous rate.

    Returns
    -------
    numpy.ndarray
        The rate, 1/s.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        if cycle_averaged:
            return model.cycle_averaged_rate(field_strength)
        return model.instantaneous_rate(field_strength)


def amplitude(intensity: ArrayLike) -> Float[np.ndarray, " *shape"]:
    r"""Return the amplitude of the field of a cycle-averaged intensity.

    Parameters
    ----------
    intensity : array_like
        Intensity, W/m\ :sup:`2`.

    Returns
    -------
    numpy.ndarray
        :math:`F = (2I/\varepsilon_0 c)^{1/2}`, V/m.

    Examples
    --------
    1e14 W/cm\ :sup:`2`:

    >>> f"{float(amplitude(1e18)):.4e}"
    '2.7449e+10'
    """
    return np.asarray(np.sqrt(2 * np.asarray(intensity, dtype=float) / (EPS0 * C)))


def gaussian_pulse(
    peak_intensity: float,
    fwhm: float,
    wavelength: float,
    *,
    cep: float = 0.0,
) -> tuple[Float[np.ndarray, " n"], Float[np.ndarray, " n"], Float[np.ndarray, " n"]]:
    r"""Return the times, the envelope and the field of a Gaussian pulse.

    The vector potential is :math:`A(t) = -(F_0/\omega)\, e(t) \sin(\omega t +
    \varphi)` with the envelope :math:`e(t) = \exp(-2 \ln 2\, t^2/\tau^2)`, and the
    field is :math:`-dA/dt`: the pulse then has no static component however short it
    is. :math:`\tau` is the full width at half maximum of the intensity
    :math:`e^2`. For :math:`\varphi = 0` the field has a crest at the peak of the
    envelope. This is the pulse of tRecX (``gauss``; Scrinzi 2022) and of the
    calculations of ``travers_2026_trecx``.

    Parameters
    ----------
    peak_intensity : float
        Peak of the cycle-averaged intensity, W/m\ :sup:`2`.
    fwhm : float
        Full width at half maximum of the intensity, s.
    wavelength : float
        Carrier wavelength, m.
    cep : float, optional
        Carrier-envelope phase :math:`\varphi`, rad.

    Returns
    -------
    tuple of numpy.ndarray
        The times (s), the amplitude of the envelope (V/m) and the field (V/m), on a
        uniform grid of 96 points per cycle over three FWHM
        either side of the peak.

    Examples
    --------
    >>> times, envelope, field = gaussian_pulse(1e18, 6e-15, 800e-9)
    >>> f"{float(np.max(np.abs(field)) / np.max(envelope)):.4f}"
    '1.0000'
    """
    omega = 2 * math.pi * C / wavelength
    period = wavelength / C
    half = _PULSE_HALF_WIDTH * fwhm
    count = 2 * math.ceil(half / period * _SAMPLES_PER_CYCLE) + 1
    times = np.linspace(-half, half, count)
    peak = float(amplitude(peak_intensity))
    shape = np.exp(-2 * math.log(2) * times**2 / fwhm**2)
    slope = -4 * math.log(2) * times / fwhm**2 * shape
    phase = omega * times + cep
    # -dA/dt of A = -(F0 / omega) e(t) sin(phase).
    field = peak * (shape * np.cos(phase) + slope / omega * np.sin(phase))
    return times, peak * shape, field


def cycle_pulse(
    peak_intensity: float,
    wavelength: float,
    ramp_cycles: float,
    plateau_cycles: float = 0.0,
    *,
    cep: float = 0.0,
) -> tuple[Float[np.ndarray, " n"], Float[np.ndarray, " n"], Float[np.ndarray, " n"]]:
    r"""Return the times, the envelope and the field of a pulse counted in cycles.

    The vector potential is :math:`A(t) = -(F_0/\omega)\, e(t) \sin(\omega t +
    \varphi)`, the field :math:`-dA/dt`, with the envelope rising as
    :math:`\sin^2` over ``ramp_cycles``, flat over ``plateau_cycles`` and falling as
    :math:`\sin^2` over ``ramp_cycles``: with no plateau a :math:`\sin^2` pulse of
    :math:`2 \times` ``ramp_cycles`` cycles in all, with one the flat-top pulse of
    Bauer and Deng (2018; 2-cycle ramps and a 10-cycle plateau). The time origin is
    the pulse's centre.

    Parameters
    ----------
    peak_intensity : float
        Peak of the cycle-averaged intensity, W/m\ :sup:`2`.
    wavelength : float
        Carrier wavelength, m.
    ramp_cycles : float
        Cycles of each :math:`\sin^2` ramp, > 0.
    plateau_cycles : float, optional
        Cycles of the plateau, >= 0.
    cep : float, optional
        Carrier-envelope phase :math:`\varphi`, rad, relative to the centre.

    Returns
    -------
    tuple of numpy.ndarray
        The times (s), the amplitude of the envelope (V/m) and the field (V/m), on a
        uniform grid of 96 points per cycle.

    Examples
    --------
    A 2-10-2 flat top at 800 nm lasts 14 cycles, 37.4 fs:

    >>> times, envelope, field = cycle_pulse(1e18, 800e-9, 2, 10)
    >>> f"{(times[-1] - times[0]) * 1e15:.1f}"
    '37.4'
    """
    omega = 2 * math.pi * C / wavelength
    period = wavelength / C
    total = 2 * ramp_cycles + plateau_cycles
    count = math.ceil(total * _SAMPLES_PER_CYCLE) + 1
    times = np.linspace(-total / 2, total / 2, count) * period
    ramp = ramp_cycles * period
    edge = np.clip((total / 2 * period - np.abs(times)) / ramp, 0.0, 1.0)
    shape = np.sin(math.pi / 2 * edge) ** 2
    inside = (edge > 0) & (edge < 1)
    # de/dt: rising before the plateau, falling after it.
    slope = np.where(
        inside, np.sin(math.pi * edge) * math.pi / (2 * ramp) * -np.sign(times), 0.0
    )
    peak = float(amplitude(peak_intensity))
    phase = omega * times + cep
    field = peak * (shape * np.cos(phase) + slope / omega * np.sin(phase))
    return times, peak * shape, field


def ionised_fraction(
    rate: Float[np.ndarray, " n"], times: Float[np.ndarray, " n"]
) -> float:
    r"""Return the fraction ionised by the end of a pulse.

    :math:`1 - \exp(-\int W\, dt)` (Tong and Lin 2005, eq. 6), by the trapezoid rule.

    Parameters
    ----------
    rate : numpy.ndarray
        Rate per atom along the pulse, 1/s, shape ``(n,)``.
    times : numpy.ndarray
        The times, s, shape ``(n,)``.

    Returns
    -------
    float
        The fraction, in [0, 1].

    Examples
    --------
    >>> times = np.linspace(0, 1e-14, 11)
    >>> f"{ionised_fraction(np.full(11, 1e13), times):.4f}"
    '0.0952'
    """
    return float(-np.expm1(-np.trapezoid(rate, times)))


def pulse_yield(
    model: IonisationRateModel | MolecularRateModel,
    peak_intensity: float,
    fwhm: float,
    *,
    carrier: bool = False,
    cep: float = 0.0,
    cycles: tuple[float, float] | None = None,
) -> float:
    """Return the fraction a pulse at the model's wavelength ionises.

    A Gaussian pulse (:func:`gaussian_pulse`), or with ``cycles`` a pulse counted in
    cycles (:func:`cycle_pulse`).

    Parameters
    ----------
    model : IonisationRateModel or MolecularRateModel
        The rate.
    peak_intensity : float
        Peak of the cycle-averaged intensity, W/m^2.
    fwhm : float
        Full width at half maximum of the intensity, s.
    carrier : bool, optional
        True: the instantaneous rate on the field of :func:`gaussian_pulse`. False:
        the cycle-averaged rate on its envelope.
    cep : float, optional
        Carrier-envelope phase, rad; it matters with ``carrier`` only.
    cycles : tuple of float or None, optional
        ``(ramp, plateau)`` cycles of :func:`cycle_pulse`, in place of the Gaussian;
        ``fwhm`` is then ignored.

    Returns
    -------
    float
        The ionised fraction after the pulse.

    Examples
    --------
    Argon, 800 nm, 42 fs, 1e14 W/cm^2:

    >>> model = ionisation_rate_model("Ar", 800e-9)
    >>> f"{pulse_yield(model, 1e18, 42e-15):.4f}"
    '0.0115'
    """
    if cycles is None:
        times, envelope, field = gaussian_pulse(
            peak_intensity, fwhm, model.wavelength, cep=cep
        )
    else:
        times, envelope, field = cycle_pulse(
            peak_intensity, model.wavelength, *cycles, cep=cep
        )
    if carrier:
        return ionised_fraction(_evaluate(model, field, cycle_averaged=False), times)
    return ionised_fraction(_evaluate(model, envelope, cycle_averaged=True), times)


def focal_average(
    intensity: Float[np.ndarray, " n"], local: Float[np.ndarray, " n"]
) -> Float[np.ndarray, " n"]:
    r"""Return the yield of a thin target across a Gaussian beam, up to a constant.

    Atoms at the radius :math:`r` of a Gaussian beam of peak intensity :math:`I_0`
    see :math:`I = I_0 e^{-2r^2/w^2}`, so that :math:`2\pi r\, dr = (\pi w^2/2)\,
    dI/I` and the number ionised in a target much thinner than the Rayleigh range is

    .. math::

       S(I_0) = \int_0^{I_0} P(I)\, \frac{dI}{I}

    times :math:`\pi w^2/2` and the target's areal density: Wallace et al.'s (2016)
    two-dimensional focal-volume average (their eq. 3).

    Parameters
    ----------
    intensity : numpy.ndarray
        Peak intensities, ascending, shape ``(n,)``; the first must be low enough that
        the yield below it is negligible.
    local : numpy.ndarray
        The ionised fraction :math:`P` at each, shape ``(n,)``.

    Returns
    -------
    numpy.ndarray
        :math:`S` at each intensity, shape ``(n,)``.

    Examples
    --------
    A yield proportional to :math:`I^4` averages to a quarter of itself:

    >>> intensity = np.geomspace(1e-3, 1, 400)
    >>> f"{float(focal_average(intensity, intensity**4)[-1]):.3f}"
    '0.250'
    """
    logarithm = np.log(intensity)
    steps = 0.5 * (local[1:] + local[:-1]) * np.diff(logarithm)
    # Below the first intensity the yield is taken to follow the power law of its
    # first two points, whose integral is P / slope.
    slope = math.log(local[1] / local[0]) / (logarithm[1] - logarithm[0])
    return np.concatenate([[0.0], np.cumsum(steps)]) + local[0] / slope


def _by_key() -> dict[str, Dataset]:
    """Return the curated datasets by key.

    Returns
    -------
    dict of str to Dataset
        Every dataset.
    """
    return {d.key: d for d in datasets()}


def _numbers(
    dataset: Dataset,
) -> tuple[Float[np.ndarray, " n"], Float[np.ndarray, " n"]]:
    """Return the abscissae (SI) and values of a dataset, without its excluded rows.

    Parameters
    ----------
    dataset : Dataset
        The dataset.

    Returns
    -------
    tuple of numpy.ndarray
        The abscissae and the values, scaled.
    """
    return dataset.abscissa_values(), dataset.scaled_values()


def static_rates(gas: Gas | str) -> list[RateComparison]:
    """Compare the static rate of a gas with its curated references.

    Parameters
    ----------
    gas : Gas or str
        He, Ne or Ar: the gases with tabulated static rates.

    Returns
    -------
    list of RateComparison
        One comparison per reference of :data:`STATIC_REFERENCES`, in atomic units of
        rate against the field in V/m, at the fields up to twice the
        barrier-suppression field.

    Raises
    ------
    DataError
        If the gas has no static reference.

    Examples
    --------
    >>> (first, *_) = static_rates("He")
    >>> first.label
    'Parker et al. 2009, two electrons'
    >>> bool(np.all(np.abs(first.ratio[5:20] - 1) < 0.15))
    True
    """
    species = as_gas(gas)
    if species.value not in STATIC_REFERENCES:
        msg = f"no static ionisation rates are curated for {species.value}"
        raise DataError(msg)
    by_key = _by_key()
    comparisons = []
    for key, name, label in STATIC_REFERENCES[species.value]:
        model = ionisation_rate_model(species, math.inf, name)
        field, expected = _numbers(by_key[key])
        keep = field <= _FIELD_LIMIT * model.barrier_suppression_field
        predicted = (
            _evaluate(model, field[keep], cycle_averaged=False) / AU_ANGULAR_FREQUENCY
        )
        comparisons.append(
            RateComparison(
                species,
                key,
                label,
                name,
                "static rate",
                math.inf,
                field[keep],
                expected[keep],
                predicted,
            )
        )
    return comparisons


def _coefficients(source: str, formula: str) -> dict[str, float]:
    """Return the coefficients of a published formula.

    Parameters
    ----------
    source : str
        The source's id.
    formula : str
        The formula's id in its record.

    Returns
    -------
    dict of str to float
        The coefficients, from their printed strings.

    Raises
    ------
    DataError
        If the record has no such formula.
    """
    for entry in sources()[source].formulas:
        if entry["id"] == formula:
            table = entry["coefficients"]
            assert isinstance(table, Mapping)
            return {str(k): float(str(v)) for k, v in table.items()}  # pyright: ignore[reportUnknownVariableType]
    msg = f"{source} has no formula {formula!r}"
    raise DataError(msg)


def parker_rate_au(
    intensity: ArrayLike, wavelength_nm: int
) -> Float[np.ndarray, " *shape"]:
    r"""Return Parker et al.'s fit of the ionisation rate of helium in a laser field.

    Parker, Meharg, McKenna and Taylor (2007), eq. (5), fitted "through the centre"
    of the
    rates of their integration of the two-electron time-dependent Schrödinger
    equation:

    .. math::

       \Gamma = \exp\left(-\frac{R}{\sqrt{a_0 + a_1 \bar I + a_2 \bar I^2 + a_3
           \bar I^3}}\right), \qquad \bar I = I R^2, \quad R = \lambda / 390\
           \text{nm},

    with :math:`I` in :math:`10^{14}` W/cm\ :sup:`2` and :math:`\Gamma` in atomic
    units. The rate is that of the atom in a field of constant amplitude, the
    resonances smoothed away. The fits cover 0.3 to 34 (390 nm) and 0.3 to 14.5
    (780 nm); no accuracy is stated for them, and the calculation's is 2 to 4 %.

    Parameters
    ----------
    intensity : array_like
        Intensity, W/m\ :sup:`2`.
    wavelength_nm : int
        390 or 780.

    Returns
    -------
    numpy.ndarray
        The rate, atomic units.

    Raises
    ------
    DataError
        If the wavelength is not one of the two.

    Examples
    --------
    The paper's check, 0.0107 au at 30e14 W/cm\ :sup:`2` and 390 nm:

    >>> f"{float(parker_rate_au(30e18, 390)):.4f}"
    '0.0107'
    """
    if wavelength_nm not in {390, 780}:
        msg = f"Parker et al. fitted 390 and 780 nm, not {wavelength_nm} nm"
        raise DataError(msg)
    c = _coefficients("parker_2007", "eq5_intense_field")
    ratio = wavelength_nm / 390
    scaled = np.asarray(intensity, dtype=float) / 1e18 * ratio**2
    polynomial = sum(
        c[f"a{power}_{wavelength_nm}"] * scaled**power for power in range(4)
    )
    return np.asarray(np.exp(-ratio / np.sqrt(polynomial)))


def parker_pulse_yield(peak_intensity: float, fwhm: float, wavelength_nm: int) -> float:
    r"""Return the fraction of helium that Parker et al.'s rate ionises over a pulse.

    :func:`parker_rate_au` at the intensity of the envelope, :math:`I_0 e(t)^2`, of
    :func:`gaussian_pulse`, integrated by :func:`ionised_fraction`: the rate of a
    field of constant amplitude followed adiabatically (TC-1050).

    Parameters
    ----------
    peak_intensity : float
        Peak intensity, W/m\ :sup:`2`.
    fwhm : float
        Full width at half maximum of the intensity, s.
    wavelength_nm : int
        390 or 780.

    Returns
    -------
    float
        The fraction ionised.

    Examples
    --------
    >>> f"{parker_pulse_yield(5.68981e18, 6e-15, 390):.3e}"
    '7.000e-03'
    """
    times, envelope, _ = gaussian_pulse(peak_intensity, fwhm, wavelength_nm * 1e-9)
    intensity = peak_intensity * (envelope / np.max(envelope)) ** 2
    rate = parker_rate_au(intensity, wavelength_nm) * AU_ANGULAR_FREQUENCY
    return ionised_fraction(rate, times)


def helium_laser_rates(model: str = "default") -> list[RateComparison]:
    """Compare helium's rate at 390 and 780 nm with Parker et al.'s.

    Parameters
    ----------
    model : str, optional
        The rate model.

    Returns
    -------
    list of RateComparison
        The cycle-averaged rate in atomic units against the intensity, at 780 nm
        (1e14 to 1.4e15 W/cm^2) and 390 nm (1e14 to 3e15 W/cm^2): the ranges of the
        fits above the intensity at which the rate is 1e-12 au.

    Examples
    --------
    >>> infrared, _ = helium_laser_rates()
    >>> bool(np.all(np.abs(infrared.ratio - 1) < 0.2))
    True
    """
    comparisons = []
    for wavelength_nm, upper in ((780, 1.4e19), (390, 3e19)):
        rate = ionisation_rate_model("He", wavelength_nm * 1e-9, model)
        intensity = np.geomspace(1e18, upper, 25)
        predicted = _evaluate(rate, amplitude(intensity), cycle_averaged=True)
        comparisons.append(
            RateComparison(
                Gas.He,
                "parker_2007/eq5_intense_field",
                f"Parker et al. 2007, two electrons, {wavelength_nm} nm",
                model,
                "rate",
                wavelength_nm * 1e-9,
                intensity,
                parker_rate_au(intensity, wavelength_nm),
                predicted / AU_ANGULAR_FREQUENCY,
            )
        )
    return comparisons


RMATRIX_RATES: Mapping[str, tuple[tuple[str, str], ...]] = {
    "He": (
        ("vanderhart_2006/fig6_floquet_He", "R-matrix Floquet (van der Hart 2006)"),
    ),
    "Ne": (
        ("vanderhart_2006/fig6_floquet_Ne", "R-matrix Floquet (van der Hart 2006)"),
    ),
    "Ar": (
        ("vanderhart_2006/fig6_floquet_Ar", "R-matrix Floquet (van der Hart 2006)"),
        (
            "vanderhart_2007/table1_one_state_Ar",
            "time-dependent R-matrix (van der Hart 2007, Table I)",
        ),
        (
            "vanderhart_2007/fig4_td_rmatrix_Ar",
            "time-dependent R-matrix (van der Hart 2007, Fig. 4)",
        ),
    ),
}
"""The all-electron rates at 390 nm, by gas: dataset and label (TC-700, TC-701)."""


def rmatrix_rates(gas: Gas | str, model: str = "default") -> list[RateComparison]:
    """Compare a gas's rate at 390 nm with van der Hart's all-electron rates.

    The R-matrix Floquet rates of He, Ne and Ar (van der Hart 2006, read from the image
    of its Fig. 6) and the time-dependent R-matrix rates of Ar (van der Hart 2007,
    Table I and its Fig. 4's vector drawing), each for the atom with one ionic target
    state, at 390 nm with linear polarisation: Keldysh parameters from 1.4 to 4.7. The
    references resolve the resonances of the atom, which no model has.

    Parameters
    ----------
    gas : Gas or str
        He, Ne or Ar.
    model : str, optional
        The rate model.

    Returns
    -------
    list of RateComparison
        The cycle-averaged rate in atomic units against the intensity, W/m^2.

    Raises
    ------
    DataError
        If the gas has no such reference.
    """
    species = as_gas(gas)
    if species.value not in RMATRIX_RATES:
        msg = f"no R-matrix rates at 390 nm are curated for {species.value}"
        raise DataError(msg)
    rate = ionisation_rate_model(species, 390e-9, model)
    by_key = _by_key()
    comparisons = []
    for key, label in RMATRIX_RATES[species.value]:
        intensity, values = _numbers(by_key[key])
        predicted = _evaluate(rate, amplitude(intensity), cycle_averaged=True)
        comparisons.append(
            RateComparison(
                species,
                key,
                label,
                model,
                "rate",
                390e-9,
                intensity,
                values / AU_ANGULAR_FREQUENCY,
                predicted / AU_ANGULAR_FREQUENCY,
            )
        )
    return comparisons


TDSE_YIELDS: Mapping[str, tuple[tuple[str, str], ...]] = {
    "Ar": (
        ("travers_2026_trecx/ionised_gauss6_Ar_m0", "the ionised fraction"),
        ("travers_2026_trecx/ionised_gauss10_400nm_Ar_m0", "the ionised fraction"),
        ("travers_2026_trecx/ionised_gauss10_1030nm_Ar_m0", "the ionised fraction"),
        ("travers_2026_trecx/ionised_gauss10_800nm_Ar_m0", "the ionised fraction"),
        ("travers_2026_trecx/ionised_wallace_cep0_Ar_m0", "the ionised fraction"),
        ("travers_2026_trecx/ionised_wallace_cep90_Ar_m0", "the ionised fraction"),
        ("travers_2026_trecx/ionised_gauss30_400nm_Ar_m0", "the ionised fraction"),
        ("travers_2026_trecx/ionised_gauss30_1030nm_Ar_m0", "the ionised fraction"),
    ),
    "He": (("travers_2026_trecx/ionised_gauss6_He", "the ionised fraction"),),
    "Kr": (
        ("travers_2026_trecx/ionised_gauss6_Kr_m0", "the ionised fraction"),
        ("travers_2026_trecx/ionised_wallace_cep0_Kr_m0", "the ionised fraction"),
        ("travers_2026_trecx/ionised_wallace_cep90_Kr_m0", "the ionised fraction"),
    ),
    "Xe": (
        ("travers_2026_trecx/ionised_gauss6_Xe_m0", "the ionised fraction"),
        ("travers_2026_trecx/ionised_wallace_cep0_Xe_m0", "the ionised fraction"),
        ("travers_2026_trecx/ionised_wallace_cep90_Xe_m0", "the ionised fraction"),
    ),
}
"""The yields of the one-electron calculations, by gas: dataset and measure.

The ionised fraction, the bound states up to n = 10 subtracted (TC-812); for argon at
800 nm (6 fs), 400 nm (10 fs, TC-824), 1030 nm (10 fs, TC-826) and 800 nm with the
same 10 fs pulse (TC-830); Wallace et al.'s 5.5 fs pulse at 800 nm with
carrier-envelope phases 0 and pi/2 (TC-836), each compared at its own phase; 30 fs at
400 and 1030 nm (TC-844); krypton and xenon in A2's 6 fs pulse (TC-845). The first
runs' measures (``chi_gauss_*``, ``norm_gauss_*``, ``isurff_cos8_*``) count Rydberg
states as ionised or miss part of the continuum (TC-703), and are not compared."""


def tdse_yields(
    gas: Gas | str, model: str = "one_electron", *, carrier: bool = True
) -> list[RateComparison]:
    """Compare the yields of short pulses with the one-electron calculations.

    Parameters
    ----------
    gas : Gas or str
        He or Ar.
    model : str, optional
        The rate model; the calculations are of one electron in the m = 0 orbital, to
        which ``one_electron`` corresponds.
    carrier : bool, optional
        The instantaneous rate on the field (the calculations' pulse,
        :func:`gaussian_pulse` with a carrier-envelope phase of 0), or the
        cycle-averaged rate on the envelope.

    Returns
    -------
    list of RateComparison
        One comparison per measure of the yield (:data:`TDSE_YIELDS`), against the
        peak intensity.

    Raises
    ------
    DataError
        If the gas has no such calculation.
    """
    species = as_gas(gas)
    if species.value not in TDSE_YIELDS:
        msg = f"no TDSE yields are curated for {species.value}"
        raise DataError(msg)
    by_key = _by_key()
    comparisons = []
    for key, measure in TDSE_YIELDS[species.value]:
        dataset = by_key[key]
        wavelength = float(str(dataset.parameters["wavelength_nm"])) * 1e-9
        fwhm = float(str(dataset.parameters["pulse_fs"])) * 1e-15
        cep = float(str(dataset.parameters.get("cep", 0.0)))
        rate = ionisation_rate_model(species, wavelength, model)
        intensity, expected = _numbers(dataset)
        keep = amplitude(intensity) <= _FIELD_LIMIT * rate.barrier_suppression_field
        predicted = np.array(
            [
                pulse_yield(rate, peak, fwhm, carrier=carrier, cep=cep)
                for peak in intensity[keep]
            ]
        )
        comparisons.append(
            RateComparison(
                species,
                key,
                f"one-electron TDSE, {measure}",
                model,
                "yield",
                wavelength,
                intensity[keep],
                expected[keep],
                predicted,
            )
        )
    return comparisons


HYDROGEN_ALPHA = "tong_2005/alpha_H"
"""Tong and Lin's constant of the static correction for atomic hydrogen (TC-833)."""


def hydrogen_model(wavelength: float, coulomb: str = "static") -> IonisationRateModel:
    r"""Return the rate of atomic hydrogen built as optmat's ``one_electron`` model.

    Atomic hydrogen is not one of optmat's gases, but the cleanest test of the rate's
    Coulomb factor: its one-electron calculations are exact. The model is the
    ``one_electron`` form for the 1s electron: :math:`I_p = 1/2` hartree, :math:`l = 0`,
    one electron, the asymptotic coefficient :math:`C_0 = 2` (the 1s orbital
    :math:`\pi^{-1/2} e^{-r}` is :math:`C_0 \kappa^{3/2} (\kappa r)^{n^* - 1}
    e^{-\kappa r} Y_{00}` with :math:`\kappa = n^* = 1`), Tong and Lin's static
    correction with their :math:`\alpha` for hydrogen (:data:`HYDROGEN_ALPHA`), and the
    default's channel sum and averaging (TC-827).

    Parameters
    ----------
    wavelength : float
        Wavelength, m.
    coulomb : str, optional
        The Coulomb factor (:func:`optmat.gas.ionisation_rates.coulomb_factor`).

    Returns
    -------
    IonisationRateModel
        The model, status ``user``.

    Examples
    --------
    >>> model = hydrogen_model(800e-9)
    >>> f"{float(model.keldysh(amplitude(1e18))):.3f}"
    '1.067'
    """
    alpha = float(_by_key()[HYDROGEN_ALPHA].scaled_values()[0])
    return IonisationRateModel.from_parameters(
        "H",
        wavelength,
        ionization_energy=0.5 * AU_ENERGY,
        orbital_l=0,
        electrons=1.0,
        structure="tong_2002",
        structure_coefficient=2.0,
        coulomb=coulomb,  # pyright: ignore[reportArgumentType] - checked on construction
        static_correction="tong-lin",
        static_constants=(alpha,),
        citation="atomic hydrogen, C_0 = 2, Tong and Lin's alpha (tong_2005)",
        model=f"hydrogen_{coulomb}",
        sources={"static_correction": HYDROGEN_ALPHA},
    )


HYDROGEN_YIELDS: Mapping[str, tuple[float, float] | None] = {
    "zhao_2016/fig2a_tdse_H_600nm": None,
    "zhao_2016/fig2b_tdse_H_800nm": None,
    "zhao_2016/fig2c_tdse_H_1000nm": None,
    "zhao_2016/fig2d_tdse_H_1200nm": None,
    "travers_2026_trecx/ionised_gauss15_400nm_H": None,
    "travers_2026_trecx/ionised_gauss15_800nm_H": None,
    "le_2016/fig12a_tdse_H_400nm": (6.0, 0.0),
    "le_2016/fig12b_tdse_H_800nm": (3.0, 0.0),
    "le_2016/fig12c_tdse_H_1600nm": (3.0, 0.0),
    "bauer_2018/fig1_tdse_H_800nm": (2.0, 10.0),
    "bauer_2018/fig5a_tdse_H_800nm": (2.0, 10.0),
    "bauer_2018/fig6a_tdse_H_400nm": (2.0, 10.0),
}
r"""The one-electron calculations of hydrogen against intensity, and their pulses.

None: the dataset's Gaussian (``pulse_fs`` the full width at half maximum of the
intensity, as Zhao et al.'s own PPT curves are reproduced with). Otherwise the
``(ramp, plateau)`` cycles of :func:`cycle_pulse`: Le et al.'s :math:`\sin^2` pulses of
12 and 6 cycles, taken on the field (which their PPT curves favour; on the intensity
the model's yield is about 1.35 times higher), and Bauer and Deng's 2-10-2 flat top
(TC-833 to TC-835, TC-843)."""

HYDROGEN_SPECTRA: tuple[str, ...] = (
    "bauer_2018/fig7_tdse_H_spectrum_F002",
    "bauer_2018/fig8_tdse_H_spectrum_F0063",
)
"""Bauer and Deng's calculations of hydrogen against the frequency at fixed field."""

_HYDROGEN_PRECISION = 0.1
"""Largest stated uncertainty, relative, of a hydrogen yield compared."""

_SPECTRUM_STEP = 0.005
"""Spacing, hartree, of the frequencies at which a spectrum is compared."""

_SPECTRUM_MAXIMUM = 0.25
"""Highest frequency, hartree, compared: above it the yield is resonant or
single-photon, outside the domain of a tunnelling rate."""


def hydrogen_yields(coulomb: str = "static") -> list[RateComparison]:
    """Compare the hydrogen model with the one-electron calculations of hydrogen.

    The cycle-averaged rate of :func:`hydrogen_model` on each calculation's envelope
    (:data:`HYDROGEN_YIELDS`; no calculation states its carrier-envelope phase), at
    the points whose stated uncertainty is at most 0.1 of the value: the low yields of
    Bauer and Deng's Fig. 1, read on a linear axis to 0.001, are left out.

    Parameters
    ----------
    coulomb : str, optional
        The Coulomb factor.

    Returns
    -------
    list of RateComparison
        One per calculation, against the peak intensity.
    """
    by_key = _by_key()
    comparisons = []
    for key, cycles in HYDROGEN_YIELDS.items():
        dataset = by_key[key]
        wavelength = _parameter(dataset, "wavelength_nm") * 1e-9
        fwhm = _parameter(dataset, "pulse_fs") * 1e-15
        rate = hydrogen_model(wavelength, coulomb)
        intensity, expected = _numbers(dataset)
        precise = dataset.standard_uncertainty() <= _HYDROGEN_PRECISION * expected
        intensity, expected = intensity[precise], expected[precise]
        predicted = np.array(
            [pulse_yield(rate, peak, fwhm, cycles=cycles) for peak in intensity]
        )
        comparisons.append(
            RateComparison(
                "H",
                key,
                "one-electron TDSE of hydrogen",
                rate.model,
                "yield",
                wavelength,
                intensity,
                expected,
                predicted,
            )
        )
    return comparisons


def hydrogen_spectra(coulomb: str = "static") -> list[RateComparison]:
    """Compare the hydrogen model with Bauer and Deng's yields against the frequency.

    At their fixed peak intensity, for their 2-10-2 flat-top pulse, whose duration
    scales with the period; at the vertices of their curve nearest to frequencies
    0.005 hartree apart, up to 0.25 hartree. Above 0.13 hartree
    the wavelength is below the model's tested range (350 nm): the comparison is made
    there on purpose, without the warning.

    Parameters
    ----------
    coulomb : str, optional
        The Coulomb factor.

    Returns
    -------
    list of RateComparison
        One per spectrum, against the angular frequency (hartree); ``wavelength`` is
        NaN.
    """
    by_key = _by_key()
    comparisons = []
    for key in HYDROGEN_SPECTRA:
        dataset = by_key[key]
        peak = _parameter(dataset, "intensity_w_cm2") * 1e4
        # The abscissa is the angular frequency in hartree, as printed.
        omega = np.asarray(dataset.abscissa_values(), dtype=float)
        values = dataset.scaled_values()
        grid = np.arange(omega.min(), _SPECTRUM_MAXIMUM, _SPECTRUM_STEP)
        chosen = np.unique([int(np.argmin(np.abs(omega - w))) for w in grid])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ModelRangeWarning)
            predicted = np.array(
                [
                    pulse_yield(
                        hydrogen_model(
                            2 * math.pi * C / (float(omega[k]) * AU_ANGULAR_FREQUENCY),
                            coulomb,
                        ),
                        peak,
                        0.0,
                        cycles=(2.0, 10.0),
                    )
                    for k in chosen
                ]
            )
        comparisons.append(
            RateComparison(
                "H",
                key,
                "one-electron TDSE of hydrogen against the frequency",
                f"hydrogen_{coulomb}",
                "yield",
                math.nan,
                omega[chosen],
                values[chosen],
                predicted,
            )
        )
    return comparisons


HACC_YIELDS: tuple[str, ...] = (
    "travers_2026_trecx/hacc_normloss_gauss6_Ar_c4",
    "travers_2026_trecx/hacc_normloss_gauss6_Ar_c3",
    "travers_2026_trecx/hacc_normloss_gauss6_Ar_c3m1",
    "travers_2026_trecx/hacc_normloss_gauss20_Ar_c3",
)
"""haCC's all-electron argon, the norm-loss 1000 au after the pulse: A2's 6 fs pulse
with 4 and 3 ion channels and with mMax 1, and a 20 fs pulse (TC-847, TC-861)."""

HACC_REFERENCES: Mapping[str, str] = {
    "travers_2026_trecx/hacc_normloss_gauss6_Ar_c3": (
        "travers_2026_trecx/normloss_ipmatched_c3_Ar_m0"
    ),
    "travers_2026_trecx/hacc_normloss_gauss20_Ar_c3": (
        "travers_2026_trecx/normloss_ipmatched_c3_20fs_Ar_m0"
    ),
}
"""haCC's yields with one electron run at exactly their ionisation energy and pulse
(TC-861); the others interpolate between :data:`IP_MATCHED_YIELDS`."""

IP_MATCHED_YIELDS: tuple[str, ...] = (
    "travers_2026_trecx/normloss_ipmatched_ci_Ar_m0",
    "travers_2026_trecx/normloss_ipmatched_c3m1_Ar_m0",
)
"""One electron (m = 0) at two of haCC's ionisation energies, the same pulse and
measure: the reference haCC's yields are divided by (TC-847)."""


@dataclass(frozen=True)
class AtomRatio:
    """An atom's all-electron yield over one electron's at the same ionisation energy.

    Attributes
    ----------
    reference : str
        The all-electron dataset.
    ionization_energy : float
        Its effective ionisation energy, hartree.
    intensity : numpy.ndarray
        The peak intensities, W/m^2.
    atom : numpy.ndarray
        The atom's norm-loss.
    one_electron : numpy.ndarray
        One m = 0 electron's at the same energy: the logarithm of the two calculations
        of :data:`IP_MATCHED_YIELDS` interpolated linearly in the ionisation energy.
    """

    reference: str
    ionization_energy: float
    intensity: Float[np.ndarray, " n"]
    atom: Float[np.ndarray, " n"]
    one_electron: Float[np.ndarray, " n"]

    @property
    def ratio(self) -> Float[np.ndarray, " n"]:
        """The atom's yield over one electron's.

        Returns
        -------
        numpy.ndarray
            The ratio at each intensity.
        """
        return self.atom / self.one_electron


def hacc_ratios() -> list[AtomRatio]:
    """Compare haCC's argon atom with one electron at its ionisation energy.

    Each set-up of :data:`HACC_YIELDS` has its own effective ionisation energy (the
    initial state's, against the ion's); one electron's norm-loss there is the run at
    exactly that energy and pulse where there is one (:data:`HACC_REFERENCES`;
    TC-861), and otherwise interpolated, in its logarithm and linearly in the energy,
    between the two calculations of :data:`IP_MATCHED_YIELDS` at the same intensities.
    Both measures are the norm-loss 1000 au after the pulse (TC-847).

    Returns
    -------
    list of AtomRatio
        One per haCC set-up.
    """
    by_key = _by_key()
    references = []
    for key in IP_MATCHED_YIELDS:
        dataset = by_key[key]
        intensity, values = _numbers(dataset)
        energy = _parameter(dataset, "model_ionization_energy_au")
        references.append((energy, dict(zip(intensity, values, strict=True))))
    (low, at_low), (high, at_high) = sorted(references, key=lambda r: r[0])
    ratios = []
    for key in HACC_YIELDS:
        dataset = by_key[key]
        intensity, atom = _numbers(dataset)
        energy = _parameter(dataset, "model_ionization_energy_au")
        if key in HACC_REFERENCES:
            own_intensity, own = _numbers(by_key[HACC_REFERENCES[key]])
            exact = dict(zip(own_intensity, own, strict=True))
            one = np.array([exact[i] for i in intensity])
        else:
            weight = (energy - low) / (high - low)
            one = np.array(
                [at_low[i] * (at_high[i] / at_low[i]) ** weight for i in intensity]
            )
        ratios.append(AtomRatio(key, energy, intensity, atom, one))
    return ratios


SATURATION_INTENSITIES: Mapping[str, str] = {
    "He": "uiterwaal_1998/isat_He",
    "Ar": "uiterwaal_1998/isat_Ar",
    "Kr": "uiterwaal_1998/isat_Kr",
}
"""Uiterwaal et al.'s saturation intensities at 248.6 nm, 500 fs (TC-853)."""


@dataclass(frozen=True)
class SaturationTest:
    """A model's saturation of an atom against a measured saturation intensity.

    Attributes
    ----------
    gas : Gas
        The atom.
    model : str
        The rate model.
    measured : float
        The measured saturation intensity, W/m^2.
    sigma : float
        Its stated uncertainty, W/m^2.
    saturation : float
        The model's: the peak intensity at which the integral of its rate over the
        pulse is 1, W/m^2.
    exponent : float
        The integral of the model's rate over the pulse at the measured saturation
        intensity: 1 if the model matches, the factor by which its rate is high
        otherwise.
    keldysh : float
        The Keldysh parameter at the measured saturation intensity.
    """

    gas: Gas
    model: str
    measured: float
    sigma: float
    saturation: float
    exponent: float
    keldysh: float

    @property
    def scale(self) -> float:
        """The model's saturation intensity over the measured one.

        Returns
        -------
        float
            The scale.
        """
        return self.saturation / self.measured


def saturation_test(gas: Gas | str, model: str = "default") -> SaturationTest:
    r"""Test a model against Uiterwaal et al.'s saturation intensity of an atom.

    Their saturation intensity (TC-853) is the peak intensity at which a single rate
    ionises :math:`1 - 1/e` of the atoms over the pulse: the integral of the rate over
    the pulse is 1 there. The model's cycle-averaged rate is integrated over their
    pulse, a hyperbolic secant squared of 500 fs full width at half maximum of the
    intensity at 248.6 nm (the wavelength is below the models' range, whose warning
    is suppressed), and the peak intensity at which the integral is 1 is found by
    bisection in its logarithm.

    Parameters
    ----------
    gas : Gas or str
        He, Ar or Kr.
    model : str, optional
        The rate model.

    Returns
    -------
    SaturationTest
        The comparison.

    Raises
    ------
    DataError
        If no saturation intensity is curated for the gas.
    """
    species = as_gas(gas)
    if species.value not in SATURATION_INTENSITIES:
        msg = f"no measured saturation intensity of {species.value}"
        raise DataError(msg)
    dataset = _by_key()[SATURATION_INTENSITIES[species.value]]
    wavelength = float(dataset.vacuum_wavelength()[0])
    to_si = {"W/cm2": 1e4, "W/m2": 1.0}[dataset.unit]
    measured = float(dataset.scaled_values()[0]) * to_si
    sigma = float(dataset.standard_uncertainty()[0]) * to_si
    fwhm = _parameter(dataset, "pulse_fs") * 1e-15
    width = fwhm / (2 * math.acosh(math.sqrt(2)))
    time = np.linspace(-4 * fwhm, 4 * fwhm, 4001)
    envelope = 1 / np.cosh(time / width) ** 2
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        rate = ionisation_rate_model(species, wavelength, model)

        def exponent(peak: float) -> float:
            values = rate.cycle_averaged_rate(amplitude(peak * envelope))
            return float(np.trapezoid(values, time))

        low, high = math.log(measured / 100), math.log(measured * 100)
        for _ in range(60):
            middle = 0.5 * (low + high)
            if exponent(math.exp(middle)) > 1:
                high = middle
            else:
                low = middle
        keldysh = float(rate.keldysh(amplitude(measured)))
        return SaturationTest(
            species,
            model,
            measured,
            sigma,
            math.exp(0.5 * (low + high)),
            exponent(measured),
            keldysh,
        )


WEAK_FIELD_ATOMS: tuple[str, ...] = ("He", "Ne", "Ar", "Kr", "Xe")
"""The atoms of Tolstikhina et al.'s Table I that optmat curates (TC-852)."""


@dataclass(frozen=True)
class WeakFieldFactor:
    """An atom's weak-field rate over the model's one electron.

    Attributes
    ----------
    gas : Gas
        The atom.
    kappa : float
        The Hartree-Fock orbital exponent at which both coefficients are compared, au.
    atom : float
        The atom's coefficient g of the many-electron weak-field asymptotic theory
        (Hartree-Fock, with relaxation and the occupation factor).
    one_electron : float
        The model's one electron's on the same scale, from its structure coefficient.
    """

    gas: Gas
    kappa: float
    atom: float
    one_electron: float

    @property
    def factor(self) -> float:
        """The atom's weak-field rate over the model's one electron, (g/g_1)^2.

        Returns
        -------
        float
            The factor.
        """
        return (self.atom / self.one_electron) ** 2


def weak_field_factors() -> list[WeakFieldFactor]:
    r"""Compare the atoms' weak-field coefficients with the model's one electron.

    Tolstikhina et al. (2014) give the rate of an atom in a weak static field as
    :math:`|g|^2 W(F)`, with :math:`g = \sqrt{p} R_N g_N` from Hartree-Fock orbitals
    (Table I). For one electron of the model, whose orbital's tail is Tong, Zhao &
    Lin's :math:`C_l r^{1/\kappa - 1} e^{-\kappa r}`, their Eq. (11) gives
    :math:`g_1 = \sqrt{2(2l+1)/\kappa}\, 2^{-1/\kappa} C_l`, with which their
    Eqs. (1)-(2) are the static rate of :func:`~optmat.gas.ionisation_rates` without
    Tong & Lin's correction, which tends to 1 at weak field. Both at the Hartree-Fock
    :math:`\kappa` (TC-852).

    Returns
    -------
    list of WeakFieldFactor
        One per atom of :data:`WEAK_FIELD_ATOMS`.
    """
    by_key = _by_key()
    factors = []
    for name in WEAK_FIELD_ATOMS:
        gas = as_gas(name)
        kappa = float(by_key[f"tolstikhina_2014/kappa_{name}"].scaled_values()[0])
        atom = float(by_key[f"tolstikhina_2014/g_{name}"].scaled_values()[0])
        model = ionisation_rate_model(gas, math.inf, "one_electron")
        coefficient = model.structure_coefficient
        if coefficient is None:
            msg = f"the one-electron model of {name} has no structure coefficient"
            raise DataError(msg)
        one = (
            math.sqrt(2 * (2 * model.orbital_l + 1) / kappa)
            * 2 ** (-1 / kappa)
            * coefficient
        )
        factors.append(WeakFieldFactor(gas, kappa, atom, one))
    return factors


PUBLISHED_ATOMS: Mapping[str, str] = {
    "pathak_2020/fig2_tdcasscf_Ar": "travers_2026_trecx/outside20_sin2_m1_Ar_m0",
    "pathak_2020/fig2_sae_Ar": "travers_2026_trecx/outside20_sin2_m1mk_ArMK_m0",
    "pathak_2022/fig2_tdcasscf_Kr": "travers_2026_trecx/outside20_sin2_m2_Kr_m0",
    "sato_2018/fig3a_casscf_Ar_800nm": "travers_2026_trecx/outside20_sin2_m3_Ar_m0",
    "sato_2016/fig6_casscf_He_800nm": "travers_2026_trecx/outside20_sin2_m4_He",
}
"""Published many-electron calculations and one electron in their pulses (TC-858).

Each published curve (the authors' best many-electron calculation, TD-CASSCF; and
Pathak et al. 2020's own single-active-electron curve, against one electron in their
potential) with tRecX's one m = 0 electron in the same sin^2 pulse, read with the same
measure: the probability outside 20 au at the end of the pulse."""


def published_atom_ratios() -> list[AtomRatio]:
    """Compare published many-electron calculations with one electron in their pulses.

    For each pair of :data:`PUBLISHED_ATOMS`, at the published calculation's
    intensities: the published probability outside 20 au over one electron's (TC-858).

    Returns
    -------
    list of AtomRatio
        One per published curve; ``ionization_energy`` is the one electron's.
    """
    by_key = _by_key()
    ratios = []
    for key, reference in PUBLISHED_ATOMS.items():
        intensity, atom = _numbers(by_key[key])
        own_intensity, own = _numbers(by_key[reference])
        one = {float(i): float(v) for i, v in zip(own_intensity, own, strict=True)}
        energy = _parameter(by_key[reference], "model_ionization_energy_au")
        ratios.append(
            AtomRatio(
                key,
                energy,
                intensity,
                atom,
                np.array([one[float(i)] for i in intensity]),
            )
        )
    return ratios


HELIUM_TWO_ELECTRON: Mapping[str, str] = {
    "travers_2026_trecx/twoelectron_normloss_gauss6_390nm_He": (
        "travers_2026_trecx/normloss_ipmatched_390nm_base_He"
    ),
    "travers_2026_trecx/twoelectron_normloss_gauss6_390nm_He_m5": (
        "travers_2026_trecx/normloss_ipmatched_390nm_m5_He"
    ),
    "travers_2026_trecx/hacc_normloss_gauss6_390nm_He": (
        "travers_2026_trecx/normloss_ipmatched_390nm_base_He"
    ),
}
"""Helium at 390 nm with two electrons (and haCC's) against one electron at the same
ionisation energy, pulse and measure, the norm-loss 1000 au after the pulse
(TC-1050)."""


def helium_two_electron_ratios() -> list[AtomRatio]:
    """Compare helium's two-electron and haCC yields at 390 nm with one electron's.

    For each pair of :data:`HELIUM_TWO_ELECTRON`, at the two-electron (or haCC)
    calculation's intensities (TC-1050).

    Returns
    -------
    list of AtomRatio
        One per pair; ``ionization_energy`` is the one electron's.
    """
    by_key = _by_key()
    ratios = []
    for key, reference in HELIUM_TWO_ELECTRON.items():
        intensity, atom = _numbers(by_key[key])
        own_intensity, own = _numbers(by_key[reference])
        one = {float(i): float(v) for i, v in zip(own_intensity, own, strict=True)}
        energy = _parameter(by_key[reference], "model_ionization_energy_au")
        ratios.append(
            AtomRatio(
                key,
                energy,
                intensity,
                atom,
                np.array([one[float(i)] for i in intensity]),
            )
        )
    return ratios


@dataclass(frozen=True)
class YieldRatio:
    """One calculation's yield over another's, in one pulse and with one measure.

    Attributes
    ----------
    reference : str
        The dataset compared.
    against : str
        The dataset it is divided by.
    intensity : numpy.ndarray
        The peak intensities, W/m^2.
    value : numpy.ndarray
        ``reference``'s yields.
    against_value : numpy.ndarray
        ``against``'s at the same intensities.
    """

    reference: str
    against: str
    intensity: Float[np.ndarray, " n"]
    value: Float[np.ndarray, " n"]
    against_value: Float[np.ndarray, " n"]

    @property
    def ratio(self) -> Float[np.ndarray, " n"]:
        """The yield over the other's.

        Returns
        -------
        numpy.ndarray
            ``value / against_value``.
        """
        return self.value / self.against_value

    @property
    def rate_ratio(self) -> Float[np.ndarray, " n"]:
        r"""The ratio of the two in rate form.

        Returns
        -------
        numpy.ndarray
            :math:`\ln(1 - Y) / \ln(1 - Y')`: the ratio of the rates integrated over
            the pulse, if the yield is :math:`1 - \exp(-\int W dt)`; for one
            calculation over one electron's, an effective number of electrons.
        """
        return np.log1p(-self.value) / np.log1p(-self.against_value)


HELIUM_MATCHED: Mapping[str, str] = {
    "travers_2026_pyroraptor/twoelectron_above4_gauss6_390nm_He": (
        "travers_2026_pyroraptor/sae_above4_gauss6_390nm_He"
    ),
    "travers_2026_pyroraptor/hacc1_above4_gauss6_390nm_He": (
        "travers_2026_pyroraptor/twoelectron_above4_gauss6_390nm_He"
    ),
    "travers_2026_pyroraptor/hacc5_above4_gauss6_390nm_He": (
        "travers_2026_pyroraptor/twoelectron_above4_gauss6_390nm_He"
    ),
    "travers_2026_pyroraptor/hacc1_above5_gauss6_390nm_He": (
        "travers_2026_pyroraptor/twoelectron_above5_gauss6_390nm_He"
    ),
    "travers_2026_pyroraptor/hacc5_above5_gauss6_390nm_He": (
        "travers_2026_pyroraptor/twoelectron_above5_gauss6_390nm_He"
    ),
    "travers_2026_pyroraptor/hacc1_outside20_sin2_800nm_He": (
        "sato_2016/fig6_casscf_He_800nm"
    ),
    "travers_2026_pyroraptor/hacc5_outside20_sin2_800nm_He": (
        "sato_2016/fig6_casscf_He_800nm"
    ),
}
"""Helium's calculations against each other on measures converged in the box (TC-1092).

At 390 nm, in A2's 6 fs pulse, the population above the shells n <= 4 or n <= 5 1000
au after the pulse: two electrons over one, and haCC over two electrons; at 800 nm, in
Sato et al.'s sin^2 pulse, haCC's probability outside 20 au over their MCTDHF."""


def helium_matched_ratios() -> list[YieldRatio]:
    """Compare helium's calculations on matched, box-converged measures.

    For each pair of :data:`HELIUM_MATCHED`, at the first dataset's intensities
    (TC-1092).

    Returns
    -------
    list of YieldRatio
        One per pair.
    """
    by_key = _by_key()
    ratios = []
    for key, against in HELIUM_MATCHED.items():
        intensity, value = _numbers(by_key[key])
        own_intensity, own = _numbers(by_key[against])
        other = {float(i): float(v) for i, v in zip(own_intensity, own, strict=True)}
        ratios.append(
            YieldRatio(
                key,
                against,
                intensity,
                value,
                np.array([other[float(i)] for i in intensity]),
            )
        )
    return ratios


@dataclass(frozen=True)
class SurvivalRatio:
    """An atom's ground-state depletion against one electron's in the same pulse.

    Attributes
    ----------
    reference : str
        The atom's ground-state survival.
    intensity : numpy.ndarray
        The peak intensities, W/m^2, as the reference labels them.
    atom : numpy.ndarray
        The atom's survival.
    m0 : numpy.ndarray
        One electron's survival in the orbital with m = 0, at the same ionisation
        energy, pulse and measure.
    m1 : numpy.ndarray
        The same with m = 1.
    """

    reference: str
    intensity: Float[np.ndarray, " n"]
    atom: Float[np.ndarray, " n"]
    m0: Float[np.ndarray, " n"]
    m1: Float[np.ndarray, " n"]

    @property
    def per_m0(self) -> Float[np.ndarray, " n"]:
        r"""The atom's depletion in electrons with m = 0.

        Returns
        -------
        numpy.ndarray
            :math:`\ln P_{atom} / \ln P_0`.
        """
        return np.log(self.atom) / np.log(self.m0)

    @property
    def per_shell(self) -> Float[np.ndarray, " n"]:
        r"""The atom's depletion over that of the shell's six independent electrons.

        Returns
        -------
        numpy.ndarray
            :math:`\ln P_{atom} / (2 \ln P_0 + 4 \ln P_1)`: two electrons with
            m = 0 and four with m = +1 or -1.
        """
        return np.log(self.atom) / (2 * np.log(self.m0) + 4 * np.log(self.m1))


XENON_RMT: Mapping[str, tuple[str, str]] = {
    "hartmann_2022/rmt_ground_survival_Xe": (
        "travers_2026_trecx/ground_survival_rmtpulse_m0_Xe",
        "travers_2026_trecx/ground_survival_rmtpulse_m1_Xe",
    ),
    "hartmann_2022/rmt_ground_survival_Xe@unscaled": (
        "travers_2026_trecx/ground_survival_rmtpulse_m0_unscaled_Xe",
        "travers_2026_trecx/ground_survival_rmtpulse_m1_unscaled_Xe",
    ),
}
"""Xenon's all-electron RMT calculation against one 5p electron in its field (TC-1093).

The ground-state survival at the end of RMT's run, against one electron with m = 0
and m = 1 read at the same time, at RMT's ionisation energy and (``@unscaled``) at the
potential's own, 0.18 % deeper."""


def xenon_rmt_ratios() -> list[SurvivalRatio]:
    """Compare RMT's xenon with one electron in its field and with its measure.

    For each entry of :data:`XENON_RMT`, at the one-electron runs' intensities
    (TC-1093).

    Returns
    -------
    list of SurvivalRatio
        One per entry.
    """
    by_key = _by_key()
    ratios = []
    for key, (m0, m1) in XENON_RMT.items():
        atom_intensity, atom = _numbers(by_key[key.split("@")[0]])
        values = {float(i): float(v) for i, v in zip(atom_intensity, atom, strict=True)}
        intensity, zero = _numbers(by_key[m0])
        _, one = _numbers(by_key[m1])
        ratios.append(
            SurvivalRatio(
                key,
                intensity,
                np.array([values[float(i)] for i in intensity]),
                zero,
                one,
            )
        )
    return ratios


_RMT_FIELD_PER_ROOT_INTENSITY = 0.05336
"""Hartmann et al.'s peak field per sqrt(1e14 W/cm^2), atomic units: their intensities
are labels of E0 = 0.05336 sqrt(I/1e14 W/cm^2) au (hartmann_2022; TC-1059)."""

XENON_SURVIVAL: Mapping[str, tuple[float, float, bool]] = {
    "hartmann_2022/rmt_ground_survival_Xe": (759.4e-9, 5.297e-15, True),
    "hartmann_2022/measured_neutral_abundance_Xe": (760e-9, 4.5e-15, False),
}
"""Xenon's ground-state survival after a few-cycle pulse: the dataset, its pulse's
wavelength (m) and full width at half maximum of the intensity (s), and whether its
intensities are RMT's labels of the field (TC-1093, TC-1094)."""


def xenon_survival(model: str = "default") -> list[RateComparison]:
    r"""Predict xenon's ground-state depletion after Hartmann et al.'s pulses.

    The rate's yield (:func:`pulse_yield`, the cycle-averaged rate on a Gaussian
    envelope) at each dataset's peak intensity; RMT's intensities are converted from
    its field labels (E0 = 0.05336 sqrt(I/1e14 W/cm^2) au) to the intensity whose
    field has that amplitude. The model ionises; it does not excite (TC-1093, TC-1094).

    Parameters
    ----------
    model : str, optional
        The rate model.

    Returns
    -------
    list of RateComparison
        One per entry of :data:`XENON_SURVIVAL`, quantity ``depletion``:
        :math:`-\ln P`, the survival's logarithm, of the dataset and of the model.
    """
    by_key = _by_key()
    comparisons = []
    for key, (wavelength, fwhm, labelled) in XENON_SURVIVAL.items():
        intensity, survival = _numbers(by_key[key])
        peak = intensity
        if labelled:
            field = _RMT_FIELD_PER_ROOT_INTENSITY * AU_FIELD * np.sqrt(intensity / 1e18)
            peak = field**2 * EPS0 * C / 2
        rate = ionisation_rate_model(Gas.Xe, wavelength, model)
        predicted = np.array(
            [-math.log1p(-pulse_yield(rate, float(i), fwhm)) for i in peak]
        )
        comparisons.append(
            RateComparison(
                Gas.Xe,
                key,
                "Hartmann et al. 2022",
                model,
                "depletion",
                wavelength,
                intensity,
                -np.log(survival),
                predicted,
            )
        )
    return comparisons


LONG_PULSE_YIELD = "travers_2026_trecx/ionised_gauss42_Ar_m0"
"""The one-electron calculation of argon in Wahlstrand et al.'s 42 fs pulse (TC-702)."""


def long_pulse_yield(
    model: str = "one_electron", *, carrier: bool = True
) -> RateComparison:
    """Compare argon's yield of a 42 fs pulse at 800 nm with the one-electron TDSE.

    The calculation is of one electron in the 3p, m = 0 orbital of Tong and Lin's
    potential, in the pulse of Wahlstrand et al.'s measurement (a Gaussian of 42 fs
    full width at half maximum of the intensity), and its yield is the ionised
    fraction, the bound excited states left out.

    Parameters
    ----------
    model : str, optional
        The rate model; ``one_electron`` corresponds to the calculation.
    carrier : bool, optional
        The instantaneous rate on the field (the calculation's pulse, carrier-envelope
        phase 0), or the cycle-averaged rate on the envelope.

    Returns
    -------
    RateComparison
        The yield against the peak intensity, W/m^2.
    """
    dataset = _by_key()[LONG_PULSE_YIELD]
    wavelength = float(str(dataset.parameters["wavelength_nm"])) * 1e-9
    fwhm = float(str(dataset.parameters["pulse_fs"])) * 1e-15
    rate = ionisation_rate_model(Gas.Ar, wavelength, model)
    intensity, expected = _numbers(dataset)
    predicted = np.array(
        [pulse_yield(rate, peak, fwhm, carrier=carrier) for peak in intensity]
    )
    return RateComparison(
        Gas.Ar,
        LONG_PULSE_YIELD,
        "one-electron TDSE, the ionised fraction, 42 fs",
        model,
        "yield",
        wavelength,
        intensity,
        expected,
        predicted,
    )


def _constants(dataset: Dataset) -> dict[str, float]:
    """Return the labelled constants of a model-parameter dataset.

    Parameters
    ----------
    dataset : Dataset
        An ``ionization_model_parameter`` dataset whose rows are named constants.

    Returns
    -------
    dict of str to float
        The constants by row label.
    """
    values = dataset.scaled_values()
    return {
        str(label): float(value)
        for label, value in zip(dataset.abscissae, values, strict=True)
    }


def measured_yields(
    gas: Gas | str,
    model: str = "default",
    intensities: Sequence[float] | None = None,
) -> RateComparison:
    r"""Compare the yield of a 42 fs pulse at 800 nm with Wahlstrand et al.'s.

    Wahlstrand et al. (2018) measured the electron density a pulse leaves in a thin
    gas target, from the refractive index of the plasma, against the local peak
    intensity, which they calibrated through the Kerr phase shift with a measured
    nonlinear index of argon (13 %). Their Supplemental Material gives the fits
    :math:`Y = c_1 I_0^m`, :math:`I_0` in TW/cm\ :sup:`2`.

    Parameters
    ----------
    gas : Gas or str
        Ar, Kr, Xe (:func:`optmat.gas.ionisation_rates.ionisation_rate_model`), N2 or
        O2 (:func:`optmat.gas.molecular_ionisation_rates.molecular_rate_model`).
    model : str, optional
        The rate model.
    intensities : sequence of float or None, optional
        Peak intensities, W/m\ :sup:`2`. None: ten over the range of the
        measurement as read from the paper's Fig. 3 (Ar 45-135, Kr 35-105, Xe
        18-78, N2 75-125, O2 40-105 TW/cm\ :sup:`2`; the fit's range is not
        stated).

    Returns
    -------
    RateComparison
        The yields, with the cycle-averaged rate on the envelope of a Gaussian pulse
        of 42 fs.

    Raises
    ------
    DataError
        If the gas was not measured.
    """
    species = as_gas(gas)
    ranges = {
        "Ar": (45e16, 135e16),
        "Kr": (35e16, 105e16),
        "Xe": (18e16, 78e16),
        "N2": (75e16, 125e16),
        "O2": (40e16, 105e16),
    }
    if species.value not in ranges:
        msg = f"Wahlstrand et al. (2018) fitted no yield of {species.value}"
        raise DataError(msg)
    by_key = _by_key()
    key = f"wahlstrand_2018/yield_fit_c1_{species.value}"
    dataset = by_key[key]
    # c_1 is printed in (cm^2/TW)^m: the fit takes the intensity in TW/cm^2.
    (scale,) = dataset.scaled_values()
    (power,) = by_key[f"wahlstrand_2018/yield_fit_m_{species.value}"].scaled_values()
    if intensities is None:
        peaks = np.geomspace(*ranges[species.value], 10)
    else:
        peaks = np.asarray(intensities, dtype=float)
    wavelength = float(str(dataset.parameters["wavelength_nm"])) * 1e-9
    fwhm = float(str(dataset.parameters["pulse_fs"])) * 1e-15
    rate = (
        molecular_rate_model(species, wavelength, model)
        if species.value in {"N2", "O2"}
        else ionisation_rate_model(species, wavelength, model)
    )
    predicted = np.array([pulse_yield(rate, peak, fwhm) for peak in peaks])
    expected = float(scale) * (peaks / 1e16) ** float(power)
    return RateComparison(
        species,
        key,
        "Wahlstrand et al. 2018, measured, 42 fs",
        model,
        "yield",
        wavelength,
        peaks,
        expected,
        predicted,
    )


def wallace_yield(gas: Gas | str, intensity: ArrayLike) -> Float[np.ndarray, " *shape"]:
    r"""Return Wallace et al.'s fit of their measured yield, on its free scale.

    Wallace et al. (2016), eq. (4):

    .. math::

       S(I) = \frac{\exp\left[-\alpha\, (I/I_c)^{-1/2}\right]}{1 + (I/I_c)^{\gamma}},

    for a 5.5 fs pulse at 800 nm and a target thin against the Rayleigh range
    (:func:`focal_average`). :math:`\alpha` and :math:`\gamma` were fitted to the
    authors' one-electron calculations and :math:`I_c` to the measurement, whose
    intensities are fixed by the yield of atomic hydrogen to 1.3 % (Ar), 1.5 % (Kr)
    and 2.5 % (Xe). The scale of the yield is free. Valid from 1e14 to 5e14
    W/cm\ :sup:`2`.

    Parameters
    ----------
    gas : Gas or str
        Ar, Kr or Xe.
    intensity : array_like
        Peak intensity, W/m\ :sup:`2`.

    Returns
    -------
    numpy.ndarray
        :math:`S`, with the shape of ``intensity``.

    Raises
    ------
    DataError
        If the gas was not measured.
    """
    species = as_gas(gas)
    by_key = _by_key()
    shape_key = f"wallace_2016/tableI_shape_{species.value}"
    if shape_key not in by_key:
        msg = f"Wallace et al. (2016) measured no yield of {species.value}"
        raise DataError(msg)
    shape = _constants(by_key[shape_key])
    # I_c is curated in W/cm^2.
    (critical,) = by_key[f"wallace_2016/tableI_Ic_{species.value}"].scaled_values()
    x = np.asarray(intensity, dtype=float) / (float(critical) * 1e4)
    return np.asarray(np.exp(-shape["alpha"] / np.sqrt(x)) / (1 + x ** shape["gamma"]))


def _focal_yield(
    gas: Gas, model: str, wavelength: float, fwhm: float
) -> tuple[Float[np.ndarray, " n"], Float[np.ndarray, " n"]]:
    """Return a model's yield of a short pulse, averaged over phase and focus.

    Parameters
    ----------
    gas : Gas
        The gas.
    model : str
        The rate model.
    wavelength : float
        Vacuum wavelength, m.
    fwhm : float
        Full width at half maximum of the intensity, s.

    Returns
    -------
    tuple of numpy.ndarray
        The peak intensities, W/m^2, from 2e16 to 2e19, and the decimal logarithm of
        the yield there: the instantaneous rate on the carrier, averaged over
        carrier-envelope phases of 0 and pi/2 and over the two transverse dimensions
        of the focus (:func:`focal_average`).
    """
    rate = ionisation_rate_model(gas, wavelength, model)
    # Into saturation, where the ionised fraction is 1 whatever the rate.
    grid = np.geomspace(2e16, 2e19, 160)
    local = np.array(
        [
            0.5
            * (
                pulse_yield(rate, peak, fwhm, carrier=True)
                + pulse_yield(rate, peak, fwhm, carrier=True, cep=math.pi / 2)
            )
            for peak in grid
        ]
    )
    return grid, np.log10(focal_average(grid, local))


@dataclass(frozen=True)
class ScaleTest:
    """How far a model's intensities must be scaled to match a measured yield.

    Attributes
    ----------
    gas : Gas
        The gas.
    model : str
        The rate model.
    scale : float
        The factor on the measurement's intensities at which the model's yield has
        the measurement's dependence on intensity. Below 1: the model ionises too
        much.
    calibration : float
        The relative standard uncertainty of the measurement's intensities.
    residual : float
        Root mean square of the decimal logarithm of model over measurement at that
        scale, after the free amplitude is removed.
    intensity : numpy.ndarray
        The peak intensities compared at, W/m^2.
    measured : numpy.ndarray
        The measured yield there, on its free scale.
    predicted : numpy.ndarray
        The model's focal-averaged yield at the scaled intensities, times the
        amplitude that fits it to the measurement.
    unscaled : numpy.ndarray
        The same at the measurement's own intensities (a scale of 1).
    """

    gas: Gas
    model: str
    scale: float
    calibration: float
    residual: float
    intensity: Float[np.ndarray, " n"]
    measured: Float[np.ndarray, " n"]
    predicted: Float[np.ndarray, " n"]
    unscaled: Float[np.ndarray, " n"]


def intensity_scale_test(gas: Gas | str, model: str = "default") -> ScaleTest:
    r"""Test a model on Wallace et al.'s yields, whose scale is free.

    The model's yield for a 5.5 fs pulse at 800 nm (the instantaneous rate on the
    carrier, averaged over carrier-envelope phases of 0 and :math:`\pi/2`, the
    measurement's phase being random) is averaged over the focus
    (:func:`focal_average`). With the amplitude free, what the measurement tests is
    where the yield bends over towards saturation, which depends on the absolute
    rate: the scale :math:`\eta` is found at which :math:`S_\text{model}(\eta I)`
    has the shape of the measured :math:`S(I)` between 1e14 and 5e14
    W/cm\ :sup:`2`. One per cent in :math:`\eta` is about 5 % in the rate.

    Parameters
    ----------
    gas : Gas or str
        Ar, Kr or Xe.
    model : str, optional
        The rate model.

    Returns
    -------
    ScaleTest
        The scale and the comparison at it.
    """
    species = as_gas(gas)
    by_key = _by_key()
    critical = by_key[f"wallace_2016/tableI_Ic_{species.value}"]
    calibration = float(
        critical.standard_uncertainty()[0] / critical.scaled_values()[0]
    )
    wavelength = float(str(critical.parameters["wavelength_nm"])) * 1e-9
    fwhm = float(str(critical.parameters["pulse_fs"])) * 1e-15
    grid, averaged = _focal_yield(species, model, wavelength, fwhm)
    return _scale_test(species, model, calibration, grid, averaged)


def _scale_test(
    species: Gas,
    model: str,
    calibration: float,
    grid: Float[np.ndarray, " m"],
    averaged: Float[np.ndarray, " m"],
) -> ScaleTest:
    """Fit the intensity scale of a focal-averaged yield to Wallace et al.'s.

    Parameters
    ----------
    species : Gas
        The gas.
    model : str
        What the yield is of.
    calibration : float
        The measurement's relative uncertainty in the intensity.
    grid : numpy.ndarray
        Peak intensities, W/m^2.
    averaged : numpy.ndarray
        The decimal logarithm of the focal-averaged yield there.

    Returns
    -------
    ScaleTest
        The scale and the comparison at it.
    """
    intensity = np.geomspace(1e18, 5e18, 40)
    measured = wallace_yield(species, intensity)

    def residual(scale: float) -> Float[np.ndarray, " n"]:
        """Return the decimal logarithm of model over measurement at a scale.

        Parameters
        ----------
        scale : float
            Factor on the measurement's intensities.

        Returns
        -------
        numpy.ndarray
            The logarithm at each intensity, the free amplitude not yet removed.
        """
        model_yield = np.interp(np.log(scale * intensity), np.log(grid), averaged)
        return model_yield - np.log10(measured)

    scales = np.geomspace(0.5, 2.0, 601)
    spread = [float(np.std(residual(scale))) for scale in scales]
    best = int(np.argmin(spread))
    difference = residual(float(scales[best]))
    return ScaleTest(
        species,
        model,
        float(scales[best]),
        calibration,
        spread[best],
        intensity,
        measured,
        measured * 10 ** (difference - np.mean(difference)),
        measured * 10 ** (residual(1.0) - np.mean(residual(1.0))),
    )


WALLACE_TDSE: Mapping[str, tuple[str, str]] = {
    "Ar": (
        "travers_2026_trecx/ionised_wallace_cep0_Ar_m0",
        "travers_2026_trecx/ionised_wallace_cep90_Ar_m0",
    ),
    "Kr": (
        "travers_2026_trecx/ionised_wallace_cep0_Kr_m0",
        "travers_2026_trecx/ionised_wallace_cep90_Kr_m0",
    ),
    "Xe": (
        "travers_2026_trecx/ionised_wallace_cep0_Xe_m0",
        "travers_2026_trecx/ionised_wallace_cep90_Xe_m0",
    ),
}
"""The one-electron calculations in Wallace et al.'s pulse, at phases 0 and pi/2."""


def tdse_intensity_scale_test(gas: Gas | str = "Ar", factor: float = 1.0) -> ScaleTest:
    r"""Put the one-electron calculation in Wallace et al.'s pulse through their test.

    As :func:`intensity_scale_test`, with the calculation's ionised fraction
    (:data:`WALLACE_TDSE`, the :math:`m = 0` electron, averaged over its phases 0 and
    :math:`\pi/2`) in place of a model's: :math:`-\ln(1 - P)` is interpolated
    linearly in the logarithms between the calculated intensities, extrapolated with
    the end slopes beyond them, and multiplied by ``factor``, an atom's rate over one
    electron's (TC-836).

    Parameters
    ----------
    gas : Gas or str, optional
        A gas with such a calculation.
    factor : float, optional
        The factor on the calculation's rate.

    Returns
    -------
    ScaleTest
        The scale and the comparison at it; its ``model`` names the calculation.

    Raises
    ------
    DataError
        If no such calculation is curated for the gas.
    """
    species = as_gas(gas)
    if species.value not in WALLACE_TDSE:
        msg = (
            f"no one-electron calculation in Wallace et al.'s pulse for {species.value}"
        )
        raise DataError(msg)
    by_key = _by_key()
    critical = by_key[f"wallace_2016/tableI_Ic_{species.value}"]
    calibration = float(
        critical.standard_uncertainty()[0] / critical.scaled_values()[0]
    )
    grid = np.geomspace(2e16, 2e19, 160)
    local = np.zeros_like(grid)
    for key in WALLACE_TDSE[species.value]:
        intensity, ionised = _numbers(by_key[key])
        x, y = np.log(intensity), np.log(-np.log1p(-ionised))
        at = np.log(grid)
        inside = np.interp(at, x, y)
        low = y[0] + (y[1] - y[0]) / (x[1] - x[0]) * (at - x[0])
        high = y[-1] + (y[-1] - y[-2]) / (x[-1] - x[-2]) * (at - x[-1])
        log_depletion = np.where(at < x[0], low, np.where(at > x[-1], high, inside))
        local += 0.5 * -np.expm1(-factor * np.exp(log_depletion))
    averaged = np.log10(focal_average(grid, local))
    label = f"one-electron TDSE x {factor:g}"
    return _scale_test(species, label, calibration, grid, averaged)


@dataclass(frozen=True)
class PointsTest:
    """A model against the measured points of a yield whose scale is free.

    Attributes
    ----------
    gas : Gas
        The gas.
    model : str
        The rate model.
    scale : float
        The factor on the measurement's intensities at which the model fits the
        points best. Below 1: the model ionises too much.
    calibration : float
        The relative standard uncertainty of the measurement's intensities.
    chi_square : float
        The sum of the squared residuals over the printed errors at that scale, the
        amplitude fitted.
    chi_square_unscaled : float
        The same at the measurement's own intensities (a scale of 1).
    spread : float
        Root mean square of model over measurement less 1, at the best scale.
    intensity : numpy.ndarray
        The measurement's peak intensities, W/m^2.
    measured : numpy.ndarray
        The measured yields, on their free scale.
    sigma : numpy.ndarray
        Their printed errors.
    predicted : numpy.ndarray
        The model's yield at the scaled intensities, times the fitted amplitude.
    unscaled : numpy.ndarray
        The same at a scale of 1.
    """

    gas: Gas
    model: str
    scale: float
    calibration: float
    chi_square: float
    chi_square_unscaled: float
    spread: float
    intensity: Float[np.ndarray, " n"]
    measured: Float[np.ndarray, " n"]
    sigma: Float[np.ndarray, " n"]
    predicted: Float[np.ndarray, " n"]
    unscaled: Float[np.ndarray, " n"]


def measured_points_test(gas: Gas | str = "Kr", model: str = "default") -> PointsTest:
    r"""Test a model on the measured points of Wallace et al.'s krypton yield.

    The supplement of Wallace et al. (2016) tabulates the 37 measured yields of
    Kr\ :sup:`+` behind the fit of :func:`wallace_yield`, with their errors. Unlike
    the fit, whose shape is that of the authors' own calculation, the points are the
    measurement alone. The model's yield (as in :func:`intensity_scale_test`) is
    fitted to them with its amplitude free, the points weighted by the inverse
    square of their errors, at each scale :math:`\eta` of the intensities; the scale
    of the least :math:`\chi^2` is returned.

    Parameters
    ----------
    gas : Gas or str, optional
        Kr: the one gas whose points are tabulated.
    model : str, optional
        The rate model.

    Returns
    -------
    PointsTest
        The scale and the comparison at it.

    Raises
    ------
    DataError
        If no measured points are curated for the gas.

    Examples
    --------
    >>> test = measured_points_test("Kr")
    >>> len(test.measured)
    37
    >>> bool(abs(test.scale - 1) < 0.03)
    True
    """
    species = as_gas(gas)
    by_key = _by_key()
    key = f"wallace_2016/supplement_tableI_yield_{species.value}"
    if key not in by_key:
        msg = f"no measured points of the yield of {species.value} are curated"
        raise DataError(msg)
    points = by_key[key]
    critical = by_key[f"wallace_2016/tableI_Ic_{species.value}"]
    calibration = float(
        critical.standard_uncertainty()[0] / critical.scaled_values()[0]
    )
    wavelength = float(str(points.parameters["wavelength_nm"])) * 1e-9
    fwhm = float(str(points.parameters["pulse_fs"])) * 1e-15
    grid, averaged = _focal_yield(species, model, wavelength, fwhm)
    intensity, measured = _numbers(points)
    sigma = np.asarray(points.standard_uncertainty(), dtype=float)
    weight = 1 / sigma**2

    def fitted(scale: float) -> Float[np.ndarray, " n"]:
        """Return the model at a scale, with the amplitude of least squares.

        Parameters
        ----------
        scale : float
            Factor on the measurement's intensities.

        Returns
        -------
        numpy.ndarray
            The model's yield at each point, times the fitted amplitude.
        """
        shape = 10 ** np.interp(np.log(scale * intensity), np.log(grid), averaged)
        amplitude = np.sum(weight * measured * shape) / np.sum(weight * shape**2)
        return np.asarray(amplitude * shape)

    def chi_square(scale: float) -> float:
        """Return the weighted sum of squared residuals at a scale.

        Parameters
        ----------
        scale : float
            Factor on the measurement's intensities.

        Returns
        -------
        float
            The sum.
        """
        return float(np.sum(weight * (measured - fitted(scale)) ** 2))

    scales = np.geomspace(0.7, 1.4, 701)
    best = float(scales[int(np.argmin([chi_square(scale) for scale in scales]))])
    predicted = fitted(best)
    return PointsTest(
        species,
        model,
        best,
        calibration,
        chi_square(best),
        chi_square(1.0),
        float(np.sqrt(np.mean((predicted / measured - 1) ** 2))),
        intensity,
        measured,
        sigma,
        predicted,
        fitted(1.0),
    )


# Plasma interferometry: the phase a probe picks up from the electrons a pump leaves.

_NOBLE_GASES = frozenset({"He", "Ne", "Ar", "Kr", "Xe"})

_DENSITY_PER_ATM = 2.7e25
"""Atoms per m^3 at the 1 atm the phase shifts are normalised to, as Nie et al. (2022)
state it (Fig. 4: '500 rad corresponds to 100% ionization ... (plasma density of 2.7 x
10^19 cm^-3)'): the density their comparison assumes (TC-817)."""

_LAGUERRE = np.polynomial.laguerre.laggauss(64)
"""Nodes and weights of the radial average over the probe's Gaussian profile."""

_AXIAL = np.sinh(np.linspace(-math.asinh(40.0), math.asinh(40.0), 801))
"""Positions along the axis in Rayleigh lengths of the pump, dense near the focus."""


@dataclass(frozen=True)
class PhaseShiftTest:
    """A model against a measured phase shift of a probe through a laser plasma.

    Attributes
    ----------
    gas : Gas
        The gas.
    reference : str
        The key of the measured dataset.
    model : str
        The rate model.
    polarization : str
        The pump's polarisation, linear or circular.
    intensity : numpy.ndarray
        The measurement's peak intensities, W/m^2.
    measured : numpy.ndarray
        The measured phase shift, rad, at the pressure the dataset is normalised to.
    sigma : numpy.ndarray
        Its uncertainty, rad; NaN where none is drawn.
    predicted : numpy.ndarray
        The model's phase shift at the measurement's intensities, rad.
    scale : float
        The factor on the measurement's intensities at which the model's phase shift
        fits the measured one best (in the logarithm, weighted by the error bars).
        Below 1: the model ionises too much.
    rescaled : numpy.ndarray
        The model's phase shift at the scaled intensities, rad.
    """

    gas: Gas
    reference: str
    model: str
    polarization: str
    intensity: Float[np.ndarray, " n"]
    measured: Float[np.ndarray, " n"]
    sigma: Float[np.ndarray, " n"]
    predicted: Float[np.ndarray, " n"]
    scale: float
    rescaled: Float[np.ndarray, " n"]

    @property
    def ratio(self) -> Float[np.ndarray, " n"]:
        """Model over measurement at the measurement's own intensities.

        Returns
        -------
        numpy.ndarray
            The ratio of the phase shifts at each intensity.
        """
        return self.predicted / self.measured


def _parameter(dataset: Dataset, name: str) -> float:
    """Return a numeric parameter of a dataset.

    Parameters
    ----------
    dataset : Dataset
        The dataset.
    name : str
        The parameter.

    Returns
    -------
    float
        Its value.
    """
    return float(str(dataset.parameters[name]))


def circular_yield(
    model: IonisationRateModel,
    peak_intensity: float,
    fwhm: float,
) -> float:
    r"""Return the fraction a circularly polarised Gaussian pulse ionises.

    The field of circular polarisation turns at a constant magnitude
    :math:`F = \sqrt{I/\varepsilon_0 c}` within a cycle, so the rate,
    :meth:`~optmat.gas.ionisation_rates.IonisationRateModel.circular_rate`, is
    followed along the envelope (TC-820).

    Parameters
    ----------
    model : IonisationRateModel
        The rate, at the pulse's wavelength.
    peak_intensity : float
        Peak intensity, W/m^2.
    fwhm : float
        Full width at half maximum of the intensity, s.

    Returns
    -------
    float
        The ionised fraction after the pulse.
    """
    (fraction,) = _envelope_yields(
        model, np.array([peak_intensity]), fwhm, circular=True
    )
    return float(fraction)


def _envelope_yields(
    model: IonisationRateModel | MolecularRateModel,
    peak_intensity: Float[np.ndarray, " n"],
    fwhm: float,
    *,
    circular: bool = False,
) -> Float[np.ndarray, " n"]:
    """Return the ionised fraction of Gaussian pulses, the rate taken on the envelope.

    Linear polarisation: the cycle-averaged rate at the envelope's amplitude, as
    :func:`pulse_yield` without ``carrier``. Circular: the rate at the magnitude of
    the turning field, the amplitude over the square root of 2. The envelope is
    sampled at 1201 points over three widths either side of the peak.

    Parameters
    ----------
    model : IonisationRateModel or MolecularRateModel
        The rate.
    peak_intensity : numpy.ndarray
        Peak intensities, W/m^2.
    fwhm : float
        Full width at half maximum of the intensity, s.
    circular : bool, optional
        Circular polarisation.

    Returns
    -------
    numpy.ndarray
        The ionised fraction after each pulse.
    """
    times = np.linspace(-_PULSE_HALF_WIDTH, _PULSE_HALF_WIDTH, 1201) * fwhm
    shape = np.exp(-2 * math.log(2) * (times / fwhm) ** 2)
    field = amplitude(peak_intensity)[:, None] * shape[None, :]
    if circular:
        if not isinstance(model, IonisationRateModel):
            msg = "a circular-polarisation rate is given for atoms only"
            raise DataError(msg)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ModelRangeWarning)
            rate = model.circular_rate(field / math.sqrt(2.0))
    else:
        rate = _evaluate(model, field, cycle_averaged=True)
    return -np.expm1(-np.trapezoid(rate, times, axis=1))


def plasma_phase_shift(
    ionised: Float[np.ndarray, " m"],
    grid: Float[np.ndarray, " m"],
    peak_intensity: ArrayLike,
    parameters: Mapping[str, float],
) -> Float[np.ndarray, " *shape"]:
    r"""Return the phase a probe picks up from the electrons a focused pump leaves.

    Both beams are Gaussian with the same focus and beam quality :math:`M^2`, the pump
    :math:`I(r, z) = I_0 (w_0/w)^2 \exp(-2 r^2/w^2)`, :math:`w = w_0 \sqrt{1 +
    (z/z_R)^2}`, :math:`z_R = \pi w_0^2 / M^2 \lambda`. The ionised fraction at each
    point is that of the local peak intensity. The phase is Nie et al.'s (2022) Eq.
    (3), the change of index :math:`n_e/2 n_c` averaged over the probe's profile at
    each z and integrated along the axis:

    .. math:: \Delta\phi = \frac{2\pi}{\lambda_p} \frac{n}{2 n_c} \int dz\,
       \frac{\int Y(I(r, z))\, I_p(r, z)\, r\, dr}{\int I_p(r, z)\, r\, dr},

    with :math:`n` the density of atoms and :math:`n_c = \varepsilon_0 m_e
    \omega_p^2 / e^2` the critical density at the probe's frequency. Self-focusing is
    left out, as the authors do at 800 and 400 nm.

    Parameters
    ----------
    ionised : numpy.ndarray
        The ionised fraction after the pulse, at the peak intensities ``grid``.
    grid : numpy.ndarray
        Increasing peak intensities, W/m^2; the fraction is interpolated in the
        logarithm between them and is 0 below the first.
    peak_intensity : array_like
        The pump's peak intensity at the focus, W/m^2.
    parameters : mapping
        ``wavelength_nm`` and ``probe_wavelength_nm``; ``pump_waist_um`` and
        ``probe_waist_um``, the 1/e^2 radii of the intensity at the focus; ``m2``;
        and ``density_m3``, the density of atoms.

    Returns
    -------
    numpy.ndarray
        The phase shift, rad, with the shape of ``peak_intensity``.

    Examples
    --------
    Every atom ionised everywhere: the phase of a column of electrons 80 Rayleigh
    lengths long (the axis the average covers):

    >>> beams = {"wavelength_nm": 800, "probe_wavelength_nm": 400, "m2": 1.25,
    ...          "pump_waist_um": 60.0, "probe_waist_um": 42.0, "density_m3": 2.7e25}
    >>> phase = plasma_phase_shift(np.ones(2), np.array([1e10, 1e30]), 1e18, beams)
    >>> length = 80 * np.pi * 60e-6**2 / (1.25 * 800e-9)
    >>> f"{float(phase) / length:.4e}"
    '3.0434e+04'
    """
    probe = parameters["probe_wavelength_nm"] * 1e-9
    m2 = parameters["m2"]
    pump_waist = parameters["pump_waist_um"] * 1e-6
    probe_waist = parameters["probe_waist_um"] * 1e-6
    rayleigh = math.pi * pump_waist**2 / (m2 * parameters["wavelength_nm"] * 1e-9)
    probe_rayleigh = math.pi * probe_waist**2 / (m2 * probe)
    z = _AXIAL * rayleigh
    pump = pump_waist**2 * (1 + (z / rayleigh) ** 2)  # w^2
    beam = probe_waist**2 * (1 + (z / probe_rayleigh) ** 2)
    nodes, weights = _LAGUERRE
    # With s = 2 r^2 / w_probe^2 the probe's weight is exp(-s) ds.
    local = (pump_waist**2 / pump)[:, None] * np.exp(
        -nodes[None, :] * (beam / pump)[:, None]
    )
    log_grid = np.log(grid)
    log_ionised = np.log(np.maximum(ionised, 1e-300))

    def column(peak: float) -> float:
        """Return the probe-averaged ionised fraction integrated along the axis.

        Parameters
        ----------
        peak : float
            The pump's peak intensity at the focus, W/m^2.

        Returns
        -------
        float
            The integral, m.
        """
        intensity = peak * local
        fraction = np.where(
            intensity > grid[0],
            np.exp(
                np.interp(np.log(np.maximum(intensity, grid[0])), log_grid, log_ionised)
            ),
            0.0,
        )
        return float(np.trapezoid(fraction @ weights, z))

    omega = 2 * math.pi * C / probe
    critical = EPS0 * M_E * omega**2 / E_CHARGE**2
    factor = 2 * math.pi / probe * parameters["density_m3"] / (2 * critical)
    peaks = np.asarray(peak_intensity, dtype=float)
    columns = np.vectorize(column, otypes=[float])(peaks)
    return np.asarray(factor * columns)


def phase_shift_test(reference: str, model: str = "default") -> PhaseShiftTest:
    """Compare a model with a measured phase shift through a laser plasma.

    The model's ionised fraction for the measurement's pulse (a Gaussian of the
    stated duration; linear polarisation: the cycle-averaged rate on the envelope;
    circular: :func:`circular_yield`) is integrated over the
    focus by :func:`plasma_phase_shift` with the beams the dataset states.

    Parameters
    ----------
    reference : str
        The key of an ``ionization_phase_shift`` dataset.
    model : str, optional
        The rate model.

    Returns
    -------
    PhaseShiftTest
        The comparison.

    Raises
    ------
    DataError
        If the dataset is not a phase shift, or does not state its beams.
    """
    dataset = _by_key()[reference]
    if dataset.quantity != "ionization_phase_shift":
        msg = f"{reference} is not an ionization_phase_shift dataset"
        raise DataError(msg)
    species = as_gas(dataset.material)
    wavelength = _parameter(dataset, "wavelength_nm") * 1e-9
    fwhm = _parameter(dataset, "pulse_fs") * 1e-15
    polarization = str(dataset.parameters["polarization"])
    factory = (
        ionisation_rate_model if species.value in _NOBLE_GASES else molecular_rate_model
    )
    names = (
        "wavelength_nm",
        "probe_wavelength_nm",
        "pump_waist_um",
        "probe_waist_um",
        "m2",
    )
    missing = [name for name in names if name not in dataset.parameters]
    if missing:
        msg = f"{reference} does not state {', '.join(missing)}: the focus is unknown"
        raise DataError(msg)
    beams = {name: _parameter(dataset, name) for name in names}
    pressure = _parameter(dataset, "pressure_atm")
    stated = dataset.parameters.get("density_per_atm_cm3")
    per_atm = _DENSITY_PER_ATM if stated is None else float(str(stated)) * 1e6
    beams["density_m3"] = per_atm * pressure
    intensity = dataset.abscissa_values()
    measured = dataset.scaled_values()
    sigma = (
        np.array(
            [float(s) if s.strip() else math.nan for s in dataset.sigmas], dtype=float
        )
        * dataset.scale
    )
    grid = np.geomspace(
        float(np.min(intensity)) / 1e3, float(np.max(intensity)) * 2, 160
    )
    if polarization == "circular":
        ionised = _envelope_yields(
            factory(species, wavelength, model), grid, fwhm, circular=True
        )
    else:
        ionised = _envelope_yields(factory(species, wavelength, model), grid, fwhm)
    # The phase on a grid of peak intensities, interpolated in the logarithm: the
    # search for the scale then costs no further integrals over the focus.
    peaks = np.geomspace(
        float(np.min(intensity)) / 2.5, float(np.max(intensity)) * 2, 241
    )
    log_phase = np.log(plasma_phase_shift(ionised, grid, peaks, beams))

    def phase(at: Float[np.ndarray, " n"]) -> Float[np.ndarray, " n"]:
        """Return the model's phase shift at peak intensities.

        Parameters
        ----------
        at : numpy.ndarray
            Peak intensities, W/m^2.

        Returns
        -------
        numpy.ndarray
            The phase shift, rad.
        """
        return np.exp(np.interp(np.log(at), np.log(peaks), log_phase))

    predicted = phase(intensity)
    log_sigma = np.where(np.isfinite(sigma), sigma / measured, 0.1)
    scales = np.geomspace(0.5, 2.0, 601)
    misfit = [
        float(np.sum((np.log(phase(s * intensity) / measured) / log_sigma) ** 2))
        for s in scales
    ]
    best = float(scales[int(np.argmin(misfit))])
    return PhaseShiftTest(
        species,
        reference,
        model,
        polarization,
        intensity,
        measured,
        sigma,
        predicted,
        best,
        phase(best * intensity),
    )
