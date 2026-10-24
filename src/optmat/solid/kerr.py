# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The Kerr response of a solid: the nonlinear index, its dispersion, its Raman part.

A solid's intensity-dependent index n = n0 + n2 I has an instantaneous electronic part
and a delayed nuclear (Raman) part. A measurement with pulses of length tau sees the
Raman part only in so far as the response can follow the pulse:

.. math::

    n_{2,\mathrm{meas}}(\lambda, \tau) = n_{2}(\lambda)
        \left[(1 - f_R) + f_R\,\eta(\tau)\right],
    \qquad
    \eta(\tau) = \frac{\int I(t) \int_0^\infty h_R(s)\, I(t - s)\, ds\, dt}
                        {\int I(t)^2\, dt},

where n2(lambda) is the long-pulse (total) index, f_R its Raman fraction and h_R the
normalised Raman response (Hollenbeck & Cantrell 2002, ``hollenbeck_2002``). eta is 1
for pulses much longer than the response (about 100 fs), 0 in the limit of very short
ones. Electrostriction (nanosecond pulses and longer) is not modelled; such
measurements enter only as validators.

The Raman fraction is constrained by measured steady-state Raman gains, predicted as
g(Omega) = 2 omega_s n2(lambda_p) f_R Im H_R(Omega) / c with Hollenbeck & Cantrell's
H_R (``optmat.solid.raman_response``): the gain fixes f_R n2 at the pump, and so pairs
f_R with that shape (TC-696). A prior on f_R may be set instead, or as well.

The dispersion of n2 below the two-photon absorption edge is modelled as an even
polynomial in the photon energy, n2(lambda) = N [1 + a E^2 + b E^4] (E in eV), the form
of the gases' Taylor-series dispersion (``optmat.gas.nonlinear``). The fit
(``data/solid/kerr/<name>.toml``, ``tools/fit/kerr_solid.py``) gives each source a
calibration scale with a prior. Its result is stored with its covariance (TC-393).

Values printed in esu, n = n0 + n2 <E^2> (``nonlinear_index_field``), are converted
with I = c n0 <E^2> / (4 pi) (Gaussian units): n2 [cm^2/W] = 4 pi 1e7 n2 [esu] /
(c [cm/s] n0), with n0 from the material's recommended linear model.
"""

import functools
import hashlib
import json
import math
import tomllib
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar, Final

import numpy as np
from jaxtyping import Float
from numpy.typing import ArrayLike

from optmat import _catalogue
from optmat.constants import E_CHARGE, HBAR, C
from optmat.datasets import Dataset
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.identity import Identified
from optmat.solid import SOLID_DATA_DIRECTORY
from optmat.solid.models.registry import default_model

if TYPE_CHECKING:
    from scipy.optimize import OptimizeResult

__all__ = [
    "KERR_DIRECTORY",
    "GainOptions",
    "GainPoints",
    "KerrConfig",
    "KerrModel",
    "KerrOptions",
    "KerrResult",
    "Points",
    "dumps",
    "fit",
    "gain_points",
    "kerr_models",
    "load_config",
    "points",
    "raman_response",
    "raman_weight",
]

KERR_DIRECTORY: Final = SOLID_DATA_DIRECTORY / "kerr"
"""The Kerr fit configurations and their stored results."""

_CM_S: Final = C * 100.0

# The Raman response's transform is taken on the frequency grid of a 4 ps span (it has
# decayed to below 1e-6 of its peak by 2 ps) sampled at 0.1 fs, a hundredth of its
# fastest period.
_RESPONSE_STEP: Final = 0.1e-15
_RESPONSE_SPAN: Final = 4e-12


def raman_response(time: Float[np.ndarray, " n"]) -> Float[np.ndarray, " n"]:
    """Return fused silica's Raman response h_R(t), normalised to unit integral.

    Hollenbeck & Cantrell's intermediate-broadening model (J. Opt. Soc. Am. B 19, 2886
    (2002), Eq. (9) and Table 1; ``hollenbeck_2002``), whose normalisation the paper
    does not state: :meth:`optmat.solid.raman_response.SilicaRamanResponse.response`
    of its ``hollenbeck_2002`` model, normalised analytically (TC-691).

    Parameters
    ----------
    time : numpy.ndarray
        Time after the impulse, s (the response is zero before it).

    Returns
    -------
    numpy.ndarray
        The response h_R, 1/s.

    Examples
    --------
    >>> t = np.arange(0, 4e-12, 0.1e-15)
    >>> f"{np.trapezoid(raman_response(t), t):.4f}"
    '1.0000'
    """
    from optmat.solid.raman_response import silica_raman_responses  # noqa: PLC0415

    return silica_raman_responses()["hollenbeck_2002"].response(time)


@functools.cache
def _response_spectrum() -> tuple[np.ndarray, np.ndarray]:
    """Return Re H(omega), the real part of the Raman response's Fourier transform.

    Returns
    -------
    tuple of numpy.ndarray
        Angular frequencies (rad/s, from 0) and Re H, with H(0) = 1.
    """
    from optmat.solid.raman_response import silica_raman_responses  # noqa: PLC0415

    # The analytic transform (Faddeeva), on the grid of a 4 ps span: exact at every
    # frequency, and Re H(0) = 1 (TC-691; the sampled transform it replaces was 4e-6
    # low at zero frequency, the trapezoid's error on 0.1 fs).
    n = round(_RESPONSE_SPAN / _RESPONSE_STEP)
    omega = 2 * math.pi * np.fft.rfftfreq(n, _RESPONSE_STEP)
    spectrum = silica_raman_responses()["hollenbeck_2002"].spectrum(omega)
    return omega, np.real(spectrum)


@functools.cache
def raman_weight(pulse_fwhm: float) -> float:
    """Return eta, the share of the Raman response a Gaussian pulse's index sees.

    The intensity-weighted mean of the delayed response over the pulse, relative to a
    response that follows the pulse instantly (module docstring). For a Gaussian pulse
    of rms duration s it is, exactly,
    eta = int exp(-s^2 w^2) Re H(w) dw / int exp(-s^2 w^2) dw, with H the Fourier
    transform of h_R (H(0) = 1).

    Parameters
    ----------
    pulse_fwhm : float
        The intensity FWHM of a Gaussian pulse, s; ``inf`` for a long pulse.

    Returns
    -------
    float
        The weight eta, between 0 and 1.

    Examples
    --------
    >>> f"{raman_weight(math.inf):.1f}"
    '1.0'
    >>> f"{raman_weight(10e-12):.3f}"
    '1.000'
    """
    if math.isinf(pulse_fwhm):
        return 1.0
    sigma = pulse_fwhm / (2 * math.sqrt(2 * math.log(2)))
    omega, real = _response_spectrum()
    weight = np.exp(-((sigma * omega) ** 2))
    if weight[1] < 1e-12 or np.sum(weight > 1e-12) < 8:
        # The pulse's spectrum is narrower than the grid resolves: eta = Re H(0) = 1.
        return float(real[0])
    # Re H and the weight are even in omega: the half-line sums suffice.
    return float(np.trapezoid(weight * real, omega) / np.trapezoid(weight, omega))


@dataclass(frozen=True)
class KerrOptions:
    """How one curated n2 dataset enters a Kerr fit.

    Attributes
    ----------
    key : str
        ``source/dataset``.
    group : str
        The calibration group (its own source by default): the group's values share
        one scale factor.
    scale_sigma : float
        The prior of the group's scale (relative, one sigma).
    raman : str
        ``pulse`` (eta from the stated pulse length), ``none`` (the source isolates the
        electronic part), or ``full`` (pulses much longer than the response).
    pulse_fs : float or None
        The pulse FWHM where the dataset states none, fs.
    sigma : float or None
        The relative uncertainty used where the dataset states none.
    extra : float
        A relative scatter added in quadrature.
    range_nm : tuple of float
        The wavelength window used, nm.
    """

    key: str
    group: str = ""
    scale_sigma: float = 0.1
    raman: str = "pulse"
    pulse_fs: float | None = None
    sigma: float | None = None
    extra: float = 0.0
    range_nm: tuple[float, float] = (0.0, math.inf)


@dataclass(frozen=True)
class GainOptions:
    """How one curated Raman gain enters a Kerr fit.

    Attributes
    ----------
    key : str
        ``source/dataset``: a ``raman_gain`` of the material, with its pump
        (``parameters.pump_nm``) and its Stokes shift (``parameters.shift_cm1``).
    sigma : float or None
        The relative uncertainty used where the dataset states none.
    """

    key: str
    sigma: float | None = None


@dataclass(frozen=True)
class KerrConfig:
    """A Kerr fit's configuration.

    Attributes
    ----------
    name : str
        The configuration's name.
    material : str
        The material.
    description : str
        What it is.
    range_nm : tuple of float
        The range the fitted dispersion supports, nm.
    raman_fraction : tuple of float or None
        The prior of f_R, value and sigma; None for none (the gains then set it).
    quartic : bool
        Whether the E^4 term of the dispersion is fitted (else it is zero).
    fitted, validators : tuple of KerrOptions
        The fitted datasets and those only compared.
    raman_gains : tuple of GainOptions
        The measured Raman gains fitted.
    digest : str
        The configuration's SHA-256.
    """

    name: str
    material: str
    description: str
    range_nm: tuple[float, float]
    raman_fraction: tuple[float, float] | None
    quartic: bool = True
    fitted: tuple[KerrOptions, ...] = ()
    validators: tuple[KerrOptions, ...] = ()
    raman_gains: tuple[GainOptions, ...] = ()
    digest: str = ""


def _kerr_options(entry: Mapping[str, object], where: str) -> KerrOptions:
    """Read one dataset's options.

    Parameters
    ----------
    entry : mapping
        The entry.
    where : str
        The configuration, for messages.

    Returns
    -------
    KerrOptions
        The options.

    Raises
    ------
    DataError
        If the entry has unknown keys or an unknown Raman treatment.
    """
    allowed = {
        "key",
        "group",
        "scale_sigma",
        "raman",
        "pulse_fs",
        "sigma",
        "extra",
        "range_nm",
        "note",
    }
    if set(entry) - allowed:
        msg = f"{where}: unknown keys {sorted(set(entry) - allowed)}"
        raise DataError(msg)
    raman = str(entry.get("raman", "pulse"))
    if raman not in {"pulse", "none", "full"}:
        msg = f"{where}: unknown raman treatment {raman!r}"
        raise DataError(msg)
    window = entry.get("range_nm", [0.0, math.inf])
    key = str(entry["key"])
    return KerrOptions(
        key=key,
        group=str(entry.get("group", key.split("/", maxsplit=1)[0])),
        scale_sigma=float(entry.get("scale_sigma", 0.1)),  # pyright: ignore[reportArgumentType]
        raman=raman,
        pulse_fs=None if "pulse_fs" not in entry else float(entry["pulse_fs"]),  # pyright: ignore[reportArgumentType]
        sigma=None if "sigma" not in entry else float(entry["sigma"]),  # pyright: ignore[reportArgumentType]
        extra=float(entry.get("extra", 0.0)),  # pyright: ignore[reportArgumentType]
        range_nm=(float(window[0]), float(window[1])),  # pyright: ignore[reportIndexIssue, reportArgumentType]
    )


def load_config(name: str | Path) -> KerrConfig:
    """Load a Kerr fit configuration.

    Parameters
    ----------
    name : str or pathlib.Path
        The name (``data/solid/kerr/<name>.toml``) or a path.

    Returns
    -------
    KerrConfig
        The configuration.

    Raises
    ------
    DataError
        If it is malformed or names an unknown dataset.
    """
    path = Path(name)
    if path.suffix != ".toml":
        path = KERR_DIRECTORY / f"{name}.toml"
    raw = path.read_bytes()
    record = tomllib.loads(raw.decode("utf-8"))
    allowed = {
        "material",
        "description",
        "range_nm",
        "raman_fraction",
        "quartic",
        "fitted",
        "validators",
        "raman_gains",
    }
    if set(record) - allowed:
        msg = f"{path.name}: unknown keys {sorted(set(record) - allowed)}"
        raise DataError(msg)
    fitted = tuple(_kerr_options(e, path.name) for e in record.get("fitted", []))
    validators = tuple(
        _kerr_options(e, path.name) for e in record.get("validators", [])
    )
    gains = []
    for entry in record.get("raman_gains", []):
        unknown = sorted(set(entry) - {"key", "sigma", "note"})
        if unknown:
            msg = f"{path.name}: unknown keys {unknown}"
            raise DataError(msg)
        gains.append(
            GainOptions(
                key=str(entry["key"]),
                sigma=None if "sigma" not in entry else float(entry["sigma"]),
            )
        )
    known = _catalogue.by_key()
    for options in (*fitted, *validators, *gains):
        if options.key not in known:
            msg = f"{path.name}: unknown dataset {options.key}"
            raise DataError(msg)
    low, high = record["range_nm"]
    fraction = record.get("raman_fraction")
    if fraction is None and not gains:
        msg = f"{path.name}: f_R needs a prior (raman_fraction) or raman_gains"
        raise DataError(msg)
    return KerrConfig(
        name=path.stem,
        material=str(record["material"]),
        description=str(record.get("description", "")),
        range_nm=(float(low), float(high)),
        raman_fraction=None
        if fraction is None
        else (float(fraction[0]), float(fraction[1])),
        quartic=bool(record.get("quartic", True)),
        fitted=fitted,
        validators=validators,
        raman_gains=tuple(gains),
        digest=hashlib.sha256(raw).hexdigest(),
    )


@dataclass(frozen=True)
class Points:
    """The n2 values of one dataset, in m^2/W, with what the model needs.

    Attributes
    ----------
    key : str
        The dataset.
    group : str
        Its calibration group.
    wavelength : numpy.ndarray
        Vacuum wavelengths, m.
    value : numpy.ndarray
        n2, m^2/W (self coefficients; esu values converted).
    sigma : numpy.ndarray
        Their uncertainties, m^2/W.
    eta : numpy.ndarray
        The Raman weight of each value's pulse.
    """

    key: str
    group: str
    wavelength: Float[np.ndarray, " n"]
    value: Float[np.ndarray, " n"]
    sigma: Float[np.ndarray, " n"]
    eta: Float[np.ndarray, " n"]


def _in_m2_per_w(dataset: Dataset, wavelength: np.ndarray) -> np.ndarray:
    """Return a dataset's n2 in m^2/W.

    Parameters
    ----------
    dataset : Dataset
        A `nonlinear_index` or `nonlinear_index_field` dataset.
    wavelength : numpy.ndarray
        Its vacuum wavelengths, m (for n0).

    Returns
    -------
    numpy.ndarray
        The index n2, m^2/W.

    Raises
    ------
    DataError
        If the unit or the coefficient is not one this module converts.
    """
    if dataset.parameters.get("coefficient") != "self":
        msg = f"{dataset.key}: only self coefficients are fitted"
        raise DataError(msg)
    values = dataset.scaled_values()
    if dataset.quantity == "nonlinear_index":
        return values * {"m2/W": 1.0, "cm2/W": 1e-4}[dataset.unit]
    if dataset.unit != "esu":
        msg = f"{dataset.key}: a field n2 in {dataset.unit} is not converted"
        raise DataError(msg)
    n0 = np.real(
        default_model(dataset.material).refractive_index(wavelength, warn=False)
    )
    return 4 * math.pi * 1e7 * values / (_CM_S * n0) * 1e-4


def points(options: KerrOptions) -> Points:
    """Build one dataset's points: its n2 in m^2/W with their pulses' Raman weights.

    Parameters
    ----------
    options : KerrOptions
        How it enters.

    Returns
    -------
    Points
        The values inside the window, in m^2/W.

    Raises
    ------
    DataError
        If a point lacks an uncertainty or a pulse length it needs.
    """
    dataset = _catalogue.dataset(options.key)
    wavelength = dataset.vacuum_wavelength()
    value = _in_m2_per_w(dataset, wavelength)
    stated = dataset.standard_uncertainty() / np.abs(dataset.scaled_values())
    if options.sigma is not None:
        stated = np.where(np.isfinite(stated), stated, options.sigma)
    relative = np.sqrt(stated**2 + options.extra**2)
    low, high = options.range_nm
    keep = (wavelength >= low * 1e-9) & (wavelength <= high * 1e-9)
    if options.raman == "none":
        eta = np.zeros_like(value)
    elif options.raman == "full":
        eta = np.ones_like(value)
    else:
        if "pulse_fs" in dataset.columns:
            pulse = dataset.column_values("pulse_fs")
        elif "pulse_fs" in dataset.parameters:
            pulse = np.full(len(value), float(dataset.parameters["pulse_fs"]))  # pyright: ignore[reportArgumentType]
        elif options.pulse_fs is not None:
            pulse = np.full(len(value), options.pulse_fs)
        else:
            msg = f"{dataset.key}: no pulse length; set pulse_fs or raman"
            raise DataError(msg)
        eta = np.array([raman_weight(p * 1e-15) for p in pulse])
    if not np.all(np.isfinite(relative[keep])):
        msg = f"{dataset.key}: points without an uncertainty; set sigma in the config"
        raise DataError(msg)
    return Points(
        key=dataset.key,
        group=options.group,
        wavelength=wavelength[keep],
        value=value[keep],
        sigma=(relative * np.abs(value))[keep],
        eta=eta[keep],
    )


@dataclass(frozen=True)
class GainPoints:
    """A measured Raman gain, with what the model needs to predict it.

    Attributes
    ----------
    key : str
        The dataset.
    wavelength : float
        The pump's wavelength, m.
    shift : float
        The Stokes shift, rad/s.
    value : float
        The gain, m/W.
    sigma : float
        Its uncertainty, m/W.
    response : float
        Im H_R at the shift, Hollenbeck & Cantrell's normalised response.
    """

    key: str
    wavelength: float
    shift: float
    value: float
    sigma: float
    response: float

    def factor(self) -> float:
        """Return ``2 omega_s Im H_R / c``: the gain over f_R n2(pump).

        Returns
        -------
        float
            The factor, 1/m.
        """
        omega_s = 2 * math.pi * C / self.wavelength - self.shift
        return 2 * omega_s * self.response / C


def gain_points(options: GainOptions) -> GainPoints:
    """Build a measured Raman gain for a Kerr fit.

    Parameters
    ----------
    options : GainOptions
        How it enters.

    Returns
    -------
    GainPoints
        The gain in m/W, with its pump, shift and Im H_R.

    Raises
    ------
    DataError
        If the dataset is not one gain in cm/W or m/W with its pump and shift, or has no
        uncertainty.
    """
    from optmat.solid.raman_response import silica_raman_responses  # noqa: PLC0415

    dataset = _catalogue.dataset(options.key)
    unit = {"cm/W": 1e-2, "m/W": 1.0}.get(dataset.unit)
    if dataset.quantity != "raman_gain" or unit is None or len(dataset) != 1:
        msg = f"{dataset.key}: one raman_gain in cm/W or m/W is expected"
        raise DataError(msg)
    for name in ("pump_nm", "shift_cm1"):
        if name not in dataset.parameters:
            msg = f"{dataset.key}: a fitted gain needs parameters.{name}"
            raise DataError(msg)
    value = float(dataset.scaled_values()[0]) * unit
    stated = float(dataset.standard_uncertainty()[0]) * unit
    if not math.isfinite(stated):
        if options.sigma is None:
            msg = f"{dataset.key}: no uncertainty; set sigma in the config"
            raise DataError(msg)
        stated = options.sigma * value
    shift = 2 * math.pi * _CM_S * float(str(dataset.parameters["shift_cm1"]))
    shape = silica_raman_responses()["hollenbeck_2002"]
    return GainPoints(
        key=dataset.key,
        wavelength=float(str(dataset.parameters["pump_nm"])) * 1e-9,
        shift=shift,
        value=value,
        sigma=stated,
        response=float(np.imag(shape.spectrum(shift))),
    )


def _photon_energy(wavelength: np.ndarray) -> np.ndarray:
    """Return the photon energy in eV.

    Parameters
    ----------
    wavelength : numpy.ndarray
        Vacuum wavelength, m.

    Returns
    -------
    numpy.ndarray
        The energy, eV.
    """
    return np.asarray(2 * math.pi * C * HBAR / (E_CHARGE * wavelength))


@dataclass(frozen=True)
class KerrModel(Identified):
    """A fitted Kerr response: n2(lambda) with its Raman fraction and covariance.

    Attributes
    ----------
    material : str
        The material.
    name : str
        ``kerr/<name>``.
    citation : str
        The fit.
    wavelength_range : tuple of float
        Where the dispersion is fitted, m.
    parameters : tuple of float
        N (m^2/W), a (eV^-2), b (eV^-4) and f_R.
    covariance : numpy.ndarray
        Their covariance (Birge-scaled).
    scales : dict of str to float
        The calibration scale of each group.
    status : str
        ``validated`` for the stored fit; ``user`` for a model a caller built.
    """

    IDENTITY_KIND: ClassVar[str] = "kerr"
    EVALUATION_VERSION: ClassVar[int] = 1

    material: str
    name: str
    citation: str
    wavelength_range: tuple[float, float]
    parameters: tuple[float, float, float, float]
    covariance: Float[np.ndarray, "4 4"] = field(repr=False)
    scales: Mapping[str, float] = field(default_factory=dict)

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the model's name, material and the sources of its calibrations.

        Returns
        -------
        tuple
            ``(name, material, sources)``.
        """
        return self.name, self.material, tuple(sorted(self.scales))

    @property
    def raman_fraction(self) -> float:
        """The Raman fraction f_R of the long-pulse n2.

        Returns
        -------
        float
            The fraction f_R.
        """
        return self.parameters[3]

    def _check(self, wavelength: np.ndarray) -> None:
        """Warn outside the fitted range.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.
        """
        low, high = self.wavelength_range
        if np.any(wavelength < low) or np.any(wavelength > high):
            import warnings  # noqa: PLC0415

            msg = f"{self.name} evaluated outside {low * 1e9:.4g}-{high * 1e9:.4g} nm"
            warnings.warn(msg, ModelRangeWarning, stacklevel=3)

    def n2(self, wavelength: ArrayLike | float) -> Float[np.ndarray, "*shape"]:
        """Return the long-pulse (total, electronic plus Raman) n2.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            The index n2, m^2/W.
        """
        wavelength = np.asarray(wavelength, dtype=float)
        self._check(wavelength)
        n, a, b, _ = self.parameters
        e2 = _photon_energy(wavelength) ** 2
        return np.asarray(n * (1 + a * e2 + b * e2**2))

    def n2_electronic(
        self, wavelength: ArrayLike | float
    ) -> Float[np.ndarray, "*shape"]:
        """Return the instantaneous (electronic) part, (1 - f_R) n2.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            The index n2, m^2/W.
        """
        return np.asarray((1 - self.raman_fraction) * self.n2(wavelength))

    def n2_effective(
        self, wavelength: ArrayLike | float, pulse_fwhm: float
    ) -> Float[np.ndarray, "*shape"]:
        """Return the n2 a Gaussian pulse of the given length measures.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        pulse_fwhm : float
            The pulse's intensity FWHM, s.

        Returns
        -------
        numpy.ndarray
            The index n2, m^2/W.
        """
        fraction = self.raman_fraction
        eta = raman_weight(pulse_fwhm)
        return np.asarray(self.n2(wavelength) * ((1 - fraction) + fraction * eta))

    def n2_sigma(self, wavelength: ArrayLike | float) -> Float[np.ndarray, "*shape"]:
        """Return the one-sigma uncertainty of the long-pulse n2 (from the fit).

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            The uncertainty sigma, m^2/W.
        """
        wavelength = np.asarray(wavelength, dtype=float)
        n, a, b, _ = self.parameters
        e2 = _photon_energy(wavelength) ** 2
        jac = np.stack(
            [1 + a * e2 + b * e2**2, n * e2, n * e2**2, np.zeros_like(e2)], axis=-1
        )
        return np.asarray(
            np.sqrt(np.einsum("...i,ij,...j->...", jac, self.covariance, jac))
        )


@dataclass(frozen=True)
class KerrResult:
    """A Kerr fit and its statistics.

    Attributes
    ----------
    config : KerrConfig
        The configuration.
    model : KerrModel
        The fitted model.
    chi2 : float
        The chi-squared, priors included.
    dof : int
        Degrees of freedom.
    dataset_chi2 : dict of str to tuple
        Points and chi-squared per dataset.
    validation : dict of str to list of float
        Per validator, the ratio of each value to the model's prediction for its pulse.
    """

    config: KerrConfig
    model: KerrModel
    chi2: float
    dof: int
    dataset_chi2: Mapping[str, tuple[int, float]]
    validation: Mapping[str, list[float]]


def _gain_residuals(
    gains: list[GainPoints], parameters: tuple[float, float, float, float]
) -> Float[np.ndarray, " n"]:
    """Return the weighted residuals of the measured Raman gains.

    Parameters
    ----------
    gains : list of GainPoints
        The gains.
    parameters : tuple of float
        N (m^2/W), a (eV^-2), b (eV^-4) and f_R.

    Returns
    -------
    numpy.ndarray
        ``(2 omega_s f_R n2(pump) Im H_R / c - g) / sigma`` per gain.
    """
    n, a, b, fraction = parameters
    out = []
    for g in gains:
        e2 = float(_photon_energy(np.array(g.wavelength)) ** 2)
        predicted = g.factor() * fraction * n * (1 + a * e2 + b * e2**2)
        out.append((predicted - g.value) / g.sigma)
    return np.array(out, dtype=float)


def _least_squares(
    residuals: Callable[[np.ndarray], np.ndarray], x0: np.ndarray
) -> "OptimizeResult":
    """Minimise the residuals by Levenberg-Marquardt from a start.

    scipy.optimize is imported here, on the first fit, not when the module loads.

    Parameters
    ----------
    residuals : callable
        The residual vector of the parameters.
    x0 : numpy.ndarray
        The start.

    Returns
    -------
    scipy.optimize.OptimizeResult
        The solution.
    """
    from scipy import optimize  # noqa: PLC0415 - slow to import

    return optimize.least_squares(residuals, x0, method="lm")


def fit(config: KerrConfig) -> KerrResult:
    """Fit n2(lambda), f_R and the calibration scales to the configured data.

    Parameters
    ----------
    config : KerrConfig
        The configuration.

    Returns
    -------
    KerrResult
        The fit.
    """
    observed = [points(o) for o in config.fitted]
    gains = [gain_points(o) for o in config.raman_gains]
    groups: dict[str, float] = {}
    for options in config.fitted:
        groups.setdefault(options.group, options.scale_sigma)
    names = list(groups)
    # Without a prior the start is Stolen et al.'s 0.18 (TC-393).
    f_value = 0.18 if config.raman_fraction is None else config.raman_fraction[0]
    first = np.median(np.concatenate([p.value for p in observed]))
    # The vector: log N, a, (b if fitted), f_R, then the groups' log scales.
    core = 4 if config.quartic else 3
    x0 = np.array(
        [math.log(first), 0.0] + [0.0] * (core - 3) + [f_value] + [0.0] * len(names)
    )

    def unpack(x: np.ndarray) -> tuple[float, float, float, float]:
        """Return N, a, b and f_R from the vector.

        Parameters
        ----------
        x : numpy.ndarray
            The parameters.

        Returns
        -------
        tuple of float
            The parameters N, a, b and f_R.
        """
        b = float(x[2]) if config.quartic else 0.0
        return math.exp(x[0]), float(x[1]), b, float(x[core - 1])

    def predict(x: np.ndarray, p: Points) -> np.ndarray:
        """Return the model's n2 for a dataset's points, before its group's scale.

        Parameters
        ----------
        x : numpy.ndarray
            The parameters.
        p : Points
            The points.

        Returns
        -------
        numpy.ndarray
            The index n2, m^2/W.
        """
        n, a, b, fr = unpack(x)
        e2 = _photon_energy(p.wavelength) ** 2
        return n * (1 + a * e2 + b * e2**2) * ((1 - fr) + fr * p.eta)

    def residuals(x: np.ndarray) -> np.ndarray:
        """Return the weighted residuals and the priors' terms.

        Parameters
        ----------
        x : numpy.ndarray
            The parameters log N, a, (b), f_R, then the groups' log scales.

        Returns
        -------
        numpy.ndarray
            The residuals.
        """
        scale = dict(zip(names, np.exp(x[core:]), strict=True))
        parts = [(scale[p.group] * predict(x, p) - p.value) / p.sigma for p in observed]
        parts.append(_gain_residuals(gains, unpack(x)))
        if config.raman_fraction is not None:
            prior = config.raman_fraction
            parts.append(np.array([(x[core - 1] - prior[0]) / prior[1]]))
        parts.append(np.array([x[core + i] / groups[g] for i, g in enumerate(names)]))
        return np.concatenate(parts)

    solution = _least_squares(residuals, x0)
    r = solution.fun
    chi2 = float(r @ r)
    dof = max(len(r) - len(x0), 1)
    jac = np.asarray(solution.jac, dtype=float)
    covariance = np.linalg.pinv(jac.T @ jac)
    if chi2 / dof > 1:
        covariance *= chi2 / dof
    parameters = unpack(solution.x)
    # The covariance of (N, a, b, f_R), b's row and column zero when it is not fitted.
    rows = [0, 1, 2, 3] if config.quartic else [0, 1, 3]
    full = np.zeros((4, 4))
    full[np.ix_(rows, rows)] = covariance[:core, :core]
    chain = np.diag([parameters[0], 1.0, 1.0, 1.0])
    model = KerrModel(
        material=config.material,
        name=f"kerr/{config.name}",
        citation=f"optmat Kerr fit {config.name} (TC-393)",
        wavelength_range=(config.range_nm[0] * 1e-9, config.range_nm[1] * 1e-9),
        parameters=parameters,
        covariance=chain @ full @ chain,
        scales={
            g: float(math.exp(v)) for g, v in zip(names, solution.x[core:], strict=True)
        },
    )
    dataset_chi2 = {}
    position = 0
    for p in observed:
        part = r[position : position + len(p.value)]
        dataset_chi2[p.key] = (len(p.value), float(part @ part))
        position += len(p.value)
    for g in gains:
        dataset_chi2[g.key] = (1, float(r[position] ** 2))
        position += 1
    validation = {}
    for options in config.validators:
        p = points(options)
        validation[p.key] = [float(v) for v in p.value / predict(solution.x, p)]
    return KerrResult(config, model, chi2, dof, dataset_chi2, validation)


def dumps(result: KerrResult) -> str:
    """Serialise a Kerr fit as JSON text.

    Parameters
    ----------
    result : KerrResult
        The fit.

    Returns
    -------
    str
        Indented JSON.
    """
    model = result.model
    record = {
        "name": result.config.name,
        "material": model.material,
        "config_sha256": result.config.digest,
        "range_nm": list(result.config.range_nm),
        "parameters": {
            "N_m2_per_W": model.parameters[0],
            "a_per_eV2": model.parameters[1],
            "b_per_eV4": model.parameters[2],
            "raman_fraction": model.parameters[3],
        },
        "covariance": model.covariance.tolist(),
        "scales": dict(model.scales),
        "chi2": result.chi2,
        "dof": result.dof,
        "dataset_chi2": {k: list(v) for k, v in result.dataset_chi2.items()},
        "validation": dict(result.validation),
    }
    return json.dumps(record, indent=1) + "\n"


@functools.cache
def kerr_models() -> dict[str, KerrModel]:
    """Return every stored Kerr fit, named ``kerr/<name>``.

    Returns
    -------
    dict of str to KerrModel
        The models.
    """
    out = {}
    for path in sorted(KERR_DIRECTORY.glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        p = record["parameters"]
        low, high = record["range_nm"]
        name = f"kerr/{record['name']}"
        out[name] = KerrModel(
            material=record["material"],
            name=name,
            citation=f"optmat Kerr fit {record['name']} (TC-393)",
            wavelength_range=(low * 1e-9, high * 1e-9),
            parameters=(
                p["N_m2_per_W"],
                p["a_per_eV2"],
                p["b_per_eV4"],
                p["raman_fraction"],
            ),
            covariance=np.asarray(record["covariance"], dtype=float),
            scales=record["scales"],
            status="validated",
        )
    return out
