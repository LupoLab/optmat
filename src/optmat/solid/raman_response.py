# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The Raman response of fused silica (TC-691, TC-692).

A glass's delayed third-order response is a normalised shape :math:`h_R` (unit
integral) and a size: the Raman share :math:`f_R` of the long-pulse nonlinear index
:math:`n_2`. For an envelope, the intensity-dependent index is

.. math::

    \Delta n(t) = n_2 \left[(1 - f_R) I(t) + f_R \int_0^\infty h_R(\tau)
        I(t - \tau)\, d\tau\right],

the form of Stolen et al. (1989), Blow & Wood (1989, their envelope equation), Lin &
Agrawal (2006) and optmat's Kerr fit (:mod:`optmat.solid.kerr`). For a real field, in
the convention of the gases (:mod:`optmat.gas.raman_response`),

.. math::

    P(t) = \varepsilon_0 \chi_e E^3(t) + \varepsilon_0 K E(t)
        \int_0^\infty h_R(\tau) E^2(t - \tau)\, d\tau, \qquad
    \chi_e = \tfrac{4}{3} n_0^2 \varepsilon_0 c (1 - f_R) n_2, \quad
    K = 2 n_0^2 \varepsilon_0 c f_R n_2,

because a kernel slower than the optical cycle follows only the slowly varying half
of :math:`E^2`, while the instantaneous term also has the :math:`3\omega` part: the
Raman term of the field carries 2/3 of the envelope's share (Blow & Wood's Eqs. (2)
and (8), TC-692). :math:`K h_R` is the silica counterpart of a gas's :math:`N h`.

The steady-state Stokes gain, pump and Stokes polarised alike, is
:math:`g(\Omega) = 2 \omega_s n_2 f_R \operatorname{Im} H_R(\Omega) / c`, with
:math:`H_R(\Omega) = \int_0^\infty h_R(\tau) e^{i\Omega\tau} d\tau` (Lin & Agrawal's
:math:`2\gamma f_R \operatorname{Im}\tilde R`).

**The models** (``data/solid/raman/SiO2.toml``):

- ``hollenbeck_2002``, the default: Hollenbeck & Cantrell's 13 modes,
  :math:`\sum_i A_i e^{-\gamma_i t} e^{-\Gamma_i^2 t^2 / 4} \sin(\omega_i t)`, with the
  Raman fraction of optmat's Kerr fit, which uses this shape. Its normalisation is
  analytic: :math:`\int_0^\infty e^{-\gamma t - \Gamma^2 t^2/4} \sin(\omega t)\, dt
  = (\sqrt\pi / \Gamma) \operatorname{Im} w[(\omega + i\gamma)/\Gamma]`, with
  :math:`w` the Faddeeva function. It is copolarised only.
- ``lin_2006``: Lin & Agrawal's isotropic part :math:`R_a = f_a h_a` and anisotropic
  part :math:`R_b = f_b h_b + f_c h_a`, the functions :math:`a` and :math:`b` of
  Hellwarth's form (:mod:`optmat.gas.raman_response`), with their :math:`f_R = 0.245`.
- ``blow_1989``: Blow & Wood's single oscillator with its share 0.222.

The literature forms' :math:`f_R` are theirs, not the Kerr fit's: with them the
electronic part of the fitted :math:`n_2` and the Raman part do not add up to the
fitted :math:`n_2` (``paired`` is False).

References
----------
- K. J. Blow and D. Wood, "Theoretical description of transient stimulated Raman
  scattering in optical fibers", IEEE J. Quantum Electron. 25, 2665 (1989).
- R. H. Stolen, J. P. Gordon, W. J. Tomlinson and H. A. Haus, "Raman response
  function of silica-core fibers", J. Opt. Soc. Am. B 6, 1159 (1989).
- D. Hollenbeck and C. D. Cantrell, "Multiple-vibrational-mode model for fiber-optic
  Raman gain spectrum and response function", J. Opt. Soc. Am. B 19, 2886 (2002).
- Q. Lin and G. P. Agrawal, "Raman response function for silica fibers", Opt. Lett.
  31, 3086 (2006).
"""

import functools
import math
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Any, ClassVar, Final, Literal, get_args

import numpy as np
from jaxtyping import Complex, Float
from numpy.typing import ArrayLike
from scipy import special

from optmat import _catalogue
from optmat._catalogue import records as sources_of
from optmat.constants import EPS0, C
from optmat.exceptions import DataError
from optmat.identity import Identified, source_ids
from optmat.solid import SOLID_DATA_DIRECTORY
from optmat.solid.models.base import Model as SolidModel
from optmat.solid.models.registry import default_model

if TYPE_CHECKING:
    from optmat.solid.kerr import KerrModel

__all__ = [
    "EVALUATION_VERSION",
    "RAMAN_DIRECTORY",
    "Component",
    "Form",
    "SilicaRamanResponse",
    "silica_raman_response",
    "silica_raman_responses",
]

EVALUATION_VERSION: Final = 1
"""The version of this module's evaluation code, in the identity of every model that
evaluates through it; raised whenever the code changes its output (OEP 1, section 4.8).
"""

RAMAN_DIRECTORY: Final = SOLID_DATA_DIRECTORY / "raman"
"""The solids' Raman response models."""

Form = Literal["hollenbeck", "lin_agrawal", "single_oscillator"]
"""The functional form of a model (the module docstring)."""

Component = Literal["parallel", "a", "b"]
"""A part of the normalised response: ``parallel`` (:math:`h_R = a + b`), or
Hellwarth's ``a`` (isotropic) and ``b`` (anisotropic) for a model that has them."""

_CM_S: Final = C * 100.0
"""The speed of light in cm/s: wavenumbers in cm^-1 to frequencies."""

_TIME_UNITS: Final = {"fs": 1e-15, "ps": 1e-12, "s": 1.0}
"""Units of the formula coefficients that are times."""


@dataclass(frozen=True)
class SilicaRamanResponse(Identified):
    r"""A Raman response model of fused silica.

    Attributes
    ----------
    name : str
        The model (``hollenbeck_2002``, ``lin_2006``, ``blow_1989``).
    form : str
        Its functional form (:data:`Form`).
    raman_fraction : float
        :math:`f_R`, the Raman share of the long-pulse :math:`n_2`.
    paired : bool
        True where :math:`f_R` is the one optmat's Kerr fit found with this shape.
    parameters : dict of str to numpy.ndarray
        The shape's constants in SI: for ``hollenbeck`` the amplitudes ``A``, angular
        frequencies ``omega`` (rad/s), Gaussian ``Gamma`` and Lorentzian ``gamma``
        widths (1/s); for ``lin_agrawal`` ``tau1``, ``tau2``, ``tau_b`` (s) and
        ``f_a``, ``f_b``, ``f_c``; for ``single_oscillator`` ``tau1``, ``tau2`` (s).
    description : str
        What the model is.
    sources : dict of str to str
        The datasets and formulas behind it.
    reasons : dict of str to str
        Why.
    status : str
        ``validated`` for the default, ``comparison`` for the other shapes.

    Its overrides are the shape's constants (``parameters``, any of them, merged
    into the model's) and ``raman_fraction``. Any override sets :attr:`paired` to
    False: the fraction was found with the curated shape.
    """

    IDENTITY_KIND: ClassVar[str] = "raman_response"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ("optmat.solid.raman_response",)
    INPUTS: ClassVar[tuple[str, ...]] = ("parameters", "raman_fraction")
    DERIVED: ClassVar[tuple[str, ...]] = ("paired",)

    name: str
    form: str
    raman_fraction: float
    paired: bool
    parameters: Mapping[str, Float[np.ndarray, " m"]]
    description: str = ""
    sources: Mapping[str, str] = field(default_factory=dict)
    reasons: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Check the form.

        Raises
        ------
        DataError
            If the form is not one of :data:`Form`.
        """
        super().__post_init__()
        if self.form not in get_args(Form):
            msg = f"unknown Raman form {self.form!r}; one of {get_args(Form)}"
            raise DataError(msg)

    def _override(self, changes: Mapping[str, Any]) -> "SilicaRamanResponse":
        """Return the response with its shape or fraction changed, unpaired.

        Parameters
        ----------
        changes : Mapping
            ``parameters`` (a subset of the shape's) and ``raman_fraction``.

        Returns
        -------
        SilicaRamanResponse
            The rebuilt response, ``paired`` False.

        Raises
        ------
        DataError
            If ``parameters`` names a constant the shape does not have.
        """
        changed = dict(changes)
        if "parameters" in changed:
            unknown = sorted(set(changed["parameters"]) - set(self.parameters))
            if unknown:
                msg = f"{self.name} has no shape constant {', '.join(unknown)}"
                raise DataError(msg)
            changed["parameters"] = dict(self.parameters) | {
                k: np.atleast_1d(np.asarray(v, dtype=float))
                for k, v in changed["parameters"].items()
            }
        return replace(self, **changed, paired=False)

    @property
    def has_tensor(self) -> bool:
        """Whether the model has an isotropic and an anisotropic part.

        Returns
        -------
        bool
            True for Lin & Agrawal's form.
        """
        return self.form == "lin_agrawal"

    def _weights(self, component: Component) -> tuple[float, float]:
        """Return the weights of Lin & Agrawal's h_a and h_b in a component.

        Parameters
        ----------
        component : {"parallel", "a", "b"}
            The part.

        Returns
        -------
        tuple of float
            The weights of h_a and h_b.
        """
        p = self.parameters
        f_a, f_b, f_c = (float(p[k][0]) for k in ("f_a", "f_b", "f_c"))
        return {
            "parallel": (f_a + f_c, f_b),
            "a": (f_a, 0.0),
            "b": (f_c, f_b),
        }[component]

    def _check_component(self, component: Component) -> None:
        """Refuse a tensor part of a model that has none.

        Parameters
        ----------
        component : {"parallel", "a", "b"}
            The part.

        Raises
        ------
        DataError
            If the part is not ``parallel`` and the model has no tensor form.
        """
        if component not in get_args(Component):
            msg = f"unknown component {component!r}; one of {get_args(Component)}"
            raise DataError(msg)
        if component != "parallel" and not self.has_tensor:
            msg = f"{self.name} is copolarised only: it has no parts a and b"
            raise DataError(msg)

    def response(
        self, time: ArrayLike, component: Component = "parallel"
    ) -> Float[np.ndarray, " *shape"]:
        r"""Return the normalised response :math:`h_R(\tau)` (or one of its parts).

        Parameters
        ----------
        time : array_like
            Delay after the impulse, s (zero before it).
        component : {"parallel", "a", "b"}, optional
            The part (:data:`Component`).

        Returns
        -------
        numpy.ndarray
            :math:`h_R`, 1/s; ``parallel`` integrates to one.

        Raises
        ------
        DataError
            If a tensor part is asked of a copolarised model.

        Examples
        --------
        >>> response = silica_raman_response()
        >>> t = np.arange(0.0, 3e-12, 0.05e-15)
        >>> f"{np.trapezoid(response.response(t), t):.4f}"
        '1.0000'
        """
        self._check_component(component)
        tau = np.asarray(time, dtype=float)
        t = np.maximum(tau, 0.0)
        p = self.parameters
        if self.form == "hollenbeck":
            terms = (
                p["A"]
                * np.exp(-p["gamma"] * t[..., None])
                * np.exp(-(p["Gamma"] ** 2) * t[..., None] ** 2 / 4)
                * np.sin(p["omega"] * t[..., None])
            )
            h = np.sum(terms, axis=-1) / self._hollenbeck_integral()
        elif self.form == "single_oscillator":
            h = _oscillator(t, float(p["tau1"][0]), float(p["tau2"][0]))
        else:
            weight_a, weight_b = self._weights(component)
            tau_b = float(p["tau_b"][0])
            h = weight_a * _oscillator(t, float(p["tau1"][0]), float(p["tau2"][0]))
            h = h + weight_b * (2 * tau_b - t) / tau_b**2 * np.exp(-t / tau_b)
        return np.where(tau > 0, h, 0.0)

    def _hollenbeck_integral(self) -> float:
        """Return the integral of Eq. (9) as printed, from the Faddeeva function.

        Returns
        -------
        float
            The integral, s.
        """
        p = self.parameters
        z = (p["omega"] + 1j * p["gamma"]) / p["Gamma"]
        return float(
            np.sum(p["A"] * math.sqrt(math.pi) / p["Gamma"] * np.imag(special.wofz(z)))
        )

    def spectrum(
        self, omega: ArrayLike, component: Component = "parallel"
    ) -> Complex[np.ndarray, " *shape"]:
        r"""Return :math:`H_R(\Omega) = \int_0^\infty h_R(\tau) e^{i\Omega\tau} d\tau`.

        Parameters
        ----------
        omega : array_like
            Angular frequency, rad/s.
        component : {"parallel", "a", "b"}, optional
            The part (:data:`Component`).

        Returns
        -------
        numpy.ndarray
            :math:`H_R`, dimensionless, complex; ``parallel`` is 1 at zero and its
            imaginary part is positive at a Stokes shift (gain).

        Raises
        ------
        DataError
            If a tensor part is asked of a copolarised model.

        Examples
        --------
        >>> response = silica_raman_response()
        >>> f"{response.spectrum(0.0).real:.6f}"
        '1.000000'
        """
        self._check_component(component)
        w = np.asarray(omega, dtype=float)
        p = self.parameters
        if self.form == "hollenbeck":
            x = w[..., None]
            gamma, width, centre = p["gamma"], p["Gamma"], p["omega"]
            plus = special.wofz((centre + x + 1j * gamma) / width)
            minus = special.wofz((x - centre + 1j * gamma) / width)
            terms = p["A"] * math.sqrt(math.pi) / width * (plus - minus) / 2j
            return np.asarray(np.sum(terms, axis=-1) / self._hollenbeck_integral())
        if self.form == "single_oscillator":
            return _oscillator_spectrum(w, float(p["tau1"][0]), float(p["tau2"][0]))
        weight_a, weight_b = self._weights(component)
        tau_b = float(p["tau_b"][0])
        s = 1 / tau_b - 1j * w
        boson = (2 * s * tau_b - 1) / (tau_b**2 * s**2)
        oscillator = _oscillator_spectrum(w, float(p["tau1"][0]), float(p["tau2"][0]))
        return np.asarray(weight_a * oscillator + weight_b * boson)

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the model's name, material and source records.

        Returns
        -------
        tuple
            ``("silica_raman/<name>", "SiO2", sources)``.
        """
        return f"silica_raman/{self.name}", "SiO2", source_ids(*self.sources.values())

    def _identity_extras(self) -> Mapping[str, Any]:
        """Return the models :meth:`raman_coefficient` reads, by their identities.

        Returns
        -------
        Mapping
            ``index_model`` and ``kerr_model``.
        """
        return {"index_model": self.index_model, "kerr_model": self.kerr_model}

    @functools.cached_property
    def index_model(self) -> SolidModel:
        """The linear model of silica the kernel's size reads, bound once.

        Returns
        -------
        optmat.solid.models.base.Model
            The recommended model of SiO2.
        """
        return default_model("SiO2")

    @functools.cached_property
    def kerr_model(self) -> "KerrModel":
        """The Kerr fit whose n2 the kernel's size reads, bound once.

        Returns
        -------
        optmat.solid.kerr.KerrModel
            ``kerr/SiO2``.
        """
        from optmat.solid.kerr import kerr_models  # noqa: PLC0415 - kerr imports this

        return kerr_models()["kerr/SiO2"]

    def _index_and_n2(self, wavelength: float) -> tuple[float, float]:
        """Return the linear index and the long-pulse n2 of silica at a wavelength.

        Parameters
        ----------
        wavelength : float
            Vacuum wavelength, m.

        Returns
        -------
        tuple of float
            ``n0`` (the recommended linear model's) and ``n2`` (m^2/W, the Kerr fit's).
        """
        n0 = float(np.real(self.index_model.refractive_index(wavelength)))
        n2 = float(self.kerr_model.n2(wavelength))
        return n0, n2

    def raman_coefficient(self, wavelength: float) -> float:
        r"""Return :math:`K = 2 n_0^2 \varepsilon_0 c f_R n_2`, the size of the kernel.

        In the real-field form of the module docstring, the delayed polarisation is
        :math:`\varepsilon_0 K E \int h_R E^2`. It includes :math:`n_0^2`: a code in
        the :math:`n_0 = 1` convention of the nonlinear Schrödinger equation takes
        :math:`n_2`, :math:`f_R` and :math:`h_R` from a silica bundle instead
        (:func:`optmat.solid.bundle.silica_bundle`).

        Parameters
        ----------
        wavelength : float
            Vacuum wavelength of the field, m (for :math:`n_0` and :math:`n_2`).

        Returns
        -------
        float
            :math:`K`, m^2/V^2.

        Warns
        -----
        ModelRangeWarning
            Outside the Kerr fit's wavelengths.
        """
        n0, n2 = self._index_and_n2(wavelength)
        return 2 * n0**2 * EPS0 * C * self.raman_fraction * n2

    def electronic_susceptibility(self, wavelength: float) -> float:
        r"""Return :math:`\chi_e = \tfrac43 n_0^2 \varepsilon_0 c (1 - f_R) n_2`.

        Parameters
        ----------
        wavelength : float
            Vacuum wavelength of the field, m.

        Returns
        -------
        float
            The instantaneous third-order susceptibility, m^2/V^2.

        Warns
        -----
        ModelRangeWarning
            Outside the Kerr fit's wavelengths.
        """
        n0, n2 = self._index_and_n2(wavelength)
        return 4 / 3 * n0**2 * EPS0 * C * (1 - self.raman_fraction) * n2

    def gain(
        self,
        pump_wavelength: float,
        shift: ArrayLike,
        *,
        n2: float | None = None,
    ) -> Float[np.ndarray, " *shape"]:
        r"""Return the copolarised steady-state Stokes gain.

        :math:`g = 2 \omega_s n_2 f_R \operatorname{Im} H_R(\Omega) / c`.

        Parameters
        ----------
        pump_wavelength : float
            Vacuum wavelength of the pump, m.
        shift : array_like
            The Stokes shift :math:`\Omega`, rad/s.
        n2 : float, optional
            The long-pulse :math:`n_2`, m^2/W; by default optmat's Kerr fit at the
            pump wavelength.

        Returns
        -------
        numpy.ndarray
            The gain, m/W.

        Warns
        -----
        ModelRangeWarning
            Outside the Kerr fit's wavelengths (without ``n2``).
        """
        index = self._index_and_n2(pump_wavelength)[1] if n2 is None else n2
        big_omega = np.asarray(shift, dtype=float)
        omega_s = 2 * math.pi * C / pump_wavelength - big_omega
        return np.asarray(
            2
            * omega_s
            * index
            * self.raman_fraction
            * np.imag(self.spectrum(big_omega))
            / C
        )


def _oscillator(
    t: Float[np.ndarray, " *shape"], tau1: float, tau2: float
) -> Float[np.ndarray, " *shape"]:
    """Return Blow & Wood's normalised damped oscillator.

    Parameters
    ----------
    t : numpy.ndarray
        Delays, s, >= 0.
    tau1, tau2 : float
        The period over 2 pi and the decay time, s.

    Returns
    -------
    numpy.ndarray
        ``(tau1^2 + tau2^2) / (tau1 tau2^2) exp(-t/tau2) sin(t/tau1)``, 1/s.
    """
    return np.asarray(
        (tau1**2 + tau2**2) / (tau1 * tau2**2) * np.exp(-t / tau2) * np.sin(t / tau1)
    )


def _oscillator_spectrum(
    w: Float[np.ndarray, " *shape"], tau1: float, tau2: float
) -> Complex[np.ndarray, " *shape"]:
    """Return the transform of :func:`_oscillator`.

    Parameters
    ----------
    w : numpy.ndarray
        Angular frequency, rad/s.
    tau1, tau2 : float
        As :func:`_oscillator`.

    Returns
    -------
    numpy.ndarray
        ``(Omega^2 + Gamma^2) / (Omega^2 + Gamma^2 - w^2 - 2 i Gamma w)`` with
        ``Omega = 1/tau1``, ``Gamma = 1/tau2``: 1 at zero.
    """
    big_omega, rate = 1 / tau1, 1 / tau2
    norm = big_omega**2 + rate**2
    return np.asarray(norm / (norm - w**2 - 2j * rate * w))


def _seconds(value: str) -> float:
    """Return a time written with its unit (``"12.2 fs"``) in s.

    Parameters
    ----------
    value : str
        A number and a unit.

    Returns
    -------
    float
        The time, s.

    Raises
    ------
    DataError
        If the unit is not a time.
    """
    number, _, unit = value.partition(" ")
    if unit not in _TIME_UNITS:
        msg = f"{value!r} is not a time in {sorted(_TIME_UNITS)}"
        raise DataError(msg)
    return float(number) * _TIME_UNITS[unit]


def _formula(reference: str) -> dict[str, str]:
    """Return a formula record's coefficients.

    Parameters
    ----------
    reference : str
        ``source/formula``.

    Returns
    -------
    dict of str to str
        The coefficients as printed.

    Raises
    ------
    DataError
        If the formula does not exist.
    """
    source_id, _, formula_id = reference.partition("/")
    matches = [f for f in sources_of()[source_id].formulas if f["id"] == formula_id]
    if len(matches) != 1:
        msg = f"SiO2.toml names an unknown formula {reference!r}"
        raise DataError(msg)
    return _strings(matches[0]["coefficients"])


def _strings(value: object) -> dict[str, str]:
    """Return a TOML table as strings (nothing for a value that is not a table).

    Parameters
    ----------
    value : object
        The table.

    Returns
    -------
    dict of str to str
        Its keys and values.
    """
    if not isinstance(value, Mapping):
        return {}
    return {str(k): str(v) for k, v in value.items()}  # pyright: ignore[reportUnknownVariableType]


def _hollenbeck_parameters(reference: str) -> dict[str, np.ndarray]:
    """Return Hollenbeck & Cantrell's modes in SI.

    Table 1's footnote: ``omega = 2 pi c`` x position, ``Gamma = pi c`` x the
    "Gaussian FWHM" column, ``gamma = pi c`` x the "Lorentzian FWHM" column; the
    "Peak Intensity" column is Eq. (9)'s ``A_i'/omega_i`` (TC-382, TC-691).

    Parameters
    ----------
    reference : str
        The ``raman_mode`` dataset.

    Returns
    -------
    dict of str to numpy.ndarray
        ``A``, ``omega``, ``Gamma``, ``gamma``.
    """
    key = reference
    dataset = _catalogue.by_key("SiO2", quantity="raman_mode")[key]
    return {
        "A": dataset.scaled_values(),
        "omega": 2 * math.pi * _CM_S * dataset.column_values("shift"),
        "Gamma": math.pi * _CM_S * dataset.column_values("gaussian_width"),
        "gamma": math.pi * _CM_S * dataset.column_values("lorentzian_width"),
    }


def _fraction(entry: Mapping[str, object], form: str) -> float:
    """Return a model's Raman fraction from its entry.

    Parameters
    ----------
    entry : mapping
        The model's table in ``SiO2.toml``.
    form : str
        Its form.

    Returns
    -------
    float
        :math:`f_R`, the share of the long-pulse n2.
    """
    reference = str(entry["raman_fraction"])
    if reference.startswith("kerr/"):
        from optmat.solid.kerr import kerr_models  # noqa: PLC0415 - kerr imports this

        return kerr_models()[reference].raman_fraction
    coefficients = _formula(reference)
    if form == "single_oscillator":
        # alpha acts on E^2; the envelope's Raman term carries 2/3 of 1 - alpha.
        alpha = float(coefficients["alpha"])
        numerator, _, denominator = coefficients["envelope_factor"].partition("/")
        envelope = float(numerator) / float(denominator)
        return envelope * (1 - alpha) / (alpha + envelope * (1 - alpha))
    return float(coefficients["f_R"])


def _build(name: str, entry: Mapping[str, object]) -> SilicaRamanResponse:
    """Build a model from its entry in ``SiO2.toml``.

    Parameters
    ----------
    name : str
        The model.
    entry : mapping
        Its table.

    Returns
    -------
    SilicaRamanResponse
        The model.
    """
    form = str(entry["form"])
    sources = {
        "formula": str(entry["formula"]),
        "raman_fraction": str(entry["raman_fraction"]),
    }
    if form == "hollenbeck":
        parameters = _hollenbeck_parameters(str(entry["modes"]))
        sources["modes"] = str(entry["modes"])
    else:
        c = _formula(str(entry["formula"]))
        parameters = {
            "tau1": np.array([_seconds(c["tau_1"])]),
            "tau2": np.array([_seconds(c["tau_2"])]),
        }
        if form == "lin_agrawal":
            parameters["tau_b"] = np.array([_seconds(c["tau_b"])])
            for key in ("f_a", "f_b", "f_c"):
                parameters[key] = np.array([float(c[key])])
    reasons = entry.get("reasons", {})
    return SilicaRamanResponse(
        name=name,
        form=form,
        raman_fraction=_fraction(entry, form),
        paired=bool(entry["paired"]),
        parameters=parameters,
        description=str(entry.get("description", "")),
        sources=sources,
        reasons=_strings(reasons),
    )


@functools.cache
def silica_raman_responses() -> dict[str, SilicaRamanResponse]:
    """Return every Raman response model of fused silica.

    Returns
    -------
    dict of str to SilicaRamanResponse
        By name.
    """
    with (RAMAN_DIRECTORY / "SiO2.toml").open("rb") as handle:
        table = tomllib.load(handle)
    default = str(table["default"])
    return {
        str(k): replace(
            _build(str(k), v),
            status="validated" if k == default else "comparison",
            status_reason=""
            if k == default
            else f"a shape kept for comparison; the default is {default}",
        )
        for k, v in dict(table["models"]).items()
    }


def silica_raman_response(model: str = "default") -> SilicaRamanResponse:
    """Return a Raman response model of fused silica.

    Parameters
    ----------
    model : str, optional
        ``default`` (Hollenbeck & Cantrell with the Kerr fit's f_R), ``lin_2006`` or
        ``blow_1989``.

    Returns
    -------
    SilicaRamanResponse
        The model.

    Raises
    ------
    DataError
        If the model is unknown.

    Examples
    --------
    >>> silica_raman_response().name
    'hollenbeck_2002'
    >>> f"{silica_raman_response('blow_1989').raman_fraction:.3f}"
    '0.222'
    """
    models = silica_raman_responses()
    if model == "default":
        with (RAMAN_DIRECTORY / "SiO2.toml").open("rb") as handle:
            model = str(tomllib.load(handle)["default"])
    if model not in models:
        msg = f"no silica Raman model {model!r}; one of {sorted(models)}"
        raise DataError(msg)
    return models[model]
