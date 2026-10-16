# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Two-photon ionisation of the noble gases: the continuum model (OEP 6, Phase 1).

The generalised two-photon ionisation cross section :math:`\sigma^{(2)}` of a gas's
ground state is defined by the rate per atom

.. math:: W = \sigma^{(2)} \Phi^2, \qquad \Phi = I / \hbar\omega,

with :math:`\Phi` the photon flux; in SI, :math:`\sigma^{(2)}` is in m^4 s (1 cm^4 s
= 1e-8 m^4 s). A quasi-monochromatic wave in a gas of number density :math:`N` then
loses intensity as :math:`dI/dz = -\beta I^2` with :math:`\beta = 2 N \sigma^{(2)} /
\hbar\omega`, two photons per event.

Each model is a curated reference curve that the author approved under data rule 8
(TC-943), interpolated linearly in :math:`\log\sigma^{(2)}` between the curve's own
drawn vertices and never fitted or extrapolated:

- helium: Mihelič and Horvat (2021), the 1S and 1D partial cross sections summed;
- neon and argon: McKenna and van der Hart (2004), R-matrix Floquet;
- xenon: Lagutin et al. (2017), linear and circular polarisation.

Krypton has none: no vector-drawn or tabulated cross section is curated (OEP 6, Q2),
and asking for it raises :class:`~optmat.exceptions.UnavailableModelError`.

A model covers the window from the lowest two-photon threshold :math:`I_p/2` to the
first one-photon resonance :math:`E_1`, and states its regions
(:attr:`~TwoPhotonIonisationModel.regions`):

- ``below_threshold``: the direct continuum cross section is zero. Bound two-photon
  excitation and (2+1) ionisation there are not modelled (OEP 6, Phase 2).
- ``supported``: the reference's values.
- ``unsupported_detail``: the reference's values, with
  :class:`~optmat.exceptions.ModelRangeWarning`. The source does not qualify them: an
  LS-coupled calculation between the ion's two fine-structure thresholds, or the wing
  of the first one-photon resonance, which the sources place differently.
- ``unavailable``: no value. The reference does not reach it (between xenon's two
  thresholds, for instance), or it is at or above :math:`E_1`.

Evaluation never raises on a frequency grid. It returns NaN where the model is
unavailable, and :meth:`~TwoPhotonIonisationModel.available` gives the mask (the
author's decision, TC-943). A propagation code decides what to do there.

The model supplies the degenerate cross section. It does not supply a time-resolved
loss for a broadband pulse, nor the mixed pairs of two bands. The pair-spectral yield
of OEP 6, section 4, is an integral identity in a stated approximation, and its time
dependence is a documented gap (TC-903). Nor does it double as a strong-field rate:
applying it and :func:`~optmat.gas.ionisation_rates.ionisation_rate_model` to the same
spectral component counts that component twice.

Examples
--------
>>> import math
>>> from optmat.constants import C
>>> argon = two_photon_ionisation_model("Ar")
>>> [round(2 * math.pi * C / t.omega * 1e9, 2) for t in argon.thresholds]
[157.34, 155.59]
>>> omega = 2 * math.pi * C / 120e-9  # 120 nm
>>> f"{float(argon.cross_section(omega)):.3e}"  # m^4 s
'1.146e-58'
>>> bool(argon.available(2 * math.pi * C / 100e-9))  # beyond E1 (106.7 nm)
False
"""

import functools
import itertools
import json
import math
import tomllib
import warnings
from collections.abc import Mapping
from dataclasses import dataclass, field, fields, is_dataclass
from typing import Any, ClassVar, Final

import numpy as np
from jaxtyping import Float
from numpy.typing import ArrayLike

from optmat.constants import HBAR
from optmat.exceptions import (
    ArgumentError,
    DataError,
    ModelRangeWarning,
    UnavailableModelError,
)
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.gases import Gas, as_gas
from optmat.identity import Identified, parse_status, source_ids

__all__ = [
    "DERIVATION_VERSION",
    "EVALUATION_VERSION",
    "POLARIZATIONS",
    "REGION_LABELS",
    "TwoPhotonIonisationModel",
    "TwoPhotonReference",
    "TwoPhotonRegion",
    "TwoPhotonSpread",
    "TwoPhotonThreshold",
    "TwoPhotonUncertainty",
    "two_photon_ionisation_model",
    "two_photon_ionisation_models",
]

EVALUATION_VERSION: Final = 1
"""The version of the evaluation: interpolation, regions and conversions."""
DERIVATION_VERSION: Final = 1
"""The version of the stored models' derivation, raised whenever it writes otherwise."""
POLARIZATIONS: Final = ("linear", "circular")
REGION_LABELS: Final = (
    "below_threshold",
    "supported",
    "unsupported_detail",
    "unavailable",
)
"""The kinds of region a model's frequency axis is divided into."""

_SELECTION_FILE: Final = GAS_DATA_DIRECTORY / "two_photon.toml"
"""The approved references and their choices (TC-943)."""
_DIRECTORY: Final = GAS_DATA_DIRECTORY / "models" / "two_photon"
_DEFINED: Final = frozenset({"below_threshold", "supported", "unsupported_detail"})
_STANDING: Final = frozenset({"status", "status_reason", "modified", "base"})


@dataclass(frozen=True)
class TwoPhotonThreshold:
    """A threshold or resonance as a photon's angular frequency, with its source.

    Attributes
    ----------
    name : str
        What it is, e.g. ``"Ar+ 2P3/2"`` or ``"E1"``.
    omega : float
        The photon's angular frequency, rad/s: the level's energy over the number of
        photons that reach it, over hbar.
    photons : int
        2 for an ionisation threshold, 1 for :math:`E_1`.
    energy : float
        The energy above the ground state, J.
    datasets : tuple of str
        The curated datasets it is computed from.
    """

    name: str
    omega: float
    photons: int
    energy: float
    datasets: tuple[str, ...]


@dataclass(frozen=True)
class TwoPhotonRegion:
    """A half-open band of photon angular frequency, ``[lower, upper)``, and its label.

    Attributes
    ----------
    label : str
        One of :data:`REGION_LABELS`.
    lower : float
        The lower edge, rad/s.
    upper : float
        The upper edge, rad/s; inf for the last region.
    reason : str
        Why the band is so labelled.
    """

    label: str
    lower: float
    upper: float
    reason: str


@dataclass(frozen=True)
class TwoPhotonReference:
    """The approved reference curve of a model, and why it was chosen.

    Attributes
    ----------
    datasets : tuple of str
        Its datasets, summed where they are partial cross sections.
    gauge : str
        The gauge, as the record states it.
    method : str
        The calculation.
    intensity_w_cm2 : float or None
        The intensity of a non-perturbative calculation, W/cm^2; None for a
        lowest-order one.
    reason : str
        Why it is the reference (the author's approval, TC-943).
    registers : tuple of str
        The register entries.
    """

    datasets: tuple[str, ...]
    gauge: str
    method: str
    intensity_w_cm2: float | None
    reason: str
    registers: tuple[str, ...]


@dataclass(frozen=True)
class TwoPhotonSpread:
    """Another calculation against the reference over their supported overlap.

    Attributes
    ----------
    label : str
        The calculation.
    datasets : tuple of str
        Its datasets, summed where they are partial cross sections.
    form : str
        ``"vector"`` or ``"raster"``: how its values were read.
    role : str
        ``"independent"``, another calculation, which enters the envelope; or
        ``"gauge"``, the reference's other gauge, a diagnostic only.
    lower : float
        The lowest photon angular frequency compared, rad/s.
    upper : float
        The highest, rad/s.
    points : int
        The number of reference vertices compared.
    lower_ratio : float
        The smallest ratio of it to the reference.
    upper_ratio : float
        The largest.
    """

    label: str
    datasets: tuple[str, ...]
    form: str
    role: str
    lower: float
    upper: float
    points: int
    lower_ratio: float
    upper_ratio: float


@dataclass(frozen=True)
class TwoPhotonUncertainty:
    """What is known of a model's uncertainty, kept as three separate kinds (Q6).

    Attributes
    ----------
    source : float or None
        The source's own relative accuracy, None where it states none (never zero).
    source_level : str
        Its level, as stated, or ``"not stated"``.
    source_provenance : str
        Where it is stated.
    digitising : float
        The relative digitising uncertainty at a drawn vertex.
    digitising_provenance : str
        Its record and register entry.
    spreads : tuple of TwoPhotonSpread
        The other calculations against the reference, over the supported region.
    """

    source: float | None
    source_level: str
    source_provenance: str
    digitising: float
    digitising_provenance: str
    spreads: tuple[TwoPhotonSpread, ...]

    @property
    def envelope(self) -> tuple[float, float] | None:
        """The smallest and largest ratio of the independent calculations.

        It is an empirical envelope of other calculations, not a confidence interval
        or a bound on the true cross section.

        Returns
        -------
        tuple of float or None
            ``(lower, upper)`` ratios to the reference; None with no independent
            calculation.
        """
        independent = [s for s in self.spreads if s.role == "independent"]
        if not independent:
            return None
        return (
            min(s.lower_ratio for s in independent),
            max(s.upper_ratio for s in independent),
        )


@dataclass(frozen=True, eq=False)
class TwoPhotonIonisationModel(Identified):
    r"""The two-photon ionisation cross section of a noble gas's ground state.

    Built by :func:`two_photon_ionisation_model` from its stored derivation. SI
    throughout: photon angular frequency in rad/s, cross section in m^4 s.

    Attributes
    ----------
    gas : Gas
        The species.
    polarization : str
        ``"linear"`` or ``"circular"``.
    derivation_version : int
        :data:`DERIVATION_VERSION` when it was derived.
    omega : numpy.ndarray
        The reference's vertices, photon angular frequency, rad/s, increasing.
    sigma : numpy.ndarray
        The reference's cross section at them, m^4 s.
    thresholds : tuple of TwoPhotonThreshold
        The lowest ionic channel's two-photon threshold first; the ion's
        :math:`{}^2P_{1/2}` channel's second, for Ne to Xe.
    e1 : TwoPhotonThreshold
        The first one-photon resonance, where the window ends.
    regions : tuple of TwoPhotonRegion
        Bands covering every positive frequency, in order.
    reference : TwoPhotonReference
        The approved reference.
    uncertainty : TwoPhotonUncertainty
        Its uncertainty, by kind.
    input_digests : Mapping of str to str
        SHA-256 digests of what it was derived from: the selection and each dataset.
    """

    gas: Gas
    polarization: str
    derivation_version: int
    omega: Float[np.ndarray, " n"]
    sigma: Float[np.ndarray, " n"]
    thresholds: tuple[TwoPhotonThreshold, ...]
    e1: TwoPhotonThreshold
    regions: tuple[TwoPhotonRegion, ...]
    reference: TwoPhotonReference
    uncertainty: TwoPhotonUncertainty
    input_digests: Mapping[str, str] = field(default_factory=dict, repr=False)

    IDENTITY_KIND: ClassVar[str] = "two_photon_ionisation"
    EVALUATION_VERSION: ClassVar[int] = EVALUATION_VERSION

    def __post_init__(self) -> None:
        """Take ownership, and check the grid and the regions.

        Raises
        ------
        DataError
            If the grid is not increasing and positive, or the regions do not tile the
            positive axis with known labels.
        """
        super().__post_init__()
        where = f"two-photon model {self.gas.value} ({self.polarization})"
        if (
            self.omega.shape != self.sigma.shape
            or len(self.omega) < 2
            or np.any(np.diff(self.omega) <= 0)
            or np.any(self.sigma <= 0)
        ):
            msg = f"{where}: its grid must be increasing with positive cross sections"
            raise DataError(msg)
        edges = [(r.lower, r.upper) for r in self.regions]
        tiled = bool(edges) and edges[0][0] == 0 and math.isinf(edges[-1][1])
        tiled = tiled and all(a[1] == b[0] for a, b in itertools.pairwise(edges))
        if not tiled or any(r.label not in REGION_LABELS for r in self.regions):
            msg = f"{where}: its regions must tile [0, inf) with known labels"
            raise DataError(msg)

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the model's name, gas and source records.

        Returns
        -------
        tuple
            ``("two_photon/<gas>/<polarization>", gas, sources)``.
        """
        keys = [*self.reference.datasets, *self.e1.datasets]
        keys += [k for t in self.thresholds for k in t.datasets]
        keys += [k for s in self.uncertainty.spreads for k in s.datasets]
        return (
            f"two_photon/{self.gas.value}/{self.polarization}",
            self.gas.value,
            source_ids(*(f"{key}/" for key in keys)),
        )

    @property
    def digest(self) -> str:
        """The model's digest, ``identity.digest``.

        Returns
        -------
        str
            64 hexadecimal digits.
        """
        return self.identity.digest

    def _region_indices(self, omega: np.ndarray) -> np.ndarray:
        """Return the index of each frequency's region.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon angular frequencies, rad/s, positive.

        Returns
        -------
        numpy.ndarray
            Indices into :attr:`regions`.
        """
        lowers = np.array([r.lower for r in self.regions])
        return np.asarray(np.searchsorted(lowers, omega, side="right") - 1)

    def region(self, omega: ArrayLike) -> np.ndarray:
        """Return the label of the region each frequency lies in.

        Parameters
        ----------
        omega : array_like
            Photon angular frequency, rad/s, finite and positive.

        Returns
        -------
        numpy.ndarray
            Labels (:data:`REGION_LABELS`), with the shape of ``omega``.

        Raises
        ------
        ArgumentError
            If a frequency is not finite and positive.
        """
        w = _frequencies(omega)
        labels = np.array([r.label for r in self.regions], dtype=object)
        return np.asarray(labels[self._region_indices(w)])

    def available(self, omega: ArrayLike) -> np.ndarray:
        """Return where the model has a value: below threshold, or on its reference.

        Parameters
        ----------
        omega : array_like
            Photon angular frequency, rad/s, finite and positive.

        Returns
        -------
        numpy.ndarray
            Booleans with the shape of ``omega``: False where :meth:`cross_section` is
            NaN.

        Raises
        ------
        ArgumentError
            If a frequency is not finite and positive.
        """
        labels = self.region(omega)
        return np.isin(labels, list(_DEFINED))

    def cross_section(self, omega: ArrayLike) -> Float[np.ndarray, " *shape"]:
        r"""Return the two-photon ionisation cross section at photon frequencies.

        Parameters
        ----------
        omega : array_like
            Photon angular frequency, rad/s, finite and positive.

        Returns
        -------
        numpy.ndarray
            :math:`\sigma^{(2)}`, m^4 s, with the shape of ``omega``: exactly 0 below
            the lowest threshold, NaN where :meth:`available` is False.

        Raises
        ------
        ArgumentError
            If a frequency is not finite and positive.

        Warns
        -----
        ModelRangeWarning
            If a frequency lies in an ``unsupported_detail`` region.
        """
        w = _frequencies(omega)
        labels = np.array([r.label for r in self.regions], dtype=object)
        where = labels[self._region_indices(w)]
        out = np.full(w.shape, np.nan)
        out[where == "below_threshold"] = 0.0
        on_curve = (where == "supported") | (where == "unsupported_detail")
        out[on_curve] = np.exp(np.interp(w[on_curve], self.omega, np.log(self.sigma)))
        if np.any(where == "unsupported_detail"):
            reasons = {
                r.reason for r in self.regions if r.label == "unsupported_detail"
            }
            warnings.warn(
                f"{self.identity.name}: values in an unsupported-detail region: "
                + " ".join(sorted(reasons)),
                ModelRangeWarning,
                stacklevel=2,
            )
        return out

    def pair_cross_section(self, pair_omega: ArrayLike) -> Float[np.ndarray, " *shape"]:
        r"""Return the cross section indexed by a photon pair's summed frequency.

        Parameters
        ----------
        pair_omega : array_like
            The sum of the two photons' angular frequencies, rad/s, finite and
            positive.

        Returns
        -------
        numpy.ndarray
            :math:`\sigma^{(2)}(\Omega/2)`, m^4 s, as :meth:`cross_section`.

        Raises
        ------
        ArgumentError
            If a frequency is not finite and positive.

        Warns
        -----
        ModelRangeWarning
            If a frequency lies in an ``unsupported_detail`` region.
        """
        return self.cross_section(_frequencies(pair_omega) / 2)

    def rate(
        self, omega: ArrayLike, intensity: ArrayLike
    ) -> Float[np.ndarray, " *shape"]:
        r"""Return the two-photon ionisation rate per atom of a monochromatic wave.

        Parameters
        ----------
        omega : array_like
            Photon angular frequency, rad/s, finite and positive.
        intensity : array_like
            Intensity, W/m^2, finite and non-negative; broadcast with ``omega``.

        Returns
        -------
        numpy.ndarray
            :math:`W = \sigma^{(2)} (I / \hbar\omega)^2`, 1/s; NaN where the
            model is unavailable.

        Raises
        ------
        ArgumentError
            If a frequency is not finite and positive, or an intensity is not finite
            and non-negative.

        Warns
        -----
        ModelRangeWarning
            If a frequency lies in an ``unsupported_detail`` region.
        """
        w = _frequencies(omega)
        i = _non_negative(intensity, "intensity")
        return np.asarray(self.cross_section(w) * (i / (HBAR * w)) ** 2)

    def loss_coefficient(
        self, omega: ArrayLike, number_density: ArrayLike
    ) -> Float[np.ndarray, " *shape"]:
        r"""Return the two-photon loss coefficient beta, dI/dz = -beta I^2.

        Parameters
        ----------
        omega : array_like
            Photon angular frequency, rad/s, finite and positive.
        number_density : array_like
            The density of ground-state atoms, m^-3, finite and non-negative;
            broadcast with ``omega``.

        Returns
        -------
        numpy.ndarray
            :math:`\beta = 2 N \sigma^{(2)} / \hbar\omega`, m/W: two photons per
            event; NaN where the model is unavailable.

        Raises
        ------
        ArgumentError
            If a frequency is not finite and positive, or a density is not finite and
            non-negative.

        Warns
        -----
        ModelRangeWarning
            If a frequency lies in an ``unsupported_detail`` region.
        """
        w = _frequencies(omega)
        n = _non_negative(number_density, "number_density")
        return np.asarray(2 * n * self.cross_section(w) / (HBAR * w))


def _frequencies(value: ArrayLike) -> np.ndarray:
    """Return frequencies as a float array, checked finite and positive.

    Parameters
    ----------
    value : array_like
        Frequencies, rad/s.

    Returns
    -------
    numpy.ndarray
        The frequencies.

    Raises
    ------
    ArgumentError
        If one is not finite and positive.
    """
    w = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(w) & (w > 0)):
        msg = "frequencies must be finite and positive, rad/s"
        raise ArgumentError(msg)
    return w


def _non_negative(value: ArrayLike, name: str) -> np.ndarray:
    """Return a quantity as a float array, checked finite and non-negative.

    Parameters
    ----------
    value : array_like
        The quantity.
    name : str
        Its name, for the error.

    Returns
    -------
    numpy.ndarray
        The quantity.

    Raises
    ------
    ArgumentError
        If a value is not finite and non-negative.
    """
    x = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(x) & (x >= 0)):
        msg = f"{name} must be finite and non-negative"
        raise ArgumentError(msg)
    return x


@functools.cache
def _selection() -> dict[str, Any]:
    """Return the selection file, parsed.

    Returns
    -------
    dict
        Its tables, by gas.
    """
    with _SELECTION_FILE.open("rb") as f:
        return tomllib.load(f)


def _plain(value: object) -> Any:
    """Return a value as plain JSON values (an infinite edge as null).

    Parameters
    ----------
    value : object
        A field's value.

    Returns
    -------
    object
        Dicts, lists, strings, numbers, booleans or None.
    """
    if is_dataclass(value) and not isinstance(value, type):
        value = {
            f.name: getattr(value, f.name)
            for f in fields(value)
            if f.name not in _STANDING
        }
    if isinstance(value, Mapping):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list, np.ndarray)):
        return [_plain(v) for v in value]
    if isinstance(value, Gas):
        return value.value
    if isinstance(value, (float, np.floating)):
        return None if math.isinf(value) else float(value)
    return value


def _optional(value: Any) -> float | None:
    """Return a stored optional number.

    Parameters
    ----------
    value : object
        A number or None.

    Returns
    -------
    float or None
        The number.
    """
    return None if value is None else float(value)


def _dumps(model: TwoPhotonIonisationModel) -> str:
    """Serialise a model as JSON text, as stored.

    Parameters
    ----------
    model : TwoPhotonIonisationModel
        The model.

    Returns
    -------
    str
        Indented JSON.
    """
    return json.dumps(_plain(model), indent=1, allow_nan=False) + "\n"


def _object(cls: type, raw: Any, where: str) -> Mapping[str, Any]:
    """Return a stored object's fields, checking they are exactly the class's.

    Parameters
    ----------
    cls : type
        The dataclass.
    raw : object
        The JSON object.
    where : str
        For error messages.

    Returns
    -------
    Mapping
        The fields.

    Raises
    ------
    DataError
        If it is not an object with exactly the class's fields.
    """
    names = {f.name for f in fields(cls)} - _STANDING
    if not isinstance(raw, Mapping) or set(raw) != names:
        msg = f"{where}: a {cls.__name__} needs exactly the fields {sorted(names)}"
        raise DataError(msg)
    return raw


def _threshold(raw: Any, where: str) -> TwoPhotonThreshold:
    """Read a stored threshold.

    Parameters
    ----------
    raw : object
        The JSON object.
    where : str
        For error messages.

    Returns
    -------
    TwoPhotonThreshold
        The threshold.
    """
    d = _object(TwoPhotonThreshold, raw, where)
    return TwoPhotonThreshold(
        name=str(d["name"]),
        omega=float(d["omega"]),
        photons=int(d["photons"]),
        energy=float(d["energy"]),
        datasets=tuple(str(k) for k in d["datasets"]),
    )


def _loads(text: str, where: str) -> TwoPhotonIonisationModel:
    """Read a stored model from its JSON text.

    Parameters
    ----------
    text : str
        The JSON, as :func:`_dumps` writes it.
    where : str
        For error messages.

    Returns
    -------
    TwoPhotonIonisationModel
        The model, with the status its selection declares.

    Raises
    ------
    DataError
        If the text is not a model of the current derivation version.
    """
    d = _object(TwoPhotonIonisationModel, json.loads(text), where)
    if d["derivation_version"] != DERIVATION_VERSION:
        msg = (
            f"{where}: derivation version {d['derivation_version']}, not "
            f"{DERIVATION_VERSION}; re-derive it"
        )
        raise DataError(msg)
    ref = _object(TwoPhotonReference, d["reference"], where)
    unc = _object(TwoPhotonUncertainty, d["uncertainty"], where)
    gas = as_gas(d["gas"])
    polarization = str(d["polarization"])
    status, reason = parse_status(
        _selection()[gas.value][polarization], where, default="validated"
    )
    return TwoPhotonIonisationModel(
        gas=gas,
        polarization=polarization,
        derivation_version=int(d["derivation_version"]),
        omega=np.array(d["omega"], dtype=float),
        sigma=np.array(d["sigma"], dtype=float),
        thresholds=tuple(_threshold(t, where) for t in d["thresholds"]),
        e1=_threshold(d["e1"], where),
        regions=tuple(
            TwoPhotonRegion(
                label=str(r["label"]),
                lower=float(r["lower"]),
                upper=math.inf if r["upper"] is None else float(r["upper"]),
                reason=str(r["reason"]),
            )
            for r in (_object(TwoPhotonRegion, r, where) for r in d["regions"])
        ),
        reference=TwoPhotonReference(
            datasets=tuple(str(k) for k in ref["datasets"]),
            gauge=str(ref["gauge"]),
            method=str(ref["method"]),
            intensity_w_cm2=(
                None
                if ref["intensity_w_cm2"] is None
                else float(ref["intensity_w_cm2"])
            ),
            reason=str(ref["reason"]),
            registers=tuple(str(r) for r in ref["registers"]),
        ),
        uncertainty=TwoPhotonUncertainty(
            source=_optional(unc["source"]),
            source_level=str(unc["source_level"]),
            source_provenance=str(unc["source_provenance"]),
            digitising=float(unc["digitising"]),
            digitising_provenance=str(unc["digitising_provenance"]),
            spreads=tuple(
                TwoPhotonSpread(
                    label=str(s["label"]),
                    datasets=tuple(str(k) for k in s["datasets"]),
                    form=str(s["form"]),
                    role=str(s["role"]),
                    lower=float(s["lower"]),
                    upper=float(s["upper"]),
                    points=int(s["points"]),
                    lower_ratio=float(s["lower_ratio"]),
                    upper_ratio=float(s["upper_ratio"]),
                )
                for s in (_object(TwoPhotonSpread, s, where) for s in unc["spreads"])
            ),
        ),
        input_digests={str(k): str(v) for k, v in d["input_digests"].items()},
        status=status,
        status_reason=reason,
    )


@functools.cache
def _stored(gas: Gas, polarization: str) -> TwoPhotonIonisationModel:
    """Return a stored model.

    Parameters
    ----------
    gas : Gas
        The species.
    polarization : str
        ``"linear"`` or ``"circular"``.

    Returns
    -------
    TwoPhotonIonisationModel
        The model of ``data/gas/models/two_photon/<gas>_<polarization>.json``.

    Raises
    ------
    DataError
        If the file is missing, malformed, or holds another model.
    """
    path = _DIRECTORY / f"{gas.value}_{polarization}.json"
    if not path.is_file():
        msg = f"{path.name}: the stored two-photon model is missing; derive it"
        raise DataError(msg)
    model = _loads(path.read_text(encoding="utf-8"), path.name)
    if (model.gas, model.polarization) != (gas, polarization):
        msg = f"{path.name} holds the model of {model.gas.value} ({model.polarization})"
        raise DataError(msg)
    return model


def two_photon_ionisation_model(
    gas: Gas | str, polarization: str = "linear"
) -> TwoPhotonIonisationModel:
    """Return a noble gas's two-photon ionisation model (OEP 6, Phase 1).

    Parameters
    ----------
    gas : Gas or str
        The gas: He, Ne, Ar or Xe (Kr has none).
    polarization : str, optional
        ``"linear"`` (default), or ``"circular"`` where an approved reference gives it
        (Xe).

    Returns
    -------
    TwoPhotonIonisationModel
        The model.

    Raises
    ------
    UnavailableModelError
        For a gas or polarisation with no approved reference, with the reason.
    ArgumentError
        For a polarisation other than linear or circular.

    Examples
    --------
    >>> xenon = two_photon_ionisation_model("Xe", "circular")
    >>> xenon.reference.datasets
    ('lagutin_2017/fig8_cipfcp_circular',)
    """
    if polarization not in POLARIZATIONS:
        msg = f"polarization must be one of {POLARIZATIONS}, not {polarization!r}"
        raise ArgumentError(msg)
    species = as_gas(gas)
    entry = _selection().get(species.value)
    if entry is None:
        msg = f"{species.value}: two-photon ionisation is curated for the noble gases"
        raise UnavailableModelError(msg)
    if "missing" in entry:
        raise UnavailableModelError(f"{species.value}: {entry['missing']}")
    if polarization not in entry:
        msg = (
            f"{species.value}: no approved {polarization} reference (OEP 6, Q4); a "
            "factor read from a raster figure is not a model"
        )
        raise UnavailableModelError(msg)
    return _stored(species, polarization)


def two_photon_ionisation_models() -> dict[tuple[Gas, str], TwoPhotonIonisationModel]:
    """Return every two-photon ionisation model, by gas and polarisation.

    Returns
    -------
    dict
        ``{(gas, polarization): model}`` for every approved reference.

    Examples
    --------
    >>> sorted(f"{g.value}/{p}" for g, p in two_photon_ionisation_models())
    ['Ar/linear', 'He/linear', 'Ne/linear', 'Xe/circular', 'Xe/linear']
    """
    return {
        (as_gas(gas), polarization): _stored(as_gas(gas), polarization)
        for gas, entry in _selection().items()
        for polarization in POLARIZATIONS
        if polarization in entry
    }
