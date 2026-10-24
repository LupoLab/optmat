# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The models fitted here, loaded from ``data/gas/models/fitted/<name>.json``.

Each file is written by ``tools/fit/run.py`` from a fit configuration
(``data/gas/fits/<name>.toml``) and records the model's terms, their covariance, the
fitted dataset scales, the chi-squared per dataset and the configuration's SHA-256.
"""

import functools
import json
import math
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, ClassVar

import numpy as np
from jaxtyping import Float

from optmat import _validity
from optmat.constants import AU_ANGULAR_FREQUENCY, C
from optmat.exceptions import DataError
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.absorption import MeasuredAbsorption, measured_absorption
from optmat.gas.gases import as_gas
from optmat.gas.infrared import (  # pyright: ignore[reportPrivateUsage]
    _HARTREE_WAVENUMBER,
    bound_cells,
)
from optmat.gas.models._defaults import curated_status
from optmat.gas.models.oscillator import OscillatorModel

__all__ = ["FittedModel", "TemperatureCoefficient", "fitted_model", "fitted_models"]


@dataclass(frozen=True)
class TemperatureCoefficient:
    """A measured temperature coefficient of a gas polarizability (OEP 1, section 4.9).

    Attributes
    ----------
    coefficient_per_k : float
        ``A_theta``, per K, in ``A_R(T) = A_R(T0) [1 + A_theta (T - T0)]``.
    sigma_per_k : float
        Its standard uncertainty, per K.
    wavelength : float
        The vacuum wavelength it is measured at, m.
    dataset : str
        The curated dataset it comes with.
    reason : str
        Where it is printed.
    """

    coefficient_per_k: float
    sigma_per_k: float
    wavelength: float
    dataset: str
    reason: str


@functools.cache
def _temperatures() -> dict[str, Any]:
    """Return ``data/gas/temperatures.toml``.

    Returns
    -------
    dict
        Its tables.
    """
    with (GAS_DATA_DIRECTORY / "temperatures.toml").open("rb") as handle:
        return tomllib.load(handle)


@dataclass(frozen=True)
class FittedModel(OscillatorModel):
    """An oscillator model fitted here, with its covariance and statistics.

    Its overrides add two inputs to the oscillator model's: ``infrared``, a pair of
    sequences (strengths, energies in hartree) that replaces the fixed infrared terms
    whole, and ``infrared_scale``, a factor on their strengths. The term counts are
    fixed otherwise. An override drops the covariance and the chi-squared (NaN): the
    fit's uncertainty does not carry over.

    Attributes
    ----------
    n_lines : int
        How many leading terms are resonance lines at fixed energies.
    n_fixed : int
        How many trailing terms are fixed (the infrared bands), outside the covariance.
    covariance : numpy.ndarray
        Covariance of the strengths and of the effective terms' energies.
    reduced_chi2 : float
        Chi-squared per degree of freedom of the fit.
    record : dict
        The stored record, for reports.
    """

    n_lines: int = 0
    n_fixed: int = 0
    covariance: Float[np.ndarray, "q q"] = field(
        default_factory=lambda: np.zeros((0, 0)), repr=False
    )
    reduced_chi2: float = 0.0
    record: Mapping[str, object] = field(default_factory=dict, repr=False)

    EVALUATION_VERSION: ClassVar[int] = 1
    INPUTS: ClassVar[tuple[str, ...]] = (
        *OscillatorModel.INPUTS,
        "infrared",
        "infrared_scale",
    )
    DERIVED: ClassVar[tuple[str, ...]] = (
        *OscillatorModel.DERIVED,
        "covariance",
        "n_fixed",
        "n_lines",
        "reduced_chi2",
    )

    @property
    def gap_reasons(self) -> tuple[str, ...]:
        """The reason of each gap between two qualified intervals (OEP 3).

        Returns
        -------
        tuple of str
            ``absorption``, ``representation``, ``representation: kept band`` or
            ``width``; empty for a fit with one interval.
        """
        return tuple(str(r) for r in self.record.get("gap_reasons", ()))  # pyright: ignore[reportGeneralTypeIssues]

    @property
    def n_rotational(self) -> int:
        """How many of the fixed terms are the rotational system (OEP 3).

        Returns
        -------
        int
            They follow the fitted terms; 0 for a fit without the derived rules.
        """
        return int(str(self.record.get("n_rotational", 0)))

    @property
    def derived_sha256(self) -> str:
        """The SHA-256 of the derived tables the fit read, or empty (OEP 3).

        Returns
        -------
        str
            :func:`optmat.gas.infrared.derived_digest` when the fit was made.
        """
        return str(self.record.get("derived_sha256", ""))

    def _gap_note(self, wavelength: Float[np.ndarray, "*shape"]) -> str:
        """Return the reason of the gap the first wavelength between intervals is in.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.

        Returns
        -------
        str
            ``"; <low>-<high> nm is a gap (<reason>)"``, or empty.
        """
        reasons = self.gap_reasons
        values = np.ravel(np.asarray(wavelength, dtype=float))
        for (_, low), (high, _), reason in zip(
            self.intervals[:-1], self.intervals[1:], reasons, strict=False
        ):
            if np.any((values > low) & (values < high)):
                return f"; {low * 1e9:.5g}-{high * 1e9:.5g} nm is a gap ({reason})"
        return ""

    def representation_bound_au(
        self, omega: Float[np.ndarray, "*shape"], order: int
    ) -> Float[np.ndarray, "*shape"]:
        """Return the bound on the representation error of the fixed infrared terms.

        ``B_m(omega) >= |D_m[model] - D_m[reference]|`` (OEP 3, section 4) from the
        gas's derived tables: order 0 bounds the polarizability, 2 and 3 the dispersion
        ``d^m (omega alpha) / d omega^m``. It is defined in every qualified interval.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energies, hartree.
        order : int
            0, 2 or 3.

        Returns
        -------
        numpy.ndarray
            The bound, atomic units (those of alpha times E_h^(1 - m) for order m),
            NaN outside every qualified interval.

        Raises
        ------
        ValueError
            For another order.
        DataError
            For a fit that read no derived tables.

        Warns
        -----
        ModelRangeWarning
            Outside every qualified interval.
        """
        if order not in (0, 2, 3):
            msg = f"the representation bound is of order 0, 2 or 3, not {order}"
            raise ValueError(msg)
        if not self.derived_sha256:
            msg = f"{self.name} read no derived tables: it has no representation bound"
            raise DataError(msg)
        omega = np.asarray(omega, dtype=float)
        wavelength = np.asarray(2 * np.pi * C / (omega * AU_ANGULAR_FREQUENCY))
        self._check_range(wavelength)
        cells = bound_cells(self.gas)
        wavenumber = np.ravel(omega) * _HARTREE_WAVENUMBER
        index = np.searchsorted(cells["high_cm"], wavenumber, side="left")
        index = np.minimum(index, len(cells["high_cm"]) - 1)
        inside = (wavenumber >= cells["low_cm"][index]) & (
            wavenumber <= cells["high_cm"][index]
        )
        qualified = ~_validity.outside(np.ravel(wavelength), self.intervals)
        values = np.where(inside & qualified, cells[f"b{order}"][index], np.nan)
        return values.reshape(omega.shape)

    def _override(self, changes: Mapping[str, Any]) -> "FittedModel":
        """Return the fit with some inputs changed.

        Parameters
        ----------
        changes : Mapping
            Inputs to new values.

        Returns
        -------
        FittedModel
            The rebuilt model.

        Raises
        ------
        ValueError
            If the terms change in number other than through ``infrared``.
        """
        changed = dict(changes)
        n = len(self.strengths) - self.n_fixed
        strengths = list(changed.pop("strengths", self.strengths))
        energies = list(changed.pop("energies", self.energies))
        labels = list(changed.pop("labels", self.labels))
        if not len(strengths) == len(energies) == len(labels) == len(self.strengths):
            msg = f"{self.name}: the number of terms is fixed; use infrared="
            raise ValueError(msg)
        n_fixed = self.n_fixed
        if "infrared" in changed:
            band_strengths, band_energies = changed.pop("infrared")
            n_fixed = len(band_strengths)
            strengths = strengths[:n] + [float(f) for f in band_strengths]
            energies = energies[:n] + [float(e) for e in band_energies]
            labels = labels[:n] + [f"infrared {k + 1}" for k in range(n_fixed)]
        scale = float(changed.pop("infrared_scale", 1.0))
        strengths = strengths[:n] + [f * scale for f in strengths[n:]]
        return super()._override(
            changed
            | {
                "strengths": tuple(strengths),
                "energies": tuple(energies),
                "labels": tuple(labels),
                "n_fixed": n_fixed,
            }
        )

    def _stated_uncertainty(self) -> Mapping[str, Any]:
        """Return the fit's uncertainty, which an override drops.

        Returns
        -------
        Mapping
            The covariance and the reduced chi-squared, NaN.
        """
        return {
            "covariance": np.full(np.shape(self.covariance), math.nan),
            "reduced_chi2": math.nan,
        }

    @property
    def reference_temperature(self) -> float | None:
        """The temperature the fit's polarizability refers to, K, where it declares one.

        Returns
        -------
        float or None
            293.15 K for N2, O2, CO2, N2O, H2 and D2; None for the atoms and for CH4,
            CO and SF6, whose levels rest on data at other temperatures without a
            curated coefficient (TC-730). Evaluation does not depend on it.
        """
        value = self.record.get("reference_temperature_k")
        return None if value is None else float(value)  # pyright: ignore[reportArgumentType]

    @property
    def infrared_temperature(self) -> float | None:
        """The temperature of the fit's fixed infrared terms, K.

        Returns
        -------
        float or None
            296 K for HITRAN's line intensities (298.1 K for SF6's cross-section);
            None for a fit without infrared terms.
        """
        if self.n_fixed == 0:
            return None
        table = _temperatures()["infrared"]
        return float(table.get(str(self.gas), table["default_k"]))

    @property
    def temperature_coefficient(self) -> TemperatureCoefficient | None:
        """The measured temperature coefficient of the polarizability, if curated.

        Returns
        -------
        TemperatureCoefficient or None
            For N2, O2 and CO2 (TC-730); None otherwise.
        """
        entry = _temperatures()["coefficients"].get(str(self.gas))
        if entry is None:
            return None
        return TemperatureCoefficient(
            coefficient_per_k=float(entry["coefficient_per_k"]),
            sigma_per_k=float(entry["sigma_per_k"]),
            wavelength=float(entry["wavelength_nm"]) * 1e-9,
            dataset=str(entry["dataset"]),
            reason=str(entry["reason"]),
        )

    def _identity_temperatures(self) -> tuple[float | None, float | None]:
        """Return the reference temperature and the infrared terms' temperature.

        Returns
        -------
        tuple
            K, or None where there is none.
        """
        return self.reference_temperature, self.infrared_temperature

    def _identity_extras(self) -> Mapping[str, Any]:
        """Return the virial ratio and the temperatures the fit refers to.

        Returns
        -------
        Mapping
            ``refractivity_virial``, ``reference_temperature``,
            ``infrared_temperature`` and ``temperature_coefficient``.
        """
        extras = dict(super()._identity_extras())
        extras["reference_temperature"] = self.reference_temperature
        extras["infrared_temperature"] = self.infrared_temperature
        extras["temperature_coefficient"] = self.temperature_coefficient
        return extras

    def _sources(self) -> tuple[str, ...]:
        """Return the sources of the datasets the fit was fitted to.

        Returns
        -------
        tuple of str
            Sorted source ids.
        """
        fitted = self.record.get("dataset_chi2", {})
        assert isinstance(fitted, Mapping)
        return tuple(sorted({str(key).split("/", 1)[0] for key in fitted}))

    def polarizability_sigma_au(
        self, omega: Float[np.ndarray, " n"], *, birge: bool = True
    ) -> Float[np.ndarray, " n"]:
        """Return the one-sigma uncertainty of the polarizability.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energies, hartree.
        birge : bool, optional
            Scale the covariance by the reduced chi-squared when it exceeds one.

        Returns
        -------
        numpy.ndarray
            One-sigma uncertainty, atomic units.
        """
        n_fitted = len(self.strengths) - self.n_fixed
        strengths = np.asarray(self.strengths)[:n_fitted]
        energies = np.asarray(self.energies)[:n_fitted]
        denominator = energies[:, None] ** 2 - np.asarray(omega) ** 2
        gradient = np.vstack(
            [
                1.0 / denominator,
                (-2.0 * strengths[:, None] * energies[:, None] / denominator**2)[
                    self.n_lines :
                ],
            ]
        )
        covariance = self.covariance
        if birge and self.reduced_chi2 > 1:
            covariance = covariance * self.reduced_chi2
        return np.sqrt(np.einsum("in,ij,jn->n", gradient, covariance, gradient))

    def dispersion_sigma_au(
        self, omega: Float[np.ndarray, " n"], order: int, *, birge: bool = True
    ) -> Float[np.ndarray, " n"]:
        """Return the absolute one-sigma uncertainty of a dispersion derivative.

        ``sqrt(g^T Sigma g)``, with ``g`` the gradient of
        ``d^k (omega alpha) / d omega^k``
        (:meth:`~optmat.gas.models.oscillator.OscillatorModel.dispersion_au`) with
        respect to the fitted parameters and ``Sigma`` their covariance. It stays
        finite where the derivative crosses zero (OEP 2, section 4). The fixed terms
        carry none.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energies, hartree.
        order : int
            The order ``k >= 2``.
        birge : bool, optional
            Scale the covariance by the reduced chi-squared when it exceeds one.

        Returns
        -------
        numpy.ndarray
            One-sigma uncertainty, atomic units (those of alpha times
            ``E_h^(1 - k)``).

        Raises
        ------
        ValueError
            If ``order`` is below 2.
        """
        if order < 2:
            msg = f"order must be at least 2, got {order}"
            raise ValueError(msg)
        n_fitted = len(self.strengths) - self.n_fixed
        strengths = np.asarray(self.strengths)[:n_fitted, None]
        energies = np.asarray(self.energies)[:n_fitted, None]
        omega = np.asarray(omega)[None, :]
        # omega f / (E^2 - omega^2) = (f / 2) [1 / (E - omega) - 1 / (E + omega)].
        k, sign = order, (-1) ** order
        factor = 0.5 * math.factorial(k)
        per_strength = factor * (
            (energies - omega) ** -(k + 1) - sign * (energies + omega) ** -(k + 1)
        )
        per_energy = (
            factor
            * strengths
            * (k + 1)
            * (
                -((energies - omega) ** -(k + 2))
                + sign * (energies + omega) ** -(k + 2)
            )
        )
        gradient = np.vstack([per_strength, per_energy[self.n_lines :]])
        covariance = self.covariance
        if birge and self.reduced_chi2 > 1:
            covariance = covariance * self.reduced_chi2
        return np.sqrt(np.einsum("in,ij,jn->n", gradient, covariance, gradient))

    def dispersion_sigma_relative(
        self, omega: Float[np.ndarray, " n"], order: int, *, birge: bool = True
    ) -> Float[np.ndarray, " n"]:
        """Return the relative one-sigma uncertainty of a dispersion coefficient.

        In the dilute limit the propagation constant is
        ``beta = (omega / c) [1 + N alpha / (2 epsilon_0)]``, so for ``k >= 2`` each
        ``beta_k = d^k beta / d omega^k`` is proportional to
        ``d^k (omega alpha) / d omega^k`` and its relative uncertainty is that of the
        derivative, independent of the density (``k = 2`` is the group-velocity
        dispersion, ``k = 3`` the third-order dispersion). Lorentz-Lorenz changes the
        relative uncertainty only at the order of the susceptibility (below 1e-3 of
        itself at a few bar). The magnetic term is frequency-independent and drops out.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energies, hartree, below the lowest term's energy.
        order : int
            The order ``k >= 2`` of the coefficient.
        birge : bool, optional
            Scale the covariance by the reduced chi-squared when it exceeds one.

        Returns
        -------
        numpy.ndarray
            ``sigma(beta_k) / |beta_k|``: :meth:`dispersion_sigma_au` over
            ``|dispersion_au|``.

        Raises
        ------
        ValueError
            If ``order`` is below 2.
        """
        sigma = self.dispersion_sigma_au(omega, order, birge=birge)
        return sigma / np.abs(self.dispersion_au(np.asarray(omega), order))


def _absorption(record: Mapping[str, object]) -> MeasuredAbsorption | None:
    """Rebuild a stored fit's measured-absorption term.

    Parameters
    ----------
    record : mapping
        The stored record.

    Returns
    -------
    MeasuredAbsorption or None
        The term at the recorded smoothing, or None.

    Raises
    ------
    DataError
        If the rebuilt term's pieces differ from those the fit used.
    """
    stored = record.get("absorption")
    if not isinstance(stored, Mapping):
        return None
    absorption = measured_absorption(
        str(record["gas"]),
        fwhm_cm=float(stored["fwhm_cm"]),  # pyright: ignore[reportArgumentType]
    )
    pieces = [[p.source, p.from_nm, p.to_nm] for p in absorption.pieces]
    if pieces != stored["pieces"]:
        msg = f"{record['name']}: the absorption pieces changed since the fit; refit"
        raise DataError(msg)
    return absorption


@functools.cache
def fitted_model(name: str) -> FittedModel:
    """Return one stored fitted model, reading only its record.

    Parameters
    ----------
    name : str
        The model's name, ``fit/<name>``.

    Returns
    -------
    FittedModel
        The model.

    Raises
    ------
    KeyError
        If no fit is stored under that name.
    DataError
        If the stored record names another fit.
    """
    stem = name.removeprefix("fit/")
    path = GAS_DATA_DIRECTORY / "models" / "fitted" / f"{stem}.json"
    if stem == name or not path.is_file():
        raise KeyError(name)
    record = json.loads(path.read_text(encoding="utf-8"))
    if record["name"] != stem:
        msg = f"{path.name}: the fit is named {record['name']!r}, not its file name"
        raise DataError(msg)
    # A fit with the derived windows (OEP 3) records its qualified intervals;
    # otherwise its one interval is range_nm.
    spans = record.get("intervals_nm", [record["range_nm"]])
    intervals = tuple((float(a) * 1e-9, float(b) * 1e-9) for a, b in spans)
    # The covariance covers every fitted strength and each effective term's
    # energy; the fixed (infrared) terms come last and are outside it.
    n_fixed = int(record.get("n_fixed", 0))
    n_terms = len(record["labels"]) - n_fixed
    n_lines = n_terms - (len(record["covariance"]) - n_terms)
    gas = as_gas(record["gas"])
    status, reason = curated_status(gas, name)
    return FittedModel(
        gas=gas,
        name=name,
        status=status,
        status_reason=reason,
        citation=f"fitted here ({record['name']}.toml, sha256 "
        f"{record['config_sha256'][:12]})",
        wavelength_range=(intervals[0][0], intervals[-1][1]),
        intervals=intervals,
        kind="fitted",
        strengths=tuple(record["strengths"]),
        energies=tuple(record["energies_hartree"]),
        labels=tuple(record["labels"]),
        offset_au=float(record.get("chi0_au", 0.0)),
        retardation_au=float(record.get("retardation_au", 0.0)),
        absorption=_absorption(record),
        n_lines=n_lines,
        n_fixed=n_fixed,
        covariance=np.array(record["covariance"]),
        reduced_chi2=record["chi2"] / max(record["dof"], 1),
        record=record,
    )


@functools.cache
def fitted_models() -> dict[str, FittedModel]:
    """Return every stored fitted model, by name.

    Returns
    -------
    dict of str to FittedModel
        Models named ``fit/<name>``.
    """
    paths = sorted((GAS_DATA_DIRECTORY / "models" / "fitted").glob("*.json"))
    return {f"fit/{p.stem}": fitted_model(f"fit/{p.stem}") for p in paths}
