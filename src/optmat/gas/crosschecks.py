# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Test the recommended models against dispersion data that no single fit uses.

Two kinds of measurement constrain dispersion without giving a refractive index:

- **Phase-matching density ratios** (Shelton 1986, 1990): two gases give the same
  electric-field-induced second-harmonic wave-vector mismatch at densities in the ratio
  ``rho_gas / rho_ref = Delta alpha_ref / Delta alpha_gas`` with
  ``Delta alpha = alpha(2 omega) - alpha(omega)`` (to first order in density). Each
  ratio ties the near-infrared dispersion of two gases together to about 0.1 %, which is
  a direct test of the neon, krypton and xenon models where they have no near-infrared
  data of their own. The fits use most of them (TC-96); here every one is compared.
- **Group-delay dispersion** (Borzsonyi et al. 2008): ``beta_2 = (lambda^3 / 2 pi c^2)
  d^2 n / d lambda^2`` per metre at the measured state.

This module predicts each from the recommended models (``models.registry``) and returns
the comparison.
"""

from dataclasses import dataclass

import numpy as np

from optmat.constants import C
from optmat.datasets import datasets
from optmat.exceptions import OptmatError
from optmat.gas.models.registry import default_model

__all__ = ["Comparison", "group_delay_dispersions", "phase_match_ratios"]


@dataclass(frozen=True)
class Comparison:
    """A measured value against the models' prediction.

    Attributes
    ----------
    key : str
        The dataset.
    row : str
        The row within it.
    wavelength : float
        The (fundamental) vacuum wavelength, m.
    measured : float
        The measured value.
    sigma : float
        Its one-sigma uncertainty (NaN where none is stated).
    predicted : float
        The models' prediction.
    """

    key: str
    row: str
    wavelength: float
    measured: float
    sigma: float
    predicted: float

    @property
    def deviation(self) -> float:
        """The relative deviation of the prediction from the measurement.

        Returns
        -------
        float
            ``predicted / measured - 1``.
        """
        return self.predicted / self.measured - 1


def phase_match_ratios() -> list[Comparison]:
    """Compare every curated phase-matching density ratio with the models.

    Returns
    -------
    list of Comparison
        One per row of every ``phase_match_density_ratio`` dataset whose two gases have
        recommended models; the prediction is the ratio of the models' polarizability
        differences between the fundamental and its second harmonic.
    """
    comparisons = []
    for dataset in datasets(quantity="phase_match_density_ratio"):
        if dataset.use not in {"primary", "check"}:
            continue
        try:
            gas_model = default_model(dataset.material)
            reference_model = default_model(
                str(dataset.parameters["reference_material"])
            )
        except OptmatError:
            continue
        harmonic = int(dataset.parameters.get("harmonic", 2))  # pyright: ignore[reportArgumentType]
        wavelength = dataset.vacuum_wavelength()
        sigma = dataset.standard_uncertainty()
        for i, row in enumerate(dataset.rows):
            points = np.array([wavelength[i], wavelength[i] / harmonic])
            gas_alpha = gas_model.polarizability(points, warn=False)
            reference_alpha = reference_model.polarizability(points, warn=False)
            predicted = (reference_alpha[1] - reference_alpha[0]) / (
                gas_alpha[1] - gas_alpha[0]
            )
            comparisons.append(
                Comparison(
                    key=dataset.key,
                    row=row,
                    wavelength=float(wavelength[i]),
                    measured=float(dataset.scaled_values()[i]),
                    sigma=float(sigma[i]),
                    predicted=float(predicted),
                )
            )
    return comparisons


def group_delay_dispersions() -> list[Comparison]:
    """Compare every curated group-delay dispersion with the models.

    Returns
    -------
    list of Comparison
        One per row of every ``group_delay_dispersion`` dataset (fs^2/m at the
        dataset's state), predicted by a centred second difference of the recommended
        model's ``n`` in wavelength.
    """
    comparisons = []
    for dataset in datasets(quantity="group_delay_dispersion"):
        try:
            model = default_model(dataset.material)
        except OptmatError:
            continue
        if dataset.state is None:  # pragma: no cover - the format requires a state
            continue
        wavelength = dataset.vacuum_wavelength()
        sigma = dataset.standard_uncertainty()
        for i, row in enumerate(dataset.rows):
            lam = float(wavelength[i])
            step = lam * 1e-3
            n = model.refractivity(
                np.array([lam - step, lam, lam + step]), dataset.state, warn=False
            )
            second = (n[0] - 2 * n[1] + n[2]) / step**2
            predicted = lam**3 / (2 * np.pi * C**2) * second * 1e30  # fs^2/m
            comparisons.append(
                Comparison(
                    key=dataset.key,
                    row=row,
                    wavelength=lam,
                    measured=float(dataset.scaled_values()[i]),
                    sigma=float(sigma[i]),
                    predicted=float(predicted),
                )
            )
    return comparisons
