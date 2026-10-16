# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The fitted models of solids, loaded from their stored records.

A fit (:mod:`optmat.solid.fitting`, ``tools/fit/run_solid.py``) stores its full UDM
parameter set and the covariance of its free parameters in
``data/solid/models/fitted/<name>.json``. The model is named ``fit/<name>``, and its
uncertainty in n and k follows by linear propagation.
"""

import functools
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, ClassVar

import numpy as np
from jaxtyping import Float
from numpy.typing import ArrayLike

from optmat.solid import SOLID_DATA_DIRECTORY
from optmat.solid.models._defaults import curated
from optmat.solid.models.base import complex_index
from optmat.solid.models.formulas import UdmModel
from optmat.solid.udm import from_parameters

__all__ = ["FittedModel", "fitted_models"]


@dataclass(frozen=True)
class FittedModel(UdmModel):
    """A UDM fitted here, with the covariance of its free parameters.

    Attributes
    ----------
    attributes : tuple of str
        The UDM's attributes.
    values : dict of str to float
        Every parameter.
    free : tuple of str
        The fitted parameters.
    scales : dict of str to tuple of str
        The fitted groups of parameters, each scaled by one factor.
    covariance : numpy.ndarray
        Their covariance (with the datasets' offsets after them).
    reduced_chi2 : float
        The fit's reduced chi-squared.

    Its overrides are the parameters (``values``, any subset, merged into the fit's)
    and the intervals; the UDM is rebuilt from them, and the covariance and the
    chi-squared are dropped (NaN).
    """

    attributes: tuple[str, ...] = ()
    values: Mapping[str, float] = field(default_factory=dict, repr=False)
    free: tuple[str, ...] = ()
    scales: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    covariance: Float[np.ndarray, "p p"] = field(
        default_factory=lambda: np.zeros((0, 0)), repr=False
    )
    reduced_chi2: float = math.nan

    EVALUATION_VERSION: ClassVar[int] = 1
    INPUTS: ClassVar[tuple[str, ...]] = (*UdmModel.INPUTS, "values")
    DERIVED: ClassVar[tuple[str, ...]] = (
        *UdmModel.DERIVED,
        "covariance",
        "free",
        "reduced_chi2",
        "scales",
        "udm",
    )

    def _override(self, changes: Mapping[str, Any]) -> "FittedModel":
        """Return the fit with some parameters changed, its UDM rebuilt.

        Parameters
        ----------
        changes : Mapping
            Inputs to new values; ``values`` names any of the parameters.

        Returns
        -------
        FittedModel
            The rebuilt model.

        Raises
        ------
        ValueError
            If ``values`` names a parameter the UDM does not have.
        """
        changed = dict(changes)
        if "values" in changed:
            unknown = sorted(set(changed["values"]) - set(self.values))
            if unknown:
                msg = f"{self.name} has no parameter {', '.join(unknown)}"
                raise ValueError(msg)
            values = dict(self.values) | {
                k: float(v) for k, v in changed["values"].items()
            }
            changed |= {
                "values": values,
                "udm": from_parameters(self.attributes, values),
            }
        return super()._override(changed)

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

    def index_sigma(
        self, wavelength: ArrayLike | float
    ) -> tuple[Float[np.ndarray, "*shape"], Float[np.ndarray, "*shape"]]:
        """Return the one-sigma uncertainties of n and k from the fit's covariance.

        Linear propagation with a central-difference Jacobian (step 1e-6 of each
        parameter's value): the parameters' own uncertainty only, not the model's
        form.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.

        Returns
        -------
        tuple of numpy.ndarray
            The uncertainties sigma_n and sigma_k.
        """
        omega = self._omega(wavelength, warn=False)
        energy = np.atleast_1d(self._energy(omega))

        def index(values: Mapping[str, float]) -> np.ndarray:
            """Return n + ik at ``energy`` for a parameter set.

            Parameters
            ----------
            values : mapping
                The parameters.

            Returns
            -------
            numpy.ndarray
                The complex index.
            """
            return complex_index(
                from_parameters(self.attributes, values).permittivity(energy)
            )

        # One column per free parameter (a step of 1e-6 of its value), then one per
        # group scale (a relative step of 1e-6 of all its members at once): the
        # covariance's order and units.
        columns = []
        for name in self.free:
            step = 1e-6 * max(abs(self.values[name]), 1e-12)
            high, low = dict(self.values), dict(self.values)
            high[name] += step
            low[name] -= step
            columns.append((index(high) - index(low)) / (2 * step))
        for members in self.scales.values():
            high, low = dict(self.values), dict(self.values)
            for name in members:
                high[name] *= 1 + 1e-6
                low[name] *= 1 - 1e-6
            columns.append((index(high) - index(low)) / 2e-6)
        jac = np.array(columns).T
        size = len(columns)
        cov = self.covariance[:size, :size]
        sigma_n = np.sqrt(np.einsum("ij,jk,ik->i", np.real(jac), cov, np.real(jac)))
        sigma_k = np.sqrt(np.einsum("ij,jk,ik->i", np.imag(jac), cov, np.imag(jac)))
        shape = np.shape(omega)
        return sigma_n.reshape(shape), sigma_k.reshape(shape)


@functools.cache
def fitted_models() -> dict[str, FittedModel]:
    """Return every stored solid fit as a model, named ``fit/<name>``.

    Returns
    -------
    dict of str to FittedModel
        The models.
    """
    out = {}
    for path in sorted((SOLID_DATA_DIRECTORY / "models" / "fitted").glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        attributes = tuple(record["attributes"])
        values = {str(k): float(v) for k, v in record["values"].items()}
        low, high = record["range_nm"]
        name = f"fit/{record['name']}"
        out[name] = FittedModel(
            material=record["material"],
            name=name,
            citation=(f"optmat fit {record['name']} (TC-322), from {record['start']}"),
            wavelength_range=(low * 1e-9, high * 1e-9),
            kind="fitted",
            udm=from_parameters(attributes, values),
            attributes=attributes,
            values=values,
            free=tuple(record["free"]),
            scales={k: tuple(v) for k, v in record.get("scales", {}).items()},
            covariance=np.asarray(record["covariance"], dtype=float),
            reduced_chi2=record["chi2"] / record["dof"],
        )
    return curated(out)
