# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The silica nonlinear bundle: n2, f_R and h_R as one object (OEP 1, section 4.7).

A propagation code in the :math:`n_0 = 1` convention of the nonlinear Schrödinger
equation takes the nonlinear response

.. math::

   R(t) = (1 - f_R)\,\delta(t) + f_R\,h_R(t), \qquad \int_0^\infty h_R\,dt = 1,

with the total long-pulse :math:`n_2`: the triple :math:`(n_2, f_R, h_R)`. A bundle is
that triple, so its parts cannot be mismatched. Every bundle holds the **total**
:math:`n_2(\lambda)` of the Kerr fit (``kerr/SiO2``), the quantity the data measure,
and its shape's own :math:`f_R`: its electronic part is :math:`(1 - f_R)\,n_2`, and a
scalar and a vector run agree in the long-pulse limit.

- ``"default"``: Hollenbeck and Cantrell's shape with the Kerr fit's :math:`f_R`,
  ``paired``, copolarised only;
- ``"tensor"``: Lin and Agrawal's shape with its own :math:`f_R` (TC-734), which has
  the isotropic and anisotropic parts;
- one bundle per silica response, under the response's name.

Each carries its gain check against Stolen and Ippen's measured Stokes gain (526.0 nm,
330 cm^-1; ``stolen_1973``), a consistency check for a shape other than the one the
Kerr fit's :math:`f_R` was set with. The silica response's ``raman_coefficient`` and
``electronic_susceptibility`` include :math:`n_0^2`; a code in the :math:`n_0 = 1`
convention takes the bundle's :math:`n_2`, :math:`f_R` and :math:`h_R` instead.
"""

import functools
import math
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from typing import ClassVar, Final

import numpy as np
from jaxtyping import Float
from numpy.typing import ArrayLike

from optmat import _catalogue
from optmat.constants import C
from optmat.exceptions import DataError
from optmat.identity import Identified, source_ids
from optmat.solid.kerr import KerrModel, kerr_models
from optmat.solid.raman_response import (
    RAMAN_DIRECTORY,
    Component,
    SilicaRamanResponse,
    silica_raman_responses,
)

__all__ = [
    "EVALUATION_VERSION",
    "GainCheck",
    "SilicaNonlinearBundle",
    "silica_bundle",
    "silica_bundles",
]

EVALUATION_VERSION: Final = 1
"""The version of this module's evaluation code, in every bundle's identity."""

_GAIN_RECORD: Final = "stolen_1973/gain_amplifier"


@dataclass(frozen=True)
class GainCheck:
    """A bundle's Stokes gain against Stolen and Ippen's measurement.

    Attributes
    ----------
    pump_wavelength : float
        m.
    shift : float
        The Stokes shift, cm^-1.
    predicted : float
        The bundle's copolarised gain there, with its total n2 and f_R, m/W.
    measured : float
        Stolen and Ippen's, m/W.
    sigma : float
        Its standard uncertainty, m/W.
    reference : str
        The curated dataset.
    """

    pump_wavelength: float
    shift: float
    predicted: float
    measured: float
    sigma: float
    reference: str

    @property
    def deviation(self) -> float:
        """The predicted gain's relative departure from the measured one.

        Returns
        -------
        float
            ``predicted / measured - 1``.
        """
        return self.predicted / self.measured - 1


@dataclass(frozen=True)
class SilicaNonlinearBundle(Identified):
    """The consistent triple (n2, f_R, h_R) of fused silica.

    Attributes
    ----------
    name : str
        The bundle's name: ``default``, ``tensor`` or a silica response's.
    kerr : KerrModel
        The Kerr fit, whose total long-pulse n2 the bundle holds.
    response : SilicaRamanResponse
        The Raman shape, with its f_R.
    status : str
        ``validated`` for the default and the tensor default and their responses'
        bundles; ``comparison`` for the other shapes.
    """

    IDENTITY_KIND: ClassVar[str] = "silica_bundle"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ("optmat.solid.bundle",)

    name: str
    kerr: KerrModel
    response: SilicaRamanResponse

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the bundle's name, material and source records.

        Returns
        -------
        tuple
            ``("silica_bundle/<name>", "SiO2", sources)``.
        """
        found = source_ids(*self.response.sources.values(), *self.kerr.scales)
        return f"silica_bundle/{self.name}", "SiO2", found

    @property
    def raman_fraction(self) -> float:
        """The shape's f_R, the Raman share of the total n2.

        Returns
        -------
        float
            The fraction f_R.
        """
        return self.response.raman_fraction

    @property
    def paired(self) -> bool:
        """Whether f_R is the one the Kerr fit found with this shape.

        Returns
        -------
        bool
            The response's ``paired``.
        """
        return self.response.paired

    def n2(self, wavelength: ArrayLike | float) -> Float[np.ndarray, "*shape"]:
        """Return the total long-pulse n2, the Kerr fit's.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            The n2, m^2/W.

        Warns
        -----
        ModelRangeWarning
            Outside the Kerr fit's wavelengths.
        """
        return self.kerr.n2(wavelength)

    def n2_electronic(
        self, wavelength: ArrayLike | float
    ) -> Float[np.ndarray, "*shape"]:
        """Return the instantaneous part, ``(1 - f_R) n2``, with the shape's f_R.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            The electronic n2, m^2/W.
        """
        return np.asarray((1 - self.raman_fraction) * self.n2(wavelength))

    def n2_raman(self, wavelength: ArrayLike | float) -> Float[np.ndarray, "*shape"]:
        """Return the delayed part, ``f_R n2``.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            The Raman n2, m^2/W.
        """
        return np.asarray(self.raman_fraction * self.n2(wavelength))

    def h_r(
        self, time: ArrayLike, component: Component = "parallel"
    ) -> Float[np.ndarray, " *shape"]:
        """Return the normalised Raman response h_R (or one of its tensor parts).

        Parameters
        ----------
        time : array_like
            Delay after the impulse, s.
        component : {"parallel", "a", "b"}, optional
            The part; the tensor parts only for a shape that has them.

        Returns
        -------
        numpy.ndarray
            The response h_R, 1/s; ``parallel`` integrates to one.
        """
        return self.response.response(time, component)

    @functools.cached_property
    def gain_check(self) -> GainCheck:
        """The bundle's copolarised Stokes gain against Stolen and Ippen's.

        Returns
        -------
        GainCheck
            At their pump and shift, with the bundle's total n2 and f_R.
        """
        record = _catalogue.dataset(_GAIN_RECORD)
        pump = float(str(record.parameters["pump_nm"])) * 1e-9
        shift = float(str(record.parameters["shift_cm1"]))
        n2 = float(self.n2(pump))
        predicted = float(
            self.response.gain(pump, 2 * math.pi * C * 100 * shift, n2=n2)
        )
        # The record's gain is in cm/W: to m/W.
        return GainCheck(
            pump_wavelength=pump,
            shift=shift,
            predicted=predicted,
            measured=float(record.scaled_values()[0]) / 100,
            sigma=float(record.standard_uncertainty()[0]) / 100,
            reference=record.key,
        )


def _names() -> tuple[str, str]:
    """Return the names of the default and tensor-default responses.

    Returns
    -------
    tuple of str
        From ``data/solid/raman/SiO2.toml``.
    """
    with (RAMAN_DIRECTORY / "SiO2.toml").open("rb") as handle:
        table = tomllib.load(handle)
    return str(table["default"]), str(table["tensor"])


@functools.cache
def silica_bundles() -> Mapping[str, SilicaNonlinearBundle]:
    """Return every silica bundle, by name.

    Returns
    -------
    Mapping
        ``default``, ``tensor`` and one per silica response.
    """
    kerr = kerr_models()["kerr/SiO2"]
    responses = silica_raman_responses()
    default, tensor = _names()
    trusted = {default, tensor}
    chosen = {"default": default, "tensor": tensor} | {k: k for k in responses}
    out: dict[str, SilicaNonlinearBundle] = {}
    for name, response_name in chosen.items():
        validated = response_name in trusted
        out[name] = SilicaNonlinearBundle(
            name=name,
            kerr=kerr,
            response=responses[response_name],
            status="validated" if validated else "comparison",
            status_reason=""
            if validated
            else f"a shape kept for comparison; the default is {default}",
        )
    return out


def silica_bundle(name: str = "default") -> SilicaNonlinearBundle:
    """Return a silica bundle.

    Parameters
    ----------
    name : str, optional
        ``default`` (Hollenbeck, paired), ``tensor`` (Lin and Agrawal, their f_R) or a
        silica response's name.

    Returns
    -------
    SilicaNonlinearBundle
        The bundle.

    Raises
    ------
    DataError
        For an unknown name.

    Examples
    --------
    >>> bundle = silica_bundle()
    >>> bundle.paired, f"{bundle.raman_fraction:.3f}"
    (True, '0.219')
    >>> f"{silica_bundle('tensor').gain_check.deviation:+.3f}"
    '-0.079'
    """
    bundles = silica_bundles()
    if name not in bundles:
        msg = f"no silica bundle {name!r}; one of {sorted(bundles)}"
        raise DataError(msg)
    return bundles[name]
