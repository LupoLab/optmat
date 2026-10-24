# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Published dispersion formulas, built from the source records.

A source record's ``[[formulas]]`` entry that names a ``form`` is a model. Its
coefficients are kept **as printed**, as strings; a misprinted coefficient is corrected
by a ``corrections`` entry that names the term, the printed and the corrected value,
the reason and the evidence, exactly like a dataset erratum, and the printed value is
checked on loading.

The forms (``lambda`` the vacuum wavelength, ``s = 1/lambda`` in um^-1, ``nu`` the
vacuum wavenumber in cm^-1, ``omega`` the angular frequency in atomic units):

``dn_sigma2``
    ``n - 1 = scale (A + sum_i B_i / (C_i - s^2))``, the Peck and Bideau-Mehu form and
    refractiveindex.info's formula 6.
``dn_wavenumber2``
    ``n - 1 = scale (A + sum_i B_i / (C_i - nu^2))``, the atmospheric-science form (He
    et al. 2021, Wilmouth & Sayres 2020).
``n2_sellmeier``
    ``n^2 - 1 = A + sum_i B_i lambda^2 / (lambda^2 - C_i)`` with lambda in um
    (Borzsonyi et al. 2008, Ermolov et al. 2015), refractiveindex.info's formula 2.
``alpha_series``
    ``alpha = sum_k c_k omega^(2k)`` in atomic units, a Cauchy series of the
    polarizability itself (Bhatia & Drachman as quoted by Stone & Stejskal 2004).
``molar_refractivity_series``
    ``A_R = sum_k a_k lambda^(-2k)`` in cm^3/mol with lambda in nm (Stone & Stejskal).
``oscillator``
    ``alpha = sum_i f_i / (E_i^2 - omega^2)`` in atomic units, a pseudo-spectrum (Kumar
    and Meath); the entry's ``energy_unit`` ("hartree", the default, "eV" or
    "rydberg") says in what the ``E_i`` are printed.

The first three give ``n - 1`` (or ``n^2 - 1``) at the formula's reference ``state``;
the polarizability follows by inverse Lorentz-Lorenz with the density the source meant
(``density``: real or ideal gas). The last two are per particle and need no state.
"""

import functools
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, ClassVar, Final

import numpy as np
from jaxtyping import Complex, Float

from optmat._catalogue import formula_sources, source
from optmat.conditions import State
from optmat.constants import (
    AU_ANGULAR_FREQUENCY,
    AU_POLARIZABILITY,
    EPS0,
    N_A,
    C,
)
from optmat.datasets import parse_state
from optmat.exceptions import DataError
from optmat.gas.eos import number_density
from optmat.gas.gases import Gas, as_gas
from optmat.gas.lorentz_lorenz import (
    dpolarizability_dsusceptibility,
    polarizability_from_refractivity,
    polarizability_from_susceptibility,
)
from optmat.gas.models._defaults import curated_status
from optmat.gas.models.base import Model

__all__ = ["FORMS", "FormulaModel", "published_model", "published_models"]

FORMS: Final = frozenset(
    {
        "dn_sigma2",
        "dn_wavenumber2",
        "n2_sellmeier",
        "alpha_series",
        "molar_refractivity_series",
        "oscillator",
    }
)
"""The formula forms a source record may name."""

_PER_PARTICLE: Final = frozenset(
    {"alpha_series", "molar_refractivity_series", "oscillator"}
)

_ENERGY_UNITS: Final[dict[str, float]] = {
    "hartree": 1.0,
    "rydberg": 0.5,
    # eV per hartree, CODATA: the Hartree energy over the elementary charge.
    "eV": 1.0 / 27.211386245988,
}


@dataclass(frozen=True)
class FormulaModel(Model):
    """A published formula as a polarizability model.

    Attributes
    ----------
    form : str
        One of :data:`FORMS`.
    terms : dict of str to tuple of float
        The coefficients after corrections: ``A``, ``B``, ``C``, ``scale`` for the
        refractivity forms, ``c`` or ``a`` for the series.
    state : State or None
        The reference state of a refractivity form.
    ideal : bool
        True if the source related its refractivity to the state with the ideal gas.
    energy_scale : float
        Hartree per printed energy unit, for the ``oscillator`` form.

    Its overrides are the coefficients (``terms``) and the intervals.
    """

    form: str = ""
    terms: Mapping[str, tuple[float, ...]] = field(default_factory=dict, repr=False)
    state: State | None = None
    ideal: bool = False
    energy_scale: float = 1.0

    EVALUATION_VERSION: ClassVar[int] = 1
    INPUTS: ClassVar[tuple[str, ...]] = ("intervals", "terms")

    def _term(self, name: str, default: float = 0.0) -> tuple[float, ...]:
        """Return a coefficient list, or ``(default,)`` if the formula has none.

        Parameters
        ----------
        name : str
            The term.
        default : float, optional
            The value of an absent scalar term.

        Returns
        -------
        tuple of float
            The coefficients.
        """
        return self.terms.get(name, (default,))

    def _polarizability(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Float[np.ndarray, "*shape"]:
        """Evaluate the formula and convert to a polarizability.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            Polarizability, C m^2/V.
        """
        if self.form in _PER_PARTICLE:
            return self._per_particle(wavelength)
        if self.state is None:  # pragma: no cover - rejected when the model is built
            msg = f"{self.name}: a refractivity formula needs a state"
            raise DataError(msg)
        density = self.reference_density
        (a,), b, c = self._term("A"), self._term("B"), self._term("C")
        (scale,) = self._term("scale", 1.0)
        if self.form == "n2_sellmeier":
            lam2 = (wavelength * 1e6) ** 2
            with np.errstate(invalid="ignore"):  # inf / inf in the static limit
                ratio = [np.where(np.isinf(lam2), 1.0, lam2 / (lam2 - ci)) for ci in c]
            chi = a + sum(bi * r for bi, r in zip(b, ratio, strict=True))
            return polarizability_from_susceptibility(
                chi, density, self.refractivity_virial
            )
        if self.form == "dn_sigma2":
            x2 = (1e-6 / wavelength) ** 2
        else:  # dn_wavenumber2
            x2 = (1e-2 / wavelength) ** 2
        dn = scale * (a + sum(bi / (ci - x2) for bi, ci in zip(b, c, strict=True)))
        return polarizability_from_refractivity(dn, density, self.refractivity_virial)

    def _dcomplex_polarizability_domega(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the formula's analytic d alpha / d omega.

        A refractivity form is differentiated in its own variable, proportional to
        ``omega`` (``1/lambda`` in um^-1 or cm^-1), and reduced through the Jacobian
        of inverse Lorentz-Lorenz at the reference density: ``d alpha / d chi`` and
        ``d chi / d(n - 1) = 2 n`` (OEP 1, section 4.1).

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m (``inf`` for the static limit, where it is 0).

        Returns
        -------
        numpy.ndarray
            C m^2 V^-1 s (real, as complex).
        """
        if self.form in _PER_PARTICLE:
            return np.asarray(self._dper_particle(wavelength), dtype=complex)
        (a,), b, c = self._term("A"), self._term("B"), self._term("C")
        (scale,) = self._term("scale", 1.0)
        if self.form == "n2_sellmeier":
            u = 1e-6 / wavelength  # um^-1, proportional to omega
            chi = a + sum(bi / (1 - ci * u**2) for bi, ci in zip(b, c, strict=True))
            dchi = sum(
                bi * 2 * ci * u / (1 - ci * u**2) ** 2
                for bi, ci in zip(b, c, strict=True)
            )
            per_omega = 1e-6 / (2 * math.pi * C)
        else:
            unit = 1e-6 if self.form == "dn_sigma2" else 1e-2  # um^-1 or cm^-1
            x = unit / wavelength
            dn = scale * (
                a + sum(bi / (ci - x**2) for bi, ci in zip(b, c, strict=True))
            )
            ddn = scale * sum(
                bi * 2 * x / (ci - x**2) ** 2 for bi, ci in zip(b, c, strict=True)
            )
            chi = dn * (dn + 2)
            dchi = 2 * (1 + dn) * ddn
            per_omega = unit / (2 * math.pi * C)
        jacobian = dpolarizability_dsusceptibility(
            chi, self.reference_density, self.refractivity_virial
        )
        return np.asarray(jacobian * dchi * per_omega, dtype=complex)

    def _dper_particle(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Float[np.ndarray, "*shape"]:
        """Return d alpha / d omega of a per-particle series.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            C m^2 V^-1 s.
        """
        omega = 2 * math.pi * C / wavelength / AU_ANGULAR_FREQUENCY  # hartree
        au = AU_POLARIZABILITY / AU_ANGULAR_FREQUENCY
        if self.form == "oscillator":
            energies = np.asarray(self._term("E")) * self.energy_scale
            slope = sum(
                2 * f * omega / (e**2 - omega**2) ** 2
                for f, e in zip(self._term("f"), energies, strict=True)
            )
            return np.asarray(slope * au)
        if self.form == "alpha_series":
            slope = sum(
                2 * k * ck * omega ** (2 * k - 1)
                for k, ck in enumerate(self._term("c"))
                if k > 0
            )
            return np.asarray(slope * au + 0 * omega)
        s = 1e-9 / wavelength  # molar_refractivity_series, nm^-1
        slope = sum(
            2 * k * ak * s ** (2 * k - 1)
            for k, ak in enumerate(self._term("a"))
            if k > 0
        )
        return np.asarray(
            3 * EPS0 * slope * 1e-6 / N_A * 1e-9 / (2 * math.pi * C) + 0 * s
        )

    @functools.cached_property
    def reference_density(self) -> float:
        """The number density of the formula's reference state, resolved once.

        Returns
        -------
        float
            The density, m^-3 (ideal where the source used the ideal gas); NaN for a
            per-particle series, which has no reference state.
        """
        if self.state is None:
            return math.nan
        return number_density(self.gas, self.state, ideal=self.ideal)

    def _sources(self) -> tuple[str, ...]:
        """Return the formula's source record.

        Returns
        -------
        tuple of str
            The source, the first part of the name.
        """
        return (self.name.split("/", 1)[0],)

    def _identity_extras(self) -> Mapping[str, Any]:
        """Return the values resolved once: the reference density and the virial.

        Returns
        -------
        Mapping
            ``refractivity_virial`` and ``reference_density``.
        """
        return {
            **super()._identity_extras(),
            "reference_density": self.reference_density,
        }

    def _per_particle(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Float[np.ndarray, "*shape"]:
        """Evaluate a per-particle series.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            Polarizability, C m^2/V.
        """
        if self.form == "oscillator":
            omega2 = (2 * math.pi * C / wavelength / AU_ANGULAR_FREQUENCY) ** 2
            energies = np.asarray(self._term("E")) * self.energy_scale
            alpha = sum(
                f / (e**2 - omega2)
                for f, e in zip(self._term("f"), energies, strict=True)
            )
            return np.asarray(alpha * AU_POLARIZABILITY)
        if self.form == "alpha_series":
            omega2 = (2 * math.pi * C / wavelength / AU_ANGULAR_FREQUENCY) ** 2
            alpha = sum(ck * omega2**k for k, ck in enumerate(self._term("c")))
            return np.asarray(alpha * AU_POLARIZABILITY)
        inverse2 = (1e-9 / wavelength) ** 2  # molar_refractivity_series, lambda in nm
        molar = sum(ak * inverse2**k for k, ak in enumerate(self._term("a")))
        return np.asarray(3 * EPS0 * molar * 1e-6 / N_A)


def _corrected(
    terms: dict[str, list[str]], corrections: Sequence[Mapping[str, object]], where: str
) -> None:
    """Apply a formula's corrections in place, checking each printed coefficient.

    Parameters
    ----------
    terms : dict of str to list of str
        The printed coefficients, changed in place.
    corrections : sequence of mapping
        Entries with ``term``, ``index``, ``printed``, ``corrected``, ``reason`` and
        ``evidence``.
    where : str
        The formula, for error messages.

    Raises
    ------
    DataError
        If a correction names a missing coefficient, lacks its justification, or its
        printed value is not the recorded one.
    """
    for correction in corrections:
        missing = {"term", "index", "printed", "corrected", "reason", "evidence"}
        missing -= set(correction)
        if missing:
            msg = f"{where}: correction lacks {sorted(missing)}"
            raise DataError(msg)
        term, index = str(correction["term"]), int(correction["index"])  # pyright: ignore[reportArgumentType]
        if term not in terms or not 0 <= index < len(terms[term]):
            msg = f"{where}: correction of missing coefficient {term}[{index}]"
            raise DataError(msg)
        if terms[term][index] != correction["printed"]:
            msg = (
                f"{where}: correction's printed {term}[{index}] "
                f"{correction['printed']!r} is not the recorded {terms[term][index]!r}"
            )
            raise DataError(msg)
        terms[term][index] = str(correction["corrected"])


def _range(entry: Mapping[str, object], where: str) -> tuple[float, float]:
    """Return a formula's stated range in m.

    Parameters
    ----------
    entry : mapping
        The formula entry, with ``range_nm = [low, high]`` (``high`` may be "inf").
    where : str
        The formula, for error messages.

    Returns
    -------
    tuple of float
        The range, m.

    Raises
    ------
    DataError
        If the range is missing or malformed.
    """
    value = entry.get("range_nm")
    if not isinstance(value, Sequence) or len(value) != 2:
        msg = f"{where}: a formula with a form needs range_nm = [low, high]"
        raise DataError(msg)
    low, high = (float(v) * 1e-9 for v in value)  # pyright: ignore[reportArgumentType]
    return low, high


def _model(source_id: str, citation: str, entry: Mapping[str, object]) -> FormulaModel:
    """Build one formula model from its record entry.

    Parameters
    ----------
    source_id : str
        The source.
    citation : str
        The source's citation.
    entry : mapping
        The ``[[formulas]]`` entry.

    Returns
    -------
    FormulaModel
        The model.

    Raises
    ------
    DataError
        If the entry is malformed.
    """
    where = f"{source_id}/{entry.get('id')}"
    form = str(entry["form"])
    if form not in FORMS:
        msg = f"{where}: unknown form {form!r}"
        raise DataError(msg)
    raw = entry.get("terms")
    if not isinstance(raw, Mapping):
        msg = f"{where}: a formula with a form needs a terms table"
        raise DataError(msg)
    terms = {
        str(k): [str(x) for x in (v if isinstance(v, list) else [v])]
        for k, v in raw.items()
    }
    corrections = entry.get("corrections", [])
    _corrected(terms, corrections, where)  # pyright: ignore[reportArgumentType]
    state_spec = entry.get("state")
    state = None
    if form not in _PER_PARTICLE:
        state = parse_state(state_spec, where)
        if state is None:
            msg = f"{where}: a refractivity formula needs its state"
            raise DataError(msg)
    gas = as_gas(str(entry["material"]))
    status, reason = curated_status(gas, where)
    return FormulaModel(
        gas=gas,
        name=where,
        status=status,
        status_reason=reason,
        citation=citation,
        wavelength_range=_range(entry, where),
        kind=str(entry.get("kind", "published")),
        form=form,
        terms={k: tuple(float(x) for x in v) for k, v in terms.items()},
        state=state,
        ideal=entry.get("density") == "ideal",
        energy_scale=_ENERGY_UNITS[str(entry.get("energy_unit", "hartree"))],
    )


def published_model(name: str) -> FormulaModel:
    """Return one gas formula model, reading only its source's record.

    :func:`published_models` reads every record with formulas; a lookup by name needs
    one.

    Parameters
    ----------
    name : str
        ``<source>/<formula>``.

    Returns
    -------
    FormulaModel
        The model, the same as :func:`published_models` holds under that name.

    Raises
    ------
    KeyError
        If no gas formula has that name.

    Examples
    --------
    >>> published_model("bideau_mehu_1981/sellmeier_Ar").gas.value
    'Ar'
    """
    source_id, _, formula_id = name.partition("/")
    try:
        record = source(source_id)
    except KeyError:
        raise KeyError(name) from None
    for entry in record.formulas:
        if (
            entry.get("id") == formula_id
            and "form" in entry
            and str(entry.get("material")) in set(Gas)
        ):
            return _model(record.id, record.citation, entry)
    raise KeyError(name)


@functools.cache
def published_models() -> dict[str, FormulaModel]:
    """Return every formula model defined in the source records, by name.

    Returns
    -------
    dict of str to FormulaModel
        The models, named ``<source>/<formula>``.

    Examples
    --------
    >>> model = published_models()["bideau_mehu_1981/sellmeier_Ar"]
    >>> f"{float(model.refractivity(253.7e-9)) * 1e6:.1f}"
    '303.2'
    """
    return {
        f"{source.id}/{entry['id']}": _model(source.id, source.citation, entry)
        for source in formula_sources(frozenset(g.value for g in Gas))
        for entry in source.formulas
        # The gas forms only: a solid's formulas are optmat.solid's (TC-318).
        if "form" in entry and str(entry.get("material")) in set(Gas)
    }
