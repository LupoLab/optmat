# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The refractivity virial coefficient of each gas, from the curated data.

The refractivity of a real gas departs from Lorentz-Lorenz through the pair term of the
virial expansion, ``(n^2 - 1)/(n^2 + 2) = A_R rho + B_R rho^2 + ...``: a relative
correction ``B_R rho / A_R`` of order 1e-6 (helium) to 1e-4 (SF6) at 1 atm, and of order
1e-3 to 1e-2 at the tens of bar of gas-filled fibres. Which measured ``B_R`` each gas
uses is chosen in ``data/gas/refractivity_virial.toml``, which points at curated
datasets; the values themselves come from those datasets.
"""

import functools
import math
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Final, cast

from optmat import _catalogue
from optmat.exceptions import DataError
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.gases import Gas, as_gas

__all__ = [
    "VIRIAL_FILE",
    "LikeTerm",
    "ThirdTerm",
    "like_term",
    "pair_rules",
    "third_term",
    "virial_ratio",
    "virial_sources",
]

VIRIAL_FILE: Final = GAS_DATA_DIRECTORY / "refractivity_virial.toml"
"""The per-gas choice of refractivity virial coefficient."""


@functools.cache
def virial_sources() -> dict[Gas, tuple[str, str]]:
    """Return, per gas, the dataset keys of its ``B_R`` and ``A_R``.

    Returns
    -------
    dict of Gas to tuple of str
        ``(B_R key, A_R key)`` for every gas with a pair term.

    Raises
    ------
    DataError
        If an entry names a dataset that does not exist.
    """
    with VIRIAL_FILE.open("rb") as handle:
        record = tomllib.load(handle)
    known = _catalogue.by_key()
    choices = {}
    for name, entry in record.items():
        # A formula's b (an estimate, `b_formula`) has no B_R and A_R datasets.
        if name == "pairs" or "b_formula" in entry:
            continue
        for key in (entry["B_R"], entry["A_R"]):
            if key not in known:
                msg = f"refractivity_virial.toml [{name}]: unknown dataset {key!r}"
                raise DataError(msg)
        choices[as_gas(name)] = (entry["B_R"], entry["A_R"])
    return choices


@functools.cache
def _formula_b(gas: Gas | str) -> tuple[float, float, str, float, str] | None:
    """Return a gas's b from a published formula, where its entry names one.

    The formula ``b = (b_eps + b_R2 sigma^2) [1 + b_theta (T - 303 K)]`` (cm^3/mol,
    sigma = 1/lambda in um^-1) is evaluated at the entry's ``wavelength_nm`` and
    ``temperature_k``.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    tuple or None
        ``(b, sigma, formula key, temperature, wavelength)``, b and sigma in m^3/mol;
        None where the gas's entry names no formula.

    Raises
    ------
    DataError
        If the formula is not a curated one.
    """
    species = as_gas(gas)
    entry = _record().get(species.value, {})
    if "b_formula" not in entry:
        return None
    source, formula_id = str(entry["b_formula"]).split("/", 1)
    found = [f for f in _catalogue.source(source).formulas if f["id"] == formula_id]
    if not found:
        msg = (
            f"refractivity_virial.toml [{species.value}]: no formula "
            f"{entry['b_formula']!r}"
        )
        raise DataError(msg)
    coefficients = cast("Mapping[str, str]", found[0]["coefficients"])
    wavenumber = 1e3 / float(entry["wavelength_nm"])  # um^-1
    temperature = float(entry["temperature_k"])
    b = (float(coefficients["b_eps"]) + float(coefficients["b_R2"]) * wavenumber**2) * (
        1 + float(coefficients["b_theta"]) * (temperature - 303.0)
    )
    b *= 1e-6  # cm^3/mol -> m^3/mol
    sigma = abs(b) * float(entry["relative_uncertainty"])
    return (
        b,
        sigma,
        str(entry["b_formula"]),
        temperature,
        f"{float(entry['wavelength_nm']):g}",
    )


@functools.cache
def virial_ratio(gas: Gas | str) -> float:
    """Return ``B_R / A_R`` for a gas, m^3/mol (0 where none is curated).

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    float
        The ratio the Lorentz-Lorenz functions take as ``virial``.

    Examples
    --------
    >>> f"{virial_ratio('Ar') * 1e6:.3f}"
    '0.419'
    >>> virial_ratio("D2")
    0.0
    """
    formula = _formula_b(gas)
    if formula is not None:
        return formula[0]
    choice = virial_sources().get(as_gas(gas))
    if choice is None:
        return 0.0
    by_key = _catalogue.by_key(gas)
    b_r = float(by_key[choice[0]].scaled_values()[0]) * 1e-12  # cm^6/mol^2 -> m^6/mol^2
    a_r = float(by_key[choice[1]].scaled_values()[0]) * 1e-6  # cm^3/mol -> m^3/mol
    return b_r / a_r


@functools.cache
def _record() -> dict[str, Any]:
    """Return ``refractivity_virial.toml``.

    Returns
    -------
    dict
        The file's tables.
    """
    with VIRIAL_FILE.open("rb") as handle:
        return tomllib.load(handle)


@functools.cache
def pair_rules() -> Mapping[str, Any]:
    """Return the rules for the unlike pairs of a mixture (OEP 4, section 7.2).

    Returns
    -------
    Mapping
        ``rule`` (``"mean"``), ``reason``, ``uncertainty``, ``unknown_b``
        (cm^3/mol), ``third_ratio`` and ``density_fraction``.
    """
    return dict(_record()["pairs"])


@dataclass(frozen=True)
class LikeTerm:
    """A gas's own pair term, ``B_R = b A_R``, as a mixture takes it.

    Attributes
    ----------
    gas : str
        The species.
    b : float
        ``B_R / A_R``, m^3/mol; 0 for an uncurated term.
    sigma : float
        Its standard uncertainty, m^3/mol; for an uncurated term the declared scale.
    status : str
        ``curated``, ``provisional`` (a published estimate, with its stated
        uncertainty: O2's, egan_2026) or ``unknown``.
    source : str
        The ``B_R`` dataset or formula, or empty.
    temperature : float
        The temperature it was measured at, K; NaN where there is none.
    frequency : str
        ``static``, or the vacuum wavelength in nm.
    contribution : str
        ``induced``: an optical value, or a static one whose induced part is stated.
    """

    gas: str
    b: float
    sigma: float
    status: str
    source: str
    temperature: float
    frequency: str
    contribution: str


@functools.cache
def like_term(gas: Gas | str) -> LikeTerm:
    """Return a gas's pair term with its uncertainty, or its declared scale.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    LikeTerm
        The curated ``b`` of :func:`virial_ratio`, or ``unknown`` with ``b = 0`` and
        :func:`pair_rules`' ``unknown_b`` as its scale.

    Raises
    ------
    DataError
        If a static ``B_R`` does not state that it is induced.

    Examples
    --------
    >>> like_term("D2").status
    'unknown'
    >>> like_term("O2").status
    'provisional'
    >>> f"{like_term('Ar').b * 1e6:.3f}"
    '0.419'
    """
    species = as_gas(gas)
    formula = _formula_b(species)
    if formula is not None:
        b, sigma, key, temperature, frequency = formula
        return LikeTerm(
            species.value,
            b,
            sigma,
            "provisional",
            key,
            temperature,
            frequency,
            "induced",
        )
    choice = virial_sources().get(species)
    if choice is None:
        scale = float(pair_rules()["unknown_b"]) * 1e-6
        return LikeTerm(species.value, 0.0, scale, "unknown", "", math.nan, "", "")
    by_key = _catalogue.by_key(species)
    b_r, a_r = by_key[choice[0]], by_key[choice[1]]
    a_value = float(a_r.scaled_values()[0]) * 1e-6
    sigma = float(b_r.standard_uncertainty()[0]) * 1e-12 / a_value
    wavelength = float(b_r.abscissae[0])
    contribution = _record()[species.value].get("contribution", "")
    if math.isinf(wavelength) and contribution != "induced":
        msg = (
            f"refractivity_virial.toml [{species.value}]: a static B_R is an optical "
            "term only where its induced part is stated (contribution = 'induced')"
        )
        raise DataError(msg)
    return LikeTerm(
        gas=species.value,
        b=virial_ratio(species),
        sigma=sigma,
        status="curated",
        source=choice[0],
        temperature=float(str(b_r.parameters.get("temperature_k", math.nan))),
        frequency="static" if math.isinf(wavelength) else f"{wavelength:g}",
        contribution="induced",
    )


@dataclass(frozen=True)
class ThirdTerm:
    """The bound on a gas's third refractivity virial term, for the optical domain.

    Attributes
    ----------
    gas : str
        The species.
    ratio : float
        ``|C_R| / A_R``, m^6/mol^2; NaN for a gas with neither ``B_R`` nor ``C_R``,
        which has no third-term bound (OEP 4, section 4.3).
    status : str
        ``measured`` (a curated C_R), ``bound`` (a measured bound), ``scale`` (the
        declared ``third_ratio B_R^2 / A_R``) or ``none``.
    source : str
        The C_R dataset, or what gives the bound.
    """

    gas: str
    ratio: float
    status: str
    source: str


@functools.cache
def third_term(gas: Gas | str) -> ThirdTerm:
    """Return the bound on a gas's third refractivity virial term (OEP 4, 4.3).

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    ThirdTerm
        From its measured C_R, a measured bound, the declared scale on its B_R, or
        none.

    Examples
    --------
    >>> term = third_term("Ar")
    >>> term.status, f"{term.ratio * 1e12:.2f}"
    ('measured', '20.64')
    """
    species = as_gas(gas)
    entry = _record().get(species.value, {})
    by_key = _catalogue.by_key(species)
    formula = _formula_b(species)
    if formula is not None:
        scale = float(pair_rules()["third_ratio"]) * formula[0] ** 2
        return ThirdTerm(species.value, scale, "scale", formula[2])
    if species not in virial_sources():
        return ThirdTerm(species.value, math.nan, "none", "")
    b_key, a_key = virial_sources()[species]
    a_value = float(by_key[a_key].scaled_values()[0]) * 1e-6
    if "C_R" in entry:
        c_value = float(by_key[entry["C_R"]].scaled_values()[0]) * 1e-18
        return ThirdTerm(
            species.value, abs(c_value) / a_value, "measured", entry["C_R"]
        )
    if "C_R_bound" in entry:
        bound = float(entry["C_R_bound"]) * 1e-18
        return ThirdTerm(
            species.value, bound / a_value, "bound", "refractivity_virial.toml"
        )
    b_value = float(by_key[b_key].scaled_values()[0]) * 1e-12
    scale = float(pair_rules()["third_ratio"]) * b_value**2 / a_value**2
    return ThirdTerm(species.value, scale, "scale", b_key)
