# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Golden values of every model of the public surface, with their identities.

    uv run python tools/golden_models.py

writes ``tests/data/golden_models.json``: for each model of the registries (the gas
dispersion, air, transport, hyperpolarizability, rate and Raman-response models; the
equations of state, pure, air and mixed; the solid, Kerr and silica Raman models; the
hyperpolarizability support records; the two-photon ionisation models and bound-level
records), a few values at fixed arguments and the model's identity digest.
Gas Raman identities also go into ``tests/data/golden_raman_identities.json``, once
per response. Their derived widths depend on numerical arithmetic: the tests compare
their canonical records with a rounding tolerance, retaining exact checks on versions
and discrete metadata. Other models retain the exact digest check.
``tests/test_golden_models.py`` independently compares the evaluated values. A change
beyond these arithmetic tolerances requires review and re-recording, with a changelog
line and an ``EVALUATION_VERSION`` bump where evaluation changed.
"""

import json
import math
import warnings
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import numpy as np

from optmat.conditions import State
from optmat.exceptions import ModelAbsorptionWarning, ModelRangeWarning
from optmat.gas.air import AIR_COMPOSITIONS, air_model
from optmat.gas.eos import air_equation_of_state, equation_of_state
from optmat.gas.gases import Gas, info
from optmat.gas.ionisation_rates import ionisation_rate_models
from optmat.gas.mixture import mixture_optics
from optmat.gas.mixture_eos import mixture_equation_of_state
from optmat.gas.models.registry import _all_models as gas_models
from optmat.gas.molecular_ionisation_rates import molecular_rate_models
from optmat.gas.nonlinear import HyperpolarizabilitySupport, hyperpolarizability_models
from optmat.gas.raman_response import GasRamanResponse, raman_response
from optmat.gas.transport import TRANSPORT_GASES, transport_model
from optmat.gas.two_photon import two_photon_ionisation_models
from optmat.gas.two_photon_bound import two_photon_bound_levels
from optmat.solid.bundle import silica_bundles
from optmat.solid.kerr import kerr_models
from optmat.solid.models.registry import _all_models as solid_models
from optmat.solid.raman_response import silica_raman_responses

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "tests" / "data" / "golden_models.json"
RAMAN_OUTPUT = ROOT / "tests" / "data" / "golden_raman_identities.json"
POINTS = 5
# The longest wavelength sampled where a model's range is open to infinity, m.
LONGEST = 20e-6
RATE_WAVELENGTH = 800e-9
FIELDS = (2e10, 5e10, 1e11, 2e11)  # V/m
RAMAN_GASES = ("N2", "O2", "H2", "D2", "CH4", "SF6", "CO2", "N2O", "dry air")
RAMAN_TEMPERATURE = 295.0
RAMAN_DENSITY = 2.5e25  # m^-3, about 1 atm
TIMES = (1e-14, 1e-13, 1e-12, 1e-11)  # s
OMEGAS = (1e12, 1e13, 5e13, 1e14)  # rad/s
KERR_WAVELENGTHS = (400e-9, 800e-9, 1064e-9, 1550e-9)
NONLINEAR_WAVELENGTHS = (800e-9, 1064e-9, 1550e-9)
AIR_WAVELENGTHS = (400e-9, 633e-9, 800e-9, 1.2e-6)
TRANSPORT_TEMPERATURES = (250.0, 294.0, 350.0, 500.0)  # K, inside every gas's ranges
# A support record's numbers (OEP 5), m^-1 (max_w2 m^-2); NaN where a limit is missing.
SUPPORT_NUMBERS = (
    "max_w2",
    "leg_max",
    "pair_max",
    "e1",
    "e2",
    "e2_any",
    "ionisation",
    "leg_max_in_w2",
    "band_edge/linear",
    "band_edge/any",
)
MIXTURES = (("He", "Ar"), ("Ne", "Ar"), ("He", "Ne"), ("N2", "O2", "Ar"))
MIXTURE_DENSITY = 1e26  # m^-3 of each species: inside every mixture's optical domain
MIXTURE_COMPOSITIONS = (
    ("He+Ar", {"He": 0.5, "Ar": 0.5}),
    ("He+Ne", {"He": 0.5, "Ne": 0.5}),
    ("Ne+Ar", {"Ne": 0.5, "Ar": 0.5}),
    ("dry air", "dry air (N2, O2, Ar)"),
    ("CH4+N2", {"CH4": 0.5, "N2": 0.5}),
)
MIXTURE_PRESSURES = (1e5, 1e6, 2e6)  # Pa: within every default's truncation bound


def _encode(values: Any) -> list[Any]:
    """Return values as JSON lists, a complex number as ``[re, im]``."""
    array = np.asarray(values)
    if np.iscomplexobj(array):
        return np.stack([array.real, array.imag], axis=-1).tolist()
    return array.astype(float).tolist()


def _grid(low: float, high: float) -> list[float]:
    """Return ``POINTS`` wavelengths spread geometrically inside a range."""
    high = min(high, max(LONGEST, 2 * low))
    return list(np.geomspace(low * (1 + 1e-3), high * (1 - 1e-3), POINTS))


def _support_numbers(record: HyperpolarizabilitySupport) -> list[float]:
    """Return a support record's numbers in the order of ``SUPPORT_NUMBERS``."""
    limits = [
        record.leg_max,
        record.pair_max,
        record.e1,
        record.e2,
        record.e2_any,
        record.ionisation,
    ]
    edges = [record.band_edge("linear"), record.band_edge("any")]
    return [
        record.max_w2,
        *(math.nan if x is None else x.wavenumber for x in limits),
        record.leg_max_in_w2,
        *(math.nan if x is None else x for x in edges),
    ]


def _entry(
    identifier: str, obj: object, quantity: str, arguments: dict[str, Any], values: Any
) -> dict[str, Any]:
    """Return one golden record."""
    return {
        "id": identifier,
        "class": f"{type(obj).__module__}.{type(obj).__qualname__}",
        "quantity": quantity,
        "arguments": arguments,
        "values": _encode(values),
    }


def _mixture_cases() -> Iterator[
    tuple[str, object, str, dict[str, Any], Callable[[], Any]]
]:
    """Yield the golden values of the mixtures' optics and equations of state."""
    for species in MIXTURES:
        optics = mixture_optics(species)
        densities = dict.fromkeys(optics.species, MIXTURE_DENSITY)
        for quantity in ("susceptibility", "dsusceptibility_domega"):
            yield (
                f"mixture/{'+'.join(optics.species)}",
                optics,
                quantity,
                {"wavelength": list(AIR_WAVELENGTHS), "density": MIXTURE_DENSITY},
                lambda optics=optics, q=quantity, d=densities: getattr(optics, q)(
                    np.asarray(AIR_WAVELENGTHS), d
                ),
            )
    for name, composition in MIXTURE_COMPOSITIONS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ModelRangeWarning)
            eos = mixture_equation_of_state(composition)
        yield (
            f"mixture_eos/{name}",
            eos,
            "number_density",
            {"pressure": list(MIXTURE_PRESSURES), "temperature": 294.0},
            lambda eos=eos: [
                eos.number_density(State(pressure=p, temperature=294.0))
                for p in MIXTURE_PRESSURES
            ],
        )
        yield (
            f"mixture_eos/{name}",
            eos,
            "molar_mass",
            {},
            lambda eos=eos: eos.molar_mass,
        )


def _eos_cases() -> Iterator[
    tuple[str, object, str, dict[str, Any], Callable[[], Any]]
]:
    """Yield the molar masses of the pure gases' equations of state and of air's."""
    equations = [equation_of_state(g) for g in Gas if info(g).fluid is not None]
    for eos in [*equations, air_equation_of_state()]:
        yield (
            f"eos/{eos.material}",
            eos,
            "molar_mass",
            {},
            lambda eos=eos: eos.molar_mass,
        )


def _cases() -> Iterator[tuple[str, object, str, dict[str, Any], Callable[[], Any]]]:
    """Yield ``(id, model, quantity, arguments, evaluate)`` for every golden value."""
    for name, gas in sorted(gas_models().items()):
        # The first qualified interval: a molecule with OEP 3's windows has poles
        # (its bands) between its intervals; every other model has one, its range.
        grid = _grid(*gas.intervals[0])
        yield (
            f"gas/{name}",
            gas,
            "complex_polarizability",
            {"wavelength": grid},
            lambda gas=gas, grid=grid: gas.complex_polarizability(grid, warn=False),
        )
        yield (
            f"gas/{name}",
            gas,
            "dpolarizability_domega",
            {"wavelength": grid},
            lambda gas=gas, grid=grid: gas.dpolarizability_domega(grid, warn=False),
        )
    for name in AIR_COMPOSITIONS:
        air = air_model(name)
        for quantity in ("complex_polarizability", "dpolarizability_domega"):
            yield (
                f"air/{name}",
                air,
                quantity,
                {"wavelength": list(AIR_WAVELENGTHS)},
                lambda air=air, q=quantity: getattr(air, q)(
                    np.asarray(AIR_WAVELENGTHS), warn=False
                ),
            )
    yield from _eos_cases()
    yield from _mixture_cases()
    for gas in TRANSPORT_GASES:
        model = transport_model(gas)
        for quantity in ("viscosity", "heat_capacity"):
            yield (
                f"transport/{gas}",
                model,
                quantity,
                {"temperature": list(TRANSPORT_TEMPERATURES)},
                lambda model=model, q=quantity: getattr(model, q)(
                    np.asarray(TRANSPORT_TEMPERATURES)
                ),
            )
    for name, bundle in sorted(silica_bundles().items()):
        yield (
            f"silica_bundle/{name}",
            bundle,
            "n2_electronic",
            {"wavelength": list(KERR_WAVELENGTHS)},
            lambda bundle=bundle: bundle.n2_electronic(np.asarray(KERR_WAVELENGTHS)),
        )
    for name, solid in sorted(solid_models().items()):
        grid = _grid(*solid.wavelength_range)
        for quantity in ("permittivity", "dpermittivity_domega"):
            yield (
                f"solid/{name}",
                solid,
                quantity,
                {"wavelength": grid},
                lambda solid=solid, grid=grid, q=quantity: getattr(solid, q)(
                    grid, warn=False
                ),
            )
    for gas, hyper in sorted(
        hyperpolarizability_models().items(), key=lambda kv: kv[0].value
    ):
        for process in ("static", "kerr", "thg"):
            yield (
                f"hyperpolarizability/{gas.value}",
                hyper,
                f"gamma_au/{process}",
                {"wavelength": list(NONLINEAR_WAVELENGTHS)},
                lambda hyper=hyper, process=process: hyper.gamma_au(
                    process, np.asarray(NONLINEAR_WAVELENGTHS)
                ),
            )
        support = hyper.support
        assert support is not None
        yield (
            f"hyperpolarizability_support/{gas.value}",
            support,
            "numbers",
            {"names": list(SUPPORT_NUMBERS)},
            lambda support=support: _support_numbers(support),
        )
    for (gas, polarization), model in sorted(
        two_photon_ionisation_models().items(),
        key=lambda kv: (kv[0][0].value, kv[0][1]),
    ):
        # Five frequencies across the first supported band (OEP 6), rad/s.
        band = next(r for r in model.regions if r.label == "supported")
        omegas = list(np.linspace(band.lower, band.upper, POINTS + 2)[1:-1])
        yield (
            f"two_photon/{gas.value}/{polarization}",
            model,
            "cross_section",
            {"omega": omegas},
            lambda model=model, omegas=omegas: model.cross_section(np.asarray(omegas)),
        )
    for gas in ("He", "Ne", "Ar", "Kr", "Xe"):
        record = two_photon_bound_levels(gas)
        yield (
            f"two_photon_bound/{gas}",
            record,
            "sigma0",
            {"strengths": [f"{s.dataset}#{s.row}" for s in record.strengths]},
            lambda record=record: [
                math.nan if s.sigma0 is None else s.sigma0 for s in record.strengths
            ],
        )
    yield from _response_cases()


def _response_cases() -> Iterator[
    tuple[str, object, str, dict[str, Any], Callable[[], Any]]
]:
    """Yield the golden values of the rates, the Raman responses and silica's models."""
    for family, registry in (
        ("rate", ionisation_rate_models()),
        ("molecular_rate", molecular_rate_models()),
    ):
        for form, by_gas in sorted(registry.items()):
            for gas, base in sorted(by_gas.items(), key=lambda kv: kv[0].value):
                rate = base.at_wavelength(RATE_WAVELENGTH)
                for quantity in ("instantaneous_rate", "cycle_averaged_rate"):
                    yield (
                        f"{family}/{form}/{gas.value}@800nm",
                        rate,
                        quantity,
                        {"field": list(FIELDS)},
                        lambda rate=rate, q=quantity: getattr(rate, q)(
                            np.asarray(FIELDS)
                        ),
                    )
    for species in RAMAN_GASES:
        response = raman_response(
            species, RAMAN_TEMPERATURE, wavelength=RATE_WAVELENGTH
        )
        yield (
            f"raman_response/{species}@295K,800nm",
            response,
            "response",
            {"time": list(TIMES), "density": RAMAN_DENSITY},
            lambda r=response: r.response(np.asarray(TIMES), RAMAN_DENSITY),
        )
        yield (
            f"raman_response/{species}@295K,800nm",
            response,
            "decay_rates",
            {"density": RAMAN_DENSITY},
            lambda r=response: r.decay_rates(RAMAN_DENSITY),
        )
        yield (
            f"raman_response/{species}@295K,800nm",
            response,
            "frequencies",
            {"density": RAMAN_DENSITY},
            lambda r=response: r.frequencies(RAMAN_DENSITY),
        )
    for name, kerr in sorted(kerr_models().items()):
        yield (
            name,
            kerr,
            "n2",
            {"wavelength": list(KERR_WAVELENGTHS)},
            lambda kerr=kerr: kerr.n2(np.asarray(KERR_WAVELENGTHS)),
        )
    for name, silica in sorted(silica_raman_responses().items()):
        yield (
            f"silica_raman/{name}",
            silica,
            "response",
            {"time": list(TIMES)},
            lambda s=silica: s.response(np.asarray(TIMES)),
        )
        yield (
            f"silica_raman/{name}",
            silica,
            "spectrum",
            {"omega": list(OMEGAS)},
            lambda s=silica: s.spectrum(np.asarray(OMEGAS)),
        )
        yield (
            f"silica_raman/{name}",
            silica,
            "raman_coefficient",
            {"wavelength": RATE_WAVELENGTH},
            lambda s=silica: s.raman_coefficient(RATE_WAVELENGTH),
        )


def records() -> list[dict[str, Any]]:
    """Return every golden record, evaluated now.

    Returns
    -------
    list of dict
        ``{id, class, quantity, arguments, values}`` per value, and ``digest`` where
        the model has an identity. Gas Raman records also include ``identity_record``.
    """
    out = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        warnings.simplefilter("ignore", ModelAbsorptionWarning)
        for identifier, obj, quantity, arguments, evaluate in _cases():
            entry = _entry(identifier, obj, quantity, arguments, evaluate())
            identity = getattr(obj, "identity", None)
            if identity is not None:
                entry["digest"] = identity.digest
            if isinstance(obj, GasRamanResponse):
                entry["identity_record"] = obj.identity_record()
            out.append(entry)
    return out


def _finite(value: Any) -> Any:
    """Return JSON-safe values: non-finite floats as strings."""
    if isinstance(value, list):
        return [_finite(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return repr(value)
    return value


def main() -> None:
    """Write the golden values and canonical gas Raman identities."""
    data = records()
    identities = {}
    for entry in data:
        if "identity_record" in entry:
            identities[entry["id"]] = entry.pop("identity_record")
        entry["values"] = _finite(entry["values"])
    OUTPUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    RAMAN_OUTPUT.write_text(json.dumps(identities, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {len(data)} records to {OUTPUT.relative_to(ROOT)}")
    print(f"wrote {len(identities)} identities to {RAMAN_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
