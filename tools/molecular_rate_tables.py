# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The tables of the molecular-rates page (``docs/molecular_ionisation_rates.md``).

    uv run python tools/molecular_rate_tables.py

prints the tables; ``tests/test_molecular_rates_page.py`` checks that the page holds
exactly what this script prints:

- ``molecular-rate-constants``: the constants of the default model per molecule;
- ``molecular-rate-models``: what each named model is and for which molecules;
- ``molecular-rate-forms``: the models against one another at 800 nm;
- ``molecular-rate-static``: the models at fixed nuclei against the orientation
  average of the all-electron static rates of N2, CO, H2, O2 and CH4 (Majety &
  Scrinzi; Jagau);
- ``molecular-rate-vibration``: H2's vibrational factor against Saenz's;
- ``molecular-rate-wahlstrand``: the models against Wahlstrand et al.'s yields;
- ``molecular-rate-solprop``: what optmat's ionisation energies change against
  solprop's.
"""

import json
import math
import warnings
from pathlib import Path

import numpy as np

from optmat.constants import AU_FIELD, E_CHARGE, EPS0, C
from optmat.exceptions import ModelRangeWarning
from optmat.gas.ionisation_validation import measured_yields
from optmat.gas.molecular_ionisation_rates import (
    molecular_rate_model,
    molecular_rate_models,
)
from optmat.gas.molecular_ionisation_validation import (
    orientation_averaged_rates,
    vibrational_factors,
)

ROOT = Path(__file__).resolve().parents[1]
MOLECULES = ("H2", "D2", "N2", "CO", "O2", "CO2", "N2O", "CH4", "SF6")
INTENSITIES = (5e17, 1e18, 2e18)
STATIC_MODELS = (
    "default",
    "moadk_zhao_2010",
    "moadk_zhao_2017",
    "moadk_tong_2002",
    "adk_hydrogenic",
)
STATIC_FIELDS = (0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14)
STATIC_REFERENCES = (
    ("N2", "majety_2015"),
    ("CO", "majety_2015"),
    ("CO", "jagau_2018"),
    ("H2", "majety_2015"),
    ("H2", "jagau_2018"),
    ("O2", "jagau_2018"),
    ("CH4", "jagau_2018"),
)


def _table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def constants_table() -> str:
    """Return the constants of the default model per molecule.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for gas in MOLECULES:
        models = molecular_rate_models()["default"]
        if not any(k.value == gas for k in models):
            rows.append([gas, "", "", "", "", "", "", "no default: see below"])
            continue
        model = molecular_rate_model(gas, 800e-9)
        rows.append(
            [
                gas,
                f"{model.electrons:g}",
                f"{model.tong_lin_alpha:g}",
                f"{model.ionization_energy / E_CHARGE:.3f}",
                f"{model.adiabatic_ionization_energy / E_CHARGE:.3f}",
                model.sources["structure"].split("/")[0],
                f"{model.vibrational_constant:g}" if model.vibrational_constant else "",
                f"{model.uncertainty_factor:g}",
            ]
        )
    return _table(
        [
            "Molecule",
            "N_eff",
            "alpha",
            "I_p, rate (eV)",
            "I_p, adiabatic (eV)",
            "Coefficients",
            "c (au)",
            "Factor",
        ],
        rows,
    )


def models_table() -> str:
    """Return what each named model is and the molecules it is given for.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for name, by_gas in molecular_rate_models().items():
        model = next(iter(by_gas.values()))
        gases = [g for g in MOLECULES if any(k.value == g for k in by_gas)]
        rows.append(
            [
                f"`{name}`",
                model.form,
                "static"
                if model.form in {"moadk", "adk"}
                else f"{model.coulomb} Coulomb factor",
                ", ".join(gases),
            ]
        )
    return _table(["Model", "Form", "Frequency", "Molecules"], rows)


def forms_table() -> str:
    """Return the cycle-averaged rates of the forms at 800 nm, one electron.

    Returns
    -------
    str
        A Markdown table, the rates in 1/s.
    """
    names = list(molecular_rate_models())
    rows = []
    for gas in MOLECULES:
        for intensity in INTENSITIES:
            amplitude = math.sqrt(2 * intensity / (EPS0 * C))
            cells = []
            for name in names:
                by_gas = molecular_rate_models()[name]
                if not any(k.value == gas for k in by_gas):
                    cells.append("")
                    continue
                model = molecular_rate_model(gas, 800e-9, name)
                cells.append(f"{float(model.cycle_averaged_rate(amplitude)):.2e}")
            rows.append([gas, f"{intensity / 1e4:.0e}", *cells])
    return _table(["Molecule", "I (W/cm2)", *(f"`{name}`" for name in names)], rows)


def static_table() -> str:
    """Return the static forms over the orientation-averaged all-electron rates.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for gas, reference in STATIC_REFERENCES:
        first = orientation_averaged_rates(gas, STATIC_MODELS[0], reference)
        field = first.field / AU_FIELD
        for label, values in (
            ("reference: rate at the first angle over its average", first.anisotropy),
            *(
                (f"`{name}`", orientation_averaged_rates(gas, name, reference).ratio)
                for name in STATIC_MODELS
                if gas in {g.value for g in molecular_rate_models()[name]}
            ),
        ):
            cells = []
            for wanted in STATIC_FIELDS:
                near = np.flatnonzero(np.abs(field / wanted - 1) < 1e-6)
                cells.append("" if near.size == 0 else f"{values[near[0]]:.2f}")
            rows.append([gas, reference, label, *cells])
    return _table(
        ["Molecule", "Reference", "", *(f"{f:g} au" for f in STATIC_FIELDS)], rows
    )


def vibration_table() -> str:
    """Return H2's vibrational factor against Saenz's.

    Returns
    -------
    str
        A Markdown table.
    """
    comparison = vibrational_factors()
    rows = [
        [
            f"{field / AU_FIELD:g}",
            f"{expected:.3f}",
            f"{predicted:.3f}",
            f"{predicted / expected:.3f}",
        ]
        for field, expected, predicted in zip(
            comparison.field, comparison.expected, comparison.predicted, strict=True
        )
    ]
    return _table(["F (au)", "Saenz", "exp(c/F^2)", "Model / Saenz"], rows)


def wahlstrand_table() -> str:
    """Return the models over Wahlstrand et al.'s measured yields of N2 and O2.

    Returns
    -------
    str
        A Markdown table: the models' yields of a 42 fs pulse at 800 nm over the
        measured fit, and the default's over the measurement divided by argon's (N2)
        and xenon's (O2) at the same intensities.
    """
    names = ("default", "moadk_zhao_2010", "ppt_talebpour", "ppt_talebpour_ilkov")
    rows = []
    for gas, atom in (("N2", "Ar"), ("O2", "Xe")):
        comparisons = [measured_yields(gas, name) for name in names]
        intensities = comparisons[0].x
        # The atom's model with the static Coulomb factor, as the molecules' (TC-848).
        reference = measured_yields(
            atom, "default_static", intensities=list(intensities)
        )
        for k in range(0, len(intensities), 3):
            rows.append(
                [
                    gas,
                    f"{intensities[k] / 1e16:.0f}",
                    f"{comparisons[0].expected[k]:.2e}",
                    *(f"{c.ratio[k]:.2f}" for c in comparisons),
                    f"{comparisons[0].ratio[k] / reference.ratio[k]:.2f} ({atom})",
                ]
            )
    return _table(
        [
            "Molecule",
            "I (TW/cm2)",
            "Measured yield",
            *(f"`{name}`" for name in names),
            "`default`, over the atom's `default_static`",
        ],
        rows,
    )


def solprop_table() -> str:
    """Return what optmat's ionisation energies change against solprop's.

    Returns
    -------
    str
        A Markdown table.
    """
    fixture = json.loads(
        (ROOT / "tests" / "data" / "solprop_molecular_rates.json").read_text()
    )
    rows = []
    for gas, entry in fixture["molecules"].items():
        model = molecular_rate_model(gas, 800e-9, "adk_hydrogenic")
        other = model.with_ionization_energy(entry["potential_j"], "solprop's")
        ratios = [
            float(model.instantaneous_rate(f) / other.instantaneous_rate(f))
            for f in (2e10, 4e10)
        ]
        shift = (model.ionization_energy - other.ionization_energy) / E_CHARGE
        rows.append(
            [
                gas,
                f"{model.ionization_energy / E_CHARGE:.4f}",
                f"{other.ionization_energy / E_CHARGE:.4f}",
                f"{1e3 * shift:+.1f}",
                f"{ratios[0]:.3f}",
                f"{ratios[1]:.3f}",
            ]
        )
    return _table(
        [
            "Molecule",
            "I_p, optmat (eV)",
            "I_p, solprop (eV)",
            "Difference (meV)",
            "Rate ratio at 2e10 V/m",
            "at 4e10 V/m",
        ],
        rows,
    )


TABLES = {
    "molecular-rate-constants": constants_table,
    "molecular-rate-models": models_table,
    "molecular-rate-forms": forms_table,
    "molecular-rate-static": static_table,
    "molecular-rate-vibration": vibration_table,
    "molecular-rate-wahlstrand": wahlstrand_table,
    "molecular-rate-solprop": solprop_table,
}


def main() -> None:
    """Print the tables."""
    with warnings.catch_warnings():
        # The forms have no correction above the barrier-suppression field, and are
        # shown there as the codes evaluate them.
        warnings.simplefilter("ignore", ModelRangeWarning)
        for name, function in TABLES.items():
            print(f"<!-- {name}:start -->\n{function()}\n<!-- {name}:end -->\n")


if __name__ == "__main__":
    main()
