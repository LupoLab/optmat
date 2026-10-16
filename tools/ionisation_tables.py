# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The tables of the ionisation page (``docs/ionisation.md``), from the curated data.

    uv run python tools/ionisation_tables.py

prints the three tables; ``tests/test_ionisation_page.py`` checks that the page holds
exactly what this script prints, so the page is regenerated whenever the parameters or
their sources change:

- ``ionisation-parameters``: the parameters of each gas and their sources;
- ``ionisation-validation``: the static tunnelling rate of Tong, Zhao & Lin (2002),
  evaluated with the curated parameters, against the ab initio static-field rates;
- ``ionisation-solprop``: solprop's values against the curated ones.
"""

import json
from pathlib import Path

import numpy as np

from optmat.constants import AU_FIELD, AU_POLARIZABILITY, E_CHARGE
from optmat.datasets import datasets
from optmat.gas import _tunnelling as tunnelling  # pyright: ignore[reportPrivateUsage]
from optmat.gas.ionisation import IonisationParameters, ionisation_parameters

ROOT = Path(__file__).resolve().parents[1]
GASES = (
    "He",
    "Ne",
    "Ar",
    "Kr",
    "Xe",
    "H2",
    "D2",
    "N2",
    "O2",
    "CH4",
    "SF6",
    "N2O",
    "CO2",
)

# The ab initio references (tier 1, TC-439) and the fields at which the page compares:
# the lowest the reference gives, and fields up to about the barrier-suppression field.
REFERENCES = (
    ("He", "parker_2009/table1_present", "Parker et al. 2009, full TDSE"),
    ("He", "majety_2015/table1_hacc5_He", "Majety & Scrinzi 2015, haCC(5)"),
    ("H2", "saenz_2002/tableI_a_v0", "Saenz 2002, frozen nuclei, parallel"),
    ("H2", "majety_2015/table2_parallel_hacc5_H2", "Majety & Scrinzi 2015, parallel"),
    ("N2", "majety_2015/tableA1_N2_0deg", "Majety & Scrinzi 2015, parallel"),
    ("Ar", "zhang_2014/tableIX_gsz_Ar", "Zhang et al. 2014, SAE (GSZ potential)"),
)
FIELDS = (0.03, 0.045, 0.05, 0.053, 0.06, 0.08, 0.094, 0.1, 0.12, 0.14)


def _table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def _ev(energy: float | None) -> str:
    return "" if energy is None else f"{energy / E_CHARGE:.4f}"


def _au(alpha: float | None) -> str:
    return "" if alpha is None else f"{alpha / AU_POLARIZABILITY:.4g}"


def parameter_table() -> str:
    """Return the table of each gas's parameters.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for gas in GASES:
        p = ionisation_parameters(gas)
        structure = "" if p.default_structure is None else p.default_structure
        if p.default_structure is not None:
            c = p.structure[p.default_structure].coefficients
            structure += ": " + ", ".join(
                f"C{ell} = {v:g}" for ell, v in sorted(c.items())
            )
        rows.append(
            [
                gas,
                _ev(p.ionization_energy),
                p.ionization_energy_source,
                _ev(p.vertical_ionization_energy),
                f"{p.orbital} (m = {p.orbital_m})",
                f"{p.n_star:.4f}",
                structure,
                _au(p.polarizability),
                _au(p.ion_polarizability),
                ""
                if p.barrier_suppression_alpha is None
                else f"{p.barrier_suppression_alpha:g}",
                ""
                if p.barrier_suppression_zhang is None
                else ", ".join(f"{a:g}" for a in p.barrier_suppression_zhang),
                "" if p.effective_charge is None else f"{p.effective_charge:g}",
            ]
        )
    return _table(
        [
            "Gas",
            "I_p (eV)",
            "Source",
            "Vertical (eV)",
            "Orbital",
            "n*",
            "Structure coefficients",
            "alpha (au)",
            "alpha+ (au)",
            "alpha_BSI",
            "Zhang a1, a2, a3",
            "Z_eff (PPT, l = 0)",
        ],
        rows,
    )


def _rate(p: IonisationParameters, field: float) -> float:
    """Return the one-electron static rate for the orbital aligned with the field."""
    assert p.default_structure is not None
    s = p.structure[p.default_structure]
    factors = tunnelling.aligned_structure_factors(s.coefficients, s.m)
    return tunnelling.static_rate_au(field, p.ionization_energy_au, factors)


def validation_table() -> str:
    """Return the ratios of the tunnelling rate to the ab initio rates.

    Returns
    -------
    str
        A Markdown table: per reference, the ratio at each field it gives.
    """
    by_key = {d.key: d for d in datasets()}
    rows = []
    for gas, key, label in REFERENCES:
        dataset = by_key[key]
        fields = dataset.abscissa_values() / AU_FIELD
        rates = dataset.scaled_values()
        p = ionisation_parameters(gas)
        cells = []
        for field in FIELDS:
            match = np.flatnonzero(np.isclose(fields, field, rtol=0, atol=5e-4))
            if match.size == 0:
                cells.append("")
                continue
            i = int(match[0])
            cells.append(f"{_rate(p, float(fields[i])) / rates[i]:.2f}")
        rows.append([gas, label, *cells])
    return _table(["Gas", "Reference", *(f"F = {f:g}" for f in FIELDS)], rows)


def solprop_table() -> str:
    """Return solprop's values against the curated ones.

    Returns
    -------
    str
        A Markdown table.
    """
    path = ROOT / "tests" / "data" / "solprop_ionisation.json"
    fixture = json.loads(path.read_text())
    rows = []
    for gas in GASES:
        if gas not in fixture["gases"]:
            continue  # CO2 is optmat's, not solprop's
        theirs = fixture["gases"][gas]
        p = ionisation_parameters(gas)
        ip = theirs["potential_j"]
        d_ip = 1e3 * (ip - p.ionization_energy) / E_CHARGE
        alpha = theirs.get("polarizability_si")
        ion = theirs.get("ion_polarizability_si")
        rows.append(
            [
                gas,
                f"{ip / E_CHARGE:.4f}",
                f"{d_ip:+.1f}",
                "" if alpha is None else f"{alpha / p.polarizability - 1:+.1%}",
                ""
                if ion is None or p.ion_polarizability is None
                else f"{ion / p.ion_polarizability - 1:+.1%}",
                ""
                if "quantum_numbers" not in theirs
                else _quantum(theirs["quantum_numbers"]),
            ]
        )
    return _table(
        [
            "Gas",
            "solprop I_p (eV)",
            "- optmat (meV)",
            "alpha: solprop / optmat - 1",
            "alpha+: solprop / optmat - 1",
            "solprop (n, l, Z)",
        ],
        rows,
    )


def _quantum(numbers: dict[str, float]) -> str:
    return f"({numbers['principal']}, {numbers['orbital']}, {numbers['charge']:g})"


def main() -> None:
    """Print the three tables."""
    for name, table in (
        ("ionisation-parameters", parameter_table()),
        ("ionisation-validation", validation_table()),
        ("ionisation-solprop", solprop_table()),
    ):
        print(f"<!-- {name}:start -->\n{table}\n<!-- {name}:end -->\n")


if __name__ == "__main__":
    main()
