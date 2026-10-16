# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The tables of the Raman page (``docs/raman.md``), from the curated parameters.

    uv run python tools/raman_tables.py

prints the tables; ``tests/test_raman_page.py`` checks that the page holds exactly
what this script prints:

- ``raman-parameters``: each gas's constants, strengths and widths at 295 K and 1 atm;
- ``raman-gain``: the steady-state Stokes gain the parameters give, against the
  measured gains;
- ``raman-rotational``: the long-pulse rotational nonlinear index against the
  published values;
- ``raman-solprop``: solprop's values against the curated ones.
"""

import json
import math
from pathlib import Path

import numpy as np

from optmat.constants import AU_POLARIZABILITY, EPS0, HBAR, K_B, C
from optmat.datasets import datasets
from optmat.gas.linewidths import AMAGAT
from optmat.gas.raman import (
    AnisotropyTable,
    raman_lines,
    raman_parameters,
    rotational_nonlinear_index,
    stokes_gain,
)

ROOT = Path(__file__).resolve().parents[1]
GASES = ("N2", "O2", "H2", "D2", "CH4", "SF6", "CO2", "N2O")
ONE_ATM_295 = AMAGAT * 273.15 / 295.0
# The measured gains (TC-477, TC-478): the dataset, the gas, the line's J (None for an
# unresolved band), the pump wavelength (nm) and the density (amagat).
GAINS = (
    ("ottusch_1988/gain_H2", "H2", 1, 532.0, 60.0),
    ("bischel_1986b/gain_H2", "H2", 1, None, 20.0),
    ("ottusch_1988/gain_D2", "D2", 2, 532.0, 60.0),
    ("ottusch_1988/gain_CH4", "CH4", None, 532.0, 115.0),
)


def _table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def _range(values: np.ndarray, digits: str) -> str:
    low, high = float(np.min(values)), float(np.max(values))
    if f"{low:{digits}}" == f"{high:{digits}}":
        return f"{low:{digits}}"
    return f"{low:{digits}}-{high:{digits}}"


def _anisotropy_au(value: float | AnisotropyTable) -> str:
    if isinstance(value, AnisotropyTable):
        low, high = value.value[0], value.value[-1]
        return f"{low / AU_POLARIZABILITY:.3f}-{high / AU_POLARIZABILITY:.3f} (per J)"
    return f"{value / AU_POLARIZABILITY:.3f}"


def parameter_table() -> str:
    """Return each gas's constants, strengths and widths.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for gas in GASES:
        p = raman_parameters(gas)
        lines = raman_lines(gas, 295.0, ONE_ATM_295)
        branch = np.array(lines.branch)
        assert lines.fwhm is not None
        s = branch == "S"
        q = branch == "Q"
        rot = p.rotational
        band = p.vibrational[0]
        element = band.transition_polarizability
        element_cell = (
            _range(np.asarray(element.value) / AU_POLARIZABILITY, ".4f") + " (per J)"
            if isinstance(element, AnisotropyTable)
            else f"{element / AU_POLARIZABILITY:.4f}"
        )
        rows.append(
            [
                gas,
                "" if rot is None else f"{rot.constants.b:.5f}",
                "" if rot is None else f"{rot.constants.d:.3e}",
                "" if rot is None else _anisotropy_au(rot.anisotropy),
                "" if not s.any() else _range(lines.fwhm[s] / 1e9, ".2f"),
                f"{band.band_origin:.3f}",
                element_cell,
                _range(lines.fwhm[q] / 1e9, ".3g"),
            ]
        )
    return _table(
        [
            "Gas",
            "B (cm-1)",
            "D (cm-1)",
            "Delta alpha (au)",
            "S FWHM (GHz)",
            "nu_0 (cm-1)",
            "alpha_01 (au)",
            "Q FWHM (GHz)",
        ],
        rows,
    )


def _gain(gas: str, j: int | None, pump_nm: float, amagat: float) -> float:
    density = amagat * AMAGAT
    lines = raman_lines(gas, 295.0, density, wavelength=pump_nm * 1e-9)
    gain = stokes_gain(lines, pump_nm * 1e-9, density)
    select = np.array(lines.branch) == "Q"
    if j is not None:
        select &= lines.j == j
    return float(np.max(gain[select])) * 1e11  # m/W to cm/GW


def gain_table() -> str:
    """Return the model's steady-state gains against the measured ones.

    Returns
    -------
    str
        A Markdown table.
    """
    by_key = {d.key: d for d in datasets()}
    rows = []
    for key, gas, j, pump, amagat in GAINS:
        dataset = by_key[key]
        values = dataset.scaled_values() * (1e9 if dataset.unit == "cm/W" else 1.0)
        sigmas = dataset.standard_uncertainty() * (
            1e9 if dataset.unit == "cm/W" else 1.0
        )
        pumps = (
            [pump] * len(dataset)
            if pump is not None
            else [float(r) for r in dataset.rows]
        )
        for measured, sigma, nm in zip(values, sigmas, pumps, strict=True):
            model = _gain(gas, j, nm, amagat)
            rows.append(
                [
                    gas,
                    "Q" if j is None else f"Q({j})",
                    f"{nm:g}",
                    f"{amagat:g}",
                    key,
                    f"{measured:.2f} +- {sigma:.2f}",
                    f"{model:.2f}",
                    f"{model / measured - 1:+.0%}",
                ]
            )
    return _table(
        [
            "Gas",
            "Line",
            "Pump (nm)",
            "Density (amagat)",
            "Measurement",
            "Measured (cm/GW)",
            "Model (cm/GW)",
            "Model / measured - 1",
        ],
        rows,
    )


# The published rotational n2 (self coefficients at 1 atm, TC-493).
ROTATIONAL = (
    ("zahedpour_2015/n2_rot_N2", "N2"),
    ("zahedpour_2019/n2_rot_N2", "N2"),
    ("zahedpour_2015/n2_rot_O2", "O2"),
    ("zahedpour_2019/n2_rot_O2", "O2"),
    ("wahlstrand_2015/n2_rot_H2", "H2"),
    ("wahlstrand_2015/n2_rot_D2", "D2"),
)


def rotational_table() -> str:
    """Return the long-pulse rotational n2 against the published values.

    Returns
    -------
    str
        A Markdown table.
    """
    by_key = {d.key: d for d in datasets()}
    rows = []
    for key, gas in ROTATIONAL:
        dataset = by_key[key]
        assert dataset.state is not None
        temperature = dataset.state.temperature
        density = dataset.state.pressure / (K_B * temperature)  # ideal gas
        lines = raman_lines(gas, temperature, density)
        model = rotational_nonlinear_index(lines, density) * 1e4 / 1e-20
        published = float(dataset.scaled_values()[0]) / 1e-20
        rows.append(
            [
                gas,
                key,
                f"{temperature:g}",
                f"{published:g}",
                f"{model:.1f}",
                f"{model / published - 1:+.0%}",
            ]
        )
    return _table(
        [
            "Gas",
            "Published",
            "T (K)",
            "n2,rot published (1e-20 cm2/W)",
            "Model",
            "Model / published - 1",
        ],
        rows,
    )


def solprop_table() -> str:
    """Return solprop's values against the curated ones.

    Returns
    -------
    str
        A Markdown table.
    """
    fixture = json.loads(
        (ROOT / "tests" / "data" / "solprop_raman.json").read_text(encoding="utf-8")
    )
    rows = []
    for gas in GASES:
        theirs = fixture["raman"].get(gas, {})
        if theirs.get("status") != "available":
            continue  # CO2 is optmat's, not solprop's; solprop has no N2O lines
        p = raman_parameters(gas)
        cells = [gas]
        rot, band = theirs["rotational"], theirs["vibrational"]
        if rot is not None and p.rotational is not None:
            ours = p.rotational
            anisotropy = ours.anisotropy
            value = (
                float(anisotropy(np.array([1]))[0])
                if isinstance(anisotropy, AnisotropyTable)
                else anisotropy
            )
            cells += [
                f"{rot['b'] / 100:.4f} / {ours.constants.b:.4f}",
                f"{rot['d'] / 100:.3e} / {ours.constants.d:.3e}",
                f"{rot['polarisability_anisotropy'] * 4 * math.pi * EPS0 / value - 1:+.1%}",  # noqa: E501 - one table cell
            ]
        else:
            cells += ["", "", ""]
        ours_band = p.vibrational[0]
        k_theirs = (4 * math.pi * EPS0 * band["polarisability_derivative"]) ** 2 / (
            2 * band["reduced_mass"] * band["frequency"] * EPS0
        )
        element = ours_band.transition_polarizability
        if isinstance(element, AnisotropyTable):
            element = float(element(np.array([1 if gas == "H2" else 2]))[0])
        k_ours = element**2 / (HBAR * EPS0)
        cells += [
            f"{band['frequency'] / (2 * math.pi * C * 100):.2f}"
            f" / {ours_band.band_origin:.2f}",
            f"{k_theirs / k_ours - 1:+.1%}",
        ]
        rows.append(cells)
    return _table(
        [
            "Gas",
            "B: solprop / optmat (cm-1)",
            "D: solprop / optmat (cm-1)",
            "Delta alpha: solprop / optmat - 1",
            "nu: solprop / optmat (cm-1)",
            "Vibrational A: solprop / optmat - 1",
        ],
        rows,
    )


def main() -> None:
    """Print the tables."""
    for name, table in (
        ("raman-parameters", parameter_table()),
        ("raman-gain", gain_table()),
        ("raman-rotational", rotational_table()),
        ("raman-solprop", solprop_table()),
    ):
        print(f"<!-- {name}:start -->\n{table}\n<!-- {name}:end -->\n")


if __name__ == "__main__":
    main()
