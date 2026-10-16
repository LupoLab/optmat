# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The curated two-photon data of the noble gases on one footing: the manual's tables.

Writes the marker tables of the manual page *Two-photon ionisation and absorption*
(``docs/two_photon.md``); ``tests/test_two_photon_page.py`` checks that the page holds
this script's output. Nothing is fitted: the tables state what the curated records
hold, compared with one another (as reconciled in the register entries of the
curation).

* **Thresholds:** the two-photon ionisation thresholds of each gas, half the neutral's
  ionisation energy (the ion's ²P₃/₂ ground level) and half of it plus the ion's ²P₁/₂
  level, from ``nist_asd_2026``; the first one-photon resonance E₁ and half the lowest
  two-photon level, E₂/2, from the same record (OEP 5's selections).
* **Coverage:** every two-photon ionisation dataset of each gas: its source, form
  (vector, raster, table, text), polarisation, gauge and the range of photon energies
  it covers inside the window from I_p/2 to E₁.
* **Comparison:** the total cross section of each independent calculation or
  measurement at a few photon energies off resonance, on one footing (cm⁴ s,
  Gamma = sigma Phi^2, photon energy), and its ratio to a named comparator. Partial
  cross sections drawn for the final symmetries are summed; values are interpolated
  linearly in log sigma between a curve's own vertices and never extrapolated.

McCown et al.'s coefficient alpha (cm⁴ W⁻¹, rate alpha I²/ħω; their Eq. 9) is converted
to sigma = alpha ħω at the printed wavelength.

Usage::

    uv run python tools/two_photon_tables.py [--update-docs]
"""

import argparse
import functools
import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from optmat.constants import E_CHARGE, HBAR, C
from optmat.datasets import Dataset, sources

ROOT = Path(__file__).resolve().parents[1]
PAGES = (ROOT / "docs/two_photon.md", ROOT / "docs/user-guide/two-photon-data.md")
GASES = ("He", "Ne", "Ar", "Kr", "Xe")
HC_EV_CM = 2 * math.pi * HBAR * C / E_CHARGE * 100
"""hc in eV cm, from the package's CODATA constants."""

E1 = {
    "He": "level_He_1s2p_1P1",
    "Ne": "level_Ne_odd_J1",
    "Ar": "level_Ar_odd_J1",
    "Kr": "level_Kr_odd_J1",
    "Xe": "level_Xe_odd_J1",
}
E2 = {
    "He": "level_He_1s2s_1S0",
    "Ne": "level_Ne_even_J2",
    "Ar": "level_Ar_even_J2",
    "Kr": "level_Kr_even_J2",
    "Xe": "level_Xe_even_J2",
}
"""The levels OEP 5's support records select (``data/gas/nonlinear_support.toml``)."""


@dataclass(frozen=True)
class Curve:
    """One calculation or measurement: its datasets, summed if partial."""

    label: str
    keys: tuple[str, ...]
    note: str = ""


COMPARISONS: dict[str, tuple[tuple[float, ...], tuple[Curve, ...]]] = {
    "He": (
        (13.0, 14.0, 16.0, 18.0, 20.0),
        (
            Curve(
                "Saenz & Lambropoulos 1999, ¹S + ¹D (velocity)",
                ("saenz_1999/fig1a_sigma2_1S", "saenz_1999/fig1b_sigma2_1D"),
            ),
            Curve(
                "Mihelič & Horvat 2021, L = 0 + 2",
                ("mihelic_2021/fig3a_sigma2_L0", "mihelic_2021/fig3a_sigma2_L2"),
            ),
            Curve(
                "Nikolopoulos & Lambropoulos 2001, S + D (Fig. 4)",
                ("nikolopoulos_2001/fig4_sigma2_S", "nikolopoulos_2001/fig4_sigma2_D"),
            ),
            Curve(
                "van der Hart & Bingham 2005, R-matrix Floquet (1e13 W/cm²)",
                ("vanderhart_2005/fig2_sigma2_rmf",),
            ),
        ),
    ),
    "Ne": (
        (11.0, 12.0, 14.0, 16.0),
        (
            Curve("McKenna & van der Hart 2004", ("mckenna_2004/fig1_sigma2_Ne",)),
            Curve(
                "Moccia, Rahman & Rizzo 1983, velocity (raster)",
                ("moccia_1983/fig1_sigma2_V_Ne",),
            ),
            Curve(
                "Moccia, Rahman & Rizzo 1983, length (raster)",
                ("moccia_1983/fig1_sigma2_L_Ne",),
            ),
            Curve(
                "McGuire 1981, circular (raster)",
                ("mcguire_1981/fig1_sigma2_Ne",),
                "circular",
            ),
        ),
    ),
    "Ar": (
        (8.0, 9.0, 10.0, 11.0),
        (
            Curve("McKenna & van der Hart 2004", ("mckenna_2004/fig5_sigma2_Ar",)),
            Curve(
                "Petrov et al. 2016, corr + CP, length", ("petrov_2016/fig10_corrcp_L",)
            ),
            Curve(
                "Petrov et al. 2016, corr + CP, velocity",
                ("petrov_2016/fig10_corrcp_V",),
            ),
            Curve("Petrov et al. 2016, corr, length", ("petrov_2016/fig10_corr_L",)),
            Curve(
                "Starace & Jiang 1987, IIC, length (raster)",
                ("starace_1987/fig7_iic_linear_Ar",),
            ),
            Curve(
                "Starace & Jiang 1987, ESI, length (raster)",
                ("starace_1987/fig7_esi_linear_Ar",),
            ),
            Curve(
                "Pan, Gao & Starace 1990, length (raster)",
                ("pan_1990b/fig4c_all_L_Ar",),
            ),
            Curve(
                "Pan, Gao & Starace 1990, velocity (raster)",
                ("pan_1990b/fig4c_all_V_Ar",),
            ),
            Curve(
                "Pindzola & Kelly 1975, second order, length (raster)",
                (
                    "pindzola_1975/fig2_L_p1D_Ar",
                    "pindzola_1975/fig3_L_p1S_Ar",
                    "pindzola_1975/fig4_L_f1D_Ar",
                ),
            ),
            Curve(
                "Moccia, Rahman & Rizzo 1983, velocity (raster)",
                ("moccia_1983/fig2_sigma2_V_Ar",),
            ),
            Curve(
                "Moccia, Rahman & Rizzo 1983, length (raster)",
                ("moccia_1983/fig2_sigma2_L_Ar",),
            ),
            Curve("Fink & Johnson 1990, linear (raster)", ("fink_1990/fig3_linear",)),
            Curve(
                "McGuire 1981, circular (raster)",
                ("mcguire_1981/fig2_sigma2_Ar",),
                "circular",
            ),
        ),
    ),
    "Kr": (
        (7.2, 7.5, 8.0, 8.5, 9.5),
        (
            Curve(
                "L'Huillier, Tang & Lambropoulos 1989, MQDT (raster)",
                ("lhuillier_1989/fig6_total_Kr",),
            ),
            Curve(
                "L'Huillier, Tang & Lambropoulos 1989, above ²P₁/₂ (raster)",
                ("lhuillier_1989/fig10_total_Kr",),
            ),
            Curve("Fink & Johnson 1990, linear (raster)", ("fink_1990/fig4_linear",)),
            Curve(
                "Fink & Johnson 1990, circular (raster)",
                ("fink_1990/fig4_circular",),
                "circular",
            ),
            Curve(
                "McGuire 1981, circular (raster)",
                ("mcguire_1981/fig3_sigma2_Kr",),
                "circular",
            ),
        ),
    ),
    "Xe": (
        (6.42, 7.0, 7.5, 8.0),
        (
            Curve(
                "Lagutin et al. 2017, CIPFCP, linear",
                ("lagutin_2017/fig8_cipfcp_linear",),
            ),
            Curve(
                "Lagutin et al. 2017, CIPFCP, circular",
                ("lagutin_2017/fig8_cipfcp_circular",),
                "circular",
            ),
            Curve(
                "Karamatskou & Santra 2017, TDCIS (Fig. 1)",
                ("karamatskou_2017/fig1_interchannel",),
            ),
            Curve("Fink & Johnson 1990, linear (raster)", ("fink_1990/fig5_linear",)),
            Curve(
                "Gangopadhyay et al. 1986, MQDT, Table I (sum of cores)",
                ("gangopadhyay_1986/tab1_jc32", "gangopadhyay_1986/tab1_jc12"),
            ),
            Curve(
                "Gangopadhyay et al. 1986, MQDT, between the thresholds (raster)",
                ("gangopadhyay_1986/fig1",),
            ),
            Curve(
                "L'Huillier & Wendin 1987 (J. Phys. B), RPAE (raster)",
                ("lhuillier_1987b/fig4_rpae",),
            ),
            Curve(
                "L'Huillier & Wendin 1987, RPAE (raster)",
                ("lhuillier_1987/fig6_rpae_Xe",),
            ),
            Curve(
                "L'Huillier, Tang & Lambropoulos 1989, MQDT (raster)",
                ("lhuillier_1989/fig6_total_Xe",),
            ),
            Curve(
                "L'Huillier, Tang & Lambropoulos 1989, above ²P₁/₂ (raster)",
                ("lhuillier_1989/fig10_total_Xe",),
            ),
            Curve("McCown, Ediger & Eden 1982, measured", ("mccown_1982/alpha_Xe",)),
            Curve("Bischel et al. 1979, estimate", ("bischel_1979/alpha_Xe",)),
            Curve("Hodges, Lee & Moseley 1981, lower limit", ("hodges_1981/alpha_Xe",)),
            Curve("Santra & Greene 2004", ("santra_2004/sigma2_Xe",)),
            Curve(
                "McGuire 1981, circular (raster)",
                ("mcguire_1981/fig4_sigma2_Xe",),
                "circular",
            ),
        ),
    ),
}
"""Per gas: the photon energies (eV) and the curves compared; the first is the
comparator of the ratios."""
EDGE = 1e-9
"""Relative tolerance at a curve's ends, for the round trip through the wavelength."""
POINT_WINDOW = 0.01
"""A one-point dataset (a measurement at one wavelength) is compared at a tabulated
energy within this, eV."""


@functools.cache
def _all() -> dict[str, Dataset]:
    return {d.key: d for s in sources().values() for d in s.datasets}


def dataset(key: str) -> Dataset:
    """Return the curated dataset ``source/id``."""
    return _all()[key]


def _point(key: str) -> float:
    """Return a one-value dataset's value as printed (after scale), as a float."""
    return float(dataset(key).scaled_values()[0])


def level_ev(key: str) -> float:
    """Return a ``nist_asd_2026`` level or ionisation energy in eV."""
    d = dataset(f"nist_asd_2026/{key}")
    value = _point(d.key)
    return value if d.unit == "eV" else value * HC_EV_CM


def thresholds(gas: str) -> dict[str, float]:
    """Return the two-photon thresholds, E₁ and E₂/2 of a gas, in eV."""
    ionisation = level_ev(f"ie_{gas}_I")
    split = 0.0 if gas == "He" else level_ev(f"level_{gas}+_2P12")
    return {
        "lower": ionisation / 2,
        "upper": (ionisation + split) / 2,
        "e1": level_ev(E1[gas]),
        "e2_half": level_ev(E2[gas]) / 2,
    }


def photon_ev(d: Dataset) -> np.ndarray:
    """Return a dataset's abscissa as the photon energy, eV."""
    return 2 * math.pi * HBAR * C / E_CHARGE / d.vacuum_wavelength()


def cm4_s(d: Dataset) -> np.ndarray:
    """Return a two-photon dataset's values in cm⁴ s (Gamma = sigma Phi^2)."""
    values = d.scaled_values()
    if d.unit == "cm4 s":
        return values
    if d.unit == "cm4 W-1":
        # McCown et al.'s Eq. 9: dn/dt = alpha [Xe] I^2 / hbar omega, so sigma =
        # alpha hbar omega with hbar omega in J (TC-884).
        return values * photon_ev(d) * E_CHARGE
    msg = f"{d.key}: no conversion from {d.unit!r} to cm4 s"
    raise ValueError(msg)


def at(key: str, energy: float) -> float | None:
    """Return a dataset's value at a photon energy (eV), or None outside it."""
    d = dataset(key)
    x, y = photon_ev(d), cm4_s(d)
    if len(x) == 1:
        return float(y[0]) if abs(float(x[0]) - energy) <= POINT_WINDOW else None
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    # A table row printed at the energy itself is in range: the conversion through
    # the vacuum wavelength moves it by a rounding error.
    if not x[0] * (1 - EDGE) <= energy <= x[-1] * (1 + EDGE):
        return None
    return float(np.exp(np.interp(energy, x, np.log(y))))


def curve_at(curve: Curve, energy: float) -> float | None:
    """Return a curve's total at a photon energy: the sum of its partials."""
    parts = [at(key, energy) for key in curve.keys]
    if any(p is None for p in parts):
        return None
    return float(sum(p for p in parts if p is not None))


def _form(d: Dataset) -> str:
    if d.digitised:
        return d.digitised
    return "text" if d.table == "text" else "table"


def coverage(gas: str) -> list[dict[str, str]]:
    """Return one row per two-photon dataset of the gas: what covers the window."""
    edges = thresholds(gas)
    rows = []
    for d in _all().values():
        if not (
            d.material == gas
            and d.quantity == "multiphoton_ionization_cross_section"
            and str(d.parameters.get("photons")) == "2"
        ):
            continue
        x = photon_ev(d)
        inside = x[(x >= edges["lower"]) & (x < edges["e1"])]
        span = f"{inside.min():.3f} to {inside.max():.3f}" if len(inside) > 1 else ""
        if len(inside) == 1:
            span = f"{inside[0]:.3f}"
        rows.append(
            {
                "dataset": f"`{d.key}`",
                "form": _form(d),
                "use": d.use,
                "polarisation": str(d.parameters.get("polarization", "")),
                "channel": str(d.parameters.get("channel", "")) or "total",
                "window": span or "outside",
            }
        )
    return rows


def _number(value: float | None) -> str:
    return "—" if value is None else f"{value:.3g}"


def _table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def threshold_table() -> str:
    """Return the thresholds table."""
    nm = 1239.841984  # eV nm, only to print wavelengths (hc to 10 digits)
    rows = []
    for gas in GASES:
        t = thresholds(gas)
        upper = "—" if gas == "He" else f"{t['upper']:.4f} ({nm / t['upper']:.2f})"
        rows.append(
            [
                gas,
                f"{t['lower']:.4f} ({nm / t['lower']:.2f})",
                upper,
                f"{t['e1']:.4f} ({nm / t['e1']:.2f})",
                f"{t['e2_half']:.4f} ({nm / t['e2_half']:.2f})",
            ]
        )
    return _table(
        [
            "Gas",
            "I_p/2, ²P₃/₂ (eV (nm))",
            "²P₁/₂ threshold (eV (nm))",
            "E₁ (eV (nm))",
            "E₂/2 (eV (nm))",
        ],
        rows,
    )


def coverage_table(gas: str) -> str:
    """Return the coverage table of one gas."""
    rows = [
        [
            r["dataset"],
            r["form"],
            r["use"],
            r["polarisation"],
            r["channel"],
            r["window"],
        ]
        for r in sorted(coverage(gas), key=lambda r: r["dataset"])
    ]
    return _table(
        [
            "Dataset",
            "Form",
            "Use",
            "Polarisation",
            "Channel",
            "Photon energies in the window (eV)",
        ],
        rows,
    )


def comparison_table(gas: str) -> str:
    """Return the comparison table of one gas: values and ratios to the first curve."""
    energies, curves = COMPARISONS[gas]
    reference = curves[0]
    rows = []
    for curve in curves:
        cells = [curve.label]
        for energy in energies:
            value = curve_at(curve, energy)
            base = curve_at(reference, energy)
            if value is None or curve is reference or base is None:
                cells.append(_number(value))
            else:
                cells.append(f"{value:.3g} ({value / base:.2f})")
        rows.append(cells)
    header = [r"$\sigma^{(2)}$ (cm⁴ s) at", *(f"{e:g} eV" for e in energies)]
    return _table(header, rows)


def tables() -> dict[str, str]:
    """Return every marker table of the page, by marker name."""
    result = {"two-photon-thresholds": threshold_table()}
    for gas in GASES:
        result[f"two-photon-coverage-{gas}"] = coverage_table(gas)
        result[f"two-photon-comparison-{gas}"] = comparison_table(gas)
    return result


def main(argv: Sequence[str] | None = None) -> None:
    """Print the tables, or rewrite them in the manual page."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    parser.add_argument("--update-docs", action="store_true")
    args = parser.parse_args(argv)
    result = tables()
    if not args.update_docs:
        for name, table in result.items():
            print(f"<!-- {name}:start -->\n{table}\n<!-- {name}:end -->\n")
        return
    pages = {path: path.read_text(encoding="utf-8") for path in PAGES}
    for name, table in result.items():
        start, end = f"<!-- {name}:start -->", f"<!-- {name}:end -->"
        matches = [path for path, text in pages.items() if start in text]
        if len(matches) != 1:
            raise ValueError(f"Expected one destination for {name}, got {matches}")
        path = matches[0]
        head, rest = pages[path].split(start, 1)
        _, tail = rest.split(end, 1)
        pages[path] = f"{head}{start}\n{table}\n{end}{tail}"
    for path, text in pages.items():
        path.write_text(text, encoding="utf-8")
        print(f"updated {path}")


if __name__ == "__main__":
    main()
