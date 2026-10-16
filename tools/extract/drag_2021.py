# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Drag, Marmuse and Blondel (2021): Table 1's compiled two-photon data.

The library holds the accepted manuscript (HAL hal-03449890; IOP's 'Accepted
Manuscript' cover, PDF p. 2, then the manuscript's pages 1-33 as PDF pp. 3-35). The
paper's own results, the integrated cross sections of the Xe 6p'[3/2]_2 and
6p'[1/2]_0 levels, are printed only in the text (abstract, p. 1, and Sec. 4.1, p. 25):
they are inline points of the source record. This script extracts Table 1 (p. 5,
``tbl_e87cefb4180f``, merged cells, read from its JSON cells; image checked), 'the state
of our knowledge' of other authors' integrated two-photon excitation cross sections of
Kr and Xe, whose absolute values (m^4) are reproduced datasets: Khambatta et al. (1988)
for Kr 6p[3/2]_2, Dakka et al. (2018) for Kr 5p[3/2]_2 (dakka_2018), Kroll and Bischel
(1990) for Xe 6p[1/2]_0, 6p[3/2]_2 and 6p[5/2]_2 (kroll_1990). The ratios to H, N and
O (Niemi et al. 2001, 2005; Goehlich et al. 1998; Elliott et al. 2016) have no quantity
and are not extracted. A cell reads ``5.4 \times 10^{-45} \pm 50\%``: the mantissa is
the value, the power of ten the dataset's scale (the cells of one dataset share it,
checked), the percentage its relative uncertainty.

:func:`checks` tests the reproductions against the curated originals (the conversion
10^-35 cm^4 = 10^-43 m^4, and the uncertainties), and the paper's own numbers: the
asymmetric errors against Table 3's uncertainty budget, the ratio of the two results,
the wavelengths against the level energies the text gives, and Table 2's derived
quantities (the Doppler peak of Eq. 32 normalised on the 2 omega scale, the resonant
and effective cross sections, the Rabi frequency of Eq. 11).
"""

import math
import re

from _json_table import json_table
from _library import Point, table, write

from optmat.constants import ATOMIC_MASS, HBAR, K_B, C
from optmat.datasets import DATA_DIRECTORY, Dataset, load_source
from optmat.gas.eos import equation_of_state

DIRECTORY = "Drag_2021_MeasurementTwoPhoton"
TABLE_1 = "tbl_e87cefb4180f"
TABLE_2 = "tbl_99ff3ba44ede"
TABLE_3 = "tbl_35c5917c5d18"
REPRODUCED = {
    "Khambatta et al. (1988)": "Kr/drag_2021__khambatta_1988_Kr.csv",
    "Dakka et al. (2018)": "Kr/drag_2021__dakka_2018_Kr.csv",
    "Kröll & Bischel (1990)": "Xe/drag_2021__kroll_1990_Xe.csv",
}
"""Table 1's absolute cross sections: reference -> data file. The Kroll and Bischel rows
are three, the reference printed once on the middle row (merged cells)."""
_CELL = re.compile(r"^(\d\.\d) \\times 10\^\{(-\d+)\} \\pm (\d+)\\%$")
_LEVEL = re.compile(r"^(\d[spdf]'?\[\d/\d\])_\{(\d)\}$")

RESULTS = {"6p'[3/2]2": "sigma_6p1_32_2", "6p'[1/2]0": "sigma_6p1_12_0"}
"""The record's inline datasets of the paper's results (p. 25 and the abstract)."""
LEVELS_PER_M = {"6p'[3/2]2": 8916235.6, "6p'[1/2]0": 8986001.5}
"""The levels' energies the text gives (pp. 3, 5; Kramida et al. 2020), m^-1."""
STATISTICAL = {"6p'[3/2]2": 0.07, "6p'[1/2]0": 0.12}
"""Table 3's statistical dispersion, '+-7%/ +-12%', for the two levels respectively."""
PRINTED_RELATIVE = {"6p'[3/2]2": (0.34, 0.26), "6p'[1/2]0": (0.40, 0.29)}
"""'1.00^{+34}_{-26} or 1.00^{+40}_{-29}' (p. 26): the relative uncertainties."""
RATIO = 1.38
"""p. 26: the ratio of the 6p'[1/2]_0 to the 6p'[3/2]_2 cross section (+21/-18 %)."""
TABLE_2_TEMPERATURE_K = 300.0
"""Table 2's caption: 'the normalized Doppler broadening factor (at 300 K)'."""
TABLE_2_TOLERANCE = 0.05
"""Table 2's derived values against the printed: two printed figures on each side."""


def reproduced() -> dict[str, list[list[str]]]:
    """Return Table 1's rows of each absolute-value reference."""
    grid = json_table(DIRECTORY, TABLE_1)
    out: dict[str, list[list[str]]] = {name: [] for name in REPRODUCED}
    # The Kroll and Bischel rows: the reference on the middle one, blank around it.
    kroll = "Kröll & Bischel (1990)"
    for index, cells in enumerate(grid[1:], start=1):
        name = cells[5]
        if not name and kroll in (grid[index - 1][5], grid[index + 1][5]):
            name = kroll
        if name in out:
            out[name].append(cells)
    return out


def _parse(cells: list[str]) -> tuple[str, str, str, int, int]:
    """Return a row's level, wavelength, mantissa, power of ten and percentage."""
    level = _LEVEL.match(cells[1])
    value = _CELL.match(cells[3])
    if not level or not value:
        msg = f"Table 1: unexpected row {cells}"
        raise ValueError(msg)
    label = f"{level.group(1)}{level.group(2)}"
    return label, cells[2], value.group(1), int(value.group(2)), int(value.group(3))


def _original(source: str, dataset: str) -> dict[str, float]:
    """Return a curated dataset's values in SI (m^4), by row."""
    record = load_source(DATA_DIRECTORY / "sources" / f"{source}.toml")
    data = next(d for d in record.datasets if d.id == dataset)
    factor = {"cm4": 1e-8, "m4": 1.0}[data.unit]
    return dict(zip(data.rows, data.scaled_values() * factor, strict=True))


def _check_reproduced(out: dict[str, object], failures: list[str]) -> None:
    """Check Table 1's absolute values against the curated originals and the record."""
    rows = reproduced()
    paschen = {"6p[1/2]0": "2p5", "6p[3/2]2": "2p6", "6p[5/2]2": "2p9"}
    kroll = _original("kroll_1990", "sigma_tot")
    kroll_sigma = {"2p5": 0.30, "2p6": 0.30, "2p9": 0.50}
    for cells in rows["Kröll & Bischel (1990)"]:
        label, _, mantissa, exponent, percent = _parse(cells)
        quoted = float(mantissa) * 10**exponent
        original = kroll[paschen[label]]
        out[f"Table 1 {label} / kroll_1990 sigma_tot"] = quoted / original
        if abs(quoted / original - 1) > 1e-9:
            failures.append(f"{label}: {quoted} m4 is not kroll_1990's {original}")
        if percent / 100 != kroll_sigma[paschen[label]]:
            failures.append(f"{label}: +-{percent} % is not kroll_1990's")
    dakka = _original("dakka_2018", "sigma0_2p6")["2p6"]
    _, _, mantissa, exponent, percent = _parse(rows["Dakka et al. (2018)"][0])
    quoted = float(mantissa) * 10**exponent
    out["Table 1 Kr 5p[3/2]2 / dakka_2018 sigma0"] = quoted / dakka
    if abs(quoted / dakka - 1) > 1e-9:
        failures.append(f"Dakka: {quoted} m4 is not dakka_2018's {dakka}")
    out["Table 1 Dakka +-% against 2.1/5.2"] = (percent, 100 * 2.1 / 5.2)
    # The record's scale and relative uncertainty of each reproduced row are Table 1's.
    datasets = _record()
    for name, file in REPRODUCED.items():
        data = datasets[file.removesuffix(".csv").split("__")[1]]
        for cells in rows[name]:
            label, _, _, exponent, percent = _parse(cells)
            spec = data.uncertainty
            relative = spec.relative_rows.get(label, spec.relative)
            if data.scale != 10.0**exponent or relative != percent / 100:
                failures.append(f"{name} {label}: the record's scale or +-% differs")


def _record() -> dict[str, Dataset]:
    """Return this source's datasets, by id."""
    record = load_source(DATA_DIRECTORY / "sources" / "drag_2021.toml")
    return {d.id: d for d in record.datasets}


def _check_results(out: dict[str, object], failures: list[str]) -> None:
    """Check the record's results: errors against Table 3, ratio, wavelengths."""
    datasets = _record()
    results = {
        level: (
            datasets[dataset].abscissae[0],
            float(datasets[dataset].values[0]),
            float(str(datasets[dataset].parameters["error_plus"])),
            float(str(datasets[dataset].parameters["error_minus"])),
        )
        for level, dataset in RESULTS.items()
    }
    budget = [cells[1] for cells in table(DIRECTORY, TABLE_3)[2:]]
    common = math.prod(1 + float(cell.strip("±%")) / 100 for cell in budget)
    for level, (_, value, plus, minus) in results.items():
        factor = (1 + STATISTICAL[level]) * common
        out[f"{level}: product of (1 + u), 1 - 1/product"] = (factor, 1 - 1 / factor)
        printed_plus, printed_minus = PRINTED_RELATIVE[level]
        if abs(factor - 1 - printed_plus) > 0.01:
            failures.append(f"{level}: budget gives +{factor - 1:.3f}")
        if abs(plus / value - printed_plus) > 0.005:
            failures.append(f"{level}: +{plus} is not {printed_plus} of {value}")
        out[f"{level}: -error / value, printed relative"] = (
            minus / value,
            printed_minus,
        )
        if abs(minus / value - (1 - 1 / factor)) > 0.01:
            failures.append(f"{level}: -{minus} is not (1 - 1/{factor:.3f}) {value}")
    ratio = results["6p'[1/2]0"][1] / results["6p'[3/2]2"][1]
    out["ratio 6p'[1/2]0 / 6p'[3/2]2"] = (ratio, RATIO)
    if f"{ratio:.2f}" != f"{RATIO:.2f}":
        failures.append(f"ratio {ratio:.3f}, printed {RATIO}")
    # The wavelengths are half the levels' energies, in vacuum.
    for level, energy in LEVELS_PER_M.items():
        wavelength = 2 / energy * 1e9
        out[f"{level}: 2 / E (nm, vacuum)"] = wavelength
        if f"{wavelength:.1f}" != results[level][0]:
            failures.append(f"{level}: 2/E = {wavelength:.3f} nm")


def _table_2() -> dict[str, float]:
    """Return Table 2's printed values, by symbol."""
    out = {}
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        cell = (
            cells[2].replace("10^{", "1 \\times 10^{")
            if cells[2][:3] == "10^"
            else cells[2]
        )
        mantissa, _, power = cell.partition(" \\times 10^{")
        scale = 10 ** int(power.rstrip("}")) if power else 1
        out[cells[1]] = float(mantissa) * scale
    return out


def _check_table_2(out: dict[str, object], failures: list[str]) -> None:
    """Check Table 2 (6p'[1/2]0) against its own definitions.

    D(0) = c / (2 omega_0) sqrt(M / (2 pi k T)) (Eq. 32), the Doppler distribution
    normalised on the 2 omega scale; sigma_D(0) = sigma D(0); sigma(0) = sigma / (pi
    beta) (Eq. 10); Omega = sqrt(2 sigma / pi) Phi_0 (Eq. 11).
    """
    printed = _table_2()
    omega = 2 * math.pi * C / (printed["\\lambda"] * 1e-9)
    mass = equation_of_state("Xe").molar_mass / 1e-3 * ATOMIC_MASS
    thermal = math.sqrt(mass / (2 * math.pi * K_B * TABLE_2_TEMPERATURE_K))
    sigma = printed["\\bar{\\sigma}^{(2)}"]
    derived = {
        "D(0)": (C / (2 * omega) * thermal, "D(0)"),
        "sigma_D(0) = sigma D(0)": (
            sigma * printed["D(0)"],
            "\\sigma_{\\text{D}}^{(2)}(0)",
        ),
        "sigma(0) = sigma / (pi beta)": (
            sigma / (math.pi * printed["\\beta \\geq \\Gamma/2"]),
            "\\sigma^{(2)}(0)",
        ),
        "Phi_0 = I / (hbar omega)": (printed["I"] / (HBAR * omega), "\\Phi_{0}"),
        "Omega = sqrt(2 sigma / pi) Phi_0": (
            math.sqrt(2 * sigma / math.pi) * printed["\\Phi_{0}"],
            "\\Omega",
        ),
    }
    for name, (computed, symbol) in derived.items():
        out[f"Table 2 {name}"] = (computed, printed[symbol])
        if abs(computed / printed[symbol] - 1) > TABLE_2_TOLERANCE:
            failures.append(f"Table 2 {name}: {computed:.3e} vs {printed[symbol]}")


def checks() -> dict[str, object]:
    """Return the checks of Table 1's reproductions and of the paper's own numbers.

    Raises if one fails.
    """
    out: dict[str, object] = {}
    failures: list[str] = []
    _check_reproduced(out, failures)
    _check_results(out, failures)
    _check_table_2(out, failures)
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return Table 1's absolute cross sections, by reference."""
    files: dict[str, list[Point]] = {}
    for name, rows in reproduced().items():
        points = []
        exponents = set()
        for cells in rows:
            label, wavelength, mantissa, exponent, _ = _parse(cells)
            exponents.add(exponent)
            points.append(Point(label, wavelength, mantissa, "", " | ".join(cells)))
        if len(exponents) != 1:
            msg = f"Table 1, {name}: the powers of ten differ, {exponents}"
            raise ValueError(msg)
        files[REPRODUCED[name]] = points
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()
