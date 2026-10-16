# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Pindzola and Kelly (1975): two-photon ionisation cross sections of Ar.

The paper prints its partial cross sections of two-photon ionisation of the Ar 3p shell
(many-body perturbation theory, Hartree-Fock V^(N-1) orbitals, length and velocity
forms, linear polarisation) only in figures; its one table (Table I) is a sum-rule test
of radial matrix elements, not a cross section. The PDF is a scan: every page is one
embedded 1-bit image, 4938 x 6686 pixels (about 580 dpi). The figures are read from
those native pixels (data rule 8, TC-699). Journal page = PDF page + 1542.

* **Fig. 2** (p. 6/1548, ``fig_4c9504d1d62e``): the kp 1D partial cross section,
  sigma/F in 10^-18 cm^2 per atomic unit of flux, 0 to 10: length (L, solid) and
  velocity (V, dashed), second order.
* **Fig. 3** (p. 6/1548, ``fig_7688b2052b30``): kp 1S, the same.
* **Fig. 4** (p. 6/1548, ``fig_17415fef50fa``): kf 1D, 0 to 20+: second-order length
  (L) and velocity (V), both solid. The ground-state-correlated LG and VG (dashed, run
  within a few pixels of each other and cross twice) and the approximations AHF and AHY
  are not read.
* **Fig. 6** (p. 7/1549, ``fig_e8afb9fc4099``): kf 1D with ground-state and
  intermediate-state correlations, LGI (solid) and VGI (dashed); Fig. 4's L (solid) and
  V (dashed) drawn again, read for the checks only.

The partial cross sections are 'calculated ... at intervals of 0.5 eV starting at omega
= 8.0 eV' (p. 5/1547): the drawn curves join those points, so the values are read at
8.0, 8.5, ..., 14.0 eV only.

The reading uses the machinery of :mod:`pan_1990b` (the same kind of framed plot on a
scan of the same resolution): the tick marks inside each side of the frame, found in
the image and checked to lie on one linear scale per side (:class:`pan_1990b.Figure`,
:class:`pan_1990b.Calibration`); each curve traced from a seed pixel and its style
checked (:func:`pan_1990b.trace`, :func:`pan_1990b.style`); a value where the line is
clear of every other dark pixel (:func:`pan_1990b._read`). The samples left out are
listed in :data:`LEFT_OUT`, those read on two curves within a line width of each other
in :data:`TOGETHER`; the script fails if either changes. The ``printed`` cell of each
value gives its pixel (column, row).

**Checks** (:func:`checks`, run by :func:`main` and :func:`extract`; it fails on a
miss): the vertical lines A-F drawn at the thresholds and resonances against the
energies the text prints ('omega = |epsilon_3p| / 2 = 7.88 eV', 11.62, 14.09, 14.15 and
15.755 eV; p. 5/1547) on the calibrated axes (``LINE_MARGIN``); Fig. 6's L and V against
Fig. 4's (the same curves drawn twice, ``SAME``). Reported: the pixel sizes; the length
and velocity sums of the three channels at the computed points, against McKenna and
van der Hart 2004 (TC-856), Pan, Gao and Starace 1990 (pan_1990b) and Moccia, Rahman
and Rizzo 1983 (moccia_1983).
"""

import numpy as np
import pan_1990b as pan
from _figure import image
from _library import Point, write

from optmat.datasets import DATA_DIRECTORY, Dataset, load_source

DIRECTORY = "Pindzola_1975_TwoPhotonIonization"
SOURCE = "pindzola_1975"
# sigma/F in 10^-18 cm^2 per atomic unit of flux, as cm^4 s: 1e-18 cm^2 times a0^2 t0
# (CODATA 2018), the inverse of the atomic unit of flux.
BOHR_CM = 5.29177210903e-9
TIME_AU_S = 2.4188843265857e-17
UNIT = 1e-18 * BOHR_CM**2 * TIME_AU_S
# The computed points (p. 5/1547), eV.
GRID = 0.5

_E = pan.Figure  # (page, table, unit, top, bottom, left, right, x_bottom, x_top,
# y_left, y_right, corners); the unit is that of the printed axis, 10^-18 cm^2/a.u.

FIGURES = {
    "fig2": _E(
        6,
        "fig_4c9504d1d62e",
        1.0,
        809,
        1909,
        672,
        2353,
        (
            (715.0, 7),
            (891.0, 8),
            (1067.5, 9),
            (1244.0, 10),
            (1424.5, 11),
            (1601.5, 12),
            (1779.0, 13),
            (1957.5, 14),
            (2133.5, 15),
            (2312.5, 16),
        ),
        (
            (719.0, 7),
            (894.5, 8),
            (1071.0, 9),
            (1247.0, 10),
            (1426.0, 11),
            (1601.5, 12),
            (1776.5, 13),
            (1952.5, 14),
            (2130.0, 15),
            (2307.5, 16),
        ),
        (
            (1815.5, 1),
            (1729.5, 2),
            (1640.0, 3),
            (1554.0, 4),
            (1465.0, 5),
            (1376.5, 6),
            (1286.5, 7),
            (1201.5, 8),
            (1116.5, 9),
            (1027.0, 10),
        ),
        (
            (1820.0, 1),
            (1732.5, 2),
            (1643.5, 3),
            (1556.0, 4),
            (1467.0, 5),
            (1377.5, 6),
            (1289.0, 7),
            (1202.5, 8),
            (1116.5, 9),
            (1028.0, 10),
        ),
        (None, None),
    ),
    "fig3": _E(
        6,
        "fig_7688b2052b30",
        1.0,
        4280,
        5377,
        667,
        2347,
        (
            (709.0, 7),
            (884.5, 8),
            (1061.5, 9),
            (1238.5, 10),
            (1416.5, 11),
            (1592.0, 12),
            (1770.5, 13),
            (1950.0, 14),
            (2127.5, 15),
            (2304.5, 16),
        ),
        (
            (712.5, 7),
            (887.0, 8),
            (1063.5, 9),
            (1239.5, 10),
            (1417.5, 11),
            (1591.0, 12),
            (1768.0, 13),
            (1945.5, 14),
            (2122.5, 15),
            (2300.5, 16),
        ),
        (
            (5285.5, 1),
            (5198.5, 2),
            (5113.5, 3),
            (5024.5, 4),
            (4936.5, 5),
            (4850.5, 6),
            (4759.5, 7),
            (4673.5, 8),
            (4586.0, 9),
            (4497.0, 10),
        ),
        (
            (5288.5, 1),
            (5200.5, 2),
            (5115.0, 3),
            (5026.5, 4),
            (4938.0, 5),
            (4850.5, 6),
            (4761.0, 7),
            (4673.5, 8),
            (4586.0, 9),
            (4497.5, 10),
        ),
        (None, None),
    ),
    "fig4": _E(
        6,
        "fig_17415fef50fa",
        1.0,
        812,
        1914,
        2747,
        4424,
        (
            (2786.5, 7),
            (2964.0, 8),
            (3141.5, 9),
            (3318.5, 10),
            (3495.5, 11),
            (3671.0, 12),
            (3847.5, 13),
            (4026.0, 14),
            (4203.0, 15),
            (4381.5, 16),
        ),
        (
            (2792.0, 7),
            (2968.0, 8),
            (3143.5, 9),
            (3320.0, 10),
            (3496.0, 11),
            (3670.5, 12),
            (3845.5, 13),
            (4023.5, 14),
            (4200.5, 15),
            (4378.0, 16),
        ),
        ((1029.0, 20), (1247.5, 15), (1468.0, 10), (1688.5, 5)),
        ((1034.0, 20), (1253.5, 15), (1474.5, 10), (1696.5, 5)),
        (None, None),
    ),
    "fig6": _E(
        7,
        "fig_e8afb9fc4099",
        1.0,
        3804,
        4902,
        2732,
        4412,
        (
            (2775.0, 7),
            (2951.5, 8),
            (3129.5, 9),
            (3305.0, 10),
            (3481.5, 11),
            (3659.5, 12),
            (3836.0, 13),
            (4014.5, 14),
            (4191.5, 15),
            (4369.5, 16),
        ),
        (
            (2777.5, 7),
            (2952.5, 8),
            (3130.5, 9),
            (3306.0, 10),
            (3481.5, 11),
            (3658.0, 12),
            (3833.0, 13),
            (4009.5, 14),
            (4187.0, 15),
            (4363.5, 16),
        ),
        ((4023.0, 20), (4241.0, 15), (4462.0, 10), (4680.5, 5)),
        ((4026.0, 20), (4244.5, 15), (4465.0, 10), (4684.0, 5)),
        (None, None),
    ),
}


def _file(name: str) -> str:
    """Return a dataset's data file."""
    return f"Ar/{SOURCE}__{name}_Ar.csv"


_S = pan.Segment
SEGMENTS: tuple[pan.Segment, ...] = (
    _S(_file("fig2_L_p1D"), "fig2", "solid", (935, 1656), (8.0, 11.0)),
    _S(_file("fig2_L_p1D"), "fig2", "solid", (1644, 1737), (12.0, 14.0)),
    _S(_file("fig2_V_p1D"), "fig2", "dashed", (935, 1702), (8.0, 11.0), (8.0, 9.5)),
    _S(_file("fig2_V_p1D"), "fig2", "dashed", (1644, 1789), (12.0, 14.0), (12.0, 13.2)),
    _S(_file("fig3_L_p1S"), "fig3", "solid", (929, 4720), (8.0, 11.0)),
    _S(_file("fig3_L_p1S"), "fig3", "solid", (1638, 5301), (12.0, 14.0)),
    _S(_file("fig3_V_p1S"), "fig3", "dashed", (929, 5071), (8.0, 11.0)),
    _S(_file("fig3_V_p1S"), "fig3", "dashed", (1638, 5329), (12.0, 13.0)),
    _S(_file("fig4_L_f1D"), "fig4", "solid", (3363, 1324), (8.0, 14.0), (8.0, 13.5)),
    _S(_file("fig4_V_f1D"), "fig4", "solid", (3363, 1581), (8.0, 14.0)),
    _S(_file("fig6_LGI_f1D"), "fig6", "solid", (3483, 4384), (8.0, 11.5)),
    _S(_file("fig6_LGI_f1D"), "fig6", "solid", (3704, 4487), (12.0, 14.0)),
    _S(_file("fig6_VGI_f1D"), "fig6", "dashed", (3483, 4432), (8.0, 11.5), (9.5, 11.5)),
    _S(
        _file("fig6_VGI_f1D"),
        "fig6",
        "dashed",
        (3704, 4446),
        (12.0, 14.0),
        (12.0, 13.0),
    ),
    _S("check:fig6_L_f1D", "fig6", "solid", (3483, 4294), (8.0, 14.0), (8.0, 13.5)),
    _S("check:fig6_V_f1D", "fig6", "dashed", (3483, 4553), (8.0, 14.0)),
)

# The computed points that cannot be read cleanly (another curve within a line width,
# the axis, or the line too steep), by file.
LEFT_OUT: dict[str, tuple[float, ...]] = {
    _file("fig2_L_p1D"): (9.5, 10.5, 11.0, 12.0, 13.5, 14.0),
    _file("fig2_V_p1D"): (9.5, 10.5, 11.0, 12.0, 13.0, 13.5, 14.0),
    _file("fig3_L_p1S"): (12.0, 14.0),
    _file("fig3_V_p1S"): (11.0, 12.0, 13.0),
    _file("fig4_L_f1D"): (9.5, 14.0),
    _file("fig4_V_f1D"): (14.0,),
    _file("fig6_LGI_f1D"): (8.5, 11.5, 14.0),
    _file("fig6_VGI_f1D"): (8.5, 14.0),
    "check:fig6_L_f1D": (8.0, 11.5, 14.0),
}
# The points read on two curves of a figure within a line width: left out of both.
TOGETHER: dict[str, tuple[float, ...]] = {
    _file("fig2_L_p1D"): (10.0,),
    _file("fig2_V_p1D"): (10.0,),
    _file("fig6_LGI_f1D"): (8.0, 9.0, 13.0, 13.5),
    _file("fig6_VGI_f1D"): (8.0, 9.0, 13.0, 13.5),
}


def _grid(span: tuple[float, float]) -> list[float]:
    """Return the computed points within a span, every ``GRID`` eV."""
    first, last = (round(x / GRID) for x in span)
    return [round(k * GRID, 2) for k in range(first, last + 1)]


def read_all(strict: bool = True) -> dict[str, list[pan.Reading]]:
    """Return every segment's readings, by data file, in photon-energy order.

    Raises if a segment's line style is not the one written, or (if ``strict``) the
    samples left out are not exactly those of :data:`LEFT_OUT` and :data:`TOGETHER`.
    """
    pages: dict[int, np.ndarray] = {}
    out: dict[str, list[pan.Reading]] = {}
    missing: dict[str, list[float]] = {}
    scales = calibrations(pages)
    for segment in SEGMENTS:
        figure = FIGURES[segment.figure]
        dark = pages[figure.page]
        calibration = scales[segment.figure]
        line = pan._Line(dark, figure, pan._follow(dark, calibration, segment, strict))
        for wanted in _grid(segment.span):
            reading = pan._read(calibration, line, wanted)
            if reading is None:
                missing.setdefault(segment.file, []).append(wanted)
            else:
                out.setdefault(segment.file, []).append(reading)
    together = _together(out)
    for name, samples in together.items():
        out[name] = [r for r in out[name] if r.wanted not in samples]
    for name in set(missing) | set(LEFT_OUT) | set(together) | set(TOGETHER):
        if strict and tuple(missing.get(name, ())) != LEFT_OUT.get(name, ()):
            msg = f"{name}: left out {missing.get(name)}, not {LEFT_OUT.get(name)}"
            raise ValueError(msg)
        if strict and tuple(together.get(name, ())) != TOGETHER.get(name, ()):
            msg = f"{name}: together {together.get(name)}, not {TOGETHER.get(name)}"
            raise ValueError(msg)
    if not strict:
        print("LEFT_OUT", {k: tuple(v) for k, v in missing.items()})
        print("TOGETHER", {k: tuple(v) for k, v in together.items()})
    return {name: sorted(rows, key=lambda r: r.wanted) for name, rows in out.items()}


def _together(out: dict[str, list[pan.Reading]]) -> dict[str, list[float]]:
    """Return the samples where two curves of a figure are read within a line width."""
    figure_of = {segment.file: segment.figure for segment in SEGMENTS}
    found: dict[str, set[float]] = {}
    names = sorted(out)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            if figure_of[a] != figure_of[b]:
                continue
            rows_b = {r.wanted: r.row for r in out[b]}
            for r in out[a]:
                if r.wanted in rows_b and abs(r.row - rows_b[r.wanted]) < pan.LINE:
                    found.setdefault(a, set()).add(r.wanted)
                    found.setdefault(b, set()).add(r.wanted)
    return {name: sorted(samples) for name, samples in found.items()}


def calibrations(
    pages: dict[int, np.ndarray] | None = None,
) -> dict[str, pan.Calibration]:
    """Return every figure's calibration, after checking its ticks in the image."""
    pages = {} if pages is None else pages
    out = {}
    for name, figure in FIGURES.items():
        if figure.page not in pages:
            pages[figure.page] = image(DIRECTORY, figure.page) < pan.DARK
        pan._check_ticks(pages[figure.page], figure)
        out[name] = pan.Calibration(figure, pages[figure.page])
    return out


# The energies (eV) the text prints for the lines drawn (p. 5/1547): the two-photon
# threshold A, the 4s, 5s and 3d 1P resonances, the single-photon threshold.
PRINTED = {
    "threshold 7.88": 7.88,
    "4s 1P 11.62": 11.62,
    "5s 1P 14.09": 14.09,
    "3d 1P 14.15": 14.15,
    "single-photon threshold 15.755": 15.755,
}
LINE_MARGIN = 0.03
# Fig. 6's line C, the 3d 1D autoionising resonance, 'obtained semiempirically from
# Moore's tables' (p. 7/1549): its energy is not printed; where it is drawn (eV).
UNPRINTED = {"fig6": 13.82}
# A column dark over this share of the frame's height is a drawn vertical line.
LINE_SHARE = 0.6
# How far Fig. 6's L and V may stray from Fig. 4's.
SAME = 0.03
COMPARE_AT = (8.0, 9.0, 10.0, 11.0)
OTHERS = {
    "mckenna_2004": "mckenna_2004/fig5_sigma2_Ar",
    "pan_1990b (all processes, L)": "pan_1990b/fig4c_all_L_Ar",
    "moccia_1983 L": "moccia_1983/fig2_sigma2_L_Ar",
    "moccia_1983 V": "moccia_1983/fig2_sigma2_V_Ar",
}


def lines(pages: dict[int, np.ndarray]) -> dict[str, list[float]]:
    """Return the photon energies (eV) of each figure's drawn vertical lines."""
    out = {}
    for name, calibration in calibrations(pages).items():
        figure = FIGURES[name]
        dark = pages[figure.page]
        columns = [
            c
            for c in range(figure.left + 8, figure.right - 8)
            if dark[figure.top + 5 : figure.bottom - 5, c].mean() > LINE_SHARE
        ]
        groups = np.split(np.array(columns), np.flatnonzero(np.diff(columns) > 2) + 1)
        middle = (figure.top + figure.bottom) / 2
        out[name] = [
            round(calibration.energy(float(g.mean()), middle), 3)
            for g in groups
            if g.size
        ]
    return out


def _values(readings: dict[str, list[pan.Reading]]) -> dict[str, dict[float, float]]:
    """Return each dataset's values, cm^4 s, by computed point, keyed by short name."""
    return {
        name.removeprefix("check:")
        .removeprefix(f"Ar/{SOURCE}__")
        .removesuffix("_Ar.csv"): {r.wanted: r.value * UNIT for r in rows}
        for name, rows in readings.items()
    }


def _dataset(key: str) -> Dataset:
    """Return another record's dataset by its key, ``<source>/<id>``."""
    source, _ = key.split("/")
    record = load_source(DATA_DIRECTORY / "sources" / f"{source}.toml")
    return next(d for d in record.datasets if d.key == key)


def _interpolate(dataset: Dataset, energy: float) -> float | None:
    """Return a dataset's value at a photon energy (eV), or None outside its range."""
    x = dataset.abscissa_values()
    if dataset.abscissa_unit == "hartree":
        x = x * 27.211386245988
    if not x.min() <= energy <= x.max():
        return None
    return float(np.interp(energy, x, dataset.scaled_values()))


def checks(readings: dict[str, list[pan.Reading]]) -> dict[str, object]:
    """Return the checks of the reading; raise if one fails its margin.

    Every drawn vertical line within ``LINE_MARGIN`` of an energy the text prints, or of
    :data:`UNPRINTED`; Fig. 6's L and V against Fig. 4's (``SAME``). Reported: the
    pixel sizes; the totals (the sum of the three channels at the computed points where
    all are read) against the other calculations of the library.
    """
    out: dict[str, object] = {}
    failures: list[str] = []
    pages: dict[int, np.ndarray] = {}
    for name, found in lines(pages).items():
        matched = []
        for energy in found:
            label = min(PRINTED, key=lambda k: abs(PRINTED[k] - energy))
            offset = round(energy - PRINTED[label], 3)
            if abs(offset) > LINE_MARGIN:
                unprinted = UNPRINTED.get(name)
                if unprinted is None or abs(energy - unprinted) > LINE_MARGIN:
                    failures.append(f"{name}: a line at {energy} eV")
                label, offset = "3d 1D autoionising (no energy printed)", 0.0
            matched.append((energy, label, offset))
        out[f"{name} lines (eV, printed, offset)"] = matched
    for name, calibration in calibrations(pages).items():
        x = abs(calibration.scales["bottom"][0])
        y = abs(calibration.scales["left"][0])
        out[f"{name} pixel (eV, printed unit), worst tick residual"] = (
            round(x, 4),
            round(y, 4),
            round(calibration.residual(), 1),
        )
    values = _values(readings)
    for gauge in ("L", "V"):
        again, first = values[f"fig6_{gauge}_f1D"], values[f"fig4_{gauge}_f1D"]
        shared = sorted(set(again) & set(first))
        ratios = [again[e] / first[e] for e in shared]
        out[f"fig6 {gauge} / fig4 {gauge}"] = (
            shared[0],
            shared[-1],
            round(min(ratios), 3),
            round(max(ratios), 3),
            len(shared),
        )
        if max(abs(r - 1) for r in ratios) > SAME:
            failures.append(
                f"fig6 {gauge}: {min(ratios):.3f}-{max(ratios):.3f} of fig4"
            )
    totals = {
        "L (second order)": ("fig2_L_p1D", "fig3_L_p1S", "fig4_L_f1D"),
        "V (second order)": ("fig2_V_p1D", "fig3_V_p1S", "fig4_V_f1D"),
        "L with LGI f 1D": ("fig2_L_p1D", "fig3_L_p1S", "fig6_LGI_f1D"),
    }
    others = {label: _dataset(key) for label, key in OTHERS.items()}
    for label, parts in totals.items():
        shared = sorted(set.intersection(*(set(values[p]) for p in parts)))
        total = {e: sum(values[p][e] for p in parts) for e in shared}
        out[f"total {label}, cm4 s"] = [(e, float(f"{total[e]:.3e}")) for e in shared]
        for other, dataset in others.items():
            ratios = []
            for e in COMPARE_AT:
                theirs = _interpolate(dataset, e) if e in total else None
                if theirs:
                    ratios.append((e, round(total[e] / theirs, 2)))
            out[f"total {label} / {other}"] = ratios
    if failures:
        raise ValueError("; ".join(failures))
    return out


def _points(readings: dict[str, list[pan.Reading]]) -> dict[str, list[Point]]:
    """Return the data files' points: sigma/F as printed against photon energy (eV)."""
    return {
        name: [
            Point(
                f"{r.wanted:.1f}",
                f"{r.energy:.3f}",
                f"{r.value:.3f}",
                "",
                r.printed,
            )
            for r in rows
        ]
        for name, rows in readings.items()
        if not name.startswith("check:")
    }


def extract() -> dict[str, list[Point]]:
    """Return the partial cross sections read from Figs. 2, 3, 4 and 6."""
    readings = read_all()
    checks(readings)
    return _points(readings)


def main() -> None:
    """Write the data files and print the checks."""
    readings = read_all()
    results = checks(readings)
    write(_points(readings))
    for label, result in results.items():
        print(f"{label}: {result}")


if __name__ == "__main__":
    main()
