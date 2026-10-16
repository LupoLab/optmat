# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The molecular fits' representation error and absorption in their gaps (OEP 3, Q13).

    uv run python tools/hitran/gap_bound.py [CO CO2 CH4 N2O]

A mixture whose constituent absorbs only as a trace can be transparent where the
constituent's own model is not qualified. Widening the mixture's intervals there needs,
for that constituent, its absorption and its representation error in its gaps (OEP 3,
Q13): ``windows.py`` stores both only in its windows. This tool measures them in the
gaps, from 12.5 um (800 cm^-1) to the top of the absorption grid, against the same
line-by-line reference smoothed at sigma = 1 cm^-1, for the stored fit (its electronic
terms and its fixed infrared terms as fitted). It writes
``src/optmat/data/external/hitran/<gas>_gap_bound.csv``: geometric cells (ratio 1.001)
with the largest |D_m[fit] - D_m[reference]| (m = 0, 2, 3, au) and the largest smoothed
absorption of the pure gas at 101325 Pa and 296 K (m^-1) on the 0.25 cm^-1 grid and on
the grid offset by half a step; a cell holding a pole of the fit has ``pole = 1``, and
its bound is unbounded. :func:`optmat.gas.air.air_model` reads the tables.
"""

import itertools
import sys
from types import SimpleNamespace

import numpy as np

from optmat.gas.models.compact import parent_digest
from optmat.gas.models.fitted import fitted_models

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import windows

VERSION = 1


def _values(reference: windows.Reference, fit: object, nu: np.ndarray) -> dict:
    """Return |D_m[fit] - D_m[reference]| and the absorption at points of the grid."""
    sample = reference.sample(nu)
    omega = nu / windows.HCM
    model = {
        0: fit.polarizability_au(omega),  # pyright: ignore[reportAttributeAccessIssue]
        2: fit.dispersion_au(omega, 2),  # pyright: ignore[reportAttributeAccessIssue]
        3: fit.dispersion_au(omega, 3),  # pyright: ignore[reportAttributeAccessIssue]
    }
    out: dict[int | str, np.ndarray] = {
        m: np.abs(model[m] - sample[f"d{m}"]) for m in (0, 2, 3)
    }
    out["absorption"] = sample["absorption"]
    return out


def gaps(fit: object, settings: windows.Settings) -> list[tuple[float, float]]:
    """Return the fit's gaps in wavenumber, inside the absorption grid."""
    low_cm, high_cm = settings.low_cm, settings.top_cm - 50
    spans = sorted(
        (1e-2 / b, 1e-2 / a)
        for a, b in fit.intervals  # pyright: ignore[reportAttributeAccessIssue]
    )
    found, start = [], low_cm
    for a, b in spans:
        if a > start:
            found.append((start, min(a, high_cm)))
        start = max(start, b)
    if start < high_cm:
        found.append((start, high_cm))
    return [(a, b) for a, b in found if b > a]


def main(gases: list[str]) -> None:
    """Derive and write the gap tables of each gas."""
    settings = windows.Settings()
    for gas in gases:
        fit = fitted_models()[f"fit/{gas}"]
        lines = windows.load_lines(gas)
        parent = windows.load_parent(gas)
        reference = windows.build_reference(lines, parent, settings)
        poles = np.asarray(fit.energies) * windows.HCM
        h = settings.grid_cm
        rows = []
        for a, b in gaps(fit, settings):
            nu = np.arange(a, b + h / 2, h)
            nu = nu[nu <= reference.fine_nu[-1]]
            if nu.size == 0:
                continue
            on = _values(reference, fit, nu)
            off = _values(reference, fit, np.clip(nu + h / 2, None, nu[-1]))
            edges = [a]
            while edges[-1] < b:
                edges.append(min(b, edges[-1] * settings.cell_ratio))
            for low, high in itertools.pairwise(edges):
                inside = (nu >= low) & (nu <= high)
                if not inside.any():  # a cell narrower than the grid: its neighbours
                    j = int(np.searchsorted(nu, low))
                    inside[max(j - 1, 0) : j + 1] = True
                pole = bool(np.any((poles >= low) & (poles <= high)))
                row: dict[str, object] = {
                    "low_cm": f"{low:.3f}",
                    "high_cm": f"{high:.3f}",
                }
                for key in (0, 2, 3):
                    worst = max(
                        float(on[key][inside].max()), float(off[key][inside].max())
                    )
                    row[f"b{key}"] = "inf" if pole else windows._up(worst)
                absorption = max(
                    float(on["absorption"][inside].max()),
                    float(off["absorption"][inside].max()),
                )
                row["absorption_m"] = (
                    windows._up(absorption) if absorption > 0 else "0.0"
                )
                row["pole"] = int(pole)
                rows.append(row)
        # The provenance reads the gas, the parent and the settings alone.
        analysis = SimpleNamespace(gas=gas, parent=parent, settings=settings)
        header = [
            line.replace(
                "# Derived by tools/hitran/windows.py (version 1",
                f"# Derived by tools/hitran/gap_bound.py (version {VERSION}",
            )
            for line in windows.provenance(analysis)  # pyright: ignore[reportArgumentType]
        ]
        header.insert(4, f"# fit {fit.name}: stored sha256 {parent_digest(fit.name)}")
        header.append(
            "# cells: geometric, in the fit's gaps; the largest |D_m - D_m ref| (au) "
            "and the largest smoothed absorption (m^-1) in each; pole = 1 where a pole "
            "of the fit lies in the cell (unbounded)"
        )
        fields = ["low_cm", "high_cm", "b0", "b2", "b3", "absorption_m", "pole"]
        path = windows.OUT / f"{gas}_gap_bound.csv"
        windows._write(path, header, fields, rows)
        count = len(gaps(fit, settings))
        print(f"{gas}: {len(rows)} cells in {count} gaps -> {path.name}")


if __name__ == "__main__":
    main(sys.argv[1:] or list(windows.GASES))
