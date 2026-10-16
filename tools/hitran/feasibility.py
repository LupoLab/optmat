# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The feasibility tables of OEP 3, section 8, recomputed from the cached line lists.

    uv run python tools/hitran/feasibility.py CO2 4.65 12.2 [--keep 1e-7] [--sigma 1]

For one gas and candidate window (vacuum wavelengths in um) this prints, per number of
nodes per band (1, 2, 4, 8):

- **against the lines outside the window** (the first table of section 8): every band
  split at the window and replaced by its rule, compared with the undamped sum over the
  same lines; the error in n - 1 (core maximum, the core excluding 2 % of the window at
  each edge) and in the total D_2 over the electronic |D_2| (median, and the fraction
  below 1 %);
- **against the complete reference** (the second table): the smoothed sum over every
  line, the in-window bands folded (or, with ``--keep``, those above that fraction of
  the molecule's M_0 kept), the rotational system as a 4-node rule; the fraction of
  the window passing the tolerance of section 7.

These are 1500 geometric points per window and a fixed rule per band, as section 8
states; ``windows.py`` derives the stored windows on its own grids and adaptive rules.
"""

import argparse
from collections import defaultdict

import numpy as np
import windows as w

POINTS = 1500


def _rules(
    f: np.ndarray, e: np.ndarray, groups: list[np.ndarray], n: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return the poles of every group's n-node rule."""
    parts = [w.gauss_rule(f[g], e[g], n) for g in groups]
    return np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts])


def table(  # noqa: PLR0915 - the two tables of section 8, step by step
    gas: str, low_um: float, high_um: float, keep: float | None, sigma: float
) -> None:
    """Print both feasibility tables of one window."""
    lines = w.load_lines(gas)
    parent = w.load_parent(gas)
    f, e = w.oscillator(lines.nu, lines.strength)
    rotational = lines.rotational[lines.band]
    a_cm, b_cm = 1e4 / high_um, 1e4 / low_um
    inside = ~rotational & (lines.nu >= a_cm) & (lines.nu <= b_cm)
    outside = ~rotational & ~inside
    groups = defaultdict(list)
    for i in np.flatnonzero(outside):
        groups[(int(lines.band[i]), bool(lines.nu[i] > b_cm))].append(i)
    group_index = [np.array(g) for g in groups.values()]
    lam = np.geomspace(low_um, high_um, POINTS)
    nu = 1e4 / lam
    omega = nu / w.HCM
    core = (lam > low_um * 1.02) & (lam < high_um / 1.02)
    electronic = parent.electronic(omega)
    exact_outside = w.pole_dispersion(f[outside], e[outside], omega)
    settings = w.Settings(sigma_cm=sigma)
    reference_lines = w.dispersion(
        w.direct_derivatives(lines.nu, lines.strength, sigma, nu),
        omega,
    )
    reference = {m: electronic[m] + reference_lines[m] for m in (0, 2, 3)}
    tol0 = settings.alpha_fraction * parent.sigma(omega)
    tol2 = settings.gvd_fraction * np.maximum(
        np.abs(reference[2]), np.abs(electronic[2])
    )
    total_m0 = float(f[~rotational].sum())
    bands_inside = defaultdict(list)
    for i in np.flatnonzero(inside):
        bands_inside[int(lines.band[i])].append(i)
    kept = [
        np.array(ix)
        for ix in bands_inside.values()
        if keep is not None and f[ix].sum() / total_m0 > keep
    ]
    folded = [
        np.array(ix)
        for ix in bands_inside.values()
        if not (keep is not None and f[ix].sum() / total_m0 > keep)
    ]
    print(
        f"{gas}, {low_um}-{high_um} um: {inside.sum()} lines inside in "
        f"{len(bands_inside)} bands ({f[inside].sum() / total_m0:.1e} of M_0); "
        f"{len(kept)} kept, {len(folded)} folded; sigma {sigma} cm^-1"
    )
    for n in (1, 2, 4, 8):
        pf, pe = _rules(f, e, group_index, n)
        rule = w.pole_dispersion(pf, pe, omega)
        err0 = np.abs(rule[0] - exact_outside[0]) / np.abs(
            electronic[0] + exact_outside[0]
        )
        err2 = np.abs(rule[2] - exact_outside[2]) / np.abs(electronic[2])
        # The complete reference: the rules, the kept bands' rules, the fold, rotation.
        kf, ke = _rules(f, e, kept, n) if kept else (np.zeros(0), np.zeros(0))
        poles_f, poles_e = np.concatenate([pf, kf]), np.concatenate([pe, ke])
        if folded:
            index = np.concatenate(folded)
            below = poles_e[poles_e < a_cm / w.HCM]
            above = poles_e[poles_e > b_cm / w.HCM]
            upper = above.min() if len(above) else float(parent.energies.min())
            low, high = w.fold_poles(
                float(f[index].sum()),
                float(np.sum(f[index] / e[index] ** 2)),
                float(below.max()),
                float(upper),
            )
            poles_f = np.concatenate([poles_f, [low, high]])
            poles_e = np.concatenate([poles_e, [below.max(), upper]])
        if rotational.any():
            rf, re_ = w.gauss_rule(f[rotational], e[rotational], 4)
            poles_f, poles_e = (
                np.concatenate([poles_f, rf]),
                np.concatenate([poles_e, re_]),
            )
        model = w.pole_dispersion(poles_f, poles_e, omega)
        passing = (np.abs(electronic[0] + model[0] - reference[0]) <= tol0) & (
            np.abs(electronic[2] + model[2] - reference[2]) <= tol2
        )
        if kept:
            for ix in kept:
                span = (lines.nu[ix].min() - sigma, lines.nu[ix].max() + sigma)
                passing &= (nu < span[0]) | (nu > span[1])
        print(
            f"  {n} nodes: outside lines: n - 1 core max {err0[core].max():.1e}, "
            f"GVD median {np.median(err2):.1e}, < 1 % over {np.mean(err2 < 1e-2):.0%} "
            f"(core {np.mean(err2[core] < 1e-2):.0%}); complete reference: passing "
            f"{np.mean(passing):.0%} (core {np.mean(passing[core]):.0%})"
            + (
                f", {lam[passing].min():.2f}-{lam[passing].max():.2f} um"
                if passing.any()
                else ""
            )
        )


def main() -> None:
    """Parse the arguments and print the tables."""
    parser = argparse.ArgumentParser(description="The feasibility tables of OEP 3.")
    parser.add_argument("gas")
    parser.add_argument("low_um", type=float)
    parser.add_argument("high_um", type=float)
    parser.add_argument("--keep", type=float, default=None)
    parser.add_argument("--sigma", type=float, default=1.0)
    args = parser.parse_args()
    table(args.gas, args.low_um, args.high_um, args.keep, args.sigma)


if __name__ == "__main__":
    main()
