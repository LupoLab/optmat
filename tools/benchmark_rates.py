# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The time of a rate evaluation against the number of fields, and against another tree.

    uv run python tools/benchmark_rates.py
    uv run python tools/benchmark_rates.py --against <checkout>/src --rounds 5

times ``instantaneous_rate`` of the default argon model and of the default N2 and H2
molecular models at 800 nm, at 1, 64 and 4096 fields: the best of several repeats, in
ms. With ``--against``, the same timings of the package under that source directory
(an older commit, extracted with ``git archive <commit> src | tar -x -C <dir>``) are
taken in alternation with this tree's, round by round in fresh processes, so that a
loaded machine slows both alike; the table gives each one's best round and their
ratio. The rates of the two are also compared, and any difference is reported: an
optimisation that keeps the numbers keeps them to the last bit (TC-837).
"""

import argparse
import json
import os
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np

SIZES = (1, 64, 4096)
CASES = (("Ar", "atom"), ("N2", "molecule"), ("H2", "molecule"))
WAVELENGTH = 800e-9


def _model(gas: str, kind: str) -> object:
    # Imported here so that a subprocess times the tree on its own path.
    if kind == "atom":
        from optmat.gas.ionisation_rates import ionisation_rate_model  # noqa: PLC0415

        return ionisation_rate_model(gas, WAVELENGTH)
    from optmat.gas.molecular_ionisation_rates import (  # noqa: PLC0415
        molecular_rate_model,
    )

    return molecular_rate_model(gas, WAVELENGTH)


def measure(repeats: int = 5) -> dict[str, dict[str, object]]:
    """Return the best time (ms) and the rates of each case, in this process.

    Parameters
    ----------
    repeats : int, optional
        Timed repeats of each case; the best is kept.

    Returns
    -------
    dict
        ``{"<gas>/<size>": {"ms": float, "rates": list of str}}``, the rates as
        ``float.hex`` strings.
    """
    warnings.simplefilter("ignore")
    out: dict[str, dict[str, object]] = {}
    for gas, kind in CASES:
        model = _model(gas, kind)
        rate = model.instantaneous_rate  # pyright: ignore[reportAttributeAccessIssue]
        for size in SIZES:
            field = np.linspace(1e9, 6e10, size)
            values = rate(field)
            calls = max(5, 400 // size)
            best = float("inf")
            for _ in range(repeats):
                start = time.perf_counter()
                for _ in range(calls):
                    rate(field)
                best = min(best, (time.perf_counter() - start) / calls)
            out[f"{gas}/{size}"] = {
                "ms": best * 1e3,
                "rates": [float(v).hex() for v in np.ravel(values)],
            }
    return out


def _run(source: str | None) -> dict[str, dict[str, object]]:
    env = dict(os.environ)
    if source is not None:
        env["PYTHONPATH"] = source
    result = subprocess.run(
        [sys.executable, __file__, "--measure"],
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(result.stdout)


def main() -> None:
    """Print the timings, and the comparison with another tree where asked."""
    parser = argparse.ArgumentParser(description="Time the rate evaluations.")
    parser.add_argument("--against", help="another tree's src directory")
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--measure", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.measure:
        print(json.dumps(measure()))
        return
    trees: dict[str, str | None] = {"this": None}
    if args.against:
        trees["against"] = str(Path(args.against).resolve())
    best: dict[str, dict[str, float]] = {name: {} for name in trees}
    rates: dict[str, dict[str, object]] = {}
    for _ in range(args.rounds):
        for name, source in trees.items():
            for case, entry in _run(source).items():
                ms = float(entry["ms"])  # pyright: ignore[reportArgumentType]
                best[name][case] = min(best[name].get(case, ms), ms)
                rates.setdefault(name, {})[case] = entry["rates"]
    header = f"{'case':12s}" + "".join(f"{name:>12s}" for name in trees)
    print(header + ("       ratio" if args.against else ""))
    for case in best["this"]:
        row = f"{case:12s}" + "".join(f"{best[n][case]:10.3f}ms" for n in trees)
        if args.against:
            row += f"{best['against'][case] / best['this'][case]:11.2f}x"
        print(row)
    if args.against:
        differ = [c for c in rates["this"] if rates["this"][c] != rates["against"][c]]
        print("rates identical to the last bit" if not differ else f"DIFFER: {differ}")


if __name__ == "__main__":
    main()
