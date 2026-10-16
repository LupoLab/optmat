# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Vogl & Hall (1972): the He-Ne system by the Burnett method at 40, 0, -40 C.

Table I (``tbl_b47536f1a05a``, p. 530, cells file; image checked) prints, for 100.00 %
helium and for the 49.4 mole % helium-neon mixture, P (atm), Z and sigma_z (x 10^5) at
+40, 0 and -40 C, in three column triples. Each isotherm is two Burnett runs, printed
one after the other with a blank line between them on the page; the cells file drops
the blank line, so a run is recognised here where the pressure rises again. Only the
mixture is extracted (the pure gases' Z are left to the reference equations of state),
one file per isotherm, rows ``run<r>_<k>`` in the printed order; sigma_z is written in
units of Z (the printed integer times 10^-5), the raw cells kept in ``printed``.

Table II (``tbl_00584e89d11b``, p. 531) prints B_2 (sigma_B2) per composition and
temperature: one file each for pure He, pure Ne and the mixture, rows by the
temperature in C. Table III (``tbl_aa7d07823f8a``, p. 532) prints B_12 at the three
temperatures, without an uncertainty. The B_3 rows of Table II are not extracted (no
third-virial quantity in the vocabulary).
"""

import re
from decimal import Decimal

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Vogl_1972_CompressibilityDataVirial"
TABLE_I = "tbl_b47536f1a05a"
TABLE_II = "tbl_00584e89d11b"
TABLE_III = "tbl_aa7d07823f8a"
TEMPERATURES = ("40", "0", "-40")
"""The isotherms of Tables I and II, in C, in the order of their column triples."""
NAMES = {"40": "40C", "0": "0C", "-40": "m40C"}
MIXTURE_TITLE = "Compressibility data for 49.4% helium"
_PARENTHESISED = re.compile(r"^\(([\d.]+)\)$")


def _mixture_isotherms() -> dict[str, list[Point]]:
    """Return the mixture's Z of Table I, one file per isotherm."""
    grid = json_table(DIRECTORY, TABLE_I)
    # The header's spans are uneven (the first P (atm) spans two rows): the
    # temperatures sit at columns 1, 3 and 6, the triples at 0-2, 3-5 and 6-8.
    if [grid[1][i].replace("\N{MINUS SIGN}", "-") for i in (1, 3, 6)] != [
        "+40°C",
        "0°C",
        "-40°C",
    ]:
        msg = f"unexpected Table I temperatures {grid[1]!r}"
        raise ValueError(msg)
    start = next(i for i, cells in enumerate(grid) if cells[0] == MIXTURE_TITLE)
    files: dict[str, list[Point]] = {}
    for k, temperature in enumerate(TEMPERATURES):
        points: list[Point] = []
        run, index, previous = 1, 0, float("inf")
        for cells in grid[start + 1 :]:
            pressure, z, sigma = cells[3 * k : 3 * k + 3]
            if not pressure:
                continue
            if float(pressure) > previous:
                # The second Burnett run begins at a higher pressure again.
                run, index = run + 1, 0
            previous, index = float(pressure), index + 1
            points.append(
                Point(
                    f"run{run}_{index}",
                    number(pressure),
                    number(z),
                    f"{Decimal(number(sigma)).scaleb(-5):f}",
                    f"{pressure} | {z} | {sigma}",
                )
            )
        if run != 2:
            msg = f"Table I, 49.4 % He at {temperature} C: {run} runs, not two"
            raise ValueError(msg)
        files[f"He/vogl_1972__Z_x0494_{NAMES[temperature]}.csv"] = points
    return files


def _second_virials() -> dict[str, list[Point]]:
    """Return Table II's B_2 of He, Ne and the mixture against temperature."""
    grid = json_table(DIRECTORY, TABLE_II)
    blocks = {
        "100% He": "He/vogl_1972__B_He.csv",
        "49.4% He": "He/vogl_1972__B_mix_x0494.csv",
        "100% Ne": "Ne/vogl_1972__B_Ne.csv",
    }
    files: dict[str, list[Point]] = {}
    for i, cells in enumerate(grid):
        if cells[3] not in blocks:
            continue
        row = grid[i + 1]
        if not row[0].startswith("B_{2}"):
            msg = f"Table II: no B_2 row under {cells[3]!r}"
            raise ValueError(msg)
        points = []
        for k, temperature in enumerate(TEMPERATURES):
            value, sigma = row[1 + 2 * k], row[2 + 2 * k]
            match = _PARENTHESISED.match(sigma)
            if match is None:
                msg = f"Table II: unexpected sigma cell {sigma!r}"
                raise ValueError(msg)
            points.append(
                Point(
                    temperature,
                    temperature,
                    number(value),
                    match.group(1),
                    f"{grid[1][1 + 2 * k]} | {value} | {sigma}",
                )
            )
        files[blocks[cells[3]]] = points
    if len(files) != len(blocks):
        msg = f"Table II: found {sorted(files)}"
        raise ValueError(msg)
    return files


def _interaction() -> dict[str, list[Point]]:
    """Return Table III's B_12 against temperature."""
    grid = json_table(DIRECTORY, TABLE_III)
    if grid[1] != ["T(°C)", "B_{12}(cm3/mole)"]:
        msg = f"unexpected Table III header {grid[1]!r}"
        raise ValueError(msg)
    points = []
    for cells in grid[2:]:
        temperature, value = number(cells[0]), number(cells[1])
        points.append(Point(temperature, temperature, value, "", " | ".join(cells)))
    if [p.row for p in points] != list(TEMPERATURES):
        msg = f"Table III temperatures {[p.row for p in points]!r}"
        raise ValueError(msg)
    return {"He/vogl_1972__B12_Ne.csv": points}


def extract() -> dict[str, list[Point]]:
    """Return the mixture's Z isotherms, the B_2 of Table II and the B_12."""
    return _mixture_isotherms() | _second_virials() | _interaction()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
