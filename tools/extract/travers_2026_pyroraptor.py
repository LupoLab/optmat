# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract J. C. Travers's pyroraptor helium yields from the archived result tables.

Like ``travers_2026_trecx`` this script does not read the paper library: the source is
a set of calculations, archived in ``data/external/pyroraptor`` (``PROVENANCE.md``
there), whose four tables hold the runs' values as written by ``repr``. Kept as
strings:

* at 390 nm, in A2's 6 fs Gaussian (rows ``3p68e14`` and ``5p69e14``), 1000 au after
  the pulse, the populations above the shells n <= 4 (``above4``: ionisation plus
  excitation to n >= 5) and n <= 5 (``above5``) of
  - both electrons (``twoelectron_``), the largest box, R0 = 97 bohr
    (:data:`TWO_ELECTRON_BOX`);
  - haCC with one ion channel (``hacc1_``) and with the n <= 2 channels (``hacc5_``),
    neutral square 6;
  - one electron (``sae_``, ``above4`` only), the run :data:`SAE_RUN`;
* at 800 nm, in Sato et al.'s six-cycle sin^2 pulse (rows ``I4e14``, ``I8e14`` and
  ``I14e14``, as ``sato_2016``), 1 - P(both electrons inside 20 bohr) at the end of the
  pulse of haCC(1) (square 6) and haCC(5) (``hacc1_`` and ``hacc5_outside20_``).

:func:`uncertainties` gives each dataset's relative uncertainty by the rule its record
states; :func:`checks` asserts what the record says of the runs.
"""

import csv
import math

from _library import DATA, Point, write

DIRECTORY = DATA / "external" / "pyroraptor"
SOURCE = "travers_2026_pyroraptor"
ROWS_390 = {"3p68e14": "3.684031e14", "5p69e14": "5.68981e14"}
SAE_LABELS = {"3.684e14": "3p68e14", "5.690e14": "5p69e14"}  # the driver's K2 labels
ROWS_800 = {"I4e14": "4e14", "I8e14": "8e14", "I14e14": "14e14"}
TWO_ELECTRON = "he_two_electron_390nm.csv"
TWO_ELECTRON_BOX = "97"
TWO_ELECTRON_END = 1744.144720033084  # au, 1000 au after the pulse
TWO_ELECTRON_TIME = 1e-6  # au: the two intensities' exports differ in the last digit
HACC_390 = "he_hacc_390nm.csv"
HACC_390_RUNS = {
    "hacc1": {"3p68e14": "l_sq6", "5p69e14": "l_hi_sq6"},
    "hacc5": {"3p68e14": "h5_sq6", "5p69e14": "h5_sq6_hi"},
}
HACC_390_BASE = "l_base"  # one channel, square 3: the convergence variations' base
HACC_390_LADDER = ("l_sq5", "l_sq6")  # the last step of the neutral-square ladder
HACC_390_BOXES = ("l_base", "l_r60", "l_r80")
SAE = "he_sae_390nm.csv"
SAE_RUN = "o13"
HACC_800 = "he_hacc_800nm.csv"
HACC_800_RUNS = {"hacc1": "h1/{}-sq6", "hacc5": "h5/{}-best"}
HACC_800_BASE = "h1/{}-base"  # square 3
HACC_800_VARIATIONS = (
    "dt2", "dt3", "ord12", "ord14", "r0hi", "r0lo", "lhi", "llo", "cut100", "cut400",
    "rg12", "rg28",
)  # fmt: skip
HACC_800_LADDER = ("h1/{}-base", "h1/{}-sq6")
HACC5_800_LMAX = ("h5/8e14-best", "h5/8e14-lhi")
TAGS_800 = {"I4e14": "4e14", "I8e14": "8e14", "I14e14": "1.4e15"}
# The two-electron calculation's angular basis: the inner-arm ladder of the absorbed
# norm-loss at 3.684e14 W/cm2, BOX 40, inner = 2 -> 3 -> 4, as printed in pyroraptor's
# two-electron report (-0.97 % and -0.67 %); its geometric tail puts the limit 3.1 %
# below inner = 2. Not run at 5.69e14 (TC-1092).
ANGULAR_STEPS = (-0.0097, -0.0067)


def _rows(name: str) -> list[dict[str, str]]:
    with (DIRECTORY / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _two_electron() -> dict[tuple[str, str], dict[str, str]]:
    """Return the two-electron rows 1000 au after the pulse, by (row, box)."""
    out = {}
    for row in _rows(TWO_ELECTRON):
        if abs(float(row["time_au"]) - TWO_ELECTRON_END) > TWO_ELECTRON_TIME:
            continue
        label = "3p68e14" if row["record"].startswith("k2b-3p68e14") else "5p69e14"
        out[label, row["box_r0"]] = row
    return out


def _point(row: str, x: str, value: str, *cells: str) -> Point:
    float(value)  # a plain number, or the extraction stops
    return Point(row, x, value, "", " | ".join(cells))


def extract() -> dict[str, list[Point]]:
    """Return the helium yields of every calculation and measure."""
    files: dict[str, list[Point]] = {}
    two = _two_electron()
    for cut in ("4", "5"):
        column = f"above_n{cut}"
        files[f"He/{SOURCE}__twoelectron_above{cut}_gauss6_390nm_He.csv"] = [
            _point(
                label,
                ROWS_390[label],
                (r := two[label, TWO_ELECTRON_BOX])[column],
                TWO_ELECTRON,
                r["record"],
                f"time_au {r['time_au']}",
                f"{column} {r[column]}",
            )
            for label in ROWS_390
        ]
        hacc = {row["record"]: row for row in _rows(HACC_390)}
        for model, runs in HACC_390_RUNS.items():
            files[f"He/{SOURCE}__{model}_above{cut}_gauss6_390nm_He.csv"] = [
                _point(
                    label,
                    ROWS_390[label],
                    hacc[name][column],
                    HACC_390,
                    name,
                    f"{column} {hacc[name][column]}",
                )
                for label, name in runs.items()
            ]
    sae = {
        SAE_LABELS[row["intensity_label"]]: row
        for row in _rows(SAE)
        if row["record"] == SAE_RUN
    }
    files[f"He/{SOURCE}__sae_above4_gauss6_390nm_He.csv"] = [
        _point(
            label,
            ROWS_390[label],
            sae[label]["above_n4"],
            SAE,
            SAE_RUN,
            f"intensity_label {sae[label]['intensity_label']}",
            f"above_n4 {sae[label]['above_n4']}",
        )
        for label in ROWS_390
    ]
    sato = {row["record"]: row for row in _rows(HACC_800)}
    for model, pattern in HACC_800_RUNS.items():
        files[f"He/{SOURCE}__{model}_outside20_sin2_800nm_He.csv"] = [
            _point(
                label,
                intensity,
                (r := sato[pattern.format(TAGS_800[label])])["beyond_rc"],
                HACC_800,
                r["record"],
                f"intensity_wcm2 {r['intensity_wcm2']}",
                f"beyond_rc {r['beyond_rc']}",
            )
            for label, intensity in ROWS_800.items()
        ]
    return files


def angular_excess() -> float:
    """Return how far the two-electron yield at inner = 2 lies above its angular limit.

    The geometric extrapolation of :data:`ANGULAR_STEPS`: the two steps and the tail
    beyond the second, d2 r / (1 - r) with r = d2 / d1.
    """
    first, second = ANGULAR_STEPS
    ratio = second / first
    return -(first + second + second * ratio / (1 - ratio))


def uncertainties() -> dict[str, dict[str, float]]:
    """Return each dataset's relative uncertainty, by row.

    * Two electrons: the angular excess (:func:`angular_excess`), which the box adds
      to (the change from R0 = 78 to 97 bohr).
    * haCC at 390 nm: the largest change of the measure among the one-channel
      variations at 3.684e14 W/cm2 (time step, order, cut, angle, lMax, rMax, the
      gauge radius) and the last step of the neutral-square ladder, for every row and
      both channel sets.
    * One electron at 390 nm: the largest change among the four discretisations.
    * haCC at 800 nm: per intensity, the largest change among haCC(1)'s variations
      and the last step of its square ladder; haCC(5) adds its lMax step at 8e14 to
      every row.
    """
    out: dict[str, dict[str, float]] = {}
    two = _two_electron()
    hacc = {row["record"]: row for row in _rows(HACC_390)}
    base = hacc[HACC_390_BASE]
    variations = [
        name
        for name, row in hacc.items()
        if row["ion_principal"] == "1"
        and row["intensity_wcm2"] == ROWS_390["3p68e14"]
        and row["square"] == base["square"]
        and name != HACC_390_BASE
    ]
    for cut in ("4", "5"):
        column = f"above_n{cut}"
        out[f"twoelectron_above{cut}_gauss6_390nm_He"] = {
            label: angular_excess()
            + abs(
                float(two[label, TWO_ELECTRON_BOX][column])
                / float(two[label, "78"][column])
                - 1
            )
            for label in ROWS_390
        }
        spread = max(
            abs(float(hacc[name][column]) / float(base[column]) - 1)
            for name in variations
        )
        old, new = HACC_390_LADDER
        square = abs(float(hacc[new][column]) / float(hacc[old][column]) - 1)
        for model in HACC_390_RUNS:
            out[f"{model}_above{cut}_gauss6_390nm_He"] = dict.fromkeys(
                ROWS_390, spread + square
            )
    sae = _rows(SAE)
    out["sae_above4_gauss6_390nm_He"] = {
        label: max(
            abs(float(row["above_n4"]) / float(own["above_n4"]) - 1)
            for row in sae
            if row["intensity_label"] == own["intensity_label"]
        )
        for own in sae
        if own["record"] == SAE_RUN
        for label in [SAE_LABELS[own["intensity_label"]]]
    }
    sato = {row["record"]: row for row in _rows(HACC_800)}
    one: dict[str, float] = {}
    for label, tag in TAGS_800.items():
        reference = float(sato[HACC_800_BASE.format(tag)]["beyond_rc"])
        spread = max(
            abs(float(sato[f"h1/{tag}-{v}"]["beyond_rc"]) / reference - 1)
            for v in HACC_800_VARIATIONS
        )
        old, new = (name.format(tag) for name in HACC_800_LADDER)
        square = abs(float(sato[new]["beyond_rc"]) / float(sato[old]["beyond_rc"]) - 1)
        one[label] = spread + square
    best, wider = (sato[name]["beyond_rc"] for name in HACC5_800_LMAX)
    lmax = abs(float(wider) / float(best) - 1)
    out["hacc1_outside20_sin2_800nm_He"] = one
    out["hacc5_outside20_sin2_800nm_He"] = {k: v + lmax for k, v in one.items()}
    return out


def checks() -> None:
    """Assert what the source record states about the runs.

    The two-electron ``above4`` is box-converged from R0 = 59 bohr (within 0.1 %), while
    the norm-loss-like ``above6`` is not; haCC's ``above4`` and ``above5`` do not change
    with rMax 40, 60 and 80 (within 0.01 %); the angular excess is 3.1 %.
    """
    two = _two_electron()
    for label in ROWS_390:
        values = [float(two[label, box]["above_n4"]) for box in ("59", "78", "97")]
        assert max(values) / min(values) - 1 < 1e-3, label
        far = [float(two[label, box]["above_n6"]) for box in ("59", "97")]
        assert far[0] / far[1] - 1 > 0.05, label
    hacc = {row["record"]: row for row in _rows(HACC_390)}
    for column in ("above_n4", "above_n5"):
        values = [float(hacc[name][column]) for name in HACC_390_BOXES]
        assert max(values) / min(values) - 1 < 1e-4, column
    assert math.isclose(angular_excess(), 0.031, abs_tol=5e-4)


def main() -> None:
    """Write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()
