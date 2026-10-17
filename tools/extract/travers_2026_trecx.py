# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract J. C. Travers's tRecX yields from the archived result tables.

Unlike the other scripts this one does not read the paper library: the source is a
set of calculations, archived with their inputs in ``data/external/trecx``
(``PROVENANCE.md`` there). ``ar_results.csv`` and ``he_results.csv`` hold, per peak
intensity (rows ``index`` 0 to 9), the yields of each family of runs; the numbers are
kept as the strings of those files.

Per target (Ar 3p m = 0, Ar 3p m = 1, He 1s), three files:

* ``chi_gauss_<target>``: 1 - <chi[0,20]> at the end of the Gaussian-pulse run;
* ``norm_gauss_<target>``: 1 - <norm> at the end of the same run;
* ``isurff_cos8_<target>``: the integral of the iSurff photoelectron spectrum of the
  cos8-pulse run (indices 0 to 6 only; the cells of 7 to 9 are empty).

The rows are labelled by the index of the intensity, which is also the two-digit
suffix of the run (``ar_scan_04`` is row 4). The ``printed`` cell names the run. Rows
at the numerical noise floor and the two unstable He runs are written like the others;
the source record excludes them (``excluded_rows``), with the reasons.

The 42 fs argon runs (``ar_42fs/results.csv``, TC-702) give one file,
``ionised_gauss42_Ar_m0``: the ionised fraction 1 - (ground + bound states up to
n = 10) at the end of the pulse (column ``ion_proj_le10_pe``) of the 17 base runs
(``base_*``), the row the run's intensity tag (``4p00e13``), in order of intensity.

The 6 fs runs again (``scan_6fs_v2/results.csv``, TC-812) give
``ionised_gauss6_<target>``: the same ionised fraction of the base runs
(``<family>_NN``), the row the index of the intensity as in round 1.
:func:`six_fs_uncertainty` gives their relative uncertainty by the rule of the runs'
report (the source record's ``relative_rows``).

The static rates (``static_rates/results.csv``, TC-813) give ``static_Ar_m0`` and
``static_Ar_m1``: Gamma = -2 Im E of the traced resonance (column ``Gamma``, atomic
units) at each field of the base runs, the row the field as written (``0.05329``).
:func:`static_uncertainty` gives their relative uncertainty from the direct solves of
``static_rates/results_nearest.csv``.
"""

import csv
import itertools

from _library import DATA, Point, write

DIRECTORY = DATA / "external" / "trecx"
SOURCE = "travers_2026_trecx"
# target: (material, results file, {dataset stem: (column, family of runs)})
TARGETS = {
    "Ar_m0": (
        "Ar",
        "ar_results.csv",
        {
            "chi_gauss": ("m0_gauss_chi_loss", "ar_scan"),
            "norm_gauss": ("m0_gauss_norm_loss", "ar_scan"),
            "isurff_cos8": ("m0_cos8_iSurff", "ar_isurff"),
        },
    ),
    "Ar_m1": (
        "Ar",
        "ar_results.csv",
        {
            "chi_gauss": ("m1_gauss_chi_loss", "ar_m1"),
            "norm_gauss": ("m1_gauss_norm_loss", "ar_m1"),
            "isurff_cos8": ("m1_cos8_iSurff", "ar_m1isurff"),
        },
    ),
    "He": (
        "He",
        "he_results.csv",
        {
            "chi_gauss": ("gauss_chi_loss", "he_scan"),
            "norm_gauss": ("gauss_norm_loss", "he_scan"),
            "isurff_cos8": ("cos8_iSurff", "he_isurff"),
        },
    ),
}


def _rows(name: str) -> list[dict[str, str]]:
    """Return the rows of a result table, cells as written."""
    with (DIRECTORY / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def extract() -> dict[str, list[Point]]:
    """Return the yields of every target and measure."""
    files: dict[str, list[Point]] = {}
    for target, (material, name, measures) in TARGETS.items():
        rows = _rows(name)
        for stem, (column, family) in measures.items():
            points = []
            for row in rows:
                value = row[column].strip()
                if not value:  # no run at this intensity (iSurff above index 6)
                    continue
                float(value)  # a plain number, or the extraction stops
                intensity = row["I_Wcm2"].strip()
                run = f"{family}_{int(row['index']):02d}"
                printed = f"{name} | {run} | I_Wcm2 {intensity} | {column} {value}"
                points.append(Point(row["index"], intensity, value, "", printed))
            files[f"{material}/{SOURCE}__{stem}_{target}.csv"] = points
    files[f"Ar/{SOURCE}__ionised_gauss42_Ar_m0.csv"] = _long_pulse()
    for family, (material, target) in SIX_FS.items():
        stem = f"{material}/{SOURCE}__ionised_gauss6_{target}.csv"
        files[stem] = _six_fs(family)
    for m in (0, 1):
        files[f"Ar/{SOURCE}__static_Ar_m{m}.csv"] = _static(m)
    return files | _scans() | _xe_rmt()


def _scans() -> dict[str, list[Point]]:
    """Return the ionised fractions of the scans run by A2's rule (round 2's parts)."""
    files: dict[str, list[Point]] = {}
    tables = {
        "gauss10_{}nm": WAVELENGTH_TABLES,
        "gauss30_{}nm": LONG_TABLES,
        "wallace_{}": WALLACE_TABLES,
    }
    for pattern, by_name in tables.items():
        for name, table in by_name.items():
            for family, (material, target) in WAVELENGTH.items():
                label = pattern.format(name)
                stem = f"{material}/{SOURCE}__ionised_{label}_{target}.csv"
                files[stem] = _six_fs(family, table)
    for nm, table in HYDROGEN_TABLES.items():
        files[f"H/{SOURCE}__ionised_gauss15_{nm}nm_H.csv"] = _six_fs("h_m0", table)
    for energy, table in IP_MATCHED_TABLES.items():
        for (measure, column), (family, (_, target)) in itertools.product(
            IP_MATCHED_MEASURES.items(), WAVELENGTH.items()
        ):
            stem = f"Ar/{SOURCE}__{measure}_ipmatched_{energy}_{target}.csv"
            files[stem] = _six_fs(family, table, column)
    for setup in HACC_RUNS:
        files[f"Ar/{SOURCE}__{hacc_dataset(setup)}.csv"] = _hacc(setup)
    for phase, table in WALLACE_KR_TABLES.items():
        stem = f"Kr/{SOURCE}__ionised_wallace_{phase}_Kr_m0.csv"
        files[stem] = _six_fs("kr_m0", table)
    for phase, table in WALLACE_XE_TABLES.items():
        stem = f"Xe/{SOURCE}__ionised_wallace_{phase}_Xe_m0.csv"
        files[stem] = _six_fs("xe_m0", table)
    files |= _helium_390()
    for (part, family), (material, target) in SIN2_FAMILIES.items():
        for measure, column in SIN2_MEASURES.items():
            stem = f"{material}/{SOURCE}__{measure}_sin2_{part}_{target}.csv"
            files[stem] = _six_fs(family, SIN2_TABLES[part], column)
    for family, (material, target) in KR_XE.items():
        stem = f"{material}/{SOURCE}__ionised_gauss6_{target}.csv"
        files[stem] = _six_fs(family, KR_XE_TABLE)
        files[f"{material}/{SOURCE}__static_{target}.csv"] = _kr_xe_static(material)
    return files


def _kr_xe_static(gas: str) -> list[Point]:
    """Return the static rate of Kr 4p or Xe 5p (m = 0), by field."""
    points = []
    for row in _rows(KR_XE_STATIC):
        if row["group"] != "base" or row["gas"] != gas.lower():
            continue
        if row["exit"] != "0" or row["warnings"].strip():
            msg = f"{row['name']}: exit {row['exit']}, warnings {row['warnings']!r}"
            raise ValueError(msg)
        value = row["Gamma"].strip()
        float(value)
        cells = (KR_XE_STATIC, row["name"], f"F {row['F']}", f"Gamma {value}")
        points.append(Point(row["F"], row["F"], value, "", " | ".join(cells)))
    return points


def kr_xe_static_uncertainty(gas: str) -> dict[str, float]:
    """Return the relative uncertainty of each static rate of Kr or Xe.

    The largest relative change of Gamma among the traces of the varied
    discretisations at the field (angular functions, radial elements and order, ECS
    radius 30 and 20 bohr, ECS angle 0.3 and 0.6, A2's inner grid) against the base.
    """
    rows = [row for row in _rows(KR_XE_STATIC) if row["gas"] == gas.lower()]
    base = {row["F"]: float(row["Gamma"]) for row in rows if row["group"] == "base"}
    changes = dict.fromkeys(base, 0.0)
    for row in rows:
        if row["group"] == "conv":
            change = abs(float(row["Gamma"]) / base[row["F"]] - 1)
            changes[row["F"]] = max(changes[row["F"]], change)
    return changes


LONG_PULSE = "ar_42fs/results.csv"
LONG_PULSE_COLUMN = "ion_proj_le10_pe"
LONG_PULSE_COMMIT = "a987d7c790e52b6d7f44a5615a665e0abc3e7f4a"


def _long_pulse() -> list[Point]:
    """Return the ionised fraction of the 42 fs argon scan, by intensity."""
    rows = [r for r in _rows(LONG_PULSE) if r["name"].startswith("base_")]
    points = []
    for row in sorted(rows, key=lambda r: float(r["I"])):
        commit = row["commit"].strip()
        if (
            not commit
            or not LONG_PULSE_COMMIT.startswith(commit)
            or (row["completed"] != "True")
        ):
            msg = f"{row['name']}: not a completed run of {LONG_PULSE_COMMIT[:8]}"
            raise ValueError(msg)
        value = row[LONG_PULSE_COLUMN].strip()
        float(value)
        tag = row["name"].removeprefix("base_")
        cells = (
            LONG_PULSE,
            row["name"],
            f"I {row['I']}",
            f"{LONG_PULSE_COLUMN} {value}",
        )
        printed = " | ".join(cells)
        points.append(Point(tag, row["I"], value, "", printed))
    return points


SIX_FS_TABLE = "scan_6fs_v2/results.csv"
# family of runs: (material, target of the dataset)
SIX_FS = {"ar_m0": ("Ar", "Ar_m0"), "ar_m1": ("Ar", "Ar_m1"), "he": ("He", "He")}
SIX_FS_VARIATIONS = (
    "eta40",
    "eta50",
    "nel20",
    "ord14",
    "boxp40",
    "dt23",
    "th23",
    "opthr11",
    "cut150",
    "dt12",
    "cut200",
    "th12",
    "a2grid",
)
# The 10 fs runs at 400, 800 and 1030 nm (TC-824, TC-830, TC-826): the same families
# and rule.
WAVELENGTH_TABLE = "scan_wavelength/results_400.csv"
WAVELENGTH_TABLES = {
    400: WAVELENGTH_TABLE,
    800: "scan_wavelength/results_800.csv",
    1030: "scan_wavelength/results_1030.csv",
}
WAVELENGTH = {"ar_m0": ("Ar", "Ar_m0"), "ar_m1": ("Ar", "Ar_m1")}
# The same at 400 and 1030 nm with a 30 fs pulse (TC-844).
LONG_TABLES = {
    400: "scan_wavelength/results_400_30fs.csv",
    1030: "scan_wavelength/results_1030_30fs.csv",
}
# Krypton and xenon in A2's 6 fs pulse, one electron each (TC-845): the same rule.
KR_XE_TABLE = "scan_kr_xe/results.csv"
KR_XE = {"xe_m0": ("Xe", "Xe_m0"), "kr_m0": ("Kr", "Kr_m0")}
KR_XE_STATIC = "scan_kr_xe/static_results.csv"
# One electron at haCC's ionisation energies (part J; TC-847): Tong and Lin's argon
# potential scaled to 0.560275 au (the configuration interaction's) and 0.563940 au (the
# 3-channel mMax 1 set-up's); the ionised fraction and the norm-loss at +1000 au.
IP_MATCHED_TABLES = {
    "ci": "ip_matched/results_ci.csv",
    "c3m1": "ip_matched/results_c3m1.csv",
    # E3 (TC-861): at the 3-channel mMax 0 set-up's own 0.562699 au, A2's 6 fs pulse
    # and the 20 fs one.
    "c3": "ip_matched/results_c3_6fs.csv",
    "c3_20fs": "ip_matched/results_c3_20fs.csv",
}
# A set without convergence runs of its own takes those of the set with the same pulse
# and discretisation (the ionisation energies 0.2 % apart).
IP_MATCHED_CONVERGENCE = {"c3": "c3m1"}
IP_MATCHED_MEASURES = {"ionised": LONG_PULSE_COLUMN, "normloss": "norm_loss_end"}
# haCC's argon atom (part E; TC-847): the norm-loss at +1000 au, by set-up (number of
# ion channels; mMax 0), at E2's intensities, with A2's 6 fs pulse.
HACC_TABLE = "ar_hacc/results.csv"
# The runs of each set-up by intensity (the run's own label), and the one-channel base
# of each intensity, against which the variations of the discretisation are taken. At
# 9.6e13 the one- and four-channel runs are the MPI checks' (complete, the full window).
HACC_RUNS = {
    "c4": {
        "4p03e13": "e2_4p03e13_c4",
        "9p61e13": "cost_l20_c4_9p61e13",
        "2p29e14": "e2_2p29e14_c4",
    },
    "c3": {"9p61e13": "e2_9p61e13_c3"},
    "c3m1": {"9p61e13": "e2_9p61e13_c3m1"},
    "c3_20fs": {"9p61e13": "e3_9p61e13_c3_20fs"},
}
# Each set-up's pulse (E3, the 20 fs Gaussian: TC-861) and its own lMax 30 run.
HACC_PULSE = {"c4": "gauss6", "c3": "gauss6", "c3m1": "gauss6", "c3_20fs": "gauss20"}
HACC_LMAX30 = {
    "e2_9p61e13_c3": "e2_9p61e13_c3l30",
    "e3_9p61e13_c3_20fs": "e3_9p61e13_c3l30_20fs",
}
HACC_BASE = {
    "4p03e13": "e2_4p03e13_base",
    "9p61e13": "cost_l20_c1_9p61e13",
    "2p29e14": "e2_2p29e14_base",
}
# The variations of the discretisation (one channel), and the lMax 30 and 40 runs.
HACC_VARIATIONS = ("acc8", "r120", "d025", "d04", "l30", "l40")
# Atomic hydrogen, a 15 fs Gaussian at 400 and 800 nm (TC-843): the same rule.
HYDROGEN_TABLES = {
    400: "scan_hydrogen/results_400.csv",
    800: "scan_hydrogen/results_800.csv",
}
# Wallace et al.'s pulse, 5.5 fs, at two carrier-envelope phases (TC-836): the same
# families; the pi/2 runs have no convergence runs and take phase 0's.
WALLACE_TABLES = {
    "cep0": "scan_wallace/results_cep0.csv",
    "cep90": "scan_wallace/results_cep90.csv",
}
# One electron in the sin^2 pulses of published many-electron calculations (part M;
# TC-858): E(t) = E0 sin(omega t) sin^2(pi t / tau) on the field, read at tau, the end
# of the pulse; the probability outside 20 au (1 - <chi[0,20]>, their measure) and the
# ionised fraction.
SIN2_TABLES = {
    "m1": "scan_sin2/results_m1_ar.csv",
    "m1mk": "scan_sin2/results_m1_mk.csv",
    "m2": "scan_sin2/results_m2_kr.csv",
    "m3": "scan_sin2/results_m3_ar.csv",
    "m4": "scan_sin2/results_m4_he.csv",
}
SIN2_FAMILIES = {
    ("m1", "ar_m0"): ("Ar", "Ar_m0"),
    ("m1", "ar_m1"): ("Ar", "Ar_m1"),
    ("m1mk", "armk_m0"): ("Ar", "ArMK_m0"),
    ("m2", "kr_m0"): ("Kr", "Kr_m0"),
    ("m2", "kr_m1"): ("Kr", "Kr_m1"),
    ("m3", "ar_m0"): ("Ar", "Ar_m0"),
    ("m3", "ar_m1"): ("Ar", "Ar_m1"),
    ("m4", "he"): ("He", "He"),
}
SIN2_MEASURES = {"outside20": "chi20_loss_pe", "ionised": LONG_PULSE_COLUMN}
# Helium with two active electrons (part K, step 2; TC-1050): A2's 6 fs Gaussian at 390
# nm, the norm-loss at +1000 au, by set-up (the base basis, BOX 60, |m| <= 2) and
# intensity; the base's variations at 5.69e14; haCC's helium (He19, one channel,
# lMax 40, rMax 40) and its variations; one electron at the bases' ionisation energies.
HE_TWO_ELECTRON_TABLE = "he_two_electron/results_2e.csv"
HE_TWO_ELECTRON_RUNS = {
    "base": {"5p69e14": "k2_5p69e14", "3p68e14": "k2_3p68e14"},
    "box60": {"5p69e14": "k2_5p69e14_box60"},
    "m5": {"5p69e14": "k2_5p69e14_m5ls4"},
}
HE_TWO_ELECTRON_BASE = "k2_5p69e14"
HE_TWO_ELECTRON_VARIATIONS = (
    "k2_5p69e14_dt23",
    "k2_5p69e14_A20",
    "k2_5p69e14_m5ls4",
    "k2_5p69e14_box60",
)
HE_HACC_TABLE = "he_two_electron/results_hacc.csv"
HE_HACC_RUNS = {"5p69e14": "h2_5p69e14_r40", "3p68e14": "h2_3p68e14_r40"}
# (base, variation): rMax 60, lMax 50, and the n = 2 ion channels with mMax 1, the last
# only at lMax 20 (out of memory at 40).
HE_HACC_VARIATIONS = (
    ("h2_5p69e14_r40", "h2_5p69e14_r60"),
    ("h2_5p69e14_r40", "h2_5p69e14_r40_l50"),
    ("h2_5p69e14_r40_l20", "h2_5p69e14_r40_l20_c5m1"),
)
HE_ONE_ELECTRON_TABLES = {
    "base": "he_two_electron/results_sae_base.csv",
    "m5": "he_two_electron/results_sae_m5.csv",
    "tl": "he_two_electron/results_sae_tl.csv",
}
# Krypton in the same pulse, m = 0, with part G's potential and grid (TC-846).
WALLACE_KR_TABLES = {
    "cep0": "scan_wallace/results_kr_cep0.csv",
    "cep90": "scan_wallace/results_kr_cep90.csv",
}
# Xenon in the same pulse, m = 0, with part G's potential and grid (TC-849).
WALLACE_XE_TABLES = {
    "cep0": "scan_wallace/results_xe_cep0.csv",
    "cep90": "scan_wallace/results_xe_cep90.csv",
}


XE_RMT_TABLE = "xe_rmt/results.csv"
XE_RMT_COLUMN = "P_ground_3000"  # the 5p population at 3000 au, RMT's end of run
XE_RMT_SETS = {"": "base", "_unscaled": "unscaled"}  # dataset suffix: group of runs
XE_RMT_ENERGIES = {"base": -0.444, "unscaled": -0.4448128882}  # au, the 5p
XE_RMT_SETTLED = 1.1e-5  # P at 1500 au against 3000 au, relative (at most 1.005e-5)
# The convergence runs at 1.5 and 2.5e14 W/cm2 (m = 0); A2's inner grid moves the 5p
# energy (to -0.443953 au), so it is a check of the energy, not of the convergence.
XE_RMT_CONVERGENCE = ("1.5", "2.5")
XE_RMT_NOT_CONVERGENCE = ("_a2grid",)


def _six_fs_runs(table: str = SIX_FS_TABLE) -> dict[str, dict[str, str]]:
    """Return the completed runs of a scan's table by name."""
    rows = _rows(table)
    if any(row["completed"] != "True" for row in rows):
        msg = f"{table}: a run did not complete"
        raise ValueError(msg)
    return {row["name"]: row for row in rows}


def _six_fs(
    family: str, table: str = SIX_FS_TABLE, column: str = LONG_PULSE_COLUMN
) -> list[Point]:
    """Return a measure (the ionised fraction) of one family of runs, by intensity."""
    base = sorted(
        (
            row
            for row in _six_fs_runs(table).values()
            if row["family"] == family and row["group"] == "base"
        ),
        key=lambda row: int(row["row"]),
    )
    points = []
    for row in base:
        value = row[column].strip()
        float(value)
        cells = (table, row["name"], f"I {row['I']}", f"{column} {value}")
        points.append(Point(row["row"], row["I"], value, "", " | ".join(cells)))
    return points


def six_fs_uncertainty(
    family: str, table: str = SIX_FS_TABLE, column: str = LONG_PULSE_COLUMN
) -> dict[str, float | None]:
    """Return the relative uncertainty of each row of a family of a scan's runs.

    The rule of the runs' report (its ``tables.py``): the largest relative change of
    the ionised fraction among the convergence runs of the nearest convergence rows
    below and above, counting only convergence rows above the floor (10 times the loss
    of the run without a field), and at least that loss over the yield. None below the
    floor. A family with no convergence runs of its own (m = 1 at 400 nm) takes the same
    atom's m = 0 changes at the same row. ``column`` is the measure (the ionised
    fraction, the probability outside 20 au, or the norm-loss of the runs at haCC's
    ionisation energies).
    """
    runs = _six_fs_runs(table)

    def value(name: str) -> float:
        return float(runs[name][column])

    drift = float(runs[f"{family}_nofield"]["norm_loss_end"])
    base = {
        int(row["row"]): row["name"]
        for row in runs.values()
        if row["family"] == family and row["group"] == "base"
    }
    convergence = sorted(
        index
        for index, name in base.items()
        if f"{name}_eta40" in runs and value(name) > 10 * drift
    )
    if not convergence and family.endswith("_m1"):
        reference = six_fs_uncertainty(
            family.removesuffix("_m1") + "_m0", table, column
        )
        return {
            row: None
            if not value(base[int(row)]) > 10 * drift
            else max(reference[row] or 0.0, drift / value(base[int(row)]))
            for row in map(str, sorted(base))
        }
    uncertainty: dict[str, float | None] = {}
    for index, name in sorted(base.items()):
        if not value(name) > 10 * drift:
            uncertainty[str(index)] = None
            continue
        below = [c for c in convergence if c <= index][-1:]
        above = [c for c in convergence if c >= index][:1]
        changes = [
            abs(value(f"{base[c]}_{v}") / value(base[c]) - 1)
            for c in {*below, *above}
            for v in SIX_FS_VARIATIONS
            if f"{base[c]}_{v}" in runs
        ]
        uncertainty[str(index)] = max([*changes, drift / value(name)])
    return uncertainty


def wallace_uncertainty(
    family: str, phase: str, tables: dict[str, str] | None = None
) -> dict[str, float | None]:
    """Return the relative uncertainty of each row of the Wallace-pulse runs.

    Phase 0 by :func:`six_fs_uncertainty`; phase pi/2, run without convergence runs of
    its own on the same discretisation, takes phase 0's at the same row, and at least
    its own loss without a field over the yield. ``tables`` are the phases' tables,
    argon's (:data:`WALLACE_TABLES`) by default.
    """
    tables = WALLACE_TABLES if tables is None else tables
    reference = six_fs_uncertainty(family, tables["cep0"])
    if phase == "cep0":
        return reference
    runs = _six_fs_runs(tables[phase])
    drift = float(runs[f"{family}_nofield"]["norm_loss_end"])
    result: dict[str, float | None] = {}
    for row, value in reference.items():
        name = f"{family}_{int(row):02d}"
        own = float(runs[name][LONG_PULSE_COLUMN])
        result[row] = (
            None if value is None or not own > 10 * drift else max(value, drift / own)
        )
    return result


def xe_rmt_dataset(m: int, suffix: str) -> str:
    """Return the dataset id of part N's xenon 5p survival for one m and potential."""
    return f"ground_survival_rmtpulse_m{m}{suffix}_Xe"


def _xe_rmt_rows() -> list[dict[str, str]]:
    """Return part N's rows, after checking each run is done and settled."""
    rows = _rows(XE_RMT_TABLE)
    for row in rows:
        if row["done"] != "True" or row["exit"] != "0":
            msg = f"{row['name']}: not complete"
            raise ValueError(msg)
        if row["group"] in XE_RMT_ENERGIES:
            energy = XE_RMT_ENERGIES[row["group"]]
            if abs(float(row["E_5p"]) - energy) > 1e-9:
                msg = f"{row['name']}: the 5p at {row['E_5p']} au, not {energy}"
                raise ValueError(msg)
            if abs(float(row["settle_1500_3000_rel"])) > XE_RMT_SETTLED:
                msg = f"{row['name']}: not settled by 1500 au"
                raise ValueError(msg)
    return rows


def _xe_rmt() -> dict[str, list[Point]]:
    """Return part N's 5p survival at 3000 au, by m and potential, rows by label."""
    rows = _xe_rmt_rows()
    files = {}
    for suffix, group in XE_RMT_SETS.items():
        for m in (0, 1):
            points = []
            for row in rows:
                if row["group"] != group or row["m"] != str(m):
                    continue
                value = row[XE_RMT_COLUMN].strip()
                float(value)
                cells = (XE_RMT_TABLE, row["name"], f"label {row['label']}")
                printed = " | ".join((*cells, f"{XE_RMT_COLUMN} {value}"))
                x = f"{row['label']}e14"
                points.append(Point(row["label"], x, value, "", printed))
            files[f"Xe/{SOURCE}__{xe_rmt_dataset(m, suffix)}.csv"] = points
    return files


def xe_rmt_uncertainty(m: int, suffix: str) -> dict[str, float]:
    """Return the relative uncertainty of part N's survival P, by row.

    The largest relative change of the depletion 1 - P among the convergence runs at
    the convergence intensity nearer the row's (:data:`XE_RMT_CONVERGENCE`, m = 0),
    carried to P: u_D (1 - P) / P. The m = 1 runs and the unscaled potential take the
    same u_D.
    """
    rows = _xe_rmt_rows()
    change = {}
    for label in XE_RMT_CONVERGENCE:
        base = next(
            r
            for r in rows
            if r["group"] == "base" and r["m"] == "0" and r["label"] == label
        )
        depletion = 1 - float(base[XE_RMT_COLUMN])
        change[label] = max(
            abs((1 - float(r[XE_RMT_COLUMN])) / depletion - 1)
            for r in rows
            if r["group"] == "conv"
            and r["label"] == label
            and not r["name"].endswith(XE_RMT_NOT_CONVERGENCE)
        )
    out = {}
    for row in rows:
        if row["group"] != XE_RMT_SETS[suffix] or row["m"] != str(m):
            continue
        nearest = min(
            XE_RMT_CONVERGENCE, key=lambda c: abs(float(c) - float(row["label"]))
        )
        survival = float(row[XE_RMT_COLUMN])
        out[row["label"]] = change[nearest] * (1 - survival) / survival
    return out


def he_two_electron_dataset(setup: str) -> str:
    """Return the dataset id of one of the two-electron helium set-ups."""
    suffix = "" if setup == "base" else f"_{setup}"
    return f"twoelectron_normloss_gauss6_390nm_He{suffix}"


def _he_table(table: str, runs: dict[str, str]) -> list[Point]:
    """Return the norm-loss at +1000 au of helium runs, by intensity label."""
    rows = {row["name"]: row for row in _rows(table)}
    points = []
    for label, name in runs.items():
        row = rows[name]
        if row["done"] != "True":
            msg = f"{name}: not complete"
            raise ValueError(msg)
        value = row["norm_loss_end"].strip()
        float(value)
        cells = (table, name, f"I_Wcm2 {row['I_Wcm2']}", f"norm_loss_end {value}")
        points.append(Point(label, row["I_Wcm2"], value, "", " | ".join(cells)))
    return points


def _helium_390() -> dict[str, list[Point]]:
    """Return part K's helium yields at 390 nm: two electrons, haCC, one electron."""
    files = {
        f"He/{SOURCE}__{he_two_electron_dataset(setup)}.csv": _he_table(
            HE_TWO_ELECTRON_TABLE, runs
        )
        for setup, runs in HE_TWO_ELECTRON_RUNS.items()
    }
    files[f"He/{SOURCE}__hacc_normloss_gauss6_390nm_He.csv"] = _he_table(
        HE_HACC_TABLE, HE_HACC_RUNS
    )
    for energy, table in HE_ONE_ELECTRON_TABLES.items():
        for measure, column in IP_MATCHED_MEASURES.items():
            stem = f"He/{SOURCE}__{measure}_ipmatched_390nm_{energy}_He.csv"
            files[stem] = _six_fs("he", table, column)
    return files


def he_two_electron_uncertainty() -> float:
    """Return the relative uncertainty of the two-electron helium yields.

    The largest change at +1000 au among the base's variations at 5.69e14 W/cm2 (step
    x 2/3, l <= 19, |m| <= 2, BOX 60), applied to every set-up and intensity: the box,
    -5.9 %, which BOX 80 (not run) would settle.
    """
    rows = {row["name"]: row for row in _rows(HE_TWO_ELECTRON_TABLE)}
    base = float(rows[HE_TWO_ELECTRON_BASE]["norm_loss_end"])
    return max(
        abs(float(rows[name]["norm_loss_end"]) / base - 1)
        for name in HE_TWO_ELECTRON_VARIATIONS
    )


def he_hacc_uncertainty() -> float:
    """Return the relative uncertainty of haCC's helium yields at 390 nm.

    The largest change at +1000 au among rMax 60, lMax 50 and the n = 2 ion channels
    (with mMax 1, against the one-channel run at the same lMax 20).
    """
    rows = {row["name"]: row for row in _rows(HE_HACC_TABLE)}
    changes = []
    for old, new in HE_HACC_VARIATIONS:
        if rows[new]["done"] != "True":
            msg = f"{new}: not complete in {HE_HACC_TABLE}"
            raise ValueError(msg)
        old_value = float(rows[old]["norm_loss_end"])
        changes.append(abs(float(rows[new]["norm_loss_end"]) / old_value - 1))
    return max(changes)


def he_one_electron_uncertainty(energy: str, column: str) -> dict[str, float | None]:
    """Return the uncertainty of one electron of helium at 390 nm.

    As :func:`six_fs_uncertainty`; the sets without convergence runs of their own (the
    |m| <= 2 basis's and Tong and Lin's ionisation energies) take the base's.
    """
    base = six_fs_uncertainty("he", HE_ONE_ELECTRON_TABLES["base"], column)
    own = six_fs_uncertainty("he", HE_ONE_ELECTRON_TABLES[energy], column)
    return {
        row: None if value is None else max(value, base[row] or 0.0)
        for row, value in own.items()
    }


def hacc_dataset(setup: str) -> str:
    """Return the dataset id of one of haCC's set-ups."""
    return f"hacc_normloss_{HACC_PULSE[setup]}_Ar_{setup.removesuffix('_20fs')}"


def ip_matched_uncertainty(
    energy: str, family: str, column: str
) -> dict[str, float | None]:
    """Return the uncertainty of one electron at one of haCC's ionisation energies.

    As :func:`six_fs_uncertainty`; a set without convergence runs of its own takes, row
    by row, those of :data:`IP_MATCHED_CONVERGENCE`'s set, and at least its own loss
    without a field over the yield.
    """
    own = six_fs_uncertainty(family, IP_MATCHED_TABLES[energy], column)
    if energy not in IP_MATCHED_CONVERGENCE:
        return own
    other = six_fs_uncertainty(
        family, IP_MATCHED_TABLES[IP_MATCHED_CONVERGENCE[energy]], column
    )
    return {
        row: None if value is None else max(value, other[row] or 0.0)
        for row, value in own.items()
    }


def _hacc_runs() -> dict[str, dict[str, str]]:
    """Return haCC's argon runs by name."""
    return {row["name"]: row for row in _rows(HACC_TABLE)}


def _hacc(setup: str) -> list[Point]:
    """Return haCC's argon norm-loss at +1000 au of one set-up, by intensity."""
    runs = _hacc_runs()
    points = []
    for label, name in HACC_RUNS[setup].items():
        row = runs[name]
        if row["done"] != "True":
            msg = f"{name}: not complete"
            raise ValueError(msg)
        value = row["norm_loss_end"].strip()
        float(value)
        cells = (HACC_TABLE, name, f"I_Wcm2 {row['I_Wcm2']}", f"norm_loss_end {value}")
        points.append(Point(label, row["I_Wcm2"], value, "", " | ".join(cells)))
    return points


def hacc_uncertainty(setup: str) -> dict[str, float]:
    """Return the relative uncertainty of each of haCC's argon yields of one set-up.

    The largest relative change among the one-channel variations at the same
    intensity (accuracy, rMax 120, radial order 19 and 12, lMax 30 and 40), each against
    the one-channel base, and the set-up's own lMax 30 run where there is one. The
    angular momentum dominates: lMax 20 is high by 1 to 7 %.
    """
    runs = _hacc_runs()
    uncertainty: dict[str, float] = {}
    for label, name in HACC_RUNS[setup].items():
        base = float(runs[HACC_BASE[label]]["norm_loss_end"])
        changes = [
            abs(float(runs[f"e2_{label}_{v}"]["norm_loss_end"]) / base - 1)
            for v in HACC_VARIATIONS
            if f"e2_{label}_{v}" in runs
        ]
        lmax30 = HACC_LMAX30.get(name, f"{name}l30")
        if lmax30 in runs:
            own = float(runs[name]["norm_loss_end"])
            changes.append(abs(float(runs[lmax30]["norm_loss_end"]) / own - 1))
        uncertainty[label] = max(changes)
    return uncertainty


STATIC_TABLE = "static_rates/results.csv"
STATIC_NEAREST = "static_rates/results_nearest.csv"
# The ECS radius moved outwards loses the resonance (the runs' report): not a variation.
STATIC_LEFT_OUT = ("box60",)


def _static(m: int) -> list[Point]:
    """Return the static rate of Ar 3p with magnetic quantum number m, by field."""
    points = []
    for row in _rows(STATIC_TABLE):
        if row["group"] != "base" or row["m"] != str(m):
            continue
        if row["exit"] != "0" or row["warnings"].strip():
            msg = f"{row['name']}: exit {row['exit']}, warnings {row['warnings']!r}"
            raise ValueError(msg)
        value = row["Gamma"].strip()
        float(value)
        cells = (STATIC_TABLE, row["name"], f"F {row['F']}", f"Gamma {value}")
        points.append(Point(row["F"], row["F"], value, "", " | ".join(cells)))
    return points


def static_uncertainty(m: int) -> dict[str, float]:
    """Return the relative uncertainty of each static rate of Ar 3p, m.

    The largest relative change of Gamma among the direct solves of the varied
    discretisations at the field (angular functions, radial elements and order, ECS
    angle 0.3 and 0.6, ECS radius 30 and 20 bohr) against the traced base value.
    """
    base = {
        row["name"]: float(row["Gamma"])
        for row in _rows(STATIC_TABLE)
        if row["group"] == "base" and row["m"] == str(m)
    }
    fields = _base_fields(m)
    changes: dict[str, float] = {}
    for row in _rows(STATIC_NEAREST):
        name, _, variation = row["name"].rpartition("_")
        if name not in base or variation in STATIC_LEFT_OUT:
            continue
        change = abs(float(row["Gamma"]) / base[name] - 1)
        changes[fields[name]] = max(changes.get(fields[name], 0.0), change)
    return changes


def _base_fields(m: int) -> dict[str, str]:
    """Return the field, as written, of each base static run of m."""
    return {
        row["name"]: row["F"]
        for row in _rows(STATIC_TABLE)
        if row["group"] == "base" and row["m"] == str(m)
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
