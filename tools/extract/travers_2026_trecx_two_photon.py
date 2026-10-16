# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract J. C. Travers's tRecX two-photon cross sections from the archived tables.

Like ``travers_2026_trecx.py`` this script does not read the paper library: the source
is a set of calculations (the two-photon TDSE campaign, TC-1048, run on a Slurm cluster
on 9-11 October 2026), archived with their inputs and scripts in
``data/external/trecx/two_photon`` (``PROVENANCE.md`` there). Three files, each the
degenerate generalised two-photon ionisation cross section sigma2 = P / int Phi^2 dt at
1e12 W/cm2 in a 20 fs Gaussian (the full width at half maximum of the intensity), Phi =
I(t) / (hbar omega), I(t) = I0 exp(-4 ln2 t^2 / T^2), so that int Phi^2 dt = (I0 / hbar
omega)^2 T sqrt(pi / (8 ln2)) (the runs' ``extract_run1.py``):

* ``sigma2_hacc_He``: helium in haCC (He19, one channel), run 1; P the two-photon peak
  of the iSurff spectrum (``results_hacc.csv``, column ``sigma2_peak2_cm4s``, kept as
  the table's string);
* ``sigma2_sae_He``: one-electron helium (Tong and Lin's potential), run 1; P the
  projection 1 - (bound states) at the end of the window, ``P_ion_proj_end``, less
  the same column of the no-field run of the same discretisation (the drift);
* ``sigma2_sae_Ar``: one-electron argon (Tong and Lin's 3p), run 3; P = 2 P(m = 0) +
  4 P(m = 1), each the two-photon peak of its run's spectrum (``P_p11``), the
  combination over the six 3p electrons agreed for the runs.

The rows are the photon energy as the tables write it (``13.0``, ``8p5`` becomes
``8.5``). :func:`uncertainties` gives each row's relative uncertainty by the rule the
record states; :func:`checks` compares every value with the numbers the runs'
summaries print and raises on a miss.
"""

import csv
import math
import re

from _library import DATA, Point, write

DIRECTORY = DATA / "external" / "trecx" / "two_photon"
SOURCE = "travers_2026_trecx_two_photon"
EV_J = 1.602176634e-19  # J per eV (exact, SI 2019)
PULSE_FS = 20.0
INTENSITY = "1000000000000.0"  # the 1e12 W/cm2 rows, as the tables write them
# The summaries' printed sigma2, in 1e-52 cm4 s (run 1, summary_run1.md and REPORT.md)
# and 1e-50 cm4 s (run 3, summary_runs2to6.md); rows not printed are absent.
PRINTED_HACC = {
    "13.0": "11.890",
    "14.0": "12.804",
    "16.0": "10.723",
    "18.0": "9.063",
    "20.0": "10.175",
}
PRINTED_SAE_HE = {
    "13.0": "7.456",
    "14.0": "7.328",
    "16.0": "6.524",
    "18.0": "6.228",
    "20.0": "9.346",
}
PRINTED_SAE_AR = {
    "8.0": "1.921",
    "8.5": "1.916",
    "9.0": "2.094",
    "10.0": "2.503",
    "11.0": "3.052",
}
# Run 4's lowest-order amplitudes as the runs' extract_runs.py holds them (au).
PRINTED_M = {"0": "8.914061", "1": "2.418379"}


def _rows(name: str) -> list[dict[str, str]]:
    """Return the rows of an archived table, cells as written."""
    with (DIRECTORY / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def phi2_integral(intensity: float, photon_ev: float) -> float:
    """Return int Phi^2 dt (cm^-4 s^-1) of the 20 fs Gaussian, as the runs define it."""
    flux = intensity / (photon_ev * EV_J)
    return flux**2 * PULSE_FS * 1e-15 * math.sqrt(math.pi / (8 * math.log(2)))


def _run1(model: str) -> dict[str, dict[str, str]]:
    """Return run 1's rows of a model by name."""
    return {row["name"]: row for row in _rows(f"run1/results_{model}.csv")}


def _hacc() -> list[Point]:
    rows = _run1("hacc")
    points = []
    for row in rows.values():
        if row["variant"] != "base" or row["I_Wcm2"] != INTENSITY:
            continue
        value = row["sigma2_peak2_cm4s"].strip()
        recomputed = float(row["P_peak2"]) / phi2_integral(
            float(row["I_Wcm2"]), float(row["E_ph_eV"])
        )
        if not math.isclose(float(value), recomputed, rel_tol=1e-12):
            msg = f"{row['name']}: sigma2 {value} against {recomputed:.6e}"
            raise ValueError(msg)
        printed = f"run1/results_hacc.csv | {row['name']} | sigma2_peak2_cm4s {value}"
        points.append(Point(row["E_ph_eV"], row["E_ph_eV"], value, "", printed))
    return points


def _sae_he() -> list[Point]:
    rows = _run1("sae")
    drift = rows["he_nofield"]["P_ion_proj_end"]
    points = []
    for row in rows.values():
        if row["variant"] != "base" or row["I_Wcm2"] != INTENSITY:
            continue
        ionised = float(row["P_ion_proj_end"]) - float(drift)
        sigma = ionised / phi2_integral(float(row["I_Wcm2"]), float(row["E_ph_eV"]))
        printed = (
            f"run1/results_sae.csv | {row['name']} | P_ion_proj_end "
            f"{row['P_ion_proj_end']} - he_nofield P_ion_proj_end {drift}"
        )
        energy = row["E_ph_eV"]
        points.append(Point(energy, energy, f"{sigma:.6e}", "", printed))
    return points


def _ar_degenerate() -> dict[str, tuple[dict[str, str], dict[str, str]]]:
    """Return run 3's degenerate 1e12 base runs by photon energy: (m = 0, m = 1)."""
    rows = {row["name"]: row for row in _rows("run3/results.csv")}
    pattern = re.compile(r"ar_m0_r3_(\d+(?:p\d+)?)eV_I1e12$")
    out = {}
    for name, row in rows.items():
        found = pattern.match(name)
        if found:
            energy = str(float(found.group(1).replace("p", ".")))
            out[energy] = (row, rows[name.replace("_m0_", "_m1_")])
    return dict(sorted(out.items(), key=lambda item: float(item[0])))


def _sae_ar() -> list[Point]:
    points = []
    for energy, (m0, m1) in _ar_degenerate().items():
        for row in (m0, m1):
            if row["I_per_colour"] != INTENSITY or row["fwhm_fs"] != "20.0":
                msg = f"{row['name']}: not a 1e12 W/cm2, 20 fs run"
                raise ValueError(msg)
        atom = 2 * float(m0["P_p11"]) + 4 * float(m1["P_p11"])
        sigma = atom / phi2_integral(float(INTENSITY), float(energy))
        printed = (
            f"run3/results.csv | 2 x {m0['name']} P_p11 {m0['P_p11']} + 4 x "
            f"{m1['name']} P_p11 {m1['P_p11']}"
        )
        points.append(Point(energy, energy, f"{sigma:.6e}", "", printed))
    return points


def extract() -> dict[str, list[Point]]:
    """Return the three cross sections."""
    return {
        f"He/{SOURCE}__sigma2_hacc_He.csv": _hacc(),
        f"He/{SOURCE}__sigma2_sae_He.csv": _sae_he(),
        f"Ar/{SOURCE}__sigma2_sae_Ar.csv": _sae_ar(),
    }


def _change(value: float, base: float) -> float:
    return abs(value / base - 1)


def uncertainties() -> dict[str, dict[str, float]]:
    """Return each dataset's relative uncertainty per row, by the record's rule.

    The larger of (a) the largest change of the convergence runs of the family (one
    photon energy, 1e12 W/cm2, each one-electron variant less its own no-field drift)
    and (b) the row's disagreement between two readings of the same run's ionisation:
    haCC's spectrum peak against its norm-loss; one-electron helium's projection
    against the spectrum peak of the BOX 80 and BOX 120 runs; one-electron argon's
    spectrum peak against its drift-corrected projection (both combined as the atom).

    Returns
    -------
    dict
        Dataset id to {row: relative uncertainty}.
    """
    hacc = _run1("hacc")
    base = float(hacc["he_16eV_I1e12"]["sigma2_peak2_cm4s"])
    convergence = max(
        _change(float(hacc[f"he_16eV_I1e12_{v}"]["sigma2_peak2_cm4s"]), base)
        for v in ("l15", "r60")
    )
    out: dict[str, dict[str, float]] = {"sigma2_hacc_He": {}}
    for row in hacc.values():
        if row["variant"] == "base" and row["I_Wcm2"] == INTENSITY:
            spread = _change(float(row["P_peak2"]), float(row["norm_loss_end"]))
            out["sigma2_hacc_He"][row["E_ph_eV"]] = max(convergence, spread)
    sae = _run1("sae")

    def ionised(name: str, variant: str = "") -> float:
        drift = sae[f"he_nofield_{variant}" if variant else "he_nofield"]
        return float(sae[name]["P_ion_proj_end"]) - float(drift["P_ion_proj_end"])

    base = ionised("he_16eV_I1e12")
    convergence = max(
        _change(ionised(f"he_16eV_I1e12_{v}", v), base)
        for v in ("eta15", "th23", "dt23", "cut150", "box80", "box120")
    )
    out["sigma2_sae_He"] = {}
    for row in sae.values():
        if row["variant"] == "base" and row["I_Wcm2"] == INTENSITY:
            name = row["name"]
            spread = max(
                _change(float(sae[f"{name}_{box}"]["P_peak2"]), ionised(name))
                for box in ("box80", "box120")
            )
            out["sigma2_sae_He"][row["E_ph_eV"]] = max(convergence, spread)
    ar = {row["name"]: row for row in _rows("run3/results.csv")}
    base = float(ar["ar_m0_r3_9eV_I1e12"]["P_p11"])
    convergence = max(
        _change(float(ar[f"ar_m0_r3_9eV_I1e12_{v}"]["P_p11"]), base)
        for v in ("eta15", "box120", "dt23", "cut150", "th23")
    )
    out["sigma2_sae_Ar"] = {}
    for energy, (m0, m1) in _ar_degenerate().items():
        peak = 2 * float(m0["P_p11"]) + 4 * float(m1["P_p11"])
        projection = 2 * float(m0["P_proj_corr"]) + 4 * float(m1["P_proj_corr"])
        out["sigma2_sae_Ar"][energy] = max(convergence, _change(peak, projection))
    return out


def checks() -> dict[str, object]:
    """Compare every value with the runs' printed summaries; raise on a miss.

    Returns
    -------
    dict
        What was compared, and the largest relative uncertainty of each dataset.
    """
    files = extract()
    failures = []
    expected = {
        f"He/{SOURCE}__sigma2_hacc_He.csv": (PRINTED_HACC, 1e52),
        f"He/{SOURCE}__sigma2_sae_He.csv": (PRINTED_SAE_HE, 1e52),
        f"Ar/{SOURCE}__sigma2_sae_Ar.csv": (PRINTED_SAE_AR, 1e50),
    }
    compared = 0
    for relative, (printed, factor) in expected.items():
        by_row = {p.row: p for p in files[relative]}
        if not set(printed) <= set(by_row):
            failures.append(f"{relative}: rows {sorted(set(printed) - set(by_row))}")
        for row, text in printed.items():
            digits = len(text.split(".")[1])
            got = f"{float(by_row[row].value) * factor:.{digits}f}"
            compared += 1
            if got != text:
                failures.append(f"{relative} {row}: {got} against printed {text}")
    source = (DIRECTORY / "extract_runs.py.txt").read_text(encoding="utf-8")
    held = re.search(r"M_PT = \{0: ([\d.]+), 1: ([\d.]+)\}", source)
    if held is None or (held.group(1), held.group(2)) != tuple(PRINTED_M.values()):
        failures.append("run 4's amplitudes are not those of extract_runs.py")
    for row in _rows("run4/results.csv"):
        ratio = float(row["P_4p_TDSE_over_PT"])
        limit = 0.002 if row["I_per_colour"] == "10000000000.0" else 0.015
        if abs(ratio - 1) > limit:
            failures.append(f"{row['name']}: TDSE / PT {ratio}")
    if failures:
        raise ValueError("\n".join(failures))
    worst = {k: max(v.values()) for k, v in uncertainties().items()}
    return {"printed values compared": compared, "largest uncertainty": worst}


def main() -> None:
    """Write the data files and report the checks and uncertainties."""
    write(extract())
    print(checks())
    for dataset, rows in uncertainties().items():
        print(dataset, {row: round(u, 5) for row, u in rows.items()})


if __name__ == "__main__":
    main()
