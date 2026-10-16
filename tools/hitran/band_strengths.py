# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Sum HITRAN line intensities into vibrational band strengths.

    uv run python tools/hitran/band_strengths.py CH4 N2O SF6

For a gas with a HITRAN line list (CH4, N2O), downloads every line of every
isotopologue between 0 and ``NU_MAX`` cm^-1 through the HITRAN API (version 2), keeps
the raw ``.par`` records in a cache outside the repository
(``~/.cache/optmat/hitran``), groups the lines by isotopologue and
vibrational band (the global upper and lower quanta) and writes
``data/external/hitran/<gas>_bands.csv``: per band the number of lines, the
intensity-weighted wavenumber, the range, the summed intensity ``S`` (cm/molecule at
296 K, per molecule of the natural isotopic mixture, as HITRAN defines it), the sum of
``S / nu^2`` (cm^3/molecule) and the HITRAN reference codes of the intensities. A gas
with cross-sections only (SF6) is integrated from the one file that covers every band
(``CROSS_SECTIONS``) over the windows given there, into the same columns, and the
strongest window is checked against an independent measurement. Only these derived
totals are committed; see ``data/external/hitran/PROVENANCE.md`` for the policy and
the citations.

The API key is read from ``$HITRAN_API_KEY``, or from the file named by
``$HITRAN_API_KEY_FILE``. It is never printed or written.
"""

import csv
import itertools
import json
import os
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

HOST = "https://hitran.org"
NU_MAX = 11000.0  # cm^-1: the fundamentals and the main combination bands and overtones
MOLECULES = {"CH4": 6, "N2O": 4, "SF6": 30, "CO2": 2, "O3": 3, "CO": 5}
CACHE = Path("~/.cache/optmat/hitran").expanduser()
# SF6 has no line list in HITRAN: Ref. [1164] (PNNL, 298.1 K, N2-broadened, 0.112 cm^-1)
# covers 560-6500 cm^-1, every band; Ref. [1638] (296.1 K, air) is the check. The
# windows separate the bands at the minima of the spectrum.
CROSS_SECTIONS = {
    "SF6": {
        "primary": "SF6_298.1K-760.0Torr_560.0-6500.0_0.11_N2_220_43.xsc",
        "windows": [(560.0, 700.0), (700.0, 1100.0), (1100.0, 6500.0)],
        "check": (
            "SF6_296.1K_745.7Torr_800.0-1075.0_0.040_air_140.xsc",
            (800.0, 1075.0),
        ),
    }
}
OUT = Path(__file__).resolve().parents[2] / "src/optmat/data/external/hitran"


def _key() -> str:
    """Return the API key from the environment or the key file."""
    if key := os.environ.get("HITRAN_API_KEY"):
        return key.strip()
    if file := os.environ.get("HITRAN_API_KEY_FILE"):
        return Path(file).expanduser().read_text().strip()
    msg = "set HITRAN_API_KEY or HITRAN_API_KEY_FILE (a HITRAN API key)"
    raise SystemExit(msg)


def _get(section: str, params: dict[str, object]) -> dict:
    """Query an API section and return the decoded JSON."""
    key = _key()
    url = f"{HOST}/api/v2/{key}/{section}?" + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=600) as response:
            return json.loads(response.read())
    except Exception as error:  # never let the key reach a traceback or a log
        msg = f"HITRAN API request to {section!r} failed: {type(error).__name__}"
        raise SystemExit(msg) from None


def _download(gas: str) -> Path:
    """Fetch (once) the lines of every isotopologue of a gas; return the cache file."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{gas}_0_{NU_MAX:.0f}.par"
    if path.exists():
        return path
    molecule = MOLECULES[gas]
    isotopologues = _get("isotopologues", {"molecule_id__in": molecule})["content"][
        "data"
    ]
    ids = ",".join(str(iso["id"]) for iso in isotopologues)
    header = _get(
        "transitions",
        {
            "iso_ids_list": ids,
            "numin": 0,
            "numax": NU_MAX,
            "head": "False",
            "fixwidth": 0,
            "request_params": "par_line,trans_id,global_iso_id",
        },
    )
    info = _get("info", {})["content"]["data"]
    url = f"{HOST}/{info['results_dir']}/{header['content']['data']}"
    with urllib.request.urlopen(url, timeout=3600) as response:
        path.write_bytes(response.read())
    (CACHE / f"{gas}_0_{NU_MAX:.0f}.json").write_text(
        json.dumps({"message": header["message"], "timestamp": header["timestamp"]})
    )
    return path


def _cross_section(name: str) -> tuple[list[float], list[float]]:
    """Fetch (once) a HITRAN cross-section file and return wavenumbers and sigma."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / name
    if not path.exists():
        info = _get("info", {})["content"]["data"]
        url = f"{HOST}/{info['xsec_dir']}/{name}"
        with urllib.request.urlopen(url, timeout=600) as response:
            path.write_bytes(response.read())
    header, *body = path.read_text().splitlines()
    fields = header.split()
    numin, numax, count = float(fields[1]), float(fields[2]), int(fields[3])
    sigma = [float(v) for line in body for v in line.split()][:count]
    step = (numax - numin) / (count - 1)
    return [numin + i * step for i in range(count)], sigma


def _integrate(nu: list[float], sigma: list[float], lo: float, hi: float) -> tuple:
    """Return the trapezoidal integrals of sigma and sigma / nu^2 over a window."""
    pairs = [(n, s) for n, s in zip(nu, sigma, strict=True) if lo <= n < hi]
    strength = weighted = moment = 0.0
    for (n0, s0), (n1, s1) in itertools.pairwise(pairs):
        width = n1 - n0
        strength += 0.5 * (s0 + s1) * width
        weighted += 0.5 * (s0 * n0 + s1 * n1) * width
        moment += 0.5 * (s0 / n0**2 + s1 / n1**2) * width
    return strength, weighted / strength, moment


def cross_section_bands(gas: str) -> None:
    """Write the band table of a gas known from cross-sections only."""
    spec = CROSS_SECTIONS[gas]
    nu, sigma = _cross_section(spec["primary"])
    rows = []
    for lo, hi in spec["windows"]:
        strength, centre, moment = _integrate(nu, sigma, lo, hi)
        rows.append(
            [
                "xsec",
                f"{lo:.0f}-{hi:.0f} cm-1",
                "",
                "",
                centre,
                lo,
                hi,
                strength,
                moment,
                "HITRAN xsec Ref. [1164]",
            ]
        )
    _write(gas, rows)
    check_name, (lo, hi) = spec["check"]
    check = _integrate(*_cross_section(check_name), lo, hi)[0]
    primary = _integrate(nu, sigma, lo, hi)[0]
    print(
        f"{gas}: {lo:.0f}-{hi:.0f} cm-1, Ref. [1164] {primary:.4e} against "
        f"Ref. [1638] {check:.4e} cm/molecule ({primary / check - 1:+.2%})"
    )


def _write(gas: str, rows: list[list]) -> None:
    """Write a band table, strongest band first."""
    OUT.mkdir(parents=True, exist_ok=True)
    rows = sorted(rows, key=lambda row: -row[7])
    with (OUT / f"{gas}_bands.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "global_iso_id",
                "upper",
                "lower",
                "n_lines",
                "wavenumber_cm",
                "min_cm",
                "max_cm",
                "strength_cm_per_molecule",
                "strength_over_nu2_cm3_per_molecule",
                "intensity_refs",
            ]
        )
        for row in rows:
            iso_id, upper, lower, n, centre, lo, hi, strength, moment, refs = row
            writer.writerow(
                [
                    iso_id,
                    upper,
                    lower,
                    n,
                    f"{centre:.4f}",
                    f"{lo:.4f}",
                    f"{hi:.4f}",
                    f"{strength:.6e}",
                    f"{moment:.6e}",
                    refs,
                ]
            )


def main(gas: str) -> None:
    """Write the band table of one gas."""
    if gas in CROSS_SECTIONS:
        cross_section_bands(gas)
        return
    path = _download(gas)
    bands: dict[tuple[str, str, str], dict] = defaultdict(
        lambda: {
            "n": 0,
            "s": 0.0,
            "snu": 0.0,
            "s2": 0.0,
            "lo": 1e9,
            "hi": 0.0,
            "refs": set(),
        }
    )
    with path.open() as lines:
        for line in lines:
            record, *extra = line.rstrip("\r\n").split(",")
            if len(record) < 160:
                continue
            iso_id = extra[-1] if extra else record[2]
            nu = float(record[3:15])
            strength = float(record[15:25])
            key = (iso_id, record[67:82].strip(), record[82:97].strip())
            band = bands[key]
            band["n"] += 1
            band["s"] += strength
            band["snu"] += strength * nu
            band["s2"] += strength / nu**2
            band["lo"] = min(band["lo"], nu)
            band["hi"] = max(band["hi"], nu)
            band["refs"].add(record[135:137].strip())
    _write(
        gas,
        [
            [
                iso_id,
                upper,
                lower,
                band["n"],
                band["snu"] / band["s"],
                band["lo"],
                band["hi"],
                band["s"],
                band["s2"],
                " ".join(sorted(band["refs"])),
            ]
            for (iso_id, upper, lower), band in bands.items()
        ],
    )
    total = sum(band["s"] for band in bands.values())
    count = sum(band["n"] for band in bands.values())
    print(f"{gas}: {count} lines, {len(bands)} bands, total {total:.4e} cm/molecule")


if __name__ == "__main__":
    for name in sys.argv[1:]:
        main(name)
