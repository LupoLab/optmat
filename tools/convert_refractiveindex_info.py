# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Convert the refractiveindex.info YAML snapshot to one JSON file the package reads.

The snapshot (``data/external/refractiveindex_info/``, see its PROVENANCE.md) is YAML;
the package reads JSON with the standard library, so this script, run by hand when the
snapshot is refreshed, needs PyYAML only for itself::

    uv run --with pyyaml python tools/convert_refractiveindex_info.py

It writes two files, one record per entry with its references and comments as text, its
conditions, and its DATA blocks unchanged (formula type, wavelength range in um,
coefficients; or tabulated rows):

- ``data/external/refractiveindex_info.json``: the gases' ``nk`` entries, with the gas;
- ``data/external/refractiveindex_info_solid.json``: the solids' ``nk`` and ``n2``
  entries, with the optmat material id, the kind (``nk`` or ``n2``) and, for a uniaxial
  crystal's ``-o``/``-e`` page, the axis (TC-320).
"""

import json
import re
from pathlib import Path
from typing import Any

import yaml  # pyright: ignore[reportMissingModuleSource] - run with `uv run --with pyyaml`

DATA = Path(__file__).resolve().parents[1] / "src" / "optmat" / "data"
SNAPSHOT = DATA / "external" / "refractiveindex_info"


GASES = {
    "He",
    "Ne",
    "Ar",
    "Kr",
    "Xe",
    "H2",
    "D2",
    "N2",
    "O2",
    "SF6",
    "N2O",
    "CO2",
    "CH4",
}

# The database's book (directory) of each solid, and the optmat material it is. The
# database files Yb-doped KGW under KGW's n2 only; the dopant is recorded per entry.
SOLIDS = {
    "SiO2": "SiO2",
    "MgF2": "MgF2",
    "CaF2": "CaF2",
    "BaF2": "BaF2",
    "LiF": "LiF",
    "Al2O3": "Al2O3",
    "Y3Al5O12": "YAG",
    "Si": "Si",
    "BK7": "BK7",
    "optical": "N-BK7",  # specs/schott/optical/N-BK7.yml, Schott's catalogue formula
    "Yb-KGd(WO4)2": "KGW",
}


# Pages whose axis is named only in the catalogue (catalog-nk.yml, catalog-n2.yml at the
# pinned commit): "Malitson 1962: alpha-Al2O3 (Sapphire); n(o)" and
# "DeSalvo et al. 1996: (e)".
AXES = {"main/Al2O3/nk/Malitson.yml": "o", "main/Al2O3/n2/DeSalvo.yml": "e"}


def _blocks(content: dict[str, Any]) -> list[dict[str, object]]:
    blocks = []
    for block in content.get("DATA", []):
        kind = block["type"]
        if kind.startswith("formula"):
            blocks.append(
                {
                    "type": kind,
                    "wavelength_range_um": [
                        float(x) for x in str(block["wavelength_range"]).split()
                    ],
                    "coefficients": [
                        float(x) for x in str(block["coefficients"]).split()
                    ],
                }
            )
        else:
            rows = [
                [float(x) for x in line.split()]
                for line in str(block["data"]).strip().splitlines()
            ]
            blocks.append({"type": kind, "rows": rows})
    return blocks


def _solid_entries() -> list[dict[str, object]]:
    records = []
    for path in sorted(SNAPSHOT.rglob("*.yml")):
        kind = path.parent.name if path.parent.name in {"nk", "n2"} else "nk"
        book = path.parent.parent.name if kind == path.parent.name else path.parent.name
        if book in GASES or book == "CH4 - methane" or book not in SOLIDS:
            continue
        content = yaml.safe_load(path.read_text(encoding="utf-8"))
        match = re.search(r"-([oe])(-|$)", path.stem)
        axis = AXES.get(
            str(path.relative_to(SNAPSHOT)), match.group(1) if match else ""
        )
        records.append(
            {
                "material": SOLIDS[book],
                "kind": kind,
                "axis": axis,
                "dopant": "Yb" if book.startswith("Yb-") else "",
                "entry": path.stem,
                "file": str(path.relative_to(SNAPSHOT)),
                "references": content.get("REFERENCES", ""),
                "comments": content.get("COMMENTS", ""),
                "conditions": content.get("CONDITIONS", {}),
                "data": _blocks(content),
            }
        )
    return records


def _entries() -> list[dict[str, object]]:
    records = []
    for path in sorted(SNAPSHOT.glob("*/*/nk/*.yml")):
        gas = path.parent.parent.name
        if gas not in GASES | {"CH4 - methane"}:
            continue
        content = yaml.safe_load(path.read_text(encoding="utf-8"))
        blocks = []
        for block in content.get("DATA", []):
            kind = block["type"]
            if kind.startswith("formula"):
                blocks.append(
                    {
                        "type": kind,
                        "wavelength_range_um": [
                            float(x) for x in str(block["wavelength_range"]).split()
                        ],
                        "coefficients": [
                            float(x) for x in str(block["coefficients"]).split()
                        ],
                    }
                )
            else:
                rows = [
                    [float(x) for x in line.split()]
                    for line in str(block["data"]).strip().splitlines()
                ]
                blocks.append({"type": kind, "rows": rows})
        records.append(
            {
                "gas": gas,
                "entry": path.stem,
                "file": str(path.relative_to(SNAPSHOT)),
                "references": content.get("REFERENCES", ""),
                "comments": content.get("COMMENTS", ""),
                "conditions": content.get("CONDITIONS", {}),
                "data": blocks,
            }
        )
    return records


def main() -> None:
    """Write the JSON files."""
    for name, records in (
        ("refractiveindex_info.json", _entries()),
        ("refractiveindex_info_solid.json", _solid_entries()),
    ):
        out = DATA / "external" / name
        out.write_text(
            json.dumps(records, indent=1, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        print(f"wrote {out.relative_to(DATA)} ({len(records)} entries)")


if __name__ == "__main__":
    main()
