# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The line identifiers of the curated Raman responses (OEP 1, section 6.7).

    uv run python tools/raman_line_ids.py

writes ``tests/data/raman_line_ids.json``: for each curated gas and mixture at 295 K,
the ``line_id`` of every line, in order. ``tests/test_raman_groups.py`` compares the
responses with it, so an identifier cannot change unnoticed; a line that is added or
removed is re-recorded deliberately, with a changelog line.
"""

import json
from pathlib import Path

from optmat.gas.raman_response import raman_response

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "tests" / "data" / "raman_line_ids.json"
GASES = ("N2", "O2", "H2", "D2", "CH4", "SF6", "CO2", "N2O", "dry air")
TEMPERATURE = 295.0


def snapshot() -> dict[str, list[str]]:
    """Return every curated response's line identifiers."""
    return {gas: list(raman_response(gas, TEMPERATURE).line_id) for gas in GASES}


def main() -> None:
    """Write the snapshot."""
    OUTPUT.write_text(json.dumps(snapshot(), indent=1) + "\n", encoding="utf-8")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
