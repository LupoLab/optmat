# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Every register entry cited anywhere in the package exists in ``REGISTER.md``.

The review register is the evidence behind the data: code, data records and the
manual cite its entries as ``TC-n``, and a citation of an entry that does not exist is a
broken trail of provenance.
"""

import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).parents[1]
REGISTER = (ROOT / "REGISTER.md").read_text(encoding="utf-8")
ENTRIES = set(re.findall(r"^## (TC-\d+)\b", REGISTER, re.MULTILINE))
CITATION = re.compile(r"\bTC-(\d+)\b")
TEXT = {".py", ".md", ".toml", ".rst", ".txt", ".json"}


def _tracked() -> list[Path]:
    listed = subprocess.run(
        ["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    return [ROOT / name for name in listed if Path(name).suffix in TEXT]


def test_the_register_numbers_its_entries_once_and_in_order():
    numbers = [int(n) for n in re.findall(r"^## TC-(\d+)\b", REGISTER, re.MULTILINE)]
    assert numbers == sorted(numbers)
    assert len(numbers) == len(set(numbers))


def test_every_cited_entry_exists():
    missing = {}
    for path in _tracked():
        if not path.is_file():
            continue
        for number in CITATION.findall(path.read_text(encoding="utf-8")):
            if f"TC-{number}" not in ENTRIES:
                missing.setdefault(f"TC-{number}", []).append(path.name)
    assert not missing, missing


def test_every_entry_has_a_status():
    for block in re.split(r"(?m)^(?=## TC-)", REGISTER)[1:]:
        heading = block.splitlines()[0]
        assert "**Status:**" in block or "Status:" in block, heading
