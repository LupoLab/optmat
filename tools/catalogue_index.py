# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Write the index of the curated catalogue.

    uv run python tools/catalogue_index.py

writes ``src/optmat/data/catalogue_index.json`` from the source records: which record
holds the datasets and formulas of each material and quantity (:mod:`optmat._index`).
Run it after adding or editing a source record; ``tests/test_catalogue_index.py`` fails
until the committed index matches the records.
"""

from optmat._index import INDEX_FILE, build, dumps

if __name__ == "__main__":
    INDEX_FILE.write_text(dumps(build()), encoding="utf-8")
    print(f"wrote {INDEX_FILE}")
