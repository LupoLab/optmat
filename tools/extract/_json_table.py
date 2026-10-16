# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Read a library table that has no CSV rendering, from its ``tables/<id>.json`` cells.

Some library tables (those with merged header cells, mostly) are published only as
``.json``/``.html``. This helper returns them as the same list of rows of cell strings
as ``_library.table``: each cell's text at its anchor (first row and column), and ""
under a row or column span, so that column indices match the printed layout.
"""

import json

from _library import LIBRARY


def json_table(directory: str, table_id: str) -> list[list[str]]:
    """Return a library table's rows (headers first) from its JSON cells."""
    path = LIBRARY / directory / "tables" / f"{table_id}.json"
    record = json.loads(path.read_text(encoding="utf-8"))
    grid = [["" for _ in range(record["columns"])] for _ in range(record["rows"])]
    for cell in record["cells"]:
        grid[cell["row"]][cell["column"]] = cell["text"].strip()
    return grid
