# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Every text file is read and written as UTF-8, whatever the platform's default.

The curated records are UTF-8 (author names, degree signs, micro signs), and Python
3.12 to 3.14 on Windows decodes a file without a stated encoding with the system code
page (cp1252), which fails on them. A text-mode ``open``, ``Path.open``, ``gzip.open``,
``read_text`` or ``write_text`` therefore names its encoding.
"""

import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
CHECKED = ("src", "tests", "tools", "docs")
SKIPPED = {"_build", "auto_gallery", "generated", ".venv", "build"}
TEXT_METHODS = {"read_text", "write_text"}
OTHER_OPENS = {"tarfile", "zipfile", "h5py", "os"}


def _files():
    for directory in CHECKED:
        for path in sorted((ROOT / directory).rglob("*.py")):
            if not SKIPPED.intersection(path.relative_to(ROOT).parts):
                yield path


def _mode(call: ast.Call) -> object:
    for keyword in call.keywords:
        if keyword.arg == "mode":
            value = keyword.value
            return value.value if isinstance(value, ast.Constant) else None
    # The mode follows the path in open() and gzip.open(); it is first in Path.open().
    method = isinstance(call.func, ast.Attribute) and not _is_gzip(call)
    position = 0 if method else 1
    if len(call.args) > position:
        argument = call.args[position]
        return argument.value if isinstance(argument, ast.Constant) else None
    return "rb" if _is_gzip(call) else "r"


def _is_gzip(call: ast.Call) -> bool:
    function = call.func
    return (
        isinstance(function, ast.Attribute)
        and isinstance(function.value, ast.Name)
        and function.value.id == "gzip"
    )


def _missing_encoding(call: ast.Call) -> bool:
    if any(keyword.arg in {"encoding", None} for keyword in call.keywords):
        return False
    function = call.func
    if isinstance(function, ast.Attribute) and function.attr in TEXT_METHODS:
        return True
    is_open = (isinstance(function, ast.Name) and function.id == "open") or (
        isinstance(function, ast.Attribute)
        and function.attr == "open"
        and not (
            isinstance(function.value, ast.Name) and function.value.id in OTHER_OPENS
        )
    )
    if not is_open:
        return False
    mode = _mode(call)
    # A mode built at run time is not judged here.
    return isinstance(mode, str) and "b" not in mode


@pytest.mark.parametrize("path", list(_files()), ids=lambda p: str(p.relative_to(ROOT)))
def test_text_files_name_their_encoding(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    missing = [
        f"{path.relative_to(ROOT)}:{node.lineno}: {ast.unparse(node)[:80]}"
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and _missing_encoding(node)
    ]
    assert not missing, "\n".join(missing)


def test_the_check_sees_each_kind_of_call():
    calls = {
        "path.read_text()": True,
        'path.read_text(encoding="utf-8")': False,
        'path.write_text("x")': True,
        "open(path)": True,
        'open(path, "rb")': False,
        'open(path, "w", encoding="utf-8")': False,
        'path.open("r")': True,
        'path.open("rb")': False,
        'gzip.open(path, "rt")': True,
        "gzip.open(path)": False,
        'gzip.open(path, "rt", encoding="utf-8")': False,
        "tarfile.open(path)": False,
    }
    for source, expected in calls.items():
        call = ast.parse(source, mode="eval").body
        assert isinstance(call, ast.Call)
        assert _missing_encoding(call) is expected, source
