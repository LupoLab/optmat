# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""A relative tolerance states its absolute floor.

``pytest.approx(x, rel=r)`` keeps a default absolute tolerance of 1e-12. In SI a
polarizability is near 1e-40 and a cross section near 1e-56, so that floor would accept
any value of their size, zero included, whatever ``rel`` says. Every ``approx`` that
gives ``rel`` therefore gives ``abs`` too (``abs=0.0`` unless a floor is meant).
"""

import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]


def _is_approx(function: ast.expr) -> bool:
    return (isinstance(function, ast.Attribute) and function.attr == "approx") or (
        isinstance(function, ast.Name) and function.id == "approx"
    )


def _floorless(call: ast.Call) -> bool:
    keywords = {keyword.arg for keyword in call.keywords}
    return _is_approx(call.func) and "rel" in keywords and not {"abs", None} & keywords


@pytest.mark.parametrize(
    "path", sorted((ROOT / "tests").rglob("*.py")), ids=lambda p: p.name
)
def test_every_relative_tolerance_states_its_floor(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    missing = [
        f"{path.relative_to(ROOT)}:{node.lineno}"
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and _floorless(node)
    ]
    assert not missing, "\n".join(missing)


def test_the_check_sees_a_missing_floor():
    def call(source: str) -> ast.Call:
        node = ast.parse(source, mode="eval").body
        assert isinstance(node, ast.Call)
        return node

    assert _floorless(call("pytest.approx(x, rel=1e-9)"))
    assert _floorless(call("approx(x, rel=1e-9)"))
    assert not _floorless(call("pytest.approx(x, rel=1e-9, abs=0.0)"))
    assert not _floorless(call("pytest.approx(x)"))
    assert not _floorless(call("pytest.approx(x, abs=1e-3)"))
