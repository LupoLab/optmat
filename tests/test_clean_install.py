# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The release test: the wheel and the sdist install and run cleanly (OEP 1, 6.12).

Each artifact is built from this checkout and installed with the ``gas`` extra into a
fresh environment, offline (from uv's cache), away from the source tree and without the
paper library. In that environment the quickstart's snippets and the propagation page's
run, and every name of the public surface's snapshot resolves. That proves the data
files ship. The lane is ``poe test-release``; it takes minutes.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
UV = shutil.which("uv")

pytestmark = pytest.mark.release

CHECK = """
import importlib, json, re, sys
from pathlib import Path

import optmat

source = Path(sys.argv[1])
assert source not in Path(optmat.__file__).parents, optmat.__file__
blocks = re.compile(r"^```python\\n(.*?)^```", re.MULTILINE | re.DOTALL)
for page in (
    "docs/introduction/installation.md", "docs/quickstart.md", "docs/propagation.md"
):
    namespace = {}
    for number, block in enumerate(blocks.findall((source / page).read_text()), 1):
        exec(compile(block, f"{page}, block {number}", "exec"), namespace)
api = json.loads((source / "tests" / "data" / "public_api.json").read_text())
for module, names in api.items():
    loaded = importlib.import_module(module)
    for name in names:
        assert hasattr(loaded, name), (module, name)
print(len(api), sum(len(names) for names in api.values()))
"""


def _run(*args: str, cwd: Path) -> str:
    """Run a command offline and return its output, failing on an error."""
    env = {k: v for k, v in os.environ.items() if k != "PAPER_LIBRARY"}
    env |= {"UV_OFFLINE": "1", "PYTHONNOUSERSITE": "1"}
    result = subprocess.run(
        args, cwd=cwd, env=env, capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


@pytest.fixture(scope="module")
def artifacts(tmp_path_factory):
    """Build the wheel and the sdist from this checkout."""
    if UV is None:  # pragma: no cover - the lane runs where uv is installed
        pytest.skip("the release test needs uv")
    out = tmp_path_factory.mktemp("dist")
    _run(UV, "build", "--offline", "--out-dir", str(out), cwd=ROOT)
    return {"wheel": next(out.glob("*.whl")), "sdist": next(out.glob("*.tar.gz"))}


def test_the_sdist_holds_the_tracked_tree_only(artifacts):
    with tarfile.open(artifacts["sdist"]) as sdist:
        names = sdist.getnames()
    top = names[0].split("/")[0]
    assert not [n for n in names if n.startswith(f"{top}/.claude/")]
    assert not [n for n in names if "/docs/_build/" in n]
    tracked = subprocess.run(
        ["git", "ls-files", "-z", "src/optmat/data"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split("\0")[:-1]
    assert {f"{top}/{n}" for n in tracked} <= set(names)


def test_the_wheel_ships_every_data_file(artifacts):
    with zipfile.ZipFile(artifacts["wheel"]) as wheel:
        names = set(wheel.namelist())
    tracked = subprocess.run(
        ["git", "ls-files", "-z", "src/optmat/data"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split("\0")[:-1]
    assert {n.removeprefix("src/") for n in tracked} <= names


@pytest.mark.parametrize("kind", ["wheel", "sdist"])
def test_a_clean_install_runs_the_quickstart_and_the_surface(artifacts, kind, tmp_path):
    assert UV is not None
    venv = tmp_path / "venv"
    _run(
        UV,
        "venv",
        "--offline",
        "-q",
        str(venv),
        "--python",
        sys.executable,
        cwd=tmp_path,
    )
    python = venv / "bin" / "python"
    _run(
        UV, "pip", "install", "--offline", "-q", "--python", str(python),
        f"{artifacts[kind]}[gas]", cwd=tmp_path,
    )  # fmt: skip
    script = tmp_path / "check.py"
    script.write_text(CHECK, encoding="utf-8")
    out = _run(str(python), str(script), str(ROOT), cwd=tmp_path)
    modules, names = (int(x) for x in out.split())
    assert modules == len(
        json.loads((ROOT / "tests/data/public_api.json").read_text(encoding="utf-8"))
    )
    assert names > 0
    assert re.fullmatch(r"\d+ \d+\n", out)
