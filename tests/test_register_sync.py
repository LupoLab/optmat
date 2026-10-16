# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Concurrent branches retain their evidence and citations through integration."""

import importlib.util
import re
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parents[1] / "tools" / "register_sync.py"
SPEC = importlib.util.spec_from_file_location("register_sync", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SYNC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SYNC)


def git(root, *args, check=True):
    return subprocess.run(
        [
            "git",
            "-c",
            "user.name=Register test",
            "-c",
            "user.email=register@example.invalid",
            "-c",
            "commit.gpgsign=false",
            "-c",
            "core.hooksPath=/dev/null",
            *args,
        ],
        cwd=root,
        text=True,
        capture_output=True,
        check=check,
    )


def commit(root):
    git(root, "add", ".")
    git(root, "commit", "-m", "Synthetic register fixture")


def entry(number, subject):
    return f"\n## TC-{number} {subject}\n\n- **Status:** decided (2026-09-30).\n"


@pytest.fixture
def branches(tmp_path):
    git(tmp_path, "init", "-b", "main")
    register = tmp_path / "REGISTER.md"
    common = "# Synthetic register\n" + entry(1, "Established evidence")
    register.write_text(common)
    commit(tmp_path)
    git(tmp_path, "switch", "-c", "topic")
    register.write_text(common + entry(2, "Branch first") + entry(3, "Branch second"))
    (tmp_path / "branch notes.md").write_text(
        "TC-1 stays; TC-2 and TC-3 move once; TC-20 stays. Link: #tc-2-branch-first\n"
    )
    (tmp_path / "source.toml").write_text('evidence = "TC-2, TC-3"\n')
    commit(tmp_path)
    git(tmp_path, "switch", "main")
    register.write_text(common + entry(2, "Incoming evidence"))
    (tmp_path / "incoming.md").write_text("Incoming evidence is TC-2.\n")
    commit(tmp_path)
    git(tmp_path, "switch", "topic")
    return tmp_path


@pytest.mark.parametrize("style", ["merge", "diff3", "zdiff3"])
def test_preview_renumber_and_merge_preserve_both_branches(branches, style, capsys):
    root = branches
    git(root, "config", "merge.conflictStyle", style)
    before = (root / "REGISTER.md").read_bytes()
    SYNC.prepare(root)
    assert "TC-2 -> TC-3" in capsys.readouterr().out
    assert (root / "REGISTER.md").read_bytes() == before
    assert git(root, "status", "--porcelain").stdout == ""
    SYNC.prepare(root, write=True)
    assert (root / "branch notes.md").read_text() == (
        "TC-1 stays; TC-3 and TC-4 move once; TC-20 stays. Link: #tc-3-branch-first\n"
    )
    assert (root / "source.toml").read_text() == 'evidence = "TC-3, TC-4"\n'
    assert not (root / "incoming.md").exists()
    commit(root)
    # Repeating preparation against the same main must not move entries again.
    SYNC.prepare(root, write=True)
    assert git(root, "status", "--porcelain").stdout == ""
    merged = git(root, "merge", "--no-edit", "main", check=False)
    assert merged.returncode == 1
    assert "CONFLICT (content): Merge conflict in REGISTER.md" in merged.stdout
    conflicted = (root / "REGISTER.md").read_bytes()
    SYNC.resolve(root)
    assert (root / "REGISTER.md").read_bytes() == conflicted
    SYNC.resolve(root, write=True)
    resolved = (root / "REGISTER.md").read_text()
    incoming = git(root, "show", "main:REGISTER.md").stdout
    assert resolved.startswith(incoming.rstrip() + "\n\n")
    assert re.findall(r"^## TC-(\d+)", resolved, re.MULTILINE) == ["1", "2", "3", "4"]
    assert "TC-3 Branch first" in resolved
    assert "TC-4 Branch second" in resolved
    assert (root / "incoming.md").read_text() == "Incoming evidence is TC-2.\n"
    SYNC.resolve(root, write=True)
    assert (root / "REGISTER.md").read_text() == resolved
    # The tool deliberately leaves review, staging and the merge commit to its caller.
    assert git(root, "ls-files", "-u").stdout
    commit(root)
    SYNC.prepare(root, write=True)
    assert git(root, "status", "--porcelain").stdout == ""


def test_prepare_refuses_uncommitted_files_and_edits_to_established_entries(branches):
    root = branches
    notes = root / "branch notes.md"
    notes.write_text(notes.read_text() + "Work in progress.\n")
    register = root / "REGISTER.md"
    before = register.read_bytes()
    with pytest.raises(ValueError, match="clean worktree"):
        SYNC.prepare(root, write=True)
    assert register.read_bytes() == before
    commit(root)
    register.write_text(
        register.read_text().replace("Established evidence", "Edited evidence")
    )
    commit(root)
    before = register.read_bytes()
    with pytest.raises(ValueError, match="existing register text changed"):
        SYNC.prepare(root, write=True)
    assert register.read_bytes() == before


@pytest.mark.parametrize("prepared", [False, True])
def test_resolve_refuses_collisions_including_main_advancing_after_preparation(
    branches, prepared
):
    root = branches
    if prepared:
        SYNC.prepare(root, write=True)
        commit(root)
        git(root, "switch", "main")
        register = root / "REGISTER.md"
        register.write_text(register.read_text() + entry(3, "Main advanced"))
        commit(root)
        git(root, "switch", "topic")
    assert git(root, "merge", "--no-edit", "main", check=False).returncode == 1
    before = (root / "REGISTER.md").read_bytes()
    with pytest.raises(ValueError, match="overlap or precede main"):
        SYNC.resolve(root, write=True)
    assert (root / "REGISTER.md").read_bytes() == before


def test_resolve_preserves_manual_conflict_edits(branches):
    root = branches
    SYNC.prepare(root, write=True)
    commit(root)
    assert git(root, "merge", "--no-edit", "main", check=False).returncode == 1
    register = root / "REGISTER.md"
    register.write_text(register.read_text() + "Manual evidence added during review.\n")
    before = register.read_bytes()
    with pytest.raises(ValueError, match="conflict was edited"):
        SYNC.resolve(root, write=True)
    assert register.read_bytes() == before


def test_no_new_entries_need_no_renumbering_even_when_the_preamble_changes(branches):
    root = branches
    git(root, "switch", "main")
    git(root, "switch", "-c", "wording")
    register = root / "REGISTER.md"
    register.write_text(
        register.read_text().replace("# Synthetic register", "# Wording")
    )
    commit(root)
    before = register.read_bytes()
    SYNC.prepare(root, write=True)
    assert register.read_bytes() == before
    assert git(root, "status", "--porcelain").stdout == ""
