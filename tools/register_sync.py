# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Reconcile provisional register numbers before merging main into a topic branch.

Preview with ``poe register-sync``; use ``--write`` to apply the displayed mapping.
Commit it, then merge main. If Git conflicts on the appended register entries,
``poe register-sync --resolve --write`` combines the two appendices. Neither mode
stages, commits or merges. Existing entry edits require manual reconciliation.
"""

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

__all__ = ["prepare", "resolve"]

HEADING = re.compile(r"^## TC-(\d+)\b", re.MULTILINE)
CITATION = re.compile(r"\b(TC-)(\d+)\b", re.IGNORECASE)


def _git(root: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=root, check=True, capture_output=True, text=True
    ).stdout


def _numbers(text: str) -> list[int]:
    numbers = [int(number) for number in HEADING.findall(text)]
    if not numbers or numbers != sorted(set(numbers)):
        msg = "register entries must have unique numbers in increasing order"
        raise ValueError(msg)
    return numbers


def _appendix(base: str, branch: str) -> str:
    _numbers(base)
    _numbers(branch)
    prefix = base.rstrip()
    if branch.rstrip() == prefix:
        return ""
    if not branch.startswith(prefix + "\n"):
        msg = "existing register text changed; reconcile those edits manually"
        raise ValueError(msg)
    return branch[len(prefix) :].lstrip("\n")


def _replace(text: str, mapping: dict[int, int]) -> str:
    # One substitution pass avoids cascading when an old and new number coincide.
    return CITATION.sub(
        lambda match: match[1] + str(mapping.get(int(match[2]), int(match[2]))), text
    )


def _apply(root: Path, changes: dict[str, str], *, write: bool) -> None:
    for name in changes:
        if (root / name).is_symlink():
            msg = f"refusing to rewrite a symlink: {name}"
            raise ValueError(msg)
    for name, content in sorted(changes.items()):
        print(f"{'Write' if write else 'Would update'} {name}")
        if write:
            (root / name).write_bytes(content.encode("utf-8"))
    if not changes:
        print("No register changes needed.")


def _unmodified_conflict(root: Path, base: str, branch: str, incoming: str) -> bool:
    with tempfile.TemporaryDirectory() as directory:
        files = [Path(directory) / name for name in ("branch", "base", "incoming")]
        for path, content in zip(files, (branch, base, incoming), strict=True):
            path.write_text(content, encoding="utf-8")
        result = subprocess.run(
            ["git", "merge-file", "--stdout", *map(str, files)],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
        )
    if result.returncode != 1:
        return False
    # Git labels conflicts with branch names; merge-file labels them with paths.
    labels = re.compile(r"(?m)^(<<<<<<<|\|\|\|\|\|\|\||>>>>>>>) .*$")
    actual = (root / "REGISTER.md").read_text(encoding="utf-8")
    return labels.sub(r"\1", actual) == labels.sub(r"\1", result.stdout)


def prepare(root: Path, onto: str = "main", *, write: bool = False) -> None:
    """Renumber only branch additions and their tracked UTF-8 citations.

    The merge base supplies the fixed numbering boundary, so no offset file or
    shared reservation is needed. The working tree must be clean; commit all new
    source files first. This function leaves incoming entries to the Git merge.
    """
    if _git(root, "status", "--porcelain"):
        msg = "commit branch work first; register preparation needs a clean worktree"
        raise ValueError(msg)
    target = _git(root, "rev-parse", "--verify", f"{onto}^{{commit}}").strip()
    base = _git(root, "merge-base", "HEAD", target).strip()
    common = _git(root, "show", f"{base}:REGISTER.md")
    incoming = _git(root, "show", f"{target}:REGISTER.md")
    branch = (root / "REGISTER.md").read_text(encoding="utf-8")
    fixed = _numbers(common)
    if _numbers(branch) == fixed:
        print("No branch-only register entries.")
        return
    appendix = _appendix(common, branch)
    incoming_numbers = _numbers(incoming)
    if incoming_numbers[: len(fixed)] != fixed:
        msg = "incoming register changes established entry numbers; review manually"
        raise ValueError(msg)
    mapping = {
        old: incoming_numbers[-1] + offset
        for offset, old in enumerate(_numbers(appendix), start=1)
        if old != incoming_numbers[-1] + offset
    }
    print(f"Register base {base[:12]}; target {onto} at {target[:12]}.")
    for old, new in mapping.items():
        print(f"TC-{old} -> TC-{new}")
    changes = {}
    changed_paths = set(
        _git(root, "diff", "--name-only", "-z", base, "HEAD").split("\0")
    )
    for name in _git(root, "ls-files", "-z").split("\0"):
        if not name:
            continue
        path = root / name
        if path.is_symlink():
            msg = f"review tracked symlink manually before renumbering: {name}"
            raise ValueError(msg)
        data = path.read_bytes()
        if b"\0" in data:
            continue
        try:
            content = data.decode("utf-8")
        except UnicodeDecodeError:
            continue
        updated = _replace(content, mapping)
        if updated != content:
            if name not in changed_paths:
                msg = f"a provisional citation occurs in an unchanged file: {name}"
                raise ValueError(msg)
            changes[name] = updated
    _apply(root, changes, write=write)
    if changes:
        print(
            "Review all citations, including private notes, and commit before merging."
        )


def resolve(root: Path, *, write: bool = False) -> None:
    """Combine a register append conflict after numbers have been reconciled.

    Read Git's three index stages, preserving the incoming register verbatim and
    appending this branch's entries. Refuse overlapping numbers, old-entry edits on
    this branch, or a target that advanced past the prepared numbering.
    """
    base = _git(root, "show", ":1:REGISTER.md")
    branch = _git(root, "show", ":2:REGISTER.md")
    incoming = _git(root, "show", ":3:REGISTER.md")
    appendix = _appendix(base, branch)
    fixed = _numbers(base)
    numbers = _numbers(incoming)
    if numbers[: len(fixed)] != fixed:
        msg = "incoming register changes established entry numbers; review manually"
        raise ValueError(msg)
    if not appendix or _numbers(appendix)[0] <= numbers[-1]:
        msg = (
            "branch numbers overlap or precede main; reconcile manually before merging"
        )
        raise ValueError(msg)
    if (root / "REGISTER.md").read_text(
        encoding="utf-8"
    ) == incoming.rstrip() + "\n\n" + appendix:
        print("Register conflict already resolved; review and stage it.")
        return
    if not _unmodified_conflict(root, base, branch, incoming):
        msg = "the register conflict was edited or is unusual; reconcile it manually"
        raise ValueError(msg)
    _apply(root, {"REGISTER.md": incoming.rstrip() + "\n\n" + appendix}, write=write)
    print(
        "Review REGISTER.md, resolve other conflicts, then stage and commit the merge."
    )


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--onto", default="main", help="local integration ref (default: main)"
    )
    parser.add_argument(
        "--write", action="store_true", help="apply the previewed changes"
    )
    parser.add_argument(
        "--resolve", action="store_true", help="resolve the register append conflict"
    )
    args = parser.parse_args()
    try:
        root = Path(_git(Path.cwd(), "rev-parse", "--show-toplevel").strip())
        if args.resolve:
            resolve(root, write=args.write)
        else:
            prepare(root, args.onto, write=args.write)
    except (ValueError, subprocess.CalledProcessError) as error:
        detail = (
            error.stderr.strip()
            if isinstance(error, subprocess.CalledProcessError)
            else str(error)
        )
        print(f"Register sync refused: {detail}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(_main())
