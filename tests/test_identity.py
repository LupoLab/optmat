# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The identity core: the canonical record, its digest, the status and ownership.

OEP 1, section 4.8: identical inputs give identical digests in every process, any
change to a parameter changes the digest, and a model owns read-only storage.
"""

import dataclasses
import enum
import pickle
import subprocess
import sys
import types
from collections import namedtuple
from typing import ClassVar

import numpy as np
import pytest

from optmat.exceptions import DataError
from optmat.identity import (
    STATUSES,
    FrozenMapping,
    Identified,
    ModelIdentity,
    canonical,
    canonical_json,
    frozen,
    parse_status,
    sha256_digest,
)

# A shared evaluator module for the toy model, registered so importlib finds it.
EVALUATOR = types.ModuleType("optmat_test_evaluator")
setattr(EVALUATOR, "EVALUATION_VERSION", 3)  # noqa: B010 - a module has no such slot
sys.modules[EVALUATOR.__name__] = EVALUATOR


@dataclasses.dataclass(frozen=True)
class Toy(Identified):
    IDENTITY_KIND: ClassVar[str] = "toy"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = (EVALUATOR.__name__,)

    name: str
    values: np.ndarray
    table: dict[str, float] = dataclasses.field(default_factory=dict)

    def _identity_label(self):
        return f"toy/{self.name}", self.name, ("toy_2026",)


@dataclasses.dataclass(frozen=True)
class Bound(Identified):
    """A model bound to another, whose digest enters its record."""

    IDENTITY_KIND: ClassVar[str] = "toy"
    EVALUATION_VERSION: ClassVar[int] = 1

    part: Toy

    def _identity_label(self):
        return "bound", "x", ()

    def _identity_extras(self):
        return {"resolved": 1.5}


class Colour(enum.StrEnum):
    RED = "red"


@dataclasses.dataclass(frozen=True)
class Plain:
    x: float


@dataclasses.dataclass(frozen=True)
class Versioned:
    EVALUATION_VERSION: ClassVar[int] = 7
    x: float


def _toy(**changes):
    base = {"name": "a", "values": np.array([1.0, 2.0]), "table": {"k": 1.0}}
    return Toy(**(base | changes))


class TestStatus:
    def test_the_vocabulary(self):
        assert STATUSES == (
            "validated",
            "provisional",
            "untested",
            "comparison",
            "user",
        )

    def test_a_default_applies_where_none_is_declared(self):
        assert parse_status({}, "x", default="comparison") == ("comparison", "")

    @pytest.mark.parametrize(
        ("entry", "match"),
        [
            ({"status": "good"}, "not one of"),
            ({"status": "user"}, "not one of"),
            ({"status": "provisional"}, "needs a status_reason"),
            ({"status": "untested"}, "needs a status_reason"),
        ],
    )
    def test_bad_statuses_are_refused(self, entry, match):
        with pytest.raises(DataError, match=match):
            parse_status(entry, "x", default="validated")


class TestFrozenMapping:
    def test_it_reads_like_a_dict_and_copies(self):
        source = {"a": 1}
        m = FrozenMapping(source)
        source["a"] = 2
        assert m["a"] == 1
        assert list(m) == ["a"]
        assert len(m) == 1
        assert m == {"a": 1}
        assert repr(m) == "FrozenMapping({'a': 1})"

    def test_it_cannot_be_changed_or_hashed(self):
        m = FrozenMapping({"a": 1})
        with pytest.raises(TypeError):
            m["a"] = 2  # pyright: ignore[reportIndexIssue]
        with pytest.raises(TypeError):
            hash(m)

    def test_it_pickles(self):
        m = FrozenMapping({"a": (1, 2)})
        assert pickle.loads(pickle.dumps(m)) == m


class TestFrozen:
    def test_an_array_is_copied_read_only(self):
        a = np.arange(3.0)
        b = frozen(a)
        assert b is not a
        assert not b.flags.writeable
        a[0] = 9.0
        assert b[0] == 0.0

    def test_a_read_only_array_is_still_copied(self):
        # A caller could make an array it owns writeable again, so it is never shared.
        a = np.arange(3.0)
        a.flags.writeable = False
        assert frozen(a) is not a

    def test_collections(self):
        m = frozen({"x": [1, np.zeros(1)], "y": (2,)})
        assert isinstance(m, FrozenMapping)
        assert isinstance(m["x"], tuple)
        assert not m["x"][1].flags.writeable
        pair = namedtuple("pair", "a b")(1, 2)
        assert frozen(pair) is pair
        assert frozen("s") == "s"


class TestCanonical:
    @pytest.mark.parametrize(
        ("value", "form"),
        [
            (None, None),
            (True, True),
            (3, 3),
            ("s", "s"),
            (Colour.RED, "red"),
            (0.1, 0.1),
            (float("inf"), "inf"),
            (float("-inf"), "-inf"),
            (float("nan"), "nan"),
            (1 + 2j, [1.0, 2.0]),
            (np.float64(0.5), 0.5),
            (np.bool_(True), True),
            (np.array([1.0, 2.0]), [1.0, 2.0]),
            (np.array([1, 2]), [1.0, 2.0]),
            (np.array([1.0, np.inf]), [1.0, "inf"]),
            (np.array([1 + 1j]), [[1.0, 1.0]]),
            (np.array([True]), [True]),
            ({Colour.RED: (1, [2])}, {"red": [1, [2]]}),
            ({3, 1}, [1, 3]),
            (Plain(1.0), {"x": 1.0, "__type__": f"{__name__}.Plain"}),
            (
                Versioned(1.0),
                {"x": 1.0, "__type__": f"{__name__}.Versioned", "__version__": 7},
            ),
        ],
    )
    def test_forms(self, value, form):
        assert canonical(value) == form

    @pytest.mark.parametrize("value", [object(), Plain])
    def test_no_form(self, value):
        with pytest.raises(TypeError, match="no canonical form"):
            canonical(value)

    def test_json_is_sorted_and_compact(self):
        assert canonical_json({"b": [1.5], "a": "é"}) == '{"a":"é","b":[1.5]}'.encode()

    def test_the_digest_is_sha256(self):
        # The SHA-256 of b"[]", FIPS 180-4.
        assert sha256_digest([]) == (
            "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945"
        )


class TestIdentified:
    def test_the_identity(self):
        toy = _toy(status="validated", status_reason="checked")
        identity = toy.identity
        assert isinstance(identity, ModelIdentity)
        assert identity.name == "toy/a"
        assert identity.kind == "toy"
        assert identity.material == "a"
        assert identity.sources == ("toy_2026",)
        assert (identity.status, identity.status_reason) == ("validated", "checked")
        assert len(identity.digest) == 64
        assert identity.modified == ()
        assert identity.base == ""
        assert identity.reference_temperature is None
        assert identity.infrared_temperature is None
        assert toy.identity is identity

    def test_the_record(self):
        record = _toy().identity_record()
        assert record["class"] == f"{__name__}.Toy"
        assert record["versions"] == {f"{__name__}.Toy": 1, EVALUATOR.__name__: 3}
        assert set(record["fields"]) == {"name", "values", "table"}
        assert record["extras"] == {}

    def test_a_curated_status_does_not_enter_the_digest(self):
        assert _toy().identity.digest == _toy(status="validated").identity.digest
        assert _toy().status == "user"

    @pytest.mark.parametrize(
        "change",
        [
            {"name": "b"},
            {"values": np.array([1.0, 2.0000000000000004])},
            {"table": {"k": 2.0}},
            {"table": {"j": 1.0}},
        ],
    )
    def test_any_change_changes_the_digest(self, change):
        assert _toy(**change).identity.digest != _toy().identity.digest

    def test_an_evaluator_version_enters_the_digest(self, monkeypatch):
        before = _toy().identity.digest
        monkeypatch.setattr(EVALUATOR, "EVALUATION_VERSION", 4)
        assert _toy().identity.digest != before

    def test_storage_is_owned_and_read_only(self):
        values = np.array([1.0, 2.0])
        toy = _toy(values=values)
        values[0] = 5.0
        assert toy.values[0] == 1.0
        with pytest.raises(ValueError, match="read-only"):
            toy.values[0] = 3.0
        with pytest.raises(TypeError):
            toy.table["k"] = 3.0  # pyright: ignore[reportIndexIssue]
        with pytest.raises(dataclasses.FrozenInstanceError):
            toy.name = "c"  # pyright: ignore[reportAttributeAccessIssue]

    def test_replace_rebuilds_owned_storage(self):
        toy = dataclasses.replace(_toy(), values=np.array([3.0]))
        assert not toy.values.flags.writeable

    def test_pickling_keeps_storage_read_only_and_the_digest(self):
        toy = _toy()
        digest = toy.identity.digest
        copy = pickle.loads(pickle.dumps(toy))
        assert not copy.values.flags.writeable
        assert isinstance(copy.table, FrozenMapping)
        assert copy.identity.digest == digest

    def test_a_bound_components_digest_enters_the_record(self):
        bound = Bound(part=_toy())
        record = bound.identity_record()
        assert record["fields"]["part"] == {
            "__type__": f"{__name__}.Toy",
            "digest": _toy().identity.digest,
        }
        assert record["extras"] == {"resolved": 1.5}
        assert Bound(part=_toy(name="z")).identity.digest != bound.identity.digest

    def test_the_digest_is_the_same_in_another_process(self):
        code = (
            "import numpy as np\n"
            "from optmat.identity import sha256_digest\n"
            "print(sha256_digest({'a': np.array([0.1, 1e300]), 'b': {3, 1, 2}, "
            "'c': (1 + 2j, float('nan'))}))"
        )
        here = sha256_digest(
            {
                "a": np.array([0.1, 1e300]),
                "b": {3, 1, 2},
                "c": (1 + 2j, float("nan")),
            }
        )
        for seed in ("0", "1"):
            out = subprocess.run(
                [sys.executable, "-c", code],
                capture_output=True,
                text=True,
                check=True,
                env={"PYTHONHASHSEED": seed},
            )
            assert out.stdout.strip() == here
