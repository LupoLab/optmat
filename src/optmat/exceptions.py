# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The package's exception and warning hierarchy.

Every error the package raises derives from :class:`OptmatError`, so a caller can
catch the package's refusals without catching unrelated ones.
"""

__all__ = [
    "ArgumentError",
    "DataError",
    "ModelAbsorptionWarning",
    "ModelRangeWarning",
    "OptmatError",
    "StateError",
    "UnavailableModelError",
]


class OptmatError(Exception):
    """Base class of every error raised by optmat."""


class StateError(OptmatError, ValueError):
    """A thermodynamic state (pressure, temperature) is unphysical or out of domain."""


class DataError(OptmatError, ValueError):
    """A curated data file is malformed, inconsistent or refers to something missing."""


class ArgumentError(OptmatError, ValueError):
    """An argument is outside its domain: not finite, or of the wrong sign."""


class UnavailableModelError(OptmatError, LookupError):
    """No curated model exists for what was asked: the data do not support one."""


class ModelRangeWarning(UserWarning):
    """A model is evaluated outside the range over which its source supports it."""


class ModelAbsorptionWarning(UserWarning):
    """A real-part method is evaluated where the gas absorbs noticeably.

    The real polarizability or ``n - 1`` alone then leaves out the absorption: use
    :meth:`~optmat.gas.models.base.Model.complex_polarizability` or
    :meth:`~optmat.gas.models.base.Model.refractive_index`.
    """
