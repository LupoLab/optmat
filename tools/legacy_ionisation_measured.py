# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Luna.jl's and solprop's ionisation rates through optmat's tests on measured yields.

    uv run python tools/legacy_ionisation_measured.py

reads the rates of ``tests/data/legacy_ionisation_rates.json``
(``tools/legacy_ionisation_rates.jl`` and ``.py``), interpolates their logarithm
linearly in the logarithm of the field (extrapolated below the grid with its first
slope, held above it), cycle-averages them numerically, and puts them through
:func:`optmat.gas.ionisation_validation.intensity_scale_test` (Wallace et al.) and
:func:`~optmat.gas.ionisation_validation.measured_yields` (Wahlstrand et al.) in place
of an optmat model; optmat's default is printed for comparison (TC-815). It also puts
them through :func:`~optmat.gas.ionisation_validation.phase_shift_test` (the
measurements of the phase shift at 800 nm, linear polarisation: Nambu et al., Nie et
al.) and :func:`~optmat.gas.ionisation_validation.helium_laser_rates` (Parker et al.'s
two-electron helium at 780 nm, against the codes' rates at 800 nm, up to 1e15 W/cm^2,
the end of the fixture's fields) (TC-829). Luna stops above twice the
barrier-suppression field (xenon above 3.5e14 W/cm^2, inside Wallace et al.'s range):
there its last rate is held, which Luna itself does not do.
"""

import contextlib
import dataclasses
import json
import warnings
from collections.abc import Iterator
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike

from optmat.datasets import datasets
from optmat.exceptions import ModelRangeWarning
from optmat.gas import ionisation_validation as validation
from optmat.gas.ionisation_rates import IonisationRateModel, ionisation_rate_model

FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "tests"
    / "data"
    / "legacy_ionisation_rates.json"
)
GASES = ("Ar", "Kr", "Xe")
CODES = ("luna", "solprop")
PARKER_LIMIT = 1e19
"""W/m^2: the largest intensity of the fixture's fields (1e15 W/cm^2)."""


class TabulatedRate(IonisationRateModel):
    """A code's rate from the fixture, as an optmat model of the same gas.

    A subclass of :class:`~optmat.gas.ionisation_rates.IonisationRateModel`, so that
    the validators accept it where their types are checked: its two rates are the
    code's, interpolated; every other attribute (the barrier-suppression field that
    bounds some comparisons) is that of optmat's default of the gas at the fixture's
    wavelength.
    """

    @classmethod
    def from_fixture(cls, legacy: dict, code: str, gas: str) -> "TabulatedRate":
        """Return a code's rate of a gas from the fixture.

        Parameters
        ----------
        legacy : dict
            The fixture.
        code : str
            ``luna`` or ``solprop``.
        gas : str
            The gas.

        Returns
        -------
        TabulatedRate
            The rate.
        """
        base = ionisation_rate_model(gas, float(legacy["wavelength"]))
        model = cls(
            **{
                f.name: getattr(base, f.name)
                for f in dataclasses.fields(base)
                if f.init
            }
        )
        field = np.array(legacy["field_v_per_m"], dtype=float)
        rate = np.array(
            [np.nan if r is None else r for r in legacy[code]["gases"][gas]["rate"]],
            dtype=float,
        )
        kept = np.isfinite(rate) & (rate > 0)
        order = np.argsort(field[kept])
        object.__setattr__(model, "_log_field", np.log(field[kept][order]))
        object.__setattr__(model, "_log_rate", np.log(rate[kept][order]))
        return model

    def instantaneous_rate(self, field_strength: ArrayLike) -> np.ndarray:
        """Return the code's rate, 1/s, at a field, V/m.

        Parameters
        ----------
        field_strength : array_like
            The field, V/m.

        Returns
        -------
        numpy.ndarray
            The rate: log-linear in the logarithm of the field between the fixture's
            fields, extrapolated below them with the first slope, held above them.
        """
        log_field, log_rate = self.__dict__["_log_field"], self.__dict__["_log_rate"]
        field = np.abs(np.asarray(field_strength, dtype=float))
        rate = np.zeros_like(field)
        nonzero = field > 0
        x = np.log(field[nonzero])
        slope = (log_rate[1] - log_rate[0]) / (log_field[1] - log_field[0])
        below = log_rate[0] + slope * (x - log_field[0])
        inside = np.interp(x, log_field, log_rate)
        rate[nonzero] = np.exp(np.where(x < log_field[0], below, inside))
        return rate

    def cycle_averaged_rate(self, amplitude: ArrayLike) -> np.ndarray:
        """Return the code's rate averaged over a cycle of a field of this amplitude.

        Parameters
        ----------
        amplitude : array_like
            The amplitude, V/m.

        Returns
        -------
        numpy.ndarray
            The average of :meth:`instantaneous_rate` over the phase (400 points).
        """
        phase = np.linspace(0.0, np.pi / 2, 400)
        values = np.asarray(amplitude, dtype=float)
        rate = self.instantaneous_rate(values[..., None] * np.cos(phase))
        return np.trapezoid(rate, phase, axis=-1) / (np.pi / 2)


@contextlib.contextmanager
def substituted(legacy: dict, code: str, gas: str) -> Iterator[None]:
    """Make the validators' rate model a code's tabulated rate, for one gas.

    Parameters
    ----------
    legacy : dict
        The fixture.
    code : str
        ``luna``, ``solprop``, ``optmat`` (optmat's default, unchanged) or
        ``optmat_800nm`` (the default at the fixture's 800 nm, whatever the test's
        wavelength).
    gas : str
        The gas.

    Yields
    ------
    None
        While the substitution holds.
    """
    original = validation.ionisation_rate_model
    if code == "optmat_800nm":
        fixed = original(gas, float(legacy["wavelength"]))
        validation.ionisation_rate_model = (  # pyright: ignore[reportAttributeAccessIssue]
            lambda *_args, rate=fixed, **_kwargs: rate
        )
    elif code != "optmat":
        tabulated = TabulatedRate.from_fixture(legacy, code, gas)
        validation.ionisation_rate_model = (  # pyright: ignore[reportAttributeAccessIssue]
            lambda *_args, rate=tabulated, **_kwargs: rate
        )
    try:
        yield
    finally:
        validation.ionisation_rate_model = original


def phase_shift_references() -> list[str]:
    """Return the measured phase shifts of the noble gases at 800 nm, linear light.

    Returns
    -------
    list of str
        Dataset keys.
    """
    return sorted(
        d.key
        for d in datasets()
        if d.quantity == "ionization_phase_shift"
        and d.material in ("He", "Ne", "Ar", "Kr", "Xe")
        and abs(float(str(d.parameters["wavelength_nm"])) - 800) < 15
        and d.parameters["polarization"] == "linear"
    )


def parker_ratios() -> np.ndarray:
    """Return the rate over Parker et al.'s at 780 nm, to the fixture's last field.

    Returns
    -------
    numpy.ndarray
        The ratios, with the substitution in force.
    """
    infrared, _ = validation.helium_laser_rates()
    return infrared.ratio[infrared.x <= PARKER_LIMIT]


def _line(label: str, gas: str) -> str:
    scale = validation.intensity_scale_test(gas)
    measured = validation.measured_yields(gas)
    return (
        f"{label:8s} {gas}: Wallace scale {scale.scale:.3f} (calibration "
        f"{100 * scale.calibration:.1f} %, spread {scale.residual:.3f}); over "
        f"Wahlstrand {measured.ratio.min():.2f} to {measured.ratio.max():.2f}"
    )


def main() -> None:
    """Print the tests of each code's rate."""
    legacy = json.loads(FIXTURE.read_text(encoding="utf-8"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        for code in (*CODES, "optmat"):
            for gas in GASES:
                with substituted(legacy, code, gas):
                    print(_line(code, gas))
            for key in phase_shift_references():
                gas = key.split("_")[-2] if key.endswith("_LP") else key.split("_")[-1]
                with substituted(legacy, code, gas):
                    test = validation.phase_shift_test(key)
                print(f"{code:8s} {key}: intensity scale {test.scale:.3f}")
            with substituted(legacy, code, "He"):
                ratio = parker_ratios()
            low, high = ratio.min(), ratio.max()
            print(f"{code:8s} He: over Parker et al. {low:.2f} to {high:.2f}")
        with substituted(legacy, "optmat_800nm", "He"):
            ratio = parker_ratios()
        print(f"optmat at 800 nm, He: {ratio.min():.2f} to {ratio.max():.2f}")


if __name__ == "__main__":
    main()
