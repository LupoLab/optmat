# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Zahedpour, Hancock & Milchberg (2019): Table 1, electronic n2 of N2, O2, Ar.

Table 1 is on p. 846 (PDF page 4). It has two blocks of rows, each headed by a row of
pump wavelengths in um: 3.0-6.5 (measured at 1 atm) with a last column n2,rot, and
7.0-10.5 (measured at 42 atm) with a last column n2,eff(tau_peak). Every measured cell
prints ``value ± error`` in 1e-20 cm2/W, "at atmospheric pressure". The CSV matches the
table image cell for cell (checked 2026-09-26); the header cell "10.5^c" carries the
footnote mark c (the comparison with Pigeon et al. at 10.6 um), stripped here.

Written files:

- one file per gas with the 16 measured n2, rows by pump wavelength (um);
- n2,rot of N2 and O2 (first block, last column; Ar's printed 0 is not written), and
  n2,eff(tau_peak) of N2 and O2 (second block, last column; the pulse width printed in
  the cell goes to the dataset's parameters; Ar's "8.5 (average)" is not written, see
  the record). These are single values with no pump wavelength: their abscissa is
  ``not stated``.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Zahedpour_2019_UltrashortInfrared25"
TABLE = "tbl_6c9829dbaca8"
GASES = {r"\text{N}_2": "N2", r"\text{O}_2": "O2", "Ar": "Ar"}
HEADER_ROWS = (0, 4)  # the two rows of pump wavelengths
N_WAVELENGTHS = 8  # measured columns per block; the ninth is n2,rot or n2,eff
_EFFECTIVE = re.compile(r"^\s*([\d.]+)\s*\(\s*(\d+)\s*fs\s*\)\s*$")


def _value_error(cell: str) -> tuple[str, str]:
    """Split ``value ± error`` into the printed numbers."""
    value, error = cell.split("±")
    return number(value), number(error)


def _wavelength(cell: str) -> str:
    """Return a header wavelength as printed, without a footnote mark (``10.5^c``)."""
    return number(cell.split("^", maxsplit=1)[0])


def extract() -> dict[str, list]:
    """Return the n2 file per gas and the n2,rot and n2,eff files of N2 and O2."""
    rows = table(DIRECTORY, TABLE)
    files: dict[str, list] = {}
    for header in HEADER_ROWS:
        header_cells = rows[header][1 : 1 + N_WAVELENGTHS]
        wavelengths = [_wavelength(cell) for cell in header_cells]
        for cells in rows[header + 1 : header + 4]:
            gas = GASES[cells[0]]
            printed = " | ".join(cells)
            path = f"{gas}/zahedpour_2019__n2_{gas}.csv"
            for wavelength, cell in zip(
                wavelengths, cells[1 : 1 + N_WAVELENGTHS], strict=True
            ):
                value, error = _value_error(cell)
                files.setdefault(path, []).append(
                    Point(wavelength, wavelength, value, error, printed)
                )
            if gas == "Ar":
                continue  # n2,rot printed 0; n2,eff an average of the measured n2
            last = cells[1 + N_WAVELENGTHS]
            if header == 0:
                name, value = "n2_rot", number(last)
            else:
                match = _EFFECTIVE.match(last)
                if match is None:
                    msg = f"unexpected n2,eff cell {last!r}"
                    raise ValueError(msg)
                name, value = "n2_eff", number(match.group(1))
            files[f"{gas}/zahedpour_2019__{name}_{gas}.csv"] = [
                Point(name, "not stated", value, "", printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
