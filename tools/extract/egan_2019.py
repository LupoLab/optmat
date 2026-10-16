# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Egan, Stone, Scherschligt & Harvey (2019): molar refractivities at 633 nm.

Table III (``tbl_6d3145de2fcf``, CSV; image checked) compares the present A_R
(cm^3/mol, 632.9908 nm, 293.1529 K) of Ne, Ar, Xe, N2, CO2 and N2O with the literature.
The extraction merged each gas's rows into one cell per column: the value cell holds
all the values run together (``0.999(1)1.001(1)...``) and the reference cell all the
references (``Burns et al. (Ref. 45) Hohm and Kerl (Ref. 46) ... This work``). They are
split on the uncertainty parentheses and on "(Ref. n)"; the image shows one value per
reference, in the same order, and the script checks that the counts agree. The note
mark ``a`` of Egan and Stone's argon value (adjusted by +19.3e-6 A_R) is stripped. One
file per gas and reference; the parenthesised last-digit standard uncertainty goes to
the sigma column.
"""

import re

from _library import Point, table, write

DIRECTORY = "Egan_2019_MeasuredRelationshipBetween"
TABLE = "tbl_6d3145de2fcf"
GASES = {"Ne": "Ne", "Ar": "Ar", "Xe": "Xe", "N2": "N2", "CO2": "CO2", "N2O": "N2O"}
# Reference number -> file stem.
REFERENCES = {
    "12": "egan_2016",
    "25": "achtermann_1991",
    "26": "achtermann_1993",
    "27": "hohm_1994",
    "45": "burns_1986",
    "46": "hohm_kerl_1990",
    "47": "birch_1991",
    "48": "coulon_1981",
    "49": "hou_1994",
    "50": "egan_2011_adjusted",
    "51": "hohm_trumper_1996",
    "52": "montixi_1983",
    "53": "achtermann_1986",
}
_VALUE = re.compile(r"(?P<number>\d+\.\d+)\((?P<digits>\d+)\)a?")
_REFERENCE = re.compile(r"\(Ref\. (?P<ref>\d+)\)|This work")
WAVELENGTH = {"present": "632.9908", "literature": "633"}


def _split(number: str, digits: str) -> tuple[str, str]:
    """Return a value and its last-digit uncertainty in the value's unit."""
    decimals = len(number.split(".")[1])
    return number, f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and reference of Table III."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE)[2:]:
        gas = GASES[cells[0].strip()]
        values = list(_VALUE.finditer(cells[1]))
        references = list(_REFERENCE.finditer(cells[2]))
        if "".join(m[0] for m in values) != cells[1].strip() or len(values) != len(
            references
        ):
            msg = f"{gas}: {len(values)} values for {len(references)} references"
            raise ValueError(msg)
        printed = " | ".join(cells)
        for value, reference in zip(values, references, strict=True):
            number, sigma = _split(value["number"], value["digits"])
            if reference["ref"] is None:
                name, wavelength = "A_R", WAVELENGTH["present"]
            else:
                name, wavelength = (
                    REFERENCES[reference["ref"]],
                    WAVELENGTH["literature"],
                )
            path = f"{gas}/egan_2019__{name}_{gas}.csv"
            row = f"{value[0]} | {reference[0]}"
            files[path] = [
                Point(wavelength, wavelength, number, sigma, f"{row} || {printed}")
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
