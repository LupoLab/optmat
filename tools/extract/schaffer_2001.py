# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Schaffer, Brodeur & Mazur (2001): bulk damage thresholds of 4 dielectrics.

Table 1 (``tbl_dbd271b58058``, JSON cells only: merged header) gives per material the
band gap and, at 400 and 800 nm (110 fs), the fitted threshold intensity I_th and the
critical power P_cr used in the fit of Eq. (7). Checked against the table image: the
cells agree. Every I_th is printed ``a \times 10^{13}``; the script checks the exponent
and writes the mantissa (the dataset's scale is 1e13). SF11 has no 400-nm values (a
dash: no damage at 400 nm, p. 1792).

Table 2 (``tbl_ec8b126f418c``) is scrambled in the library: its one data row holds the
characters of every cell of each column interleaved (e.g. the material column
``SCFCFuoasF1ren12dinsgil0ic2a11``). STRUCTURAL FIX: the table is transcribed below
from the table image (``tables/tbl_ec8b126f418c.jpg``), and the script checks that each
scrambled cell is exactly the multiset of the characters of the transcribed column (the
em dash is extracted as "-", subscripts and spaces dropped); a mismatch stops the
extraction. The electron densities N (Keldysh photoionization only, at the Table 1
threshold intensities) become ``ionization_yield`` datasets whose abscissa is that
intensity; the Keldysh parameters gamma have no quantity and are quoted in the record.
"""

from collections import Counter

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Schaffer_2001_LaserInducedBreakdown"
TABLE1 = "tbl_dbd271b58058"
TABLE2 = "tbl_ec8b126f418c"

MATERIALS = {
    "SF11": "SF11",
    "Corning 0211": "Corning-0211",
    "Fused silica": "SiO2",
    "CaF_2": "CaF2",
    "CaF2": "CaF2",
}
STEMS = {"SF11": "SF11", "Corning-0211": "Corning0211", "SiO2": "SiO2", "CaF2": "CaF2"}
EXPONENT = r"\times 10^{13}"
DASH = "\N{EM DASH}"
TIMES = "\N{MULTIPLICATION SIGN}"

# Table 2 as printed (the image): material, gamma and N at 400 nm, gamma and N at
# 800 nm. N is (mantissa, power of ten) of the printed "a x 10^b"; a dash is no value.
Density = tuple[str, str] | None
TABLE2_IMAGE: list[tuple[str, str, Density, str, Density]] = [
    ("SF11", DASH, None, "1.2", ("3.0", "21")),
    ("Corning 0211", "2.7", ("5.2", "21"), "1.0", ("2.7", "21")),
    ("Fused silica", "2.5", ("3.7", "20"), "1.2", ("4.1", "19")),
    ("CaF2", "2.3", ("4.8", "19"), "1.2", ("2.5", "18")),
]


def _characters(text: str) -> Counter[str]:
    """Return the characters of a cell, without spaces, the em dash read as "-"."""
    return Counter(text.replace(DASH, "-").replace(" ", ""))


def _printed(density: Density) -> str:
    """Return a density as the characters the print shows, "a x 10b" or a dash."""
    if density is None:
        return DASH
    mantissa, power = density
    return f"{mantissa}{TIMES}10{power}"


def _check_table2() -> None:
    """Check the transcription of Table 2 against the scrambled library cells."""
    scrambled = table(DIRECTORY, TABLE2)[2]
    columns = [
        "".join(row[0] for row in TABLE2_IMAGE),
        "".join(row[1] for row in TABLE2_IMAGE),
        "".join(_printed(row[2]) for row in TABLE2_IMAGE),
        "".join(row[3] for row in TABLE2_IMAGE),
        "".join(_printed(row[4]) for row in TABLE2_IMAGE),
    ]
    for index, (cell, transcribed) in enumerate(zip(scrambled, columns, strict=True)):
        if _characters(cell) != _characters(transcribed):
            msg = f"Table 2 column {index}: {cell!r} does not scramble {transcribed!r}"
            raise ValueError(msg)


def _mantissa(cell: str) -> str:
    r"""Return the mantissa of an I_th cell printed ``a \times 10^{13}``."""
    if not cell.endswith(EXPONENT):
        msg = f"I_th cell {cell!r} is not printed x 10^13"
        raise ValueError(msg)
    return number(cell)


def extract() -> dict[str, list[Point]]:
    """Return Table 1's band gaps, thresholds and critical powers, and Table 2's N."""
    _check_table2()
    files: dict[str, list[Point]] = {}
    intensities: dict[tuple[str, str], str] = {}
    for cells in json_table(DIRECTORY, TABLE1)[2:]:
        material = MATERIALS[cells[0]]
        stem = STEMS[material]
        printed = " | ".join(cells)
        files[f"{material}/schaffer_2001__band_gap_{stem}.csv"] = [
            Point("1", "", number(cells[1]), "", printed)
        ]
        threshold = []
        critical = []
        for wavelength, (i_col, p_col) in (("400", (2, 3)), ("800", (4, 5))):
            if cells[i_col] == DASH:
                continue
            mantissa = _mantissa(cells[i_col])
            intensities[material, wavelength] = mantissa
            threshold.append(Point(wavelength, wavelength, mantissa, "", printed))
            power = number(cells[p_col])
            critical.append(Point(wavelength, wavelength, power, "", printed))
        files[f"{material}/schaffer_2001__threshold_{stem}.csv"] = threshold
        if material == "Corning-0211":
            # P_cr at 800 nm was a free parameter of the fit; at 400 nm it was scaled
            # from it with Eq. (4) (p. 1792): two datasets of different standing.
            fit = f"{material}/schaffer_2001__critical_power_fit_{stem}.csv"
            scaled = f"{material}/schaffer_2001__critical_power_scaled_{stem}.csv"
            files[fit] = [p for p in critical if p.row == "800"]
            files[scaled] = [p for p in critical if p.row == "400"]
        else:
            files[f"{material}/schaffer_2001__critical_power_{stem}.csv"] = critical
    for name, gamma400, n400, gamma800, n800 in TABLE2_IMAGE:
        material = MATERIALS[name]
        stem = STEMS[material]
        for wavelength, gamma, density in (
            ("400", gamma400, n400),
            ("800", gamma800, n800),
        ):
            if density is None:
                continue
            mantissa, power = density
            printed = f"gamma {gamma} | N {mantissa} x 10^{power} (from the image)"
            x = f"{intensities[material, wavelength]}e13"
            name_ = f"{material}/schaffer_2001__keldysh_density_{stem}_{wavelength}.csv"
            files[name_] = [Point("1", x, f"{mantissa}e{power}", "", printed)]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
