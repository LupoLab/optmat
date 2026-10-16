# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract the n, k tables of the films of Franta et al. (2015, Proc. SPIE 9628).

The paper shows the films' optical constants only as figures and does not print the
UDM parameters ("will be presented elsewhere", p. 8). The authors' newAD2 site
publishes, for three of its films, the parameter file and the n, k table computed from
it, "SPIE 9628 (2015) 96281U" (attached to the paper in the library): SiO2
X2551 and MgF2 X2935 (fitted 80 cm-1 to 50 eV) and amorphous Al2O3 X2890 (80 cm-1 to
10.8 eV).

Each ``.Enk`` has 4001 rows, photon energy (eV), n, k, from 0.01 to 100 eV. Rows within
the fitted range (all rows start above 80 cm-1 = 0.0099 eV) are one dataset per
quantity; rows above the fitted upper limit are the model's extrapolation and form a
separate dataset (the record marks it excluded).
"""

from _library import Point, attachment, write

DIRECTORY = "Franta_2015_DispersionModelOptical"

FILMS = {
    # material id, file stem, upper limit of the fitted range (eV), dataset prefix
    "SiO2": ("SiO2-X2551-SPIE9628", 50.0, "X2551"),
    "MgF2": ("MgF2-X2935-SPIE9628", 50.0, "X2935"),
    "Al2O3_amorphous": ("Al2O3-X2890-SPIE9628", 10.8, "X2890"),
}


def _rows(stem: str) -> list[tuple[str, str, str, str]]:
    """Return the (energy, n, k, raw line) of a .Enk file, as printed."""
    path = attachment(DIRECTORY, f"{stem}.Enk")
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        energy, n, k = line.split()
        rows.append((energy, n, k, line.strip()))
    return rows


def extract() -> dict[str, list[Point]]:
    """Return n and k, within and beyond the fitted range, for each film."""
    files: dict[str, list[Point]] = {}
    for material, (stem, upper, prefix) in FILMS.items():
        parts: dict[str, list[Point]] = {}
        for index, (energy, n, k, raw) in enumerate(_rows(stem), start=1):
            where = "fit" if float(energy) <= upper else "extrapolated"
            for quantity, value in (("n", n), ("k", k)):
                name = f"{material}/franta_2015b__{prefix}_{quantity}_{where}.csv"
                parts.setdefault(name, []).append(
                    Point(str(index), energy, value, "", raw)
                )
        files.update(parts)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
