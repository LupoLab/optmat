# refractiveindex.info database snapshot

Gas entries of the refractiveindex.info database (M. N. Polyanskiy, "Refractiveindex.info
database of optical constants", Sci. Data 11, 94 (2024), doi:10.1038/s41597-023-02898-2),
repository https://github.com/polyanskiy/refractiveindex.info-database, commit
c5c2f188e848453def5970e347399d653df2ffc2 (branch main, fetched 2026-09-24).

Licence: CC0 1.0 (`LICENSE`, `legalcode.txt`, `permission.txt` copied from the repository).

Copied: `database/data/main/<gas>/` for He, Ne, Ar, Kr, Xe, H2, D2, N2, O2, SF6, N2O
and CO2, and `database/data/organic/CH4 - methane/` (renamed `organic/CH4/`), without
the liquid and solid entries. Used **only as a comparison source**: nothing here enters
a fit, and every value here that comes from a paper is curated from the paper itself.

## Solid entries (added 2026-09-26, same commit; TC-320)

Fetched from `raw.githubusercontent.com` at the same commit, c5c2f18, and kept at their
paths under `database/data/`. The scope:

- **nk:**
  - bulk fused silica (`main/SiO2/nk/`: Malitson, Nyakuchena, Arosa, Franta, Popova;
    not the crystalline-quartz or thin-film pages);
  - MgF2, CaF2, BaF2, LiF, Al2O3 (sapphire) and Y3Al5O12 (YAG): the crystal pages, not
    the thin films;
  - Si: the crystalline pages, not amorphous silicon;
  - BK7: `glass/optical/BK7/nk/Lane.yml`, and Schott's N-BK7 catalogue formula
    `specs/schott/optical/N-BK7.yml`.
  - KGW has no nk page in the database.
- **n2:** every n2 page of the same materials:
  - fused silica, but not quartz's `Adair-o` and `Adair-e`;
  - `glass/optical/BK7/n2/`;
  - `other/doped crystals/Yb-KGd(WO4)2/n2/` (Yb-doped KGW, the database's only KGW n2
    page).

The catalogues that selected them are `database/catalog-nk.yml` and
`database/catalog-n2.yml` at the same commit. `tools/convert_refractiveindex_info.py`
writes them to `data/external/refractiveindex_info_solid.json`, with the optmat
material id, the kind (nk or n2), the axis of a uniaxial crystal's page and the dopant.

Known errors in the files, recorded in the register rather than edited here:
- `main/SiO2/nk/Popova.yml` cites Kitamura et al., Appl. Opt. **33**, 8118 (2007); the
  volume is **46**.

## Stahrenberg et al.'s silver and copper (added 2026-09-28, same commit)

`main/Ag/nk/Stahrenberg.yml` and `main/Cu/nk/Stahrenberg.yml`, fetched from
`raw.githubusercontent.com` at commit c5c2f18. They are the only publication of the
numbers behind Stahrenberg et al., Phys. Rev. B 64, 115111 (2001), whose paper gives
figures only: the files say "Numerical data kindly provided by Prof. Dr. Norbert Esser".
They are taken as the source of those numbers (decided 2026-09-28): they are
curated in `stahrenberg_2001` as data reproduced through the database, not as a
comparison. SHA-256: Ag `0919e06925c847acf0b681c9d685455d60b88d267ab726aac2dd2c23be9b6fe4`
(identical to solprop's pinned copy), Cu
`854458100a2c3a509cb14aef188a940337e883d4b40205e980c1f8196d00e0d9`.

## Franta et al.'s silicon at 20 C, 25 C and 300 K (curated 2026-10-14, same commit)

`main/Si/nk/Franta-20C.yml`, `Franta-25C.yml` and `Franta-300K.yml`, from the solid
entries above. The pages cite Franta et al., Appl. Surf. Sci. 421, 405 (2017) with
"Numerical data kindly provided by Daniel Franta", and say that the numerical
procedures were updated in J. Appl. Phys. 123, 185707 (2018). The paper prints its
model's parameters only in part, and no n, k table: the pages are the only numerical
form of the temperature-dependent model. Like Stahrenberg's, they are taken as the
source of those numbers (TC-1075): curated in `franta_2017` as data, the 20 C page as
the fit's (TC-1088), the others as checks. The other Franta pages (10-500 K) stay
comparison files. SHA-256: 20C
`2f692bcd8b141e66785192acf2b4ae004016a9a8de4051d08a692bb0a343376a`, 25C
`c3e042b71ffeff1d801c9cb16123bb07132cbb7abb762e2553b23fc4bf4cb16a`, 300K
`45e271e0e8d20aa9617da5e98522cbbf29867ab7294d382760188b4d1fa5f969`.

## Terms

The database is dedicated to the public domain (CC0 1.0; `LICENSE`, `legalcode.txt`,
`permission.txt`). The waiver is the database's: `specs/schott/optical/N-BK7.yml` holds
the Sellmeier coefficients of Schott's catalogue, which refractiveindex.info reproduces
from Schott's published data sheet, and Schott is cited for it.
