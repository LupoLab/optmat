# NIST Chemistry WebBook: fluid properties at 0.1 MPa (Ne, Kr, Xe, N2O)

NIST Chemistry WebBook, SRD 69 (NIST Standard Reference Database 69), National
Institute of Standards and Technology, fluid data ("Isobaric Properties"),
https://webbook.nist.gov (retrieved 2026-09-30). The captured pages carry no citation of
the fluid-data section; each names its equation of state and its transport model.

The isobaric tables at 0.1 MPa were fetched from the WebBook's fluid pages, in K, MPa,
mol/l, kJ/mol, m/s, uPa s and N/m, with the default reference state, as tab-separated
text (`Action=Data&Wide=on`), and are kept byte for byte:

| File | WebBook ID | Temperatures (K) | Rows | SHA-256 | Captured (UTC) |
|---|---|---|---:|---|---|
| `Ne_0.1MPa_200-400K.tsv` | C7440019 | 200-400, step 5 | 41 | `9cae296eeb1d0c00bde99c940e17f1015bca2b18205b744f17c52baeee109f06` | 2026-09-30 15:43:37 |
| `Ne_0.1MPa_293.15-298.15K.tsv` | C7440019 | 293.15, 298.15 | 2 | `7dd0fee3be703284552b68c998665a380adfe0470297720652070eac736e2500` | 2026-09-30 15:43:40 |
| `Kr_0.1MPa_200-400K.tsv` | C7439909 | 200-400, step 5 | 41 | `cfe3ae90a39c04489e5a130405e551691a1f52fdd6cdf51c2365d32ec64dd0ef` | 2026-09-30 15:43:43 |
| `Kr_0.1MPa_293.15-298.15K.tsv` | C7439909 | 293.15, 298.15 | 2 | `608975ce28ee69635ca46a99dc8408925f8e4a6ef03b4cdda6e010fa2392c347` | 2026-09-30 15:43:45 |
| `Xe_0.1MPa_200-400K.tsv` | C7440633 | 200-400, step 5 | 41 | `b8f3b17945b877840134c804881cf1fe88049081311472ec604a2a300be9e56a` | 2026-09-30 15:43:48 |
| `Xe_0.1MPa_293.15-298.15K.tsv` | C7440633 | 293.15, 298.15 | 2 | `7cb6ce1be3cf3138a85b1113aa280c8400e3707218fd9a7094d8694f82bcf303` | 2026-09-30 15:43:51 |
| `N2O_0.1MPa_200-400K.tsv` | C10024972 | 200-400, step 5 | 41 | `b922b622af65c902bdbbcd345823def7f86ef3d212d8640175d4af842553ab0c` | 2026-09-30 15:43:54 |
| `N2O_0.1MPa_293.15-298.15K.tsv` | C10024972 | 293.15, 298.15 | 2 | `351722345c6517f8a1f28f1c95c7df860d7d2a2d4cbdd92dac2287030ac53fd8` | 2026-09-30 15:43:57 |

The query string was, with the ID and temperatures substituted:
`https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C7440019&Type=IsoBar&P=0.1&THigh=400&TLow=200&TInc=5&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=mol%2Fl&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm`.
Every row is in the vapour phase.

**Role.** A check, not a source of values. The pages name the transport model as
Huber's, documented in NISTIR 8209 (`huber_2018`), whose dilute-gas form optmat
evaluates for these four gases. The tables are at 0.1 MPa, above the dilute limit by the
density term: at 298.15 K by 0.013 % (Ne), 0.13 % (Kr), 0.21 % (Xe) and 0.21 % (N2O),
at 200 K by 0.028 %, 0.33 %, 0.54 % and 0.52 %, falling with temperature, all below the
models' stated uncertainties
(`tests/test_transport.py`).

**Uncertainty as the pages state it** (gas at atmospheric pressure): Ne 3 %, Kr 3 %, Xe
1 %, N2O 10 %. Huber's report states 2 %, 3 %, 0.7 % and 5 % (k = 2) for the same
models; the report governs.

# NIST Chemistry WebBook: constants of diatomic molecules (H2, D2, N2, O2, CO)

The WebBook's "Constants of diatomic molecules" pages, compiled by K. P. Huber and G. H.
Herzberg (data collected through February 1977; Huber & Herzberg, *Molecular Spectra and
Molecular Structure IV. Constants of Diatomic Molecules*, Van Nostrand Reinhold, 1979),
fetched on 2026-10-07 as HTML from
`https://webbook.nist.gov/cgi/cbook.cgi?ID=<ID>&Mask=1000` and stored byte for byte,
gzip-compressed (`mtime = 0`):

| File | WebBook ID | Bytes (uncompressed) | SHA-256 of the uncompressed page |
|---|---|---:|---|
| `H2_diatomic_constants.html.gz` | C1333740 | 261137 | `a4e9a291fca5b46ae80648d26eab0700b0eb6969571981a5ed404a71d0c864f9` |
| `D2_diatomic_constants.html.gz` | C7782390 | 112228 | `88c0508bcfd37be39d57d9189ba99c6b6b4c6b14389cbe8f69eaf98c0a476971` |
| `N2_diatomic_constants.html.gz` | C7727379 | 277230 | `cc3e4e37c13872d3b72c03526e9c78ed51d92094f0c37d6cc79de136cf97ed7d` |
| `O2_diatomic_constants.html.gz` | C7782447 | 264662 | `2479ad7e5635a58c3ac59e2710027f18194b3dc7cfba7db5dabe6b973bf951f9` |
| `CO_diatomic_constants.html.gz` | C630080 | 225870 | `11aaeb4806adacc9f8c9b46b2dae85aac0b28e67870226b2c4b255d2056e2a2f` |

Each page's table gives, per electronic state, Te (cm^-1), the vibrational and
rotational constants, the observed transitions (Trans.) and the position of the 0-0
band (ν00). Parentheses and brackets around a value mark it as doubtful or deduced (the
page's symbol notes). The source record `nist_webbook_2026` holds the Te of the states
that bound the off-resonant third-order response (datasets `state_<gas>_<state>`), as
reproduced from Huber & Herzberg.

## Terms

NIST Standard Reference Data. In the United States, works of the federal government
are not subject to copyright; under the Standard Reference Data Act, the Secretary of
Commerce secures copyright in the compilations elsewhere ("all rights reserved"), and
NIST asks for its recommended citation (above). The query output is redistributed
unchanged, with that citation, so that the values optmat derives from it can be traced
(`REUSE.toml`: `LicenseRef-optmat-third-party`).
