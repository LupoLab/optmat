# NIST Atomic Spectra Database: noble-gas resonance lines

A. Kramida, Yu. Ralchenko, J. Reader and NIST ASD Team, NIST Atomic Spectra Database
(version 5.12), https://physics.nist.gov/asd, National Institute of Standards and
Technology. Queried 2026-09-24 with the lines form, tab-separated output (`format=3`),
vacuum wavelengths in nm, energies in cm^-1, oscillator strengths `fik` with their
accuracy grades (`Acc`), configurations, terms, J and statistical weights, and the
bibliographic keys of the transition probabilities (`tp_ref`) and lines (`line_ref`):

| File | Spectrum | Wavelength window (nm) |
|---|---|---|
| `He_I_lines.tsv` | He I | 40-60 |
| `Ne_I_lines.tsv` | Ne I | 40-76 |
| `Ar_I_lines.tsv` | Ar I | 70-110 |
| `Kr_I_lines.tsv` | Kr I | 80-125 |
| `Xe_I_lines.tsv` | Xe I | 90-150 |

The query string was, with the spectrum and window substituted:
`https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=Ar+I&limits_type=0&low_w=70&upp_w=110&unit=1&submit=Retrieve+Data&de=0&format=3&line_out=0&en_unit=0&output=0&bibrefs=1&page_size=15&show_obs_wl=1&show_calc_wl=1&unc_out=1&order_out=0&max_low_enrg=&show_av=2&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&f_out=on&intens_out=on&max_str=&allowed_out=1&forbid_out=1&min_accur=&min_intens=&conf_out=on&term_out=on&enrg_out=on&J_out=on&g_out=on`.

The files are kept as returned. They give the positions of the resonance lines (known to
1e-7 relative) that anchor the oscillator models' poles, and the oscillator strengths
that serve as priors, each with NIST's accuracy grade (AAA <= 0.3 %, AA <= 1 %, A+ <= 2 %,
A <= 3 %, B+ <= 7 %, B <= 10 %, C+ <= 18 %, C <= 25 %, D+ <= 40 %, D <= 50 %, E > 50 %).

## Ionisation energies

`ionization_energies.csv`: the ionisation energies of He I, He II, Ne I, Ne II, Ar I,
Ar II, Kr I, Kr II, Xe I and Xe II in eV, with their uncertainties, ground shells and
levels, the level of the ion left behind and NIST's bibliographic keys, queried
2026-09-26 with the ionization-energies form, CSV output (`format=2`, cells written as
`="..."`), kept as returned (sha256
`77ef199112bb9efa0563aac90e765c18b021edbc10652c5a932fab3e3e27b9d2`). The query was
`https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=He+I%3B+He+II%3B+Ne+I%3B+Ne+II%3B+Ar+I%3B+Ar+II%3B+Kr+I%3B+Kr+II%3B+Xe+I%3B+Xe+II&units=1&format=2&order=0&at_num_out=on&sp_name_out=on&ion_charge_out=on&el_name_out=on&seq_out=on&shells_out=on&level_out=on&ion_conf_out=on&e_out=0&unc_out=on&biblio=on`.
The source record is `nist_asd_2026`.

## Levels of the noble-gas ions

The energy levels of Ne II, Ar II, Kr II and Xe II in cm^-1 above the ion's ground
level, with their uncertainties, configurations, terms, J, Landé factors or leading
percentages and NIST's bibliographic keys, queried 2026-09-28 with the levels form,
tab-delimited output (`format=3`, cells in double quotes), kept as returned:

| File | Spectrum | Lines | sha256 |
|---|---|---|---|
| `Ne_II_levels.tsv` | Ne II | 386 | `8568223b0558d1ff31820db06ef7730985ee1f6ee51e0aa36b565ede70be8bfb` |
| `Ar_II_levels.tsv` | Ar II | 420 | `bcc7d15ea55def60e20414fa40cd764507ed43f99b7a13374b2252d515b1a0e0` |
| `Kr_II_levels.tsv` | Kr II | 164 | `b05bd725fb46a1e7b0f801de308a8bce321547ed49e64d95bdd13c2df691c11c` |
| `Xe_II_levels.tsv` | Xe II | 165 | `97801228090952bc85d0518fdacb4fd757ca2bb995b67c8c7d3d602cc2b00830` |

The query was, with the spectrum substituted (`Ne+II`, `Ar+II`, `Kr+II`, `Xe+II`):
`https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=Ar+II&units=0&format=3&output=0&page_size=15&multiplet_ordered=0&conf_out=on&term_out=on&level_out=on&unc_out=1&j_out=on&lande_out=on&perc_out=on&biblio=on&temp=&submit=Retrieve+Data`.

The columns are Configuration, Term, J, Prefix, Level (cm-1), Suffix, Uncertainty
(cm-1), then Leading percentages (Ne II, Kr II) or Lande (Ar II, Xe II), and Reference.
The second row of each file is the np5 2P*<1/2> level, the upper level of the ion's
ground term; the source record `nist_asd_2026` holds it (datasets `level_<ion>_2P12`).
NIST gives no uncertainty for that level of Kr II and Xe II, and no bibliographic key
for that of Kr II.

## Levels of the neutral noble gases

The energy levels of He I, Ne I, Ar I, Kr I and Xe I in cm^-1 above the atom's ground
level, queried 2026-10-07 with the same levels form and options as the ions' levels
above (tab-delimited output, `format=3`, cells in double quotes), kept as returned:

| File | Spectrum | Lines | sha256 |
|---|---|---|---|
| `He_I_levels.tsv` | He I | 844 | `bf3d600fbfa4f49792f4a885ab5bba438aa1092a334a9e54e885cb4a40e69b29` |
| `Ne_I_levels.tsv` | Ne I | 376 | `586e71f2f50b4b6d908f2fafa31bb1e6baa210a918e7f5cf2090e79f132adede` |
| `Ar_I_levels.tsv` | Ar I | 505 | `db515019bc50c5c99a214ea0e8de0c81067f13b1d48a073534783cd2ba2ca934` |
| `Kr_I_levels.tsv` | Kr I | 529 | `b9578d615a66beb870ad2159ad260167447258b55522c418abc5b66ec9fbaba2` |
| `Xe_I_levels.tsv` | Xe I | 446 | `0d6e09af74c0cb61f53dcd027151759a6d8937afaaf3cf167cf9a0a52bb223eb` |

The query was the ions' with the spectrum `He+I`, `Ne+I`, `Ar+I`, `Kr+I` or `Xe+I`. The
columns are Configuration, Term, J, Prefix, Level (cm-1), Suffix, Uncertainty (cm-1),
then Leading percentages (Ne I, Kr I) or Lande (where given), and Reference. NIST
brackets every helium level: "Square brackets indicate the energies determined by
interpolation, extrapolation, or other semi-empirical procedure relying on some known
experimental values" (levels help page, read 2026-10-07). NIST gives no uncertainty and
no bibliographic key for the argon and krypton levels used here, and no uncertainty for
the xenon levels. The source record `nist_asd_2026` holds, per atom, the lowest levels
that bound the frequencies at which the third-order response is off resonance (datasets
`level_<atom>_*`; the manual page *Third-order nonlinearity*).

## Lines of the neutral atoms beyond the resonance lines

The lines of Ne I, Ar I, Kr I and Xe I from the end of the resonance-line windows
above to 5000 nm: the transitions among excited levels (ns-np, np-nd, ...) with their
transition probabilities `Aki`, oscillator strengths `fik`, accuracy grades, energies,
configurations, terms and J, for the sums over intermediate states of two-photon
amplitudes (OEP 6, Q9; OEP 8). Queried 2026-10-09 with the lines form as above, but
with **vacuum wavelengths at every wavelength** (`show_av=3`), kept as returned. The
default `show_av=2` gives air wavelengths from 200 to 2000 nm. The first query of these
windows used `show_av=1`, which is vacuum below 1000 nm and *wavenumber* above: NIST
then writes minus the wavenumber in cm⁻¹, times 1e-9, into the nm columns. Those files
were replaced on the same day (TC-985); every other column was identical.

| File | Spectrum | Wavelength window (nm) | Lines (rows) | sha256 |
|---|---|---|---|---|
| `Ne_I_lines_nir.tsv` | Ne I | 76-5000 | 1483 | `de139b0e9bd6edd77caf6be0e0dd86e04ac374e6db84e9210b65f2317f798398` |
| `Ar_I_lines_nir.tsv` | Ar I | 110-5000 | 440 | `679c4ceb3b2d104738fb93265fec0b0e85bdae53f2f826e1c3f9965ea8361c65` |
| `Kr_I_lines_nir.tsv` | Kr I | 125-5000 | 746 | `3b67d4004eafad8f841b423cd2331220111b20d08f16da580daa1b25002513cb` |
| `Xe_I_lines_nir.tsv` | Xe I | 150-5000 | 980 | `e359a15ef9f62054e8a4a459abd22d7fcd4146d16aef39cab88036361c68aa3f` |

The query string was the resonance lines' with `show_av=3` and the spectrum and window
substituted:
`https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=Ar+I&limits_type=0&low_w=110&upp_w=5000&unit=1&submit=Retrieve+Data&de=0&format=3&line_out=0&en_unit=0&output=0&bibrefs=1&page_size=15&show_obs_wl=1&show_calc_wl=1&unc_out=1&order_out=0&max_low_enrg=&show_av=3&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&f_out=on&intens_out=on&max_str=&allowed_out=1&forbid_out=1&min_accur=&min_intens=&conf_out=on&term_out=on&enrg_out=on&J_out=on&g_out=on`.
Each window starts where the resonance-line file's ends, so the two files of a gas do
not overlap. Not every line has a transition probability. Of the ns-np lines (3s-3p
of Ne, 4s-4p of Ar, 5s-5p of Kr, 6s-6p of Xe), NIST gives Aki for 30 of 30 (Ne), 29 of
29 (Ar), 29 of 30 (Kr) and 22 of 26 (Xe). The source record is `nist_asd_2026`
(TC-983, corrected by TC-985).

## Terms

NIST Standard Reference Data. In the United States, works of the federal government
are not subject to copyright; under the Standard Reference Data Act, the Secretary of
Commerce secures copyright in the compilations elsewhere ("all rights reserved"), and
NIST asks for its recommended citation (above). The query output is redistributed
unchanged, with that citation, so that the values optmat derives from it can be traced
(`REUSE.toml`: `LicenseRef-optmat-third-party`).
