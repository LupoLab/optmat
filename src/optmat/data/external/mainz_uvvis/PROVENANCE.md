# MPI-Mainz UV/VIS Spectral Atlas: ozone absorption cross sections

Retrieved on 2026-09-25 from the open-access MPI-Mainz UV/VIS Spectral Atlas
of Gaseous Molecules of Atmospheric Interest, which asks to be cited as:

H. Keller-Rudek, G. K. Moortgat, R. Sander and R. Sörensen, "The MPI-Mainz UV/VIS
spectral atlas of gaseous molecules of atmospheric interest", Earth Syst. Sci. Data 5,
365-373 (2013), doi:10.5194/essd-5-365-2013.

The files were fetched from
`https://uv-vis-spectral-atlas-mainz.org/uvvis_data/cross_sections/Ozone/<atlas name>`.
Each atlas page (`.../uvvis/cross_sections/Ozone/<atlas name>`) gives the original
bibliography and a comment. Each file is stored byte for byte, gzip-compressed
(`mtime = 0`), under a name without parentheses or spaces. Every file has two columns,
wavelength (nm) and cross section (cm^2 per molecule), unless noted.

The original measurements must be cited with the atlas. Their curated source records
hold what each paper states: uncertainty, wavelength scale, temperature and resolution.

| File | Atlas name | SHA-256 of the uncompressed file | Lines | Role |
|---|---|---|---:|---|
| `O3_Serdyuchenko2014_293K_213-1100nm2013_version.txt.gz` | `O3_Serdyuchenko(2014)_293K_213-1100nm(2013 version).txt` | `0d4da0b946e64b20b006655385600d797e773075a406d0fe90f254498adb23b7` | 88668 | primary, 213-1100 nm |
| `O3_Serdyuchenko2014_243K_213-1100nm2013_version.txt.gz` | `O3_Serdyuchenko(2014)_243K_213-1100nm(2013 version).txt` | `afc6e73d35c81811b3460c25304dd6e825786d04ef6a94756c7b85680a6fe521` | 88668 | temperature sensitivity |
| `O3_Yoshino1993_295K_185-254nm0.0135nm.txt.gz` | `O3_Yoshino(1993)_295K_185-254nm(0.0135nm).txt` | `03923d989d79c76a07de1ac79de18feb975ec943203a97872b2e781277009cfb` | 5122 | 185-213 nm |
| `O3_Yoshino1993_295K_185-250nmabsolute.txt.gz` | `O3_Yoshino(1993)_295K_185-250nm(absolute).txt` | `9fca86087695082a3f72be167c1dceb98ff80f3faf9c825bc21c12a8231d5500` | 14 | Yoshino's absolute points (three columns: wavelength, cross section, error limit) |
| `O3_Mason1996_298K_110.2-337.4nm.txt.gz` | `O3_Mason(1996)_298K_110.2-337.4nm.txt` | `132e63efde34232331ac23fa30f8b8d16509133b2f7d033953ce3e52de3f0d36` | 351 | 110-185 nm (digitised by the atlas from Fig. 5) |
| `O3_Mason1996_298K_110.2-139.2nm.txt.gz` | `O3_Mason(1996)_298K_110.2-139.2nm.txt` | `1f00a35def670dc230325fdca4b87079dd6e7bf46d08cc0c33314ed5d8d8ef6b` | 309 | the same, 110-140 nm (digitised from Fig. 7) |
| `O3_OgawaCook1958_295K_52-131nm.txt.gz` | `O3_OgawaCook(1958)_295K_52-131nm.txt` | `79997b5be3cd04dafe51d2b6b3e3f0ca9d2ae9c5463a314a6ff25b55086c7465` | 24 | check on the curated `ogawa_1958` Table I (the atlas converts k (atm^-1 cm^-1, base e, 273 K) with 3.72e-20) |
| `O3_Brion1998_295K_195.00-830.00nm.txt.gz` | `O3_Brion(1998)_295K_195.00-830.00nm.txt` | `347e1c262e0ab3a5835504679d4ae8395208fd1552fc1245beb473e130a384ca` | 63501 | validator (BDM, 0.01 nm; J. Brion, personal communication to the atlas, 1998) |
| `O3_Bogumil2003_293K_230-1070nm.txt.gz` | `O3_Bogumil(2003)_293K_230-1070nm.txt` | `8a430cef6df1adb92a0d1248393bc69e8410e327985f0055e5aad5b825a7fd05` | 4244 | validator (SCIAMACHY PFM; vacuum wavelengths; 0.24-0.54 nm resolution; scaled to Bass & Paur) |
| `O3_Burrows1999_293K_230-794nmvac.txt.gz` | `O3_Burrows(1999)_293K_230-794nm(vac).txt` | `deff829c0552b47ef27c0a804d825a0eaa95bc27c7cd9c7e299352e1622db208` | 3338 | validator (GOME FM; vacuum wavelengths) |
| `O3_Daumont1992_295K_195.0-208.3nm.txt.gz` | `O3_Daumont(1992)_295K_195.0-208.3nm.txt` | `beedfe254feb905ef6405777635bfb09161330c8fc05d9969667dc2c4b22bce3` | 403 | validator (digitised by the atlas from Fig. 2) |
| `O3_JPL-20102011_293-298K_121.6-827.5nmrec.txt.gz` | `O3_JPL-2010(2011)_293-298K_121.6-827.5nm(rec).txt` | `d60fb5bf38135e8ba5624f930a8cd03df570407174d982c2912c9eb638b326ce` | 168 | comparison (JPL-2010 recommendation, interval averages) |
| `O3_Ackerman1971_298K_116.5-735nmint-c.txt.gz` | `O3_Ackerman(1971)_298K_116.5-735nm(int-c).txt` | `bb1382f06dbe02c1971400fa0e67c63090d61bafb7582ee52c6b7745f3f45feb` | 162 | comparison (interval averages, 116-730 nm) |

**Atlas comments that matter:**
- **Serdyuchenko 2014:** the IUP Bremen data of 22 July 2013 at 0.01 nm steps; resolution
  0.02-0.06 nm in the UV and visible, 0.12-0.24 nm in the near infrared; eleven
  temperatures, 193-293 K. The atlas does not say whether the wavelengths are vacuum or
  air: see `serdyuchenko_2014` and `gorshelev_2014`.
- **Yoshino 1993:** the 0.0135 nm relative spectrum is put on an absolute scale at 12
  lines (Cu I, Cu II, Hg I), from the CfA molecular database. It supersedes Freeman et
  al. 1984 and Yoshino et al. 1988. The paper (J. Geophys. Res. 98, 5205 (1993)) is not
  in the library.
- **Mason 1996 and Daumont 1992:** digitised from the papers' figures, not tables.

The other temperatures of Serdyuchenko 2014 (193-283 K) are available from the same place
when needed; they are not committed.

## Terms

The atlas asks to be cited (above) and states no licence or terms for reusing its files
(checked 2026-10-16). The data are the original authors' measurements, cited in each
file's source record. The files are redistributed unchanged so that the values optmat
derives from them can be traced (`REUSE.toml`: `LicenseRef-optmat-third-party`).
