# Henke atomic scattering factors

Files `<element>.nff` from the Center for X-Ray Optics, Lawrence Berkeley National
Laboratory, https://henke.lbl.gov/optical_constants/sf/<element>.nff, fetched
2026-09-24 for He, Ne, Ar, Kr, Xe, H, C, N, O, F and S, and 2026-09-28 for the metals
Al, Ag, Au and Cu.

Each file tabulates the photon energy E (eV) and the atomic scattering factors f1 and f2
(electrons per atom) from 10 eV to 30 keV; f1 = -9999 marks energies where no f1 is
given. The tables are those of B. L. Henke, E. M. Gullikson and J. C. Davis, "X-ray
interactions: photoabsorption, scattering, transmission, and reflection at
E = 50-30000 eV, Z = 1-92", At. Data Nucl. Data Tables 54, 181 (1993), as updated by
CXRO. They treat atoms as independent (no chemical or molecular effects) and are not
reliable near absorption edges or below about 30 eV; they are kept for the extreme
ultraviolet (a later phase), where n = 1 - delta - i beta = 1 - r_0 lambda^2 N (f1 + i f2) /
(2 pi) (Henke et al. 1993, Eq. 3, p. 184; source henke_1993).

## Terms

The CXRO pages carry a copyright notice (Center for X-Ray Optics, Lawrence Berkeley
National Laboratory) but state no licence or terms for reusing the tables (checked
2026-10-16). The `.nff` files are redistributed unchanged, with the citation above, so
that the values optmat derives from them can be traced (`REUSE.toml`:
`LicenseRef-optmat-third-party`).
