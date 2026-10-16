# The extreme ultraviolet

Read {doc}`conventions` and {doc}`models` for the polarizability and density
conventions. This chapter describes the absorbing response above the first
resonance; the {doc}`quickstart` introduces the below-resonance models. See
{ref}`helium <sphx_glr_auto_gallery_dispersion_plot_euv_he.py>` and
{ref}`heavier gases <sphx_glr_auto_gallery_dispersion_plot_euv_heavy.py>` in the gallery,
{mod}`optmat.gas.euv` for the interface, and {doc}`/the-data/sources` for the records.

Above its first resonance a gas absorbs, and its refractive index is complex. The package
builds it from the gas's dipole oscillator-strength distribution (`optmat.gas.euv`): the discrete
lines from the NIST database and the photoabsorption continuum from the curated
cross sections, joined to Henke's scattering factors at high energy. Both parts of the
polarizability follow:

$$
\mathrm{Re}\,\alpha(\omega) = \sum_i \frac{f_i}{E_i^2-\omega^2}
+ \mathcal{P}\!\int_I^\infty \frac{g(E)\,\mathrm{d}E}{E^2-\omega^2},\qquad
\mathrm{Im}\,\alpha(\omega) = \frac{\pi\,g(\omega)}{2\omega},
$$

with $g = \mathrm{d}f/\mathrm{d}E$ the continuum density (the photoionisation cross
section is $109.76$ Mb eV $\times\,g$). The same distribution gives the dipole sums
$S(k)$, which test it: the Thomas-Reiche-Kuhn sum $S(0)$ is the number of electrons, and
$S(-2)$ the static polarizability. The measured continuum is scaled, within its stated
uncertainty, so that $S(-2)$ equals the ab initio value (a constrained distribution in the
sense of Kumar and Meath); the higher-energy data are left alone and $S(0)$ is a check.
The index is $n = \sqrt{1+\chi}$ with the Lorentz-Lorenz $\chi$ of the complex
polarizability, so $\mathrm{Re}\,n - 1$ turns negative above threshold and
$\mathrm{Im}\,n$ is the absorption index.

A distribution is configured in `data/gas/euv/<gas>.toml`: the continuum sources, the
Rydberg continuation of the lines, and the constraint.

The continuum constraint uses a bounded local secant iteration for small adjustments
within the measured uncertainty. If the sum does not respond to scaling, or a step
would make the continuum scale non-positive or non-finite, `constrained` raises
`DataError` (TC-736). A large change can require a different root solver; an invalid
trial step does not establish that the requested target is physically impossible.

## Helium

NIST's $1snp$ lines ($n = 2$-$10$, continued by the $n^{*-3}$ law to $n = 60$), Samson et
al.'s absolute photoionisation cross section from the threshold (24.587 eV) to 120 eV
(1-2 %), and Henke's $f_2$ from 120 eV to 30 keV. The measured continuum is scaled by
1.012 to Puchalski et al.'s static polarizability. Checks:

- $S(0) = 1.979$ against 2 (the missing 1 % is at high energy);
- below threshold the transform of the photoabsorption reproduces the ab initio
  refractive polarizability (`fit/He`) to 0.01 % in the infrared and 0.1 % at 70 nm;
- the absorption agrees with Liu's ab initio complex polarizability (27-58 eV) to 1 %;
- $f_1$ agrees with Henke's to 1 % above 70 eV; below about 50 eV Henke's is not reliable
  (as Henke et al. state).

Between the first line (58.4 nm) and the threshold (50.4 nm) the Rydberg lines make the
index singular; the model gives it between the lines but carries no line widths, so it is
not an absorption model there.

## Neon

NIST's lines with oscillator strengths (the 3s-6s and 3d-5d members of both series,
0.246 in all), the unlisted higher members as one pseudo-line at 21.55 eV (0.046: Chan et
al.'s 1992 discrete total to 21.616 eV, the midpoint of the two thresholds, less the NIST
lines), Chan et al.'s continuum from 21.6 to 250 eV and Henke's $f_2$ above. Chan's neon
continuum is normalised to the Thomas-Reiche-Kuhn sum rather than absolute; the
constraint to Lesiuk et al.'s static polarizability scales it by 1.054. Checks:

- $S(0) = 9.72$ against 10;
- below the first line the transform agrees with `fit/Ne` to 0.02 % in the visible and
  infrared, 0.34 % at 120 nm (the fit's own uncertainty there is 1.5 %);
- $S(-4) = 2.826$ and $S(-6) = 4.765$, 0.8 % and 3.4 % below Lesiuk et al.'s ab initio
  2.850(7) and 4.932(14). The photoabsorption gives *less* dispersion than theory, where
  the refractometry gives more (TC-92, TC-126).

## Argon, krypton and xenon

Chan et al.'s (1992) absolute dipole (e,e) oscillator strengths: their discrete lines
below the $^2P_{3/2}$ limit (with the unlisted higher members, the difference between their
"total to ionization" and the listed lines, as one pseudo-line just below the limit), their
continuum (5 %) from above the $^2P_{1/2}$ limit to 380-500 eV, and Henke's $f_2$ above.
The autoionising region between the two limits is not tabulated in the paper (only
plotted); its density is taken as the first tabulated value down to the lower limit, an
approximation the sum rules test. The measured continuum is then scaled to the static
polarizability (Lesiuk & Jeziorski for argon, Huot & Bose for krypton and xenon):

| Gas | Scale | $S(0)$ (TRK) | Against the fit, 0.5-2 µm | Note |
|---|---|---|---|---|
| Ar | 0.981 | 18.51 (18) | 0.003 % | $S(-4) = 27.933$ against Lesiuk's ab initio 27.976(15), 2.9 sigma: with the authors' 1993 erratum, the 7s/5d pair's 0.0426 is counted once and the pseudo-line takes the rest (TC-710) |
| Kr | 0.994 | 34.51 (36) | 0.09-0.13 % | |
| Xe | 1.0005 | 50.85 (54) | 0.004-0.05 % | |

Without the autoionising-region fill the scales would be 1.006, 1.10 and 1.30: xenon's
autoionising region carries a large share of its near-threshold strength.

Below the first resonance the fitted models are to be preferred: Chan's line energies are
printed to 1 meV, which near a resonance moves the polarizability by about 1 %.
