# Third-order conventions

Read {doc}`/nonlinear` for the scope and first calculation. This chapter defines
Taylor-series gamma, units, frequency arguments and the conversions needed to compare
sources. Continue with {doc}`nonlinear-models`. The API is {mod}`optmat.gas.nonlinear`.

## Conventions

### The Taylor-series gamma

The common currency is Shelton & Rice's Taylor-series gamma (Chem. Rev. 94, 3, 1994),
written $\gamma_T$ or simply $\gamma$. For a static field the dipole moment is

$$
p = \alpha E + \frac{\beta}{2} E^2 + \frac{\gamma}{6} E^3 + \dots
$$

For fields at optical frequencies the same $\gamma$ multiplies the same $1/6$ and the
product of the real fields, with the frequency arguments of the mixing channel. For
an instantaneous (dispersionless) response the induced polarization of a gas of number
density $N$ is simply

$$
P(t) = \frac{N \gamma}{6} E(t)^3 .
$$

Written with complex amplitudes, $E(t) = \tfrac{1}{2}\left[\mathcal{E} e^{-i\omega t} +
\text{c.c.}\right]$ and likewise for $P$, this gives:
- the Kerr (self-phase) amplitude $P_\omega = N\gamma_\text{Kerr}
  |\mathcal{E}|^2\mathcal{E}/8$;
- the third harmonic $P_{3\omega} = N\gamma_\text{THG}\mathcal{E}^3/24$.

### Units

gamma is stored in atomic units:

$$
1\ \text{au} = \frac{e^4 a_0^4}{E_h^3} = 6.235377\times10^{-65}\ \text{C}^4\,\text{m}^4\,
\text{J}^{-3} = 5.0367\times10^{-40}\ \text{esu}
$$

(`optmat.gas.nonlinear.AU_HYPERPOLARIZABILITY`; the esu factor is
`optmat.gas.nonlinear.ESU_PER_SI`, $c^4\times10^{-9}$ in SI numbers, as Shelton &
Rice). Wavenumbers $\nu$ are vacuum wavenumbers in cm$^{-1}$, and wavelengths are vacuum
wavelengths.

### Susceptibility and n2

The quantities a propagation code uses follow from gamma at the gas's number density
$N$, in the dilute limit:
- **The third-order susceptibility** of an instantaneous response,
  $P = \varepsilon_0\chi^{(3)} E^3$, is $\chi^{(3)} = N\gamma/(6\varepsilon_0)$.
- **The per-molecule coefficient** is therefore $\gamma_3 = \gamma/(6\varepsilon_0)$
  in m$^5$/V$^2$, with $\chi^{(3)} = \gamma_3 N$
  ({meth}`~optmat.gas.nonlinear.HyperpolarizabilityModel.gamma3`).
- **The self Kerr coefficient** of linearly polarised light, $n = n_0 + n_2 I$, is
  $n_2 = N\gamma_\text{Kerr}/(8\varepsilon_0^2 c)$
  ({meth}`~optmat.gas.nonlinear.HyperpolarizabilityModel.n2`). One atomic unit of gamma is
  $8.91\times10^{-23}$ cm$^2$/W of $n_2$ at 1 atm and 0 °C.
- **A weak probe** sees twice the self coefficient for the same polarization (cross-phase
  modulation), with the pump-probe gamma.

`n2` uses the real-gas density (the ideal gas for O3, which has no equation of state).
The Lorentz local-field factors, $[(n^2+2)/3]^4 \approx 1 + 8(n-1)/3$, are neglected;
they are 7e-4 for argon at 1 bar.

## The sources' conventions and their conversion

The curated data keep the source's own quantity, convention and conditions; see
{doc}`/the-data/quantities`, *Third-order quantities*. {func}`~optmat.gas.nonlinear.gamma_observation`
converts them on the way into a fit. Every factor was verified on the page images and by
reproducing a printed value (TC-313).

| Source quantity | Convention | To $\gamma_T$ |
|---|---|---|
| gamma in au, SI (C$^4$ m$^4$ J$^{-3}$) or esu | `T` (Taylor) | as is |
| Ward's per-molecule chi in esu (Ward & New 1969, Ward & Miller 1979, Ward & Elliott 1984, Lundeen et al. 1983) | `ward_1969` | $\times 6$ |
| Lehmeier et al.'s gamma, from $P = 4\varepsilon_0\chi E E E$ (SI or esu) | `lehmeier_1985` | $\times 24$ |
| CARS $\chi^{(3)}$ per amagat, esu (Hahn & Lee, Farrow et al.) | `chi_cars_per_amagat` | $2\chi/N_0$, with $N_0$ Loschmidt's number: $10^{-18}$ cm$^3$/(erg amagat) is 147.8 au |
| $n_2$, cm$^2$/W or m$^2$/W, self coefficient | quantity `nonlinear_index` | $8\varepsilon_0^2 c\, n_2/N$ at the dataset's state |
| $n_2$, probe coefficient | `coefficient = "probe"` | halved first |
| $n_2$ with a probe at another wavelength | `parameters.probe_nm` | as the self coefficient, with the pump-probe frequency sums |

The dataset's `component` must be the mean gamma (`mean`, `parallel` or unstated). For
an atom the laboratory `zzzz` component is the mean. An n2's density is the one its
source meant, at the dataset's state (real gas unless the source used the ideal gas).
A dataset option `factor` covers a unit the conversion does not know. Its one use is
Bishop & Pipin's helium, tabulated in the reduced-mass atomic unit of 4He: 6.24136e-65
against 6.235377e-65, a factor of 1.00096 (their p. 3549).

Other conventions in the literature, recorded and converted in the curated data or its
notes (TC-313), include:
- Rosasco's $\sigma = N\gamma_T/3$;
- Hellwarth's $\sigma$ in esu;
- Bishop's dispersion coefficients A and A', which map to Shelton's A and B through the
  square and fourth power of 219474.63 cm$^{-1}$ per hartree.

## Processes and their frequencies

A process is named by its frequency arguments $(-\omega_\sigma; \omega_1, \omega_2,
\omega_3)$, with $\omega_\sigma = \omega_1 + \omega_2 + \omega_3$. The fundamental (or
pump) wavenumber is $\nu$; $\nu_s$ is a Stokes and $\nu_r$ a probe wavenumber. The model
uses the two symmetric sums

$$
W_2 = \nu_\sigma^2 + \nu_1^2 + \nu_2^2 + \nu_3^2, \qquad
W_4 = \nu_\sigma^4 + \nu_1^4 + \nu_2^4 + \nu_3^4 ,
$$

({func}`~optmat.gas.nonlinear.frequency_sums`). A negative frequency enters as its
square, so its sign does not matter.

| Process (`process`) | Arguments | $W_2$ | $W_4$ | Where it enters |
|---|---|---|---|---|
| `static` | $(0; 0, 0, 0)$ | 0 | 0 | the static limit, $\gamma_0$ |
| `kerr` | $(-\omega; \omega, \omega, -\omega)$ | $4\nu^2$ | $4\nu^4$ | $n_2$, SPM; DFWM |
| `pump_probe` | $(-\omega_r; \omega_r, \omega_p, -\omega_p)$ | $2\nu^2 + 2\nu_r^2$ | $2\nu^4 + 2\nu_r^4$ | XPM; a probe's index |
| `thg` | $(-3\omega; \omega, \omega, \omega)$ | $12\nu^2$ | $84\nu^4$ | the third harmonic |
| `eshg` | $(-2\omega; \omega, \omega, 0)$ | $6\nu^2$ | $18\nu^4$ | the ESHG data |
| `dc_kerr` | $(-\omega; \omega, 0, 0)$ | $2\nu^2$ | $2\nu^4$ | the dc Kerr data |
| `cars` | $(-\omega_{as}; \omega_p, \omega_p, -\omega_s)$ | $(2\nu-\nu_s)^2 + 2\nu^2 + \nu_s^2$ | likewise | CARS, four-wave mixing |

The wavelength of a stored value is always the fundamental (input) wavelength, or the
pump for CARS and pump-probe.

## Parts: electronic, vibrational, rotational

A measured third-order response contains:
- the electronic gamma, $\gamma^e$, instantaneous on the scale of the pulse;
- a vibrational part, $\gamma^v$;
- for a molecule that is not spherical, the rotational (orientational) part.

Which of these a measurement contains depends on the process and on the pulse. The
data record it (`part`: `electronic`, `vibrational`, `rotational`, `total`,
`instantaneous` and so on). The model is electronic:
- **H2, D2, N2, O2:** Shelton's ESHG includes a vibrational and rotational part,
  $\gamma^{vR}$. It is subtracted before fitting: for H2 and D2 from Shelton's tabulated
  values at the same wavelengths (option `subtract_dataset`), for N2 and O2 from his
  formula $\gamma^{vR} = G\nu^{-2} + H\nu^{-4}$ (Table V; option `subtract_formula`).
- **CO2, CH4, SF6:** the fits are to the total ESHG gamma. Its vibrational part is
  small but not separated ("we are unable to accurately determine $\langle\gamma^e
  \rangle$ or its static limit for these molecules", Shelton 1990, p. 13). It is
  process-dependent, which the model cannot represent, and it is in the 2 % floor of
  the accuracy page.
- **Atoms:** gamma is purely electronic.
- **Femtosecond n2 of molecules:** the rotational response is delayed and depends on the
  pulse. The Maryland group separate it and report the electronic n2; values called
  "instantaneous" or measured with longer pulses contain part or all of it, and are 1.3
  to 6 times higher. In propagation the rotational and vibrational Raman responses are
  separate terms, with their own data. The gamma here is the instantaneous electronic
  term they add to.
- **N2O:** its only value, Wahlstrand et al.'s n2, contains the vibrational part.
