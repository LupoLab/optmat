# The Raman response

Read {doc}`raman` for the line parameters and {doc}`/user-guide/nonlinear-conventions`
for the electronic/delayed split. This chapter gives the tensor kernel and its
propagation interfaces; {doc}`/user-guide/raman-validation` gives the comparisons
and open questions. The {doc}`quickstart` is the introductory route.

The delayed third-order response of the molecular gases and of fused silica, as a
function of the delay and of the frequency, ready for a propagation code to sample:
{mod}`optmat.gas.raman_response` and {mod}`optmat.solid.raman_response`. The gases'
responses are built from the lines of {doc}`raman`; silica's from Hollenbeck &
Cantrell's modes and the Raman fraction of the Kerr fit of {doc}`silica`. The checks
against measurements are in {mod}`optmat.gas.raman_response_validation` and
{mod}`optmat.solid.raman_response_validation`.

```python
import numpy as np

from optmat.gas.linewidths import AMAGAT
from optmat.gas.raman_response import raman_response

t = np.linspace(0.0, 10e-12, 1001)  # delay, s
omega = 2 * np.pi * np.linspace(0.1e12, 5e12, 50)  # rad/s
densities = np.array([0.5, 1.0, 2.0]) * AMAGAT  # m^-3

air = raman_response("dry air", 295.0, wavelength=800e-9)
air.response(t, 1.0 * AMAGAT)  # h(t), m^5 V^-2 s^-1, t in s
air.spectrum(omega, 1.0 * AMAGAT, "b")  # H(omega), or one of its parts
air.decay_rates(densities), air.frequencies(densities)  # (lines, densities), cheap
air.frequency, air.isotropic, air.anisotropic, air.damping_group  # the lines
air.density_range, air.sources, air.reasons

from optmat.solid.raman_response import silica_raman_response

silica = silica_raman_response()  # Hollenbeck & Cantrell, the Kerr fit's f_R
silica.response(t), silica.spectrum(omega)  # normalised: integral 1, H(0) = 1
silica.raman_coefficient(1030e-9), silica.electronic_susceptibility(1030e-9)
shift = 2 * np.pi * 13.2e12  # the gain peak, 440 cm^-1, rad/s
silica.gain(1030e-9, shift)  # m/W
```

## The convention

For a field linearly polarised along $x$, the Raman polarisation of a gas of number
density $N$ is (TC-442, TC-685)

$$
P_x(t) = \varepsilon_0 N E_x(t) \int_0^\infty h(\tau)\, E_x^2(t - \tau)\, d\tau,
\qquad
h(\tau) = \sum_k A_k \sin(\Omega_k \tau)\, e^{-\Gamma_k \tau},
$$

with $A_k$ in m⁵ V⁻² s⁻¹ per molecule of the gas (a mixture's lines carry their
constituent's mole fraction) and $\Gamma_k = \pi\,\Delta\nu_k$, $\Delta\nu_k$ the FWHM
in Hz. The transform $H(\omega) = \int_0^\infty h(\tau)\, e^{i\omega\tau} d\tau$ is
analytic, $\sum_k A_k\Omega_k / (\Omega_k^2 + \Gamma_k^2 - \omega^2 - 2i\Gamma_k\omega)$,
and its imaginary part is the Stokes gain: $g = \omega_s N \operatorname{Im}H(\Omega) /
(c^2\varepsilon_0)$ for pump and Stokes polarised alike.

**The tensor form.** For any polarisation, following Hellwarth,

$$
P_i(t) = \varepsilon_0 N \int_0^\infty R_{ijkl}(\tau)\, E_j(t)\, E_k(t-\tau)\, E_l(t-\tau)\, d\tau,
\qquad
R_{ijkl} = a(\tau)\,\delta_{ij}\delta_{kl}
  + \tfrac12\, b(\tau)\,(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}),
$$

and $h = a + b$. Each line has an isotropic share, from the mean polarizability, which
enters $a$ only, and an anisotropic share $h_2$, from the anisotropy, for which
$a = -h_2/2$ and $b = 3h_2/2$: the change of the anisotropic polarizability follows the
traceless part of $\mathbf{E}\mathbf{E}$. A weak probe polarised across the pump sees
$a$: $-1/2$ of the parallel rotational response, and all of the isotropic vibrational
one ({meth}`~optmat.gas.raman_response.GasRamanResponse.component_amplitudes`).

- The rotational S lines are anisotropic.
- A vibrational Q(J) line's anisotropic share is $(4/45)\, b_{J,J}\, \gamma'^2 P_J /
  (\hbar\varepsilon_0)$, with the Placzek-Teller coefficient $b_{J,J} = J(J+1) /
  [(2J-1)(2J+3)]$ ({func}`~optmat.gas.raman_response.placzek_teller_q`, checked by
  quadrature). The transition anisotropy $\gamma'$ is Hunt et al.'s per J for H2 and
  D2, and elsewhere follows from the Q branch's depolarisation ratio
  $\rho = 3b\gamma'^2 / (45\bar\alpha'^2 + 4b\gamma'^2)$ (TC-686).
- The vibrational O and S branches ($\Delta J = \mp 2$), anisotropic only, are left
  out: they would add 8 % of N2's Q branch, 19 % of O2's and 3 % of H2's and D2's,
  spread over some 100 cm⁻¹ about the band.

## The gases

The lines are those of {doc}`raman`: per J, with their own strengths and widths, the Q
branches resolved (TC-541, TC-542). Every curated width law has the form
$\Delta\nu = c + d/n + \beta n$ at a given temperature and composition ($n$ in
amagat), so the response keeps the three coefficients of each line and its density
dependence costs nothing to evaluate (TC-687). The amplitudes depend on the
temperature (the populations) and the pump wavelength (the vibrational strengths'
dispersion, TC-542) only.

<!-- response-terms:start -->
| Gas | S lines | Q lines | Damping groups | Densities (amagat) | Temperatures (K) | Anisotropic Q share |
|---|---|---|---|---|---|---|
| N2 | 53 | 53 | 43 | 0.054-2.09 | 293-300 | 2.9% |
| O2 | 32 | 32 | 38 | 0.0345-4.2 | 290-1322 | 6.7% |
| H2 | 10 | 10 | 3 | 0.943-3.29e+03 | 77-300 | 2.0% |
| D2 | 15 | 15 | 12 | 0.616-1.02e+03 | 291-298 | 1.5% |
| CH4 | 0 | 1 | 1 | no limit | 290-300 | 0.0% |
| SF6 | 0 | 1 | 1 | no limit | 290-300 | 0.0% |
| CO2 | 61 | 2 | 8 | 0.0108-13.1 | 283-300 | 4.7% |
| N2O | 0 | 2 | 2 | above 0.0205 | 290-300 | 0.0% |
| dry air | 85 | 85 | 82 | 0.054-2.09 | 293-300 | 3.6% |
<!-- response-terms:end -->

- **Damping groups.** Lines with the same group share their width and shift at every
  density: a propagation code can sum them into one damped kernel with no
  approximation. H2's Q lines are one group (one width law for every Q(J)); N2's S
  lines outside Herring et al.'s S(2)-S(14) share the end values' widths.
- **The composition.** A mixture's widths use each constituent's true collision
  partners: N2's lines in air are broadened by N2, O2 and (as O2) Ar, not by "the total
  density". Dry air is CIPM-2007's N2, O2 and Ar, renormalised: CO2 would add 0.3 % to
  the rotational response (TC-689).
- **The densities** (TC-688) between which the lines are independent collisional
  Lorentzians, found for each band's strongest line:
  - *below*, molecular motion sets the width. Where a line has a Dicke term (H2, D2),
    the diffusion law holds above Herring et al.'s cutoff $\rho_c = 3.33\, d /
    \Delta\nu_D$, with $\Delta\nu_D$ the Doppler width for forward scattering (their
    Eqs. (3) and (4)); below it, $d/n$ is capped at $\Delta\nu_D$, a Lorentzian stand-in
    for the Doppler limit. Otherwise the limit is where the collisional width falls to
    $\Delta\nu_D$;
  - *above*, the width of the strongest line reaches the spacing to its neighbour: the
    lines overlap, and their mixing (the collisional narrowing of a Q branch) is not
    modelled (TC-527). For N2 this is 2.1 amagat, where Koszykowski et al. see the
    branch begin to collapse.
- **Line shifts.** H2's Q lines move by $-96$ MHz/amagat at 298 K (Bischel & Dyer's
  Eq. (4)), 1.85 times the collisional width's coefficient, and the response follows
  them ({meth}`~optmat.gas.raman_response.GasRamanResponse.frequencies`). The other
  curated shifts are at most 5 % of their widths' coefficients and are not used.
- **The temperatures** are those of each band's width measurements; outside them the
  response warns when it is built.

## Fused silica

A glass's response is a normalised shape and a size: the Raman share $f_R$ of the
long-pulse $n_2$, as the Kerr fit defines it (TC-393). For an envelope,
$\Delta n = n_2[(1-f_R) I + f_R \int h_R(\tau) I(t-\tau)\, d\tau]$. For a real field,
in the gases' convention,

$$
P = \varepsilon_0\chi_e E^3 + \varepsilon_0 K E \int_0^\infty h_R(\tau) E^2(t-\tau)\, d\tau,
\qquad
\chi_e = \tfrac43 n_0^2 \varepsilon_0 c\,(1-f_R)\, n_2,
\quad
K = 2 n_0^2 \varepsilon_0 c\, f_R\, n_2,
$$

because a kernel slower than the optical cycle follows only the slowly varying half of
$E^2$, while the instantaneous term also has its $3\omega$ part. The Raman term of the
field therefore carries 2/3 of the envelope's share: Blow & Wood's $\alpha = 0.7$ on
$E^2$ is a Raman share of 0.222 of $n_2$ (their Eqs. (2) and (8); TC-692). The gain is
$g(\Omega) = 2\omega_s n_2 f_R \operatorname{Im}H_R(\Omega)/c$.

<!-- silica-models:start -->
| Model | f_R | Paired with the Kerr fit | Isotropic and anisotropic parts | Gain peak (cm-1) | Peak gain at 526 nm (1e-11 cm/W) | at 795.5 nm | at 1064 nm |
|---|---|---|---|---|---|---|---|
| hollenbeck_2002 | 0.219 | yes | no | 440 | 1.91 | 1.15 | 0.83 |
| lin_2006 | 0.245 | no | yes | 435 | 1.91 | 1.16 | 0.83 |
| blow_1989 | 0.222 | no | no | 436 | 2.10 | 1.27 | 0.91 |
<!-- silica-models:end -->

- **`hollenbeck_2002`, the default.** Hollenbeck & Cantrell's 13 modes, fitted to both
  the measured gain spectrum and the measured response. Table 1 is used as printed with
  its footnote's conversions, and normalised analytically (TC-382, TC-691). Paired
  with the Kerr fit: its $f_R$ is the one that fit found with this shape, so the
  electronic part and the Raman part add up to the fitted $n_2$. It is copolarised
  only.
- **`lin_2006`.** Lin & Agrawal's isotropic damped oscillator and anisotropic Boson
  peak: the one model with both parts of the tensor. With their $n_2 = 2.6\times
  10^{-20}$ m²/W and $f_R = 0.245$ it reproduces their peak gains, 1.2e-11 cm/W at 795.5
  nm and 1.81e-11 cm/W at 526 nm, to 0.3 % and 0.6 %, once their $\omega_0$ (the
  pump) stands in for the Stokes frequency. Their $f_R$ is not the Kerr fit's.
- **`blow_1989`.** The single oscillator, 12.2 and 32 fs, with its share 0.222.

## Using it in a propagation code

A propagation code samples the response once, when it sets up a run, as it builds its
ionisation rate tables. The responses here are functions of continuous delay and
frequency. The grid, the taper, the padding and any tabulation belong to the code.

- **The status first.** `raman_response_status(gas)` gives the status and reason the
  response would have, from the curated choices alone, without computing the lines.
  It returns None where no constituent is Raman-active: the noble gases, alone or
  mixed. A molecule without curated choices, and CIPM-2007's full air, are refused, as
  `raman_response` refuses them. Collision partners are not checked: a mixture of CO2
  with N2, whose broadening coefficients by one another are not curated, has a status
  here and is refused by `raman_response` (TC-840).

- **The cost.** The response of a gas is a sum of lines, and a code that convolves by
  FFT pays for them only when it builds the kernel's spectrum: once for a constant
  density. Along a density gradient each line's decay rate changes, and rebuilding the
  kernel costs one exponential per line and sample. With per-J widths the lines are
  many: 106 for N2 and 170 for air. At $N = 2^{14}$ that is
  about 7 ms per rebuild for N2 against 0.2 ms for one convolution (numpy, one core).
  The convolution itself costs the same whatever the number of lines.
- **Grouping.** Lines of one damping group can be summed before damping, exactly.
- **Low rank.** Over the densities a run visits, the family of kernels is close to a
  sum of a few separable terms, $H(\omega; n) \approx \sum_r c_r(n)\, U_r(\omega)$.
  With the $U_r$ precomputed, a rebuild is $r$ multiply-adds per sample, with no
  exponentials and no FFT. The table gives the smallest $r$ at which the largest error
  on a 64-density grid, over the kernel's spectrum on a propagation-like window, falls
  below each fraction of the kernel's largest value:

<!-- response-rank:start -->
| Gas | Densities (amagat) | Window | Lines | Damping groups | Rank at 1e-03 | Rank at 1e-04 | Rank at 1e-06 |
|---|---|---|---|---|---|---|---|
| N2 | 0.5-2 | 131 ps at 4 fs | 106 | 43 | 3 | 4 | 5 |
| N2 | 0.1-10 | 131 ps at 4 fs | 106 | 43 | 5 | 6 | 9 |
| O2 | 0.1-10 | 131 ps at 4 fs | 64 | 38 | 5 | 6 | 8 |
| dry air | 0.1-10 | 131 ps at 4 fs | 170 | 82 | 5 | 6 | 8 |
| H2 | 1-30 | 66 ps at 2 fs | 20 | 3 | 4 | 5 | 6 |
| D2 | 1-30 | 66 ps at 2 fs | 30 | 12 | 4 | 4 | 5 |
<!-- response-rank:end -->

### The damping groups

`response.groups()` gives each damping group once, in SI, with the conversions written
in optmat (OEP 1, sections 4.5 and 6.7). A group is a set of lines that share their
width law and shift exactly, $(c, d, \beta, \Delta\nu_D, s)$, and its share of the
response is

$$
R_g(t; N) = e^{-\Gamma_g(N)\,t} \sum_{k \in g} A_k \sin\big((\omega_k + \sigma_g N)\,t\big),
\qquad
\Gamma_g(N) = \gamma_c + \min(\gamma_d/N,\ \gamma_D) + \gamma_\beta N,
$$

with $N$ the total number density in m⁻³, $\gamma_c = \pi c$, $\gamma_d = \pi d
N_\text{amg}$, $\gamma_D = \pi \Delta\nu_D$, $\gamma_\beta = \pi\beta/N_\text{amg}$ and
$\sigma_g = 2\pi s/N_\text{amg}$. The groups' responses sum to `response` exactly.
Where $\sigma_g = 0$ the row of amplitudes and frequencies is fixed and only the decay
depends on the density. The one group with a shift is H2's Q branch
(`group.shifted`), so a response of H2 along a density gradient is exact only at a
fixed density.

`decay_rate(N)` is defined for $N \ge 0$. At $N = 0$ it is the law's limit: the capped
Dicke term tends to its cap, so $\Gamma_g(0) = \gamma_c + \gamma_D$ for a group with a
Dicke term and $\gamma_c$ otherwise. At 295 K that limit is zero for every line of N2,
O2, dry air, N2O and CO2 and some of D2's, and finite for H2 (at least 153 MHz in
$\Gamma/\pi$), CH4 and SF6. The group view does not warn below the validated
`density_range`; a code declares the densities it will use once, with
`check_densities(N)`, which does.

Each line has a stable identifier, `line_id`: `species/band/branch/J`, with the band's
short curated key (`rotation`, `v1`, `nu1`, `nu3`, `fermi_1285`, `fermi_1388`) and `-`
for an unresolved band. The identifiers are unique and are kept in a snapshot
(`tests/data/raman_line_ids.json`), so one changes only deliberately.

```python
h2 = raman_response("H2", 295.0)
q = next(g for g in h2.groups() if g.shifted)
assert q.line_ids[0] == "H2/v1/Q/0"
rates = q.decay_rate(np.array([0.0, 1.0]) * AMAGAT)  # 1/s
assert rates[0] > 0
h2.check_densities(np.array([1.0, 10.0]) * AMAGAT)
```

### The two tensor parts

A code working in a trace and traceless basis needs each part written out. The
isotropic part is $R^{(0)}_{ijkl} = h_0\,\delta_{ij}\delta_{kl}$. The anisotropic part is
$R^{(2)}_{ijkl} = h_2\left[\tfrac34(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}) -
\tfrac12\delta_{ij}\delta_{kl}\right]$, which is $a = -h_2/2$, $b = 3h_2/2$ in the
convention above. The two elements tell them apart:

- **copolarised**, a probe along $x$ after a pump along $x$: $R_{xxxx} = a + b =
  h_0 + h_2$ (`component="parallel"`);
- **cross-polarised**, a probe along $x$ after a pump along $y$: $R_{xxyy} = a = h_0 -
  \tfrac12 h_2$ (`component="a"`).

The weights in the transverse trace and traceless basis, $(1, 0)$ for an
isotropic row and $(\tfrac14, \tfrac34)$ for a rank-two row, are this decomposition for
a transverse drive $Q$ in three dimensions. $\tfrac12 h_2 (3Q - \operatorname{tr}Q\, I)$
has trace part $\tfrac14 h_2 \operatorname{tr}Q$ and traceless part $\tfrac32 h_2
(Q - \tfrac12\operatorname{tr}Q\, I)$, which gives $h_2(\tfrac14 + \tfrac34) = h_2$
copolarised and $h_2(\tfrac14 - \tfrac34) = -\tfrac12 h_2$ cross-polarised.

```python
iso, aniso = air.isotropic, air.anisotropic
assert np.allclose(
    air.component_amplitudes("parallel"), (iso + aniso / 4) + 3 * aniso / 4
)
assert np.allclose(air.component_amplitudes("a"), (iso + aniso / 4) - 3 * aniso / 4)
```

<span id="validation"></span>
<span id="open-questions"></span>
<span id="against-solprop"></span>

See {doc}`/user-guide/raman-validation` for comparisons with measured quantities.
