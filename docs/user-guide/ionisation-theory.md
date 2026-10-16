# Noble-gas ionisation: rate theory

Read the {doc}`/ionisation_rates` example and {doc}`/ionisation` parameters first. This chapter derives the rate, its cycle average and energy loss, and explains the omitted terms. Validation begins with {doc}`ionisation-static`. The callable interface is {mod}`optmat.gas.ionisation_rates`.

## The model

Atomic units: $\kappa = \sqrt{2 I_p}$, $n^* = Z/\kappa$, $F_0 = \kappa^3$, the Keldysh
parameter $\gamma = \omega\kappa/F$ and the barrier-suppression field $F_b = I_p^2/4Z$.

$$
W(F, \omega) = N_\text{eff} \sum_c p_c\, w_c(F, \omega)\, S_c(F),
\qquad
w = Q\, w^\text{sr} .
$$

| Factor | What it is | From | Register |
|---|---|---|---|
| $w^\text{sr}$ | the cycle-averaged rate of one electron in the orbital $(l, m = 0)$, in a short-range potential, at any $\gamma$ | Perelomov, Popov & Terent'ev 1966, eq. (54) | TC-656 |
| $Q$ | the Coulomb factor: the static factor $(2F_0/F)^{2n^*}$ times $(1 + 2\gamma/e)^{-2n^* c s}$, $s = n^*\omega/2I_p$, $c = 5.25$ | Perelomov & Popov 1967, eq. (42), as in Larochelle et al. 1998, eqs. (9)-(10); the second factor fitted here to one-electron calculations | TC-653, TC-670, TC-831, TC-848 |
| $S$ | the correction of the static rate at finite field, $\exp[-\alpha (Z^2/I_p)(F/\kappa^3)]$ | Tong & Lin 2005, eq. (2) | TC-655 |
| $C_l$ | the asymptotic coefficient of the orbital | Tong, Zhao & Lin 2002, Table II | TC-657 |
| $N_\text{eff}$ | the rate of the atom over that of the one electron | set here from all-electron calculations and measurements | TC-658, TC-846, TC-847 |
| $p_c$ | the weights 2/3 and 1/3 of the ion's levels $^2P_{3/2}$ and $^2P_{1/2}$, each with its own $I_p$ | the Clebsch-Gordan coefficients; levels from NIST ASD | TC-659 |

**The short-range rate** is

$$
w^\text{sr}_{lm} = I_p\, |C_{\kappa l}|^2 \sqrt{\frac{6}{\pi}}\, f(l, m)
  \left(\frac{F\sqrt{1 + \gamma^2}}{2F_0}\right)^{|m| + 3/2}
  A_m(\omega, \gamma)\,
  \exp\left[-\frac{2F_0}{3F}\, g(\gamma)\right],
$$

with $g$ of their eq. (33), $A_m$ of eq. (55) and $f(l, m) = (2l+1)(l+|m|)!/(2^{|m|}
|m|!(l-|m|)!)$. For $\gamma \to 0$, $g \to 1$ and $A_m \to 1$, and $Q\,w^\text{sr}$ is
the cycle average of the static tunnelling rate (their eqs. 4 and 7),

$$
w_\text{stat} = I_p\, |C_{\kappa l}|^2 f(l, m)
  \left(\frac{2F_0}{F}\right)^{2n^* - |m| - 1} e^{-2F_0/3F},
$$

which is also Tong, Zhao & Lin's eq. (7) and Ammosov, Delone & Krainov's eq. (1) for
$m = 0$. Tong et al.'s coefficient is $C_l = C_{\kappa l}\,\kappa^{n^* + 1/2}$.

**The sum over channels** in $A_m$ has a step wherever the number of photons needed,
$\nu = (I_p + U_p)/\omega$, passes an integer. By default the sum is replaced by its
integral over the number of photons,

$$
\sum_{n \ge \nu} e^{-\alpha(n - \nu)}\, w_m\!\left(\sqrt{\beta(n - \nu)}\right)
\;\to\; \frac{\sqrt\pi}{2}\, |m|!\, \frac{\beta^{|m|}}{(\alpha + \beta)^{|m|+1}}
\sqrt{\frac{\beta}{\alpha}},
$$

which is smooth in the field, as a quadrature over a beam profile or an interpolation
table needs. It is within 10 % of the sum for $\gamma < 1.5$ and 50 % up to $\gamma =
3$ (`channel_sum="sum"` gives the sum; TC-654). The integral stays the default; for
any model the sum is
`dataclasses.replace(ionisation_rate_model("Ar", 193e-9), channel_sum="sum")`, the
better choice deep in the multiphoton regime, where the integral gives the yield a
non-integer order ({doc}`ionisation-measured`; TC-856).

## Three parts, each tested against what isolates it

The rate of the atom is the rate of one electron times the number of electrons that
behave like it, and the rate of one electron is a static rate times what the field's
frequency does to it. The three parts are set and tested separately, each against the
references that isolate it, so that an error in one is not absorbed into another:

| Part | What sets it | Tested against | Per gas or global |
|---|---|---|---|
| **The Coulomb factor's dependence on the frequency** ($c$ in $Q$) | one electron's physics, the same form for every atom | one-electron time-dependent Schrödinger calculations, which isolate it exactly: fitted to argon at 400, 800 and 1030 nm (10 fs), tested on argon in other pulses (6, 30 fs, Wallace et al.'s 5.5 fs) and on hydrogen, helium, neon, krypton and xenon | one constant for all gases |
| **The static normalisation** ($C_l$ and $\alpha$ in $S$) | the static rate of one electron | exact static rates of one electron in model potentials (tRecX, Bahl et al., Zhang et al.) | per gas |
| **The effective number of electrons** ($N_\text{eff}$) | the atom's many-electron physics | all-electron static rates (Jagau; Parker for helium), all-electron calculations in a pulse (haCC; R-matrix), and calibrated measured yields (Wallace et al.: krypton and xenon) | per gas |

Why they are kept apart:

- **The Coulomb factor is fitted to one-electron calculations only.** A one-electron
  model compared with a one-electron calculation contains no $N_\text{eff}$, so nothing
  about the atom enters the fit. Fitting $c$ per gas against measured yields would
  absorb errors of $N_\text{eff}$ and $C_l$ into the Coulomb factor and hide them:
  that is what the static factor did. It is 1.1 to 4 times too high above $\gamma = 1$
  for one electron, and this excess made up for krypton's missing electrons, so that
  the old default matched Wallace et al.'s krypton with $N_\text{eff} = 1$ (TC-846).
- **The frequency-dependent factor tends to the static one as $\omega \to 0$,** so the
  static references, and the values of $N_\text{eff}$ set from them, do not depend on
  it.
- **$N_\text{eff}$ is set with the Coulomb factor fixed.** In a pulse, argon ionises
  like about 0.8 to 0.9 of one $m = 0$ electron (haCC at 800 nm, R-matrix at 390 nm),
  krypton like about 1.4 (Wallace et al.'s yields on the exact calculation; TC-846;
  0.91 in Pathak et al.'s time-dependent CASSCF at 2e14 W/cm2, TC-860),
  helium like about 1.5 (Parker et al.), and xenon like about 1.1 (Wallace et al.'s
  yields on the exact calculation; TC-849). These are ratios to one electron of a
  model potential, without the spin-orbit splitting. $N_\text{eff}$ itself is the
  atom's rate over *the model's own* one-electron rate, with its coefficient $C_l$ and
  the ion's two levels weighted 2/3 and 1/3: argon's against Jagau's all-electron
  rates, helium's against Parker et al.'s, krypton's and xenon's against Wallace et
  al.'s yields, which test the rate only through the shape of the saturation and so
  weakly (a factor of 0.8 moves the scale by about 3 %). The weighting
  makes the model's rate 0.71 to 0.95 times that of one level (xenon to argon), so a
  ratio to an exact one-electron calculation is divided by it, and by the model's
  ratio to the calculation, before it is compared with $N_\text{eff}$.
- **The static normalisation is checked in a static field,** where neither of the
  other two enters: argon's and xenon's are within 5 to 15 %, krypton's 1.08 to 1.17
  high.

The default before TC-848, with the static Coulomb factor, is kept as `default_static`
(and `one_electron_static`): the comparisons that the register made with it can be
reproduced, and both appear in the validation tables.

## Two rates

| Method | Argument | Use |
|---|---|---|
| `cycle_averaged_rate` | the amplitude of the field | an envelope calculation; this is the rate the theory defines |
| `instantaneous_rate` | the field itself | a calculation that resolves the carrier |

The two are made consistent: **the average of the instantaneous rate over a cycle is
the cycle-averaged rate**, to 4 % at 400 nm and 1.5 % at 800 nm up to the model's
maximum field. A pulse then ionises the same fraction whichever is used.

The usual practice is to divide the cycle-averaged rate by $(3F/\pi F_0)^{1/2}$, the
factor of the cycle average in the tunnelling limit (Perelomov et al., eq. 7), and to
evaluate the result at the instantaneous field. That factor holds for $\gamma \to 0$
only. At finite $\gamma$ the rate is less steep in the field than the tunnelling
exponential, a larger part of the cycle contributes, and the rate so obtained ionises
more over a cycle than the rate it was made from: for argon with the form of Ilkov et
al., 1.11 times as much at 800 nm and 1e14 W/cm² and 1.85 times at 400 nm and 2e13
W/cm² (TC-662).

Here the factor is the rate's own average over a cycle,
$\frac{2}{\pi}\int_0^{\pi/2} W(F\cos\varphi)\,d\varphi / W(F)$, computed by
quadrature; it tends to $(3F/\pi F_0)^{1/2}$ in the tunnelling limit. (An expansion in
the rate's slope and curvature, {func}`~optmat.gas.ionisation_rates.cycle_average`,
gives the same to 1 % up to the barrier-suppression field, but not above it, where the
static correction bends the rate: for argon it was 24 % low at twice that field, made
the rate turn over just above it, and left rounding noise of 3e-9; TC-828.) The model
is then:

1. in a **static field**, the tunnelling rate with its correction, $N_\text{eff}
   \sum_c p_c\, w_{\text{stat},c} S_c$: exact in the limit the references test best;
2. the **cycle-averaged rate** is the average of that static rate over a cycle, times
   the enhancement of the theory's rate by the frequency, $E = Q\,w^\text{sr} /
   [w_\text{stat}\,(3F/\pi F_0)^{1/2}]$, which is 1 for $\gamma \to 0$;
3. the **instantaneous rate** is the cycle-averaged rate over the average over a cycle
   of the rate with the enhancement.

Above the maximum field, twice the barrier-suppression field where Tong and Lin's
correction ends, the rate is held at its value there, with a warning: the correction
extrapolated would make it turn over at 2.7 to 3.1 times that field (TC-816,
TC-828).

What this does not give is the **timing** within the cycle at $\gamma \gtrsim 1$,
where ionisation is less confined to the crests of the field than any function of the
instantaneous field implies (Yudin & Ivanov 2001). A plasma current that depends on
that timing inherits the error; the density does not.

The models of the literature's forms keep the factor of the tunnelling limit
(`averaging = "tunnelling"`), so that they are what the codes compute.

The rate is a function of one wavelength. For a field of many colours (a
supercontinuum, a pump with a dispersive wave) it is evaluated at the wavelength
given, usually the pump's. In the vacuum ultraviolet it fails: for argon at 140 nm a
rate of this kind is an order of magnitude too high, and the models warn below 350 nm
(TC-653).

```{figure} /auto_gallery/ionisation/images/sphx_glr_plot_ionisation_rates_001.*
:target: ../auto_gallery/ionisation/plot_ionisation_rates.html
:alt: The rates of the noble gases and the consistency of the two rates

The rates of the five gases at 800 nm, argon's against the wavelength, and the average
of the instantaneous rate over a cycle against the cycle-averaged rate.
```

## The energy per event

A consumer that removes the ionisation energy from the field needs it consistent with
the rate (OEP 1, section 4.6). Each channel $c$ of the ion has its energy $I_c$ and
rate $\Gamma_c$, the channel's weight included. `channel_rates(F)` gives them: the
channel terms of the rate, each multiplied by the one averaging factor the model
applies to their sum, so they sum to `instantaneous_rate(F)` exactly.
`energy_loss_rate(F)` is $L = \sum_c I_c \Gamma_c$, in W per atom, zero at zero field:
the quantity to average over a cycle beside the rate, as $\langle L\rangle / \langle
\Gamma\rangle$. `energy_per_event(F)` is $\bar I = L/\Gamma$, computed from the
logarithms of the channel terms so that it is defined where the rates underflow, and
equal to the lowest channel's energy at zero field. For argon at 800 nm it moves from
15.797 eV at 5 × 10⁹ V/m to 15.814 eV at 5 × 10¹⁰ V/m, between channels of 15.760 and
15.937 eV. The molecular models have one channel, and their $\bar I$ is the energy the
exponent uses; `energy_kind` says which (`vertical` for every molecular default,
`adiabatic` for the atoms, where the two coincide).

The rates are for a linearly polarised field. The instantaneous rate along an
elliptical field trajectory is a quasi-static approximation, which optmat does not
validate.

```python
from optmat.constants import E_CHARGE
from optmat.gas.ionisation_rates import ionisation_rate_model

argon = ionisation_rate_model("Ar", 800e-9)
fields = [5e9, 5e10]  # V/m
rates = argon.channel_rates(fields)  # 1/s, one row per channel
assert abs(rates.sum(axis=0) / argon.instantaneous_rate(fields) - 1).max() < 1e-12
mean = argon.energy_per_event(fields) / E_CHARGE  # eV
assert 15.76 < mean[0] < mean[1] < 15.94
```

## What is left out, and why

**The factor $(1 + \gamma^2)^{-n^*}$.** Ilkov, Decker & Chin (1992, eq. 8) and Bergé
et al. (2007, eq. B14), whose form most propagation codes use, carry the Coulomb
factor as $(2F_0/F\sqrt{1 + \gamma^2})^{2n^*}$, and cite Perelomov & Popov (1967) for
it. It is not there. Their eq. (42) is $(2F_0/F)^{2n^*}$, of which they write that
"the dependence ... on the frequency has dropped out, i.e., the Coulomb correction in
the ionization probability $w(F, \omega)$ has the same form as for a static field";
Perelomov, Popov & Kuznetsov (1968) repeat that it "does not depend on $\gamma$ at
all". The other form appears when the power of the short-range rate, which does carry
$\sqrt{1 + \gamma^2}$, and that of the Coulomb factor are written under one bracket
(TC-670). Against the reference calculations of
{doc}`ionisation-laser` and {doc}`ionisation-one-electron` it is too low by 1.4 at $\gamma = 0.74$, 1.8 at $\gamma = 1$ and 7 at $\gamma
= 2.8$. Popruzhenko et al. (2008) derive $(1 + 2\gamma/e)^{-2n^*}$, which is lower still
for $\gamma < 3.3$. Both are available (`coulomb="ilkov"`, `"popruzhenko"`, and the
models `ppt_ilkov` and `ppt_popruzhenko`), to reproduce rates computed with them.

**The range of the static factor.** Perelomov, Popov & Kuznetsov (1968) derive it for
a Coulomb field that does not distort the path of the electron under the barrier:
$\gamma \ll \gamma_k = 2I_p/n^*\omega$ (their eq. 20), and more exactly $s =
1/\gamma_k$ below a critical $s_\text{cr}(\gamma)$, which their Table I gives from
$\gamma = 3$ (0.040) to 20 (0.025). At 800 nm $s$ is 0.023 for He, 0.029 for Ne, 0.046
for Ar, 0.055 for Kr and 0.068 for Xe, and twice that at 400 nm. At 800 nm helium
and neon are inside the range at every $\gamma$ of use; argon, krypton and xenon
leave it between $\gamma = 1$ and 3, and at 400 nm so does helium. That the factor
nevertheless reproduces the one-electron calculations of argon to $\gamma = 2.8$
was a finding of the earlier comparison, not of the theory (TC-670). The later
pulse calculations motivated the empirical frequency-dependent factor
({doc}`ionisation-one-electron`; TC-848).

**A Stark shift of the ionisation energy**, $I_p \to I_p + (\alpha_0 - \alpha_+)F^2/2$.
To first order in the field it multiplies the rate by $\exp(-\kappa\,\Delta\alpha\,F)$,
which is one term of the first-order correction of the weak-field asymptotic theory
(Trinh et al. 2015, eq. 3b). For helium it is $-1.48$ of a total coefficient $B =
-2.60$. The rest are of the same size. The exact rates to which $S$ was fitted contain
all of them for the active electron, and the all-electron rates on which $N_\text{eff}$
rests contain them for the atom. Adding the shift to either counts it twice (TC-661).

**A factor for the polarisation of the ion**, $\exp(-2\alpha_+ F)$, has no identified source in the papers reviewed (TC-661). Its nearest relative is Zhao & Yuan's (2014) eq. (9),
$\exp[-(2\alpha_+ F/\kappa)(1/z_0 - 1/z_1)]$ with $z_0$ the radius of the core and
$z_1$ the tunnel exit, derived for calcium. For argon at 0.1 au that is about 0.6; the
factor above is 0.26. Majety & Scrinzi (2015) find that the first excited channel of
the ion changes the static rates of neon and argon "by less than 2%" (TC-661).

**The electrons with $m = \pm 1$.** Their rate is smaller by $F/2F_0$ in the
tunnelling limit: 2 to 4 % of the $m = 0$ yield below 2e14 W/cm² in the one-electron
calculations of argon. They are inside $N_\text{eff}$.

**The second electron.** The models give the rate of the first ionisation. Beyond
saturation the ion ionises in turn (27.6 eV for Ar⁺).
