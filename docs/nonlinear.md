# Third-order nonlinearity

The Kerr terms (self- and cross-phase modulation) and third-harmonic generation of a
propagation need the gas's third-order response at the frequencies involved. optmat
provides it as the electronic second hyperpolarizability gamma of each gas, a function
of the process and of all four frequencies. It is fitted to the curated third-order
data (TC-224 to TC-314) and stored with its covariance.

This page introduces the response and its use. Continue with
{doc}`/user-guide/nonlinear-conventions`, {doc}`/user-guide/nonlinear-models`,
{doc}`/user-guide/nonlinear-support` and {doc}`nonlinear_accuracy` for the
conversions, theory, data support and uncertainty. The {doc}`quickstart` introduces
the material interfaces; {doc}`raman_response` describes the delayed response.

## What is modelled

- **The electronic second hyperpolarizability** $\gamma(-\omega_\sigma; \omega_1,
  \omega_2, \omega_3)$ of an isolated molecule, orientationally averaged (the mean gamma,
  the one an isotropic gas shows for parallel linear polarizations).
- **Below resonance:** from the static limit to a fundamental of 400 nm (700 nm for
  O3), and to 20 µm.
- **For every third-order process at once:** the Kerr effect, third-harmonic generation,
  ESHG, the dc Kerr effect, CARS-type four-wave mixing and a nondegenerate pump-probe
  index.

What is not modelled here:
- **The vibrational and rotational (Raman) response:** in propagation it is a separate,
  delayed term. The one exception is the small vibrational part of the ESHG of CO2, CH4
  and SF6, which could not be separated; see {doc}`/user-guide/nonlinear-conventions`.
- **Ionisation and plasma.**
- **Pair (collision-induced) contributions:** gamma is per molecule, and the
  third-order susceptibility is linear in density.
- **Higher-order Kerr terms:** the claims of such terms were refuted (TC-313).

## Detailed chapters

<span id="conventions"></span>
<span id="the-taylor-series-gamma"></span>
<span id="units"></span>
<span id="susceptibility-and-n2"></span>
<span id="the-sources-conventions-and-their-conversion"></span>
<span id="processes-and-their-frequencies"></span>
<span id="parts-electronic-vibrational-rotational"></span>

See {doc}`/user-guide/nonlinear-conventions`.

<span id="the-dispersion-model"></span>
<span id="the-time-domain-form"></span>
<span id="the-fits"></span>
<span id="the-fit"></span>
<span id="the-choices-gas-by-gas"></span>
<span id="values"></span>

See {doc}`/user-guide/nonlinear-models`.

<span id="configuration"></span>
<span id="running-and-storing"></span>

For configuration and regeneration procedures, see {doc}`/developers/fitting`.

<span id="where-the-fits-are-supported"></span>
<span id="the-legs-of-the-data"></span>
<span id="resonance-limits"></span>
<span id="does-gamma-increase-with-each-frequency"></span>
<span id="quantities-derived-from-the-facts"></span>
<span id="the-support-record"></span>

See {doc}`/user-guide/nonlinear-support`.

## Using the models

```python
from optmat.conditions import State
from optmat.gas.gases import Gas
from optmat.gas.nonlinear import hyperpolarizability_models

argon = hyperpolarizability_models()[Gas.Ar]
argon.gamma_au("kerr", 800e-9)  # 1258 au: SPM and XPM at 800 nm
argon.gamma_au("thg", 800e-9)  # 1480 au: the third harmonic of 800 nm
argon.gamma_au("pump_probe", 800e-9, 400e-9)  # the index a 400 nm probe sees
argon.gamma_sigma_au("kerr", 800e-9)  # statistical one sigma, au
argon.gamma3("kerr", 800e-9)  # gamma_T / (6 epsilon_0), m^5/V^2
argon.n2(800e-9, State(1e5, 293.15))  # self n2 at 1 bar, 20 C, m^2/W
argon.parameters  # gamma0, A, B, B'
argon.record["validation"]  # the validators against the fit
```

`gamma_au` takes arrays of wavelengths, `inf` for the static limit, and for `cars` and
`pump_probe` the second wavelength. For refitting commands and a Python fitting example, see {doc}`/developers/fitting`.

## Limitations

- **Resonance:** the form holds below resonance, within the fitted $W_2$. It says nothing
  near a resonance (O3's third harmonic, the vacuum-ultraviolet harmonics of short
  fundamentals, two-photon resonances), nor past one: a leg beyond $E_1$ or a pair sum
  beyond $E_2$ ({doc}`/user-guide/nonlinear-support`). Xenon's n2 at 248.4 nm, 1 nm above the
  two-photon resonance with 6p[1/2]_0 at 249.6 nm, is negative and large: -(2.06 +-
  0.14) and -(1.82 +- 0.08) x 10^-14 esu measured, -(1.54 +- 0.30) calculated, about
  -8e-21 m2/W at 1 atm (`lehmberg_1995`, TC-854); a fundamental near 250 nm in xenon is
  outside the model.
- **The level:** it rests on the ESHG scale, which Shelton himself states to about 1 %.
  The Maryland n2 lie up to 16 % lower (0.84-1.03 at 800 nm), an open question (accuracy
  page).
- **The unseparated vibrational parts** of CO2, CH4 and SF6.
- **N2O and O3** are known to 15 %; CO to 10 %.
- **Dilute limit:** no local-field factors, and no pair or density dependence.
- **Temperature:** the electronic gamma is taken as independent of temperature.

The support record is what
[OEP 5](oeps/archive/oep-0005-hyperpolarizability-support.md) proposed and implemented:
its theory, decisions and phase report are the record.
