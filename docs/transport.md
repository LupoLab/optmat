# Transport of the dilute gases

A gas-filled fibre's fill needs, besides the index, three material properties of each gas
(OEP 1, section 4.3):

- the **dilute-gas shear viscosity** $\eta(T)$, Pa s: the zero-density limit, a function
  of temperature alone, as kinetic theory and Wilke's rule assume;
- the **ideal-gas molar isobaric heat capacity** $c_p^\circ(T)$, J mol⁻¹ K⁻¹: the
  temperature derivative of the ideal-gas part of the gas's equation of state, not a
  value at a finite pressure;
- the **molar mass** $M$, kg/mol, the one the equation of state uses.

`optmat.gas.transport.transport_model(gas)` gives all three for the thirteen gases
listed below, with their temperature ranges, sources, status and identity ({doc}`identity`).
The caller combines them for a flow problem, sound speed or mean free path.
Read {doc}`conventions` for state and density definitions; the {doc}`quickstart`
introduces the material interfaces, and {doc}`mixtures` explains mixture densities.

```python
import numpy as np
from optmat.gas.transport import (
    heat_capacity,
    mixture_viscosity,
    molar_mass,
    transport_model,
    viscosity,
)

argon = transport_model("Ar")
assert argon.identity.sources == ("lemmon_2004", "tegeler_1999")
eta = viscosity("Ar", np.array([250.0, 294.0, 350.0]))  # Pa s
assert 2.2e-5 < eta[1] < 2.3e-5
assert abs(float(heat_capacity("Xe", 294.0)) / 20.786156 - 1) < 1e-7  # (5/2) R
assert molar_mass("N2") == 0.02801348
low, high = argon.temperature_range("viscosity")  # K
```

## The sources

| Gas | Viscosity (form) | K | Stated uncertainty | $c_p^\circ$, $M$ | K | Status |
|---|---|---|---|---|---|---|
| He | `arp_1998` (CoolProp) | 110-1500 | at most 10 % above 100 K, for the correlation as a whole | `harvey_2023` | 2.1768-2000 | validated |
| Ne | `huber_2018` (Huber's dilute form) | 24.56-725 | 2 % (k = 2) | Thol et al. 2019, not in hand | 24.56-725 | validated |
| Ar | `lemmon_2004` (CoolProp) | 83.806-2000 | 0.5 % in the dilute gas | `tegeler_1999` | 83.806-2000 | validated |
| Kr | `huber_2018` (Huber's dilute form) | 115.77-750 | 3 % (k = 2) | `lemmon_2006` | 115.77-750 | validated |
| Xe | `huber_2018` (Huber's dilute form) | 161.4-750 | 0.7 % (k = 2) | `lemmon_2006` | 161.4-750 | validated |
| H2 | `muzny_2013` (CoolProp) | 13.957-1000 | 0.1 % from 200 K to 400 K up to 0.1 MPa | `leachman_2009` | 13.957-1000 | validated |
| D2 | `assael_1987` (correlation; Table 2 the check) | 30-2200 | 1 % from 250 K to 350 K, 2 % below, at most 4 % at the highest temperatures | `richardson_2014` | 18.724-600 | validated |
| N2 | `lemmon_2004` (CoolProp) | 63.151-2000 | 0.5 % in the dilute gas | `span_2000` | 63.151-2000 | validated |
| O2 | `lemmon_2004` (CoolProp) | 54.361-2000 | 1 % in the dilute gas above 200 K, 5 % below | `schmidt_1985`, `stewart_1991` | 54.361-1700 | validated |
| CH4 | `quinones_cisneros_2006` (CoolProp) | 95-600 | not stated (the fit's average deviation, 0.31 %) | `setzmann_1991` | 90.6941-625 | provisional |
| SF6 | `quinones_cisneros_2012` (CoolProp) | 223.555-1000 | 0.3 % (k = 2) from 298 K to 700 K, 1.5 % at the triple point and 1000 K | `guder_2009` | 223.555-625 | validated |
| N2O | `huber_2018` (Huber's dilute form) | 182.33-525 | 5 % (k = 2) | `lemmon_2006` | 182.33-525 | validated |
| CO2 | `laesecke_2017` (CoolProp) | 200-2000 | 0.2 % from 200 K to 700 K, 1 % to 2000 K | `span_1996` | 216.592-2000 | validated |

The choices are curated in `data/gas/transport.toml`; each record says what was checked
against its paper.

**CoolProp's correlations.** For eight gases the viscosity is CoolProp's correlation at
10⁻⁶ mol/m³, where its density terms are below 10⁻¹⁰. Each reproduces the dilute-gas form
printed in its paper, with the paper's constants (recorded as formulas), to below 10⁻⁹.
The ranges are the papers': Lemmon and Jacobsen's equations hold over their equations
of state's states, Arp et al. use Maitland and Smith's dilute values above 110 K, and
Laesecke and Muzny state an uncertainty from 200 K, and Quiñones-Cisneros et al. extend
SF6's to 1000 K (TC-765). Methane's has no stated uncertainty (TC-763).

**Ne, Kr, Xe and N2O.** CoolProp has no viscosity for these. optmat evaluates Huber's
dilute-gas form (NISTIR 8209: Chapman-Enskog with Neufeld et al.'s collision integral,
his Lennard-Jones parameters and his factor $F_c$), which reproduces his check values to
2×10⁻⁷. These are the models behind the NIST WebBook's tables, kept in
`data/external/nist_webbook` as the check. At 0.1 MPa the tables lie above the dilute
values by the density term: 0.013 %, 0.13 %, 0.21 % and 0.21 % at 298.15 K, and at most
0.54 % (Xe at 200 K). For Ne, Kr and Xe, $F_c$ was adjusted to Assael et al.'s (2018)
reference values at 25 °C, so those are not an independent check.

**D2.** Assael, Mixafendi and Wakeham's (1987) correlation is evaluated, and their
Table 2 is its check. The paper's printed Eq. (4) does not reproduce the table: it is
the authors' universal corresponding-states functional (1986, Eq. 13), not their
hydrogen correlation (1986, Eq. 3). With the latter, the D2 mass and the paper's
quantum factor $1 - e^{-f T^{*1/3}}$ (Eq. 7) at every temperature, the table is
reproduced within 1.5 units of its last digit, 5.1×10⁻⁴, except for the 90 K cell,
degraded in the scan, 1.2×10⁻³ off (TC-760, TC-765). Until 0.2.1 the table itself was
served as a cubic spline in $\log\eta$ against $\log T$; its knots' four-digit
rounding, 2×10⁻⁴, made the spline wiggle at 10⁻⁴, so that no smooth representation of
it (a Chebyshev series over 60-600 K) came closer than 4×10⁻⁴. The correlation is
analytic: degree 20 in $\ln T$ reaches 10⁻¹¹ (TC-839).

**The heat capacity** is the equation's dimensionless $c_p^\circ/R$ times CODATA's $R$, the
convention of the number density ({doc}`conventions`). The noble gases give
$\tfrac52 R$ exactly, whatever gas constant their equation was published with (helium's
in CoolProp is 1.5 ppm below the one Harvey et al. give). O2's ideal part in CoolProp is
a refit of Wagner, Ewers and Schmidt's (Schmidt and Wagner 1985; Stewart et al. 1991,
Eq. 11), within 3.3 ppm to 1500 K. Above 1700 K it departs by more than that
equation's stated 0.003 J mol⁻¹ K⁻¹, so O2's heat capacity stops there (TC-765).

## Mixtures: Wilke's rule

`mixture_viscosity(mixture, T)` is Wilke's rule (`wilke_1950`, Eqs. 13-14):

$$
\eta_\text{mix} = \sum_i \frac{x_i\,\eta_i}{\sum_j x_j\,\Phi_{ij}}, \qquad
\Phi_{ij} = \frac{\left[1 + (\eta_i/\eta_j)^{1/2}(M_j/M_i)^{1/4}\right]^2}{\left[8\,(1 + M_i/M_j)\right]^{1/2}}.
$$

It is a kinetic-theory approximation: Wilke reports an average deviation of 0.97 % for
his binary mixtures (0.49 % without H2-Ar, the worst) and 1.9 % for multicomponent ones.
A mixture is a named composition (`"dry air"` is `"dry air (N2, O2, Ar)"`, as for the air
model) or a mapping of gases to mole fractions summing to one. CIPM-2007's full
composition is refused: CO has no transport model. The ideal mixture's heat capacity is
$\sum_i x_i c_{p,i}^\circ$.

```python
air = float(mixture_viscosity("dry air", 294.0))
assert 1.80e-5 < air < 1.83e-5
he_ar = mixture_viscosity({"He": 0.5, "Ar": 0.5}, np.array([250.0, 300.0]))
```

(against-solprop)=
## Finite-pressure values

The dilute-gas limit should be distinguished from a viscosity evaluated at finite
pressure. For example, the nitrogen correlation at 294 K gives 1.7612×10⁻⁵ Pa s at
1 bar and 1.7599×10⁻⁵ Pa s in the dilute limit (TC-762). The difference comes from
the density term. Likewise, a finite-pressure heat-capacity ratio need not equal
$c_p^\circ/(c_p^\circ-R)$.

The {doc}`/developers/testing` guide describes regression checks retained for
historical calculations. For current model selection, use the sources and ranges
above and the {mod}`optmat.gas.transport` reference.

This page documents what [OEP 1](oeps/archive/oep-0001-propagation-interface.md)
proposed and implemented: its theory, decisions and phase reports are the record.
