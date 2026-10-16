# Third-order extrapolation

Read {doc}`nonlinear-support` for the meanings of fitted legs, frequency sums and
resonance limits. These numerical comparisons test possible extrapolation rules;
they do not supply a universal error bound. For supported calculations and their
uncertainties, return to {doc}`/nonlinear` and {doc}`/nonlinear_accuracy`.

## Beyond the fitted data: bounds on the legs

A propagation code that evaluates every four-wave quartet must decide what to do with
the legs beyond the data, potentially acting on each leg separately. That decision is the code's (TC-868); this section gives evidence for it. It
evaluates caps $\nu_j \to \min(|\nu_j|, \nu_c)$ at four
places (the manual's *Where the fits are supported* gives the legs and the resonances;
`tools/fit/nonlinear_support.py` computes everything here, TC-866 and TC-867):

- **(a)** the Kerr diagonal of the largest fitted $W_2$, $\sqrt{W_{2,\max}/4}$;
- **(b)** the largest fitted leg;
- **(b')** $\sqrt{3W_{2,\max}/4}$, the largest leg any quartet within the fitted $W_2$
  can have;
- **(c)** the sufficient resonance bound $\min(E_1, E_2/2)$ with $E_2$ of any
  polarisation (not defined where $E_1$ or $E_2$ is not a curated level).

### On the data and within the fitted W2

Per cap: its wavelength (nm); the fitted values it changes, of all, and the largest
change; then the worst change of a quartet within the fitted $W_2$ (20 000 random
frequency-conserving quartets with $W_2 \le W_{2,\max}$, and the third harmonic at
$W_{2,\max}$, whose output leg is the extreme):

<!-- nonlinear-caps:start -->
| Gas | (a) Kerr diagonal | (b) largest fitted leg | (b') √(3 W2,max / 4) | (c) min(E1, E2/2) |
|---|---|---|---|---|
| He | 263.1: 7/52, -13.3 %; -13.3 % | 151.9: 0/52, exact; 0 | 151.9: 0/52, exact; 0 | 120.3: 0/52, exact; 0 |
| Ne | 373.9: 11/26, -6.9 %; -8.2 % | 229.0: 0/26, exact; -1.4 % | 215.9: 0/26, exact; 0 | 134.9: 0/26, exact; 0 |
| Ar | 374.0: 27/41, -15.2 %; -17.9 % | 229.0: 0/41, exact; -3.1 % | 215.9: 0/41, exact; 0 | 192.1: 0/41, exact; 0 |
| Kr | 398.5: 4/6, -18.6 %; -21.7 % | 244.1: 0/6, exact; -3.9 % | 230.1: 0/6, exact; 0 | 219.4: 0/6, exact; 0 |
| Xe | 398.5: 4/10, -25.6 %; -29.7 % | 244.1: 0/10, exact; -5.4 % | 230.1: 0/10, exact; 0 | 258.8: 2/10, -5.1 %; -10.2 % |
| H2 | 374.0: 27/56, -16.6 %; -19.5 % | 229.0: 0/56, exact; -3.4 % | 215.9: 0/56, exact; 0 | 201.7: 0/56, exact; 0 |
| D2 | 389.2: 18/46, -16.3 %; -17.7 % | 231.5: 0/46, exact; -1.5 % | 224.7: 0/46, exact; 0 | 201.1: 0/46, exact; 0 |
| N2 | 373.9: 27/37, -14.2 %; -16.7 % | 229.0: 0/37, exact; -2.9 % | 215.9: 0/37, exact; 0 | 290.1: 9/37, -8.8 %; -11.3 % |
| O2 | 389.2: 18/20, -20.7 %; -23.9 % | 238.3: 0/20, exact; -4.5 % | 224.7: 0/20, exact; 0 | 304.6: 7/20, -13.6 %; -17.3 % |
| CO2 | 398.5: 7/11, -18.0 %; -18.4 % | 231.5: 0/11, exact; -0.2 % | 230.1: 0/11, exact; 0 | n/a |
| CH4 | 389.2: 12/13, -21.6 %; -25.1 % | 238.3: 0/13, exact; -4.6 % | 224.7: 0/13, exact; 0 | n/a |
| SF6 | 374.0: 6/8, -13.3 %; -15.3 % | 229.0: 0/8, exact; -2.9 % | 215.9: 0/8, exact; 0 | n/a |
| CO | 401.0: 3/6, -16.9 %; -16.8 % | 231.5: 0/6, exact; 0 | 231.5: 0/6, exact; 0 | 309.9: 1/6, -11.6 %; -11.4 % |
| N2O | 664.1: 1/1, -2.1 %; -6.3 % | 580.0: 0/1, exact; -5.3 % | 383.4: 0/1, exact; 0 | n/a |
| O3 | 694.6: 0/4, exact; -4.5 % | 694.6: 0/4, exact; -4.4 % | 401.0: 0/4, exact; 0 | n/a |
<!-- nonlinear-caps:end -->

Only (b') is exact on the data and on every quartet within the fitted $W_2$, by
construction. (b) is exact on the data and changes in-support quartets by up to 5 %
(Xe). (a) changes up to 27 fitted values of a gas, by up to 26 %. (c) is not a leg cap
at all for the pair limit: a pair resonance is reached by two legs together, and capping
single legs at $E_2/2$ cuts the data of Xe, N2, O2 and CO, whose harmonic legs pair only
with small ones.

### Against a carrier

A constant Kerr gamma evaluated at the carrier is a useful reference for testing
a per-leg cap. This comparison uses 20 000 quartets per gas, with three
inputs in the declared range (each absorbed or emitted) and the output wherever it
falls, as four-wave mixing makes it. Where the quartet's $W_2$ is within the fitted data
the fit is the reference, and the table counts the quartets where the capped value is
further from it than the carrier's. Beyond the fitted data no truth is known; the capped
value is counted as *shown no worse* when the carrier is at or below it, on the
assumption that the truth lies at or above the polynomial there. The calculations below
support that assumption for the processes and gases they cover; it is not proven for
arbitrary mixed quartets (TC-869). The rest are counted as *not shown*. Quartets with a
leg beyond $E_1$ or a pair beyond $E_2$ are past a resonance, where the polynomial says
nothing; for CH4, SF6, N2O and O3 the column counts quartets beyond their absorption
evidence instead. The
carrier is at 400 nm, the hardest of the three (the last column: (b') with 800 and
1030 nm carriers):

<!-- nonlinear-carrier:start -->
| Gas | (a): worse in support / not shown beyond | (b) | (b') | (b'): shown no worse beyond | Past a resonance | (b'), 800 and 1030 nm carriers: worse or not shown |
|---|---|---|---|---|---|---|
| He | 3072 / 0 | 0 / 0 | 0 / 0 | 90 | 0 | 0 |
| Ne | 2789 / 2683 | 0 / 0 | 0 / 0 | 2670 | 0 | 0 |
| Ar | 2677 / 2688 | 0 / 0 | 0 / 0 | 2673 | 0 | 0 |
| Kr | 1949 / 3576 | 188 / 883 | 0 / 0 | 3565 | 0 | 0 |
| Xe | 1830 / 3537 | 145 / 1091 | 0 / 26 | 3502 | 37 | 0 |
| H2 | 2694 / 2688 | 0 / 0 | 0 / 0 | 2673 | 0 | 0 |
| D2 | 2251 / 3236 | 12 / 0 | 0 / 0 | 3241 | 0 | 0 |
| N2 | 2673 / 2683 | 0 / 0 | 0 / 0 | 2670 | 0 | 0 |
| O2 | 2172 / 1532 | 220 / 184 | 0 / 0 | 1534 | 1707 | 0 |
| CO2 | 1939 / 3576 | 4 / 9 | 0 / 0 | 3565 | 0 | 0 |
| CH4 | 2209 / 1767 | 220 / 184 | 0 / 0 | 1778 | 1463 | 0 |
| SF6 | 2633 / 1219 | 0 / 0 | 0 / 0 | 1210 | 1463 | 0 |
| CO | 1973 / 3539 | 0 / 47 | 0 / 37 | 3508 | 117 | 0 |
| N2O | 0 / 10839 | 0 / 10837 | 0 / 10730 | 0 | 2471 | 0 |
| O3 | 0 / 0 | 0 / 0 | 0 / 0 | 0 | 3628 | 0 |
<!-- nonlinear-carrier:end -->

(a) is worse than the 400 nm carrier for 11-16 % of the quartets within the fitted
$W_2$ of every gas with ESHG or THG data: the argon quartet 400 + 400 - 800 → 266.7 nm, at 0.98 of the largest fitted $W_2$,
is the example (-9.3 % against -4.3 %). (b') is never worse within the data, and
beyond them is shown no worse except for a few dozen xenon and CO quartets that mix
legs above and below a 400 nm carrier. N2O's 400 nm carrier is itself beyond its single
fitted value, and none of its quartets beyond the data can be shown no worse than it;
against 800 and 1030 nm carriers all can, as for every gas.

### The polynomial beyond the data

A capped value lies below the physical gamma only if the physical gamma lies above the
polynomial beyond the data. The calculations that the windows leave out, and the validators beyond
the largest fitted $W_2$, test that (model, with the calculation's fitted scale, over
value, at the nearest and the farthest point):

<!-- nonlinear-beyond:start -->
| Gas | Dataset | Role | Process, fundamentals | W2 / W2,max | Model / value, nearest → farthest |
|---|---|---|---|---|---|
| He | `bishop_1989/omega_kerr_zzzz` | windowed | Kerr, 253.1-227.8 nm | 1.08-1.33 | 0.997 → 0.992 |
| He | `bishop_1989/omega_eshg_zzzz` | windowed | ESHG, 303.8-227.8 nm | 1.12-2.00 | 0.998 → 0.969 |
| He | `bishop_1989/omega_thg_zzzz` | windowed | THG, 414.2-227.8 nm | 1.21-4.00 | 0.997 → 0.720 |
| Ne | `hattig_1998/ccsd_thg_Ne` | windowed | THG, 632.8-457.9 nm | 1.05-2.00 | 0.999 → 0.993 |
| Ar | `hattig_1998/ccsd_thg_Ar` | windowed | THG, 632.8-488.0 nm | 1.05-1.76 | 1.001 → 0.967 |
| H2 | `bishop_1991/thg_par_elec_H2` | windowed | THG, 632.8-488.0 nm | 1.05-1.76 | 0.997 → 0.946 |
| D2 | `bishop_1991/thg_par_elec_D2` | windowed | THG, 632.8-488.0 nm | 1.13-1.91 | 0.998 → 0.951 |
| N2 | `hattig_1998/ccsd_thg_N2` | windowed | THG, 632.8-457.9 nm | 1.05-2.00 | 1.002 → 0.964 |
| CH4 | `hattig_1998/ccsd_thg_CH4` | validator | THG, 671.5-476.5 nm | 1.01-2.00 | 1.195 → 1.179 |
<!-- nonlinear-beyond:end -->

Farther beyond the data every calculation lies above the polynomial, increasingly with
$W_2$: at 1.8 to 2 times the fitted $W_2$ by 0.7 % (Ne) to 3.6 % (N2), by 5 % for H2 and
D2, and by 28 % for helium's third harmonic at four times. Not every point does: just
beyond the window argon's THG sits at 1.001 and N2's at 1.002, within the fit's scatter.
These comparisons include the fitted scales and cover particular processes; they give no
universal lower bound and no error bar on extrapolation (TC-869). CH4's validator is off in level (Hättig &
Jørgensen's CCSD was rejected as CH4's shape, TC-314), but falls relative to the model
in the same way. No calculation reaches beyond twice the fitted $W_2$ except helium's.

### Representative four-wave-mixing quartets

Deep-ultraviolet four-wave mixing in argon and vacuum-ultraviolet mixing in helium: the polynomial, each cap and the carriers,
in au:

<!-- nonlinear-cases:start -->
| Quartet | W2 / W2,max | Largest leg / fitted | Largest leg / E1 | Largest pair / E2 | Polynomial (au) | (a) | (b) | (b') | Kerr carrier 400 / 800 / 1030 nm | Past a resonance |
|---|---|---|---|---|---|---|---|---|---|---|
| Ar 400 + 400 - 800 -> 266.7 nm | 0.98 | 0.86 | 0.40 | 0.473 | 1691.7 | 1534.7 | 1691.7 | 1691.7 | 1619.8 / 1257.7 / 1220.4 | no |
| Ar SPM at 266.7 nm | 1.97 | 0.86 | 0.40 | 0.710 | 2537.6 | 1705.4 | 2537.6 | 2537.6 | 1619.8 / 1257.7 / 1220.4 | no |
| He 250 + 250 - 1500 -> 136.4 nm | 1.49 | 1.11 | 0.43 | 0.481 | 65.9 | 53.2 | 62.6 | 62.6 | 48.5 / 44.4 / 43.9 | no |
| He 200 + 200 - 600 -> 120 nm | 2.11 | 1.27 | 0.49 | 0.601 | 79.2 | 53.9 | 69.5 | 69.5 | 48.5 / 44.4 / 43.9 | no |
| He SPM at 120 nm | 4.81 | 1.27 | 0.49 | 1.002 | 163.3 | 57.2 | 102.6 | 102.6 | 48.5 / 44.4 / 43.9 | yes |
<!-- nonlinear-cases:end -->

The argon quartet is within the fitted $W_2$ and is given exactly by (b) and (b'). Its
SPM at 266.7 nm is at twice the fitted $W_2$ but its legs are within the fitted legs:
(b') leaves it uncapped; the table above suggests the polynomial underestimates
there, without proving it. The two helium four-wave
quartets are beyond the fitted legs and $W_2$, below every resonance: (b') caps their
output legs. Helium's SPM at 120 nm is past the 1s2s ¹S two-photon resonance (its pair
is 1.002 times $E_2$): the polynomial's 163 au means nothing there.
