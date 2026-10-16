# Noble-gas ionisation: one-electron calculations

One-electron calculations isolate the Coulomb correction from the effective electron count. This chapter retains the convergence checks, pulse definitions, hydrogen tests, and evidence for the empirical frequency-dependent factor. Read {doc}`ionisation-theory` first; compare the resulting atomic rate with {doc}`ionisation-measured`.

## One-electron calculations of argon: 800 nm (6 fs, and Wallace et al.'s 5.5 fs), 400, 800 and 1030 nm (10 fs), 400 and 1030 nm (30 fs)

Time-dependent Schrödinger calculations with tRecX for one electron in Tong & Lin's
potentials (`travers_2026_trecx`, tier 3 of TC-439), against the `one_electron` model
with each Coulomb factor. The yield is the ionised fraction at the end of the pulse,
the population of the bound states up to n = 10 subtracted, with its uncertainty from
convergence runs (TC-812). The first runs' measures, the population beyond 20 bohr and
the integral of the photoelectron spectrum, are kept in the record but not compared:
the first counts the atoms left in Rydberg states, up to 30 % of the yield, and the
second is low (TC-703). Rows within a factor of 10 of the numerical loss of a run
without a field are left out.

<!-- ionisation-rate-tdse:start -->
| Gas | nm | I (W/cm2) | gamma | TDSE: ionised fraction | ± (%) | `one_electron` (the default's factor) | static | `ppt_ilkov` | `ppt_popruzhenko` | the default's, on the carrier |
|---|---|---|---|---|---|---|---|---|---|---|
| Ar | 800 | 2.61e+13 | 2.25 | 2.130e-07 | 1.7 | 1.07 | 1.67 | 0.36 | 0.37 | 1.07 |
| Ar | 800 | 4.03e+13 | 1.81 | 5.222e-06 | 1.7 | 1.07 | 1.58 | 0.50 | 0.48 | 1.07 |
| Ar | 800 | 6.22e+13 | 1.46 | 9.954e-05 | 1.7 | 1.05 | 1.46 | 0.68 | 0.60 | 1.04 |
| Ar | 800 | 9.61e+13 | 1.17 | 1.406e-03 | 0.6 | 1.01 | 1.34 | 0.91 | 0.76 | 1.01 |
| Ar | 800 | 1.48e+14 | 0.94 | 1.441e-02 | 0.6 | 0.95 | 1.21 | 1.18 | 0.94 | 0.95 |
| Ar | 800 | 2.29e+14 | 0.76 | 1.010e-01 | 0.6 | 0.90 | 1.09 | 1.46 | 1.17 | 0.90 |
| Ar | 800 | 3.54e+14 | 0.61 | 3.766e-01 | 0.6 | 0.99 | 1.13 | 1.74 | 1.51 | 0.99 |
| Ar | 800 | 5.47e+14 | 0.49 | 7.872e-01 | 0.6 | 1.03 | 1.08 | 1.26 | 1.25 | 1.03 |
| Ar | 800 | 8.44e+14 | 0.40 | 9.699e-01 | 0.1 | 1.02 | 1.02 | 1.03 | 1.03 | 1.02 |
| Ar | 400 | 3.00e+13 | 4.19 | 1.655e-04 | 0.3 | 1.07 | 3.90 | 0.28 | 0.37 | 1.04 |
| Ar | 400 | 4.63e+13 | 3.37 | 8.591e-04 | 0.3 | 1.27 | 4.01 | 0.46 | 0.56 | 1.24 |
| Ar | 400 | 7.16e+13 | 2.71 | 6.302e-03 | 0.3 | 0.98 | 2.69 | 0.49 | 0.55 | 0.95 |
| Ar | 400 | 1.11e+14 | 2.18 | 5.120e-02 | 0.3 | 0.60 | 1.44 | 0.43 | 0.44 | 0.59 |
| Ar | 400 | 1.71e+14 | 1.76 | 1.313e-01 | 0.2 | 0.99 | 1.99 | 1.00 | 0.95 | 0.97 |
| Ar | 400 | 2.64e+14 | 1.41 | 3.175e-01 | 0.2 | 1.29 | 2.03 | 1.69 | 1.55 | 1.27 |
| Ar | 400 | 4.07e+14 | 1.14 | 7.318e-01 | 0.2 | 1.11 | 1.30 | 1.33 | 1.30 | 1.10 |
| Ar | 400 | 6.29e+14 | 0.92 | 9.093e-01 | 0.2 | 1.09 | 1.10 | 1.10 | 1.10 | 1.08 |
| Ar | 400 | 9.71e+14 | 0.74 | 9.436e-01 | 0.2 | 1.06 | 1.06 | 1.06 | 1.06 | 1.06 |
| Ar | 1030 | 2.39e+13 | 1.83 | 3.309e-08 | 4.2 | 0.83 | 1.13 | 0.33 | 0.31 | 0.83 |
| Ar | 1030 | 3.68e+13 | 1.47 | 1.216e-06 | 3.3 | 1.05 | 1.36 | 0.57 | 0.50 | 1.04 |
| Ar | 1030 | 5.69e+13 | 1.18 | 3.908e-05 | 3.3 | 1.01 | 1.26 | 0.75 | 0.62 | 1.01 |
| Ar | 1030 | 8.79e+13 | 0.95 | 8.390e-04 | 0.9 | 0.96 | 1.16 | 0.95 | 0.76 | 0.96 |
| Ar | 1030 | 1.36e+14 | 0.77 | 1.059e-02 | 0.9 | 0.99 | 1.17 | 1.31 | 1.03 | 0.99 |
| Ar | 1030 | 2.10e+14 | 0.62 | 9.535e-02 | 0.9 | 0.91 | 1.04 | 1.53 | 1.23 | 0.91 |
| Ar | 1030 | 3.24e+14 | 0.50 | 4.016e-01 | 0.9 | 1.02 | 1.11 | 1.74 | 1.54 | 1.02 |
| Ar | 1030 | 5.00e+14 | 0.40 | 8.726e-01 | 0.2 | 1.00 | 1.03 | 1.14 | 1.14 | 1.00 |
| Ar | 800 | 2.33e+13 | 2.38 | 1.456e-07 | 2.3 | 1.11 | 1.76 | 0.34 | 0.35 | 1.10 |
| Ar | 800 | 3.63e+13 | 1.91 | 3.806e-06 | 2.3 | 1.17 | 1.74 | 0.50 | 0.49 | 1.16 |
| Ar | 800 | 5.65e+13 | 1.53 | 7.785e-05 | 2.3 | 1.20 | 1.69 | 0.73 | 0.65 | 1.19 |
| Ar | 800 | 8.78e+13 | 1.23 | 1.379e-03 | 2.3 | 1.03 | 1.38 | 0.87 | 0.73 | 1.02 |
| Ar | 800 | 1.37e+14 | 0.98 | 1.468e-02 | 0.8 | 1.04 | 1.33 | 1.21 | 0.98 | 1.03 |
| Ar | 800 | 2.13e+14 | 0.79 | 1.126e-01 | 0.8 | 0.97 | 1.19 | 1.51 | 1.21 | 0.97 |
| Ar | 800 | 3.31e+14 | 0.63 | 4.524e-01 | 0.8 | 1.02 | 1.16 | 1.64 | 1.44 | 1.02 |
| Ar | 800 | 5.14e+14 | 0.51 | 8.927e-01 | 0.8 | 1.01 | 1.05 | 1.12 | 1.12 | 1.01 |
| Ar | 800 | 8.00e+14 | 0.41 | 9.704e-01 | 0.1 | 1.03 | 1.03 | 1.03 | 1.03 | 1.03 |
| Ar | 800 | 5.00e+13 | 1.62 | 2.644e-05 | 1.5 | 0.86 | 1.24 | 0.48 | 0.44 | 0.86 |
| Ar | 800 | 6.27e+13 | 1.45 | 9.829e-05 | 1.5 | 1.02 | 1.42 | 0.67 | 0.59 | 1.02 |
| Ar | 800 | 7.86e+13 | 1.30 | 4.097e-04 | 1.5 | 0.99 | 1.34 | 0.77 | 0.66 | 0.99 |
| Ar | 800 | 9.85e+13 | 1.16 | 1.503e-03 | 1.5 | 0.99 | 1.32 | 0.92 | 0.76 | 0.99 |
| Ar | 800 | 1.23e+14 | 1.03 | 5.294e-03 | 1.5 | 0.95 | 1.23 | 1.03 | 0.84 | 0.95 |
| Ar | 800 | 1.55e+14 | 0.92 | 1.576e-02 | 1.5 | 0.97 | 1.24 | 1.24 | 0.99 | 0.97 |
| Ar | 800 | 1.94e+14 | 0.82 | 4.112e-02 | 0.7 | 1.03 | 1.28 | 1.53 | 1.21 | 1.03 |
| Ar | 800 | 2.43e+14 | 0.74 | 1.065e-01 | 0.7 | 0.98 | 1.19 | 1.65 | 1.32 | 0.98 |
| Ar | 800 | 3.05e+14 | 0.66 | 2.322e-01 | 0.7 | 0.97 | 1.14 | 1.76 | 1.46 | 0.97 |
| Ar | 800 | 3.82e+14 | 0.59 | 4.096e-01 | 0.7 | 1.03 | 1.16 | 1.79 | 1.58 | 1.03 |
| Ar | 800 | 4.79e+14 | 0.52 | 6.424e-01 | 0.7 | 1.02 | 1.11 | 1.48 | 1.42 | 1.02 |
| Ar | 800 | 6.00e+14 | 0.47 | 8.189e-01 | 0.3 | 1.04 | 1.09 | 1.22 | 1.21 | 1.04 |
| Ar | 800 | 5.00e+13 | 1.62 | 2.644e-05 | 1.5 | 0.86 | 1.24 | 0.48 | 0.44 | 0.86 |
| Ar | 800 | 6.27e+13 | 1.45 | 9.644e-05 | 1.5 | 1.04 | 1.45 | 0.68 | 0.60 | 1.04 |
| Ar | 800 | 7.86e+13 | 1.30 | 4.337e-04 | 1.5 | 0.93 | 1.27 | 0.73 | 0.62 | 0.93 |
| Ar | 800 | 9.85e+13 | 1.16 | 1.592e-03 | 1.5 | 0.94 | 1.25 | 0.86 | 0.72 | 0.94 |
| Ar | 800 | 1.23e+14 | 1.03 | 5.247e-03 | 1.5 | 0.96 | 1.24 | 1.04 | 0.84 | 0.96 |
| Ar | 800 | 1.55e+14 | 0.92 | 1.644e-02 | 1.5 | 0.93 | 1.19 | 1.19 | 0.95 | 0.93 |
| Ar | 800 | 1.94e+14 | 0.82 | 4.092e-02 | 0.7 | 1.04 | 1.28 | 1.53 | 1.22 | 1.03 |
| Ar | 800 | 2.43e+14 | 0.74 | 1.103e-01 | 0.7 | 0.95 | 1.14 | 1.59 | 1.28 | 0.95 |
| Ar | 800 | 3.05e+14 | 0.66 | 2.181e-01 | 0.7 | 1.04 | 1.22 | 1.88 | 1.56 | 1.04 |
| Ar | 800 | 3.82e+14 | 0.59 | 4.240e-01 | 0.7 | 0.99 | 1.12 | 1.73 | 1.52 | 0.99 |
| Ar | 800 | 4.79e+14 | 0.52 | 6.343e-01 | 0.7 | 1.03 | 1.12 | 1.50 | 1.43 | 1.03 |
| Ar | 800 | 6.00e+14 | 0.47 | 8.374e-01 | 0.3 | 1.02 | 1.06 | 1.19 | 1.19 | 1.02 |
| Ar | 400 | 3.00e+13 | 4.19 | 8.627e-04 | 0.3 | 0.61 | 2.24 | 0.16 | 0.21 | 0.60 |
| Ar | 400 | 4.63e+13 | 3.37 | 2.345e-03 | 0.3 | 1.39 | 4.39 | 0.51 | 0.62 | 1.36 |
| Ar | 400 | 7.16e+13 | 2.71 | 1.631e-02 | 0.3 | 1.12 | 3.06 | 0.57 | 0.64 | 1.10 |
| Ar | 400 | 1.11e+14 | 2.18 | 1.844e-01 | 0.3 | 0.49 | 1.11 | 0.35 | 0.36 | 0.48 |
| Ar | 400 | 1.71e+14 | 1.76 | 3.528e-01 | 0.2 | 0.97 | 1.69 | 0.98 | 0.94 | 0.95 |
| Ar | 400 | 2.64e+14 | 1.41 | 7.362e-01 | 0.2 | 1.08 | 1.30 | 1.22 | 1.18 | 1.07 |
| Ar | 400 | 4.07e+14 | 1.14 | 9.373e-01 | 0.2 | 1.06 | 1.07 | 1.07 | 1.07 | 1.06 |
| Ar | 400 | 6.29e+14 | 0.92 | 9.660e-01 | 0.2 | 1.04 | 1.04 | 1.04 | 1.04 | 1.04 |
| Ar | 400 | 9.71e+14 | 0.74 | 9.800e-01 | 0.2 | 1.02 | 1.02 | 1.02 | 1.02 | 1.02 |
| Ar | 1030 | 2.39e+13 | 1.83 | 8.234e-08 | 4.1 | 1.00 | 1.36 | 0.39 | 0.37 | 1.00 |
| Ar | 1030 | 3.68e+13 | 1.47 | 3.752e-06 | 3.9 | 1.02 | 1.32 | 0.55 | 0.49 | 1.01 |
| Ar | 1030 | 5.69e+13 | 1.18 | 1.475e-04 | 3.9 | 0.80 | 1.00 | 0.60 | 0.50 | 0.80 |
| Ar | 1030 | 8.79e+13 | 0.95 | 2.311e-03 | 0.4 | 1.04 | 1.26 | 1.04 | 0.83 | 1.04 |
| Ar | 1030 | 1.36e+14 | 0.77 | 3.184e-02 | 0.4 | 0.98 | 1.15 | 1.29 | 1.02 | 0.98 |
| Ar | 1030 | 2.10e+14 | 0.62 | 2.560e-01 | 0.4 | 0.93 | 1.05 | 1.47 | 1.22 | 0.93 |
| Ar | 1030 | 3.24e+14 | 0.50 | 8.056e-01 | 0.4 | 0.98 | 1.03 | 1.21 | 1.17 | 0.98 |
| Ar | 1030 | 5.00e+14 | 0.40 | 9.899e-01 | 0.1 | 1.01 | 1.01 | 1.01 | 1.01 | 1.01 |
| He | 800 | 2.39e+14 | 0.93 | 4.049e-06 | 0.8 | 1.10 | 1.21 | 0.88 | 0.72 | 1.10 |
| He | 800 | 3.68e+14 | 0.75 | 1.004e-04 | 0.4 | 1.13 | 1.23 | 1.09 | 0.87 | 1.13 |
| He | 800 | 5.69e+14 | 0.60 | 1.561e-03 | 0.4 | 1.13 | 1.21 | 1.29 | 1.04 | 1.14 |
| He | 800 | 8.79e+14 | 0.48 | 1.541e-02 | 1.8 | 1.12 | 1.19 | 1.49 | 1.22 | 1.13 |
| He | 800 | 1.36e+15 | 0.39 | 9.651e-02 | 1.8 | 1.12 | 1.17 | 1.68 | 1.43 | 1.12 |
| He | 800 | 2.10e+15 | 0.31 | 3.598e-01 | 1.8 | 1.13 | 1.16 | 1.70 | 1.55 | 1.13 |
| Kr | 800 | 1.69e+13 | 2.63 | 5.529e-07 | 1.5 | 0.98 | 1.83 | 0.22 | 0.30 | 0.97 |
| Kr | 800 | 2.61e+13 | 2.12 | 1.155e-05 | 1.5 | 0.90 | 1.55 | 0.29 | 0.36 | 0.90 |
| Kr | 800 | 4.03e+13 | 1.71 | 1.148e-04 | 1.5 | 1.39 | 2.23 | 0.64 | 0.74 | 1.38 |
| Kr | 800 | 6.22e+13 | 1.37 | 1.297e-03 | 1.5 | 1.46 | 2.19 | 0.96 | 1.03 | 1.45 |
| Kr | 800 | 9.61e+13 | 1.10 | 1.297e-02 | 0.5 | 1.28 | 1.81 | 1.21 | 1.22 | 1.27 |
| Kr | 800 | 1.48e+14 | 0.89 | 9.102e-02 | 0.5 | 1.13 | 1.50 | 1.48 | 1.45 | 1.13 |
| Kr | 800 | 2.29e+14 | 0.72 | 3.608e-01 | 0.5 | 1.11 | 1.34 | 1.72 | 1.69 | 1.11 |
| Kr | 800 | 3.54e+14 | 0.58 | 7.626e-01 | 0.5 | 1.09 | 1.17 | 1.30 | 1.30 | 1.09 |
| Kr | 800 | 5.47e+14 | 0.46 | 9.510e-01 | 0.5 | 1.04 | 1.05 | 1.05 | 1.05 | 1.04 |
| Kr | 800 | 5.00e+13 | 1.53 | 4.361e-04 | 0.7 | 1.19 | 1.85 | 0.66 | 0.73 | 1.19 |
| Kr | 800 | 6.27e+13 | 1.37 | 1.272e-03 | 0.7 | 1.42 | 2.13 | 0.94 | 1.01 | 1.41 |
| Kr | 800 | 7.86e+13 | 1.22 | 4.825e-03 | 0.7 | 1.20 | 1.75 | 0.96 | 1.00 | 1.20 |
| Kr | 800 | 9.85e+13 | 1.09 | 1.392e-02 | 0.7 | 1.23 | 1.73 | 1.18 | 1.19 | 1.22 |
| Kr | 800 | 1.23e+14 | 0.97 | 3.994e-02 | 0.7 | 1.15 | 1.57 | 1.32 | 1.30 | 1.15 |
| Kr | 800 | 1.55e+14 | 0.87 | 9.607e-02 | 0.7 | 1.15 | 1.52 | 1.56 | 1.52 | 1.15 |
| Kr | 800 | 1.94e+14 | 0.78 | 1.831e-01 | 0.7 | 1.29 | 1.62 | 1.96 | 1.91 | 1.29 |
| Kr | 800 | 2.43e+14 | 0.69 | 3.639e-01 | 0.7 | 1.19 | 1.42 | 1.86 | 1.82 | 1.19 |
| Kr | 800 | 3.05e+14 | 0.62 | 5.918e-01 | 0.7 | 1.13 | 1.27 | 1.57 | 1.56 | 1.13 |
| Kr | 800 | 3.82e+14 | 0.55 | 7.732e-01 | 0.7 | 1.11 | 1.18 | 1.29 | 1.29 | 1.11 |
| Kr | 800 | 4.79e+14 | 0.49 | 9.140e-01 | 0.7 | 1.05 | 1.07 | 1.09 | 1.09 | 1.05 |
| Kr | 800 | 6.00e+14 | 0.44 | 9.586e-01 | 0.3 | 1.04 | 1.04 | 1.04 | 1.04 | 1.04 |
| Kr | 800 | 5.00e+13 | 1.53 | 4.435e-04 | 0.7 | 1.17 | 1.81 | 0.64 | 0.72 | 1.17 |
| Kr | 800 | 6.27e+13 | 1.37 | 1.169e-03 | 0.7 | 1.54 | 2.32 | 1.03 | 1.10 | 1.54 |
| Kr | 800 | 7.86e+13 | 1.22 | 5.270e-03 | 0.7 | 1.10 | 1.60 | 0.88 | 0.91 | 1.10 |
| Kr | 800 | 9.85e+13 | 1.09 | 1.469e-02 | 0.7 | 1.16 | 1.64 | 1.12 | 1.13 | 1.16 |
| Kr | 800 | 1.23e+14 | 0.97 | 3.582e-02 | 0.7 | 1.28 | 1.75 | 1.47 | 1.45 | 1.28 |
| Kr | 800 | 1.55e+14 | 0.87 | 9.607e-02 | 0.7 | 1.15 | 1.52 | 1.56 | 1.52 | 1.15 |
| Kr | 800 | 1.94e+14 | 0.78 | 1.934e-01 | 0.7 | 1.22 | 1.54 | 1.86 | 1.81 | 1.22 |
| Kr | 800 | 2.43e+14 | 0.69 | 3.946e-01 | 0.7 | 1.10 | 1.31 | 1.71 | 1.68 | 1.10 |
| Kr | 800 | 3.05e+14 | 0.62 | 5.932e-01 | 0.7 | 1.13 | 1.26 | 1.56 | 1.55 | 1.13 |
| Kr | 800 | 3.82e+14 | 0.55 | 8.069e-01 | 0.7 | 1.07 | 1.13 | 1.24 | 1.23 | 1.07 |
| Kr | 800 | 4.79e+14 | 0.49 | 8.996e-01 | 0.7 | 1.07 | 1.09 | 1.11 | 1.11 | 1.07 |
| Kr | 800 | 6.00e+14 | 0.44 | 9.602e-01 | 0.3 | 1.03 | 1.04 | 1.04 | 1.04 | 1.03 |
| Xe | 800 | 1.69e+13 | 2.45 | 4.377e-05 | 1.3 | 1.03 | 2.29 | 0.22 | 0.40 | 1.03 |
| Xe | 800 | 2.61e+13 | 1.97 | 6.108e-04 | 1.3 | 0.85 | 1.70 | 0.27 | 0.45 | 0.84 |
| Xe | 800 | 4.03e+13 | 1.59 | 4.922e-03 | 1.3 | 0.98 | 1.78 | 0.46 | 0.70 | 0.97 |
| Xe | 800 | 6.22e+13 | 1.28 | 4.242e-02 | 1.3 | 0.81 | 1.35 | 0.58 | 0.81 | 0.81 |
| Xe | 800 | 9.61e+13 | 1.03 | 1.811e-01 | 0.3 | 0.98 | 1.45 | 1.02 | 1.33 | 0.97 |
| Xe | 800 | 1.48e+14 | 0.83 | 5.203e-01 | 0.3 | 1.07 | 1.34 | 1.39 | 1.56 | 1.07 |
| Xe | 800 | 2.29e+14 | 0.67 | 8.301e-01 | 0.3 | 1.11 | 1.17 | 1.20 | 1.20 | 1.11 |
| Xe | 800 | 5.00e+13 | 1.42 | 1.333e-02 | 0.8 | 0.92 | 1.60 | 0.53 | 0.77 | 0.91 |
| Xe | 800 | 6.27e+13 | 1.27 | 4.052e-02 | 1.2 | 0.81 | 1.34 | 0.57 | 0.81 | 0.80 |
| Xe | 800 | 7.86e+13 | 1.14 | 8.465e-02 | 1.2 | 0.95 | 1.49 | 0.83 | 1.13 | 0.95 |
| Xe | 800 | 9.85e+13 | 1.02 | 1.833e-01 | 1.2 | 0.97 | 1.43 | 1.03 | 1.34 | 0.97 |
| Xe | 800 | 1.23e+14 | 0.91 | 3.377e-01 | 1.2 | 1.02 | 1.39 | 1.27 | 1.54 | 1.02 |
| Xe | 800 | 1.55e+14 | 0.81 | 5.704e-01 | 1.2 | 1.00 | 1.24 | 1.31 | 1.46 | 1.00 |
| Xe | 800 | 1.94e+14 | 0.72 | 7.053e-01 | 1.2 | 1.12 | 1.26 | 1.36 | 1.39 | 1.12 |
| Xe | 800 | 2.43e+14 | 0.65 | 9.035e-01 | 1.2 | 1.03 | 1.08 | 1.11 | 1.11 | 1.03 |
| Xe | 800 | 3.05e+14 | 0.58 | 9.181e-01 | 1.2 | 1.07 | 1.09 | 1.09 | 1.09 | 1.07 |
| Xe | 800 | 5.00e+13 | 1.42 | 1.378e-02 | 0.8 | 0.89 | 1.54 | 0.51 | 0.75 | 0.88 |
| Xe | 800 | 6.27e+13 | 1.27 | 3.845e-02 | 1.2 | 0.85 | 1.41 | 0.61 | 0.85 | 0.85 |
| Xe | 800 | 7.86e+13 | 1.14 | 9.063e-02 | 1.2 | 0.89 | 1.40 | 0.78 | 1.05 | 0.88 |
| Xe | 800 | 9.85e+13 | 1.02 | 1.508e-01 | 1.2 | 1.18 | 1.74 | 1.26 | 1.62 | 1.17 |
| Xe | 800 | 1.23e+14 | 0.91 | 2.908e-01 | 1.2 | 1.19 | 1.62 | 1.47 | 1.79 | 1.18 |
| Xe | 800 | 1.55e+14 | 0.81 | 5.044e-01 | 1.2 | 1.13 | 1.40 | 1.48 | 1.65 | 1.13 |
| Xe | 800 | 1.94e+14 | 0.72 | 7.714e-01 | 1.2 | 1.03 | 1.15 | 1.24 | 1.27 | 1.02 |
| Xe | 800 | 2.43e+14 | 0.65 | 8.530e-01 | 1.2 | 1.09 | 1.14 | 1.17 | 1.17 | 1.09 |
| Xe | 800 | 3.05e+14 | 0.58 | 9.344e-01 | 1.2 | 1.06 | 1.07 | 1.07 | 1.07 | 1.06 |
<!-- ionisation-rate-tdse:end -->

```{figure} /auto_gallery/ionisation/images/sphx_glr_plot_ionisation_frequency_001.*
:target: ../auto_gallery/ionisation/plot_ionisation_frequency.html
:alt: The three Coulomb factors against the calculations of helium and argon

The three Coulomb factors against the one-electron yields of argon and helium and the
two-electron rates of helium.
```

The following comparisons explain why the former static Coulomb factor was
replaced (TC-812, TC-848). "The static factor" in them is that of
`one_electron_static`; the current `one_electron` column above includes the empirical
frequency-dependent factor, named `frequency` where a bullet compares it.

- **Earlier analytic Coulomb factors.** The static factor of `one_electron_static` is
  1.08 to 1.67 times above them for argon, the most at the largest Keldysh parameter,
  and 1.16 to 1.23 times for helium. Ilkov et al.'s and Popruzhenko et al.'s cross them
  near $\gamma = 1$: 2 to 3 times too low at $\gamma \approx 2$, up to 1.7 times too
  high at $\gamma \approx 0.6$. The static factor has the smallest error over the
  range (TC-653, TC-812).
- **At 400 nm the excess is larger.** The static factor is 1.4 to 3.9 times the
  calculation at Keldysh parameters 1.4 to 4.2, the excess Nambu et al. measured at
  400 nm, now for one electron: it is in the rate's dependence on frequency, not in the
  atom. Popruzhenko et al.'s factor alone is 3 to 4 times low there. Half of its
  dependence on the Keldysh parameter (in the logarithm) fits 400 nm (0.96 on average)
  but is low at 800 nm (0.65 to 0.73); a third of it fits the three argon calculations
  together to a factor of 1.3 (TC-824).
- **At 1030 nm the static factor is nearly right.** With the same 10 fs pulse as at
  400 nm, the static-factor model is 1.02 to 1.35 times the calculation at Keldysh parameters 0.4
  to 1.8; at $\gamma \approx 1.8$ it is 1.1 times at 1030 nm, 1.6 at 800 nm (6 fs) and
  2.0 at 400 nm. The fraction of Popruzhenko et al.'s dependence that fits each set on
  average is 0.15 at 1030 nm, 0.16 (42 fs) and 0.24 (6 fs) at 800 nm, and 0.48 at 400
  nm: it grows with the frequency, but the pulse length matters as well. A single
  fraction, 0.25 to 0.3, fits the four argon sets to a factor of 1.28 (TC-826).
- **One pulse at three wavelengths.** With the 10 fs pulse at 800 nm as well, the
  static-factor model at $\gamma = 1.8$ is 1.1 times the calculation at 1030 nm, 1.7 at 800 nm and
  1.9 at 400 nm, so the excess grows with the frequency at one pulse length; at 800
  nm the 10 fs excess (1.75 to 1.05 over $\gamma$ 2.4 to 0.5) is a little above the
  6 fs one. The fraction that fits is 0.15, 0.29 and 0.48, about five times $s =
  n^*\omega/2I_p$ at each wavelength (TC-830).
- **Wallace et al.'s pulse** (5.5 fs, carrier-envelope phases 0 and $\pi/2$): the
  phase changes the yield by up to 6 %, oscillating with the intensity. The static
  factor is 1.06 to 1.45 times the calculation, the `frequency` factor 0.86
  to 1.04. Put through Wallace et al.'s test in place of a model
  ({func}`~optmat.gas.ionisation_validation.tdse_intensity_scale_test`), the
  calculation itself matches their measured yield at an intensity scale of 1.002, within
  its calibration of 1.3 %: the measurement agrees with exact one-electron physics. The
  test hardly sees the rate's magnitude: a factor of 0.8 on the calculation's rate (an
  atom's, TC-830) moves the scale to 1.016, and the static and frequency-dependent models sit
  either side of the calculation, at 0.97 and 1.02 for one electron (TC-836).
- **A 30 fs pulse at 400 and 1030 nm** (the 10 fs sets' intensities): the pulse three
  times longer ionises 2.3 to 5.4 times more at 400 nm and 2.5 to 3.8 times at 1030 nm
  away from saturation, about the 3 a rate gives for a Gaussian three times longer:
  the dependence on the pulse length beyond that is small. The `frequency`
  factor, fitted to the 10 fs sets, is 0.87 (400 nm) and 0.95 (1030 nm) times these on
  geometric average, the static factor 2.0 and 1.18 (TC-844).
- **Krypton and xenon** (A2's 6 fs pulse, one electron in Zhang et al.'s potentials;
  TC-845): the static factor is 1.74 (Kr) and 1.61 (Xe) times the calculation on
  geometric average, the `frequency` factor 1.16 and 0.94. Krypton's 1.16 is
  its static normalisation, not the Coulomb factor: the one-electron model is 1.08 to
  1.17 times tRecX's static rates of the same potential (Bahl et al.'s: 1.10 to 1.22);
  xenon's is 0.84 to 1.04. Zhao et al.'s krypton, 2.2 times below the static factor
  and 1.5 below the frequency-dependent factor, used another potential (Garvey et al.'s); with Zhang et
  al.'s the excess is argon's and xenon's. With the frequency-dependent factor, one electron of xenon
  is not ionised more, against argon's, than the model says: the measured excess of
  xenon over argon (TC-819) is not in the one-electron physics.
- **Krypton in Wallace et al.'s pulse** (TC-846): through their test the exact
  one-electron calculation of krypton needs the intensities scaled by 1.045, beyond the
  calibration of 1.5 %, where argon's matched at 1.002: the measured krypton atom
  ionises more than one $m = 0$ electron. A factor of about 1.4 on the calculation's
  rate matches it (1.47: 0.989): the krypton atom ionises like about 1.4 electrons
  without the spin-orbit splitting, where argon's ionises like 0.8 to 1 (TC-830,
  TC-842). (Santra et al.'s krypton rate, 1.34 times tRecX's one electron at 0.10 au,
  was first read as a second, all-electron indication; it is an independent-electron
  sum in another potential; TC-851.) The phase changes krypton's yield by up to 10 %, against argon's 6 %.
- **Xenon in Wallace et al.'s pulse** (TC-849): the exact one-electron calculation of
  xenon passes their test, at 1.019 within the calibration of 2.5 %, as argon's does
  (1.002) and krypton's does not (1.045). A factor of 1.09 on its rate gives 1.000; the
  test allows about 0.93 to 1.27. The measured xenon atom ionises like about one
  electron *without the spin-orbit splitting*, at the ionisation energy of
  $^2P_{3/2}$. The model's $N_\text{eff}$ multiplies the rate weighted over the two
  levels of the ion (2/3 at $^2P_{3/2}$, 1/3 at $^2P_{1/2}$, 1.3 eV higher for xenon),
  which is 0.71 to 0.79 times that rate for xenon from 5e13 to 2e14 W/cm² (krypton
  0.78 to 0.85, argon 0.91 to 0.95). With one electron the default was at 1.133; with
  xenon's $N_\text{eff}$ set to 1.5 (about $1.09/(0.98 \times 0.75)$, the model's one
  electron being 0.96 to 1.00 times the calculation here; TC-850) it is at 1.026,
  where argon's (1.033) and krypton's (1.023) are. That moved Nambu et al.'s xenon at
  400 nm from 1.11 to 0.98, at 800 nm from 1.00 to 0.92, Nie et al.'s from 0.87 to
  0.80, and Wahlstrand et al.'s from 1.9–4.6 to 2.9–6.8. One electron of xenon in the frequency-dependent factor is
  0.80 to 1.18 times the calculation in this pulse (geometric means 0.96 and 1.00 at
  the two phases), the static factor 1.24 to 1.73: the factor holds in a second
  pulse. The phase changes xenon's yield by up to 18 %.
- **An independent code agrees.** Zhao et al.'s (2016) one-electron calculations (15
  fs) put the static factor 1.23 times high for argon at 800 nm, as tRecX's do, 1.75 to
  2.2 for krypton and xenon at 800 nm, 1.4 to 1.6 for helium and neon at 400 nm, and,
  with Tong and Lin's static correction, 1.09, 1.17, 1.28 and 1.53 for hydrogen at
  1200, 1000, 800 and 600 nm (TC-827): at equal Keldysh parameter the excess still
  changes with the wavelength. It orders with $s = n^*\omega/2I_p$, the parameter of the static
  factor's range of validity (Perelomov, Popov & Kuznetsov 1968; TC-670, TC-825).
- **The excess is in the dependence on frequency, not in the static limit.** The same
  model is 0.94 to 1.01 times tRecX's static rates of the same argon electron ({doc}`ionisation-static`),
  and Zhang et al.'s.
- The rates of the $m = 1$ electrons of argon, static and in the 6 fs pulse, are curated
  with them (`travers_2026_trecx/static_Ar_m1`, `ionised_gauss6_Ar_m1`): 2 to 22 % of the
  $m = 0$ rate in a static field from 0.05 to 0.2 au, 2 to 39 % of its yield in the
  pulse.

## One-electron calculation of argon, 800 nm, 42 fs

The same model of argon (3p, m = 0) in the pulse of Wahlstrand et al.'s measurement: a
Gaussian of 42 fs, 17 peak intensities (`travers_2026_trecx`, TC-702). The yield is
the ionised fraction, the population of the bound states up to n = 10 subtracted; the
excited states hold 1 to 30 % of it, most at the resonant intensities. Its uncertainty
from convergence runs is 0.1 to 2.6 %.

<!-- ionisation-rate-tdse42:start -->
| I (W/cm2) | gamma | TDSE: ionised fraction | `one_electron_static` | `one_electron` | `default` | TDSE / measured (the atom) |
|---|---|---|---|---|---|---|
| 4.00e+13 | 1.82 | 4.490e-05 | 1.21 | 0.82 | 0.75 | 1.28 |
| 4.50e+13 | 1.71 | 1.101e-04 | 1.11 | 0.77 | 0.69 | 1.71 |
| 5.00e+13 | 1.62 | 1.948e-04 | 1.27 | 0.89 | 0.81 | 1.76 |
| 5.50e+13 | 1.55 | 2.558e-04 | 1.80 | 1.27 | 1.16 | 1.41 |
| 6.20e+13 | 1.46 | 9.175e-04 | 1.08 | 0.77 | 0.71 | 2.74 |
| 6.60e+13 | 1.41 | 1.466e-03 | 0.99 | 0.72 | 0.66 | 3.18 |
| 7.00e+13 | 1.37 | 1.689e-03 | 1.24 | 0.90 | 0.82 | 2.70 |
| 7.40e+13 | 1.33 | 2.157e-03 | 1.35 | 0.99 | 0.91 | 2.60 |
| 7.80e+13 | 1.30 | 2.730e-03 | 1.46 | 1.07 | 0.99 | 2.51 |
| 8.20e+13 | 1.27 | 4.201e-03 | 1.27 | 0.94 | 0.87 | 2.98 |
| 8.70e+13 | 1.23 | 6.802e-03 | 1.10 | 0.82 | 0.76 | 3.56 |
| 9.20e+13 | 1.20 | 8.767e-03 | 1.17 | 0.88 | 0.81 | 3.45 |
| 9.70e+13 | 1.17 | 1.186e-02 | 1.16 | 0.87 | 0.81 | 3.55 |
| 1.05e+14 | 1.12 | 1.480e-02 | 1.42 | 1.08 | 1.00 | 2.95 |
| 1.12e+14 | 1.09 | 2.134e-02 | 1.38 | 1.06 | 0.99 | 3.05 |
| 1.20e+14 | 1.05 | 4.048e-02 | 1.04 | 0.80 | 0.75 | 4.06 |
| 1.35e+14 | 0.99 | 6.226e-02 | 1.20 | 0.94 | 0.88 | 3.41 |
<!-- ionisation-rate-tdse42:end -->

- **The models follow the calculation within its resonance structure.** The
  `one_electron_static` model is 0.99 to 1.8 times the calculation, 1.24 on geometric
  average, the most just past the channel closings: the calculated yield against
  intensity is not smooth (its local slope swings between 2.4 and 10.7), which no rate
  is.
- **One electron already exceeds the measured yield of the atom**, by 1.3 to 4.1 times,
  and excitation cannot explain it. The disagreement with Wahlstrand et al. is not in
  the rate's dependence on frequency or pulse length, nor in counting Rydberg atoms as
  ionised: it is between one active electron and the measurement (TC-660, TC-702).

## One-electron calculations of hydrogen

Atomic hydrogen is the cleanest test of the Coulomb factor: its one-electron
calculations are exact, and its static rate is known. Four sets are curated: Zhao et
al. (2016; 15 fs Gaussian, 600 to 1200 nm), Le et al. (2016; 12- and 6-cycle
$\sin^2$ pulses at 400, 800 and 1600 nm) and Bauer and Deng (2018; a 2-10-2 cycle
flat top at 400 and 800 nm, and against the frequency at two fixed fields), read from
the figures (TC-825, TC-834, TC-835), and optmat's own tRecX calculations (15 fs
Gaussian at 400 and 800 nm, converged to 0.06 %; TC-843). The first two come from one
group; Bauer and Deng's and tRecX's are independent of it. At 800 nm, with the same
pulse, Zhao et al.'s read from their figure are 0.72 to 1.17 times tRecX's, the spread
that of the resonances between tRecX's coarser intensities: the reading and the
pulse's convention are confirmed. The model is `one_electron`'s form for the 1s electron, with
Tong and Lin's $\alpha = 6$ for hydrogen
({func}`~optmat.gas.ionisation_validation.hydrogen_model`; TC-833), on each
calculation's envelope. Model over calculation, geometric means of the yields below
0.75 in three ranges of the Keldysh parameter:

<!-- ionisation-rate-hydrogen:start -->
| Calculation | Wavelength or frequency | static: gamma 0 to 1 | static: gamma 1 to 2 | static: gamma 2 and above | frequency: gamma 0 to 1 | frequency: gamma 1 to 2 | frequency: gamma 2 and above | popruzhenko: gamma 0 to 1 | popruzhenko: gamma 1 to 2 | popruzhenko: gamma 2 and above |
|---|---|---|---|---|---|---|---|---|---|---|
| `zhao_2016/fig2a_tdse_H_600nm` | 600 nm |  | 1.47 | 1.95 |  | 0.81 | 0.83 |  | 0.33 | 0.23 |
| `zhao_2016/fig2b_tdse_H_800nm` | 800 nm | 1.28 | 1.32 | 1.01 | 0.93 | 0.86 | 0.56 | 0.44 | 0.32 | 0.14 |
| `zhao_2016/fig2c_tdse_H_1000nm` | 1000 nm | 1.16 | 1.24 |  | 0.92 | 0.88 |  | 0.44 | 0.30 |  |
| `zhao_2016/fig2d_tdse_H_1200nm` | 1200 nm | 1.09 | 1.13 |  | 0.92 | 0.87 |  | 0.45 | 0.30 |  |
| `travers_2026_trecx/ionised_gauss15_400nm_H` | 400 nm |  | 2.19 | 1.65 |  | 1.07 | 0.47 |  | 0.59 | 0.20 |
| `travers_2026_trecx/ionised_gauss15_800nm_H` | 800 nm | 1.19 | 1.21 |  | 0.96 | 0.81 |  | 0.56 | 0.32 |  |
| `le_2016/fig12a_tdse_H_400nm` | 400 nm |  | 3.72 | 4.44 |  | 1.44 | 0.81 |  | 0.73 | 0.26 |
| `le_2016/fig12b_tdse_H_800nm` | 800 nm | 1.32 | 1.41 | 1.35 | 0.96 | 0.93 | 0.68 | 0.45 | 0.35 | 0.14 |
| `le_2016/fig12c_tdse_H_1600nm` | 1600 nm | 1.05 | 1.11 |  | 0.93 | 0.93 |  | 0.48 | 0.33 |  |
| `bauer_2018/fig1_tdse_H_800nm` | 800 nm | 1.27 | 1.36 |  | 1.03 | 0.95 |  | 0.59 | 0.40 |  |
| `bauer_2018/fig5a_tdse_H_800nm` | 800 nm | 1.20 | 1.17 | 0.89 | 0.96 | 0.76 | 0.46 | 0.54 | 0.28 | 0.10 |
| `bauer_2018/fig6a_tdse_H_400nm` | 400 nm |  | 2.63 | 1.67 |  | 1.42 | 0.37 |  | 0.83 | 0.13 |
| `bauer_2018/fig7_tdse_H_spectrum_F002` | 0.04 to 0.25 au |  | 0.38 | 13.10 |  | 0.26 | 0.66 |  | 0.06 | 0.38 |
| `bauer_2018/fig8_tdse_H_spectrum_F0063` | 0.04 to 0.25 au | 1.14 | 2.22 | 2.90 | 0.94 | 1.26 | 0.75 | 0.52 | 0.67 | 0.78 |
<!-- ionisation-rate-hydrogen:end -->

- **Below $\gamma = 1$** every calculation agrees: the static factor is 1.05 to 1.32
  times them, the `frequency` factor 0.92 to 1.03, Popruzhenko et al.'s 0.44
  to 0.59.
- **From $\gamma = 1$ to 2** the static factor is 1.1 to 1.5 times them at 600 nm and
  longer, 2.2 to 3.7 at 400 nm; the frequency-dependent factor 0.76 to 0.95, and 1.07 to 1.44 at 400
  nm.
- **Above $\gamma = 2$**, in the multiphoton regime, the yields carry resonances
  (Bauer and Deng's excitation reaches 34 % at 400 nm) and no form follows them: the
  static factor is 0.9 to 4.4 times them, the frequency-dependent factor 0.37 to 0.83, too low in
  every set. Against the frequency at 0.02 au (Bauer and Deng's Fig. 7, $\gamma = 2$
  and above at every frequency compared) the static factor is 13 times high and every
  form fails.
- **Fitted to argon, the frequency-dependent factor's constant overcorrects hydrogen** at $\gamma > 1$
  (TC-827: hydrogen alone would take $c$ = 2.6 to 3.7).
- **Bauer and Deng's two 800 nm sets differ**: their Fig. 5(a), whose caption says
  it shows the calculations of Fig. 1, puts its points on another intensity grid and up
  to 20 % higher; their Fig. 8 sides with Fig. 5(a). Both are curated (TC-834).

## The frequency-dependent Coulomb factor

The one-electron calculations above show the static factor too high at $\gamma > 1$ by
an amount that grows with the frequency at fixed $\gamma$; Popruzhenko et al.'s factor,
which depends on $\gamma$ alone, is too low. Every first-order Coulomb correction of
the theory depends on $\gamma$ alone (Popruzhenko et al. 2008; Karnakov et al. 2015, eq.
5.8; the analytical R-matrix theory), and none in the literature reviewed gives the
calculations' dependence on the frequency. The factor, the default's since TC-848, is
empirical:

$$
Q = \left(\frac{2F_0}{F}\right)^{2n^*}\left(1 + \frac{2\gamma}{e}\right)^{-2n^* c s},
\qquad s = \frac{n^*\omega}{2I_p},
$$

with $s$ the parameter of the static factor's range of validity (Perelomov, Popov and
Kuznetsov 1968) and $c = 5.25$ fitted to the three 10 fs calculations of argon (400,
800, 1030 nm; yields below 0.75). It tends to the static factor as $\omega \to 0$.
`default` and `one_electron` carry it (`coulomb = "frequency"`). Model over calculation
or measurement, the static factor (`default_static`, `one_electron_static`) against it
(TC-831). This table records the comparison before the Kr and Xe electron counts
were revised: it uses one electron, or the atomic channel weights with argon's
one-electron normalisation. It is not a table of the current atomic defaults:

| Test | Static | Frequency-dependent |
|---|---|---|
| tRecX argon, three 10 fs sets (fitted): rms factor | 1.78 | 1.17 |
| tRecX argon, 800 nm, 6 fs: geometric mean | 1.34 | 1.00 |
| tRecX argon, 800 nm, 42 fs | 1.24 | 0.91 |
| tRecX helium, 800 nm, 6 fs | 1.20 | 1.12 |
| Zhao et al., He and Ne at 400 nm | 1.42, 1.57 | 1.15, 1.14 |
| Zhao et al., Ar, Kr, Xe at 800 nm | 1.23, 2.22, 1.75 | 0.95, 1.53, 0.99 |
| Zhao et al., H at 600 to 1200 nm (Tong and Lin's alpha = 6) | 1.52 to 1.09 | 0.84 to 0.91 |
| Parker et al.'s helium at 390 nm | 1.27 to 1.45 | 0.85 to 1.28 |
| Nambu et al., 400 nm, Ar, Kr, Xe: intensity scale | 0.50 to 0.73 | 0.90 to 1.11 (Xe circular 1.41) |
| Nambu et al. and Nie et al., 800 nm, 16 fits: rms intensity scale | 1.11 | 1.08 |
| Wallace et al., Ar, Kr, Xe: intensity scale | 0.98, 1.00, 1.00 | 1.03, 1.07, 1.13 |
| Wahlstrand et al.: rate over measured | 1.7 to 7.9 | 1.0 to 4.6 |

It brings every calculation within 0.84 to 1.15 on average except krypton's, which
stays 1.5 times high (Zhao et al.'s potential; with Zhang et al.'s, tRecX gives 1.16, its
static normalisation, TC-845), and
makes hydrogen's ratio flat in the wavelength. It removes most
of the excess at 400 nm. Against Wallace et al.'s intensity scale it is 3 to 13 % beyond
the calibration of 1.3 to 2.5 % (Ar, Kr, Xe). For argon that is no evidence against it:
the exact one-electron calculation in their pulse matches the measurement at 1.002,
and the frequency-dependent factor is within 2 % of the calculation's scale where the static factor is
within 3 % on the other side (TC-836). For krypton the exact calculation itself is
4.5 % beyond the calibration: the measured atom ionises like about 1.35 electrons,
which the frequency-dependent factor, with argon's one electron, leaves out (TC-846). For xenon the
exact calculation is within the calibration, at 1.019 (TC-849). The fit is to one
gas, one pulse and one carrier-envelope phase; $c$ for hydrogen alone would be 2.6 to
3.7. Adopted as the default in TC-848, with krypton's $N_\text{eff}$ raised to 1.35;
with it the default matches Wallace et al.'s argon at 1.033 and krypton at 1.023, and,
with xenon's $N_\text{eff}$ of 1.5 (TC-850), xenon at 1.026.
