CoolProp is loaded without the superancillaries it builds for all of its fluids, unless
`COOLPROP_DISABLE_SUPERANCILLARIES_ENTIRELY` is already set or CoolProp was imported
first: the first real-gas density takes 0.07 s instead of 0.46 s. No density at a given
state changed and no recommended value moved, but the identity of every equation of
state, and of the mixture, transport and air models built on one, changes: the
equation-of-state `EVALUATION_VERSION` is now 2. A gas's critical temperature is now its
equation's stated value (oxygen 154.581 K instead of 154.599 K; others within 2e-5), and
saturation pressures, so where a gas's domain ends, move by at most 7e-11 more than 1 K
below the critical point and by up to 7e-8 within it (TC-1096).
