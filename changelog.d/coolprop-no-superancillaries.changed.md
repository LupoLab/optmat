CoolProp is loaded without the superancillaries it builds for all of its fluids, unless
`COOLPROP_DISABLE_SUPERANCILLARIES_ENTIRELY` is already set or CoolProp was imported
first: the first real-gas density takes 0.07 s instead of 0.46 s. No density changed and
no recommended value moved; saturation pressures, and so where a gas's domain ends, move
by at most 3e-11 (oxygen at its triple point: 2e-9) (TC-1096).
