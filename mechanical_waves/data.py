from uncertainties import ufloat

wide_unit = {
    "l": ufloat(45.6e-2, 0, "m"),
    "m": ufloat(43.2e-3, 0, "kg/bar"),
    "d": ufloat(1.27e-2, 0, "m/bar"),
    "n": ufloat(72, 0, "1"),
    "L": None,
    "f_forced": None,
    "f_free": None
}

narrow_unit = {
    "l": ufloat(22.8e-2, 0, "m"),
    "m": ufloat(21.9e-3, 0, "kg/bar"),
    "d": ufloat(1.27e-2, 0, "m/bar"),
    "n": ufloat(72, 0, "1"),
    "L": None,
    "f_forced": None,
    "f_free": None
}

freq_meas_err = 0.05

forced_modes = {
    1: ufloat(0.3, freq_meas_err, "Hz"),
    2: ufloat(0.538, freq_meas_err, "Hz"),
    3: ufloat(0.789, freq_meas_err, "Hz"),
    4: ufloat(1.119, freq_meas_err, "Hz"),
    5: ufloat(1.353, freq_meas_err, "Hz"),
    6: ufloat(1.579, freq_meas_err, "Hz")
}

free_modes = {
    2: ufloat(0.508, freq_meas_err, "Hz"),
    3: ufloat(0.779, freq_meas_err, "Hz"),
    4: ufloat(1, freq_meas_err, "Hz"),
    5: ufloat(1.275, freq_meas_err, "Hz"),
    6: ufloat(1.540, freq_meas_err, "Hz"),
    7: ufloat(1.818, freq_meas_err, "Hz")
}

transition_unit = {
    "n": ufloat(47, 0, "1")
}

k = ufloat((0.49 * 22.5e-2) / 0.11, 0, "?")



