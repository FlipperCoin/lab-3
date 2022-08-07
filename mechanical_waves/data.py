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

transition_unit = {
    "n": ufloat(47, 0, "1")
}

k = ufloat((0.49 * 22.5e-2) / 0.11, 0, "?")



