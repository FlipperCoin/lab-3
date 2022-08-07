from uncertainties.umath import sqrt, floor
import numpy as np
from data import *

def calc_bar_inertia(m, l):
    return m * (l ** 2) / 12

def calc_unit_impedance(k, i):
    return sqrt(k * i)

max_freq = ufloat(2, 0, "Hz")

for unit in [wide_unit, narrow_unit]:
    unit["L"] = unit["d"] * unit["n"]
    unit["I"] = calc_bar_inertia(unit["m"], unit["l"])
    unit["Z"] = calc_unit_impedance(k, unit["I"])

    # forced end
    freq_step = sqrt(k / unit["I"]) / (2 * unit["L"])
    max_n = floor(max_freq / freq_step)
    unit["f_forced"] = np.array([n * freq_step for n in range(1, max_n)])

    # free end
    freq = lambda n: sqrt(k/unit["I"])*(2*n-1)/(4*unit["L"])
    max_n = floor(0.5*(8*unit["L"]*sqrt(unit["L"]/k)+1))
    unit["f_free"] = np.array([freq(n) for n in range(1, max_n)])



