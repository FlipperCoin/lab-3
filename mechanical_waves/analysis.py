from uncertainties.umath import sqrt, floor
import numpy as np
import matplotlib.pyplot as plt # for plotting figures and setting their properties
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
    freq_step = unit["d"] * sqrt(k / unit["I"]) / (2 * unit["L"])
    max_n = floor(max_freq / freq_step)
    unit["f_forced"] = np.array([n * freq_step for n in range(1, max_n)])

    # free end
    freq = lambda n: unit["d"]*sqrt(k/unit["I"])*(2*n-1)/(4*unit["L"])
    max_n = floor(0.5*(8*unit["L"]*(1/unit["d"])*sqrt(unit["I"]/k)+1))
    unit["f_free"] = np.array([freq(n) for n in range(1, max_n)])


df_forced = max(np.abs(wide_unit["f_forced"] - [value.nominal_value for value in forced_modes.values()]))

plt.errorbar(forced_modes.keys(), [value.nominal_value for value in forced_modes.values()], yerr=[value.std_dev for value in forced_modes.values()], fmt="ro", label="measurements")
plt.plot(range(1, len(wide_unit["f_forced"])+1), [val.nominal_value for val in wide_unit["f_forced"]], label="prediction")
plt.xlabel("mode number $[1]$")
plt.ylabel(r"frequency $[\frac{1}{sec}]$")
plt.legend()
plt.grid()
plt.show()

df_free = max(np.abs(wide_unit["f_free"][1:] - [value.nominal_value for value in free_modes.values()]))

plt.errorbar(free_modes.keys(),  [value.nominal_value for value in free_modes.values()], yerr=[value.std_dev for value in free_modes.values()], fmt="ro", label="measurements")
plt.plot(range(1, len(wide_unit["f_free"])+1), [val.nominal_value for val in wide_unit["f_free"]], label="prediction")
plt.xlabel("mode number $[1]$")
plt.ylabel(r"frequency $[\frac{1}{sec}]$")
plt.grid()
plt.legend()
plt.show()


