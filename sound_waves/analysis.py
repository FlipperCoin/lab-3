import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from uncertainties.umath import sqrt
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.constants import R, Avogadro, proton_mass, neutron_mass

#%% part 0 - prep

f = unumpy.uarray([], 0)
A = unumpy.uarray([], 0)

plt.errorbar(noms(f), noms(A), xerr=devs(f), yerr=devs(A), fmt="bo", label="data")
plt.grid()
plt.xlabel(r"f $\left[ \frac{1}{sec} \right]$")
plt.ylabel(r"A $\left[ m \right]$")
plt.show()

res_idx = np.argmax(noms(A))
f_res, A_res = f[res_idx], A[res_idx]

#%% part 1 - sound phase velocity in air

x = unumpy.uarray([], 0)
phi = unumpy.uarray([], 0)
reg = linregress(noms(x), noms(phi))

plt.errorbar(noms(x), noms(phi), xerr=devs(x), yerr=devs(phi), fmt="bo", label="data")
plt.plot(noms(x), reg.intercept+reg.slope*noms(x), label="fit")
plt.legend()
plt.grid()
plt.xlabel(r"x $\left[ m \right]$")
plt.ylabel(r"phi $\left[ 1 \right]$")
plt.show()

m = ufloat(reg.slope, reg.stderr)

# deltaPhi/deltaX = 2pi/lambda
wave_len = 2*np.pi*1/m
# v = lambda*f
v_air_meas = wave_len*f_res

# v_theo = sqrt(gamma*R*T/M)
T_lab = ufloat(25, 0)
T_zero = 0 + 273.15
v_zero = 331.7
gammaR_over_M = v_zero**2/T_zero
v_air = sqrt(gammaR_over_M*(T_lab + 273.15))

#%% part 2 - sound velocity in gases
L = ufloat(20e-2, 0)

# helium
# TODO: change M to 4e-3 (4 g/mol)
m = 2 * proton_mass + 2 * neutron_mass
M = m * Avogadro
gamma = 1.66
v_helium = np.sqrt(gamma*R*(T_lab + 273.15)/M)

# CO2
# TODO: change M to 44e-3 (44 g/mol)
m = 22 * proton_mass + 22 * neutron_mass
M = m * Avogadro
gamma = 1.304
v_co2 = np.sqrt(gamma*R*(T_lab + 273.15)/M)

gas = {
    "air": {
        "v_theo": v_air,
        "T_max": None,
        "f_min": None
    },
    "he": {
        "v_theo": v_helium,
        "T_max": None,
        "f_min": None
    },
    "co2": {
        "v_theo": v_co2,
        "T_max": None,
        "f_min": None
    }
}

# T <= 1/2*L/v, f >= 2*v/L
for g in gas.keys():
    gas[g]["T_max"] = L/gas[g]["v"]/2
    gas[g]["f_min"] = 1/gas[g]["T_max"]

gas["air"]["T_pipe"] = ufloat(0, 0)
gas["air"]["T_mic"] = ufloat(0, 0)
gas["air"]["v_meas"] = L/gas["air"]["T_pipe"]

gas["he"]["T_pipe"] = ufloat(0, 0)
gas["he"]["T_mic"] = ufloat(0, 0)
gas["he"]["v_meas"] = L/gas["he"]["T_pipe"]

gas["co2"]["T_pipe"] = ufloat(0, 0)
gas["co2"]["T_mic"] = ufloat(0, 0)
gas["co2"]["v_meas"] = L/gas["co2"]["T_pipe"]

#%% part 3 - resonance
f_modes = unumpy.uarray([], 0)
n = np.array([])

reg2 = linregress(n, noms(f))
plt.errorbar(n, noms(f), yerr=devs(f), fmt="bo", label="data")
plt.plot(n, reg2.intercept + reg2.slope*n, label="fit")
plt.legend()
plt.grid()
plt.xlabel(r"n $\left[ 1 \right]$")
plt.ylabel(r"f $\left[ \frac{1}{sec} \right]$")
# f = v/lambda = nv/2L
m2 = ufloat(reg2.slope, reg2.stderr)
v_air_meas_2 = 2*L*m2

f_n = np.arange(0, 5, 1)*v_air/2/L
plt.errorbar(n, noms(f_n), yerr=devs(f_n), fmt="ro", label="theoretical")
plt.show()

#%% part 4 - standing waves
L2 = ufloat(20e-2, 0)
f_n_2 = np.arange(0, 5, 1)*v_air/2/L2
f_prim = f_n_2[3]

x2 = unumpy.uarray([], 0)
A2 = unumpy.uarray([], 0)

plt.errorbar(noms(x2), noms(A2), xerr=devs(x2), yerr=devs(A2), fmt="bo", label="data")
plt.grid()
plt.xlabel(r"x $\left[ m \right]$")
plt.ylabel(r"A $\left[ m \right]$")
plt.show()



