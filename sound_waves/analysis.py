import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from uncertainties.umath import sqrt
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.constants import R, Avogadro, proton_mass, neutron_mass

#%% part 0 - prep

f = unumpy.uarray([3,  3.2, 3.5, 3.8, 4.1, 4.4, 4.7, 5, 5.3, 5.6, 5.9, 6.2], 0.1e-3)  # KHz
A = unumpy.uarray([19, 20,  21.5,23.5, 26.5, 29, 29.7, 34, 25, 21.7, 20, 19.5 ], [4, 4, 4, 4, 2, 1, 1, 1, 1, 1, 2, 3 ])  # mV

plt.errorbar(noms(f), noms(A), xerr=devs(f), yerr=devs(A), fmt="bo", label="data")
plt.grid()
plt.xlabel(r"f $\left[ \frac{1}{sec} \right]$")
plt.ylabel(r"A $\left[ m \right]$")
plt.show()

res_idx = np.argmax(noms(A))
f_res, A_res = f[res_idx], A[res_idx]

f_res = 5000

#%% part 1 - sound phase velocity in air
#
# x0 = 127
# x = unumpy.uarray([127, 130.5,134, 137.4, 140.7,144.5,147.9,151.6, 155  ], 0.3)
# phi = unumpy.uarray([0, 177, 0,   177   , 0    ,177  ,0    ,177  ,0  ], [1.5, 0.5, 3, 5, 3, 2,3,5])
# reg = linregress(noms(x), noms(phi))
#
# plt.errorbar(noms(x), noms(phi), xerr=devs(x), yerr=devs(phi), fmt="bo", label="data")
# plt.plot(noms(x), reg.intercept+reg.slope*noms(x), label="fit")
# plt.legend()
# plt.grid()
# plt.xlabel(r"x $\left[ m \right]$")
# plt.ylabel(r"phi $\left[ 1 \right]$")
# plt.show()
#
# m = ufloat(reg.slope, reg.stderr)
#
# # deltaPhi/deltaX = 2pi/lambda
# wave_len = 2*np.pi*1/m
# # v = lambda*f
# v_air_meas = wave_len*f_res

# v_theo = sqrt(gamma*R*T/M)
T_lab = ufloat(21.3, 0.3)
T_zero = 0 + 273.15
v_zero = 331.7
gammaR_over_M = v_zero**2/T_zero
v_air = sqrt(gammaR_over_M*(T_lab + 273.15))

#%% part 2 - sound velocity in gases
L = ufloat(0.98, 0.001)

# helium
m = 4e-3
M = m
gamma = 1.66
v_helium = sqrt(gamma*R*(T_lab + 273.15)/M)

# CO2
m = 44e-3
M = m
gamma = 1.304
v_co2 = sqrt(gamma*R*(T_lab + 273.15)/M)

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
    gas[g]["T_min"] = L/gas[g]["v_theo"]*4
    gas[g]["f_max"] = 1/gas[g]["T_min"]

# edge up of 4th peak
T_air = ufloat(19.90e-3, 0.02e-3)
T_co2 = ufloat(25.04e-3, 0.02e-3)
T_he = ufloat(8.58e-3, 0.02e-3)

gas["air"]["T_pipe"] = ufloat(0, 0)
gas["air"]["T_mic"] = ufloat(0, 0)
gas["air"]["v_meas"] = L/gas["air"]["T_pipe"]

gas["he"]["T_pipe"] = ufloat(0, 0)
gas["he"]["T_mic"] = ufloat(0, 0)
gas["he"]["v_meas"] = L/gas["he"]["T_pipe"]

gas["co2"]["T_pipe"] = ufloat(0, 0)
gas["co2"]["T_mic"] = ufloat(0, 0)
gas["co2"]["v_meas"] = L/gas["co2"]["T_pipe"]
#
# #%% part 3 - resonance
# f_modes = unumpy.uarray([], 0)
# n = np.array([])
#
# reg2 = linregress(n, noms(f))
# plt.errorbar(n, noms(f), yerr=devs(f), fmt="bo", label="data")
# plt.plot(n, reg2.intercept + reg2.slope*n, label="fit")
# plt.legend()
# plt.grid()
# plt.xlabel(r"n $\left[ 1 \right]$")
# plt.ylabel(r"f $\left[ \frac{1}{sec} \right]$")
# # f = v/lambda = nv/2L
# m2 = ufloat(reg2.slope, reg2.stderr)
# v_air_meas_2 = 2*L*m2
#
# f_n = np.arange(0, 5, 1)*v_air/2/L
# plt.errorbar(n, noms(f_n), yerr=devs(f_n), fmt="ro", label="theoretical")
# plt.show()
#
# #%% part 4 - standing waves
# L2 = ufloat(20e-2, 0)
# f_n_2 = np.arange(0, 5, 1)*v_air/2/L2
# f_prim = f_n_2[3]
#
# x2 = unumpy.uarray([], 0)
# A2 = unumpy.uarray([], 0)
#
# plt.errorbar(noms(x2), noms(A2), xerr=devs(x2), yerr=devs(A2), fmt="bo", label="data")
# plt.grid()
# plt.xlabel(r"x $\left[ m \right]$")
# plt.ylabel(r"A $\left[ m \right]$")
# plt.show()
#
#
#
