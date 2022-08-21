import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from uncertainties.umath import sqrt
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.constants import R, Avogadro, proton_mass, neutron_mass

#%% part 0 - prep

f = unumpy.uarray([3,  3.2, 3.5, 3.8, 4.1, 4.4, 4.7, 5, 5.3, 5.6, 5.9, 6.2], 1e-3)  # KHz
A = unumpy.uarray([19, 20,  21.5, 23.5, 26.5, 29, 29.7, 34, 25, 21.7, 20, 19.5 ], [4, 4, 4, 4, 2, 1, 1, 1, 1, 1, 2, 3 ])  # mV

plt.errorbar(noms(f), noms(A), xerr=devs(f), yerr=devs(A), fmt="bo", label="data")
plt.grid()
plt.xlabel(r"f $\left[ KHz \right]$")
plt.ylabel(r"A $\left[ mV \right]$")
plt.show()

res_idx = np.argmax(noms(A))
f_res, A_res = f[res_idx], A[res_idx]

f_res = ufloat(5000, 1)

#%% part 1 - sound phase velocity in air

x0 = 127e-2
x = unumpy.uarray([127, 130.5, 134, 137.4, 140.7, 144.5, 147.9, 151.6, 155], 0.3)*1e-2
x = x - x0
phi = unumpy.uarray([0, 177, 360+0, 360+177, 720+0, 720+177, 1080+0, 1080+177, 360*4+0], [1.5, 0.5, 3, 5, 3, 2,  3, 5, 5])
phi = unumpy.uarray(np.deg2rad(noms(phi)), np.deg2rad(devs(phi)))
reg = linregress(noms(x), noms(phi))

plt.errorbar(noms(x), noms(phi), xerr=devs(x), yerr=devs(phi), fmt="bo", label="data")
plt.plot(noms(x), reg.intercept+reg.slope*noms(x), label="fit")
plt.legend()
plt.grid()
plt.xlabel(r"$\Delta$x $\left[ m \right]$")
plt.ylabel(r"$\Delta\phi$ $\left[ radians \right]$")
plt.show()

m = ufloat(reg.slope, reg.stderr)

# deltaPhi/deltaX = 2pi/lambda
wave_len = 2*np.pi/m
# v = lambda*f
v_air_meas = wave_len*f_res

# v_theo = sqrt(gamma*R*T/M)
T_lab = ufloat(21.3, 0.3)
T_zero = 0 + 273.15
v_zero = 331.7
gammaR_over_M = v_zero**2/T_zero
v_air = sqrt(gammaR_over_M*(T_lab + 273.15))

#%% part 2 - sound velocity in gases
L = ufloat(0.98, 0.01)

# helium
M = 4e-3
gamma = 1.66
v_helium = sqrt(gamma*R*(T_lab + 273.15)/M)

# CO2
M = 44e-3
gamma = 1.304
v_co2 = sqrt(gamma*R*(T_lab + 273.15)/M)

gas = {
    "air": {
        "v_theo": v_air,
        "T_min": None,
        "f_max": None
    },
    "he": {
        "v_theo": v_helium,
        "T_min": None,
        "f_max": None
    },
    "co2": {
        "v_theo": v_co2,
        "T_min": None,
        "f_max": None
    }
}

# T_pipe - edge up of 4th peak
gas["air"]["T_pipe"] = ufloat(19.90e-3, 0.02e-3)
gas["he"]["T_pipe"] = ufloat(8.58e-3, 0.02e-3)
gas["co2"]["T_pipe"] = ufloat(25.04e-3, 0.02e-3)

# T <= 1/2*L/v, f >= 2*v/L
for g in gas.keys():
    gas[g]["T_min"] = L/gas[g]["v_theo"]*4
    gas[g]["f_max"] = 1/gas[g]["T_min"]
    gas[g]["v_meas"] = 7 * L / gas[g]["T_pipe"]
    gas[g]["relative_err"] = gas[g]["v_meas"].std_dev / gas[g]["v_meas"].nominal_value

gas["he"]["air_prop"] = (7*L/gas["he"]["v_theo"] - gas["he"]["T_pipe"]) / (7*L/gas["he"]["v_theo"] - 7*L/gas["air"]["v_theo"])
gas["co2"]["air_prop"] = (7*L/gas["co2"]["v_theo"] - gas["co2"]["T_pipe"]) / (7*L/gas["co2"]["v_theo"] - 7*L/gas["air"]["v_theo"])


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

f_theory = n*v_air/2/L
plt.errorbar(n, noms(f_theory), yerr=devs(f_theory), fmt="ro", label="theoretical velocity")
plt.show()

#%% part 4 - standing waves
L2 = ufloat(20e-2, 0)
f2 = ufloat(0, 0)

x2 = unumpy.uarray([], 0)
A2 = unumpy.uarray([], 0)

plt.errorbar(noms(x2), noms(A2), xerr=devs(x2), yerr=devs(A2), fmt="bo", label="data")
plt.grid()
plt.xlabel(r"x $\left[ m \right]$")
plt.ylabel(r"A $\left[ mV \right]$")
plt.show()

standing_wavelen = ufloat(0, 0)
standing_a = ufloat(0, 0)
standing_phi = ufloat(0, 0)

x = np.linspace(0, x2[len(x2)-1].nominal_value, 1000)

def theo_curve(x, lamb, a, phi0):
    return 2*a*np.cos(2*np.pi*(L2/(2*lamb))-phi0/2)*np.cos(np.pi/lamb*(2*x-L2)+phi0/2)

plt.errorbar(noms(x2), noms(A2), xerr=devs(x2), yerr=devs(A2), fmt="bo", label="data")
plt.plot(x, theo_curve(x, standing_wavelen, standing_a, standing_phi), label="theoretical")
plt.legend()
plt.grid()
plt.xlabel(r"x $\left[ m \right]$")
plt.ylabel(r"A $\left[ mV \right]$")
plt.show()

v_air_meas_3 = standing_wavelen*f2

x3 = unumpy.uarray([], 0)
A3 = unumpy.uarray([], 0)

plt.errorbar(noms(x3), noms(A3), xerr=devs(x3), yerr=devs(A3), fmt="bo", label="data")
plt.grid()
plt.xlabel(r"x $\left[ m \right]$")
plt.ylabel(r"A $\left[ mV \right]$")
plt.show()


