import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from uncertainties.umath import sqrt
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.constants import R, Avogadro, proton_mass, neutron_mass
from scipy.optimize import curve_fit as cfit

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
phi = unumpy.uarray([0, 177, 360+0, 360+177, 720+0, 720+177, 1080+0, 1080+177, 360*4+0], [1.5, 0.5, 3, 5, 3, 2, 3, 5, 5])
phi = unumpy.uarray(np.deg2rad(noms(phi)), np.deg2rad(devs(phi)))
reg = linregress(noms(x), noms(phi))

plt.errorbar(noms(x), noms(phi), xerr=devs(x), yerr=devs(phi), fmt="bo", label="data")
plt.plot(noms(x), reg.intercept+reg.slope*noms(x), label="fit")
plt.legend()
plt.grid()
plt.xlabel(r"$\Delta$x $\left[ m \right]$")
plt.ylabel(r"$\Delta\phi$ $\left[ r \right]$")
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

T_lab2 = ufloat(21.9, 0.1)

f_modes = unumpy.uarray([178.3, 354.7, 530.7, 709.8, 884.4, 1061, 1236, 1415, 1592, 1768, 1942, 2121, 2299, 2475, 2652], [0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
n = np.arange(1, len(f_modes)+1)

reg2 = linregress(n, noms(f_modes))
plt.errorbar(n, noms(f_modes), yerr=devs(f_modes), fmt="bo", label="data")
plt.plot(n, reg2.intercept + reg2.slope*n, label="fit")
plt.grid()
plt.xlabel(r"n $\left[ 1 \right]$")
plt.ylabel(r"f $\left[ \frac{1}{sec} \right]$")
# f = v/lambda = nv/2L
m2 = ufloat(reg2.slope, reg2.stderr)
v_air_meas_2 = 2*L*m2

f_theory = n*v_air/2/L
plt.errorbar(n, noms(f_theory), yerr=devs(f_theory), fmt="ro", label="theoretical velocity")
plt.legend()
plt.title("normal modes frequency vs mode number n")
plt.savefig("1.png")
plt.show()

#%% part 4 - standing waves
L2 = ufloat(246e-3, 1e-3)
f2 = ufloat(4900, 1)

lamb_n = 2/np.arange(1, 15)*L2
f_n = gas["air"]["v_theo"]/lamb_n

chosen_frequency = 4900 #Hz

x2 = unumpy.uarray(np.arange(3, 25.5, 0.5), 0.1)*1e-2  # m
A2 = unumpy.uarray([35, 40, 32, 24, 11, 6, 20, 33, 40, 37, 23, 10, 3, 20, 30, 39, 36, 26, 14, 9, 16, 30, 37, 35, 28, 16, 12, 17, 27, 36, 37, 30, 22, 13, 6, 21, 33, 35, 33, 26, 19, 20, 27, 34, 32],
                   [2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 2, 1, 4, 1, 2, 4, 2, 3, 1, 2, 2, 1, 1, 3, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3])  # mV

x2 = x2[0:int(np.floor(len(x2)/2))]
A2 = A2[0:int(np.floor(len(A2)/2))]

# plt.errorbar(noms(x2), noms(A2), xerr=devs(x2), yerr=devs(A2), fmt="bo", label="data")
# plt.grid()
# plt.xlabel(r"x $\left[ m \right]$")
# plt.ylabel(r"A $\left[ mV \right]$")
# plt.show()

def theo_curve(x, lamb, a, phi0):
    return 2*a*np.abs(np.cos(np.pi/lamb*(2*x-L2.nominal_value)+phi0/2))

fit = cfit(theo_curve, noms(x2), noms(A2), p0=[0.08, 20, 0])

standing_wavelen = ufloat(fit[0][0], np.sqrt(np.diag(fit[1]))[0])
standing_a = ufloat(fit[0][1], np.sqrt(np.diag(fit[1]))[1])
standing_phi = ufloat(fit[0][2], np.sqrt(np.diag(fit[1]))[2])

x = np.linspace(0.03, x2[len(x2)-1].nominal_value, 1000)

plt.errorbar(noms(x2), noms(A2), xerr=devs(x2), yerr=devs(A2), fmt="bo", label="data")
plt.plot(x, noms(theo_curve(x, standing_wavelen.nominal_value, standing_a.nominal_value, standing_phi.nominal_value)), label="fit")
plt.legend()
plt.grid()
plt.xlabel(r"x $\left[ m \right]$")
plt.ylabel(r"A $\left[ mV \right]$")
plt.title("2 mic, standing wave")
plt.savefig("2.png")
plt.show()

v_air_meas_3 = standing_wavelen*f2

x3 = unumpy.uarray(np.arange(3, 25.5, 0.5), 0.1)*1e-2  # m
A3 = unumpy.uarray([13, 39, 55, 60, 57, 49, 35, 8, 30.5, 57, 62.9, 60, 51.4, 34,  11, 28.5, 54, 59, 61, 53, 43, 23, 21, 53, 57, 55, 53, 39, 17, 15, 45, 58, 61, 55, 44, 23, 11, 41, 54, 58, 54, 46, 30, 8.5, 31],
                   [1,  1,  2,  1,  2,  1,  1,  1, 1,    1,  0.2,  0.3,0.3,  0.5, 1,  1,    1,  2,  1,  2,  1,  1,  1,  1,  1,  1,  1,  2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0.5, 0.5])

plt.errorbar(noms(x3), noms(A3), xerr=devs(x3), yerr=devs(A3), fmt="bo", label="data")
plt.grid()
plt.xlabel(r"x $\left[ m \right]$")
plt.ylabel(r"A $\left[ mV \right]$")
plt.title("1 mic, standing wave")
plt.savefig("3.png")
plt.show()


