import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import uarray
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from uncertainties.umath import sqrt, tan
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.constants import k

#%% scattering

# prob no need for code...

#%% reflection

# calc n
theta_B = ufloat(np.deg2rad(57), np.deg2rad(5))
n = tan(theta_B)

#%% diffraction - single slit

L = ufloat(35+18+33*5, 0)
screw_val = ufloat(8, 1)
a_meas = ufloat(27, 1) - screw_val
lambd = ufloat(0, 0)

minimums = uarray([],
                  [])

# sin(theta) = (lambda/a)*n

minimums_sin_theta = minimums/sqrt(minimums**2+L**2)

n = np.arange(len(minimums))+1
reg = linregress(n, noms(minimums_sin_theta))
plt.errorbar(n, minimums_sin_theta, yerr=devs(minimums_sin_theta), fmt='bo', label='data')
plt.plot(n, n*reg.slope+reg.intercept, label='fit')
plt.grid()
plt.legend()
plt.show()

a_calc = lambd / ufloat(reg.slope, reg.stderr)

x = uarray([], [])
a_meas_2 = ufloat(11, 1) - screw_val
a_list = [a_meas, a_meas_2]
for a in a_list:
    X = a * unumpy.sin(unumpy.arctan2(x/L)) / lambd
    I = (unumpy.sin(np.pi*X)**2) / ((np.pi*X)**2)
    plt.plot(x, I)
    plt.show()

#%% diffraction & superposition - 2 slits

minimums = uarray([],
                  [])

# sin(theta) = (lambda/2d)*(2n+1) = (lambda/d)*n+(lambda/2d)

minimums_sin_theta = minimums/sqrt(minimums**2+L**2)

n = np.arange(len(minimums))
reg2 = linregress(n, noms(minimums_sin_theta))
plt.errorbar(n, minimums_sin_theta, yerr=devs(minimums_sin_theta), fmt='bo', label='data')
plt.plot(n, n*reg.slope+reg.intercept, label='fit')
plt.grid()
plt.legend()
plt.show()

d_slope = lambd/ufloat(reg2.slope, reg2.stderr)
d_intercept = ufloat(reg2.intercept, reg2.intercept_stderr)

d = d_slope + d_intercept / 2

#%% diffraction & superposition - diffraction grating

L = ufloat(0, 0)
maximums = uarray([],
                  [])
n = uarray([],
           [])

maximums_sin_theta = maximums/sqrt(maximums**2+L**2)

reg3 = linregress(n, noms(maximums_sin_theta))
plt.errorbar(n, minimums_sin_theta, yerr=devs(minimums_sin_theta), fmt='bo', label='data')
plt.plot(n, n*reg.slope+reg.intercept, label='fit')
plt.grid()
plt.legend()
plt.show()

d = lambd/ufloat(reg3.slope, reg3.stderr)
density = 1e-3/d
density_theo = ufloat(500, 50)

#%% interferometer

lambd0 = ufloat(0, 0)

p_lab = ufloat(0, 0)
T_lab = ufloat(0 + 273.15, 0)

gas = {'air': {}, 'He': {}, 'CO2': {}, 'He_CO2': {}}

perr = 0

gas['air']['p_min'] = ufloat(0, 0)
gas['air']['F'] = 10
gas['air']['p'] = uarray([], perr)
gas['air']['n_literature'] = ufloat(1.0002926, 1e-7)

gas['He']['p_min'] = ufloat(0, 0)
gas['He']['F'] = 10
gas['He']['p'] = uarray([], perr)
gas['He']['n_literature'] = ufloat(1.000036, 1e-6)

gas['CO2']['p_min'] = ufloat(0, 0)
gas['CO2']['F'] = 10
gas['CO2']['p'] = uarray([], perr)
gas['CO2']['n_literature'] = ufloat(1.000448, 6e-6)

gas['He_CO2']['p_min'] = ufloat(0, 0)
gas['He_CO2']['F'] = 10
gas['He_CO2']['p'] = uarray([], perr)

for g in gas.keys():
    Fn = np.arange(gas[g]['F'])+1
    reg4 = linregress(Fn, gas[g]['p'])
    gas[g]['alpha'] = (2*k*T_lab*lambd0/L) / ufloat(reg4.slope, reg4.stderr)

    # TODO: check about 76cmHg
    gas[g]['n_T0_p76'] = 1 + gas[g]['alpha']*76/(2*k*273.15*lambd0)

# N_He / N_CO2
gas['He_CO2']['ratio'] = (gas['He']['n_literature'] - gas['He_CO2']['n_T0_p76']) / \
                         (gas['He_CO2']['n_T0_p76'] - gas['CO2']['n_literature'])

