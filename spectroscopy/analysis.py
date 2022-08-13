import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from uncertainties import unumpy, ufloat
from uncertainties.umath import sqrt
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from scipy.constants import c, h, electron_volt

from data import *

#%% diffraction grating constant calculation

# d*sin(theta)=n*lambda
# 2d*sin(theta_min/2)=n*lambda

# same wave length, different order numbers
merc_theta_min = unumpy.uarray([-64, -36, -17, 0, 17, 36, 64], 5)
merc_theta_min = 2 * np.pi * merc_theta_min / 360
n = unumpy.uarray(range(-3, 4, 1), 0)

x = unumpy.sin(merc_theta_min / 2)
y = n
reg = linregress(noms(x), noms(y))

plt.errorbar(noms(x), noms(y), xerr=devs(x), yerr=devs(y), fmt="ro", label="data")
plt.plot(noms(x), reg.slope * noms(x) + reg.intercept, label="fit")
plt.xlabel(r"$\sin(\frac{\theta}{2})$ $[1]$")
plt.ylabel("order number n $[1]$")
plt.grid()
plt.legend()
plt.show()

d = ufloat(reg.slope, reg.stderr) * mercury_wavlen / 2

#%% Helium spectrum with diffraction grating

# different wavelength, order number +/-1
hel_theta_min = unumpy.uarray([], 5)
merc_theta_min = 2 * np.pi * merc_theta_min / 360

hel_wavelen = 2 * d * unumpy.sin(hel_theta_min / 2)

exp_hel_wavelen = np.array([396.47, 402.62, 412.08, 414.38, 438.79,
                            447.15, 471.31, 492.19, 501.57, 504.77,
                            587.56, 667.82, 706.52])

#%% prism calibration

# different wavelength, order number +/-1
hel_delta_min = unumpy.uarray([], 5)
hel_delta_min = 2 * np.pi * hel_delta_min / 360

x = 1/(hel_wavelen**2)
y = hel_delta_min

# delta_min = B_tag+C_tag/(lambda**2)
reg2 = linregress(x, y)
B_tag = ufloat(reg2.intercept, reg2.intercept_stderr)
C_tag = ufloat(reg2.slope, reg2.stderr)

plt.errorbar(x, y, xerr=devs(x), yerr=devs(y), fmt="ro", label="data")
plt.plot(x, B_tag+C_tag*x, label="fit")
plt.xlabel(r"$\frac{1}{\lambda^2}$ $[\frac{1}{m^2}]$")
plt.ylabel(r"$\delta_min$ $[1]$")
plt.grid()
plt.legend()
plt.show()

#%% hydrogen spectrum with prism

hyd_delta_min = unumpy.uarray([], 5)
hyd_wavelen = sqrt(C_tag/(hyd_delta_min - B_tag))

nu = c/hyd_wavelen
e_deltas = h*nu

e_deltas = np.take_along_axis(e_deltas, np.argsort(noms(e_deltas)))

exp_nf = np.array([3, 4, 5, 6, 7, 8])
Ry = 13.6*electron_volt
exp_e_deltas = h*c/(Ry*(1/4-1/(exp_nf**2)))

nf = np.array([])
calc_ry = h*c/(e_deltas*(1/4-1/(nf**2)))

calc_ry_mean = np.mean(noms(calc_ry))
