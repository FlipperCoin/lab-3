import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from uncertainties import unumpy, ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from scipy.constants import c, h, electron_volt

from data import *

#%% diffraction grating constant calculation

# d*sin(theta)=n*lambda
# 2d*sin(theta_min/2)=n*lambda

# same wave length, different order numbers
merc_theta_min = unumpy.uarray([48.735, 31.95, 15.885, 0, -15.795, -31.815, -48.735], 0.001)
merc_theta_min = 2 * np.pi * merc_theta_min / 360
n = unumpy.uarray([3, 2, 1, 0, -1, -2, -3], 0)

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
# color order: -20.43 light red, -19.305 red, -16.965 orange, -14.49 green, -14.22 cyan/green, -13.635 blue, -12.915 purple, 12.96 purple, 13.59 blue, 14.22 green/cyan, 14.535 green, 17.01 orange, 19.35 red, 20.43 light red
hel_theta_min = unumpy.uarray([-20.43, -19.305, -16.965, -14.49, -14.22, -13.635, -12.915, 12.96, 13.59, 14.22, 14.535, 17.01, 19.35, 20.43], 0.01)
hel_theta_min = 2 * np.pi * hel_theta_min / 360

hel_wavelen = 2 * d * unumpy.sin(hel_theta_min / 2)

exp_hel_wavelen = np.array([396.47, 402.62, 412.08, 414.38, 438.79,
                            447.15, 471.31, 492.19, 501.57, 504.77,
                            587.56, 667.82, 706.52])

#%% prism calibration

# different wavelength, order number +/-1
hel_delta_min = unumpy.uarray([-48.05, -47.565, -47.205, -47.07, -46.125, -45.54, -45.315], 0.001)
hel_delta_min = 2 * np.pi * hel_delta_min / 360
hel_delta_min *= -1

x = 1/(hel_wavelen[7:]**2)
y = hel_delta_min

# delta_min = B_tag+C_tag/(lambda**2)
reg2 = linregress(noms(x), noms(y))
B_tag = ufloat(reg2.intercept, reg2.intercept_stderr)
C_tag = ufloat(reg2.slope, reg2.stderr)

plt.errorbar(noms(x), noms(y), xerr=devs(x), yerr=devs(y), fmt="ro", label="data")
plt.plot(noms(x), B_tag.nominal_value+C_tag.nominal_value*noms(x), label="fit")
plt.xlabel(r"$\frac{1}{\lambda^2}$ $\left[ \frac{1}{m^2} \right]$")
plt.ylabel(r"$\delta_min$ $[1]$")
plt.grid()
plt.legend()
plt.show()

#%% hydrogen spectrum with prism
# red, light blue, purple, deep purple
hyd_delta_min = unumpy.uarray([-45.09, -46.8, -47.88, -48.555], 0.001)
hyd_delta_min = 2 * np.pi * hyd_delta_min / 360
hyd_delta_min *= -1
hyd_wavelen = unumpy.sqrt((C_tag/(hyd_delta_min - B_tag)))

nu = c/hyd_wavelen
e_deltas = h*nu

e_deltas = np.take_along_axis(e_deltas, np.argsort(noms(e_deltas)), axis=0)

exp_nf = np.array([3, 4, 5, 6])
Ry = 13.6*electron_volt
exp_e_deltas = Ry*(1/4-1/(exp_nf**2))

nf = np.array([3, 4, 5, 6])
calc_ry = e_deltas/(1/4-1/(nf**2))

calc_ry_mean = np.mean(noms(calc_ry))
calc_ry_mean_ev = calc_ry_mean / electron_volt
