import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import uarray
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from uncertainties.umath import sqrt, tan
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.constants import k

lambd = ufloat(632.8, 0.1) * 1e-9

T_lab = ufloat(25.3, 0.1)
T_20 = 20 + 273.15
v_20 = 1480
gammaR_over_M = v_20**2/T_20
v_water = sqrt(gammaR_over_M*(T_lab + 273.15))

#%% near-field

nf_f = uarray([2.1667300,    1.8914300, 1.8004300, 1.7001300, 2.6792300, 2.8692300, 2.2283300, 2.4196500, 3.7596400, 1.5792400], 0.0000001) * 1e6
# first 3 might need diff denu
nf_delta_x = uarray([575/11, 526/9,     561/9,     802/12,    510/12,    514/13,    553/11,    519/11,     418/14,    498/7], 1) * 5.2e-6
# dx = lambda/2 = v/(2f)
def fit_func(f, m):
    return m/f

popt, pcov = curve_fit(fit_func, noms(nf_f), noms(nf_delta_x), sigma=devs(nf_delta_x))
m = ufloat(popt[0], np.sqrt(pcov[0, 0]))
nf_v = 2*m
# TODO: v theoretical value in literature
# v_theo = ufloat()

plt.errorbar(noms(nf_f), noms(nf_delta_x), xerr=devs(nf_f), yerr=devs(nf_delta_x), fmt='bo', label='data')
plt.plot(noms(nf_f), fit_func(noms(nf_f), m.nominal_value), label='fit')
plt.xlabel(r'f $\left[ Hz \right]$')
plt.ylabel(r'$\overline{\Delta x}$ $\left[ m \right]$')
plt.legend()
plt.grid()
plt.show()

#%% far-field

f3 = ufloat(300, 1) * 1e-3

# should be same values as in NF, but only in case they give a clear picture
ff_f = uarray([1.5792400,   1.7001300, 1.8004300, 1.8914300, 2.1667300, 2.2283300, 2.4196500, 2.6792300, 2.8692300, 3.7596400], 0.0000001) *1e6
ff_delta_x = uarray([110/2, 119/2,     127/2,     136/2,     151/2,     149/2,     167/2,     184/2,     198/2,     262/2 ], 1) * 5.2e-6

# dx = f3*lambda/v * f
reg = linregress(noms(ff_f), noms(ff_delta_x))
ff_m = ufloat(reg.slope, reg.stderr)
ff_v = (f3*lambd)/ff_m

plt.errorbar(noms(ff_f), noms(ff_delta_x), xerr=devs(ff_f), yerr=devs(ff_delta_x), fmt='bo', label='data')
plt.plot(noms(ff_f), reg.intercept + reg.slope*noms(ff_f), label='fit')
plt.xlabel(r'f $\left[ Hz \right]$')
plt.ylabel(r'$\overline{\Delta x}$ $\left[ m \right]$')
plt.legend()
plt.grid()
plt.show()


