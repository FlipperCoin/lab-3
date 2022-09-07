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

#%% near-field

nf_f = uarray([], 1)
nf_delta_x = uarray([], [])

# dx = lambda/2 = v/(2f)
def fit_func(f, m):
    return m/f

popt, pcov = curve_fit(fit_func, noms(nf_f), noms(nf_delta_x), sigma=devs(nf_delta_x))
m = ufloat(popt[0], np.sqrt(pcov[0, 0]))
nf_v = 2*m
# TODO: v theoretical value in literature
# v_theo = ufloat()

plt.errorbar(noms(nf_f), noms(nf_delta_x), xerr=devs(nf_f), yerr=devs(nf_delta_x), fmt='bo', label='data')
plt.plot(noms(nf_f), fit_func(nf_f, m.nominal_value), label='fit')
plt.xlabel(r'f $\left[ Hz \right]$')
plt.ylabel(r'$\overline{\Delta x}$ $\left[ m \right]$')
plt.legend()
plt.grid()
plt.show()

#%% far-field

f3 = ufloat(300, 1) * 1e-3

# should be same values as in NF, but only in case they give a clear picture
ff_f = uarray([], 1)
ff_delta_x = uarray([], [])

# dx = f3*lambda/v * f
reg = linregress(noms(ff_f), noms(ff_delta_x))
ff_m = ufloat(reg.slope, reg.intercept)
ff_v = (f3*lambd)/ff_m

plt.errorbar(noms(ff_f), noms(ff_delta_x), xerr=devs(ff_f), yerr=devs(ff_delta_x), fmt='bo', label='data')
plt.plot(noms(ff_f), reg.intercept + reg.slope*noms(ff_f), label='fit')
plt.xlabel(r'f $\left[ Hz \right]$')
plt.ylabel(r'$\overline{\Delta x}$ $\left[ m \right]$')
plt.legend()
plt.grid()
plt.show()


