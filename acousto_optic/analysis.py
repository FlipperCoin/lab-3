import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import uarray
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from uncertainties.umath import sqrt, tan
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.fft import ifft, fft
from scipy import misc
import imageio
from scipy.constants import k

lambd = ufloat(632.8, 0.1) * 1e-9

T_lab = ufloat(25.3, 0.1)
T_20 = 293.204
v_20 = 1160.02
gammaR_over_M = v_20**2/T_20
v_ethanol = sqrt(gammaR_over_M*(T_lab + 273.15))

#%% near-field

nf_f = uarray([2.1667300,    1.8914300, 1.8004300, 1.7001300, 2.6792300, 2.8692300, 2.2283300, 2.4196500, 3.7596400, 1.5792400], 0.00002) * 1e6
sort = np.argsort(noms(nf_f))
nf_f = np.take(nf_f, sort)
# first 3 might need diff denu
nf_delta_x = uarray([575/11, 526/9,     561/9,     802/12,    510/12,    514/13,    553/11,    519/11,     418/14,    498/7], 2) * 5.2e-6
nf_delta_x = np.take(nf_delta_x, sort)
# dx = lambda/2 = v/(2f)
def fit_func(f, m):
    return m/f

popt, pcov = curve_fit(fit_func, noms(nf_f), noms(nf_delta_x), sigma=devs(nf_delta_x))
m = ufloat(popt[0], np.sqrt(pcov[0, 0]))
nf_v = 2*m

plt.errorbar(noms(nf_f), noms(nf_delta_x), xerr=devs(nf_f), yerr=devs(nf_delta_x), fmt='bo', label='data')
plt.plot(noms(nf_f), fit_func(noms(nf_f), m.nominal_value), label='fit')
plt.xlabel(r'f $\left[ Hz \right]$')
plt.ylabel(r'$\overline{\Delta x}$  $\left[ m \right]$')
plt.ticklabel_format(axis='y', scilimits=[-4, -4])
plt.legend()
plt.grid()
plt.show()

#%% far-field

f3 = ufloat(300, 1) * 1e-3

# should be same values as in NF, but only in case they give a clear picture
ff_f = uarray([1.5792400,   1.7001300, 1.8004300, 1.8914300, 2.1667300, 2.2283300, 2.4196500, 2.6792300, 2.8692300, 3.7596400], 0.00002) *1e6
ff_delta_x = uarray([110/2, 119/2,     127/2,     136/2,     151/2,     149/2,     167/2,     184/2,     198/2,     262/2 ], 6) * 5.2e-6
ff_delta_x_max = ff_delta_x + np.linspace(-3, 3, len(ff_delta_x)) * 5.2e-6
ff_delta_x_min = ff_delta_x - np.linspace(-3, 3, len(ff_delta_x)) * 5.2e-6

# uarray([104/2, 114/2,     125/2,     130/2,     147/2,     149/2,     165/2,     181/2,     198/2,     262/2 ], 2) * 5.2e-6
# ff_f = ff_f[:-1]
# ff_delta_x = ff_delta_x[:-1]

# dx = f3*lambda/v * f
# reg = linregress(noms(ff_f), noms(ff_delta_x))
# ff_m = ufloat(reg.slope, reg.stderr)

def lin_fit(x, inter, slope):
    return slope*x+inter

popt, pcov = curve_fit(lin_fit, noms(ff_f), noms(ff_delta_x), sigma=devs(ff_delta_x))
ff_m = ufloat(popt[1], np.sqrt(pcov[1][1]))

ff_v = (f3*lambd)/ff_m

popt_max, pcov_max = curve_fit(lin_fit, noms(ff_f), noms(ff_delta_x_max), sigma=devs(ff_delta_x_max))
ff_m_max = ufloat(popt_max[1], np.sqrt(pcov_max[1][1]))

ff_v_max = (f3*lambd)/ff_m_max

popt_min, pcov_min = curve_fit(lin_fit, noms(ff_f), noms(ff_delta_x_min), sigma=devs(ff_delta_x_min))
ff_m_min = ufloat(popt_min[1], np.sqrt(pcov_min[1][1]))

ff_v_min = (f3*lambd)/ff_m_min

plt.errorbar(noms(ff_f), noms(ff_delta_x), xerr=devs(ff_f), yerr=devs(ff_delta_x), fmt='bo', label='data')
plt.plot(noms(ff_f), popt[0] + popt[1]*noms(ff_f), label='fit')
plt.xlabel(r'f $\left[ Hz \right]$')
plt.ylabel(r'$\overline{\Delta x}$  $\left[ m \right]$')
plt.ticklabel_format(axis='y', scilimits=[-4, -4])
plt.legend()
plt.grid()
plt.show()

plt.errorbar(noms(ff_f), noms(ff_delta_x), xerr=devs(ff_f), yerr=devs(ff_delta_x), fmt='bo', label='data')
plt.plot(noms(ff_f), popt[0] + popt[1]*noms(ff_f), label='fit')
# plt.plot(noms(ff_f), popt_max[0] + popt_max[1]*noms(ff_f), label='fit_max')
plt.plot(noms(ff_f), popt_min[0] + popt_min[1]*noms(ff_f), label='fit_min')
plt.xlabel(r'f $\left[ Hz \right]$')
plt.ylabel(r'$\overline{\Delta x}$  $\left[ m \right]$')
plt.ticklabel_format(axis='y', scilimits=[-4, -4])
plt.legend()
plt.grid()
plt.show()
#
# i = imageio.imread("tst.tif")
# plt.imshow(i)
# plt.show()
#
# inverse = np.abs(ifft(i))
# plt.imshow(inverse.astype('uint8'))
# plt.show()
#
# i2 = imageio.imread("tst2.tif")
# plt.imshow(i2)
# plt.show()
#
# inverse2 = np.abs(fft(i2))
# plt.imshow(inverse2.astype('uint8'))
# plt.show()
