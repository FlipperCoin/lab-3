import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from uncertainties.umath import sqrt
import matplotlib.pyplot as plt
from scipy.stats import linregress

#%% michaleson
extrema_locations = unumpy.uarray([],)
I = unumpy.uarray([],)
plt.errorbar(noms(extrema_locations), noms(I), xerr=devs(extrema_locations), yerr=devs(I), fmt='bo')
wavelength = ufloat(,)
plt.xlabel('locations of exterma in intensity [m]')
plt.ylabel('Intensity [V]')
plt.grid()
plt.legend()
plt.show()

#%% Fabry - Perot
d1 = ufloat(,)
minima_locations = unumpy.uarray([],)
minima_count = ufloat(,)
d2 = ufloat(,)
wavelength = 2 * d1 #order 1?

d1_2 = ufloat(,)
minima_locations_2 = unumpy.uarray([],)
minima_count_2 = ufloat(,)
d2_2= ufloat(,)
wavelength_2 = 2 * d1 #order 1?

#%% Lloyd
d1 = ufloat(,)
h1 = ufloat(,)
distance_1 = unumpy.uarray([],)
I_min_1 = unumpy.uarray([],)
h2 = ufloat(,)
distance_2 = unumpy.uarray([],)
I_max_1 = unumpy.uarray([],)
wavelength_1 = 2*d1 #not sure

d1_2= ufloat(,)
h1_2 = ufloat(,)
distance_1_2 = unumpy.uarray([],)
I_min_2 = unumpy.uarray([],)
h2_2= ufloat(,)
distance_2_2 = unumpy.uarray([],)
I_max_2 = unumpy.uarray([],)
wavelength_2 = 2*d1_2 #not sure


