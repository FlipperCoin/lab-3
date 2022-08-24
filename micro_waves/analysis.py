import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from uncertainties.umath import sqrt
import matplotlib.pyplot as plt
from scipy.stats import linregress

#%% grid polarizer

theta_1 = unumpy.uarray([10], [2])
theta_1 = unumpy.uarray(np.deg2rad(noms(theta_1)), np.deg2rad(devs(theta_1)))
I_1 = unumpy.uarray([5], [1])

x1 = unumpy.cos(theta_1)**2
reg = linregress(noms(x1), noms(I_1))
plt.errorbar(noms(x1), noms(I_1), xerr=devs(x1), yerr=devs(I_1), fmt='bo', label='data')
plt.plot(noms(x1), reg.intercept + reg.slope*noms(x1), label='fit')
plt.grid()
plt.legend()
plt.show()

theta_2 = unumpy.uarray([10], [2])
theta_2 = unumpy.uarray(np.deg2rad(noms(theta_2)), np.deg2rad(devs(theta_2)))
I_2 = unumpy.uarray([5], [1])

x2 = unumpy.cos(theta_2)**4
reg = linregress(noms(x2), noms(I_2))
plt.errorbar(noms(x2), noms(I_2), xerr=devs(x2), yerr=devs(I_2), fmt='bo', label='data')
plt.plot(noms(x2), reg.intercept + reg.slope*noms(x2), label='fit')
plt.grid()
plt.legend()
plt.show()

#%% waveguide

# perpendicular
d = unumpy.uarray([], [])
I_3 = unumpy.uarray([], [])
plt.errorbar(noms(d), noms(I_3), xerr=devs(d), yerr=devs(I_3), fmt='bo', label='data')
plt.grid()
plt.show()

theta_3 = unumpy.uarray([10], [2])
theta_3 = unumpy.uarray(np.deg2rad(noms(theta_3)), np.deg2rad(devs(theta_3)))
I_4 = unumpy.uarray([5], [1])

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta_3, I_4, 'bo')
ax.grid(True)
plt.show()

# parallel
d2 = unumpy.uarray([], [])
I_5 = unumpy.uarray([], [])
plt.errorbar(noms(d2), noms(I_5), xerr=devs(d2), yerr=devs(I_5), fmt='bo', label='data')
plt.grid()
plt.show()

d_min = ufloat(2, 2)
lamb = 2*d_min

d_lab = ufloat(2, 2)
tmp = 1/(lamb**2)-1/(2*d_lab)**2
lamb_g_exp = sqrt(1/tmp)

node_deltas = unumpy.uarray([], [])
antinode_deltas = unumpy.uarray([], [])

lamb_g_2 = unumpy.uarray([], [])
d3 = unumpy.uarray([], [])
y = 1/(lamb_g_2**2)
x = (1/(2*d3))**2
reg2 = linregress(x, y)
plt.errorbar(noms(x), noms(y), xerr=devs(x), yerr=devs(y), fmt='bo', label='data')
lamb_2 = ufloat(reg2.intercept, reg2.intercept_stderr)

L = 15e-2

def calc_d(delta_phi):
    # tmp: 1/(2d)**2
    tmp = 1/(lamb_2**2)-(delta_phi/(2*np.pi*L)-1/lamb_2)**2
    d = 1/2*np.sqrt(1/tmp)
    return d

d_circ = calc_d(np.pi/2)
d_lin = calc_d(2*np.pi)
d_lin2 = calc_d(np.pi)
d_lin3 = calc_d(4*np.pi)

theta_circ = unumpy.uarray(range(0, 360, 15), 1)
theta_circ = unumpy.uarray(np.deg2rad(noms(theta_circ)), np.deg2rad(devs(theta_circ)))
I_circ = unumpy.uarray([5], [1])

theta_lin = unumpy.uarray(range(0, 360, 15), 1)
theta_lin = unumpy.uarray(np.deg2rad(noms(theta_lin)), np.deg2rad(devs(theta_lin)))
I_lin = unumpy.uarray([5], [1])

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta_circ, I_circ, 'bo')
ax.grid(True)
plt.show()

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta_lin, I_lin, 'bo')
ax.grid(True)
plt.show()


