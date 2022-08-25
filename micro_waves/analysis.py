import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from uncertainties.umath import sqrt
import matplotlib.pyplot as plt
from scipy.stats import linregress

#%% grid polarizer

# deltaLines = 4.75 degrees, dev = 0.25
angle_delta = ufloat(360/72,0.25)
# 17 lines per 1/4
theta_1 = unumpy.uarray([0, 5, 10, 15, 18, 23, 28, 33, 36, 41, 46, 51, 54, 59, 64, 69, 72 ], 1) * angle_delta
theta_1 = unumpy.uarray(np.deg2rad(noms(theta_1)), np.deg2rad(devs(theta_1)))
I_1 = unumpy.uarray([0, 0.076, 0.211, 0.272, 0.283, 0.263, 0.183 , 0.039, 0.0001, 0.098, 0.219, 0.275, 0.0282, 0.253, 0.171, 0.036, 0],
                    [0.001, 0.002, 0.002, 0.001, 0.001, 0.001, 0.002, 0.001, 0.0001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.0001 ])

x1 = unumpy.cos(theta_1)**2
reg = linregress(noms(x1), noms(I_1))
plt.errorbar(noms(x1), noms(I_1), xerr=devs(x1), yerr=devs(I_1), fmt='bo', label='data')
plt.plot(noms(x1), reg.intercept + reg.slope*noms(x1), label='fit')
plt.xlabel('cos^2(theta)')
plt.ylabel('I')
plt.grid()
plt.legend()
plt.show()

theta_2 = unumpy.uarray(list(range(0, 105, 15)), 2)
theta_2 = unumpy.uarray(np.deg2rad(noms(theta_2)), np.deg2rad(devs(theta_2)))
I_2 = unumpy.uarray([0.021, 0.009, 0.03, 0.135,  0.214, 0.253, 0.261],
                    [0.001, 0.001, 0.002, 0.001, 0.001, 0.001, 0.001])

x2 = unumpy.cos(theta_2)**4
reg = linregress(noms(x2), noms(I_2))
plt.errorbar(noms(x2), noms(I_2), xerr=devs(x2), yerr=devs(I_2), fmt='bo', label='data')
plt.plot(noms(x2), reg.intercept + reg.slope*noms(x2), label='fit')
plt.grid()
plt.legend()
plt.show()

# #%% waveguide

# perpendicular
d1 = unumpy.uarray([1.5, 1.8,    2,      2.5  , 3,      3.3,       3.6,    4,    4.2,  2.9], 0.3) * 10**(-2)
I_3 = unumpy.uarray([0.25,0.243  ,0.274, 0.235, 0.286,  0.362,    0.375,  0.376,  0.39   ,0.252], 0.002)
plt.errorbar(noms(d1), noms(I_3), xerr=devs(d1), yerr=devs(I_3), fmt='bo', label='data')
plt.grid()
plt.show()
#
theta_3 = unumpy.uarray([0, 5, 10, 15, 20, 25, 30,                 35,  40,       45,      50,      55], 1) * angle_delta
theta_3 = unumpy.uarray(np.deg2rad(noms(theta_3)), np.deg2rad(devs(theta_3)))
I_4 = unumpy.uarray([0.003,0.143, 0.33, 0.357, 0.379, 0.356, 0.172, 0.004, 0.18,  0.335,   0.365 ,  0.375], 0.002)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(noms(theta_3), noms(I_4), 'bo')
ax.grid(True)
plt.show()

# # parallel
d2 = unumpy.uarray([1.7,2, 2.5, 2.7, 3, 3.5, 4, 4.5, 5, 3.2 ], 0.3) * 10**(-2)

I_5 = unumpy.uarray([0.265,0.326, 0.39, 0.367, 0.403, 0.414, 0.441, 0.443, 0.428, 0.405 ], 0.001)

I_5 = unumpy.uarray([], 0.001)
# plt.errorbar(noms(d2), noms(I_5), xerr=devs(d2), yerr=devs(I_5), fmt='bo', label='data')
# plt.grid()
# plt.show()
#
d_min = ufloat(1.6e-2, 0.2e-2)
lamb = 2*d_min

d1 = ufloat(2.3, 0.01) * 10 ** (-2)
d1_min_array = unumpy.uarray([2.7 , 4.6, 6.4, 8.3, 10.1],0.1) * 10** (-2)
d1_min = diff(d1_min_array) / (len(d1_min_array)-1)
tmp = 1/(lamb**2)-1/(2*d1_min)**2
lamb_g_min = sqrt(1/tmp)

d1_max_array = unumpy.uarray([3.7 , 5.5, 7.4, 9.3, 11.1],0.1) * 10** (-2)
d1_max = diff(d1_max_array) / (len(d1_max_array)-1)
tmp = 1/(lamb**2)-1/(2*d1_max)**2
lamb_g_max1 = sqrt(1/tmp)

d2 = ufloat(2.7, 0.01) * 10 ** (-2)
d2_max_array = unumpy.uarray([5.7, 7.4],0.1) * 10** (-2)
d2_max = diff(d2_max_array) / (len(d2_max_array)-1)
tmp = 1/(lamb**2)-1/(2*d2_max)**2
lamb_g_max2 = sqrt(1/tmp)

d3 = ufloat(3, 0.01) * 10 ** (-2)
d3_max_array = unumpy.uarray([6.4, 7.8],0.1) * 10** (-2)
d3_max = diff(d3_max_array) / (len(d3_max_array)-1)
tmp = 1/(lamb**2)-1/(2*d3_max)**2
lamb_g_max3 = sqrt(1/tmp)

d4 = ufloat(2, 0.01) * 10 ** (-2)
d4_max_array = unumpy.uarray([7, 9.2],0.1) * 10** (-2)
d4_max = diff(d4_max_array) / (len(d4_max_array)-1)
tmp = 1/(lamb**2)-1/(2*d4_max)**2
lamb_g_max4 = sqrt(1/tmp)

d5 = ufloat(1.7, 0.01) * 10 ** (-2)
d5_min_array = unumpy.uarray([6.3, 8.9],0.1) * 10** (-2)
d5_min = diff(d5_min_array) / (len(d5_min_array)-1)
tmp = 1/(lamb**2)-1/(2*d5_min)**2
lamb_g_min5 = sqrt(1/tmp)

# node_deltas = unumpy.uarray([], [])
# antinode_deltas = unumpy.uarray([], [])
#
# lamb_g_2 = unumpy.uarray([], [])
# d3 = unumpy.uarray([], [])
# y = 1/(lamb_g_2**2)
# x = (1/(2*d3))**2
# reg2 = linregress(x, y)
# plt.errorbar(noms(x), noms(y), xerr=devs(x), yerr=devs(y), fmt='bo', label='data')
# lamb_2 = ufloat(reg2.intercept, reg2.intercept_stderr)
#
L = 15e-2
#
def calc_d(delta_phi):
    # tmp: 1/(2d)**2
    tmp = 1/(2.8e-2**2)-(delta_phi/(2*np.pi*L)-1/2.8e-2)**2
    d = 1/2*np.sqrt(1/tmp)
    return d

d_circ = calc_d(np.pi/2)
d_lin = calc_d(2*np.pi)
d_lin2 = calc_d(np.pi)
d_lin3 = calc_d(4*np.pi)
#
I_circ = unumpy.uarray([0.31, 0.347, 0.346, 0.316, 0.361, 0.323,0.311,  0.309, 0.317, 0.34, 0.34,0.334, 0.332, 0.354 ], 0.003)
theta_circ = unumpy.uarray(range(0, 4*len(noms(I_circ) + 4), 4), 1) * angle_delta
theta_circ = unumpy.uarray(np.deg2rad(noms(theta_circ)), np.deg2rad(devs(theta_circ)))

#

I_lin = unumpy.uarray([0.117, 0.09, 0.066, 0.134, 0.183, 0.203, 0.236, 0.268, 0.216, 0.176, 0.087, 0.078, 0.155, 0.222 ], 0.002)
theta_lin = unumpy.uarray(range(0, 4*len(noms(I_lin) + 4), 4), 1) * angle_delta
theta_lin = unumpy.uarray(np.deg2rad(noms(theta_lin)), np.deg2rad(devs(theta_lin)))

#
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(noms(theta_circ), noms(I_circ), 'bo')
ax.grid(True)
plt.show()

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(noms(theta_lin), noms(I_lin), 'bo')
ax.grid(True)
plt.show()

#

#%% michelson interferometer

x_mich = unumpy.uarray([], [])
I_mich = unumpy.uarray([], [])
plt.errorbar(noms(x_mich), noms(I_mich), devs(x_mich), devs(I_mich), fmt='bo')
plt.xlabel('x [m]')
plt.ylabel('I [W]')
plt.grid()
plt.show()

#%% fabry-perot interferometer

d1_fp_1 = ufloat(0, 0)
min_points_count_1 = 10
d2_fp_1 = ufloat(0, 0)

# d2-d1=min_points*(lamb/2)
lamb_fp_1 = 2 * (d2_fp_1 - d1_fp_1) / min_points_count_1

d1_fp_2 = ufloat(0, 0)
min_points_count_2 = 10
d2_fp_2 = ufloat(0, 0)

lamb_fp_2 = 2 * (d2_fp_2 - d1_fp_2) / min_points_count_2

#%% lloyd interferometer

h1 = ufloat(0, 0)
I_h1 = ufloat(0, 0)
h1_neighbours = unumpy.uarray([], [])
I_h1_neighbours = unumpy.uarray([], [])

h2 = ufloat(0, 0)
I_h2 = ufloat(0, 0)
h2_neighbours = unumpy.uarray([], [])
I_h2_neighbours = unumpy.uarray([], [])

h_max_neighbours = unumpy.uarray([], [])
I_h_max_neighbours = unumpy.uarray([], [])

h = np.concatenate(h1_neighbours, h2_neighbours, h_max_neighbours)  # unumpy.uarray([], [])
I_lloyd = np.concatenate(I_h1_neighbours, I_h2_neighbours, I_h_max_neighbours)  # unumpy.uarray([], [])

plt.errorbar(noms(h), noms(I_lloyd), devs(h), devs(I_lloyd), fmt='bo', label='data')
plt.xlabel('h [m]')
plt.ylabel('I [W]')
plt.grid()
plt.show()

# take into account 10 cm of effective source & detector
d1_lloyd = ufloat(0, 0)
AC = 2*d1_lloyd

# calculate lambda somehow? (maybe find several minimums and maximums, similar to FP)
# lamb_lloyd =

# random calculations (prob not needed)
lamb_manufacturer = 2.8e-2
wavelength_delta = (unumpy.sqrt((AC/2)**2 + h**2)*2-AC) / lamb_manufacturer
phase_delta = wavelength_delta*2*np.pi
