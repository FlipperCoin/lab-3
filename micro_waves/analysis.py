import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from uncertainties.umath import sqrt
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.integrate import cumtrapz

# %% grid polarizer

# deltaLines = 4.75 degrees, dev = 0.25
angle_delta = ufloat(360 / 72, 0.25)
# 17 lines per 1/4
theta_1 = unumpy.uarray([0, 5, 10, 15, 18, 23, 28, 33, 36, 41, 46, 51, 54, 59, 64, 69, 72], 1) * angle_delta
theta_1 = unumpy.uarray(np.deg2rad(noms(theta_1)), np.deg2rad(devs(theta_1)))
I_1 = unumpy.uarray(
    [0, 0.076, 0.211, 0.272, 0.283, 0.263, 0.183, 0.039, 0.0001, 0.098, 0.219, 0.275, 0.282, 0.253, 0.171, 0.036, 0],
    [0.001, 0.002, 0.002, 0.001, 0.001, 0.001, 0.002, 0.001, 0.0001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
     0.0001])

x1 = unumpy.cos(theta_1) ** 2
reg = linregress(noms(x1), noms(I_1))
plt.errorbar(noms(x1), noms(I_1), xerr=devs(x1), yerr=devs(I_1), fmt='bo', label='data')
plt.plot(noms(x1), reg.intercept + reg.slope * noms(x1), label='fit')
plt.xlabel(r'$cos^2(\theta)$')
plt.ylabel(r'V $\left[ V \right]$')
plt.grid()
plt.legend()
plt.show()

theta_2 = unumpy.uarray(list(range(0, 105, 15)), 2)
theta_2 = unumpy.uarray(np.deg2rad(noms(theta_2)), np.deg2rad(devs(theta_2)))
I_2 = unumpy.uarray([0.021, 0.009, 0.03, 0.135, 0.214, 0.253, 0.261],
                    [0.001, 0.001, 0.002, 0.001, 0.001, 0.001, 0.001])

x2 = unumpy.cos(theta_2) ** 4
reg2 = linregress(noms(x2), noms(I_2))
plt.errorbar(noms(x2), noms(I_2), xerr=devs(x2), yerr=devs(I_2), fmt='bo', label='data')
plt.plot(noms(x2), reg2.intercept + reg2.slope * noms(x2), label='fit')
plt.xlabel(r'$cos^4(\theta)$')
plt.ylabel(r'V $\left[ V \right]$')
plt.grid()
plt.legend()
plt.show()

# #%% waveguide

# perpendicular
d1 = unumpy.uarray([1.5, 1.8, 2, 2.5, 3, 3.3, 3.6, 4, 4.2, 2.9], 0.2) * 1e-2
I_3 = unumpy.uarray([0.25, 0.243, 0.274, 0.235, 0.286, 0.362, 0.375, 0.376, 0.39, 0.252], 0.002)
plt.errorbar(noms(d1), noms(I_3), xerr=devs(d1), yerr=devs(I_3), fmt='bo', label='data')
reg_perp = linregress(noms(d1), noms(I_3))
plt.plot(noms(d1), noms(d1)*reg_perp.slope+reg_perp.intercept, label='fit')
plt.xlabel(r'd $\left[ m \right] $')
plt.ylabel(r'V $\left[ V \right]$')
plt.legend()
plt.grid()
plt.show()
#
theta_3 = unumpy.uarray([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55], 1) * angle_delta - 90
theta_3 = unumpy.uarray(np.deg2rad(noms(theta_3)), np.deg2rad(devs(theta_3)))
I_4 = unumpy.uarray([0.003, 0.143, 0.33, 0.357, 0.379, 0.356, 0.172, 0.004, 0.18, 0.335, 0.365, 0.375], 0.002)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.errorbar(noms(theta_3), noms(I_4), xerr=devs(theta_3), yerr=devs(I_4), fmt='bo')
# ax.plot(noms(theta_3), np.max(noms(I_4))*np.abs(np.sin(noms(theta_3))))
ax.grid(True)
plt.show()

# parallel
d2 = unumpy.uarray([1.7, 2, 2.5, 2.7, 3, 3.5, 4, 4.5, 5, 3.2], 0.2) * 1e-2
I_5 = unumpy.uarray([0.265, 0.326, 0.39, 0.367, 0.403, 0.414, 0.441, 0.443, 0.428, 0.405], 0.001)
plt.errorbar(noms(d2), noms(I_5), xerr=devs(d2), yerr=devs(I_5), fmt='bo', label='data')
plt.xlabel(r'd $\left[ m \right] $')
plt.ylabel(r'V $\left[ V \right]$')
plt.grid()
plt.show()

d_min = ufloat(1.6e-2, 0.2e-2)
lamb = 2 * d_min
lamb_manu = 2.8e-2

def diff(arr):
    return arr[len(arr) - 1] - arr[0]

d1 = ufloat(2.3, 0.1) * 1e-2
d1_min_array = unumpy.uarray([2.7, 4.6, 6.4, 8.3, 10.1], 0.1) * 1e-2
lamb_g_1 = 2 * diff(d1_min_array) / (len(d1_min_array) - 1)
# lamb_g_1_squared = 1/(lamb_d1**2)-1/(2*d1)**2

d1_max_array = unumpy.uarray([3.7, 5.5, 7.4, 9.3, 11.1], 0.1) * 1e-2
lamb_d1_max = 2 * diff(d1_max_array) / (len(d1_max_array) - 1)
# lamb_g_1_max_squared = 1/(lamb_d1_max**2)-1/(2*d1)**2

d1_array = np.concatenate([d1_min_array, d1_max_array])
d1_array = d1_array[np.argsort(noms(d1_array))]
d1_array = d1_array - d1_array[0]
n = np.arange(0, len(d1_array))
reg_d1 = linregress(n, noms(d1_array))
lamb_g_1_fromreg = 4 * ufloat(reg_d1.slope, reg_d1.stderr)

plt.errorbar(n, noms(d1_array), yerr=devs(d1_array), fmt='bo', label='data')
plt.plot(n, n*reg_d1.slope+reg_d1.intercept, label='fmt')
plt.xlabel(r"n $\left[ 1 \right]$")
plt.ylabel(r"$\Delta x$ $\left[ m \right]$")
plt.legend()
plt.grid()
plt.show()

d2 = ufloat(2.7, 0.1) * 1e-2
d2_max_array = unumpy.uarray([5.7, 7.4], 0.1) * 1e-2
lamb_g_2 = 2 * diff(d2_max_array) / (len(d2_max_array) - 1)
# lamb_g_2_squared = 1/(lamb_d2**2)-1/(2*d2)**2

d3 = ufloat(3, 0.1) * 1e-2
d3_max_array = unumpy.uarray([6.4, 7.8], 0.1) * 1e-2
lamb_g_3 = 2 * diff(d3_max_array) / (len(d3_max_array) - 1)
# lamb_g_3_squared = 1/(lamb_d3**2)-1/(2*d3)**2

d4 = ufloat(2, 0.1) * 1e-2
d4_max_array = unumpy.uarray([7, 9.2], 0.1) * 1e-2
lamb_g_4 = 2 * diff(d4_max_array) / (len(d4_max_array) - 1)
# lamb_g_4_squared = 1/(lamb_d4**2)-1/(2*d4)**2

d5 = ufloat(1.7, 0.1) * 1e-2
d5_min_array = unumpy.uarray([6.3, 8.9], 0.1) * 1e-2
lamb_g_5 = 2 * diff(d5_min_array) / (len(d5_min_array) - 1)
# lamb_g_5_squared = 1/(lamb_d5**2)-1/(2*d5)**2

lamb_g = np.array([lamb_g_1_fromreg, lamb_g_2, lamb_g_3,
                   lamb_g_4, lamb_g_5])
d_lambda_g = np.array([d1, d2, d3, d4, d5])

y = 1 / (lamb_g**2)
x = 1 / ((2 * d_lambda_g) ** 2)
reg2 = linregress(noms(x), noms(y))
plt.errorbar(noms(x), noms(y), xerr=devs(x), yerr=devs(y), fmt='bo', label='data')
plt.plot(noms(x), noms(x) * reg2.slope + reg2.intercept, label='fit')
lamb_2 = sqrt(1 / ufloat(reg2.intercept, reg2.intercept_stderr))
plt.xlabel(r'$\frac{1}{(2d)^2}$ $\left[ \frac{1}{m^2} \right]$')
plt.ylabel(r'$\frac{1}{\lambda^2_g}$ $\left[ \frac{1}{m^2} \right]$')
plt.grid()
plt.legend()
plt.show()
#
L = 15e-2

#
def calc_d(delta_phi):
    # tmp: 1/(2d)**2
    tmp = 1 / (2.8e-2 ** 2) - (delta_phi / (2 * np.pi * L) - 1 / 2.8e-2) ** 2
    d = 1 / 2 * np.sqrt(1 / tmp)
    return d

d_circ = calc_d(np.pi / 2)
d_lin = calc_d(2 * np.pi)
d_lin2 = calc_d(np.pi)
d_lin3 = calc_d(4 * np.pi)
#
I_circ = unumpy.uarray([0.31, 0.347, 0.346, 0.316, 0.361, 0.323, 0.311, 0.309, 0.317, 0.34, 0.34, 0.334, 0.332, 0.354],
                       0.003)
theta_circ = unumpy.uarray(range(0, 4 * len(noms(I_circ) + 4), 4), 1) * angle_delta * -1
theta_circ = unumpy.uarray(np.deg2rad(noms(theta_circ)), np.deg2rad(devs(theta_circ)))

#

I_lin = unumpy.uarray([0.117, 0.09, 0.066, 0.134, 0.183, 0.203, 0.236, 0.268, 0.216, 0.176, 0.087, 0.078, 0.155, 0.222],
                      0.002)
theta_lin = unumpy.uarray(range(0, 4 * len(noms(I_lin) + 4), 4), 1) * angle_delta * -1
theta_lin = unumpy.uarray(np.deg2rad(noms(theta_lin)), np.deg2rad(devs(theta_lin)))

#
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.errorbar(noms(theta_circ), noms(I_circ), xerr=devs(theta_circ), yerr=devs(I_circ), fmt='bo')
ax.grid(True)
plt.show()

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.errorbar(noms(theta_lin), noms(I_lin), xerr=devs(theta_lin), yerr=devs(I_lin), fmt='bo')
ax.grid(True)
plt.show()

d_dphi = ufloat(2.4, 0.1)*1e-2
lamb_g_dphi = 1/sqrt(1/(lamb_manu**2)-1/((2*d_dphi)**2))
dphi = 2*np.pi*L*((1/lamb_manu)-(1/lamb_g_dphi))

d_dphi = ufloat(4.6, 0.1)*1e-2
lamb_g_dphi = 1/sqrt(1/(lamb_manu**2)-1/((2*d_dphi)**2))
dphi = 2*np.pi*L*((1/lamb_manu)-(1/lamb_g_dphi))

theta = np.arange(0, 2+1/16, 1/16) * np.pi

Ex = np.sin(theta)
Ey = np.sin(theta+ np.pi/2 + 0.1)
arg = np.arctan2(Ey, Ex)
amp = np.sqrt(Ex**2 + Ey**2)
I_circ_theo = np.array([1/(2*np.pi)*cumtrapz((amp*np.cos(arg-pos))**2, theta)[-1] for pos in theta])

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, np.max(noms(I_circ))/np.max(noms(I_circ_theo))*I_circ_theo, 'ro', label='theo')
ax.errorbar(noms(theta_circ), noms(I_circ), xerr=devs(theta_circ), yerr=devs(I_circ), fmt='bo', label='data')
ax.legend()
ax.grid(True)
plt.show()


Ex = np.sin(theta)
Ey = np.sin(theta+0.6)
arg = np.arctan2(Ey, Ex)
amp = np.sqrt(Ex**2 + Ey**2)
I_lin_theo = np.array([1/(2*np.pi)*cumtrapz((amp*np.cos(arg-pos))**2, theta)[-1] for pos in theta])

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, np.max(noms(I_lin))/np.max(noms(I_lin_theo))*I_lin_theo, 'ro', label='theo')
ax.errorbar(noms(theta_lin), noms(I_lin), xerr=devs(theta_lin), yerr=devs(I_lin), fmt='bo', label='data')
ax.legend()
ax.grid(True)
plt.show()


# %% michelson interferometer

x_max_first = 0.76
x_mich = unumpy.uarray([0.706, 0.714, 0.723, 0.729, 0.735, 0.740, 0.745, 0.754, 0.759, 0.768], 0.001)
I_mich = unumpy.uarray([0.22, 0.00, 0.206, 0.0003, 0.205, 0.0003, 0.202, 0.0003, 0.19, 0.0003],
                       [0.01, 0.01, 0.003, 0.0005, 0.005, 0.0003, 0.003, 0.0005, 0.01, 0.0004])
plt.errorbar(noms(x_mich), noms(I_mich), devs(x_mich), devs(I_mich), fmt='bo')
plt.xlabel('x [m]')
plt.ylabel('I [W]')
plt.grid()
plt.show()

lambda_mich = (x_mich[len(x_mich) - 1] - x_mich[0]) * 4 / 9

# %% fabry-perot interferometer

d1_fp_1_x1 = 52.9
d1_fp_1_x2 = 46.9
d1_fp_1 = ufloat(d1_fp_1_x1 - d1_fp_1_x2, 0.1) * 1e-2
min_points_count_1 = 10
d2_fp_1_x1 = 67.1
d2_fp_1_x2 = 46.9
d2_fp_1 = ufloat(d2_fp_1_x1 - d2_fp_1_x2, 0.1) * 1e-2

# d2-d1=min_points*(lamb/2)
lamb_fp_1 = 2 * (d2_fp_1 - d1_fp_1) / min_points_count_1

d1_fp_2_x1 = 56.0
d1_fp_2_x2 = 46.9
d1_fp_2 = ufloat(d1_fp_2_x1 - d1_fp_2_x2, 0.1) * 1e-2
min_points_count_2 = 10
d2_fp_2_x1 = 69.9
d2_fp_2_x2 = 46.9
d2_fp_2 = ufloat(d2_fp_2_x1 - d2_fp_2_x2, 0.1) * 1e-2

lamb_fp_2 = 2 * (d2_fp_2 - d1_fp_2) / min_points_count_2

# %% lloyd interferometer

h_base = 60e-2 - 5.5e-2
h1 = ufloat(61.5e-2, 0.1e-2) - h_base
I_h1 = ufloat(0.138, 0.003)
h1_neighbours = unumpy.uarray([61.2e-2, 61.3e-2, 61.4e-2, 61.5e-2, 61.6e-2, 61.7e-2, 61.8e-2],
                              [0.1e-2, 0.1e-2, 0.1e-2, 0.1e-2, 0.1e-2, 0.1e-2, 0.1e-2]) - h_base
I_h1_neighbours = unumpy.uarray([0.144, 0.143, 0.137, 0.136, 0.137, 0.138, 0.140],
                                [0.002, 0.001, 0.001, 0.002, 0.002, 0.001, 0.001])

plt.errorbar(noms(h1_neighbours),
             noms(I_h1_neighbours),
             xerr=devs(h1_neighbours),
             yerr=devs(I_h1_neighbours), fmt='bo', label='data')
plt.xlabel('h [m]')
plt.ylabel('I [W]')
plt.grid()
plt.show()

h2 = ufloat(66.2e-2, 0.1e-2) - h_base
I_h2 = ufloat(0.173, 0.001)
h2_neighbours = unumpy.uarray(
    [65.5e-2, 65.6e-2, 65.7e-2, 65.8e-2, 65.9e-2, 66.0e-2, 66.1e-2, 66.2e-2, 66.3e-2, 66.4e-2],
    0.1e-2) - h_base
I_h2_neighbours = unumpy.uarray([0.174, 0.173, 0.172, 0.171, 0.171, 0.171, 0.172, 0.173, 0.175, 0.177],
                                [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001])

plt.errorbar(noms(h2_neighbours),
             noms(I_h2_neighbours),
             xerr=devs(h2_neighbours),
             yerr=devs(I_h2_neighbours), fmt='bo', label='data')
plt.xlabel('h [m]')
plt.ylabel('I [W]')
plt.grid()
plt.show()

h_max_neighbours = unumpy.uarray([64.2e-2, 64.1e-2, 64e-2, 63.9e-2, 63.8e-2, 63.7e-2, 63.6e-2, 63.5e-2, 63.4e-2, ],
                                 0.1e-2) - h_base
I_h_max_neighbours = unumpy.uarray([0.211, 0.212, 0.212, 0.212, 0.213, 0.213, 0.212, 0.210, 0.207],
                                   [0.001, 0.001, 0.001, 0.001, 0.001, 0.002, 0.003, 0.001, 0.001])

plt.errorbar(noms(h_max_neighbours),
             noms(I_h_max_neighbours),
             xerr=devs(h_max_neighbours),
             yerr=devs(I_h_max_neighbours), fmt='bo', label='data')
plt.xlabel('h [m]')
plt.ylabel('I [W]')
plt.grid()
plt.show()
#
# h = np.concatenate([h1_neighbours, h2_neighbours, h_max_neighbours])  # unumpy.uarray([], [])
# I_lloyd = np.concatenate([I_h1_neighbours, I_h2_neighbours, I_h_max_neighbours])  # unumpy.uarray([], [])
#
# plt.errorbar(noms(h), noms(I_lloyd), devs(h), devs(I_lloyd), fmt='bo', label='data')
# plt.xlabel('h [m]')
# plt.ylabel('I [W]')
# plt.grid()
# plt.show()

# take into account 10 cm of effective source & detector
d1_lloyd = ufloat(36, 1) * 1e-2
AC = 2 * d1_lloyd

# calculate lambda somehow? (maybe find several minimums and maximums, similar to FP)
min1 = ufloat(0.070, 0.001)
min2 = ufloat(0.114, 0.001)
path1 = 2 * sqrt(d1_lloyd ** 2 + min1 ** 2)
path2 = 2 * sqrt(d1_lloyd ** 2 + min2 ** 2)
lamb_lloyd = path2 - path1

min1 = ufloat(0.070, 0.001)
max1 = ufloat(0.093, 0.001)
path1 = 2 * sqrt(d1_lloyd ** 2 + min1 ** 2)
path2 = 2 * sqrt(d1_lloyd ** 2 + max1 ** 2)
lamb_lloyd2 = 2 * (path2 - path1)

# asymmetric
min1 = ufloat(0.070, 0.001)
min2 = ufloat(0.114, 0.001)
path1 = sqrt(d1_lloyd ** 2 + min1 ** 2) + sqrt((d1_lloyd - 1e-2) ** 2 + min1 ** 2)
path2 = sqrt(d1_lloyd ** 2 + min2 ** 2) + sqrt((d1_lloyd - 1e-2) ** 2 + min2 ** 2)
lamb_lloyd3 = path2 - path1

# %% lloyd interferometer 2
#
# h_base = 60e-2 - 5.5e-2
# h1_neighbours = unumpy.uarray([61.8, 61.7, 61.6, 61.9, 62],
#                               0.1e-2) - h_base
# I_h1_neighbours = unumpy.uarray([0.15, 0.155, 0.156, 0.151, 0.152],
#                                 [0.003, 0.002, 0.002, 0.002, 0.002])
#
# plt.errorbar(noms(h1_neighbours),
#              noms(I_h1_neighbours),
#              xerr=devs(h1_neighbours),
#              yerr=devs(I_h1_neighbours), fmt='bo', label='data')
# plt.xlabel('h [m]')
# plt.ylabel('I [W]')
# plt.grid()
# plt.show()
#
# h2_neighbours = unumpy.uarray([63.9, 64, 64.1, 64.2, 64.3],
#                               0.1e-2) - h_base
# I_h2_neighbours = unumpy.uarray([0.215, 0.218, 0.217, 0.218, 0.217],
#                                 [0.002, 0.002, 0.002, 0.002, 0.002])
#
# plt.errorbar(noms(h2_neighbours),
#              noms(I_h2_neighbours),
#              xerr=devs(h2_neighbours),
#              yerr=devs(I_h2_neighbours), fmt='bo', label='data')
# plt.xlabel('h [m]')
# plt.ylabel('I [W]')
# plt.grid()
# plt.show()
#
# h_max_neighbours = unumpy.uarray([],
#                                  0.1e-2) - h_base
# I_h_max_neighbours = unumpy.uarray([],
#                                    [])
#
# plt.errorbar(noms(h_max_neighbours),
#              noms(I_h_max_neighbours),
#              xerr=devs(h_max_neighbours),
#              yerr=devs(I_h_max_neighbours), fmt='bo', label='data')
# plt.xlabel('h [m]')
# plt.ylabel('I [W]')
# plt.grid()
# plt.show()
#
# # take into account 10 cm of effective source & detector
# d2_lloyd = ufloat(35, 1) * 1e-2
# AC = 2 * d2_lloyd
#
# # calculate lambda somehow? (maybe find several minimums and maximums, similar to FP)
# min1 = ufloat(0.070, 0.001)
# min2 = ufloat(0.114, 0.001)
# path1 = 2 * sqrt(d2_lloyd ** 2 + min1 ** 2)
# path2 = 2 * sqrt(d2_lloyd ** 2 + min2 ** 2)
# lamb_lloyd = path2 - path1
#
# min1 = ufloat(0.070, 0.001)
# max1 = ufloat(0.093, 0.001)
# path1 = 2 * sqrt(d2_lloyd ** 2 + min1 ** 2)
# path2 = 2 * sqrt(d2_lloyd ** 2 + max1 ** 2)
# lamb_lloyd2 = 2 * (path2 - path1)
