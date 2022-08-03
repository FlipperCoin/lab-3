# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:08:57 2022

@author: student
"""

import numpy as np # math functions
import scipy # scientific functions
import scipy.stats # contains linregress (for linear regression)
import matplotlib.pyplot as plt # for plotting figures and setting their properties
import pandas as pd # handling data structures (loaded from files)
from scipy.optimize import curve_fit as cfit # non-linear curve fitting

#%% potential

# def potential(x,y,a,C):
#     return C * np.log(np.sqrt(np.power(x-a,2) + np.power(y,2))/a) - C * (np.log(np.sqrt(np.power(x+a,2) + np.power(y,2))/a))

# C = 1
# a = 1
# L = 3
# N = 100
# coord = np.linspace(-L, L , N)
# coord_x, coord_y = np.meshgrid(coord, coord)

# V_xy = potential(coord_x, coord_y, a, C)

# plt.figure()
# plt.pcolormesh(coord_x, coord_y, V_xy)
# plt.colorbar()

# plt.contour(coord_x, coord_y, V_xy, np.sort([-1 , 0 , 1]), cmap='hot')

# plt.figure()

# x = np.linspace(-a, a, 100)
# V_x = potential(x, np.array([0]*100), a, C)
# V_x = plt.plot(x, V_x, '.')

#%% capacitor
# eps0 = 8.854e-12 # F/m
# D = 18e-2 # m
# d = 0.5e-3 # m

# C_theor = (eps0*np.pi*(D**2))/(4*d)

# R_total = 38.4e3 # Ohm
# R = 977 # Ohm
# tau_theor = R_total * C_theor

# # plot voltage over time

# C_data = pd.read_csv('capacitor.csv')
# C_data = C_data.rename(columns = {"time (sec)":"t", "ch2":"V_R"})
# C_data["V_C"] = C_data["ch1"] - C_data["V_R"]

# plt.plot(C_data["t"],C_data["V_C"])

# # curve fitting V0 & tau

# def V_decay(t,a,b):
#     return a*np.exp(-t*b)

# fit2 = cfit(V_decay,C_data['t'], C_data["V_C"])
# plt.plot(C_data['t'], V_decay(C_data['t'],fit2[0][0],fit2[0][1]))

# # plot log(voltage) over time

# plt.figure()
# plt.plot(C_data["t"], np.log(C_data["V_C"]))
# plt.grid("both")

# # plot log(voltage) over time, linear time

# plt.figure()
# t1 = 0
# t2 = 45e-6
# inds = (C_data['t'] > t1) & (C_data['t'] < t2)
# plt.plot(C_data['t'][inds], np.log(C_data["V_C"])[inds],'.')
# plt.grid("both")

# # linear regress for V0 & tau

# reg2 = scipy.stats.linregress(C_data['t'][inds], np.log(C_data["V_C"][inds]))
# print(reg2)

# tau_reg = -1/reg2[0]
# V_0_reg = np.exp(reg2[1])

# # vector integral of V_R over time

# C_data["int_V_R"] = scipy.integrate.cumtrapz(C_data["V_R"], C_data["t"], initial = 0)

# plt.figure()
# plt.plot(C_data["int_V_R"], C_data["V_C"])

# t3 = 0
# t4 = 0.000168
# inds2 = (C_data['t'] > t3) & (C_data['t'] < t4)
# plt.figure()
# plt.plot(C_data["int_V_R"][inds2], C_data["V_C"][inds2])

# # linreg for 1/RC

# reg3 = scipy.stats.linregress(C_data["int_V_R"][inds2], C_data["V_C"][inds2])
# print(reg3)

# C_meas = (1/(reg3.slope))/R

# plt.plot(C_data["int_V_R"],C_data["int_V_R"]*reg3.slope + reg3.intercept)
# plt.xlabel("integral of V_R")
# plt.ylabel("V_C")
# plt.legend(["data","regression"])
# plt.grid("both")

#%% ohm

# def I_R(V2, R1):
#     return V2/R1

# def V_R(V1, V2):
#     return V1-V2

# def R_t(V_R, I_R):
#     return V_R/I_R

# def P_t(V_R, I_R):
#     return V_R*I_R

# def Energy(P_t, t):
#     return scipy.integrate.cumtrapz(P_t, t, initial = 0)

# R1 = 5.48

# R_data=pd.read_csv("ohm.csv", header=1)

# R_data = R_data.rename(columns = {"Time (s)":"t", "1 (VOLT)":"V1", "2 (VOLT)":"V2"})
# R_data["V_R"] = V_R(R_data["V1"],R_data["V2"])
# R_data["I_R"] = I_R(R_data["V2"],R1)
# R_data["R_t"] = R_t(R_data["V_R"],R_data["I_R"])
# R_data["P_t"] = P_t(R_data["V_R"],R_data["I_R"])
# R_data["Energy"] = Energy(R_data["P_t"], R_data["t"])

# plt.plot(R_data["Energy"], R_data["R_t"])
# plt.legend(["R_t"])

# t1 = 0.01
# t2 = 0.103
# inds = (R_data['t'] > t1) & (R_data['t'] < t2)

# plt.figure()
# plt.plot(R_data["Energy"][inds], R_data["R_t"][inds])
# plt.legend(["R_t"])

# reg2 = scipy.stats.linregress(R_data["Energy"][inds], R_data["R_t"][inds])

# R0 = reg2.intercept
# alpha_over_Cheat = reg2.slope/R0

#%% inductance

def magnetic_flux(V_t,t):
    return scipy.integrate.cumtrapz(V_t, t, initial = 0)

def max_flux(flux,t):
    return t[np.argmax(flux)]

def error_y_over_x(y,x,yerr,xerr):
    return np.sqrt(pow((1/x)*yerr,2) + pow((-y/pow(x,2))*xerr,2))

L = np.array([30,24,18,14,8])
L = L*1e-2
Ind_data = []
for n in range(0,5,1):
    df = pd.read_csv('Trace %d.csv'%(n,),header = 1)
    df = df.rename(columns = {"Time (s)":"t", "1 (VOLT)":"ref", "2 (VOLT)":"signal"})
    df["ref_flux"] = magnetic_flux(df["ref"],df["t"])
    df["signal_flux"] = magnetic_flux(df["signal"],df["t"])
    df["t_coil"] = max_flux(df["signal_flux"],df["t"]) - max_flux(df["ref_flux"],df["t"])
    Ind_data.append(df)

t_coil = np.array([df["t_coil"][0] for df in Ind_data])
L_over_tcoil = L/t_coil
ferrors =  error_y_over_x(L,t_coil,2e-3,10e-3)
plt.errorbar(t_coil, L_over_tcoil,ferrors)

reg = scipy.stats.linregress(t_coil,L_over_tcoil)
a = reg.slope*2
v0 = reg.intercept
plt.plot(t_coil,(t_coil*a/2)+v0)