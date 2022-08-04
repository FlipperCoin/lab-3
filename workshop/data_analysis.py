# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 09:53:31 2022

@author: Lab3
"""

import numpy as np # math functions
import scipy # scientific functions
import matplotlib.pyplot as plt # for plotting figures and setting their properties
import pandas as pd # handling data structures (loaded from files)
from scipy.stats import linregress # contains linregress (for linear regression)
from scipy.optimize import curve_fit as cfit # non-linear curve fitting

#%% potential
def potential(x,y,a,C):
    return C * np.log( np.sqrt(np.power(x-a,2) + np.power(y,2)) / a) - C * np.log( np.sqrt(np.power(x+a,2) + np.power(y,2)) / a)

C = 1
a = 1
L = 3
N = 100
coord = np.linspace(-L, L , N) # defines coordinates
coord_x, coord_y = np.meshgrid(coord, coord)

    
V_xy = potential(coord_x, coord_y, a, C) # calculate potential

plt.figure() # create empty figure
plt.grid()
plt.pcolormesh(coord_x, coord_y, V_xy) # create mesh color plot
plt.colorbar() # add color bar
plt.contour(coord_x, coord_y, V_xy, np.sort([-1 ,-0.75, -0.6,  -0.5, -0.25,  0 , 0.25, 0.5, 0.75,  1]), cmap='hot')
#cmap changes the colormap to the hot preset
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.show() # show figure


x = np.linspace(-L,L,N)
V_x = potential(x, [0] * N, a, C)
plt.figure()
plt.grid()
plt.plot(x,V_x,'.', label="calculated potential")
plt.xlabel("x [m]")
plt.ylabel("V [V]")
plt.show()

#%% capacitor
def V_decay(t,a,b):
    return a*np.exp(-t*b)
def linearCurve(x,a,b):
    return a*x + b

eps0 = 8.854e-12 # F/m
D = 18e-2 # m
d = 0.5e-3 # m
C_theoretical = eps0 * np.pi * np.power(D,2) / (4 * d)
R = 977 # Ohm
R_total = 38.4e3 # Ohm
tau_theoretical = R_total * C_theoretical

C_data = pd.read_csv('capacitor.csv')
C_data = C_data.rename(columns = {"time (sec)":"t", "ch2":"V_R"})
C_data["V_C"] = C_data["ch1"] - C_data["V_R"]
t = np.array(C_data['t'].values)
V_C = np.array(C_data['V_C'].values)

plt.figure()
plt.grid()
plt.plot(t,V_C,label="data")

# curve fitting means fitting a curve of certain type that is as close as possible to all points of a given data set.

fit2 = cfit(V_decay,C_data['t'], C_data["V_C"]) # needs to find a = V_0 and b = 1/tau

plt.plot(t, V_decay(t,fit2[0][0],fit2[0][1]),label="fitted curve")
plt.legend()
plt.xlabel("t [sec]")
plt.ylabel("V [V]")
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.show()

plt.figure()
plt.grid()
plt.plot(t,np.log(V_C))
plt.xlabel("t [sec]")
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel("log(V) [V]")
plt.show()

# log of V_C should be linear in time
plt.figure()
plt.grid()
t1 = 0
t2 = 0.000047
inds = (t > t1) & (t < t2)
plt.plot(t[inds], np.log(V_C)[inds],'.', label="data")
plt.xlabel("time [sec]")
plt.ylabel("log(V) [V]")
plt.show()

# linear regression means fitting a linear curve to a set of data points.
reg2 = linregress(t[inds], np.log(V_C[inds]))
print(reg2)
V_0_reg = np.exp(reg2.intercept)
tau_reg = -1/reg2.slope

C_data["int_V_R"] = scipy.integrate.cumtrapz(C_data["V_R"], x = t, initial = 0)
int_V_R = np.array(C_data["int_V_R"].values)

plt.figure()
plt.plot(int_V_R,V_C, label='data')
reg3 = linregress(int_V_R, V_C)
C_meas = int_V_R /(R*V_C)

plt.plot(C_data["int_V_R"], linearCurve(C_data["int_V_R"], reg3.slope, reg3.intercept), label='regression')
plt.xlabel("integral of V_R [V*sec]")
plt.ylabel("V_C [V]")
plt.legend()
plt.grid()
plt.show()

#%% Ohm
R1 = 5.48 #Ohm

def I_R(V2,R1):
    return V2/R1

def V_R(V1,V2):
    return V1-V2

def R_t(V_R,I_R):
    return V_R/I_R

def P_t(V_R, I_R):
    return V_R * I_R

def Energy(P_t, t):
    return scipy.integrate.cumtrapz(P_t, x = t, initial = 0)

R_data=pd.read_csv("ohm.csv", header=1, usecols=[3,4,5])
# Header is the number of the row in the csv file that contains the column names, it is needed since the given file contains header names that should not be used as data.
R_data = R_data.rename(columns = {"Time (s)":"t", "1 (VOLT)":"V1", "2 (VOLT)":"V2"})
V1 = R_data["V1"]
V2 = R_data["V2"]
t = R_data["t"]
VR = V_R(V1,V2)
IR = I_R(V2,R1)
Rt = R_t(VR,IR)
Et = Energy(P_t(VR,IR),t)
indices = (Et > 0.01) & (Et < 0.9)
Et = Et[indices]
Rt = Rt[indices]
plt.figure()
plt.plot(Et, Rt)
plt.xlabel("E [J]")
plt.ylabel("R [Ohm]")
plt.grid()
plt.show()
reg = linregress(Et, Rt)
R0 = reg.intercept
aC = reg.slope/R0

#%% inductance
def flux(voltage, time):
    return  scipy.integrate.cumtrapz(voltage, x = time, initial = 0)
def max_flux_index(flux):
    return np.argmax(flux)

h = np.array([30,24,18,14,8]) * 1e-2

Ind_data = []
plt.figure()
plt.grid()
for n in range(0,5):
    df = pd.read_csv('Trace %d.csv'%n, header = 1)
    df = df.rename(columns = {"Time (s)":"t", "1 (VOLT)":"ref", "2 (VOLT)":"signal"})
    Ind_data.append(df)
    plt.plot(df["t"], df["ref"], label=f"ref  {h[n]}")
    plt.plot(df["t"], df["signal"], label=f"signal  {h[n]}")
    plt.xlabel("t [sec]")
    plt.ylabel("V [sec]")
plt.legend()
plt.show()

plt.figure()
plt.grid()
t_coil = np.array([])
for n in range(0,5):
    df = Ind_data[n]
    ref_flux = flux(df["ref"],df["t"])
    signal_flux = flux(df["signal"],df["t"])
    ref_ind = max_flux_index(ref_flux)
    signal_ind = max_flux_index(signal_flux)  
    df["t"] = df["t"] - df["t"][ref_ind]     
    plt.plot(df["t"], ref_flux,  label = f"ref  {h[n]}")
    plt.plot(df["t"][ref_ind], ref_flux[ref_ind],"ro")
    plt.plot(df["t"], signal_flux , label = f"signal  {h[n]}")
    plt.plot(df["t"][signal_ind], signal_flux[signal_ind],"ro")
    t_coil = np.append(t_coil, df["t"][signal_ind] - df["t"][ref_ind])
    plt.xlabel("t [sec]")
    plt.ylabel("flux [V*sec]")
plt.legend()
plt.show()

plt.figure()
y_err = np.sqrt( np.power(1e-3/t_coil,2) + np.power(1e-3 * h/(t_coil**2),2)  )
plt.errorbar(t_coil,h/t_coil,y_err,1e-3,"ro")
reg = linregress(t_coil,h/t_coil)
v0 = reg.intercept
a = reg.slope * 2
plt.plot(t_coil, reg.slope * t_coil + reg.intercept)
plt.xlabel("t_coil [sec]")
plt.ylabel("h/t_coil [m/s]")
plt.grid()
plt.show()
R_squared = reg.rvalue ** 2

chi_squared = scipy.stats.chisquare(h/t_coil, reg.slope * t_coil + reg.intercept).statistic
p = reg.pvalue

