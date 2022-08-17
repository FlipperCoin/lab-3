import numpy as np
from scipy.constants import R, Avogadro, proton_mass, neutron_mass
import matplotlib.pyplot as plt

# 1
# v = sqrt(gamma*R*T/M)
T = 0 + 273.15
v = 331.7
gammaR_over_M = v**2/T
v_air = np.sqrt(gammaR_over_M*(20 + 273.15))

# 2
m = 2 * proton_mass + 2 * neutron_mass
M = m * Avogadro
gamma = 1.66
v_helium = np.sqrt(gamma*R*(20 + 273.15)/M)

# 3
m = 22 * proton_mass + 22 * neutron_mass
M = m * Avogadro
gamma = 1.304
v_co2 = np.sqrt(gamma*R*(20 + 273.15)/M)

# 4
L = 1
# T <= 1/2*L/v, f >= 2*v/L
T_air = L/v_air/2
T_helium = L/v_helium/2
T_co2 = L/v_co2/2
f_air = 1/T_air
f_helium = 1/T_helium
f_co2 = 1/T_co2

# 5_1
L = 20e-2
x = np.linspace(0, L, 1000)
a = 1
v = v_air
f = 5e3
P = np.abs(2*a*np.cos(2*np.pi*f*L/(2*v))*np.cos(np.pi*(f/v)*(2*x-L)))

plt.plot(x, P)
plt.grid()
plt.xlabel(r"x $\left[ m \right]$")
plt.ylabel(r"$\Delta P \left( x \right)$ $\left[ \frac{N}{m^2} \right]$")
plt.show()

# 5_2
phi0 = np.pi
P = np.abs(2*a*np.cos(2*np.pi*f*L/(2*v)-phi0/2)*np.cos(np.pi*(f/v)*(2*x-L)+phi0/2))

plt.plot(x, P)
plt.grid()
plt.xlabel(r"x $\left[ m \right]$")
plt.ylabel(r"$\Delta P \left( x \right)$ $\left[ \frac{N}{m^2} \right]$")
plt.show()
