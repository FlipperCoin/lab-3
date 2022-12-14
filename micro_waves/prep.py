import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

lamb = 2.8e-2
L = 15e-2

def calc_d(delta_phi):
    # tmp: 1/(2d)**2
    tmp = 1/(lamb**2)-(delta_phi/(2*np.pi*L)-1/lamb)**2
    d = 1/2*np.sqrt(1/tmp)
    return d

d_circ = calc_d(np.pi/2)
d_lin = calc_d(2*np.pi)
d_lin2 = calc_d(np.pi)
d_lin3 = calc_d(4*np.pi)

omega_t = np.arange(0, 2+1/16, 1/16) * np.pi

Ex = np.sin(omega_t)
Ey = np.sin(omega_t+np.pi/2)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(np.arctan2(Ey, Ex), np.sqrt(Ex**2+Ey**2), 'bo')
ax.grid(True)
plt.show()

Ex = 0
Ey = np.sin(omega_t+np.pi/2)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(np.arctan2(Ey, Ex), np.sqrt(Ex**2+Ey**2), 'bo')
ax.grid(True)
plt.show()

Ex = 2*np.sin(omega_t)
Ey = np.sin(omega_t+np.pi/2)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(np.arctan2(Ey, Ex), np.sqrt(Ex**2+Ey**2), 'bo')
ax.grid(True)
plt.show()

Ex = np.sin(omega_t)
Ey = np.sin(omega_t+0.6)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(np.arctan2(Ey, Ex), np.sqrt(Ex**2+Ey**2), 'bo')
ax.grid(True)
plt.show()


theta = np.arange(0, 2+1/16, 1/16) * np.pi

Ex = np.sin(omega_t)
Ey = np.sin(omega_t+0.6)
arg = np.arctan2(Ey, Ex)
amp = np.sqrt(Ex**2 + Ey**2)
I = np.array([1/(2*np.pi)*cumtrapz((amp*np.cos(arg-pos))**2, omega_t)[-1] for pos in theta])

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, I, 'bo')
ax.grid(True)
plt.show()

Ex = np.sin(omega_t)
Ey = np.sin(omega_t+ np.pi/2 + 0.1)
arg = np.arctan2(Ey, Ex)
amp = np.sqrt(Ex**2 + Ey**2)
I = np.array([1/(2*np.pi)*cumtrapz((amp*np.cos(arg-pos))**2, omega_t)[-1] for pos in theta])

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, I, 'bo')
ax.grid(True)
plt.show()