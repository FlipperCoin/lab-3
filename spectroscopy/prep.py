import numpy as np
from uncertainties import ufloat
from uncertainties.umath import sin

# 1

wave_len = 500e-9
slit_dens = 600  # slits per mm
d = 1e-3 / slit_dens

# d*sin(theta)=n*lambda
# sin(theta)=(n*lambda/d)
# arcsin(n*lambda/d)=theta_n

max_n = int(np.floor(d/wave_len))
theta_n_rad = np.array([np.arcsin(n*wave_len/d) for n in range(0, max_n+1)])

theta_n = 360*theta_n_rad/(2*np.pi)

# 2

theta_err = 5
theta = ufloat(20, theta_err)
n = 1

wave_len = d*sin(theta)/n

wave_len_err_manual = np.abs(theta_err*d*np.cos(theta.nominal_value)/n)

# 3

