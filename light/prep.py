import numpy as np

# 3

a = 100e-6
lamb = 5000e-10
L = 50e-2

prim_width = 2 * lamb * L / a

# 5
# F = p * (alpha*L)/(2*kB*T*lamb0)
# n = 1 + p * alpha/(2*kB*T)
def cels2kelv(c): return c + 273.15

K = 5
lamb0 = lamb
n = 1 + 76 * K * lamb0 * cels2kelv(27) / L / cels2kelv(0)



