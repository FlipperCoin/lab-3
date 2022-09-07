import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

L = 62e-3
modes_lambd = L * (2 / np.arange(30))
v = 1
x = np.linspace(0, L, 1000)
t = np.linspace(0, 5, 1000)
X, T, M_LAMBD = np.meshgrid(x, t, modes_lambd, indexing='ij')
modes_phi = np.sin(2*np.pi/M_LAMBD * X)*np.cos(v/M_LAMBD * T)

sup = np.sum(modes_phi, axis=2)

sup_t0 = sup[:, 0]

plt.plot(x, modes_phi[:, 0, 0], label='mode 0')
plt.plot(x, modes_phi[:, 0, 1], label='mode 1')
plt.plot(x, modes_phi[:, 0, 2], label='mode 2')
plt.plot(x, modes_phi[:, 0, 3], label='mode 3')
plt.xlabel('x')
plt.ylabel('phi')
plt.grid()
plt.legend()
plt.show()

fig, ax = plt.subplots()
line, = ax.plot(x, np.abs(sup[:, 0]))
def animate(i):
    line.set_ydata(np.abs(sup[:, i]))

ani = animation.FuncAnimation(fig, animate, interval=20)
ani.save('sup.mp4')

fig, ax = plt.subplots()
line1, = ax.plot(x, modes_phi[:, 0, 1], label='mode 1')
line2, = ax.plot(x, modes_phi[:, 0, 2], label='mode 2')
line3, = ax.plot(x, modes_phi[:, 0, 3], label='mode 3')
line4, = ax.plot(x, modes_phi[:, 0, 4], label='mode 4')
def animate2(i):
    line1.set_ydata(modes_phi[:, i, 1])
    line2.set_ydata(modes_phi[:, i, 2])
    line3.set_ydata(modes_phi[:, i, 3])
    line4.set_ydata(modes_phi[:, i, 4])


ani = animation.FuncAnimation(fig, animate2, interval=20)
ani.save('modes.mp4')
