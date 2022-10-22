import numpy as np
import matplotlib.pyplot as plt

# data to plot
x_FE = np.loadtxt('x_FE.txt')
y_FE = np.loadtxt('y_FE.txt')
z_FE = np.loadtxt('z_FE.txt')

x_rk4 = np.loadtxt('x_RK4.txt')
y_rk4 = np.loadtxt('y_RK4.txt')
z_rk4 = np.loadtxt('z_RK4.txt')

time = np.linspace(0, 50, 50000)

# Analytical solution
x0 = 20
z0 = 20
vy0 = 25

B0 = 9.65e1
V0 = 2.41e6
m = 40.08
d = 500
q = 1

w_z = np.sqrt(2 * q * V0 / (m * d**2))
w_0 = B0*q / m
w_plus = (w_0 + np.sqrt(w_0**2 - 2*w_z**2)) / 2
w_minus = (w_0 - np.sqrt(w_0**2 - 2*w_z**2)) / 2

A_plus = (vy0 + w_minus * x0) / (w_minus - w_plus)
A_minus = - (vy0 + w_plus * x0) / (w_minus - w_plus)

R_plus = np.abs(A_plus + A_minus)
R_minus = np.abs(A_plus - A_minus)

x = A_plus * np.cos(-w_plus * time) + A_minus * np.cos(-w_minus * time)
y = A_plus * np.sin(-w_plus * time) + A_minus * np.sin(-w_minus * time)
z = z0 * np.cos(w_z * time)

plt.figure()
plt.plot(x_rk4, y_rk4, label='RK4')
plt.plot(x, y, label='Analytical')
plt.axis('equal')
plt.legend()

plt.figure()
plt.plot(time, z_rk4, label='RK4')
plt.plot(time, z, label='Analytical')
plt.axis('equal')
plt.legend()

plt.figure()
plt.plot(time, z_FE, label='FE')
plt.plot(time, z, label='Analytical')
plt.axis('equal')
plt.legend()

plt.figure()
plt.plot(time, x_rk4, label='RK4')
plt.plot(time, x, label='Analytical')
plt.axis('equal')
plt.legend()

plt.figure()
plt.plot(time, y_rk4, label='RK4')
plt.plot(time, y, label='Analytical')
plt.axis('equal')
plt.legend()

plt.show()

