import numpy as np
import matplotlib.pyplot as plt

# data to plot
x_FE = np.loadtxt('x_FE.txt')
y_FE = np.loadtxt('y_FE.txt')
z_FE = np.loadtxt('z_FE.txt')

x_rk4 = np.loadtxt('x_RK4.txt')
y_rk4 = np.loadtxt('y_RK4.txt')
z_rk4 = np.loadtxt('z_RK4.txt')

x2_rk4 = np.loadtxt('x2_RK4.txt')
y2_rk4 = np.loadtxt('y2_RK4.txt')
z2_rk4  = np.loadtxt('z2_RK4.txt')

vx_1 = np.loadtxt('v_x_RK4.txt')
vy_1 = np.loadtxt('v_y_RK4.txt')
vz_1 = np.loadtxt('v_z_RK4.txt')

vx_2 = np.loadtxt('v_x2_RK4.txt')
vy_2 = np.loadtxt('v_y2_RK4.txt')
vz_2 = np.loadtxt('v_z2_RK4.txt')

x_no = np.loadtxt('x_RK4_no.txt')
y_no = np.loadtxt('y_RK4_no.txt')
z_no = np.loadtxt('z_RK4_no.txt')

vx_no = np.loadtxt('v_x_RK4_no.txt')
vy_no = np.loadtxt('v_y_RK4_no.txt')
vz_no = np.loadtxt('v_z_RK4_no.txt')

x2_no = np.loadtxt('x2_RK4_no.txt')
y2_no = np.loadtxt('y2_RK4_no.txt')
z2_no = np.loadtxt('z2_RK4_no.txt')

vx2_no = np.loadtxt('v_x2_RK4_no.txt')
vy2_no = np.loadtxt('v_y2_RK4_no.txt')
vz2_no = np.loadtxt('v_z2_RK4_no.txt')

xf = np.loadtxt('xf.txt')
yf = np.loadtxt('yf.txt')
zf = np.loadtxt('zf.txt')


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

circle1 = plt.Circle((0, 0), R_plus, color='k', fill=False, label=r"$|R_+ - R_-$|")
circle2 = plt.Circle((0, 0), R_minus, color='g', fill=False,  label=r"$R_+ + R_-$")
fig, ax = plt.subplots()

ax.plot(x_FE, y_FE, label='FE')
ax.plot(x, y, label='Analytical')
ax.add_patch(circle1)
ax.add_patch(circle2)
ax.axis('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend()

# plt.figure()
# plt.plot(time, z_rk4, label='RK4')
# plt.xlabel('time[micro seconds]')
# plt.ylabel('z [micro meter]')
# plt.plot(time, z, label='Analytical')
# plt.axis('equal')
# plt.legend()

plt.figure()
plt.plot(time, z_FE, label='FE')
plt.plot(time, z, label='Analytical')
plt.axis('equal')
plt.legend()

# plt.figure()
# plt.plot(time, x_rk4, label='RK4')
# plt.plot(time, x, label='Analytical')
# plt.axis('equal')
# plt.legend()

# plt.figure()
# plt.plot(time, y_rk4, label='RK4')
# plt.plot(time, y, label='Analytical')
# plt.axis('equal')
# plt.legend()


########## 2 particle interactions ##########
plt.figure()
plt.title('2 particle interactions')
plt.scatter(x_rk4, y_rk4, s=.2, label='particle 1')
plt.scatter(y_rk4, y2_rk4,  s=.2, label='particle 2')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.axis('equal')

plt.figure()
plt.title('2 particle no interactions')
plt.scatter(y2_no, y_no,  s=.2, label='x1')
plt.scatter(x2_no, y2_no, s=.2,  label='x2')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.axis('equal')


########## Phase space plots ##########
plt.figure()
plt.title('Phase space plot interaction')
plt.scatter(x_rk4, vx_1,  s=.2, label='particle 1')
plt.scatter(x2_rk4, vx_2,  s=.2, label='particle 2')
plt.xlabel('x')
plt.ylabel('vx')
plt.legend()
plt.axis('equal')

plt.figure()
plt.title('Phase space plot no interaction')
plt.scatter(x_no, vx_no, s=.2,  label='particle 1')
plt.scatter(x2_no, vx2_no, s=.2,  label='particle 2')
plt.xlabel('x')
plt.ylabel('vx')
plt.legend()
plt.axis('equal')

plt.figure()
plt.title('Phase space plot interaction z')
plt.scatter(z_rk4, vz_1,  s=.2, label='particle 1')
plt.scatter(z2_rk4, vz_2,  s=.2, label='particle 2')
plt.xlabel('z')
plt.ylabel('vz')
plt.legend()
plt.axis('equal')

plt.figure()
plt.title('Phase space plot no interaction z')
plt.scatter(z_no, vz_no,  s=.2, label='particle 1')
plt.scatter(z2_no, vz2_no, s=.2,  label='particle 2')
plt.xlabel('z')
plt.ylabel('vz')
plt.legend()
plt.axis('equal')


### 3D plots ###

fig = plt.figure()
fig.suptitle('3D plot interaction')
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_rk4, y_rk4, z_rk4, s=.2,  label='particle 1')
ax.scatter(x2_rk4, y2_rk4, z2_rk4,  s=.2, label='particle 2')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.legend()

fig = plt.figure()
fig.suptitle('3D plot no interaction')
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_no, y_no, z_no, s=.2,  label='particle 1')
ax.scatter(x2_no, y2_no, z2_no, s=.2,  label='particle 2')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.legend()


########## Frequency plots ##########

plt.figure()
plt.title('Frequency plot interaction')
plt.plot(xf, yf)
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.show()

