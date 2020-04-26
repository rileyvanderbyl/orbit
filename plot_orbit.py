#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import initpos_to_orbit as orb
from mpl_toolkits import mplot3d

r1 = np.array([1,0,0])
v1 = np.array([0,1,0])
mu = 100
[a,e,i,Omega,omega,M] = orb.initpos_to_orbit(r1,v1,mu)
print(r1,v1)
print([a,e,i,Omega,omega,M])






fig = plt.figure()
ax = plt.axes(projection='3d')

# t = np.linspace(1,10,1000)

dt = 0.01
steps = 10

xline = np.zeros(shape=steps)
yline = np.zeros(shape=steps)
zline = np.zeros(shape=steps)
for t in range(0,steps):

	y = orb.orbit_to_position(a,e,i,Omega,omega,M,0,dt*t,mu)
	outpos = y[0:3]
	xline[t] = outpos[0]
	yline[t] = outpos[1]
	zline[t] = outpos[2]
	print(t, xline[t],yline[t],zline[t])
ax.plot3D(xline, yline, zline, 'black')

plt.show()