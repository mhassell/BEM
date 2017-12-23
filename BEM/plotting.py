import numpy as np
from numpy import genfromtxt
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

my_data = genfromtxt('solution.csv', delimiter=',')
X = my_data[:,0]
Y = my_data[:,1]
inside = my_data[:,2]
Z = my_data[:,3]

X = X.reshape([100, 100]);
Y = Y.reshape([100, 100]);
inside = inside.reshape([100, 100]);
Z = Z.reshape([100, 100]);
Z[inside==0]=0

fig = plt.figure(figsize=(20,10))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, linewidth=0, antialiased=False, cmap=cm.coolwarm)
ax.set_zlim(-0.5, 1);
ax.set_xlim(0,1);
ax.set_ylim(0,1);
angle = 220
ax.view_init(30, angle)
plt.show()
