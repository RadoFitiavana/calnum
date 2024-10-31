import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = []
y = []
z = []

with open("solution.txt") as f :
	lines = f.readlines()
	for line in lines :
		container = line.split(" ")
		x.append(float(container[0]))
		y.append(float(container[1]))
		z.append(float(container[2]))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z)
plt.show() 
fig.savefig("img.png")