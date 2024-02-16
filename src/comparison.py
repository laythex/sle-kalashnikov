import numpy as np
import matplotlib.pyplot as plt

with open('out.txt') as file:
    cols_step, cols_max, density_step = map(float, file.readline().split())

data = np.loadtxt('out.txt', skiprows=1)
data = data.astype('float64')
data = np.transpose(data)
data = np.flip(data, 0)
ax = plt.gca()
im = ax.imshow(np.log(data), extent=[0, cols_max, 0, 1])
ax.set_aspect(cols_max)
ax.set_xlabel('matrix size')
ax.set_ylabel('relative matrix density')
plt.title(r'$\ln\ t_{dense}/t_{sparse}$')
plt.colorbar(im)
plt.savefig('plot.png', dpi=300)
plt.show()