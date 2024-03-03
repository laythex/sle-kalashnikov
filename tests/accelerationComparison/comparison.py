import numpy as np
import matplotlib.pyplot as plt

FPI = np.loadtxt('FPI.txt').transpose()
Jacobi = np.loadtxt('Jacobi.txt').transpose()
GS = np.loadtxt('GS.txt').transpose()
FPIA = np.loadtxt('FPIA.txt').transpose()

plt.scatter(*FPI, label='Fixed-Point Iteration')
plt.scatter(*Jacobi, label='Jacobi method')
plt.scatter(*GS, label='Gauss-Seidel method')
plt.scatter(*FPIA, label='Accelerated FPI')

plt.legend()
plt.xscale('linear')
plt.yscale('log')
plt.xlabel('Размер матрицы')
plt.ylabel('Время решения, нс')

plt.savefig('plot.png')
plt.show()