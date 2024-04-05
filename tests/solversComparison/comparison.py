import numpy as np
import matplotlib.pyplot as plt

FPI = np.loadtxt('FPI.txt').transpose()
FPIA = np.loadtxt('FPIA.txt').transpose()
Jacobi = np.loadtxt('Jacobi.txt').transpose()
GS = np.loadtxt('GS.txt').transpose()
SGS = np.loadtxt('SGS.txt').transpose()
AGS = np.loadtxt('AGS.txt').transpose()
GradDes = np.loadtxt('GradDes.txt').transpose()

plt.plot(*FPI, label='Fixed-Point Iteration')
plt.plot(*FPIA, label='Accelerated FPI')
plt.plot(*Jacobi, label='Jacobi method')
plt.plot(*GS, label='Gauss-Seidel method')
plt.plot(*SGS, label='Sym Gauss-Seidel method')
plt.plot(*AGS, label='Acc Gauss-Seidel method')
plt.plot(*GradDes, label='Steepest Gradient Descent')

plt.legend()
plt.xscale('linear')
plt.yscale('log')
plt.xlabel('Размер матрицы')
plt.ylabel('Время решения, нс')
plt.grid()

plt.savefig('plot.png')
plt.show()