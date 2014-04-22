import numpy as np
import matplotlib as plt
import sys, string, os

n_archivos = len(sys.argv)

for i in range(n_archivos)
	datos = np.loadtxt(sys.argv[i])

	N = np.shape(datos)[0]

	for i in range(N):
    	plt.scatter(datos[i,0],datos[i,1])
	plt.show()
