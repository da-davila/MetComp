import numpy as np
import matplotlib.pyplot as plt
#obtenemos los datos

data=loadtxt("volterra.dat")
t=data[:,0]
x=data[:,1]
y=data[:,2]

#Gráfica de X vs. T:

plt.plot(t,x)
plt.scatter(t,x)

#Gráfica de Y vs. T:

plt.plot(t,y)
plt.scatter(t,y)

#Gráfica de X vs. Y:

plt.plot(y,x)
plt.scatter(y,x)
