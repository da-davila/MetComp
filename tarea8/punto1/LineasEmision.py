import numpy as np
import pylab as plt
from math import *

obs_data = np.loadtxt("energy_counts.dat")
num_energy = obs_data[:,0]
n_counts = obs_data[:,1]

#Se toma el likelihood como chi cuadrado, pues si se toma la exponencial el resultado sera siempre 0.0
def Likelihood (n_counts, n_model):
    chi_squared = sum((n_counts-n_model)**2)
    return chi_squared

def model(num_energy, A, B, E_zero, sigma, alpha):
    return A*(num_energy**alpha) + B*np.exp(-((num_energy - E_zero)/(np.sqrt(2)*sigma))**2)

A_walk = np.empty((0))
B_walk = np.empty((0))
E_zero_walk = np.empty((0))
sigma_walk = np.empty((0))
alpha_walk = np.empty((0))
like_walk = np.empty((0))

#Se corrio el programa varias veces para encontrar el orden de magnitud de cada variable, reemplazando el valor inicial por el resultado del paso anterior, de acuerdo a esto, sera necesario dar pasos de diferente longitud para cada una.
A_walk = np.append(A_walk,10**16)
B_walk = np.append(B_walk,1000)
E_zero_walk = np.append(E_zero_walk,1387.65)
alpha_walk = np.append(alpha_walk, 100)
sigma_walk = np.append(sigma_walk,0)

n_steps = 100000

like_initial = Likelihood (n_counts,model(num_energy,A_walk[0], B_walk[0], E_zero_walk[0], alpha_walk[0], sigma_walk[0]))
like_walk = np.append(like_walk, like_initial)

for i in range(n_steps - 1):
    #como cada parametro tiene diferente orden de magnitud, cada uno tiene una longitud de paso diferente
    A_prime = np.random.normal(A_walk[i], 10**13)
    B_prime = np.random.normal(B_walk[i], 1)
    E_zero_prime = np.random.normal(E_zero_walk[i], 1)
    alpha_prime = np.random.normal(alpha_walk[i], 0.1)
    sigma_prime = np.random.normal(sigma_walk[i], 0.01)

    like_actual = Likelihood (n_counts,model(num_energy,A_walk[i], B_walk[i], E_zero_walk[i], alpha_walk[i], sigma_walk[i]))
    like_prime = Likelihood (n_counts,model(num_energy,A_prime, B_prime, E_zero_prime, alpha_prime, sigma_prime))

    alpha = like_actual - like_prime

    if(alpha>0.0):
        A_walk = np.append(A_walk, A_prime)
        B_walk = np.append(B_walk, B_prime)
        E_zero_walk = np.append(E_zero_walk, E_zero_prime)
        alpha_walk = np.append(alpha_walk, alpha_prime)
        sigma_walk = np.append(sigma_walk,sigma_prime)
        like_walk = np.append(like_walk,like_prime)
    else:
        beta = np.random.random()
        if(np.log(beta)<=alpha):
            A_walk = np.append(A_walk, A_prime)
            B_walk = np.append(B_walk, B_prime)
            E_zero_walk = np.append(E_zero_walk, E_zero_prime)
            alpha_walk = np.append(alpha_walk, alpha_prime)
            sigma_walk = np.append(sigma_walk,sigma_prime)
            like_walk = np.append(like_walk,like_prime)
        else:
            A_walk = np.append(A_walk, A_walk[i])
            B_walk = np.append(B_walk, B_walk[i])
            E_zero_walk = np.append(E_zero_walk, E_zero_walk[i])
            alpha_walk = np.append(alpha_walk, alpha_walk[i])
            sigma_walk = np.append(sigma_walk,sigma_walk[i])
            like_walk = np.append(like_walk,like_actual)

max_like_index = np.argmin(like_walk)
best_A = A_walk[max_like_index]
best_B = B_walk[max_like_index]
best_E_zero = E_zero_walk[max_like_index]
best_alpha = alpha_walk[max_like_index]
best_sigma = sigma_walk[max_like_index]


print "El valor de A es", best_A
print "El valor de B es", best_B
print "El valor de E_0 es", best_E_zero
print "El valor de alpha es", best_alpha
print "El valor de sigma es", best_sigma




plt.scatter(num_energy,n_counts)
plt.plot(num_energy,model(num_energy, best_A, best_B, best_E_zero, best_alpha, best_sigma))
plt.show()

count, bins, ignored =plt.hist(A_walk, 20, normed=True)
plt.show()
count, bins, ignored =plt.hist(B_walk, 20, normed=True)
plt.show()
count, bins, ignored =plt.hist(E_zero_walk, 20, normed=True)
plt.show()
count, bins, ignored =plt.hist(alpha_walk, 20, normed=True)
plt.show()
count, bins, ignored =plt.hist(sigma_walk, 20, normed=True)
plt.show()

