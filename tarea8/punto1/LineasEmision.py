import numpy as np
import pylab as plt
from math import *
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

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
B_walk = np.append(B_walk,1900)
E_zero_walk = np.append(E_zero_walk,1387.65)
alpha_walk = np.append(alpha_walk, 100)
sigma_walk = np.append(sigma_walk,-2)

n_steps = 10000

like_initial = Likelihood (n_counts,model(num_energy,A_walk[0], B_walk[0], E_zero_walk[0], alpha_walk[0], sigma_walk[0]))
like_walk = np.append(like_walk, like_initial)

for i in range(n_steps - 1):
    #como cada parametro tiene diferente orden de magnitud, cada uno tiene una longitud de paso diferente
    A_prime = np.random.normal(A_walk[i], 10**13)
    B_prime = np.random.normal(B_walk[i], 10)
    E_zero_prime = np.random.normal(E_zero_walk[i], 10)
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

dA = best_A/2.0
dB = best_B/2.0
dE_zero= best_E_zero/2.0
dalpha = best_alpha/2.0
dsigma = -best_sigma/2.0

index_A = np.where((np.abs(B_walk-best_B)<dB) & (np.abs(E_zero_walk-best_E_zero)<dE_zero) & (np.abs(alpha_walk-best_alpha)<dalpha) & (np.abs(sigma_walk-best_sigma)<dsigma))
selected_As = A_walk[index_A] - best_A
lA=(-1)*like_walk[index_A]
inc_A= np.sqrt(((best_A-selected_As)**2)/(2*lA**2))


index_B = np.where((np.abs(A_walk-best_A)<dA) & (np.abs(E_zero_walk-best_E_zero)<dE_zero) & (np.abs(alpha_walk-best_alpha)<dalpha) & (np.abs(sigma_walk-best_sigma)<dsigma))
selected_Bs=B_walk[index_B] - best_B
lB=(-1)*like_walk[index_B]
inc_B= np.sqrt(((best_B-selected_Bs)**2)/(2*lB**2))

index_Ezero = np.where((np.abs(B_walk-best_B)<dB) & (np.abs(A_walk-best_A)<dA) & (np.abs(alpha_walk-best_alpha)<dalpha) & (np.abs(sigma_walk-best_sigma)<dsigma))
selected_Ezeros=E_zero_walk[index_Ezero] - best_E_zero
lEzero=(-1)*like_walk[index_Ezero]
inc_Ezero= np.sqrt(((best_E_zero-selected_Ezeros)**2)/(2*lEzero**2))

index_alpha = np.where((np.abs(B_walk-best_B)<dB) & (np.abs(E_zero_walk-best_E_zero)<dE_zero) & (np.abs(A_walk-best_A)<dA) & (np.abs(sigma_walk-best_sigma)<dsigma))
selected_alphas=alpha_walk[index_alpha] - best_alpha
lalpha=(-1)*like_walk[index_alpha]
inc_alpha= np.sqrt(((best_alpha-selected_alphas)**2)/(2*lalpha**2))

index_sigma = np.where((np.abs(B_walk-best_B)<dB) & (np.abs(E_zero_walk-best_E_zero)<dE_zero) & (np.abs(alpha_walk-best_alpha)<dalpha) & (np.abs(A_walk-best_A)<dA))
selected_sigmas=sigma_walk[index_sigma] - best_sigma
lsigma=(-1)*like_walk[index_sigma]
inc_sigma= np.sqrt(((best_sigma-selected_sigmas)**2)/(2*lsigma**2))


print "El valor de A es", best_A
print "El valor de B es", best_B
print "El valor de E_0 es", best_E_zero
print "El valor de alpha es", best_alpha
print "El valor de sigma es", best_sigma
print "Las incertidumbres en A, B, E_0, alpha, sigma son:",np.amax(inc_A), np.amax(inc_B), np.amax(inc_Ezero), np.amax(inc_alpha), np.amax(inc_sigma), "respectivamente."



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


plt.plot(A_walk,B_walk)
plt.xlabel('$A$')
plt.ylabel('$B$')
plt.title('$A$ vs. $B$')
plt.show()
plt.plot(A_walk,E_zero_walk)
plt.xlabel('$A$')
plt.ylabel('$E_0$')
plt.title('$A$ vs. $E_0$')
plt.show()
plt.plot(A_walk,sigma_walk)
plt.xlabel('$A$')
plt.ylabel('$\sigma$')
plt.title('$A$ vs. $\sigma$')
plt.show()
plt.plot(A_walk,alpha_walk)
plt.xlabel('$A$')
plt.ylabel('$\\alpha$')
plt.title('$A$ vs. $\\alpha$')
plt.show()
plt.plot(B_walk,E_zero_walk)
plt.xlabel('$B$')
plt.ylabel('$E_0$')
plt.title('$B$ vs. $E_0$')
plt.show()
plt.plot(B_walk,sigma_walk)
plt.xlabel('$B$')
plt.ylabel('$\sigma$')
plt.title('$B$ vs. $\sigma$')
plt.show()
plt.plot(B_walk,alpha_walk)
plt.xlabel('$B$')
plt.ylabel('$\\alpha$')
plt.title('$B$ vs. $\\alpha$')
plt.show()
plt.plot(E_zero_walk,sigma_walk)
plt.xlabel('$E_0$')
plt.ylabel('$\sigma$')
plt.title('$E0$ vs. $\sigma$')
plt.show()
plt.plot(E_zero_walk,alpha_walk)
plt.xlabel('$E_0$')
plt.ylabel('$\\alpha$')
plt.title('$E_0$ vs. $\\alpha$')
plt.show()
plt.plot(sigma_walk,alpha_walk)
plt.xlabel('$\sigma$')
plt.ylabel('$\\alpha$')
plt.title('$\sigma$ vs. $\\alpha$')
plt.show()




