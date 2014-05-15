import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


delta_t=0.0001
#Con este deltat se demora como doce minutos. Debe ser tal que 0.001 sea multiplo entero de deltat 
time=0.800
time_iterations=int(time/delta_t)
datos=np.loadtxt("./lotka_volterra_obs.dat")
t_obs=datos[:,0]
x_obs=datos[:,1]
y_obs=datos[:,2]
x_0=x_obs[0]
y_0=y_obs[0]
t=np.linspace(0,0.799,time_iterations)

#Implemento el metodo de Euler para resolver la ecuacion diferencial
def euler(alfa,beta,gamma,delta,x_0,y_0,iterations,t_obs,delta_t):
  x=np.zeros(iterations)
  y=np.zeros(iterations)
  k= int(t_obs[0]*(1.0/delta_t))
  x[k]=x_0
  y[k]=y_0 
  #Solo la resuelvo desde el tiempo minimo en t_obs ya que los valores iniciales de x & y seran los valores observados en este tiempo
  for i in range(k+1,iterations):
     x[i]=x[i-1]+delta_t*(alfa*x[i-1]-beta*x[i-1]*y[i-1])
     y[i]=y[i-1]+delta_t*(-gamma*y[i-1]+delta*x[i-1]*y[i-1])
  return x,y  


#funcion para el likelihood. Solo uso los valores en X & Y en los tiempos donde hubo observaciones
def likelihood(x_obs,y_obs,x,y,t_obs,delta_t):
  x_model=np.zeros(len(x_obs))
  y_model=np.zeros(len(x_obs))
  for i in range(len(x_obs)):
     k=int(t_obs[i]*(1.0/delta_t))
     x_model[i]=x[k]
     y_model[i]=y[k]
  chi_squared=-0.5*(sum((x_obs-x_model)**2)+sum((y_obs-y_model)**2))
  return chi_squared

#Inicializo la caminata para las 4 variables
alfa_walk=np.empty((0))
beta_walk=np.empty((0))
gamma_walk=np.empty((0))
delta_walk=np.empty((0))
l_walk=np.empty((0))

#Inicializo las variables con numeros calculados a simple vista despues de correr varias veces el programa (para que no haya overflow)
alfa_walk=np.append(alfa_walk, 32)
beta_walk=np.append(beta_walk, 4)
gamma_walk=np.append(gamma_walk, 45)
delta_walk=np.append(delta_walk, 1.6)

x_init,y_init=euler(alfa_walk[0],beta_walk[0],gamma_walk[0],delta_walk[0],x_0,y_0,time_iterations,t_obs,delta_t)

l_walk=np.append(l_walk,likelihood(x_obs,y_obs,x_init,y_init,t_obs,delta_t))
print l_walk

#Implemento el metodo de Metropolis hastings para calcular las variables
n_iterations=15000
for i in range(n_iterations):
  alfa_prime=np.random.normal(alfa_walk[i],0.2)
  beta_prime=np.random.normal(beta_walk[i],0.2)
  gamma_prime=np.random.normal(gamma_walk[i],0.2)
  delta_prime=np.random.normal(delta_walk[i],0.2)
  x_init,y_init=euler(alfa_walk[i],beta_walk[i],gamma_walk[i],delta_walk[i],x_0,y_0,time_iterations,t_obs,delta_t)
  x_prime,y_prime=euler(alfa_prime,beta_prime,gamma_prime,delta_prime,x_0,y_0,time_iterations,t_obs,delta_t)
  l_prime=likelihood(x_obs,y_obs,x_prime,y_prime,t_obs,delta_t)
  l_init=likelihood(x_obs,y_obs,x_init,y_init,t_obs,delta_t)

  alpha = l_prime-l_init
  #Hago las modificaciones correspondientes para el alpha y b con respecto al ejemplo de Jaime, ya que no uso la funcion exponencial para el likelihood
  if(alpha>=0):
    alfa_walk=np.append(alfa_walk,alfa_prime)
    beta_walk=np.append(beta_walk,beta_prime)
    gamma_walk=np.append(gamma_walk,gamma_prime)
    delta_walk=np.append(delta_walk,delta_prime)
    l_walk=np.append(l_walk,l_prime)
  else:
      b=np.random.random()
      if(b<=np.exp(alpha)):
         alfa_walk=np.append(alfa_walk,alfa_prime)
         beta_walk=np.append(beta_walk,beta_prime)
         gamma_walk=np.append(gamma_walk,gamma_prime)
         delta_walk=np.append(delta_walk,delta_prime)
         l_walk=np.append(l_walk,l_prime)
      else:
         alfa_walk=np.append(alfa_walk,alfa_walk[i])
         beta_walk=np.append(beta_walk,beta_walk[i])
         gamma_walk=np.append(gamma_walk,gamma_walk[i])
         delta_walk=np.append(delta_walk,delta_walk[i])
         l_walk=np.append(l_walk,l_init)
  

#EStimacion de los valores mas probables de las variables
max_likelihood_id=np.argmax(l_walk)
best_alfa=alfa_walk[max_likelihood_id]
best_beta=beta_walk[max_likelihood_id]
best_gamma=gamma_walk[max_likelihood_id]
best_delta=delta_walk[max_likelihood_id]  

print "Los valores mas probables de alfa, beta, gamma y delta son:",best_alfa, best_beta, best_gamma, best_delta,"respectivamente."

#Calculo de incertidumbres

dalfa=best_alfa/20
dbeta=best_beta/20
dgamma=best_gamma/20
ddelta=best_delta/20

index_alfa = np.where((np.abs(beta_walk-best_beta)<dbeta) & (np.abs(gamma_walk-best_gamma)<dgamma) & (np.abs(delta_walk-best_delta)<ddelta))
selected_alfas=alfa_walk[index_alfa] - best_alfa
lalfa=(-1)*l_walk[index_alfa]
inc_alfa= np.sqrt(((best_alfa-selected_alfas)**2)/(2*lalfa**2))

index_beta = np.where((np.abs(alfa_walk-best_alfa)<dalfa) & (np.abs(gamma_walk-best_gamma)<dgamma) & (np.abs(delta_walk-best_delta)<ddelta))
selected_betas=beta_walk[index_beta] - best_beta
lbeta=(-1)*l_walk[index_beta]
inc_beta= np.sqrt(((best_beta-selected_betas)**2)/(2*lbeta**2))

index_gamma = np.where((np.abs(beta_walk-best_beta)<dbeta) & (np.abs(alfa_walk-best_alfa)<dalfa) & (np.abs(delta_walk-best_delta)<ddelta))  
selected_gammas=gamma_walk[index_gamma] - best_gamma
lgamma=(-1)*l_walk[index_gamma]
inc_gamma= np.sqrt(((best_gamma-selected_gammas)**2)/(2*lgamma**2))

index_delta = np.where((np.abs(beta_walk-best_beta)<dbeta) & (np.abs(gamma_walk-best_gamma)<dgamma) & (np.abs(alfa_walk-best_alfa)<dalfa))
selected_deltas=delta_walk[index_delta] - best_delta
ldelta=(-1)*l_walk[index_delta]
inc_delta= np.sqrt(((best_delta-selected_deltas)**2)/(2*ldelta**2))

print "Las incertidumbres en alfa, beta, gamma, delta son:", np.amax(inc_alfa), np.amax(inc_beta), np.amax(inc_gamma), np.amax(inc_delta), "respectivamente."



#Graficas bidimensionales de todas las combinaciones entre las variables

min_alfa = np.amin(alfa_walk)
min_beta = np.amin(beta_walk)
min_gamma = np.amin(gamma_walk)
min_delta = np.amin(delta_walk)

max_alfa = np.amax(alfa_walk)
max_beta = np.amax(beta_walk)
max_gamma = np.amax(gamma_walk)
max_delta = np.amax(delta_walk)

n_points=len(alfa_walk)

#Histograma de Alfa vs Beta
grid_alfa, grid_beta = np.mgrid[min_alfa:max_alfa:200j, min_beta:max_beta:200j]
points = np.ones((n_points,2))
points[:,0] = alfa_walk
points[:,1] = beta_walk
grid_l = griddata(points, -l_walk, (grid_alfa, grid_beta), method='cubic')
fig = plt.figure()
plt.imshow(grid_l.T, extent=(min_alfa, max_alfa, min_beta, max_beta), aspect='auto',origin='lower')
plt.savefig('alfa_vs_beta.png')
plt.close()

#Histograma de Alfa vs Gamma
grid_alfa, grid_gamma = np.mgrid[min_alfa:max_alfa:200j, min_gamma:max_gamma:200j]
points = np.ones((n_points,2))
points[:,0] = alfa_walk
points[:,1] = gamma_walk
grid_l = griddata(points, -l_walk, (grid_alfa, grid_gamma), method='cubic')
fig = plt.figure()
plt.imshow(grid_l.T, extent=(min_alfa, max_alfa, min_gamma, max_gamma), aspect='auto',origin='lower')
plt.savefig('alfa_vs_gamma.png')
plt.close()

#Histograma de Alfa vs Delta
grid_alfa, grid_delta = np.mgrid[min_alfa:max_alfa:200j, min_delta:max_delta:200j]
points = np.ones((n_points,2))
points[:,0] = alfa_walk
points[:,1] = delta_walk
grid_l = griddata(points, -l_walk, (grid_alfa, grid_delta), method='cubic')
fig = plt.figure()
plt.imshow(grid_l.T, extent=(min_alfa, max_alfa, min_delta, max_delta), aspect='auto',origin='lower')
plt.savefig('alfa_vs_delta.png')
plt.close()

#Histograma de Beta vs Gamma
grid_beta, grid_gamma = np.mgrid[min_beta:max_beta:200j, min_gamma:max_gamma:200j]
points = np.ones((n_points,2))
points[:,0] = beta_walk
points[:,1] = gamma_walk
grid_l = griddata(points, -l_walk, (grid_beta, grid_gamma), method='cubic')
fig = plt.figure()
plt.imshow(grid_l.T, extent=(min_beta, max_beta, min_gamma, max_gamma), aspect='auto',origin='lower')
plt.savefig('beta_vs_gamma.png')
plt.close()

#Histograma de Beta vs Delta
grid_beta, grid_delta = np.mgrid[min_beta:max_beta:200j, min_delta:max_delta:200j]
points = np.ones((n_points,2))
points[:,0] = beta_walk
points[:,1] = delta_walk
grid_l = griddata(points, -l_walk, (grid_beta, grid_delta), method='cubic')
fig = plt.figure()
plt.imshow(grid_l.T, extent=(min_beta, max_beta, min_delta, max_delta), aspect='auto',origin='lower')
plt.savefig('beta_vs_delta.png')
plt.close()


#Histograma de Gamma vs Delta
grid_gamma, grid_delta = np.mgrid[min_gamma:max_gamma:200j, min_delta:max_delta:200j]
points = np.ones((n_points,2))
points[:,0] = gamma_walk
points[:,1] = delta_walk
grid_l = griddata(points, -l_walk, (grid_gamma, grid_delta), method='cubic')
fig = plt.figure()
plt.imshow(grid_l.T, extent=(min_gamma, max_gamma, min_delta, max_delta), aspect='auto',origin='lower')
plt.savefig('gamma_vs_delta.png')
plt.close()

best_x,best_y=euler(best_alfa,best_beta,best_gamma,best_delta,x_0,y_0,time_iterations,t_obs,delta_t)
plt.scatter(x_obs,y_obs)
plt.plot(best_x,best_y)
plt.savefig('Datos observados vs modelo')
plt.close()
