#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*Este código permite resolver la ecuación de Volterra para el modelamiento de cazador y presa
  Tenemos las ecuaciones dx/dt=Ax-Bxy y dy/dt=-Cy+Dxy*/

/*PARTE PRIMERA: Hallar condiciones iniciales para condiciones de equilibrio asumiento que A=20,B=1,C=30,D=1
Para las condiciones de equilibrio se pide que dx/dt=dy/dt=0 con lo cual se obtiene que:
x(0)=30 y y(0)=20*/

/*PARTE SEGUNDA:Resolver las ecuaciones diferenciales de Volterra manteniendo fijo y(0)=20 y variando x(0)=N con N=1,2,3,...,29,30 con t entre 0 y 1
Lo vamos a resolver por el método de Euler
Definimos una separación de 0.01 e inicializamos las variables*/


float xfunction_prime(float t, float x, float y);
float yfunction_prime(float t, float x, float y);
void imprimir(float *x, float *y, float *t);

int main (){
  /*declaramos las variables necesarias*/
  float h;
  float min_t;
  float max_t;
  int n_points;
  float *t;
  float *x;
  float *y;
  int i;
  int j;
  int k;
  float x_prime;
  float y_prime;


  /*inicializamos las variables constantes*/
  h = 0.01;
  min_t = 0.0;
  max_t = 1.0;
  n_points = 100;

  /*separamos memoria para los punteros*/
  if(!(t = malloc(sizeof(float) * n_points))){
    printf("problem with memory allocation");
    exit(1);
  }

  if(!(x = malloc(sizeof(float) * n_points))){
    printf("problem with memory allocation");
    exit(1);
  }
  if(!(y = malloc(sizeof(float) * n_points))){
    printf("problem with memory allocation");
    exit(1);
  }


  /*llenamos los vectores, variando la entrada principal de x entre 30 y 1*/


  t[0] = min_t;
  y[0] = 20.0;
  for(i=0; i <30;i++){
    x[0] = 30 - i;
    for (j = 1; j < n_points; j++){
      x_prime = xfunction_prime(t[j-1], x[j-1], y[j-1]);
      y_prime = yfunction_prime(t[j-1], x[j-1], y[j-1]);
      x[j]=x[j-1]+h*x_prime;
      y[j]=y[j-1]+h*y_prime;
      t[j] = t[j-1]+h;
    }
    imprimir(x,y,t);
  }

}
float xfunction_prime(float t, float x, float y){
  float resp;
  resp = (20*x)-(x*y);
  return resp;
}
float yfunction_prime(float t, float x, float y){
  float resp;
  resp = (x*y)-(30*y);
  return resp;
}
void imprimir(float *x, float *y, float *t){
  int i;
  int puntos;
  puntos = sizeof(x);
    for (i=0; i<puntos; i++) {
        printf("%f %f %f \n", t[i], x[i], y[i]);
    }
}


