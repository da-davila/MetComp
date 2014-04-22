/*
* Autores:
* David Andres Davila
* Juan Sebastián Calderón
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.141592653589793
#define Gra 39.486 //units of ua+3 msun-1 yr-1

int main (int argc, char *argv[])
{
  /*declaramos las variables necesarias*/
  float x_0;
  float y_0;
  float z_0;
  float vx_0;
  float vy_0;
  float vz_0;
  float M;
  int N;
  float R;
  int i;
  float x_temp;
  float y_temp;
  float vx_temp;
  float vy_temp;
  float RadTemp;
  float AngTemp;

 /*iniciamos los valores conocidos*/
  x_0 = atof(argv[1]);
  y_0 = atof(argv[2]);
  z_0 = atof(argv[3]);
  vx_0 = atof(argv[4]);
  vy_0 = atof(argv[5]);
  vz_0 = atof(argv[6]);
  M = atof(argv[7]);
  N = atoi(argv[8]);
  R = atoi(argv[9]);

  /*escribimos las condiciones iniciales del centro de masa en pantalla*/
  printf( "%i %f %f %f %f %f %f \n", -1, x_0, y_0, z_0, vx_0,vy_0,vz_0);

  /*es necesario un for para ubicar las N particulas, se ubicaran aleatoriamente en el interior de un circulo de radio R*/
  for(i = 0; i<N; i++){
    RadTemp = drand48()*R;
    AngTemp= drand48()*2*PI;
    x_temp = sqrt(RadTemp)*cos(AngTemp) + x_0;
    y_temp = sqrt(RadTemp)*sin(AngTemp) + y_0;
    vx_temp = sqrt(Gra*M/RadTemp)*sin(AngTemp);
    vy_temp = sqrt(Gra*M/RadTemp)*sin(AngTemp);

    printf("%i %f %f %f %f %f %f \n", i, x_temp, y_temp, z_0, vx_temp,vy_temp,vz_0);

  }
	
  return 0;
}

