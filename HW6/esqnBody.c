/*
* Autores:
* David Andres Davila
* Juan Sebastian Calderon
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double *RungeKutta (double pix, double piy, double piz, double vix, double viy, double viz);
double GenerarPosicion();

int main (int argc, char *argv[])
{
/*
*Funcion principal del programa. llama a las otras funciones
*Entradas: -argc: tamanio del vector de parametros
*-*argv[]: es el vector de parametros, con el archivo de posicions generado en el primer programa
*Salidas: -0
*/
  int i;
  double *posiciones;
  double *velocidades;
posiciones = malloc(3 * sizeof(int));
velocidades = malloc(3 * sizeof(int));
  for(i = 0; i < 10000; i++){
    posiciones[i] = GenerarPosicion();
}
  RungeKutta(0.0,0.0,0.0,0.0,0.0,0.0);

	return 0;
}


double *RungeKutta (double pix, double piy, double piz, double vix, double viy, double viz)
{
/*
*funcion para resolver la ecuacion diferencial de segundo orden para una particula orbitando el centro de la esfera
*Entradas: -pix: posicion inicial en x (condicion inicial)
*-piy: posicion inicial en y (condicion inicial)
*-piz: posicion inicial en z (condicion inicial)
*-vix: velocidad inicial en x (condicion inicial)
*-viy: velocidad inicial en y (condicion inicial
*-viz: velocidad inicial en z (condicion inicial)
*Salidas: -*SolPos: vector con las posiciones de los puntos segun el tiempo
*-*SolVel: vector con las velocidades en x y en y de los puntosde la orbita 
*/
	double *SolPos;
	double *SolVel;
	SolPos = malloc(15 * sizeof(int));
	SolVel = malloc(15 * sizeof(int));
	return SolPos, SolVel;

}
double GenerarPosicion(){
  double posx;
  double posy;
  double posz;
  posx = 0;
  posy = 0;
  posz = 0;
  return posx,posy,posz;
}
