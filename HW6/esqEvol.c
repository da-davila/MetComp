/*
* Autores:
* David Andres Davila
* Juan Sebastian Calderon
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double *RungeKutta (double pix, double piy, double piz, double vix, double viy, double viz);

int main (int argc, char *argv[])
{
/*
*Funcion principal del programa. Genera cinco archivos de texto con los datos de evolucion en cinco tiempos diferentes.
*Entradas: -argc: tamanio del vector de parametros
*-*argv[]: es el vector de parametros, con el archivo de posicions generado en el primer programa
*Salidas: -0
*/
	FILE *Aentrada;
	FILE *archivo1;
	FILE *archivo2;
	FILE *archivo3;
	FILE *archivo4;
	FILE *archivo5;
	RungeKutta(0.0,0.0,0.0,0.0,0.0,0.0);

	return 0;
}


double *RungeKutta (double pix, double piy,double piz, double vix, double viy, double viz)
{
/*
*Metodo para resolver la ecuacion diferencial de segundo orden para una particula orbitando el centro de la galaxia
*Entradas: -pix: posicion inicial en x (condicion inicial)
*-piy: posicion inicial en y (condicion inicial)
*-piz: posicion inicial en z (condicion inicial)
*-vix: velocidad inicial en x (condicion inicial)
*-viy: velocidad inicial en y (condicion inicial)
*-viz: velocidad inicial en z (condicion inicial)
*Salidas: -*SolPos: vector con las posiciones de los puntos segun el tiempo
*-*SolVel: vector con las velocidades en x en y y en z de los puntosde la orbita 
*/
	double *SolPos;
	double *SolVel;
	SolPos = malloc(15 * sizeof(int));
	SolVel = malloc(15 * sizeof(int));
	return SolPos, SolVel;

}
