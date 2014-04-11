/*
* Autores:
* David Andres Davila
* Juan Sebastián Calderón
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double *PosicionesRadio (int r);
double *VelocidadesRadio (int r);

int main (int argc, char *argv[])
{
/*
*Funcion principal del programa, llama a las otras funciones. Genera un archivo con los datos de posiciones y velocidades de la galaxia
*Entradas: -argc: tamanio del vector de parametros
*-*argv[]: es el vector de parametros, con las posiciones y velocidades inicales del centro de la galaxia
*Salidas: -0
*/
	FILE *archivo;
	int radio;
	double *Posi;
	double *Velo;
	Posi = PosicionesRadio(radio);
	Velo = VelocidadesRadio(radio);
	return 0;
}

double *PosicionesRadio (int r)
{
/*
*Funcion que calcula las condiciones iniciales de posicion de los cuerpos que orbitan el centro de masa de la galaxia a un radio dado
*Entradas: -r: entero que indica el radio en kiloparsecs de la orbita a calcular
*Salidas: -*Posiciones: vector con las coordenadas en x y en y (asumiendo que el centro de la galaxia esta en el origen) de los puntos de la orbita (para el punto i su coordenada en x es la entrada 3*i, su coordenada en y es la entrada 3*i + 1 y su coordenada en z es la entrada 3*i+2)
*/
	double *Posiciones;
	int puntos;
	int total;
	/*este es un if que calcula cuantas masas hay en la orbita, lo cual depende unicamente del radio*/
	if(r == 10)
	{
		puntos = 10;
	}
	else
	{
		puntos = r;
	}
	total = 3*puntos;
	Posiciones = malloc(total * sizeof(int));
	return Posiciones;

}

double *VelocidadesRadio (int r)
{
/*
*Funcion que calcula las condiciones iniciales de velocidad de los cuerpos que orbitan el centro de masa de la galaxia a un radio dado
*Entradas: -r: entero que indica el radio en kiloparsecs de la orbita a calcular
*Salidas: -*Velocidades: vector con las velocidades en x y en y de los puntos de la orbita (para el punto i su velocidad en x es la entrada 3*i, su velocidad en y es la entrada 3*i + 1 y su velocidad en z es la entrada 3*i+2)
*/
	double *Velocidades;
	int puntos;
	int total;
	/*este es un if que calcula cuantas masas hay en la orbita, lo cual depende unicamente del radio*/
	if(r == 10)
	{
		puntos = 10;
	}
	else
	{
		puntos = r;
	}
	total = 3*puntos;
	Velocidades = malloc(total * sizeof(int));
	return Velocidades;

}
