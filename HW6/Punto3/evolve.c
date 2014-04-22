#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define GRA 39.486 //units of ua+3 msun-1 yr-1

FLOAT * get_memory(int n_points);
FLOAT cinetica(int n_points, FLOAT *v_x, FLOAT *v_y, FLOAT *v_z, FLOAT *masas);
FLOAT potencial(int n_points, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *masas);
void aceleracion(int n_points, int *ID, FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *masas);

int main(int argc, char **argv){
  /*declaramos las variables necesarias*/
  FILE *ArchEntrada;
  char *nombreEntrada;
  FLOAT T;
  int M;
  int n_points;
  int *ID;
  FLOAT *x;
  FLOAT *y;
  FLOAT *z;
  FLOAT *vx;
  FLOAT *vy;
  FLOAT *vz;
  FLOAT h;
  FILE *Arch1;
  FILE *Arch2;
  FILE *Arch3;
  FILE *Arch4;
  FILE *Arch5;
  FILE *energias;
  FLOAT *ax;
  FLOAT *ay;
  FLOAT *az;
  FLOAT *masas;
  int pasos;
  int i;
  int j;
  FLOAT Ecinetica;
  FLOAT Epotencial;

  FLOAT *ax_old;
  FLOAT *ay_old;
  FLOAT *az_old;
  FLOAT *x_temp;
  FLOAT *y_temp;
  FLOAT *z_temp;
  FLOAT *vx_temp;
  FLOAT *vy_temp;
  FLOAT *vz_temp;
  FLOAT *ax_temp;
  FLOAT *ay_temp;
  FLOAT *az_temp;
  FLOAT  *k_1_x;
  FLOAT  *k_1_y;
  FLOAT  *k_1_z;
  FLOAT  *k_1_vx;
  FLOAT  *k_1_vy;
  FLOAT  *k_1_vz;
  FLOAT  *k_2_x;
  FLOAT  *k_2_y;
  FLOAT  *k_2_z;
  FLOAT  *k_2_vx;
  FLOAT  *k_2_vy;
  FLOAT  *k_2_vz;   
  FLOAT  *k_3_x;
  FLOAT  *k_3_y;
  FLOAT  *k_3_z;   
  FLOAT  *k_3_vx;
  FLOAT  *k_3_vy;
  FLOAT  *k_3_vz; 
  FLOAT  *k_4_x;
  FLOAT  *k_4_y;
  FLOAT  *k_4_z;   
  FLOAT  *k_4_vx;
  FLOAT  *k_4_vy;
  FLOAT  *k_4_vz;


  /*iniciamos los valores conocidos(excepto los punteros, para los cuales no hemos separado la memoria)*/
  nombreEntrada = argv[1];
  T = atof(argv[2]);
  M=atoi(argv[3]);
  h = 0.01;
  ArchEntrada = fopen(nombreEntrada,"r");

  if(!ArchEntrada){
    printf("problema abriendo el archivo %s", nombreEntrada);
    exit(1);
  }
  for (i=0; i<n_points; i++) {
    fscanf(ArchEntrada,"%i %f %f %f %f %f %f \n", &(ID[i]),&(x[i]),&(y[i]),&(z[i]),&(vx[i]),&(vy[i]),&(vz[i]));
    }
    
  fclose(ArchEntrada);
  
  pasos = floor(T/h);

  /*separamos memoria para los punteros*/
  x = get_memory(n_points);
  y = get_memory(n_points);
  z = get_memory(n_points);
  vx = get_memory(n_points);
  vy = get_memory(n_points);
  vz = get_memory(n_points);
  ax = get_memory(n_points);
  ay = get_memory(n_points);
  az = get_memory(n_points);
  masas = get_memory(n_points);

  ax_old = get_memory(n_points);
  ay_old = get_memory(n_points);
  az_old = get_memory(n_points);
  x_temp = get_memory(n_points);
  y_temp = get_memory(n_points);
  z_temp = get_memory(n_points);
  vx_temp = get_memory(n_points);
  vy_temp = get_memory(n_points);
  vz_temp = get_memory(n_points);
  ax_temp = get_memory(n_points);
  ay_temp = get_memory(n_points);
  az_temp = get_memory(n_points);
  k_1_x = get_memory(n_points);
  k_1_y = get_memory(n_points);
  k_1_z = get_memory(n_points);
  k_1_vx = get_memory(n_points);
  k_1_vy = get_memory(n_points);
  k_1_vz = get_memory(n_points);
  k_2_x = get_memory(n_points);
  k_2_y = get_memory(n_points);
  k_2_z = get_memory(n_points);
  k_2_vx = get_memory(n_points);
  k_2_vy = get_memory(n_points);
  k_2_vz = get_memory(n_points);
  k_3_x = get_memory(n_points);
  k_3_y = get_memory(n_points);
  k_3_z = get_memory(n_points);
  k_3_vx = get_memory(n_points);
  k_3_vy = get_memory(n_points);
  k_3_vz = get_memory(n_points);
  k_4_x = get_memory(n_points);
  k_4_y = get_memory(n_points);
  k_4_z = get_memory(n_points);
  k_4_vx = get_memory(n_points);
  k_4_vy = get_memory(n_points);
  k_4_vz = get_memory(n_points);

  /*implementacion del Runge-Kutta de cuarto orden*/
  energias = fopen("energias.dat","w");
        
   ax_old = ax;
   ay_old = ay;
   az_old = az;
        
   aceleracion(n_points,ID, ax_old, ay_old, az_old, x, y, z, masas);
   for(j=0;j<n_points;j++){
            
     k_1_x[j] = vx[j];
     k_1_y[j] = vy[j];
     k_1_z[j] = vz[j];
     
     k_1_vx[j] = ax[j];
     k_1_vy[j] = ay[j];
     k_1_vz[j] = az[j];
            
     x_temp[j] = x[j] + (h/2.0)*k_1_x[j];
     y_temp[j] = y[j] + (h/2.0)*k_1_y[j];
     z_temp[j] = z[j] + (h/2.0)*k_1_z[j];
     
     vx_temp[j] = vx[j] + (h/2.0)*k_1_vx[j];
     vy_temp[j] = vy[j] + (h/2.0)*k_1_vy[j];
     vz_temp[j] = vz[j] + (h/2.0)*k_1_vz[j];
            
     k_2_x[j] = vx_temp[j];
     k_2_y[j] = vy_temp[j];
     k_2_z[j] = vz_temp[j];
        
     aceleracion(n_points,ID, ax_temp,ay_temp,az_temp,x_temp,y_temp,z_temp, masas);
            
     k_2_vx[j] = ax_temp[j];
     k_2_vy[j] = ay_temp[j];
     k_2_vz[j] = az_temp[j];
     
     x_temp[j] = x[j] + (h/2.0)*k_2_x[j];
     y_temp[j] = y[j] + (h/2.0)*k_2_y[j];
     z_temp[j] = z[j] + (h/2.0)*k_2_z[j];
     
     vx_temp[j] = vx[j] + (h/2.0)*k_2_vx[j];
     vy_temp[j] = vy[j] + (h/2.0)*k_2_vy[j];
     vz_temp[j] = vz[j] + (h/2.0)*k_2_vz[j];
            
     k_3_x[j] = vx_temp[j];
     k_3_y[j] = vy_temp[j];
     k_3_z[j] = vz_temp[j];
        
     aceleracion(n_points,ID,  ax_temp,ay_temp,az_temp,x_temp,y_temp,z_temp, masas);
            
     k_3_vx[j] = ax_temp[j];
     k_3_vy[j] = ay_temp[j];
     k_3_vz[j] = az_temp[j];
     
     x_temp[j] = x[j] + h*k_3_x[j];
     y_temp[j] = y[j] + h*k_3_y[j];
     z_temp[j] = z[j] + h*k_3_z[j];
     
     vx_temp[j] = vx[j] + h*k_3_vx[j];
     vy_temp[j] = vy[j] + h*k_3_vy[j];
     vz_temp[j] = vz[j] + h*k_3_vz[j];
            
     k_4_x[j] = vx_temp[j];
     k_4_y[j] = vy_temp[j];
     k_4_z[j] = vz_temp[j];
        
     aceleracion(n_points,ID , ax_temp,ay_temp,az_temp,x_temp,y_temp,z_temp, masas);
            
     k_4_vx[j] = ax_temp[j];
     k_4_vy[j] = ay_temp[j];
     k_4_vz[j] = az_temp[j];
            
     x[j] = x[j] + h*(1.0/6.0)*(k_1_x[j]+2*k_2_x[j]+2*k_3_x[j]+k_4_x[j]);
     y[j] = y[j] + h*(1.0/6.0)*(k_1_y[j]+2*k_2_y[j]+2*k_3_y[j]+k_4_y[j]);
     z[j] = z[j] + h*(1.0/6.0)*(k_1_z[j]+2*k_2_z[j]+2*k_3_z[j]+k_4_z[j]);
            
     vx[j] = vx[j] + h*(1.0/6.0)*(k_1_vx[j]+2*k_2_vx[j]+2*k_3_vx[j]+k_4_vx[j]);
     vy[j] = vy[j] + h*(1.0/6.0)*(k_1_vy[j]+2*k_2_vy[j]+2*k_3_vy[j]+k_4_vy[j]);
     vz[j] = vz[j] + h*(1.0/6.0)*(k_1_vz[j]+2*k_2_vz[j]+2*k_3_vz[j]+k_4_vz[j]);
    }
        
   Ecinetica = cinetica( n_points,vx, vy, vz, masas);
   Epotencial = potencial( n_points, x, y, z, masas);
        
   fprintf(energias,"%f %f \n", Ecinetica, Epotencial);
   fclose(energias);

  return 0;
}


FLOAT * get_memory(int n_points){
    FLOAT * x; 
    if(!(x = malloc(sizeof(FLOAT) * n_points))){
        printf("problem with memory allocation");
        exit(1);
    }
    return x;
}


FLOAT cinetica( int n_points, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *masas){
    int i;
    FLOAT energia;
    energia = 0.0;
    
    
    for (i=0; i<n_points; i++) {
      energia = energia + 0.5*masas[i]*(sqrt(pow(vx[i],2)+pow(vy[i],2)+pow(vz[i],2)));
    }
    return energia;
}


FLOAT potencial( int n_points, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *masas){
    int i,j;
    FLOAT r_ij;
    FLOAT energia;
    energia = 0.0;
    for (i=0; i<n_points; i++) {
        for (j=0; j<n_points; j++) {
            if (i!=j && j < 0) {
                r_ij = sqrt(pow((x[i] - x[j]),2.0) + pow((y[i] - y[j]),2.0) + pow((z[i] - z[j]),2.0));
                energia = energia -GRA * (masas[i] * masas[j]) / r_ij;
            }
        }
    }
    return energia;
}


void aceleracion( int n_points, int *ID, FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *masas){
    int i,j;
    FLOAT r_ij;
    for(i=0;i<n_points;i++){
        ax[i]=0.0;
        ay[i]=0.0;
        az[i]=0.0;
        
        for(j=0;j<n_points;j++){
            if(j!=i && ID[j]<0){
                r_ij = (pow((x[i] - x[j]),2.0) + pow((y[i] - y[j]),2.0) + pow((z[i] - z[j]),2.0));
                r_ij = sqrt(r_ij);
                ax[i] = ax[i] - GRA *masas[j]/ pow(r_ij,1.5) * (x[i] - x[j]);
                ay[i] = ay[i] - GRA *masas[j]/ pow(r_ij,1.5) * (y[i] - y[j]);
                az[i] = az[i] - GRA *masas[j] / pow(r_ij,1.5) * (z[i] - z[j]);
	    }
	}
    }
}
