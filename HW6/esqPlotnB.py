"""
Programa para graficar las condiciones de las particulas a partir de archivos de texto
Autores: David Andres Davila
         Juan Sebastian Calderon
"""
import numpy as np

def CargarDatos(archivo):
    """
    Carga los datos de un archivo.
    Entradas: -archivo:el archivo de texto cuyos datos se cargaran
    salidas: -datos: matriz con los datos del archivo
    """
    datos = archivo.readlines()
    return datos

def Graficar(condiciones):
    """
    genera una grafica a partir de las condiciones de los 10^24 cuerpos

    Entradas: -condiciones: matriz con las condiciones de la galaxia
    salidas: 
    """
    np.plot(condiciones)
    return
