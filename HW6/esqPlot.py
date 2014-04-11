"""
Programa para graficar las condiciones de las galaxias a partir de archivos de taxto
Autores: David Andres Davila
         Juan Sebastian Calderon
"""
import numpy as np

def CargarDatos(archivo):
    """
    Carga los datos de un archivo, por lo que deberá usarse un número arbitrario de veces (una por cada archivo)

    Entradas: -archivo:el archivo de texto cuyos datos se cargaran
    salidas: -datos: matriz con los datos del archivo
    """
    datos = archivo.readlines()
    return datos

def Graficar(condiciones):
    """
    genera una grafica a partir de las condiciones de una galaxia en el tiempo T

    Entradas: -condiciones: matriz con las condiciones de la galaxia
    salidas: 
    """
    np.plot(condiciones)
    return

