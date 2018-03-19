# -*- coding: utf-8 -*-

"""
Este módulo contiene las funciones con las que se ha obtenido el redshift 
de las observaciones de los catálogos H-ATLAS y GAMA. 

  Python 3.6.0    scipy 0.18    numpy 1.11.3

@author: Javier Gutiérrez Solórzano
"""

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import numpy as np


def z_phot_objeto(x_sed, y_sed, x_experimental, y_experimental_i, y_error_i, 
           abs_s=True):
    """
    Ajusta (utilizando dos parámetros) la SED teórica de la galaxia modelo
    a las medias experimentales del flujo espectroscópico.
    
    Parámetros:
    x_sed ------------ Array con los valores 'x' de la SED teórica.
    y_sed ------------ Array con los valores 'y' de la SED teórica.
    x_experimental --- Array de tres elementos con los valores de la longitud 
                       de onda a la que se realizan las medidas del flujo.
                       Nosotros siempre utilizaremos el mismo array que esta 
                       definido en el modulo 'main.py':
                           [250,350,500]*u.micron 
    y_experimental_i-- Array que contiene las medidas del flujo un objeto. En
                       En este trabajo se dispone de tres medidas para cada 
                       objeto, asi que cada array contendrá 3 elementos
    y_error_i--------- Array que contiene error del flujo asociado al array 
                       'y_experimental'.
    Return:
    popt-------------------- Valores óptimos del ajuste.
    np.sqrt(np.diag(pcov))-- Las diagonales proporcionan la varianza de la 
                             de los parámetros de ajuste.
    """
    #-------------------------- 
    #   Definición de la funcion desplazada.
    #--------------------------  
    def func_desplazada(x,K,C):
        """
        Obtiene una función mediante interpolación lineal del los puntos y 
        multiplica los valores asociados a cada uno de los ejes de coordenadas
        'x' e 'y' por un factor ('K' y 'C' respectivamente). 
        
        Parámetros:
        K -- Factor multiplicativo para el eje de abscisas.
        C -- Factor multiplicativo para el eje de ordenadas.

        Return:
        y -- Valores de la función obtenida por interpolcación lineal 
             asociados al eje de ordenadas tras el ajuste.
        
        """
        x_desplazada = x_sed*(K)
        f1 = interp1d(x_desplazada,y_sed,bounds_error=False,
                      fill_value=np.min(y_sed))
        y = f1(x)*C
        return y
        #   p0=[1.8, 1e-3] Limita el los valores en los que se realiza
        # el ajuste, aunque en principio no es necesario.
        
    popt, pcov = curve_fit(lambda x, K, C: func_desplazada(x, K, C), 
                           x_experimental, y_experimental_i, p0=[1.8, 1e-4], 
                           sigma=y_error_i)         
    
    return popt, np.sqrt(np.diag(pcov))

def z_phot_hatlas(x_sed,y_sed,x_experimental,y_experimental,y_error):  
    """
    Hace uso de la función 'objeto' para obtener el redshift de todo el 
    catálogo HATLAS. 
        
    Parámetros:
    x_sed ----------- Array con los valores 'x' de la SED teórica.
    y_sed ----------- Array con los valores 'y' de la SED teórica.
    x_experimental -- Array de tres elementos con los valores de la longitud 
                      de onda a la que se realizan las medidas del flujo.
                      Nosotros siempre utilizaremos el mismo array que está
                      definido en el modulo 'main.py':
                          [250,350,500]*u.micron 
    y_experimental -- Array que contiene las medidas del flujo de un conjunto
                      de objetos. Cada elemento del array debe de ser un array
                      con las medidas de flujo de un objeto.
    y_error --------- Array que contiene error del flujo asociado al array 
                      'y_experimental'.
                      
    Return:
    z -- Array con los desplazamientos al rojo del catálogo HATLAS completo.
    S -- Desviación estándar como medida de la calidad del ajuste.
        
    """
    z=[]
    S=[]
    # Cambiar len(y_experimental) por otro valor inferior para las pruebas si 
    # el catálogo es muy extenso.
    for i in range(0,len(y_experimental),1):
        [z_objeto,S_objeto] = z_phot_objeto(x_sed, y_sed, x_experimental,
        y_experimental[i], y_error[i])
        
        z.append(z_objeto[0]-1)
        S.append(S_objeto[0])

    return z, S


def id_z_validos(z_catalogo):
    """
    Proporciona el identificador de aquellas observaciones del catálogo
    HATLAS con z "válido".
    
    Parámetros:
    z_catalogo ------- Redshifts proporcionados por HATLAS.
    
    Return:
    contador_h ------- Posición que ocupa la observación en el catálogo HATLAS.
        
    """
    contador_h=[]
    for i in range(0,len(z_catalogo),1):
        if np.any(z_catalogo[i]!=-99): 
            contador_h.append(i)
    return contador_h

def z_hatlas(gsq_flag,z_spec,z_qual,z_ajuste,z_min,z_max):
    """
    Esta función se utiliza para asignar los valores definitivos a los
    redshifts a las observaciones realizadas por HATLAS. Además añade una
    etiqueta que permite identificar el origen de esa medida.
    
    Parámetros:

    gsq_flag ---- Clasificación del tipo de objeto proporcionada por HATLAS.
    z_spec ------ Redshift proporcionado por HATLAS.
    z_qual ------ Etiqueta de calidad asignada por HATLAS.
    z_ajuste ---- Estimación del redshift proporiconado por el ajuste de la
                  función "z_phot_hatlas".
    z_min ------- Valor mínimo aceptable para el z del ajuste.
    z_max ------- Valor máximo aceptable para el z del ajuste.
    
    Return:
        
    z_hatlas -------- Array con los valores definitivos de z asignados a las
                      observaciones de HATLAS.
    sigma_z_hatlas -- Errores asignados a las medias de "z_hatlas".
    flag ------------ Etiqueta de "calidad" sobre la medida del redshift.
                      *   11, 12, 13: medidas "espectroscópicas"
                        procedentes del catálogo HATLAS.
                      *   2: redshift fotométrico generado con ANNZ.
                      *   3: z obtenido mediante la función "z_phot_hatlas".
                      *   -99 para el resto de casos.    
    """
    z_hatlas=[]
    sigma_z_hatlas=[]
    flag=[]
    for i in range(0,len(z_ajuste),1):
        #0=gal, 1=star, 2=quasar, 3=quasar candidate based on colour selection.
        if np.any(gsq_flag[i]==0) or np.any(gsq_flag[i]==2) \
                 or np.any(gsq_flag[i]==3):
            # >=3 is OK to use (redshift espectroscópico)
            if np.any(z_qual[i]>=3):
                z_hatlas.append(z_spec[i])
                # Puede haber redshifts negativos, asi que abs()!
                sigma_z_hatlas.append(0.00011064*(1+np.abs(z_spec[i])))
                if np.any(z_qual[i]==5):
                    flag.append(11)
                elif np.any(z_qual[i]==4):
                    flag.append(12)
                else:
                    flag.append(12)
            # <3 (redshift fotométrico generado con ANNZ)
            elif (np.any(z_qual[i]==2) or np.any(z_qual[i]==1)) \
                 and (np.any(z_spec[i]>=0) and np.any(z_spec[i]<=0.7)):
                z_hatlas.append(z_spec[i])
                sigma_z_hatlas.append(0.023)
                flag.append(2)
            # Es una galaxia de la que no se ha podido obtener el redshift por 
            # ninguno de los métodos anteriores.
            elif  np.any(gsq_flag[i]==0) and np.any(z_ajuste[i] >= z_min) \
                        and np.any(z_ajuste[i] <= z_max) : 
                z_hatlas.append(z_ajuste[i])
                sigma_z_hatlas.append(0.115*(1+z_ajuste[i]))
                flag.append(3)
            else:
                z_hatlas.append(-99)
                sigma_z_hatlas.append(-99)
                flag.append(-99)
        else:
            z_hatlas.append(-99)
            sigma_z_hatlas.append(-99)
            flag.append(-99)
                
    return z_hatlas, sigma_z_hatlas, flag


def z_gama(z_HELIO,z_quality):
    """
    Esta función se utiliza para obtener los valores del redshift 
    espectroscópico que nos proporciona el proyecto GAMA. Los valores de z 
    asociados a los objetos del catálogo GAMA se encuentran en una columna
    del fichero 'GamaCoreDR1_v1.fits', pero resulta que no todos los valores 
    que aparecen aquí son válidos. El proyecto indica los z no válidos del 
    siguiente modo:
        -Si el corrimiento al rojo de objeto es desconocido, '9999'
        -Si se divulgará en una publicación posterior, ' − 2'
        -Si z se ha observado pero es una medida pobre, ' − 0.9'
    Esta función "filtra" aquellos valores no válidos y proporciona el array 
    con los valores que si podemos utilizar.
    
    Por otra parte, resulta necesario asociar un identificador a cada
    observación del catálogo GAMA; un entero que indica la fila que ocupaba 
    la observación en el fichero 'GamaCoreDR1_v1.fits'.
    
    Parámetros:
    z_HELIO  --------- Valores del redshift proporcionados por GAMA.
    z_quality -------- Factor de calidad asignado a z por GAMA.
    
    Return:
    z_gama ----------- Redshifts que se consideramos válidos.
    sigma_z_gama ----- Error que hemos asignado a cada valor de "z_gama".
    flag ------------- Etiqueta de calidad sobre la medida del redshift,
                       *   11, 12, 13 para medidas espectroscópicas
                       *   -99 resto de casos.
        
    """
    z_gama=[]
    sigma_z_gama=[]
    flag=[]
    for i in range(0,len(z_HELIO),1):
        if np.any(z_HELIO[i]!=-2) and np.any(z_HELIO[i]!=9999):
            if np.any(z_quality[i]>=3):
                z_gama.append(z_HELIO[i])
                # Puede haber redshifts negativos, asi que abs()!
                sigma_z_gama.append(0.00011064*(1+np.abs(z_HELIO[i])))
                if np.any(z_quality[i]==5):
                    flag.append(11)    
                elif np.any(z_quality[i]==4):
                    flag.append(12)
                else:
                    flag.append(13)
            else:
                z_gama.append(-99)
                sigma_z_gama.append(-99)
                flag.append(-99)
        else:
          z_gama.append(-99)
          sigma_z_gama.append(-99)
          flag.append(-99)
            
    return z_gama, sigma_z_gama, flag