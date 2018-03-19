# -*- coding: utf-8 -*-

"""
A partir de éste código se obtiene una tabla .fits que será uno de los ficheros
de entrada al ejecutable main_matching_gh.py ; El fichero de salida contendrá
todos los datos relativos al catálogo HATLAS que se van a utilizar en el 
matching. El cálculo de los desplazamientos al rojo a partir del ajuste a la 
SED empírica es una tarea que requiere mucho tiempo de cálculo (mas de 1 hora
para los 120000 objetos de HATLAS_DR1 con un procesador Intel Core2Duo). Como 
puede ser desesperante, a medida que se ejecuta el código aparecen distintos
mensajes en la consola sobre el tiempo de ejecución que permiten conocer las 
instrucciones que se están ejecutando.

  Python 3.6.0    scipy 0.18    numpy 1.11.3    astropy 1.3     datetime
  
@author: Javier Gutiérrez Solórzano
"""
# Para conocer el tiempo de ejecución.
import datetime
tiempo_inicio = datetime.datetime.now()

#==============================================================================
#    Nombres de los ficheros de entrada/salida y su localización.
#==============================================================================
#--------------------------
#   Directorio en el que se encuentran los modulos.
#--------------------------
import sys
sys.path.append("../packages/")
#--------------------------
#   Directorios de los ficheros.
#--------------------------
inroute = './inroute/'
outroute = './outroute_hatlas/'
#--------------------------
#   # Nombre de los ficheros.
#--------------------------
infile_template = 'SMM_template_norm.sed'
infile_catalogo1 = "HATLAS_DR1_CATALOGUE_V1.2.FITS"
outfile = "HAtlas.fits"

import datetime
tiempo_inicio = datetime.datetime.now()

#==============================================================================
#   Desplazamiento al rojo de los objetos HATLAS a partir del ajuste.
#============================================================================== 
#   Soporte de unidades físicas.
from astropy import units as u 
#   Soporte de magnitudes físicas. Solo se utiliza la velocidad de la luz.
from astropy.constants import c 

import get

#-------------------------- 
#   Entrada:
#   Leemos el fichero de la SED modelo. El fichero contiene dos columnas; la 
# primera representa las longitudes de onda en Angstrom y la segunda valores
# valores de la densidad espectral de flujo (F_(lambda)), en "unidades" (está
# normalizado en 5500A ) de : (erg) * (cm)**(-2) * (s)**(-1) * (A)**(-1)
#-------------------------- 
[x_sed, y_sed] = get.sed_func(inroute + infile_template) 
#   Unidades a las longitudes de onda.
x_sed=x_sed*u.AA
#   Unidades de densidad espectral de flujo (F_(lambda))
y_sed=y_sed*((u.erg)* (u.cm)**(-2)* (u.s)**(-1)* (u.AA)**(-1))
#   Cambio de F_(lambda) a F_(ípsilon) multiplicando por (lambda**2)/c
y_sed=y_sed*(x_sed**2)*(c**-1)
#   Dejo los arrays con las unidades que tienen los valores experimentales que 
# voy a utilizar.
x_sed=x_sed.to(u.micron)
y_sed=y_sed.to(u.Jy)
 
#-------------------------- 
#   Obtenemos los valores experimentales (mediciones realizadas por HATLAS) de
# todos los objetos de un catálogo a partir de la lectura de un fichero con 
# formato .fits
#-------------------------- 
x_experimental = [250,350,500]*u.micron 
[y_experimental, y_error]= get.datos_experimentales(inroute + infile_catalogo1)
#    Le asignamos unidades fisicas.
# EN ESTE CASO LAS MEDIDAS DE DENSIDAD DE FLUJO SE ENCUENTRAN EN mJy!
y_experimental = y_experimental*u.mJy
y_error = y_error*u.mJy
#    Hacemos un cambio a las unidades. 
y_experimental=y_experimental.to(u.Jy)
y_error = y_error.to(u.Jy)

#-------------------------- 
#   Salida.
#-------------------------- 
b = datetime.datetime.now()
print('Tiempo de ejecucion:',b-tiempo_inicio,'\n',
      'Comienza el calculo de los desplazamientos al rojo ')

import rojo
[z_ajuste, E_z_ajuste] = rojo.z_phot_hatlas(x_sed,y_sed,x_experimental 
                                            ,y_experimental,y_error)
b = datetime.datetime.now()
print('El calculo ha durado ',b-tiempo_inicio)

#==============================================================================
#  Desplazamiento al rojo HATLAS.
#============================================================================== 
gsq_flag=get.columna_fits(inroute + infile_catalogo1,32)
z_spec=get.columna_fits(inroute + infile_catalogo1,44)
z_qual=get.columna_fits(inroute + infile_catalogo1,45)

# Rango de z en el que "z_phot_hatlas" se considera valido (z= 1 - 3.5).
z_min=1
z_max=3.5

# Crea las columnas con los redshifts considerados válidos
[z_hatlas, E_z_hatlas, flag_hatlas]=rojo.z_hatlas(gsq_flag,z_spec,z_qual
                                                         ,z_ajuste,z_min,z_max)

#==============================================================================
# Creación del fichero .fits
#==============================================================================  
outfile_array_HATLAS_DR1_CATALOGUE = get.columna_fits(inroute 
                                                          + infile_catalogo1,0)
outfile_array_RA_h = get.columna_fits(inroute + infile_catalogo1,2)
outfile_array_DEC_h = get.columna_fits(inroute + infile_catalogo1,3)

from astropy.io import fits
#-------------------------- 
#   Creamos las columnas del fichero de salida.
#--------------------------  
col0 = fits.Column(name='HATLAS_DR1_CATALOGUE', format='23A'
                   , array=outfile_array_HATLAS_DR1_CATALOGUE)
col1 = fits.Column(name='RA_H', format='D',  array=outfile_array_RA_h
                   , unit='Degrees')
col2 = fits.Column(name='DEC_H', format='D', array=outfile_array_DEC_h
                   , unit='Degrees')

col3 = fits.Column(name='S250', format='D', array=get.columna_fits(inroute 
                                             + infile_catalogo1,4), unit='mJy')
col4 = fits.Column(name='S350', format='D', array=get.columna_fits(inroute 
                                             + infile_catalogo1,5), unit='mJy')
col5 = fits.Column(name='S500', format='D', array=get.columna_fits(inroute 
                                             + infile_catalogo1,6), unit='mJy')

col6 = fits.Column(name='E250', format='D', array=get.columna_fits(inroute 
                                             + infile_catalogo1,7), unit='mJy')
col7 = fits.Column(name='E350', format='D', array=get.columna_fits(inroute 
                                             + infile_catalogo1,8), unit='mJy')
col8 = fits.Column(name='E500', format='D', array=get.columna_fits(inroute 
                                             + infile_catalogo1,9), unit='mJy')

#    Desplazamiento al rojo calculado.
col9 = fits.Column(name='Z_AJUSTE', format='D', array=z_ajuste)
col10 = fits.Column(name='Error_Z_AJUSTE', format='D', array=E_z_ajuste)

#    Información proporcinada por HATLAS.
col11 = fits.Column(name='z_SPEC', format='D', array=z_spec)
col12 = fits.Column(name='z_QUAL', format='D', array=z_qual)
col13 = fits.Column(name='GSQ_FLAG', format='D', array=gsq_flag)

#    Desplazamientos al rojo considerados válidos.
col14 = fits.Column(name='Z_H', format='D', array=z_hatlas)
col15 = fits.Column(name='Error_Z_H', format='D', array=E_z_hatlas)
col16 = fits.Column(name='Flag_hatlas', format='D', array=flag_hatlas)

#    Ahora fabrico un objeto con esas columnas:
cols = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6, col7, col8
                      , col9, col10, col11, col12, col13, col14, col15, col16])

#    Ese objeto se mete en una tabla HDU.
tbhdu = fits.BinTableHDU.from_columns(cols)

tbhdu.writeto(outroute + outfile,overwrite=True)

#==============================================================================
b = datetime.datetime.now()
print('El tiempo de ejecución total ha sido de ',b-tiempo_inicio)