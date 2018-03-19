# -*- coding: utf-8 -*-

"""
A partir de este código se obtiene la tabla .fits con toda la información 
relativa al proyecto GAMA que se va a utilizar en el matching; el fichero de
salida será uno de los ficheros de entrada de main_matching_gh.py

  Python 3.6.0    scipy 0.18    numpy 1.11.3    astropy 1.3

@author: Javier Gutiérrez Solórzano
"""
#==============================================================================
#    Nombres de los ficheros de entrada/salida y su localización.
#==============================================================================
#--------------------------
#   Directorio en el que se encuentran los modulos.
#--------------------------
import sys
sys.path.append("../pakages/")
#--------------------------
#   Directorios de los ficheros.
#--------------------------
inroute = './inroute/'
outroute = './outroute_gama/'
#--------------------------
#   Nombre de los ficheros.
#--------------------------
infile_catalogo = "GamaCoreDR1_v1.fits"
outfile = "Gama.fits"

import get
import rojo

#==============================================================================
#  Desplazamiento al rojo disponible de los objetos de GAMA.
#============================================================================== 
z_HELIO = get.columna_fits(inroute + infile_catalogo,6)
z_QUALITY = get.columna_fits(inroute + infile_catalogo,7)

[z_gama, E_z_gama, flag_gama] = rojo.z_gama(z_HELIO,z_QUALITY)

#==============================================================================
# Creación del fichero fits de 'salida'.
#==============================================================================
     
outfile_array_GAMA_IAU_ID = get.columna_fits(inroute + infile_catalogo,0)
outfile_array_RA_g = get.columna_fits(inroute + infile_catalogo,3)
outfile_array_DEC_g = get.columna_fits(inroute + infile_catalogo,4)

from astropy.io import fits

#-------------------------- 
#   Creamos las columnas del fichero de salida.
#--------------------------  

col0 = fits.Column(name='GAMA_IAU_ID', format='23A'
                   , array=outfile_array_GAMA_IAU_ID)
col1 = fits.Column(name='RA_G', format='D',  array=outfile_array_RA_g
                   , unit='Degrees')
col2 = fits.Column(name='DEC_G', format='D', array=outfile_array_DEC_g
                   , unit='Degrees')

    # Información proporcinada por HATLAS.
col3 = fits.Column(name='z_HELIO', format='D', array=z_HELIO)
col4 = fits.Column(name='z_QUALITY', format='D', array=z_QUALITY)

    # Desplazamientos al rojo considerados válidos.
col5 = fits.Column(name='Z_G', format='D', array=z_gama)
col6 = fits.Column(name='Error_Z_G', format='D', array=E_z_gama)
col7 = fits.Column(name='Flag_gama', format='D', array=flag_gama)

cols = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6, col7])

tbhdu = fits.BinTableHDU.from_columns(cols)
   
tbhdu.writeto(outroute + outfile,overwrite=True)