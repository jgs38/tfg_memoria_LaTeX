# -*- coding: utf-8 -*-
"""
Este es el ejecutable con el que se realizan las simulaciones. Para ejecutarse
es necesario disponer de los ficheros de salida de los ejecutables
 main_gama.py y main_hatlas.py . Para realizar la simulación a cada observación
se le ha reasignado una posición aleatoria sobre una región con un área 
equivalente al área del catálogo al que pertenece. Después estas dos regiones 
se han solapado en un área igual al área de intersección de los catálogos 
originales y se ha procedido ha realizar un conteo de las contrapartidas
 con la función "get.selecion_candidatos".

  Python 3.6.0    numpy 1.11.3    astropy 1.3

@author: Javier Gutiérrez Solórzano
"""

#==============================================================================
#    Nombres de los ficheros de entrada/salida y su localización.
#==============================================================================
#--------------------------
#   Directorio en el que se encuentran los módulos.
#--------------------------
import sys
sys.path.append("../packages/")
#--------------------------
#   Directorios de los ficheros.
#--------------------------
inroute_g = './outroute_gama/'
inroute_h = './outroute_hatlas/'
outroute = './outroute_matching/'
#--------------------------
#   Nombre de los ficheros
#--------------------------
infile_catalogo1 = "Gama.fits"
infile_catalogo2 = "HAtlas.fits"
infile_objetos_interesantes="nombres_candidatos.txt"
outfile = "matching_gh"

#   Soporte de unidades físicas.
from astropy import units as u 

# Para el cálculo del factor de bayes espectroscópico.
z_max=3.501

import numpy as np

fwhm_sigma=2*np.sqrt(2*np.log(2))

#   Separación angular máxima entre los objetos que forman el "matching".
sep=50*u.arcsec

import get
import rojo

# Límites de las áreas en las que se realiza la simulación
#   GAMA
DEC_a_gama = 0*u.deg
DEC_b_gama = 12*u.deg

RA_a_gama = 0*u.deg

area_gama = 144*u.deg*u.deg
area_gama_sr = (area_gama.to(u.sr)).value
RA_b_gama_rad = (area_gama_sr/np.sin(DEC_b_gama))*u.rad
RA_b_gama = RA_b_gama_rad.to(u.deg)

#   HATLAS
DEC_hatlas = 12*u.deg

area_hatlas = 161*u.deg*u.deg
area_hatlas_sr = (area_hatlas.to(u.sr)).value
RA_hatlas_rad = (area_hatlas_sr/np.sin(DEC_hatlas))*u.rad
RA_hatlas = RA_hatlas_rad.to(u.deg)

#   Intersección
DEC_interseccion = 12*u.deg

area_interseccion = 130*u.deg*u.deg
area_interseccion_sr = (area_interseccion.to(u.sr)).value
RA_interseccion_rad = (area_interseccion_sr/np.sin(DEC_interseccion))*u.rad
RA_interseccion = RA_interseccion_rad.to(u.deg)


DEC_a_hatlas = 0*u.deg
DEC_b_hatlas = 12*u.deg

RA_a_hatlas = RA_b_gama-RA_interseccion
RA_b_hatlas = RA_a_hatlas+RA_hatlas

#-------------------------- 
#   Preparación del catálogo GAMA.
#-------------------------- 
#   Todos los valores del redshift de los objetos de GAMA.
z_g=get.columna_fits(inroute_g + infile_catalogo1,5)

#   Nos devuelve el identificador (posición en el fichero) de aquellos objetos
# de GAMA de los que disponemos de una medida del redshift válida.
idg_catalogo1 = rojo.id_z_validos(z_g)

#-------------------------- 
#   Preparación del catálogo HATLAS.
#--------------------------
#   Todos los valores del redshift de los objetos de HATLAS.
z_h=get.columna_fits(inroute_h + infile_catalogo2,14)

#   Nos devuelve el identificador de aquellos objetos de HATLAS de los que
# disponemos de una medida del redshift válida.
idh_catalogo2 = rojo.id_z_validos(z_h)   

id_halo_h_espec_media=[]
id_halo_h_fot_media=[]
id_halo_h_media=[]
id_halo_negrello_gonzalez_media=[]

from astropy.coordinates import SkyCoord
import factor_bayes
    
for i in range(0,1000,1):   
                             
    [ra_catalogo_GAMA,dec_catalogo_GAMA]= get.uniform_spherical_distribution(
        RA_a_gama.value,RA_b_gama.value,DEC_a_gama.value,DEC_b_gama.value,
        len(idg_catalogo1))
                                 
    [ra_catalogo_HATLAS,dec_catalogo_HATLAS] = \
        get.uniform_spherical_distribution(RA_a_hatlas.value,RA_b_hatlas.value,
        DEC_a_hatlas.value,DEC_b_hatlas.value,len(idh_catalogo2))
    
    #==========================================================================
    #   Matching para identificar aquellos objetos del catálogo H-Atlas que se 
    # encuentran a una distancia angular inferior  a "sep" segundos de arco de 
    # un objeto de GAMA.
    #==========================================================================     
    h = SkyCoord(ra=ra_catalogo_HATLAS*u.degree,
                 dec=dec_catalogo_HATLAS*u.degree)  
    g = SkyCoord(ra=ra_catalogo_GAMA*u.degree, dec=dec_catalogo_GAMA*u.degree)
    
    # d2d   Separación angular entre los objetos que forman el "matching".
    # _ indica que no queremos realizar ese cálculo.
    [idg_arrays, idh_arrays, d2d, _] = h.search_around_sky(g, sep)
    
    #  Posiciones de los objetos en el catalogos de entrada (idxg,idxh):
    idxg=get.ordena_Array(idg_arrays, idg_catalogo1)
    idxh=get.ordena_Array(idh_arrays, idh_catalogo2)
    
    #   El primer elemento del array en Python es un 0 mientras que el la tabla 
    # fits comienza por un 1. Las dos líneas siguientes resuelven ésto.
    idg=get.posicion_fits(idxg)
    idh=get.posicion_fits(idxh)
    
    #==========================================================================
    #   Factores de Bayes: posicional, fotométrico y conjunto
    #========================================================================== 
    #   GAMA
    # Errores posicionales del catálogo GAMA.
    sg=(0.7/fwhm_sigma)*u.arcsec
    sg=sg.to(u.radian)
    sg=sg.value
    
    zg = get.ordena_Array(idxg,z_g)
    szg = get.ordena_Array(idxg,
                           get.columna_fits(inroute_g + infile_catalogo1,6))                                                 
    flag_g = get.ordena_Array(idxg,
                              get.columna_fits(inroute_g + infile_catalogo1,7))
    #   HALTAS
    # Errores posicionales del catálogo HALTAS.
    sh=(17.98/fwhm_sigma)*u.arcsec  
    sh=sh.to(u.radian)
    sh=sh.value
    
    zh = get.ordena_Array(idxh,z_h)
    szh = get.ordena_Array(idxh,
                             get.columna_fits(inroute_h + infile_catalogo2,15))
    flag_h = get.ordena_Array(idxh,
                             get.columna_fits(inroute_h + infile_catalogo2,16))
    #   Distancias angulares en radianes (necesario para calcular el factor de 
    # bayes posicional).
    d2d = d2d.to(u.radian)
    d2d = d2d.value
    
    #   Cálculo de los factores de Bayes
    [bayes_posicional, bayes_z, bayes_conjunto] = factor_bayes.bayes_conjunto(
            sh,sg,d2d,zh,szh,zg,szg,z_max) 
    
    outfile_array_S250=get.ordena_Array(idxh
                            ,get.columna_fits(inroute_h + infile_catalogo2,3))
    outfile_array_S350=get.ordena_Array(idxh
                            ,get.columna_fits(inroute_h + infile_catalogo2,4))
    outfile_array_S500=get.ordena_Array(idxh
                            ,get.columna_fits(inroute_h + infile_catalogo2,5))
    
    
    [id_halo_h_espec,id_halo_h_fot,id_halo_negrello_gonzalez]= \
        get.selecion_candidatos(outfile_array_S250,outfile_array_S350,
        outfile_array_S500,flag_h,zh,zg,bayes_posicional,bayes_conjunto)
        
    print(len(id_halo_h_espec)+len(id_halo_h_fot),len(id_halo_h_espec),
          len(id_halo_h_fot),len(id_halo_negrello_gonzalez))

    id_halo_h_espec_media.append(len(id_halo_h_espec))
    id_halo_h_fot_media.append(len(id_halo_h_fot))
    id_halo_h_media.append(len(id_halo_h_espec)+len(id_halo_h_fot))
    id_halo_negrello_gonzalez_media.append(len(id_halo_negrello_gonzalez))

print("-Número total de candidatos:",np.mean(id_halo_h_media),chr(177),
      np.std(id_halo_h_media),"\n-HATLAS fotométrico:",
      np.mean(id_halo_h_espec_media), chr(177),np.std(id_halo_h_espec_media),
      "\n-HATLAS espectroscópico:",
      np.mean(id_halo_h_fot_media),chr(177),np.std(id_halo_h_fot_media),
      "\n-Candidatos que cumplen los criterios de Gonzalez y/o Negrello:",
      np.mean(id_halo_negrello_gonzalez_media),chr(177),
             np.std(id_halo_negrello_gonzalez_media))