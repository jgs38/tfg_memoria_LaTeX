# -*- coding: utf-8 -*-
"""
Este ejecutable realiza el emparejamiento entre las observaciones de los 
catálogos GAMA y HATLAS e identifica los candidatos que forman un sistema 
lente gravitatoria según nuestro criterio. También proporciona varias gráficas
que ayudan al análisis de los resultados.

  Python 3.6.0    scipy 0.18   matplotlib 2.0.0     numpy 1.11.3    astropy 1.3

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

#   Para el cálculo del factor de bayes espectroscópico.
z_max=3.501

import numpy as np
#   Resolucion angular de los catálogos a partir de la FWHM.
fwhm_sigma=2*np.sqrt(2*np.log(2))
sg=(0.7/fwhm_sigma)*u.arcsec     #  GAMA 
sh=(17.98/fwhm_sigma)*u.arcsec   #  H-ATLAS 

#   Separación angular máxima entre los objetos que forman el "matching".
sep=54*u.arcsec

#==============================================================================
#   Lectura de los datos necesarios para el matching.
#============================================================================== 
import get
import rojo

#-------------------------- 
#   Preparación del catálogo GAMA.
#-------------------------- 
#   Todos los valores del redshift de los objetos de GAMA.
z_g=get.columna_fits(inroute_g + infile_catalogo1,5)

#   Nos devuelve el identificador (posición en el fichero) de aquellos objetos
# de GAMA de los que disponemos de una medida del redhsift válida.
idg_catalogo1 = rojo.id_z_validos(z_g)   

# Posiciones de todos los objetos del catálogo GAMA
ra_c_GAMA=get.columna_fits(inroute_g + infile_catalogo1,1)
dec_c_GAMA=get.columna_fits(inroute_g + infile_catalogo1,2)

#   Posiciones de aquellos objetos de GAMA con z válidos.
ra_catalogo_GAMA = get.ordena_Array(idg_catalogo1,ra_c_GAMA) #  RA (0 y 360)
dec_catalogo_GAMA = get.ordena_Array(idg_catalogo1,dec_c_GAMA) # DEC (-90 y 90)

#-------------------------- 
#   Preparación del catálogo HATLAS.
#--------------------------
#   Todos los valores del redshift de los objetos de HATLAS.
z_h=get.columna_fits(inroute_h + infile_catalogo2,14)

#   Identificador de aquellos observaciones de HATLAS con z válido.
idh_catalogo2 = rojo.id_z_validos(z_h)   

# Posiciones de todos los objetos del catálogo HATLAS.
ra_c_HATLAS=get.columna_fits(inroute_h + infile_catalogo2,1)
dec_c_HATLAS=get.columna_fits(inroute_h + infile_catalogo2,2)

#   Posiciones de los objetos HATLAS.
ra_catalogo_HATLAS = get.ordena_Array(idh_catalogo2,ra_c_HATLAS)
dec_catalogo_HATLAS = get.ordena_Array(idh_catalogo2,dec_c_HATLAS)

#==============================================================================
#   Matching para identificar aquellas observaciones de H-Atlas que se encuen-
# tran a una distancia angular inferior a "sep" de otra observación en GAMA.
#============================================================================== 
from astropy.coordinates import SkyCoord

h = SkyCoord(ra=ra_catalogo_HATLAS*u.degree, dec=dec_catalogo_HATLAS*u.degree)  
g = SkyCoord(ra=ra_catalogo_GAMA*u.degree, dec=dec_catalogo_GAMA*u.degree)

# d2d   Separación angular entre los objetos que forman el "matching".
# _ indica que no queremos realizar ese cálculo.
[idg_arrays, idh_arrays, d2d, _] = h.search_around_sky(g, sep)

# Posiciones de las observacionesque forman el "matching" en los 
# catálogos de entrada (idxg,idxh):
idxg=get.ordena_Array(idg_arrays, idg_catalogo1)
idxh=get.ordena_Array(idh_arrays, idh_catalogo2)

#   El primer elemento del array en Python es un 0 mientras que el la tabla 
# fits comienza por un 1. Las dos líneas siguientes resuelven ésto.
idg=get.posicion_fits(idxg)
idh=get.posicion_fits(idxh)

#==============================================================================
#   Factores de Bayes: posicional, fotométrico y conjunto
#============================================================================== 

#   GAMA
sg=sg.to(u.radian)
sg=sg.value
zg = get.ordena_Array(idxg,z_g)
szg = get.ordena_Array(idxg,get.columna_fits(inroute_g + infile_catalogo1,6))
flag_g = get.ordena_Array(idxg,get.columna_fits(inroute_g + infile_catalogo1,
                                                7))
#   HALTAS
sh=sh.to(u.radian)
sh=sh.value
zh = get.ordena_Array(idxh,z_h)
szh = get.ordena_Array(idxh,get.columna_fits(inroute_h + infile_catalogo2,15))
flag_h = get.ordena_Array(idxh,get.columna_fits(inroute_h + infile_catalogo2,
                                                16))

#   Distancias angulares en radianes 
d2d = d2d.to(u.radian)
d2d = d2d.value

import factor_bayes
# sh     número   error posicional del catalogo1
# sg     número   error posicional en el catalogo2
# d2d    array    distancias entre los objetos
# zh     array    valores de z en el catalogo1
# szh    array    error de z de los objetos del catalogo1
# zg     array    valores de z sobre los objetos del catalogo2
# szg    array    error de z sobre los objetos del catalogo2
# z_max  número   valor maximo considerado

[bayes_posicional, bayes_z, bayes_conjunto] = factor_bayes.bayes_conjunto(sh,
                                                    sg,d2d,zh,szh,zg,szg,z_max)

#==============================================================================
#   Creación del fichero fits con los datos del matching completo
#==============================================================================   
              
#   outfile_array_IDNAME_GAMA = idxg
outfile_array_GAMA_IAU_ID = get.ordena_Array(idxg,get.columna_fits(inroute_g 
                                                         + infile_catalogo1,0))
outfile_array_RA_g = get.ordena_Array(idxg,get.columna_fits(inroute_g 
                                                         + infile_catalogo1,1))
outfile_array_DEC_g = get.ordena_Array(idxg,get.columna_fits(inroute_g 
                                                         + infile_catalogo1,2))

#   outfile_array_IDNAME_HATLAS = idxh
outfile_array_HATLAS_DR1_CATALOGUE = get.ordena_Array(idxh
                             ,get.columna_fits(inroute_h + infile_catalogo2,0))
outfile_array_RA_h = get.ordena_Array(idxh
                             ,get.columna_fits(inroute_h + infile_catalogo2,1))
outfile_array_DEC_h = get.ordena_Array(idxh
                             ,get.columna_fits(inroute_h + infile_catalogo2,2)) 

outfile_array_S250=get.ordena_Array(idxh
                             ,get.columna_fits(inroute_h + infile_catalogo2,3))
outfile_array_S350=get.ordena_Array(idxh
                             ,get.columna_fits(inroute_h + infile_catalogo2,4))
outfile_array_S500=get.ordena_Array(idxh
                             ,get.columna_fits(inroute_h + infile_catalogo2,5))

from astropy.io import fits
#-------------------------- 
#   Creamos las columnas del fichero de salida.
#--------------------------  
col0 = fits.Column(name='idxg', format='14A', array=idg)
col1 = fits.Column(name='GAMA_IAU_ID', format='23A'
                   , array=outfile_array_GAMA_IAU_ID)
col2 = fits.Column(name='RA_G', format='D', array=outfile_array_RA_g 
                   , unit='Degrees')
col3 = fits.Column(name='DEC_G', format='D'
                   , array=outfile_array_DEC_g, unit='Degrees')

col4 = fits.Column(name='idxh', format='14A', array=idh)
col5 = fits.Column(name='HATLAS_DR1_CATALOGUE', format='23A'
                   , array=outfile_array_HATLAS_DR1_CATALOGUE)
col6 = fits.Column(name='RA_H', format='D'
                   , array=outfile_array_RA_h, unit='Degrees')
col7 = fits.Column(name='DEC_H', format='D'
                   , array=outfile_array_DEC_h, unit='Degrees')

col8 = fits.Column(name='S250', format='D'
                   , array=outfile_array_S250, unit='mJy')
col9 = fits.Column(name='S350', format='D'
                   , array=outfile_array_S350, unit='mJy')
col10 = fits.Column(name='S500', format='D'
                   , array=outfile_array_S500, unit='mJy')

col11 = fits.Column(name='E250', format='D'
, array=get.ordena_Array(idxh,get.columna_fits(inroute_h + infile_catalogo2,6))
, unit='mJy')
col12 = fits.Column(name='E350', format='D'
, array=get.ordena_Array(idxh,get.columna_fits(inroute_h + infile_catalogo2,7))
, unit='mJy')
col13 = fits.Column(name='E500', format='D'
, array=get.ordena_Array(idxh,get.columna_fits(inroute_h + infile_catalogo2,8))
, unit='mJy')

col14 = fits.Column(name='Dist. angular', format='D', array=d2d,unit='Degrees') 
#    Desplazamiento al rojo.
col15 = fits.Column(name='Z_H', format='D', array=zh)
col16 = fits.Column(name='Error_Z_H', format='D', array=szh)
col17 = fits.Column(name='flag_H', format='D', array=flag_h)

col18 = fits.Column(name='Z_G', format='D', array=zg)
col19 = fits.Column(name='Error_Z_G', format='D', array=szg)
col20 = fits.Column(name='flag_G', format='D', array=flag_g)

#    Factores de bayes.
col21 = fits.Column(name='Bayes_posicional', format='D'
                    , array=bayes_posicional) 
col22 = fits.Column(name='Bayes_z', format='D', array=bayes_z) 
col23 = fits.Column(name='Bayes_conjunto', format='D', array=bayes_conjunto) 

#    Ahora fabrico un objeto con esas columnas:
cols = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6, col7, col8, 
                     col9, col10, col11, col12, col13, col14, col15, col16, 
                     col17, col18, col19, col20, col21, col22, col23])

#    Ese objeto se mete en una tabla HDU. 
tbhdu = fits.BinTableHDU.from_columns(cols)
   
#   Y al final lo escribimos en un fichero "outfile" en la ruta "outroute".
tbhdu.writeto(outroute + outfile + '.fits',overwrite=True)

#==============================================================================
#   Representación de los fatores de Bayes.
#==============================================================================

# Devolvemos las distancias angulares del matching a arcsec.
d2d = d2d*u.radian
d2d = d2d.to(u.arcsec)
d2d = d2d.value

import matplotlib.pyplot as plt

#--------------------------
#   Factor de Bayes posicional vs separación angular.
#--------------------------
figura=plt.figure(num = None, figsize = (12, 6.4), dpi = 80, facecolor = 'w',
                  edgecolor = 'k')

plt.scatter(d2d, bayes_posicional, color="blue",linewidth=0.05, s=1)

plt.xlabel('${\phi}_{ij} [\mathrm{arcsec}]$', fontsize = 14, color='blue')
plt.ylabel('$B^{p}_{12}$', fontsize = 16, color='blue')

plt.axhline(y=1, linewidth=1, color='k', linestyle='--')

plt.yscale('log', linthreshx=0.1)

#   Cambiamos los límites de los ejes. 
plt.axis([0, 55, 0, 1.5e9])

plt.show()
figura.savefig(outroute + outfile + 'bayes_pos' + outfile, dpi=150
               , transparent=True)

print("\nDistancia angular para para la cual el factor de Bayes posicional "
   +"toma el valor 1, B_p=1, fi=",
   (np.sqrt(-2*(sg**2+sh**2)*(np.log((sg**2+sh**2)/2)))*u.radian).to(u.arcsec))

#--------------------------
#   Factor de Bayes espectroscópico vs zh-zg.
#--------------------------
# Se hace una copia del array original para no modificar sus valores.
import copy 
bayes_z_representacion = copy.copy(bayes_z)

# Para poder representar aquellos valores con factor de Bayes fotométrico 
# inferior a 1e-10.
for i in range(0,len(bayes_z),1):
    if np.any(bayes_z[i]<2e-20):
        bayes_z_representacion[i]=2e-20

zh_g=[]
for i in range(0,len(zh),1):
    zh_g.append(zh[i]-zg[i])
    
figura=plt.figure(num = None, figsize = (12, 6.4), dpi = 80, facecolor = 'w',
                  edgecolor = 'k')

plt.scatter(zh_g, bayes_z_representacion, color="blue", linewidth=0.1, s=2)

plt.grid(False)
plt.xlabel('$z_{h}-z_{g}$', fontsize = 16, color='blue')
plt.ylabel('$B^{z}_{12}$', fontsize = 16, color='blue')

plt.axhline(y=1, linewidth=1, color='k', linestyle='--')

plt.yscale('log', linthreshx=0.1)

#   Cambiamos los limites de los ejes. 
plt.axis([-1, 3.5, 1e-20, 15000])

plt.show()
figura.savefig(outroute + outfile + 'bayes_z' + outfile, dpi=150
               , transparent=True)

#==============================================================================
#   Representación posición de las observaciones, con los que se va a realizar
# el matching.
#==============================================================================

#--------------------------
#   G09 Y BLOQUE 2
#--------------------------
figura=plt.figure(num = None, figsize = (15,5.3), dpi = 800, facecolor = 'w'
                  , edgecolor = 'k')

plt.scatter(ra_c_HATLAS, dec_c_HATLAS, color="blue", linewidth=0.1, s=2,
            label="HATLAS")
plt.scatter(ra_c_GAMA, dec_c_GAMA, color="green", linewidth=0.1, s=2,
            label="GAMA")

plt.grid(False)
plt.xlabel('RA [degree]', fontsize = 16, color='blue')
plt.ylabel('DEC [degree]', fontsize = 16, color='blue')

#   Límites de los ejes. 
plt.axis([127, 142, -2.2, 3.1])

plt.show()
figura.savefig(outroute + outfile + 'region_1', dpi=150, transparent=True)

[area_hatlas,area_gama,area_interseccion]=get.area_region(127,142,-2.1,3.3
,ra_c_GAMA,dec_c_GAMA,ra_c_HATLAS,dec_c_HATLAS,0.08)

print(area_hatlas,area_gama,area_interseccion)

#--------------------------
#   G12 Y BLOQUE 3
#--------------------------
figura=plt.figure(num = None, figsize = (15,5.3), dpi = 800, facecolor = 'w'
                  , edgecolor = 'k')

plt.scatter(ra_c_HATLAS, dec_c_HATLAS, color="blue", linewidth=0.1, s=2,
            label="HATLAS")
plt.scatter(ra_catalogo_GAMA, dec_catalogo_GAMA, color="green", linewidth=0.1,
            s=2, label="GAMA")

plt.grid(False)
plt.xlabel('RA [degree]', fontsize = 16, color='blue')
plt.ylabel('DEC [degree]', fontsize = 16, color='blue')

#   Limites de los ejes. 
plt.axis([172, 187, -3.1, 2.2])

plt.show()
figura.savefig(outroute + outfile + 'region_2', dpi=150, transparent=True)

[area_hatlas,area_gama,area_interseccion]=get.area_region(172, 187, -3.1, 2.2
,ra_c_GAMA,dec_c_GAMA,ra_c_HATLAS,dec_c_HATLAS,0.08)

print(area_hatlas,area_gama,area_interseccion)

#--------------------------
#   G15 Y BLOQUE 4
#--------------------------
figura=plt.figure(num = None, figsize = (15,5.3), dpi = 800, facecolor = 'w'
                  , edgecolor = 'k')

plt.scatter(ra_c_HATLAS, dec_c_HATLAS, color="blue", linewidth=0.1, s=2,
            label="HATLAS")
plt.scatter(ra_c_GAMA, dec_c_GAMA, color="green", linewidth=0.1, s=2, 
            label="GAMA")

plt.grid(False)
plt.xlabel('RA [degree]', fontsize = 16, color='blue')
plt.ylabel('DEC [degree]', fontsize = 16, color='blue')

#   Limites de los ejes. 
plt.axis([210, 225, -2.2, 3.1])

plt.show()
figura.savefig(outroute + outfile + 'region_3', dpi=150, transparent=True)

[area_hatlas,area_gama,area_interseccion]=get.area_region(210, 225, -2.2, 3.1
,ra_c_GAMA,dec_c_GAMA,ra_c_HATLAS,dec_c_HATLAS,0.08)

print(area_hatlas,area_gama,area_interseccion)           

#==============================================================================
#   Número de HALOS en HATLAS_DR1_CATALOGUE_V1
#==============================================================================

s250=get.columna_fits(inroute_h + infile_catalogo2,3)
s350=get.columna_fits(inroute_h + infile_catalogo2,4)
s500=get.columna_fits(inroute_h + infile_catalogo2,5)

print("\nNúmero de HALOs en HATLAS_DR1 es de: (Negrello,Gonzalez,Ambos)="
                                       , get.numero_halos(idxh,s250,s350,s500))

#==============================================================================
#   Histogramas de los redshifts
#==============================================================================

from astropy.visualization import hist as hist_astropy
#--------------------------
#   Histograma redshits catálogo GAMA con los que se va a realizar el matching.
#--------------------------
idg_pos_o=get.id_pos(z_g)
flag_g_o=get.columna_fits(inroute_g + infile_catalogo1,7)

# _amg (antes matching gama).
[id_spec11_amg,id_spec12_amg,id_spec13_amg,id_disabled_amg
                     ,idx_amg]=get.obj_by_flag_gama(idg_pos_o,flag_g_o,"array")

z_spec11_amg=get.ordena_Array(id_spec11_amg, z_g)
z_spec12_amg=get.ordena_Array(id_spec12_amg, z_g)
z_spec13_amg=get.ordena_Array(id_spec13_amg, z_g)
z_amg=get.ordena_Array(idx_amg, z_g)

figura=plt.figure(num = None, figsize = (9, 6), dpi = 80, facecolor = 'w'
                  , edgecolor = 'k')

#   "knuth"     "scott"     "freedman"      "blocks"
hist_astropy(z_spec12_amg, bins=100, range=(0,0.6), normed=False, alpha=0.6
        , color='#2e7d32',  histtype='stepfilled', label='$z_{g}$ Q=2')
hist_astropy(z_spec11_amg, bins=100, range=(0,0.6), normed=False, alpha=0.6
        , color='#558b2f',  histtype='stepfilled', label='$z_{g}$ Q=1')
hist_astropy(z_spec13_amg, bins=100, range=(0,0.6), normed=False, alpha=0.6
        , color='#9e9d24',  histtype='stepfilled', label='$z_{g}$ Q=3')
hist_astropy(z_amg, bins=100, range=(0,0.6), normed=False, alpha=0.7
        , color='black',  histtype='step', label='$z_{g}$')

plt.xlabel('$z_{g}$', fontsize = 16)
plt.ylabel('$\mathrm{N}(z_{g})$', fontsize = 16)

plt.legend(loc='best')

figura.savefig(outroute + outfile +  'histograma_z_gama_disponible.png'
               , dpi=150, transparent=True)
plt.show()

print("(Q=1,Q=2,Q=3,disabled,total) =",get.obj_by_flag_gama(idg_pos_o,
                                                            flag_g_o,"numero"))
#--------------------------
#   Histograma redshits catálogo GAMA que forman parte del matching.
#--------------------------
# _dmg (despues matching gama).
[id_spec11_dmg,id_spec12_dmg,id_spec13_dmg,id_disabled_dmg
                            ,idx_dmg]=get.obj_by_flag_gama(idxg,flag_g,"array")

z_spec11_dmg=get.ordena_Array(id_spec11_dmg, z_g)
z_spec12_dmg=get.ordena_Array(id_spec12_dmg, z_g)
z_spec13_dmg=get.ordena_Array(id_spec13_dmg, z_g)
z_dmg=get.ordena_Array(idx_dmg, z_g)

figura=plt.figure(num = None, figsize = (9, 6), dpi = 80, facecolor = 'w'
                  , edgecolor = 'k')

#   "knuth"     "scott"     "freedman"      "blocks"
hist_astropy(z_spec12_dmg, bins=100, range=(0,0.6), normed=False, alpha=0.6
             , color='#2e7d32',  histtype='stepfilled', label='$z_{g}$ Q=2')
hist_astropy(z_spec11_dmg, bins=100, range=(0,0.6), normed=False, alpha=0.6
             , color='#558b2f',  histtype='stepfilled', label='$z_{g}$ Q=1')
hist_astropy(z_spec13_dmg, bins=100, range=(0,0.6), normed=False, alpha=0.6
             , color='#9e9d24',  histtype='stepfilled', label='$z_{g}$ Q=3')
hist_astropy(z_dmg, bins=100, range=(0,0.6), normed=False, alpha=0.7
             , color='black',  histtype='step', label='$z_{g}$')

plt.xlabel('$z_{g}$', fontsize = 16)
plt.ylabel('$\mathrm{N}(z_{g})$', fontsize = 16)

plt.legend(loc='best')

figura.savefig(outroute + outfile + 'histograma_z_gama_matching.png', dpi=150
               , transparent=True)
plt.show()

print("(Q=1,Q=2,Q=3,disabled,total) =",get.obj_by_flag_gama(idxg,flag_g,
                                                                     "numero"))
#--------------------------
#   Histograma redshits catálogo HATLAS con los que se va a realizar el 
# matching.
#--------------------------
idh_pos=get.id_pos(z_h)
flag_h_o=get.columna_fits(inroute_h + infile_catalogo2,16)

# _amh (antes matching hatlas)
[id_spec_amh,id_annz_amh,id_phot_amh,id_disabled_amh
                     ,idx_amh]=get.obj_by_flag_hatlas(idh_pos,flag_h_o,"array")

z_spec_amh=get.ordena_Array(id_spec_amh, z_h)
z_annz_amh=get.ordena_Array(id_annz_amh, z_h)
z_phot_amh=get.ordena_Array(id_phot_amh, z_h)
z_amh=get.ordena_Array(idx_amh, z_h)

figura=plt.figure(num = None, figsize = (9, 6), dpi = 80, facecolor = 'w'
                  , edgecolor = 'k')

#   "knuth"     "scott"     "freedman"      "blocks"
hist_astropy(z_spec_amh, bins=100, range=(0,3.5), normed=False, alpha=0.6
    , color='green',  histtype='stepfilled', label='$z_{h}^{\mathrm{\;spec}}$')
hist_astropy(z_phot_amh, bins=100, range=(0,3.5), normed=False, alpha=0.3
    , color='red',  histtype='stepfilled', label='$z_{h}^{\mathrm{\;phot}}$')
hist_astropy(z_annz_amh, bins=100, range=(0,3.5), normed=False, alpha=0.6
    , color='blue',  histtype='stepfilled', label='$z_{h}^{\mathrm{\;ANNZ}}$')
hist_astropy(z_amh, bins=100, range=(0,3.5), normed=False, alpha=0.7
    , color='black',  histtype='step', label='${z}_{h}$')

plt.xlabel('$z_{h}$', fontsize = 16)
plt.ylabel('$\mathrm{N}(z_{h})$', fontsize = 16)

plt.legend(loc='best')

figura.savefig(outroute + outfile + 'histograma_z_hatlas_disponible.png'
               , dpi=150, transparent=True)
plt.show()

print("(z_spec,z_ANNZ,z_phot,disabled,total) =",get.obj_by_flag_hatlas(idh_pos,
      flag_h_o,"numero"))

#--------------------------
#   Histograma redshits catálogo HALTAS que forman parte del matching.
#--------------------------
# _dmh (después matching haltlas)
[id_spec_dmh,id_annz_dmh,id_phot_dmh,id_disabled_dmh
                          ,idx_dmh]=get.obj_by_flag_hatlas(idxh,flag_h,"array")

z_spec_dmh=get.ordena_Array(id_spec_dmh, z_h)
z_annz_dmh=get.ordena_Array(id_annz_dmh, z_h)
z_phot_dmh=get.ordena_Array(id_phot_dmh, z_h)
z_dmh=get.ordena_Array(idx_dmh, z_h)

figura=plt.figure(num = None, figsize = (9, 6), dpi = 80, facecolor = 'w'
                  , edgecolor = 'k')

#   "knuth"     "scott"     "freedman"      "blocks"
hist_astropy(z_spec_dmh, bins=100, range=(0,3.5), normed=False, alpha=0.6
    , color='green',  histtype='stepfilled', label='$z_{h}^{\mathrm{\;spec}}$')
hist_astropy(z_phot_dmh, bins=100, range=(0,3.5), normed=False, alpha=0.3
    , color='red',  histtype='stepfilled', label='$z_{h}^{\mathrm{\;phot}}$')
hist_astropy(z_annz_dmh, bins=100, range=(0,3.5), normed=False, alpha=0.6
    , color='blue',  histtype='stepfilled', label='$z_{h}^{\mathrm{\;ANNZ}}$')
hist_astropy(z_dmh, bins=100, range=(0,3.5), normed=False, alpha=0.7
    , color='black',  histtype='step', label='${z}_{h}$')

plt.xlabel('$z_{h}$', fontsize = 16)
plt.ylabel('$\mathrm{N}(z_{h})$', fontsize = 16)

plt.legend(loc='best')

figura.savefig(outroute + outfile + 'histograma_z_hatlas_matching.png'
               , dpi=150, transparent=True)
plt.show()

print("(z_spec,z_ANNZ,z_phot,disabled,total) =",get.obj_by_flag_hatlas(idxh,
      flag_h,"numero"))

#==============================================================================
#   Factor de Bayes fotométrico vs posicional .
#==============================================================================

#   Se puede incluir en un fichero el nombre de aquellos objetos de Hatlas de 
# los que queremos destacar el pareado en el que participan por algún motivo. 
# Hay que poner el nombre completo; Ej: "HATLAS J083945.1+023440" 
#
# En principio ninguno => fichero vacío.

fichero_con_nombres=[]
with open(inroute_h + infile_objetos_interesantes,'r') as fichero_nombres:
    fichero_con_nombres=fichero_nombres.readlines()

nombres_resaltar=[]
for o in range (0,len(fichero_con_nombres),1):
    #   Parte cada string del array en dos, por donde esta el elemento.
    # que le pasamos en split.
    nombres_confirmados_consaltolinea=fichero_con_nombres[o].split('\n')
    nombres_resaltar.append(nombres_confirmados_consaltolinea[0])

#    Necesitmos estos arrays para la representación de dos rectas en las dos
# figuras siguientes.
x = np.arange(-10, 11, 0.1)
y1 = -1*x
y100 = -1*x+2

#-------------------------- 
#   Representación de aquellos emparejamientos en las que participa una 
# observación de HATLAS cuyo z ha sido obtenido mediante "z_phot_hatlas".
#--------------------------
    
[x_0, y_0, x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4, x_5, y_5] = \
    get.representacion_candidatos(bayes_posicional,bayes_z, bayes_conjunto,
    outfile_array_HATLAS_DR1_CATALOGUE,nombres_resaltar,outfile_array_S250,
    outfile_array_S350,outfile_array_S500,flag_h,zh,zg,"ajuste")

figura=plt.figure(num = None, figsize = (12,3), dpi = 800, facecolor = 'w'
                  , edgecolor = 'k')

plt.scatter(x_3, y_3, color="#ffea00", linewidth=0.1, s=10,
            label="Resto de objetos")
plt.scatter(x_4, y_4, color="#64dd17", linewidth=0.1, s=10,
            label="Bayes mayor que 100")
plt.scatter(x_5, y_5, color="cyan", linewidth=0.1, s=10,
            label="Bayes mayor que 1")
plt.scatter(x_0, y_0, color="black", linewidth=0.1, s=40,
            label="Objetos del fichero")
plt.scatter(x_1, y_1, color="red", linewidth=0.1, s=30,
            label="Negrello et al.")
plt.scatter(x_2, y_2, color="blue", linewidth=0.1, s=30,
            label="Gonzalez-Nuevo et al.")
plt.axvline(x=0, linewidth=1, color='k')
plt.axhline(y=0, linewidth=1, color='k')

plt.plot(x, y1, color='k', linestyle='-.')
plt.plot(x, y100,'k--')

plt.grid(False)
plt.ylabel('$\mathrm{log_{10}}{({B}^{z}_{12})}$', fontsize = 16, color='blue')

#   Limites de los ejes.   [xmin,xmax,ymin,ymax]
plt.axis([-1.8, 9.25, -30, -1])

plt.yscale('symlog', linthreshx=0.1)

plt.show()
figura.savefig(outroute + outfile + 'factor_bayes_log_semilog_ajuste', dpi=150,
               transparent=True)

#-------------------------- 
#   Todas las observaciones
#--------------------------
[x_0, y_0, x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4, x_5
 , y_5] = get.representacion_candidatos(bayes_posicional,bayes_z,bayes_conjunto
 ,outfile_array_HATLAS_DR1_CATALOGUE,nombres_resaltar,outfile_array_S250
 ,outfile_array_S350,outfile_array_S500,flag_h,zh,zg,"todos")

figura=plt.figure(num = None, figsize = (12,9), dpi = 800, facecolor = 'w',
                  edgecolor = 'k')

plt.scatter(x_3, y_3, color="#ffea00", linewidth=0.1, s=10,
            label="Resto de objetos")
plt.scatter(x_4, y_4, color="#64dd17", linewidth=0.1, s=10,
            label="Bayes mayor que 100")
plt.scatter(x_5, y_5, color="cyan", linewidth=0.1, s=10,
            label="Bayes mayor que 1")
plt.scatter(x_0, y_0, color="black", linewidth=0.1, s=40,
            label="Objetos del fichero")
plt.scatter(x_1, y_1, color="red", linewidth=0.1, s=30,
            label="Negrello et al.")
plt.scatter(x_2, y_2, color="blue", linewidth=0.1, s=30,
            label="Gonzalez-Nuevo et al.")
plt.axvline(x=0, linewidth=1, color='k')
plt.axhline(y=0, linewidth=1, color='k')

plt.plot(x, y1, color='k', linestyle='-.')
plt.plot(x, y100,'k--')

plt.grid(False)
plt.xlabel('$\mathrm{log_{10}}({B}^{p}_{12})$', fontsize = 16, color='blue')
plt.ylabel('$\mathrm{log_{10}}{({B}^{z}_{12})}$', fontsize = 16, color='blue')
#   Limites de los ejes.   [xmin,xmax,ymin,ymax]
plt.axis([-1.8, 9.25, -110, 1.1])

plt.yscale('symlog', linthreshx=0.1)

plt.show()
figura.savefig(outroute + outfile + 'factor_bayes_log_semilog', dpi=150
               , transparent=True)

[id_halo_h_espec,id_halo_h_fot,id_halo_negrello_gonzalez]= \
    get.selecion_candidatos(outfile_array_S250,outfile_array_S350,
    outfile_array_S500,flag_h,zh,zg,bayes_posicional, bayes_conjunto)
    
print("Hay",len(id_halo_h_espec),len(id_halo_h_fot),
     "candidatos segun nuestro criterio y",len(id_halo_negrello_gonzalez),
     "cumplen ambos criterios")

#==============================================================================
#   Selección de los candidatos a lente gravitatoria .
#==============================================================================

#   Identificador de todos los emparejamientos que cumplen nuestro criterio
id_halo=id_halo_h_espec+id_halo_h_fot

idg_halo=get.ordena_Array(id_halo,idg)
idh_halo=get.ordena_Array(id_halo,idh)

#-------------------------- 
#   Creamos los arrays que faltan del fichero de salida
#--------------------------      
nombre_gama_halo=get.ordena_Array(id_halo,outfile_array_GAMA_IAU_ID)
nombre_hatlas_halo=get.ordena_Array(id_halo,outfile_array_HATLAS_DR1_CATALOGUE)

z_h_halo=get.ordena_Array(id_halo,zh)
ez_h_halo=get.ordena_Array(id_halo,szh)
flag_h_halo=get.ordena_Array(id_halo,flag_h)

# Los flujos se encuantran en mJy, 
s250_halo=get.ordena_Array(id_halo,outfile_array_S250)*u.Jy
s350_halo=get.ordena_Array(id_halo,outfile_array_S350)*u.Jy
s500_halo=get.ordena_Array(id_halo,outfile_array_S500)*u.Jy
s250_halo=s250_halo.to(u.mJy)
s350_halo=s350_halo.to(u.mJy)
s500_halo=s500_halo.to(u.mJy)

z_g_halo=get.ordena_Array(id_halo,zg)
ez_g_halo=get.ordena_Array(id_halo,szg)
flag_g_halo=get.ordena_Array(id_halo,flag_g)

bayes_posicional_halo=get.ordena_Array(id_halo,bayes_posicional)
bayes_z_halo=get.ordena_Array(id_halo,bayes_z)
bayes_conjunto_halo=get.ordena_Array(id_halo,bayes_conjunto)

#-------------------------- 
#   Creamos las columnas del fichero de salida.
#--------------------------  
cl0 = fits.Column(name='idxg', format='14A', array=idg_halo)
cl1 = fits.Column(name='GAMA_IAU_ID', format='23A', array=nombre_gama_halo)
cl2 = fits.Column(name='idxh', format='14A', array=idh_halo)
cl3 = fits.Column(name='HATLAS_DR1_CATALOGUE', format='23A'
                   , array=nombre_hatlas_halo)

cl4 = fits.Column(name='S250', format='D', array=s250_halo, unit='mJy')
cl5 = fits.Column(name='S350', format='D', array=s350_halo, unit='mJy')
cl6 = fits.Column(name='S500', format='D', array=s500_halo, unit='mJy')

#    Desplazamiento al rojo.
cl7 = fits.Column(name='Z_H', format='D', array=z_h_halo)
cl8 = fits.Column(name='Error_Z_H', format='D', array=ez_h_halo)
cl9 = fits.Column(name='flag_H', format='D', array=flag_h_halo)

cl10 = fits.Column(name='Z_G', format='D', array=z_g_halo)
cl11 = fits.Column(name='Error_Z_G', format='D', array=ez_g_halo)
cl12 = fits.Column(name='flag_G', format='D', array=flag_g_halo)

#    Factores de bayes.
cl13 = fits.Column(name='Bayes_posicional', format='D'
                    , array=bayes_posicional_halo) 
cl14 = fits.Column(name='Bayes_z', format='D', array=bayes_z_halo) 
cl15 = fits.Column(name='Bayes_conjunto', format='D'
                    , array=bayes_conjunto_halo) 

cls = fits.ColDefs([cl0, cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8, cl9, cl10
                        , cl11, cl12, cl13, cl14, cl15])

tbhdu = fits.BinTableHDU.from_columns(cls)
   
tbhdu.writeto(outroute + outfile + '_candidatos.fits',overwrite=True)