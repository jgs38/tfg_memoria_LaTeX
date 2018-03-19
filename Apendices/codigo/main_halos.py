# -*- coding: utf-8 -*-

"""
Calcula los redshifts de los objetos que aparecen en la tabla del artículo 
Gonzalez-Nuevo et al 2012 y compara estos valores con los que obtuvieron estos
autores. Por último se obtiene un ajuste con el que se calcularán los
errores asociados a las medidas obtenidas con las funciones definidas en el 
modulo "rojo".

  Python 3.6.0    scipy 0.18   matplotlib 2.0.0     numpy 1.11.3    astropy 1.3

@author: Javier Gutiérrez Solórzano
"""
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
outroute = './outroute/'
#--------------------------
#   Nombre de los ficheros.
#--------------------------
infile_template = 'SMM_template_norm.sed'
infile_catalogo = "catalogo_paper_joaquin_2012.fits"

outfile = "comparacion_Z.fits"
outpicture = "comparacion_Z.png"

#==============================================================================
#   Desplazamiento al rojo de los objetos H-Atlas.
#============================================================================== 
#   Soporte de unidades físicas.
from astropy import units as u 
#   Soporte de magnitudes físicas. Solo se utiliza la velocidad de la luz.
from astropy.constants import c 

#   Importamos el modulo get donde se definen las funciones para las funciones 
# utilizadas para la lectura de los ficheros y tablas.
import get

#-------------------------- 
#   #  Entrada:
#   Leemos el fichero de la SED modelo. El fichero contiene dos columnas; la 
# primera representa las longitudes de onda en Ångström y la segunda valores
# valores de la densidad espectral de flujo (F_(lambda)), en unidades de: 
# (erg)* (cm)**(-2)* (s)**(-1)* (Å)**(-1)
#-------------------------- 
[x_sed, y_sed] = get.sed_func(inroute + infile_template) 
#   Unidades a las longitudes de onda
x_sed=x_sed*u.AA # En Ångström
#   Unidades de densidad espectral de flujo (F_(lambda))
y_sed=y_sed*((u.erg)* (u.cm)**(-2)* (u.s)**(-1)* (u.AA)**(-1))
#   Cambio de F_(lambda) a F(ípsilon) multiplicando por (lambda**2)/c 
y_sed=y_sed*(x_sed**2)*(c**-1)
#   Dejo los arrays con las unidades que tienen los valores experimentales que 
# voy a utilizar
x_sed=x_sed.to(u.micron)
y_sed=y_sed.to(u.mJy) 
 
#-------------------------- 
#   Obtenemos los valores experimentales (mediciones realizadas por HATLAS) de
# todos los objetos de un catálogo a partir de la lectura de un fichero con 
# formato .fits
#-------------------------- 
x_experimental = [250,350,500]*u.micron 
[y_experimental, y_error]= get.datos_experimentales(inroute + infile_catalogo)

y_experimental = y_experimental*u.mJy
y_error = y_error*u.mJy

#-------------------------- 
#   Salida
#-------------------------- 
import numpy as np

import rojo
[z,E_z] = rojo.z_phot_hatlas(x_sed,y_sed,x_experimental,y_experimental,y_error)

#==============================================================================
#   Creación del fichero fits de 'salida'
#==============================================================================
#   Creamos varios arrays que contendrán los elementos de las otras filas 
# del fichero de salida.
outfile_array_HATLAS_DR1_CATALOGUE = get.columna_fits(inroute + 
                                                      infile_catalogo,0)

outfile_array_RA_h = get.columna_fits(inroute + infile_catalogo,2)
outfile_array_DEC_h = get.columna_fits(inroute + infile_catalogo,3)

outfile_array_S250 = get.columna_fits(inroute + infile_catalogo,4)
outfile_array_S350 = get.columna_fits(inroute + infile_catalogo,5)
outfile_array_S500 = get.columna_fits(inroute + infile_catalogo,6)

outfile_array_e_S250 = get.columna_fits(inroute + infile_catalogo,7)
outfile_array_e_S350 = get.columna_fits(inroute + infile_catalogo,8)
outfile_array_e_S500 = get.columna_fits(inroute + infile_catalogo,9)
#   outfile_array_z es un array con los valores de z que se encuentran en el 
# articulo.
outfile_array_z = get.columna_fits(inroute + infile_catalogo,10)
outfile_array_e_z = get.columna_fits(inroute + infile_catalogo,11)

from astropy.io import fits
#-------------------------- 
#   Creamos las columnas del fichero de salida.
#-------------------------- 
col0 = fits.Column(name='HATLAS_DR1_CATALOGUE', format='23A',
                   array=outfile_array_HATLAS_DR1_CATALOGUE)
col1 = fits.Column(name='RA_H', format='D',  array=outfile_array_RA_h,
                   unit='Degrees')
col2 = fits.Column(name='DEC_H', format='D', array=outfile_array_DEC_h,
                   unit='Degrees')

col3 = fits.Column(name='S500', format='D', array=outfile_array_S500, 
                   unit='Jy')
col4 = fits.Column(name='e_S500', format='D', array=outfile_array_e_S500,
                   unit='Jy')

col5 = fits.Column(name='S350', format='D', array=outfile_array_S350,
                   unit='Jy')
col6 = fits.Column(name='e_S350', format='D', array=outfile_array_e_S350,
                   unit='Jy')

col7 = fits.Column(name='S250', format='D', array=outfile_array_S250,
                   unit='Jy')
col8 = fits.Column(name='e_S250', format='D', array=outfile_array_e_S250,
                   unit='Jy')

    #   Desplazamiento al rojo
col9 = fits.Column(name='Z_ajuste', format='D', array=outfile_array_z)
col10 = fits.Column(name='Error_Z_ajuste', format='D', array=outfile_array_e_z)

col11 = fits.Column(name='Z_HALOS', format='D', array=z)
col12 = fits.Column(name='Error_Z_HALOS', format='D', array=E_z)
    #   Ahora fabrico un objeto con esas columnas:
cols = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6, col7, col8, 
                     col9, col10, col11, col12])
    #   Ese objeto se utiliza para crear la tabla HDU.
tbhdu = fits.BinTableHDU.from_columns(cols)   
    #   Y al final lo escribimos en un fichero 'outfile' en la ruta 'outroute'.
tbhdu.writeto(outroute + outfile,overwrite=True)
#-------------------------- 
#   Revisión de datos; z, E_z, outfile_array_z, outfile_array_e_z .
#--------------------------  
print('\n Longitud del array z=',len(z),'\n',
'\n z =  ',z,'\n',
'\n E_z =',E_z,'\n',
'\n Longitud del array outfile_array_z=',len(outfile_array_z),'\n',
'\n outfile_array_z =  ',outfile_array_z,'\n',
'\n outfile_array_e_z =',outfile_array_e_z)

#==============================================================================
#   Graficas
#==============================================================================
import matplotlib.pyplot as plt
#-------------------------- 
#   Graficas de los ajustes.
#-------------------------- 
for i in range(6,10,1): # 3,10
     K=rojo.z_phot_objeto(x_sed, y_sed, x_experimental, y_experimental[i], 
                   y_error[i])[0][0]
     C=rojo.z_phot_objeto(x_sed, y_sed, x_experimental, y_experimental[i],
                   y_error[i])[0][1]
    
     figura=plt.figure()

     plt.loglog(x_sed,y_sed, color="blue", linewidth=1.5, linestyle="--",
                label="SMM J2")
     plt.loglog(x_sed*K,y_sed*C, color="red")
     plt.loglog(x_experimental,y_experimental[i],'+', color='black')
    
     plt.axvline(500*(1+z[i]), linewidth=0.75, color = 'yellow')
     plt.axvline(500, color = 'green',linewidth=0.75, linestyle="--")
    
     plt.grid(False)
     
     plt.title(outfile_array_HATLAS_DR1_CATALOGUE[i], fontsize = 10, 
               color='blue')
     plt.xlabel('$ \mathrm{\lambda}\;[\mu \mathrm{m}]$', fontsize 
                           = 10, color='blue')
     plt.ylabel('$ {\mathrm{F}}_{ν}\;[ \mathrm{mJy} ]$',
                fontsize = 10, color='blue')
     #  Transforma 'i' en un string para utilizarlo en el nombre de la figura.
     indice = str(i)
     #------------------------------------------------------------------------
     #  a)
     figura.savefig(outroute + 'ajuste_' + indice, dpi=150, transparent=True)
     #-------------------------- 
#     #   b)
#     #   Cambiamos los limites de los ejes para hacer ampliar la zona en la 
#     # que los puntos experimentales se encuentran mas proximos a la curva de
#     # ajuste.
#     #      [xmin,xmax,ymin,ymax]
#     plt.axis([200, 700, 10, 500])
#     figura.savefig(outroute + 'ajuste_' + indice + '_zoom', dpi=150,
#                    transparent=True)     
     #-------------------------------------------------------------------------
     # Muestra la figura.
     plt.show()

#-------------------------- 
#   Grafica de Z_articulo vs Z_programa
#--------------------------
#   Transforma un tipo 'list' a 'array' 
#   Ajuste, escribe en la variable f los parámetros del ajuste y la matriz de 
# covarianza de estos parámetros.
f =  np.polyfit(outfile_array_z, z, 1, cov=True)   

C = f[1]  # matriz de covarianza.
p = f[0]  # valores del ajuste.

afit = p[0]   # pendiente ajustada.
bfit = p[1]   # ord. origen ajustada.
#--------------------------
def fsup(x,Cmat,nsigmas):  #   función que calcula el extremo superior de la
                           # banda, a NSIGMAS sigmas.
    ymean = afit*x+bfit
    sa    = Cmat[0,0]
    sb    = Cmat[1,1]
    sab   = Cmat[0,1]
    y     = ymean+nsigmas*np.sqrt(x*x*sa+sb+2*x*sab)
    return y
 
def finf(x,Cmat,nsigmas):  #    función que calcula el extremo inferior de la 
                           # banda, a NSIGMAS sigmas.
    ymean = afit*x+bfit
    sa    = Cmat[0,0]
    sb    = Cmat[1,1]
    sab   = Cmat[0,1]
    y     = ymean-nsigmas*np.sqrt(x*x*sa+sb+2*x*sab)
    return y
#-------------------------- 
#   Figura z_articulo vs z_programa.
#-------------------------- 
figura = plt.figure(num = None, figsize = (9, 9), dpi = 80, facecolor = 'w', 
                    edgecolor = 'k')
#   Cambiamos los limites de los ejes.
plt.axis([1, 3.5, 1, 3.5])
#   Genera 100 valores igualmente especiados entre 1 y 3.5.
x = np.linspace(1, 3.5, 100, endpoint=True)
#   Recta del ajuste.
plt.plot(x,afit*x+bfit,'r',linewidth=0.5)
#   Sombreado amarillo para la banda 1 sigma.
plt.fill_between(x,finf(x,C,2),fsup(x,C,2),facecolor='yellow', alpha=0.8)
#   Sombreado amarillo para la banda 2 sigma.
plt.fill_between(x,finf(x,C,1),fsup(x,C,1),facecolor='orange', alpha=0.8)
#   Dibuja los puntos y sus barras de error.
plt.errorbar(outfile_array_z, z, fmt='+', label="data", linewidth=0.15, 
             xerr=outfile_array_e_z, yerr=E_z, ecolor='black')
#   Sin rejilla de fondo.
plt.grid(False)
#   Etiqueta eje x.
plt.xlabel('${z}_{\mathrm{H}}$', fontsize = 16, color='blue')
#   Etiqueta eje y.
plt.ylabel('${z}_{\mathrm{phot}}$', fontsize = 16, color='blue')
#   Muestra la figura.
plt.show()
#   Guarda una imagen de la gráfica.
figura.savefig(outroute + outpicture, dpi=150, transparent=True)
#-------------------------- 
#   Muestra datos del ajuste por pantalla.
#--------------------------   

print('z_HALOS vs z_ajuste:','\n',
'--- Best fit    y = ({0}'.format(p[0]), chr(177) ,
'{0})'.format(np.sqrt(C[0,0])), chr(215),
'x + ({0}'.format(p[1]), chr(177) ,'{0})'.format(np.sqrt(C[1,1])))

#-------------------------- 
#   Grafica errores_z_articulo vs z_articulo.
#--------------------------  
from scipy.optimize import curve_fit

def func(x, A):
    """
        Modelo para nuestros datos. 
    
    """
    return A*(1+x)

#   pcov - Matriz covarianza.
#   popt - Vector con valores de ajuste.
popt_errores, pcov_errores = curve_fit(func, outfile_array_z
                                       , outfile_array_e_z)

def recta_sup(x,opt_par,cov_matrix,nsigmas):  #   función que calcula el 
                              # extremo superior de la banda, a NSIGMAS sigmas.
    ymean = opt_par[0]*(x+1)
    sa    = cov_matrix[0,0]
    y     = ymean+nsigmas*np.sqrt(sa)
    return y
 
def recta_inf(x,opt_par,cov_matrix,nsigmas):  #    función que calcula el 
                              # extremo inferior de la banda, a NSIGMAS sigmas.
    ymean = opt_par[0]*(x+1)
    sa    = cov_matrix[0,0]
    y     = ymean-nsigmas*np.sqrt(sa)
    return y

figura = plt.figure(num = None, figsize = (9, 6), dpi = 80, facecolor = 'w', 
                    edgecolor = 'k')
#   Cambiamos los limites de los ejes.
plt.axis([1.25, 3.1, 0.2, 0.7])
#   Genera 100 valores igualmente especiados entre 1 y 3.5
x = np.linspace(1, 3.2, 10, endpoint=True)
#   Recta del ajuste.
plt.plot(x,popt_errores[0]*(1+x),'r',linewidth=0.5)
#   Sombreado amarillo para la banda 1 sigma.
plt.fill_between(x,recta_inf(x,popt_errores,pcov_errores,2)
       ,recta_sup(x,popt_errores,pcov_errores,2),facecolor='yellow', alpha=0.8)
#   Sombreado amarillo para la banda 2 sigma.
plt.fill_between(x,recta_inf(x,popt_errores,pcov_errores,1)
       ,recta_sup(x,popt_errores,pcov_errores,1),facecolor='orange', alpha=0.8)
#   Dibuja los puntos.
plt.errorbar(outfile_array_z, outfile_array_e_z, fmt='+', label="data"
       , linewidth=0.15, ecolor='black')
#   Sin rejilla de fondo.
plt.grid(False)
#   Etiqueta eje x.
plt.xlabel('${z}_{\mathrm{H}}$', fontsize = 16, color='blue')
#   Etiqueta eje y.
plt.ylabel('${\sigma}^{z}$', fontsize = 16, color='blue')
#   Muestra la figura.
plt.show()
#   Guarda una imagen de la gráfica.
figura.savefig(outroute + "ajuste_errores.png", dpi=150, transparent=True)
#-------------------------- 
#   Muestra datos del ajuste por pantalla.
#--------------------------   

print('Errores z_HALOS vs z_HALOS:\n',
'--- Best fit    y = ({0}'.format(popt_errores[0]), chr(177) ,
'{0})'.format(np.sqrt(pcov_errores[0,0])), chr(215),'(1+x)')