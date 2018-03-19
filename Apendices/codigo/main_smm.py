# -*- coding: utf-8 -*-

"""
Representación gráfica de la SED de la galaxia SMM expresado en F_lambda y
 F_ípsilon .

  Python 3.6.0    matplotlib 2.0.0   astropy 1.3
  
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
inroute = './inroute/'
outroute = './outroute/'
#--------------------------
#   Nombre del fichero de entrada
#--------------------------
infile_template = 'SMM_template_norm.sed'

#==============================================================================
#   Desplazamiento al rojo de los objetos H-Atlas
#============================================================================== 
import get
#-------------------------- 
#   Leemos el fichero de la SED modelo. El fichero contiene dos columnas; la 
# primera representa las longitudes de onda en Ångström y la segunda valores
# valores de la densidad espectral de flujo (F_(lambda)), en unidades de: 
# (erg)* (cm)**(-2)* (s)**(-1)* (Å)**(-1) aunque normalizado en 5570 Ångström.
#-------------------------- 
[x_sed, y_sed] = get.sed_func(inroute + infile_template) 

#==============================================================================
#   Grafica F_(lambda)
#==============================================================================
import matplotlib.pyplot as plt

figura=plt.figure()
plt.loglog(x_sed,y_sed, color="blue", linewidth=1, linestyle="-",
           label="SMM J2")
plt.grid(False)
plt.title('SMM J2135-0102', fontsize = 16, color='blue')
plt.errorbar(x_sed[42],y_sed[42], fmt='+', label="data", linewidth=1, xerr=0, 
             yerr=0, ecolor='yellow')
plt.xlabel('$ \lambda\;[\mathrm{\AA}]$', fontsize = 16, color='blue')
plt.ylabel('$ \mathrm{F_{\lambda}}\; [@\,5570\,\mathrm{\AA}]$', fontsize = 16, 
           color='blue')
figura.savefig(outroute + 'grafica_SMM_F_lambda', dpi=150, transparent=True)
plt.show()

#==============================================================================
#   Grafica F(ν)
#==============================================================================
#   Soporte de unidades físicas
from astropy import units as u 

#   Soporte de magnitudes físicas. Solo se utiliza la velocidad de la luz.
from astropy.constants import c 

#   Unidades a las longitudes de onda
x_sed=x_sed*u.AA # En Ångström
#   Unidades de densidad espectral de flujo (F_(lambda))
y_sed=y_sed*((u.erg)* (u.cm)**(-2)* (u.s)**(-1)* (u.AA)**(-1))
#   Cambio de F_(lambda) a F(ípsilon) multiplicando por (lambda**2)/c 
y_sed=y_sed*(x_sed**2)*(c**-1)
#   Cambio de unidades 
x_sed=x_sed.to(u.micron)
y_sed=y_sed.to(u.mJy) 

figura=plt.figure()
plt.loglog(x_sed,y_sed, color="blue", linewidth=1, linestyle="-",
        label="SMM J2")
plt.grid(False)
plt.title('SMM J2135-0102', fontsize = 16, color='blue')
plt.xlabel('$ \lambda\;[\mu \mathrm{m}]$', fontsize = 16, color='blue')
plt.ylabel('$ \mathrm{F_{ν}}}$',fontsize = 16, color='blue')
figura.savefig(outroute + 'grafica_SMM_F_ípsilon', dpi=150, transparent=True)
plt.show()