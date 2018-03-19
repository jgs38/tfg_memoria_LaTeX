# -*- coding: utf-8 -*-

"""
Se realiza un ajuste lineal para obtener los errores asociados a las medidas 
espectroscópicas del redshift para los proyectos HALTAS y GAMA.
     
    Python 3.6.0   scipy 0.18  matplotlib 2.0.0 numpy 1.11.3
     
@author: Javier Gutiérrez Solórzano
"""
#==============================================================================
#    Nombres de los ficheros de entrada/salida y su localización.
#==============================================================================
#--------------------------
#   Directorio en el que se encuentran los módulos.
#--------------------------
import sys
sys.path.append("../pakages/")
#--------------------------
#   Directorios de los ficheros.
#--------------------------
outroute = './outroute/'

#   Redshift espectrocópico 
z=[0.312005,0.385261,0.133342,0.125497,0.126204,0.130568,0.123867,0.103692,
   0.223991,0.385721,0.374708,0.081022,0.394029,0.134499,0.476045,0.291690,
   0.305359,0.038159,0.326157,0.293476]
#   Errores asociados a las medidas espectroscópicas (desviación estándar)
ez=[0.000102,0.000103,0.000092,0.000079,0.000100,0.000133,0.000136,0.000048,
    0.000089,0.000117,0.000120,0.000130,0.000175,0.000119,0.000180,0.000098,
    0.000119,0.000065,0.000351,0.000351]

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def func(x, A):
    """
        Función que relaciona el redshift espectrocópico con su error. 
    """
    return A*(1+x)

#   pcov - Matriz covarianza
#   popt - Vector con valores de ajuste 
popt, pcov = curve_fit(func, z, ez)

def recta_sup(x,opt_par,cov_matrix,nsigmas):  #   función que calcula el 
            # extremo superior de la banda, a NSIGMAS sigmas
    ymean = opt_par[0]*(x+1)
    sa    = cov_matrix[0]
    y     = ymean+nsigmas*np.sqrt(sa)
    return y
 
def recta_inf(x,opt_par,cov_matrix,nsigmas):  #    función que calcula el 
             # extremo inferior de la banda, a NSIGMAS sigmas
    ymean = opt_par[0]*(x+1)
    sa    = cov_matrix[0]
    y     = ymean-nsigmas*np.sqrt(sa)
    return y

#-------------------------- 
#   Grafica z_articulo vs z_programa
#-------------------------- 
figura = plt.figure(num = None, figsize = (9, 6), dpi = 80, facecolor = 'w', 
                    edgecolor = 'k')
#   Cambiamos los limites de los ejes.
plt.axis([0, 0.6, 0, 0.0004])
#   Genera 100 valores igualmente especiados entre 1 y 3.5
x = np.linspace(0, 0.6, 10, endpoint=True)
#   Recta del ajuste
plt.plot(x,popt[0]*(x+1),'r',linewidth=0.5)
#   Sombreado amarillo para la banda 1 sigma
plt.fill_between(x,recta_inf(x,popt,pcov,2),recta_sup(x,popt,pcov,2),
                 facecolor='yellow', alpha=0.8)
#   Sombreado amarillo para la banda 2 sigma
plt.fill_between(x,recta_inf(x,popt,pcov,1),recta_sup(x,popt,pcov,1),
                 facecolor='orange', alpha=0.8)
#   Dibuja los puntos 
plt.errorbar(z, ez, fmt='+', label="data", linewidth=0.15, ecolor='black')
#   Sin rejilla de fondo.
plt.grid(False)
#   Etiqueta eje x.
plt.xlabel('${z}^{\mathrm{spec}}$', fontsize = 16, color='blue')
#   Etiqueta eje y.
plt.ylabel('${\sigma}^{z}$', fontsize = 16, color='blue')
#   Muestra la figura.
plt.show()
#   Guarda una imagen de la gráfica.
figura.savefig(outroute + "ajuste_errores_z", dpi=150, transparent=True)

#-------------------------- 
#   Muestra datos del ajuste por pantalla
#--------------------------   
print('Errores z_spec vs z_spec:\n',
'--- Best fit    y = ({0}'.format(popt[0]), chr(177) ,
'{0})'.format(np.sqrt(pcov[0,0])), chr(215),'(1+x)')