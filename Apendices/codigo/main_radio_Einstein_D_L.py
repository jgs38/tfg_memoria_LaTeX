# -*- coding: utf-8 -*-

"""
Gráfica de radio de Einstein frente a la distancia observador-lente.

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
outroute = './outroute/'

import numpy as np

import matplotlib.pyplot as plt

from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

print(cosmo)

from astropy import units as u

#   Genera 100 valores igualmente especiados entre 1 y 3.5
x_z= np.linspace(1e-5, 0.6, 100, endpoint=True)

x_pc=cosmo.luminosity_distance(np.array(x_z))
x_pc=x_pc.to(u.pc)
x_pc=x_pc.value

r_e_9=0.1*np.sqrt((1e9/np.array(x_pc)))#*u.arcsec
r_e_12=0.1*np.sqrt((1e12/np.array(x_pc)))#*u.arcsec

x_z_eje= np.linspace(0.1, 0.6, 6, endpoint=True)
x_pc_eje=np.round(cosmo.luminosity_distance(np.array(x_z_eje)).value, decimals=0, out=None)

#-------------------------- 
#   Figura radio de Einstein vs distancia observador-lente.
#-------------------------- 
figura = plt.figure(num = None, figsize = (9, 4), dpi = 80, facecolor = 'w', 
                    edgecolor = 'k')
#   Cambiamos los límites de los ejes.
plt.axis([0, 0.6, 4e-2, 50])

#   Graficas para cada masa
plt.plot(x_z,r_e_12, linewidth=0.7, color="red", linestyle='--')
plt.plot(x_z,r_e_9, linewidth=0.7, color="red")

#   Líneas horizontales
plt.axhline(y=7.63, linewidth=0.5, color='k', linestyle='-.')
plt.axhline(y=0.297, linewidth=0.5, color='k', linestyle=':')

#   Escala de cada eje
plt.yscale('log', linthreshx=0.1)
#plt.xscale('log', linthreshx=0.1)

#   Sin rejilla de fondo.
plt.grid(False)

ax = plt.gca()
#   Etiqueta eje y.
ax.set_ylabel(r"${\theta}_{E}\;[\mathrm{arcsec}]$", fontsize=16, color='blue')

#   Etiqueta eje x.
plt.xlabel(r"${z}_{L}$", fontsize = 16, color='blue')

ax2 = ax.twiny()
ax2.set_xlabel("${D}_{L}\;[\mathrm{Mpc}]$", fontsize=16, color='blue')
ax2.set_xlim(0, 0.6)
ax2.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax2.set_xticklabels([460,980,1553,2172,2833,3530])

plt.show()
#   Guarda gráfica.
figura.savefig(outroute + "radio_einstein", dpi=150, transparent=True)