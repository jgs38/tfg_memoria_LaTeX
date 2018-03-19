# -*- coding: utf-8 -*-
"""
Contiene las funciones que se utilizan para calcular los factores de Bayes.

  Python 3.6.0    numpy 1.11.3

@author: Javier Gutiérrez Solórzano
"""

import math
import numpy as np

def bayes_objeto(s1,s2,psi,z1,sz1,z2,sz2,z_max):
    """
    Permite el cálculo del factor de Bayes posicional, fotométrico y 
    conjunto para cada emparejamiento encontrado. Esta función contiene dos 
    funciones; una para calcular el factor de Bayes posicional y otra para el 
    factor de Bayes fotométrico. El factor de Bayes conjunto se calcula como el
    producto de ambos.
    
    Parámetros:
    Se explican en el comentario de cada sub-función.
    
    Return:
    bp ----- Factor de Bayes posicional.
    bz ----- Factor de Bayes fotométrico.
    bp*bz -- Factor de Bayes conjunto.
    """
    def bayesPosicional(s1,s2,psi):
        """
        s1 --- Desviación estándar de las medidas posicionales de uno de los
               catálogos.
        s2 --- Desviación estándar de las medidas posicionales del segundo
               catálogo.
        psi -- Distancia angular a la que se encuentran los objetos que 
               forman el emparejamiento.
        """
        x_1=(s1**2+s2**2)
        bayes_posicional=(2/x_1)*math.exp(-( (psi**2)/(2*x_1) ))
        
        return bayes_posicional
        
    def bayesZ(sz1,sz2,z_max,z1,z2):
        """
        sz1 ---- Desviación estándar asociada a z1.
        sz2 ---- Desviación estándar asociada a z2.
        z_max -- 
        z1 ----- Valor del redshift de los objetos de uno de los catálogos. 
        z2 ----- Valor del redshift del segundo catálogo.
        """
        if sz1>0 and sz2>0 and z_max>0 and z1>0 and z2>0 and z_max>z1 \
            and z_max>z2:
                
            sz12=sz1**2
            sz22=sz2**2
            
            f_error_1_numerador=math.erf((sz12*z2+sz22*z1)/(sz1*sz2*\
                        math.sqrt(2*(sz12+sz22))))
            f_error_2_numerador=math.erf((sz12*(z2-z_max)+sz22*(z1-z_max))/\
                        (sz1*sz2*math.sqrt(2*sz12+sz22)))
            numerador=math.sqrt(2/math.pi)*z_max*math.exp(-((z1-z2)**2)/\
                        (sz12+sz22))*(f_error_1_numerador-f_error_2_numerador)
            
            dif_error1_denominador=math.erf(z1/sz1*math.sqrt(2))-\
                                        math.erf((z1-z_max)/(sz1*math.sqrt(2)))
            dif_error2_denominador=math.erf(z2/sz2*math.sqrt(2))-\
                                        math.erf((z2-z_max)/(sz2*math.sqrt(2)))
            denominador=math.sqrt(sz12+sz22)*\
                                 dif_error1_denominador*dif_error2_denominador
                                 
            if np.any(denominador==0):
                denominador=1e-80
            
            bayes_z=numerador/denominador
        
            return bayes_z
        else:
            # Si no se cumplen las condiciones, Bz no nos aporta información
            return 1
    
    bp=bayesPosicional(s1,s2,psi)
    bz=bayesZ(sz1,sz2,z_max,z1,z2)
    
    return bp, bz, bp*bz

    
def bayes_conjunto(s1,s2,psi,z1,sz1,z2,sz2,z_max):
    """
    Aplica la función bayes_objeto a un conjunto de emparejados. 
    
    Parámetros:
    Se explican en el comentario de las funciones "bayesPosicional" y "bayesZ".
    
    Return:
    bp ----- Array que contiene los factores de Bayes posicionales.
    bz ----- Array que contiene los factores de Bayes fotométricos.
    bp*bz -- Array que contiene los factores de Bayes conjuntos.
    
    """
    array_bayes_posicional=[]
    array_bayes_z=[]
    array_bayes=[]
    for i in range(0,len(psi),1):
        [bayes_posicional, bayes_z, conjunto] = bayes_objeto(s1,s2,psi[i]
        ,z1[i],sz1[i],z2[i],sz2[i],z_max)
        array_bayes_posicional.append(bayes_posicional)
        array_bayes_z.append(bayes_z)
        array_bayes.append(conjunto)
    
    return array_bayes_posicional, array_bayes_z, array_bayes