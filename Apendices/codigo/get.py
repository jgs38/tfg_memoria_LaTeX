# -*- coding: utf-8 -*-

"""
Módulo del programa que contiene la mayoria de las funciones que utilizan en 
programas (salvo las relativas a la obtención del redshift y de la obtención
de los factores de Bayes). Se ha decidido hacer así porque muchas de estas 
funciones se utilizan en varios modulos; además, se facilita la comprensión de 
los otros módulos y se evita que el programador se pierda en los detalles.

  Python 3.6.0    numpy 1.11.3    astropy 1.3

@author: Javier Gutiérrez Solórzano
"""

from astropy.io import fits
import numpy as np 

def sed_func(fichero_sed):
    """
    Función para leer el fichero que contiene los puntos SED teórica de la 
    galaxia 'SMMJ2135-0102' que consiste en dos columnas: la primera con los 
    valores de longitudes de onda y la segunda valores de la densidad espectral
    de flujo F_(lambda).
    
    Parámetros:
    fichero_sed -- Nombre del fichero.
    
    Return:
    x_teorica ---- Array, primera columna del fichero.
    y_teorica ---- Array, segunda columna del fichero.
    
    """
    x_teorica=[]
    y_teorica=[]
    
    file=open(fichero_sed,'r')
    contenido=file.readlines()
    #   Obtenemos el número de líneas que tiene el fichero (¡cuidado con 
    # ficheros muy grandes, pueden dar problemas con el uso de la memoria!)
    lineas = len(open(fichero_sed).readlines())
        
    for i in range (0,lineas,1): 
        dato1=contenido[i].find(' ')# (separador espacio ' ')
        dato2=contenido[i].find('\n',dato1+1) # (salto de línea '\n')
        #   Cambia el tipo de cada elemento.
        valor1=float(contenido[i][0:dato1]) 
        valor2=float(contenido[i][dato1+1:dato2])
        
        x_teorica.append(valor1) 
        y_teorica.append(valor2)
    
    file.close()
    
    return x_teorica, y_teorica


def columna_fits(catalogo,n_columna):
    """
    La lectura de las columnas individuales de un fichero .fits 
    
    Parámetros:
    catalogo ---- Nombre del fichero.
    n_columnas -- Posición de la columna 
    
    Return:
    Array con el contenido de la columna 'n_columna' del fichero 'catalogo'.        
    """
        
    return fits.open(catalogo, memmap=True)[1].data.field(n_columna)

          
def datos_experimentales(inroute_infile_datos_experimentales):
    """
    Lee las columnas que contienen las medidas del flujo del
    catálogo 'HATLAS_DR1_CATALOGUE_V1.2.FITS' y devuelve los datos en un 
    formato adecuado para usarse por las funciones definidas en el módulo 
    'rojo.py'.
    
    Cada objeto astronómico del catálogo contiene tres medidas del flujo, en 
    las longitudes de onda 250, 350 y 500 micrómetros (columnas 4, 5 y 6) y con
    un error asociado a cada una de ellas (columnas 7, 8 y 9).
    
    Parámetros:
    inroute_infile_datos_experimentales -- Nombre del fichero.
    
    Return:
    y_experimental -- Array cuyos elementos son un array de 3 elementos con
                      las medidas del flujo de la observación.
    y_error --------- Array cuyos elementos son un array de 3 elementos con el
                      "error" del flujo.        
    """

    #   Carga archivo fits.
    fitsFile = fits.open(inroute_infile_datos_experimentales, memmap=True) 
        
    #   Seleccionamos toda la tabla.
    tbdata = fitsFile[1].data 
    
    #   Transformamos una columna del archivo de entrada, infile, en un 
    # array. Éste nos proporcina las condiciones de ejecucion (a partir del
    # tamaño del array) al bucle for.
    longitud_infile_array=len(np.array(tbdata.field(5)))
    
    #-------------------------- 
    # ARRAYS DE DATOS EXPERIMENTALES.
    #-------------------------- 
    
    y_experimental=[]
    for i in range(0,longitud_infile_array,1):
        y_experimental_objeto=[]
        y_experimental_objeto.append(tbdata[i].field(4))
        y_experimental_objeto.append(tbdata[i].field(5))
        y_experimental_objeto.append(tbdata[i].field(6))
    
        y_experimental.append(y_experimental_objeto)
       
    y_error=[]
    for i in range(0,longitud_infile_array,1):       
        y_error_objeto=[]
        y_error_objeto.append(tbdata[i].field(7))
        y_error_objeto.append(tbdata[i].field(8))
        y_error_objeto.append(tbdata[i].field(9))
    
        y_error.append(y_error_objeto)
        
        
    return y_experimental, y_error  


def ordena_Array(idxa,array_datos):
    """
    Reordena un array a partir de los valores numéricos de otro que contiene 
    números enteros. El elemento que ocupa la posición [i] en el array toma la 
    posición array_datos[idxa[i]] en el array de salida. La longitud de
    array_datos debe ser mayor o igual que el mayor de los idxa[i].
    
    Debe cumplirse siempre que: len(array_datos)-1 <= max(idxa)
    
    Ejemplo: 
    idxa = [0, 6, 3, 4, 5, 5, 4, 5]           
    array_datos =  ["a","b","c","d","e","f","g"] 
    ordena_Array(idxa,array_datos) = ['a', 'g', 'd', 'e', 'f', 'f', 'e', 'f']
    
    Parámetros:
    idxa --------- Array que contiene números enteros.
    array_datos -- Array con los elementos que queremos ordenar.
    
    Return:
    nuevo_orden -- Array que contiene los elementos de array_datos ordenados
                   según el orden del array idxa.        
    """
    nuevo_orden=[]
    for i in range(0,len(idxa),1):
        nuevo_orden.append(array_datos[idxa[i]])
        
    return nuevo_orden


def elimina_elementos_repetidos(idx):
    """
    Ordena los valores numéricos de un array de menor a mayor y deja un solo
    valor de aquellos que se repiten.
    
    Ejemplo:
        
    idx=[0, 9, 3, 4, 5, 5, 4, 9] 
    elimina_elementos_repetidos(idx)=[0, 3, 4, 5, 9] 
    
    Parámetros:
    idx ---------- Array con números enteros.
    
    Return:
    b  ----------- Array que contiene los mismos números que idx ordenado de 
                   menor a mayor y sin elementos repetidos.
    """
    idx_sort=sorted(idx)
    b=[]
    if len(idx)==0:
        return b
    else:
        b.append(idx[0])
        for i in range(0,len(idx),1):
            if idx_sort[i] != b[len(b)-1]:
                b.append(idx_sort[i])
        return b


def id_pos(lista):
    """
    Permite obtener un array de número naturales de la misma longitud que el 
    array de entrada.
    
    Parámetros:
    lista ------- Array cualquiera.
    
    Return:
    contador_i -- Array de números naturales de longitud len(lista).
        
    """
    contador_i=[]
    for i in range(0,len(lista),1):
        contador_i.append(i)
    return contador_i

    
def posicion_fits(idx):
    """
    Sirve para encontrar la fila que ocupa un objeto en uno de los catálogos 
    astronómicos. Esta función es necesaria porque el primer elemento de un 
    array en Python ocupa la posición 0 mientras que la primera fila del 
    catálogo, tiene asignado el valor 1.
    
    Parámetros:
    idx ------ Array de posiciones, siendo el primer elemento 0.    
    
    Return:
    id_fits -- Array de posiciones, siendo el primer elemento 1.       
    
    """
    id_fits=[]
    for i in range(0,len(idx),1):
        id_fits.append(idx[i]+1)
        
    return id_fits


def frange(start, end, step):
    """
    El comportamiento de esta función es el mismo que el de la función "range"
    pero, a diferencia de ésta, permite que cualquiera de los parámetros de
    entrada no sea un entero.
    
    Parámetros:        
    start -- Limite inferior de la secuencia (incluído en la misma).
    end ---- Límite superior de la secuencia (no está incluído).
    step --- Diferencia entre cada número de la secuencia.
        
    """
    tmp = start
    while(tmp < end):
        yield tmp
        tmp += step
        

from astropy.coordinates import SkyCoord
from astropy import units as u 

def area_region(x_inf,x_sup,y_inf,y_sup,ra_catalogo_GAMA,dec_catalogo_GAMA
                                ,ra_catalogo_HATLAS,dec_catalogo_HATLAS,paso):
    """
    Permite obtener el ángulo sólido de las regiones en las que se encuentran
    los objetos del catálogo HATLAS, GAMA y la zona de intersección. Para ello
    se crea una cuadrícula cuyos elementos están separados a una distancia del 
    orden de la separación angular de los objetos de los catálogos; después se 
    realiza un conteo de aquellos puntos de la cuadrícula que tienen objetos de
    un determinado catálogo (o ambos, para determinar la zona de intersección)
    a una distancia la mitad de la diagonal de cuadrado que conforma la celda 
    (paso*(np.sqrt(2))**(-1)) que es la distancia mínima posible para que todo 
    el espacio quede cubierto. Para obtener el ángulo sólido de la región, 
    se multiplica el número de objetos por el área de cada celda.    
    
    Parámetros:
    x_inf ---------------- Límite horizontal inferior de la cuadrícula.
    x_sup ---------------- Límite horizontal superior de la cuadrícula.
    y_inf ---------------- Límite vertical inferior de la cuadrícula.
    y_sup ---------------- Límite vertical superior de la cuadrícula.
    ra_catalogo_GAMA ----- Ascensión recta de los objetos de GAMA.
    dec_catalogo_GAMA ---- Declinación de los objetos de GAMA.
    ra_catalogo_HATLAS --- Ascensión recta de los objetos de HATLAS.
    dec_catalogo_HATLAS -- Declinación de los objetos de HATLAS.
    paso ----------------- Longitud del lado de las celdas que forman la 
                           cuadrícula.    
    Return:
    area_hatlas ---------- Ángulo sólido que cubren los objetos de HATLAS.
    area_gama ------------ Ángulo sólido que cubren los objetos de GAMA.
    area_interseccion ---- Ángulo sólido que cubren los objetos que forman la 
                           intesección de ambos catálogos.        
    """
    
    def regilla(x_inf,x_sup,y_inf,y_sup,paso):
        # Creo una regilla.
        px=[]
        py=[]
        for y in frange(y_inf,y_sup,paso):
            for x in frange(x_inf,x_sup,paso):
                py.append(y)
                px.append(x)
        
        return px, py
    [px,py]=regilla(x_inf,x_sup,y_inf,y_sup,paso)
    # Área de cada de cada elemento del grid.
    area=(paso*u.degree)**2
    # De esta forma nos aseguramos de cubrir todo el espacio del grid.
    sep=paso*((np.sqrt(2))**(-1))*u.degree
    # Asignamos unidades a la regilla.
    grid = SkyCoord(ra=px*u.degree, dec=py*u.degree)
             
    h = SkyCoord(ra=ra_catalogo_HATLAS*u.degree
                                            , dec=dec_catalogo_HATLAS*u.degree)
    [id_grid_h, _, _, _] = h.search_around_sky(grid, sep)
    area_hatlas=len(elimina_elementos_repetidos(id_grid_h))*area
                   
    g = SkyCoord(ra=ra_catalogo_GAMA*u.degree, dec=dec_catalogo_GAMA*u.degree)
    [id_grid_g, _, _, _] = g.search_around_sky(grid, sep)
    area_gama=len(elimina_elementos_repetidos(id_grid_g))*area
    
    area_interseccion=len(list(set(id_grid_g).intersection(id_grid_h)))*area
    
    return area_hatlas,area_gama,area_interseccion 


def numero_halos(idxh,s250,s350,s500):
    """
    Realiza un conteo de todos aquellos objetos del catálogo HATLAS que cumplen
    los criterios de Gonzalez-Nuevo et al. y de Negrello et al. para ser 
    lente gravitatoria partiendo de las medidas del instrumento SPIRE.
    
    Parámetros:
    idxh -- Array con las posiciones de los objetos de HATLAS.
    s250 -- Medidas del flujo en 250 microm.
    s350 -- Medidas del flujo en 350 microm.
    s500 -- Medidas del flujo en 500 microm.
    
    Return:
    len(id_halo_negrello) -- Número de objetos tipo Negrello et al.
    len(id_halo_gonzalez) -- Número de objetos tipo Gonzalez-Nuevo et al.
    len(id_halo_ambos) ----- Número de objetos que cumplen ambas condiciones.
    
    """
    id_halo_negrello=[]
    id_halo_gonzalez=[]
    id_halo_ambos=[]
    for i in range(0,len(idxh),1): 
        negrello=False
        gonzalez=False
        # Para evitar divisiones por 0 (menor que la precisión de la 
        # medida del instrumento SPIRE)
        if np.any(s250[i] == 0): 
            s250[i]=1e-6
        if np.any(s250[i] == 0): 
            s350[i]=1e-6
        s350_s250 = s350[i]/s250[i]
        s500_s350 = s500[i]/s350[i]
       
        if np.any(s500[i] >= 0.1):
            negrello=True
        if np.any(s350[i] >= 0.085) and np.any(s250[i] >= 0.035) \
                 and np.any(s350_s250 > 0.6) and np.any(s500_s350 > 0.4):
            gonzalez=True
        
        if negrello==True:
            id_halo_negrello.append(idxh[i])
        if gonzalez==True:
            id_halo_gonzalez.append(idxh[i])
        if negrello==True and gonzalez==True:
            id_halo_ambos.append(idxh[i])
           
    return len(id_halo_negrello),len(id_halo_gonzalez),len(id_halo_ambos)


def obj_by_flag_hatlas(idx,flag,tipo):
    """
    Realiza un conteo de los objetos del catálogo HATLAS en función de la 
    calidad de la medida del redshift de que se dispone y nos devuelve un array
    o un número dependiendo de la variable "tipo" que hemos escogido. (la 
    medida de la calidad, para cada objeto, se indica por el valor de flag[i]).
    
    Parámetros:
    idx --- Array con las posiciones de los objetos de HATLAS.
    flag -- Array con la etiqueta sobre la calidad del redshift.
    tipo -- Toma el valor "array" si queremos que nos devuelva los 
            identificadores de los objetos o bien el valor "numero" si estamos
            interesados en el número de objetos con una determinada medida del
            factor de calidad.
    
    Return:
    return 1 -- Devuelve las longitude de los arrays "return 2".
    return 2 -- Contiene los identificadores de los objetos, en 4 arrays 
                diferentes, dependiendo de si su redshift ha sido obtenido
                a partir de una medida espectroscópica (id_spec), si se ha
                obtenido a partir de ANNZ (id_annz), ha sido obtenida mediante
                el ajuste propuesto en este trabajo (id_phot)o bien no está 
                disponible (id_disabled). El último elemnto que devuelve es
                el array con los identificadores de todos los objetos 
                considerados.
    return 3 -- En caso de haber introducido un valor incorrecto en el 
                parámetro "tipo".
    """
    id_spec=[]
    id_annz=[]
    id_phot=[] 
    id_disabled=[]
   
    for i in range(0,len(flag),1):
        if flag[i]==11 or flag[i]==12 or flag[i]==13:
            id_spec.append(idx[i])
        elif flag[i]==2:
            id_annz.append(idx[i])
        elif flag[i]==3:
            id_phot.append(idx[i])
        else:
            id_disabled.append(idx[i])
  
    if tipo=="numero":
        return len(elimina_elementos_repetidos(id_spec)), \
               len(elimina_elementos_repetidos(id_annz)), \
               len(elimina_elementos_repetidos(id_phot)), \
               len(elimina_elementos_repetidos(id_disabled)), \
               len(elimina_elementos_repetidos(idx))
    elif tipo=="array":
        return elimina_elementos_repetidos(id_spec), \
               elimina_elementos_repetidos(id_annz), \
               elimina_elementos_repetidos(id_phot), \
               elimina_elementos_repetidos(id_disabled), \
               elimina_elementos_repetidos(idx)
    else:
        return "la variable tipo debe ser 'numero' o 'array' "     
       

def obj_by_flag_gama(idx,flag,tipo):
    """
    Realiza un conteo de los objetos del catálogo GAMA en función de la 
    calidad de la medida del redshift de que se dispone y nos devuelve un array
    o un número dependiendo de la variable "tipo" que hemos escogido. 
    
    Parámetros:
    idx --- Array con las posiciones de los objetos de GAMA.
    flag -- Array con la etiqueta sobre la calidad del redshift.
    tipo -- Toma el valor "array" si queremos que nos devuelva los 
            identificadores de los objetos o bien el número de observaciones
            con una determinada medida del factor de calidad.
    
    Return:
    return 1 -- Devuelve las longitude de los arrays "return 2".
    return 2 -- Contiene los identificadores de los objetos, en 4 arrays 
                diferentes, dependiendo de si su redshift ha sido obtenido
                a partir de una medida espectroscópica (id_spec11), si se ha
                obtenido a partir de ANNZ (id_spec12), ha sido obtenida mediante
                el ajuste propuesto en este trabajo (id_spec13)o bien no está 
                disponible (id_disabled). El último elemento que devuelve es
                el array con los identificadores de todos los objetos 
                considerados.
    return 3 -- En caso de haber introducido un valor incorrecto en el 
                parámetro "tipo", nos indica los valores posibles que puede 
                tomar.
    """
    id_spec11=[]
    id_spec12=[]
    id_spec13=[]
    id_disabled=[]
   
    for i in range(0,len(flag),1):
        if flag[i]==11:
            id_spec11.append(idx[i])
        elif flag[i]==12:
            id_spec12.append(idx[i])
        elif flag[i]==13:
            id_spec13.append(idx[i])
        else:
            id_disabled.append(idx[i])
  
    if tipo=="numero":
        return len(elimina_elementos_repetidos(id_spec11)), \
               len(elimina_elementos_repetidos(id_spec12)), \
               len(elimina_elementos_repetidos(id_spec13)), \
               len(elimina_elementos_repetidos(id_disabled)), \
               len(elimina_elementos_repetidos(idx))
    elif tipo=="array":
        return elimina_elementos_repetidos(id_spec11), \
               elimina_elementos_repetidos(id_spec12), \
               elimina_elementos_repetidos(id_spec13), \
               elimina_elementos_repetidos(id_disabled), \
               elimina_elementos_repetidos(idx)
    else:
        return "tipo debe ser 'numero' o 'array' "  
   
import copy

def representacion_candidatos(bp,bz_original,bc,nombre_hatlas,
    nombres_confirmados,s250_original,s350_original,s500,flag,zh,zg,seleccion):
    """
    Esta función clasifica los emparejamientos y agrupa sus factores de Bayes
    en varios arrays para realizar una representación adecuada de los mismos.
    
    Parámetros:
    bp ------------------- Factor de Bayes posicional.
    bz_original ---------- Factor de Bayes fotométrico.
    bc ------------------- Factor de Bayes conjunto.
    nombre_hatlas -------- Nombres de los objetos HATLAS que han participado en
                           el matching.
    nombres_confirmados -- Nombres que de los objetos de los que queremos 
                           destacar el pareado en el que participan.
    s250_original -------- Medidas del fujo en 250 microm.
    s350_original -------- Medidas del fujo en 350 microm.
    s500 ----------------- Medidas del fujo en 500 microm.
    flag ----------------- Array con la etiqueta sobre la calidad del redshift
                           de los objetos de HATLAS.
    zh ------------------- Redshift de la observación de HATLAS.
    zg ------------------- Redshift de la observación de GAMA.
    seleccion ------------ Nos permite seleccionar solo aquellos objetos del
                           catálogo HATLAS cuyo redshift ha sido obtenido a 
                           partir del ajuste propuesto en este trabajo 
                           (flag[i]==3).
                           
    Return:
    return --  Devuelve varios arrays que contienen el logaritmo en base 10 de 
               los factores de Bayes posicional y fotométrico. La casificación
               está basada en el valor de los factores de Bayes y si se cumplen
               o no los criterios de Negrello y Gonzalez-Nuevo.
    """
    #   Se van a modificar estos arrays. Para no sobreescribir sobre los
    # originales se hace una copia. Sin estas lineas de código, el fits se 
    # genera a partir de bayes_z modificado (bz).
    bz = copy.copy(bz_original)
    s250 = copy.copy(s250_original)
    s350 = copy.copy(s350_original)
    
    # Creo los arrays vacios.
    [x_0, y_0, x_1, y_1, x_2, y_2] = [[],[],[],[],[],[]]
    [x_3, y_3, x_4, y_4, x_5, y_5] = [[],[],[],[],[],[]]
    
    # Si seleccion=ajuste representamos las emparejados en los que participa 
    # una observación de HATLAS con flag=3
    if seleccion=="ajuste":
        (a,b,c,d,e)=(3,3,3,3,3)    
    else:
        (a,b,c,d,e)=(11,12,13,2,3)
    
    for i in range(0,len(nombre_hatlas),1):
        if (flag[i]==a or flag[i]==b or flag[i]==c or flag[i]==d or \
            flag[i]==e) and zh[i]>=1 and zh[i]>zg[i]:

            # Para evitar log(0) posteriormente.
            if np.any(bp[i]<1e-100):
                bp[i]=1e-100
            if np.any(bz[i]<1e-100):
                bz[i]=1e-100
            #--------------------------
            # Selecciona los objetos de cuyos nombres están en el fichero.
            #--------------------------
            for o in range(0,len(nombres_confirmados),1): 
                nombre_candidato=False
                if nombre_hatlas[i]==nombres_confirmados[o]:
                    nombre_candidato=True
                    break
            if nombre_candidato==True: 
                x_0.append(bp[i])
                y_0.append(bz[i])
            else:
                #--------------------------
                # Clasificación de los emparejamientos. 
                #--------------------------
                if np.any(s500[i] >= 0.1):  # Negrello et al.
                    x_1.append(bp[i]) 
                    y_1.append(bz[i]) 
                else:
                    # Para evitar divisiones por 0 (menor que la precisión de
                    # la medida del instrumento SPIRE)
                    if np.any(s250[i] < 1e-6): 
                        s250[i]=1e-6
                    if np.any(s350[i] < 1e-6): 
                        s350[i]=1e-6                    
                    s350_s250 = s350[i]/s250[i]
                    s500_s350 = s500[i]/s350[i]
                    
                    if np.any(s350[i] >= 0.085) and np.any(s250[i] >= 0.035) \
                             and np.any(s350_s250 > 0.6) \
                             and np.any(s500_s350 > 0.4):
                        x_2.append(bp[i])   # Condición (Gonzalez-Nuevo et al)
                        y_2.append(bz[i]) 
                            
                    elif np.any(bc[i]>100): # Bayes conjunto mayor que 100.
                        x_5.append(bp[i])
                        y_5.append(bz[i])
                    elif np.any(bc[i]>1):   # Bayes conjunto mayor que 1.
                        x_4.append(bp[i])
                        y_4.append(bz[i])
                    else:
                        x_3.append(bp[i])   # Resto.
                        y_3.append(bz[i])
            
    return np.log10(x_0), np.log10(y_0), np.log10(x_1), np.log10(y_1), \
           np.log10(x_2), np.log10(y_2), np.log10(x_3), np.log10(y_3), \
           np.log10(x_4), np.log10(y_4), np.log10(x_5), np.log10(y_5)

def selecion_candidatos(s250_original,s350_original,s500,flag_h,zh,zg,
                        bayes_posicional,bayes_conjunto): 
    """
    Permite seleccionar aquellos emparejamientos que cumplen nuestro criterio 
    para ser considerados candidatos a lente gravitatoria. Los emparejamientos
    que cumplen nuestro criterio se agrupan en dos grupos dependiendo de la 
    calidad de la medida del redhsift. Por último se cuenta el número de 
    amparejamientos que además contienen observaciones en HATLAS que cumplen 
    los criterios de Negrello y González-Nuevo para ser candidatos a SLGs.
    
    Parámetros:
    s250_original ----- Medidas del flujo en 250 microm.
    s350_original ----- Medidas del flujo en 350 microm.
    s500 -------------- Medidas del flujo en 500 microm.
    flag_h ------------ Etiqueta de calidad sobre la medida del redshift de la
                        observación perteneciente a HATLAS.
    zh ---------------- Redshift perteneceinte a HATLAS
    zg ---------------- Redshift perteneceinte a GAMA
    bayes_posicional -- Factor de Bayes posicional
    bayes_conjunto ---- Factor de Bayes conjunto
    
    Return:
    id_halo_h_espec ------------ Identificadores de los emparejamientos que 
                                 cumplen nuestros criterios para ser candidatos
                                 y estan formados por observaciones que tienen
                                 medidas espectroscópicas del redhsift.
    id_halo_h_fot -------------- Identificadores de los emparejamientos que 
                                 cumplen nuestros criterios para ser candidatos
                                 pero la observación de HATLAS tiene un 
                                 redshift fotométrico.
    id_halo_negrello_gonzalez -- Identificadores de los emparejamientos que 
                                 cumplen nuestros criterios y la observación de
                                 HATLAS cumple además los criterios de Negrello
                                 y González-Nuevo.
    
    """
    s250 = copy.copy(s250_original)
    s350 = copy.copy(s350_original)
   
    id_halo_h_espec=[]
    id_halo_h_fot=[]
    id_halo_negrello_gonzalez=[]
    
    for i in range(0,len(s250_original),1): 
        if zh[i]>=1 and zh[i]>zg[i]:
            negrello=False
            gonzalez=False
            bayes=False
            # Para evitar divisiones por 0 (menor que la precisión de la 
            # medida del instrumento SPIRE)
            if np.any(s250[i] < 1e-6): 
                s250[i]=1e-6
            if np.any(s350[i] < 1e-6): 
                s350[i]=1e-6  
            s350_s250 = s350[i]/s250[i]
            s500_s350 = s500[i]/s350[i]
           
            if np.any(s500[i] >= 0.1):
                negrello=True
            if np.any(s350[i] >= 0.085) and np.any(s250[i] >= 0.035) \
                     and np.any(s350_s250 > 0.6) and np.any(s500_s350 > 0.4):
                gonzalez=True
            if np.any(bayes_conjunto[i]<1) and np.any(bayes_posicional[i]>1):
                bayes=True
           
            if bayes==True:
                if flag_h[i]==3:
                    id_halo_h_fot.append(i)
                else:
                    id_halo_h_espec.append(i)
                
            if bayes==True and (negrello==True or gonzalez==True):
                id_halo_negrello_gonzalez.append(i)

            
    return id_halo_h_espec,id_halo_h_fot,id_halo_negrello_gonzalez

import math
import random 

def uniform_spherical_distribution(RA_a,RA_b,DEC_a,DEC_b,n_vectors):
    """
    Genera una distribución esférica uniforme de puntos aleatorios comprendida 
    entre los meridianos "RA_a" y "RA_b" y los paralelos "DEC_a" y "DEC_b". La 
    variable "n_vectors" indica a la función en número de vectores que se 
    debe generar.
    
    Parámetros:
        
    RA_a y RA_b -------- Meridianos límite.
    DEC_a y DEC_b ------ Paralelos límite.
    n_vectors -- Número de vectores que queremos generar.
    
    Return:
    RA_random --- Array con la ascensión recta de los vectores generados.
    DEC_random -- Array con la declinación de los vectores generados.
    
    
    """
    def random_3D_unit_vector(RA_a,RA_b,DEC_a,DEC_b):
        """
        Genera un vector unitario en un área definida entre "RA_a" y "RA_b" y 
        "DEC_a" y "DEC_b" de forma aleatoria
        """
        RA = random.uniform(RA_a,RA_b) #0,360
        
        theta_a = 90-DEC_a     #-90,90
        theta_a = math.radians(theta_a)
        cos_theta_a = math.cos(theta_a)
        
        theta_b = 90-DEC_b
        theta_b = math.radians(theta_b)
        cos_theta_b=math.cos(theta_b)
        
        cos_theta = random.uniform(cos_theta_a,cos_theta_b) 
                
        theta = math.acos(cos_theta)    
        theta = math.degrees(theta)
        DEC = 90-theta
        
        return RA,DEC
    
    RA_random=[]
    DEC_random=[]
    for i in range(0,n_vectors,1):
        [RA,DEC]=random_3D_unit_vector(RA_a,RA_b,DEC_a,DEC_b)
        RA_random.append(RA)
        DEC_random.append(DEC)
        
    return RA_random,DEC_random