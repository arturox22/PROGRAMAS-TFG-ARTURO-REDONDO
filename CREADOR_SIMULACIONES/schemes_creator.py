# -*- coding: utf-8 -*-

import string, random

def schemes_creator(n_strains, n_interactions, positive_interactions, av_force, n_matrixes):
    """
    DEFINICIÓN: Función que construye matrices de interacciones aleatorias
               según los parámetros que se le pasan como argumentos. Las 
               interacciones entre especies de cada matriz se representan
               en forma de esquemas metabólicos (diccionarios de python). 
               La función devuelve una lista de esquemas metabólicos, es 
               decir, una lista de diccionarios, uno por cada matriz.
               
    
    ARGUMENTOS:
        - n_strains --> Nº entero. Nº de especies. Este parámetro determinará el tamaño de 
          (INTEGER)      la matriz de interacciones (n_strains x n_strains)
                        
        - n_interactions --> Nº de interacciones entre especies distintas. 
          (INTEGER)          LAS AUTOINTERACCIONES (Una célula consigo misma)
                             QUE OCUPARÍAN LA DIAGONAL DE LA MATRIZ, NO SE 
                             ESTUDIAN, SE CONSIDERAN 0.
                             
        - positive_interactions (INTEGER) --> Nº de interacciones positivas que 
                        habrá en la matriz de interacciones. Las interacciones negativas 
                        se obtendrán de diferencia entre n_interactions y positive_interactions
                             
        - av_force --> (REAL) Fuerza media de las interacciones  ##### **** 
                        TODO ESTO ES DE AARON:
                            sera la desviación típica de la 
                            Distribución normal de la cual extraer las fuerzas
                            de cada interacción
        
        - n_matrixes --> Nª de matrices aleatorias que se quieren generar. Se 
             (INTEGER)   corresponde con el nº de archivos GRO que se crearán. 
    
    RESULTADO: La funcion devuelve la lista metabolic_schemes
        
        - metabolic_schemes --> lista de diccionarios que contendrá los esquemas 
                           metabólicos asociados a cada una de las matrices 
                           generadas. Tiene la siguiente forma:
                             
Ejemplo de metabolic_schemes para:
    n_strains = 3 
    n_interactions = 4 
    positive_interactions = 3
    av_force = 0.1 
    n_matrixes = 3
                                                    
[
    {    1: {'out': ['A', 'B'], 'in': ['C'], 'force': [-0.1]},
         2: {'out': ['C', 'D'], 'in': ['A'], 'force': [0.1]},
         3: {'out': [], 'in': ['B', 'D'], 'force': [0.1, 0.1]}     },
    
    {    1: {'out': ['A', 'B'], 'in': ['D'], 'force': [0.1]},
         2: {'out': ['C'], 'in': ['A'], 'force': [0.1]},
         3: {'out': ['D'], 'in': ['B', 'C'], 'force': [0.1, -0.1]}  },
    
    {    1: {'out': ['A', 'B'], 'in': ['D'], 'force': [0.1]},
         2: {'out': ['C'], 'in': ['A'], 'force': [-0.1]},
         3: {'out': ['D'], 'in': ['B', 'C'], 'force': [0.1, -0.1]}  }
]       


    Cada elemento de la lista se corresponde con el esquema metabólico 
    (diccionario) asignado a cada una de las matrices creadas. 
    
    Dentro de cada esquema (diccionario), los nº se corresponden con 
    especies distintas. A su vez, cada especie (cada nº) tiene dentro
    un diccionario con las interacciones que ésta presenta con el resto 
    de especies.
    
    Este ejemplo en concreto se correspondería con la matriz: 
                        
                      0      0.1     0.1
                    -0.1      0      0.1
                      0       0       0
        
        * Los elementos de diagonal principal siempre serán ceros, 
          no se consideran las autointeracciones
        
    Las matrices van a ser tratadas como vectores. Para referirnos al 
    elemento ij, usaremos la posicion que ocuparía dicho elemento si 
    todos los elementos de la matriz se representasen seguidos en línea
    (de izq a decha) 
    
        Una matriz 3 x 3 se representaría de la siguiente manera:
            
                [a11, a12, a13, a21, a22, a23, a31, a32, a33]
                
     Posiciones:   0   1    2    3    4    5    6    7    8
    
    
    """

########################################################################   
# CREACIÓN Y DEFINICIÓN DE LOS ELEMENTOS DE UNA MATRIZ NO CUADARADA DE #
#            DIMENSIONES     N * (N-1) --> "PREMATRIZ"                 # 
#       (Los elementos de la diagonal se añadirán más adelante)        # 
######################################################################## 
    
    
    # 1) Definición del numero de 0 de la matriz:
    #(nº total de elementos de la matriz NO cuadrada - nº total interacciones)
    
    n_zeros=n_strains*(n_strains-1)-n_interactions
    zeros = [0 for x in range(n_zeros)]                             
    
    
    # 2) Definición de las interacciones positivas y negativas
    
    negative_forces = [-av_force for x in range(int(round(n_interactions-positive_interactions)))]      
    positive_forces= [av_force for x in range(int(round(positive_interactions)))]
    
    
    # 3) Definicion de la lista interactions que contendrá los elementos de 
    # la matriz N x N-1 (*** En este paso el orden de los elementos no importa)
    
    interactions = zeros+negative_forces+positive_forces     


    # 4) Definición de la lista que almacenará todas las prematrices que se 
    #    generarán en el siguiente paso   
                         
    prematrixes_list=[]     
    
    
    # 5) Creación de las prematrices que serán guardadas en prematrixes_list
    # Las matrices se generan por aleatorización los elementos de las lista
    # "interactions" anteriormente creada
    
    for i in range(n_matrixes):                 
        random.shuffle(interactions)                                  
        prematrixes_list.append(list(interactions))
        
        
    
########################################################################   
#       DEFINICIÓN DE LA MATRIZ DE INTERACCIONES DE DIMENSIONES        #
#                    DIMENSIONES   N * (N-1)                           # 
#          (Ahora se añaden los elementos de la diagonal               # 
######################################################################## 
    
    # 1) Creación de matrixes_list, una lista que almacenerá todas las matrices
    #    generadas
    
    matrixes_list=[]
    
    
    # 2) Adición de la diagonal de 0 para conseguir la matriz N x N
    # Para ello se recorre prematrixes_list y en cada prematriz se 
    # añade 0 en las posiciones que ocuparía la diagonal en un vector
    
    # Las posiciones de la diagonal se consiguen con la fórmula: 
    #          range(0,n_strains**2,n_strains+1)
    
    for element in prematrixes_list:                                
        for i in range(0,n_strains**2,n_strains+1):             
            element.insert(i,0)   
        # Una vez tenemos la matriz N x N la añadimos a matrixes_list    
        # NO OLVIDAR QUE:   len (matrixes_list) == n_matrixes           
        matrixes_list.append(element) 

    #print(matrixes_list)


########################################################################   
#       DEFINICIÓN DE LOS ESQUEMAS METABOLICOS ASOCIADOS A CADA        #
#             MATRIZ (se representan como diccionarios)                #
#      Estos esquemas interpretan la posicion Ajk como el efecto de    #
#             de la especie j sobre k (al contrario que May)           # 
######################################################################## 
    
    # 1) Definicion de la lista que almacenará todos los esquemas metabólicos
    metabolic_schemes=[]     

    # 2) Para cada matriz se realiza:
    for matrix in matrixes_list:
        
        # a) Creación del esquema
        strains_dict={}
        
        # b) Creación de una lista con las letras del abecedario que se usará
        #    para nombrar a las señales que emiten las cepas
        letters=list(string.ascii_uppercase)                    

        # c) Definición de los flujos de entrada y salida de cada cepa
        for i in range(1,n_strains+1):
            strains_dict[i]={"out":[],"in":[],"force":[]}   

        #d) Relleno del esquema  --> proceso bastante complejo
        #      - Las variables j y k son los indices de la matriz con 
        #        los que establecer el metabolismo de cada especie
        
        for j in range(1,n_strains+1):                          
            for k in range(1,n_strains+1):                  
                
                # Asignación de la letra correspondiente a cada cepa
                value=matrix.pop(0)                     
                if value!=0:
                    metabolite=letters.pop(0)                   
                    
                    # Relleno de los flujos
                    strains_dict[j]["out"].append(metabolite)   
                    strains_dict[k]["in"].append(metabolite)
                    strains_dict[k]["force"].append(value)

        
        #e) Añadimos el esquema metabolico relleno a la lista de esquemas 
        metabolic_schemes.append(strains_dict)
        

    #Se devuelve el conjunto de diccionarios (esquemas metabólicos), con los que hacer los archivos .gro
    return metabolic_schemes                                              
