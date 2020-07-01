# -*- coding: utf-8 -*-

import experiment_generator
import schemes_creator


#Generador automático de archivos .gro

# Las listas Fi, Ci y Pos contienen los valores de los parámetros de Fuerza de 
# interacción media, Conectividad y Proporción de interacciones positivas y negativas
# que se estudien. Los archivos . gro generados, se crean siguiendo un diseño factorial
# completo --> se estudian crean 10 archivos .gro para cada una de las posibles combinaciones
# de Ci y Pos   y   por otro lado 10 archivos . gro por cada combinación de Fi y Pos


Fi=[0.05,0.2,0.5] 
Ci=[4,8,16]
Pos=[0.25,0.5,0.75]

name1= "grupo1_numPos-"

for i in range(0,len(Ci)):
    for j in range(0,len(Pos)):
        experiment_generator(name1+str(Pos[j])+"_Ci-"+str(Ci[i]), schemes_creator(5,int(Ci[i]),int(Ci[i]*Pos[j]),0.05,10),0.05)


name2= "grupo2_numPos-"

for i in range(0,len(Fi)):
    for j in range(0,len(Pos)):
        experiment_generator(name2+str(Pos[j])+"_Fi-"+str(Fi[i]), schemes_creator(5,16,int(16*Pos[j]),Fi[i],10),Fi[i])
    