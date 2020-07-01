# -*- coding: utf-8 -*-
import gro_writer

def experiment_generator(experiment_type, metabolic_schemes,av_force):     
    """
    DEFINICIÓN: Dada una lista de esquemas metabólicos (output del script
                "schemes_creator.py"), genera un archivo con formato .GRO
                para cada esquema
        
    ARGUMENTOS:
        
        - experiment_type: String con el nombre o tipo del experimento que se 
                           quiere realizar
        
        - metabolic_schemes: Lista de esquemas metabólicos (diccionarios). Se 
                             corresponde con el ouput generado con el script
                             "schemes_creator.py"
        
        - av_force: Nº real que se corresponde con el valor medio de fuerza de
                    interacción
        
    RESULTADO: 
        - La función no devuelve nada, directamente genera archivos .GRO, uno
          por cada esquema metabolico de la lista pasada como argumento
    
    """                
                                                                   
    # Lectura de los archivos initial_parameters.GRO y end_parameters .GRO,
    # contienen los parámetros iniciales y finales necesarios para cualquier 
    #archivo .GRO, y almacenamiento de las líneas leidas en la listas 
    # "initial_parameters" y "end_parameters"
    
    parameters_1 = open ("initial_parameters.GRO", "r")                               
    initial_parameters = parameters_1.readlines()
    parameters_1.close()
    
    parameters_2 = open ("end_parameters.GRO", "r")                               
    end_parameters = parameters_2.readlines()
    parameters_2.close()
    
    
    for i in range(len(metabolic_schemes)):                          
        
        # Variable que contiene cada uno de los esquemas metabolicos
        scheme = metabolic_schemes[i]
        
        #Creacion de un archivo.GRO para cada esquema metabolico
        name_experiment=experiment_type+"_experiment"+str(i)   
        file = open(name_experiment+".GRO", "w")                       

        # Escritura del archivo Gro llamando invocando a la función gro_writer
        gro_writer(file, initial_parameters,end_parameters, name_experiment, scheme, av_force)
        file.close()
        
    return

