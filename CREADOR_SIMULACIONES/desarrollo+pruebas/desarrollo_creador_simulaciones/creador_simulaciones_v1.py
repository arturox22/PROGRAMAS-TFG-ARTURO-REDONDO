
import pprint, random, string


def schemes_creator(n_strains, n_interactions, positive_interactions, av_force, n_matrixes):
    
    """
    DEFNICIÓN: Función que construye matrices de interacciones aleatorias
               según los parámetros que se le pasan como argumentos. Las 
               interacciones entre especies de cada matriz se representan
               en forma de esquemas metabólicos (diccionarios de python). 
               La función devuelve una lista de esquemas metabólicos, es 
               decir, una lista de diccionarios, uno por cada matriz.
               
    
    ARGUMENTOS:
        - n_strains --> Nº de especies. Este parámetro determinará el tamaño de 
                        la matriz de interacciones (n_strains x n_strains)
                        
        - n_interactions --> Nº de interacciones entre especies distintas. 
                             LAS AUTOINTERACCIONES (Una célula consigo misma)
                             QUE OCUPARÍAN LA DIAGONAL DE LA MATRIZ, NO SE 
                             ESTUDIAN, SE CONSIDERAN 0.
                             
        - positive_interactions --> Nº de interacciones positivas que habrá en la matriz de 
                    interacciones. Las interacciones negativas se obtendrán de 
                    diferencia entre n_interactions y positive_interactions
                             
        - av_force --> Fuerza media de las interacciones  ##### **** 
                        TODO ESTO ES DE AARON:
                            sera la desviación típica de la distribución
                            Distribución normal de la cual extraer las fuerzas
                            de cada interacción
        
        - n_matrixes --> Nª de matrices aleatorias que se quieren generar. Se 
                         corresponde con el nº de archivos GRO que se crearán. 
    
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
               
    """
#El input es el n de especies que va a simular, el n de interacciones y la fuerza media, que sera la desviación típica de la distribución
#Distribución normal de la cual extraer las fuerzas de cada interacción. También el número de matrices, que sera el n de archivos .gro

    n_zeros=n_strains*(n_strains-1)-n_interactions                  #Matriz n*n-1 porque no queremos la diagonal, realmente es un vector
    zeros = [0 for x in range(n_zeros)]                             #Esa matriz es realmente un vector zeros + non-zeros, estos ultimos salen de una normal

#DISEÑO DE LOS ELEMENTOS QUE FORMARAN PARTE DE LA MATRIZ N x N-1

    ##no zeros ahora son el mismo valor con distintos signos
    ##positive_interactions es el n de interacciones con signo positivo
    negative_forces = [-av_force for x in range(int(round(n_interactions-positive_interactions)))]      # nonz --> no zeros 
    positive_forces= [-av_force for x in range(int(round(positive_interactions)))]
    prematrix = zeros+negative_forces+positive_forces                             # Realmente la matriz es una lista de numeros 
    prematrixes_list=[]                                             #Vamos a generar algunas matrices posibles, para ello mezclamos los elementos del la lista matriz
    for i in range(n_matrixes):
                                                                    # con random shuffle --> reorganizas el interior de la lista creada
        random.shuffle(prematrix)                                   #Generamos algunas posibles prematrices por aleatorización. Podría salir repetidas pero es muy improbable. Se podrían hacer permutaciones pero es inviable computacionalmente (Memory Error).
        prematrixes_list.append(list(prematrix))


# LAS MATRICES SE VERAN COMO VECTORES --> A11=0, A12=1, A13=2 ....
# LOS VECTORES/MATRICES SE ALMACENAN EN UNA LISTA DE LISTAS
# EN LAS LINEAS ANTERIORES SE MEZCLA DE FORMA ALEATORIA LOS ELEMENTOS DE LA MATRIZ NxN-1
# EN LAS LINEAS DE A CONTINUACION SE INTRODUCE LA DIAGONAL DE 0 --> NO QUEREMOS AUTOINTERACCIONES
    matrixes_list=[]
    for element in prematrixes_list:                                #Modificamos cada vector/lista para añadir los 0 de la diagonal

        for i in range(0,n_strains**2,n_strains+1):             #Añadir ceros a la diagonal de la matriz (realmente una lista/vector) 
            element.insert(i,0)                             #Dentro de la lista de elementos insertas 0 en las posiciones de lo que seria la diagonal en la matriz
            #pprint.pprint(element)#QUITAR
            matrixes_list.append(element)                           #Se obtiene el conjunto de las matrices completas. tamaño=n_matrices

#GENERACION DE DICCIONARIOS CON LA INFO DEL METABOLISMO QUE SE ASIGNARA A CADA BACTERIA 
# TODOS LOS DICCIONARIOS SE GUARDAN EN UNA LISTA, metabolic_schemes
    metabolic_schemes=[]                                                 # creacion de una lista que almacenara los diccionarios que se crean mas abajo
                                                                    # estos diccionarios almacenaran los metabolitos con los que interactua cada especie
    for matrix in matrixes_list:
            strains_dict={}
            
            # LETTERS ES UN ABECEDARIO QUE SE EMPLEARA PARA DAR NOMBRES A LAS DISTINTAS MOLEC Y SEÑALES
            letters=list(string.ascii_uppercase)                    #Esto permite nombrar luego cada metabolito de interacción con una letra mayúscula
                                                                    # basicamente estamos creando un abecedario de la A a la Z en mayusculas
         
# CON ESTE BUCLE GENERAMOS UNA LISTA DE DICCIONARIOS, DE FORMA QUE 
#CADA DICCIONARIO CORRESPONDE A UNA MATRIZ
# CADA UNO DE ESTOS DICCIONARIOS POSEE A SU VEZ TANTOS DICCIONARIOS 
# COMO ESPECIES HAYA
             

# metabolic_schemes                                                       
#[{1: {'out': ['A', 'B'], 'in': ['C'], 'force': [-0.1]},
#  2: {'out': ['C', 'D'], 'in': ['A'], 'force': [0.1]},
#  3: {'out': [], 'in': ['B', 'D'], 'force': [0.1, 0.1]}},
# {1: {'out': ['A', 'B'], 'in': ['D'], 'force': [0.1]},
#  2: {'out': ['C'], 'in': ['A'], 'force': [0.1]},
#  3: {'out': ['D'], 'in': ['B', 'C'], 'force': [0.1, -0.1]}},
# {1: {'out': ['A', 'B'], 'in': ['D'], 'force': [0.1]},
#  2: {'out': ['C'], 'in': ['A'], 'force': [-0.1]},
#  3: {'out': ['D'], 'in': ['B', 'C'], 'force': [0.1, -0.1]}}]


# MECANISMO MUY COMPLICADO PARA INTERPRETAR MATRIZ, PERO EFECTIVO           
            for i in range(1,n_strains+1):
                    strains_dict[i]={"out":[],"in":[],"force":[]}   #Cada especie tiene una serie de metabolitos que produce(out) y que absorbe(in)


            for j in range(1,n_strains+1):                          #Con cada matriz se procede de este modo para elaborar el esquema metabolico de cada especie
                    for k in range(1,n_strains+1):                  #Las variables j y k son los indices de la matriz con los que establecer el metabolismo de cada especie
                            value=matrix.pop(0)                     #Extraemos el valor correpondiente y si no es 0 se establece interacción y se modifica el metabolismo de las especies corresponientes
                            if value!=0:
                                    metabolite=letters.pop(0)                   #Definimos un metabolito con el instanciador para que medie la interaccion
                                                                                # pop(0) va avanzando en el abecedario uno a uno --> primera iteracion saca A, luego B, luego C
                                    strains_dict[j]["out"].append(metabolite)   #Despues asignamos este metabolito a la lista correspondienta para cada especie
                                    strains_dict[k]["in"].append(metabolite)
                                    strains_dict[k]["force"].append(value)#Se asigna el valor en la fuerza de interacción que modificara la tasa de crecimiento

            metabolic_schemes.append(strains_dict)
    pprint.pprint(metabolic_schemes)
    return metabolic_schemes                                             #Se devuelve el conjunto de diccionarios (esquemas metabólicos), con los que hacer los archivos .gro 

#ejemplo=encoder(3,3,1,3)

#pprint.pprint(ejemplo)


#def encoder(file, initial_parameters, end_parameters, scheme, av_force):
def gro_writer(file, initial_parameters, end_parameters, name_experiment, scheme, av_force):   
 
    # Escritura de los parámetros iniciales del archivo .GRO
    file.writelines(initial_parameters)
    # Escritura del valor de las fuerzas de interacción en el archivo .GRO
    file.write('\nf_pos := '+str(av_force)+';')
    file.write('\nf_neg := -'+str(av_force)+';')
    
    
    ###########################################################################
    
    file.write("""


//-------------------------------------------------------
               
//#################################################
//##################   SIGNAL   ###################
//#################################################


""")


    for strain in scheme.keys():                                   #Cada especie esta asociada a un diccionario, esquema de su metabolismo

        metabolism = scheme[strain]
        # scheme[strain] --> cada una de las especies de cada diccionario--> con los valores out, in, force
        # siendo scheme:  --> cada uno de los elementos de la lista metabolic_schemes, cada scheme es un metabolismo para un GRO distinto
        #   1: {'out': ['A', 'B'], 'in': ['C'], 'force': [-0.1]},
        #   2: {'out': ['C', 'D'], 'in': ['A'], 'force': [0.1]},
        #   3: {'out': [], 'in': ['B', 'D'], 'force': [0.1, 0.1]}
        
        #scheme [1] --> {'out': ['A', 'B'], 'in': ['C'], 'force': [-0.1]}

        #Comenzando por los metabolitos en out, por cada uno hay que codificar una señal
        for letter in metabolism ['out']:
            #### EJEMPLO DE SINTAXIS DE SEÑALES EN GRO 1.2.1g rel
            #signal([ name := "sRock", kdiff := 0.1, kdeg := 0.01, pad := 10.0]);
            signal='signal([ name := "'+letter+'", kdiff := diff, kdeg := deg, pad := k_pad ]);\n'
            file.write(signal)
            
            

    ###########################################################################
    
    
    file.write("""


//-------------------------------------------------------
               
//#################################################
//##################   GATES   ####################
//#################################################


""")
  
    file.write('gate([ name := "gaTrue", function := "TRUE", input := {} ]);\n')
    
    for strain in scheme.keys():
        gate='gate([ name := "gaStrain'+str(strain)+'", function := "YES", input := {"strain'+str(strain)+'"}  ]);\n'
        file.write(gate)
    
    
    ###########################################################################
    
    file.write("""
    
    
//-------------------------------------------------------           
    

//#################################################
//###############   METABOLISM   ##################
//#################################################


""")

    for strain in scheme.keys():
        name_line='\nmetabolism ([name := "metStrain'+str(strain)+'"'
        gate_line=',\n\t\t\tgate := "gaStrain'+str(strain)+'"'                     #IMPORTANTE
        
       
        metabolism = scheme[strain]   ### metabolimo de cada cepa
        fluxes=''
        fluxes_in=''
        fluxes_w=''
        fluxes_lines=[]
        metabolites=''
        
        
        for metabolite in metabolism['out']:
            metabolites='"'+metabolite+'", '+metabolites
            flux=metabolite+'out'
            fluxes='"'+flux+'", '+fluxes
            flux_lines=',\n\n\t\t\t'+flux+' := [metabolite := "'+metabolite+'", f := [ bias := signal_emission ] ]'
            fluxes_lines.append(flux_lines)

        for metabolite in metabolism['in']:
            metabolites='"'+metabolite+'", '+metabolites            #Rellenar flujos
            flux=metabolite+'in'                                    #Nombrar flujos
            fluxes_in='"'+flux+'", '+fluxes_in 
            fluxes='"'+flux+'", '+fluxes                            #Rellenar metabolitos de entrada. Hay un umbral para que comience a absorber.
            flux_lines=',\n\n\t\t\t'+flux+' := [metabolite := "'+metabolite+'",\n\t\t\t\t\t functions := {"ins", "max"},\n\t\t\t\t\t ins := [ metabolites := {"'+metabolite+'"}, metabolites_w := {-0.5}, bias := 0.0 ],\n\t\t\t\t\t max := [ bias := -0.02 ],\n\t\t\t\t\t tree := [ '+metabolite+' := threshold_up, up := "max", low := "ins"]  ]'
            fluxes_lines.append(flux_lines)


        for w in metabolism['force']:                               #No es necesario el redondeo de las fuerzas de interacción
            w=round(w,3)
            if w>0:
                fluxes_w='f_pos'+', '+fluxes_w
            else:
                fluxes_w='f_neg'+', '+fluxes_w
                
                

        if metabolites:
            metabolites=metabolites[:-2]
            
        if fluxes_in:
            fluxes_in=fluxes_in[:-2]
            fluxes_w=fluxes_w[:-2]
            biomass_line=',\n\n\t\t\tbiomass := [ metabolite := "biomass", f := [ fluxes := {'+fluxes_in+'}, fluxes_w := {'+fluxes_w+'}, bias := growth_rate ] ]\n]);\n\n'
            fluxes= fluxes+', "biomass"'        
    
        else:
            biomass_line ='\n]);\n\n'
   
    
        metab_line=',\n\t\t\tmetabolites := {'+metabolites+'}'
        flux_line=',\n\t\t\tfluxes := {'+fluxes+'}'
        
        
    

        file.writelines([name_line,gate_line, metab_line,flux_line])
        file.writelines(fluxes_lines+[biomass_line])
        
        
        
    
    ###########################################################################
    
    file.write("""


//-------------------------------------------------------
               
//#################################################
//##################   COLOR   ####################
//#################################################


""")
    #### Combinaciones posibles de colores --> Posibilidad de pintar 15 cepas 
    #### de 15 colores distintos
    diferent_colors=["0,0,0,color_intensity",                                           #0001
                     "0,0,color_intensity,0",                                           #0010
                     "0,0,color_intensity,color_intensity",                             #0011
                     "0,color_intensity,0,0",                                           #0100
                     "0,color_intensity,0,color_intensity",                             #0101
                     "0,color_intensity,color_intensity,0",                             #0110
                     "0,color_intensity,color_intensity,color_intensity",               #0111
                     "color_intensity,0,0,0",                                           #1000
                     "color_intensity,0,0,color_intensity",                             #1001
                     "color_intensity,0,color_intensity,0",                             #1010
                     "color_intensity,0,color_intensity,color_intensity",               #1011
                     "color_intensity,color_intensity,0,0",                             #1100
                     "color_intensity,color_intensity,0,color_intensity",               #1101
                     "color_intensity,color_intensity,color_intensity,0",               #1110
                     "color_intensity,color_intensity,color_intensity,color_intensity"] #1111
                     
    
    for strain in scheme.keys():  
        color= 'color([ gate := "gaStrain'+str(strain)+'", channels := {'+diferent_colors[strain - 1]+'}, delta := delta_on, reverse := reverse_on ]);\n'
        file.write(color)
        
        
    
    
    
    ###########################################################################

    file.write("""


//-------------------------------------------------------
               
//#################################################
//#################   STRAIN   ####################
//#################################################


""")
    for strain in scheme.keys():  
        
        cepa= 'strain([ name := "Strain'+str(strain)+'", cell_growth := [ base_growth_rate := base_growth, metabolism_growth := met_growth_on ] ]);\n'
        file.write(cepa)
    
    
        
    
    ###########################################################################
    
    file.write("""


//-------------------------------------------------------
               
//#################################################
//#############  REPLATING FUNCTION  ##############
//#################################################


""")
    
    replating= 'replating([ alive_fraction := survivors, alive_fraction_var := survivors_variability, timing := replate_timing ]);\n'
    file.write(replating)
        
            
    
    
    ###########################################################################
    
    file.write("""


//-------------------------------------------------------
               
//#################################################
//##################   OUTPUT   ###################
//#################################################


""")
    #Creamos una string con los nombres de las cepas
    strains=''
    for strain in scheme.keys():  
        strains+='"Strain'+str(strain)+'", '
    
    #Del string anterior, eliminamos la coma y el espacio final
    strains=strains[:-2]
    
    #Generamos el output con el conteo de todass las cepas
    output='output([ path := output_path, file_name := '+name_experiment+ ',\n\t\t\t\t\t\tgate := "gaTrue",\n\t\t\t\t\t\tpopulation_level := 1,\n\t\t\t\t\t\ttiming := output_timing,\n\t\t\t\t\t\tdecimal_places := 2,\n\t\t\t\t\t\tfields := {"molecule"}, molecule := {'+strains+'}]);\n'
    file.write(output)
        

    
        
    ###########################################################################
    ########   ESCRITURA DE LOS PARÁMETROS FINALES DEL ARCHIVO .GRO     #######
    ###########################################################################
    
    file.writelines(end_parameters)
    
    
    ###########################################################################
    
    file.write("""


//-------------------------------------------------------
               
//#################################################
//##################   ECOLIS   ###################
//#################################################


""")
    for strain in scheme.keys():  
        ecolis= 'ecolis( [ num := num_cells, x := x_position, y := y_position, r := radious , strain := "Strain'+str(strain)+'", plasmids := {}, molecules := {}, mode := "default" ], program p());\n'
        file.write(ecolis)
        
        
    return


#def programer(strain_schemes,av_force):
def experiment_generator(experiment_type, metabolic_schemes,av_force):                     
                                                                   

    parameters_1 = open ("initial_parameters.GRO", "r")                               #Leemos el modelo del cual sacar el comienzo y el final del archivo .gro
    initial_parameters = parameters_1.readlines()
    parameters_1.close()
    
    parameters_2 = open ("end_parameters.GRO", "r")                               #Leemos el modelo del cual sacar el comienzo y el final del archivo .gro
    end_parameters = parameters_2.readlines()
    parameters_2.close()
    
#    initial_parameters=str(initial_parameters)
#    end_parameters=str(end_parameters)
    
    for i in range(len(metabolic_schemes)):                           #Cada esquema es un diccionario con n especies
        
        # Variable que contiene cada uno de los esquemas metabolicos
        scheme = metabolic_schemes[i]
        
        #Creacion de un archivo.GRO para casa esquema metabolico
#        exact_time = datetime.today().strftime("_%d-%m-%y_%I:%m_%p")
        name_experiment=experiment_type+"_experiment"+str(i)   #+exact_time
        file = open(name_experiment+".GRO", "w")                       #Se genera un archivo .gro para cada diseño

        # Escritura del archivo Gro 
        gro_writer(file, initial_parameters,end_parameters, name_experiment, scheme, av_force)
        file.close()
        
# codigo de aaron que hay que modificar
#        model_endmod=list(model_end)                                #Generamos una nueva lista (evitar el ) para añadir el fopen y escribirlo después
#        outputfile='fp := fopen (route1 <> "'+name_experiment+'.dat", "w" );\n'
#        model_endmod.insert(15,outputfile)
        
        
        
        

    return

    





# HAY UN ARCHIVO MODELO ..> MODELO2.GRO 
# ESTE ARCHIVO SE DIVIDE EN DOS --> PRIMERA PARTE CON PARAMETROS COMUNES
# SEGUNDA PARTE CON METABOLISMOS



# metabolic_schemes --> lista que contiene los diccionarios para cada matriz
# scheme son los elementos de metabolic_schemes --> scheme == metabolic_schemes[i]

#gro_creator(schemes_creator(5,14,0.5,10,14),0.5)
                        











