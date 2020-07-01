# -*- coding: utf-8 -*-

def gro_writer(file, initial_parameters, end_parameters, name_experiment, scheme, av_force):   
 
    """
    DEFINICIÓN: Función que a partir de un esquema metabólico (generado con 
                "schemes_creator.py") escribe un archivo con formato .GRO
    
    **** Realmente esta función, es una subfunción del script  ****
         experiment_generator.py, pero dada su gran extensión, 
                  se ha escrito en un módulo a parte.
                
    ARGUMENTOS:
        
        - file: Archivo en el que se escribirá la traducción del esquema
                metabólico a lenguaje GRO.
                
        - initial_parameters: Lista que contendrá las parámetros iniciales 
                              necesarios para cualquier simulacion con GRO.
                              
        - end_parameters: Lista que contendrá las parámetros finales 
                          necesarios para cualquier simulacion con GRO.
                
        *** Estas dos listas de parámetros se generan automáticamente ***
                              al correr el script
                              
        - name_experiment: String que designa al nombre del experimento 
                           (necesario para diferenciar los distintos ouputs
                             que se generarán al correr el archivo.GRO)
        
        - scheme: Esquema metabólico (diccionario) --> cada uno de los 
                  elemento de la lista de diccionarios generada con el
                  script "schemes_creator"
        
        - av_force: Nº real que designa el valor absoluto de la fuerza de
                    interacción entre las distintas cepas
        
        
    RESULTADOS:
        
        - Esta función no devuelve nada, directamente escribe en el archivo
          que se le pasa como parámetro
    
    """
    
    
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
        # Siendo scheme cada uno de los elementos de la lista metabolic_schemes,
        # con la siguiente forma: 
        
        #   1: {'out': ['A', 'B'], 'in': ['C'], 'force': [-0.1]},
        #   2: {'out': ['C', 'D'], 'in': ['A'], 'force': [0.1]},
        #   3: {'out': [], 'in': ['B', 'D'], 'force': [0.1, 0.1]}
        
        # scheme[strain] se corresponde con cada una de las cepas 
        # representadas con los nº 1 ,2, 3 en el ejemplo
        
        # Por ejemplo, la cepa 1 tendría el siguiente metabolismo
        #scheme [1] --> {'out': ['A', 'B'], 'in': ['C'], 'force': [-0.1]}
        
        # *** FORCE siempre se refiere a la fuerza de interacción, ***
        #      es decir, se relaciona con el flujo de entrada (IN)

        
        for letter in metabolism ['out']:
        
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
        gate='gate([ name := "gaStrain'+str(strain)+'", function := "YES", input := {"Strain'+str(strain)+'"}  ]);\n'
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
        gate_line=',\n\t\t\tgate := "gaStrain'+str(strain)+'"'                     
        
       
        # metabolismo de cada cepa
        metabolism = scheme[strain]   
       
        # String compuesto por los nombres de TODOS los flujos metabólicos de 
        # la cepa que se este analizando
        fluxes_names=''
        
        # String compuesto por los nombres de todos los flujos de entrada de 
        # la cepa que se este analizando
        fluxes_in=''
        
        # String compuesto por los coeficientes con los que se definirá la 
        # función que relacione la absorción de señal con el crecimiento de 
        # la cepa que estemos estudiando
        fluxes_w=''
        
        # Lista que almacenará las líneas en lenguaje .GRO asociadas 
        # a los flujos metabólicos de la cepa analizada
        fluxes_lines=[]
        
        # String compuesto por los nombres de los metabolitos con los que la 
        # cepa este relacionada
        metabolites=''
        
        
        for metabolite in metabolism['out']:
            metabolites=metabolites+', "'+metabolite+'"'  
            actual_flux_name=metabolite+'out'
            fluxes_names=fluxes_names+', "'+actual_flux_name+'"'
            actual_flux_line=',\n\n\t\t\t'+actual_flux_name+' := [metabolite := "'+metabolite+'", f := [ bias := signal_emission ] ]'
            fluxes_lines.append(actual_flux_line)

        for metabolite in metabolism['in']:
            metabolites=metabolites+', "'+metabolite+'"'            
            actual_flux_name=metabolite+'in'                                    
            fluxes_in=fluxes_in+', "'+actual_flux_name+'"'
            fluxes_names=fluxes_names+', "'+actual_flux_name+'"'                            
            actual_flux_line=',\n\n\t\t\t'+actual_flux_name+' := [metabolite := "'+metabolite+'",\n\t\t\t\t\t functions := {"ins", "max"},\n\t\t\t\t\t ins := [ metabolites := {"'+metabolite+'"}, metabolites_w := {-0.5}, bias := 0.0 ],\n\t\t\t\t\t max := [ bias := -threshold_up ],\n\t\t\t\t\t tree := [ '+metabolite+' := threshold_up, up := "max", low := "ins"]  ]'
            fluxes_lines.append(actual_flux_line)


        # Redondeo de las fuerzas de interacción
        for w in metabolism['force']:                               
            w=round(w,3)
            if w>0:
                fluxes_w=fluxes_w+', -f_pos'
            else:
                fluxes_w=fluxes_w+', -f_neg'
                
                
        # Con las sentencias [2:] eliminamos los dos primeros elementos de los
        # strings a los que se refiere, pues dichos strings han sido generados
        # por concatenación y tienen un ",/n" extra que daría un error al. 
        
        fluxes_names=fluxes_names[2:]
        
        if metabolites:
            metabolites=metabolites[2:]
            
        if fluxes_in:
            fluxes_in=fluxes_in[2:]
            fluxes_w=fluxes_w[2:]
            biomass_line=',\n\n\t\t\tbiomass := [ metabolite := "biomass", f := [ fluxes := {'+fluxes_in+'}, fluxes_w := {'+fluxes_w+'}, bias := 0.0 ] ]\n]);\n\n'
            fluxes_names= fluxes_names+', "biomass"'        
        
        # Si no hay flujos de entrada, no hace falta definir el flujo biomass,
        # pues en las lineas de "strain", ya se define una tasa de crecimiento
        # basal
        else:
            biomass_line ='\n]);\n\n'
   
        
        all_metabolites=',\n\t\t\tmetabolites := {'+metabolites+'}'
        all_fluxes=',\n\t\t\tfluxes := {'+fluxes_names+'}'
        
        
        # Escrbimos en ORDEN todas lineas generadas
        file.writelines([name_line,gate_line, all_metabolites,all_fluxes])
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
        
        cepa= 'strain([ name := "Strain'+str(strain)+'", cell_growth := [ base_growth_rate_rdn := [ dist_params := {base_growth, 0.0} ], metabolism_growth := met_growth_on ] ]);\n'
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
    output='output([ path := output_path, file_name := "'+name_experiment+'",\n\t\t\t\t\t\tgate := "gaTrue",\n\t\t\t\t\t\tpopulation_level := 1,\n\t\t\t\t\t\ttiming := output_timing,\n\t\t\t\t\t\tdecimal_places := 2,\n\t\t\t\t\t\tfields := {"molecule"}, molecule := {'+strains+'}]);\n'
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