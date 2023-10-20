import csv 
import numpy as np
import scipy.stats as scp
import networkx as nx
import matplotlib.pyplot as plt
import os
import time

tiempo_inicio = time.time()

# arista que contiene correlacion entre dos nodos
class Arista:
    def __init__(self, genA, genB):
        self.genA = genA
        self.genB = genB

# grafica - red genetica
def showGraph(graph,title=None):
    pos = nx.spring_layout(graph,k=100)
    nx.draw_networkx_nodes(graph, pos, cmap=plt.get_cmap('jet'), node_size = 500)
    nx.draw_networkx_labels(graph, pos)
    nx.draw_networkx_edges(graph, pos, arrows=False)
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=nx.get_edge_attributes(graph, 'weight'))
    plt.title(title)
    plt.show()
 
path = "/home/maria/Documentos/PhD_bired/Preprocesamiento_pruebas/data_paper3/scriptPython/data_forPython2/" # ruta que contiene los archivos de los DEG  
archivos = os.listdir(path)

for archivo in archivos:
    tiempo_inicio2 = time.time()
    print('archivo: {}'.format(archivo))
    if archivo.endswith(".csv"):
        experimento = archivo.split(".")  
     
        # Leer el fichero csv
        with open(path + archivo) as File:
            data = list(csv.reader(File, delimiter=',', quotechar=',', quoting=csv.QUOTE_MINIMAL)) # datos del fichero se obtienen como una lista
            # quitar fila 0 y guardarla en una lista
            ids = data.pop(0)
            threshold_1 = 0.7
            threshold_2 = 0.75
            threshold_3 = 0.8
            threshold_4 = 0.85
            threshold_5 = 0.9
            threshold_6 = 0.95
            correlation = [] # lista que guarda los nombres de los genes cuya correlacion supera el umbral y su valor
            fichero07 = [("Origen", "Destino")] #lista para volcar en el fichero
            fichero075 = [("Origen", "Destino")]
            fichero08 = [("Origen", "Destino")]
            fichero085 = [("Origen", "Destino")]
            fichero09 = [("Origen", "Destino")]
            fichero095 = [("Origen", "Destino")]
            comparedGens = [] # lista que guarda los nombres de los genes que ya se han comparado
            for gen_actual in data: # para cada gen        
                for gen_comparado in data: # se compara con el resto de genes 
                    nombre_gen_actual = gen_actual[0] # el nombre del gen es el primer elemento de la lista
                    nombre_gen_comparado = gen_comparado[0]
                    
                    # para que no se compare consigo mismo
                    if nombre_gen_actual == nombre_gen_comparado: continue
                    # no se haga una comparación ya hecha
                    if (nombre_gen_comparado, nombre_gen_actual) in comparedGens: continue
        
                    # cálculo de correlacion de Spearman entre dos genes 
                    corr_Spearman, pv = scp.spearmanr(gen_actual, gen_comparado)
                    corr_Spearman = abs(corr_Spearman) 
        
                    # cálculo de la correlacion de Kendall entre dos genes
                    corr_Kendall, pv = scp.kendalltau(gen_actual, gen_comparado)
                    corr_Kendall = abs(corr_Kendall)
        
                    # cálculo de la correlación de Pearson
                    x = [float(value) for value in gen_actual[1:]]
                    y = [float(value) for value in gen_comparado[1:]]
                    corr_Pearson, pv = scp.pearsonr(x,y)
                    corr_Pearson = abs(corr_Pearson)
        
                    # guardar en la lista comparedGens la comparación actual --> evitar que dos pares de genes se comparen mas de una vez
                    comparedGens.append((nombre_gen_actual, nombre_gen_comparado))
        
                    # major voting = si dos de las correlaciones son mayor o igual que el umbral, se guarda en la lista correlation
                    if corr_Spearman >= threshold_1 and corr_Kendall >= threshold_1 or corr_Spearman >= threshold_1 and corr_Pearson >= threshold_1 or corr_Kendall >= threshold_1 and corr_Pearson >= threshold_1:
                        correlation.append(
                            Arista(
                                    nombre_gen_actual, nombre_gen_comparado, # arista = genes que componen los nodos de una relacion
                                )
                            )
                        fichero07.append((nombre_gen_actual, nombre_gen_comparado)) #lista que contiene los datos para volcar en el fichero
                    
                    if corr_Spearman >= threshold_2 and corr_Kendall >= threshold_2 or corr_Spearman >= threshold_2 and corr_Pearson >= threshold_2 or corr_Kendall >= threshold_2 and corr_Pearson >= threshold_2:
                        correlation.append(
                            Arista(
                                    nombre_gen_actual, nombre_gen_comparado, # arista = genes que componen los nodos de una relacion
                                )
                            )
                        fichero075.append((nombre_gen_actual, nombre_gen_comparado)) #lista que contiene los datos para volcar en el fichero
                
                    if corr_Spearman >= threshold_3 and corr_Kendall >= threshold_3 or corr_Spearman >= threshold_3 and corr_Pearson >= threshold_3 or corr_Kendall >= threshold_3 and corr_Pearson >= threshold_3:
                        correlation.append(
                            Arista(
                                    nombre_gen_actual, nombre_gen_comparado, # arista = genes que componen los nodos de una relacion
                                )
                            )
                        fichero08.append((nombre_gen_actual, nombre_gen_comparado)) #lista que contiene los datos para volcar en el fichero
                
                    if corr_Spearman >= threshold_4 and corr_Kendall >= threshold_4 or corr_Spearman >= threshold_4 and corr_Pearson >= threshold_4 or corr_Kendall >= threshold_4 and corr_Pearson >= threshold_4:
                        correlation.append(
                            Arista(
                                    nombre_gen_actual, nombre_gen_comparado, # arista = genes que componen los nodos de una relacion
                                )
                            )
                        fichero085.append((nombre_gen_actual, nombre_gen_comparado)) #lista que contiene los datos para volcar en el fichero
                
                    if corr_Spearman >= threshold_5 and corr_Kendall >= threshold_5 or corr_Spearman >= threshold_5 and corr_Pearson >= threshold_5 or corr_Kendall >= threshold_5 and corr_Pearson >= threshold_5:
                        correlation.append(
                            Arista(
                                    nombre_gen_actual, nombre_gen_comparado, # arista = genes que componen los nodos de una relacion
                                )
                            )
                        fichero09.append((nombre_gen_actual, nombre_gen_comparado)) #lista que contiene los datos para volcar en el fichero
                
                    if corr_Spearman >= threshold_6 and corr_Kendall >= threshold_6 or corr_Spearman >= threshold_6 and corr_Pearson >= threshold_6 or corr_Kendall >= threshold_6 and corr_Pearson >= threshold_6:
                        correlation.append(
                            Arista(
                                    nombre_gen_actual, nombre_gen_comparado, # arista = genes que componen los nodos de una relacion
                                )
                            )
                        fichero095.append((nombre_gen_actual, nombre_gen_comparado)) #lista que contiene los datos para volcar en el ficher
                
                    #print('{} se está comparando con {}'.format(gen_actual[0], gen_comparado[0]))
                    #print('La correlacion de S es {}, la de K es {} y la de P es {}\n'.format(corr_Spearman, corr_Kendall, corr_Pearson))
            
            tiempo_fin2 = time.time()
            tiempo_total2 = tiempo_fin2 - tiempo_inicio2
            print('fin {} tardo {}'.format(archivo, tiempo_total2))
            
            # fichero que contenga origen-destino
            with open(path + '{}threshold_{}.csv'.format(experimento[0], threshold_1), 'w', newline='') as f:
                writer = csv.writer(f)
                for elemento in fichero07:
                    writer.writerow([elemento])
        
            with open(path + '{}threshold_{}.csv'.format(experimento[0], threshold_2), 'w', newline='') as f:
                writer = csv.writer(f)
                for elemento in fichero075:
                    writer.writerow([elemento])
        
            with open(path + '{}threshold_{}.csv'.format(experimento[0], threshold_3), 'w', newline='') as f:
                writer = csv.writer(f)
                for elemento in fichero08:
                    writer.writerow([elemento])
        
            with open(path + '{}threshold_{}.csv'.format(experimento[0], threshold_4), 'w', newline='') as f:
                writer = csv.writer(f)
                for elemento in fichero085:
                    writer.writerow([elemento])
        
            with open(path + '{}threshold_{}.csv'.format(experimento[0], threshold_5), 'w', newline='') as f:
                writer = csv.writer(f)
                for elemento in fichero09:
                    writer.writerow([elemento])
        
            with open(path + '{}threshold_{}.csv'.format(experimento[0], threshold_6), 'w', newline='') as f:
                writer = csv.writer(f)
                for elemento in fichero095:
                    writer.writerow([elemento])

tiempo_fin = time.time()
tiempo_total = tiempo_fin - tiempo_inicio
print('El script tardo {}'.format(tiempo_total))    
    
    