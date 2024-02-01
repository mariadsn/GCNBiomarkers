#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 16:07:42 2024

@author: maria
"""

# pruebas con grafos

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

normalBreast = nx.read_edgelist('/home/maria/Documentos/PhD_bired/Preprocesamiento_pruebas/data_paper3/MajorRevision_Enero2024/data_DEGs/final_networks_selected/normalBreast_DEGs_table_MajorRevisionthreshold_0.75.csv', delimiter=(','), nodetype=str)
tumorBreast = nx.read_edgelist('/home/maria/Documentos/PhD_bired/Preprocesamiento_pruebas/data_paper3/MajorRevision_Enero2024/data_DEGs/final_networks_selected/tumorBreast_DEGs_table_MajorRevisionthreshold_0.75.csv', delimiter=(','), nodetype=str)

interaction_common_breast = set(normalBreast.edges()) & set(tumorBreast.edges())  # interacciones compartidas
BreastCommon = nx.Graph(list(interaction_common_breast))

nodoshub_normal = ['DOC2B', 'SH2D3C', 'ARSA','MYH9', 'ANKRD20A11P', 'RGMA', 'FAM110D', 'GRP', 'ADCK2', 'C1orf122', 'LY86', 'FDX2', 'DACT3', 'SPP1', 'SHISA2',
                   'RAPGEF3', 'EVC', 'CCN4', 'FAM110C', 'CCM2L']

nodoshub_tumor = ['LONP2', 'ST6GAL2', 'SUN1', 'GPR137B', 'RIPOR3', 'DYRK2', 'MAP7D1', 'SNRK', 'SYNDIG1', 'MATN3', 'MCTS1', 'SH3BGRL2', 'CLK1', 'COL5A1',
                  'LGALS1', 'ALDH1L2', 'DEPDC7', 'PROS1', 'NAP1L5', 'LRRC59']


interacciones_hub_breast =[]
for nodo in nodoshub_normal:
    if nodo in BreastCommon.nodes:
        for vecino in BreastCommon.neighbors(nodo):
            interacciones_hub_breast.append((nodo,vecino))

interactionNormalBreastExclussive = set(normalBreast.edges()) - set(tumorBreast.edges()) # interacciones exclusivas normales pecho
normalBreastExclussive = nx.Graph(list(interactionNormalBreastExclussive))

interactionTumorBreastExclussive = set(tumorBreast.edges()) - set(normalBreast.edges()) # interacciones exclusivas tumorales pecho
tumorBreastExclussive = nx.Graph(list(interactionTumorBreastExclussive))


normalBreast_gene = set(normalBreast.nodes()) - set(tumorBreast.nodes())


normalProstate = nx.read_edgelist('/home/maria/Documentos/PhD_bired/Preprocesamiento_pruebas/data_paper3/MajorRevision_Enero2024/data_DEGs/final_networks_selected/normalProstate_DEGs_table_MajorRevisionthreshold_0.75.csv', delimiter=(','), nodetype=str)
tumorProstate = nx.read_edgelist('/home/maria/Documentos/PhD_bired/Preprocesamiento_pruebas/data_paper3/MajorRevision_Enero2024/data_DEGs/final_networks_selected/tumorProstate_DEGs_table_MajorRevisionthreshold_0.75.csv', delimiter=(','), nodetype=str)

interaction_common_prostate = set(normalProstate.edges()) & set(tumorProstate.edges()) # interacciones compartidas
ProstateCommon = nx.Graph(list(interaction_common_prostate))


nodoshub_normal = ['TIMP4', 'SMTNL2', 'RBP4', 'NAT2', 'CCDC85A', 'KLHL14', 'IGSF1', 'PDE3B','ARHGAP28','RSPH9','TRHDE','ODAD3','CFD','LOC101927668','LINC01082',
                   'FBXO2','GRTP1-AS1','GATA6-AS1','CLGN','MLC1']
nodoshub_tumor = ['CDH1','ARFGEF3','PRR15L','EHF','SPON2','CRNDE','HOXB13','PDLIM5','CXADR','FOXA1','GATA6-AS1','LEF1-AS1','HOXC6','GPR160','APBA2','SORD',
                  'RAB25','DUSP6','LINC03026','TRPM8']

interacciones_hub_prostate =[]
for nodo in nodoshub_tumor:
    if nodo in ProstateCommon.nodes:
        for vecino in ProstateCommon.neighbors(nodo):
            interacciones_hub_prostate.append((nodo,vecino))

interactionNormalProstateExclussive = set(normalProstate.edges()) - set(tumorProstate.edges()) # interacciones exclusivas normales pecho
normalProstateExclussive = nx.Graph(list(interactionNormalProstateExclussive))

interactionTumorProstateExclussive = set(tumorProstate.edges()) - set(normalProstate.edges()) # interacciones exclusivas tumorales pecho
tumorProstateExclussive = nx.Graph(list(interactionTumorProstateExclussive))


tumorProstate_gene = set(tumorProstate.nodes()) - set(normalProstate.nodes())

# guardar datos del grafos
aristas_normalPecho = list(normalBreastExclussive.edges())
df_aristas_normal = pd.DataFrame(aristas_normalPecho, columns=['Origen', 'Destino'])
df_aristas_normal.to_csv('/home/maria/Documentos/PhD_bired/Preprocesamiento_pruebas/data_paper3/MajorRevision_Enero2024/data_DEGs/final_networks_selected/edges_unico_normal_breast.csv', index=False)

aristas_pecho = list(tumorBreastExclussive.edges())
df_aristas = pd.DataFrame(aristas_pecho, columns=['Origen', 'Destino'])
df_aristas.to_csv('/home/maria/Documentos/PhD_bired/Preprocesamiento_pruebas/data_paper3/MajorRevision_Enero2024/data_DEGs/final_networks_selected/edges_unico_tumor_breast.csv', index=False)

aristas_normalprostate = list(normalProstateExclussive.edges())
df_aristas_normalprostate = pd.DataFrame(aristas_normalprostate, columns=['Origen', 'Destino'])
df_aristas_normalprostate.to_csv('/home/maria/Documentos/PhD_bired/Preprocesamiento_pruebas/data_paper3/MajorRevision_Enero2024/data_DEGs/final_networks_selected/edges_unicos_normal_prostate.csv', index=False)

aristas_prostate = list(tumorProstateExclussive.edges())
df_aristas_prostate = pd.DataFrame(aristas_prostate, columns=['Origen', 'Destino'])
df_aristas_prostate.to_csv('/home/maria/Documentos/PhD_bired/Preprocesamiento_pruebas/data_paper3/MajorRevision_Enero2024/data_DEGs/final_networks_selected/edges_unicos_tumor_prostate.csv', index=False)

#######

nodosInteres = ['TGFB3', 'PTGS1', 'DUSP6', 'SASH1', 'FIGN']
interacciones_normal_breast =[]
for nodo in nodosInteres:
    for vecino in normalBreast.neighbors(nodo):
        interacciones_normal_breast.append((nodo,vecino))

nodosInteres = ['TGFB3', 'PTGS1', 'DUSP6', 'SASH1', 'FIGN']
interacciones_tumor_breast =[]
for nodo in nodosInteres:
    for vecino in tumorBreast.neighbors(nodo):
        interacciones_tumor_breast.append((nodo,vecino))
        
interacciones_comunes_genesCommPros = set(interacciones_normal_breast) & set(interacciones_tumor_breast)  

nodosInteres = ['TGFB3', 'PTGS1', 'DUSP6', 'SASH1', 'FIGN']
interacciones_normal_prostate =[]
for nodo in nodosInteres:
    for vecino in normalProstate.neighbors(nodo):
        interacciones_normal_prostate.append((nodo,vecino))

nodosInteres = ['TGFB3', 'PTGS1', 'DUSP6', 'SASH1', 'FIGN']
interacciones_tumor_prostate =[]
for nodo in nodosInteres:
    for vecino in tumorProstate.neighbors(nodo):
        interacciones_tumor_prostate.append((nodo,vecino))
        
interacciones_comunes_genesCommBres = set(interacciones_normal_prostate) & set(interacciones_tumor_prostate) 


# combinacion de grafos
breastTogether = nx.compose(normalBreast, tumorBreast)
prostateTogether = nx.compose(normalProstate, tumorProstate)

common_breastProstate = set(breastTogether.edges()) & set(prostateTogether.edges())

      

# Aplicar el algoritmo de Girvan-Newman
communities = list(nx.community.girvan_newman(tumorProstateExclussive))

# Obtener la partición de comunidades
best_communities = communities[0]

# Dibujar el grafo con colores según las comunidades
colors = [i for i, comm in enumerate(best_communities) for _ in comm]
nx.draw(tumorProstateExclussive, node_color=colors, with_labels=True, cmap=plt.cm.rainbow)
plt.show()


pos = {node: np.random.rand(3) for node in tumorProstateExclussive.nodes()}        
# Crear un gráfico en 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Obtener posiciones de los nodos
pos = nx.spring_layout(tumorProstateExclussive)

# Dibujar nodos y aristas
nx.draw_networkx_nodes(tumorProstateExclussive, pos, ax=ax, node_color=colors, cmap=plt.cm.rainbow)
nx.draw_networkx_edges(tumorProstateExclussive, pos, ax=ax)

# Ajustar la vista en 3D
ax.view_init(elev=20, azim=30)

# Mostrar el gráfico
plt.show()


# Obtener posiciones aleatorias en 3D
pos = nx.spring_layout(tumorProstateExclussive, dim=3)

# Crear un gráfico en 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Dibujar nodos
nx.draw_networkx_nodes(tumorProstateExclussive, pos, ax=ax, node_color=colors, cmap=plt.cm.rainbow)

# Dibujar aristas en 3D
for edge in tumorProstateExclussive.edges():
    ax.plot([pos[edge[0]][0], pos[edge[1]][0]],
            [pos[edge[0]][1], pos[edge[1]][1]],
            [pos[edge[0]][2], pos[edge[1]][2]], color='gray')

# Ajustar la vista en 3D
ax.view_init(elev=20, azim=30)

# Mostrar el gráfico
plt.show()



# Aplicar el algoritmo de Girvan-Newman
comp = nx.community.girvan_newman(tumorProstateExclussive)

# Obtener el dendrograma de comunidades
dendrogram = list(comp)

# Visualizar el dendrograma
#plt.figure(figsize=(10, 5))
nx.draw(tumorProstateExclussive, with_labels=True, font_weight='bold')
plt.title('Dendrograma de Girvan-Newman')
plt.show()

# Imprimir el dendrograma (esto te dará información sobre la estructura jerárquica)
print("Dendrograma de Girvan-Newman:")
for level, communities in enumerate(dendrogram):
    print(f"Level {level}: {communities}")
        