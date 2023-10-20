#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 11:38:12 2023

@author: maria
"""

# carga de netwoks

import networkx as nx

networkProstate = nx.read_edgelist('/home/maria/Documentos/PhD_bired/Preprocesamiento_pruebas/data_paper3/uniqueNetwork_tumorProstate_010_08.csv', delimiter=(','), nodetype=str)
networkBreast = nx.read_edgelist('/home/maria/Documentos/PhD_bired/Preprocesamiento_pruebas/data_paper3/cancer_network085/breast_cancer_network.csv', delimiter=(','), nodetype=str)


# calculo de la densidad
densityProstate = nx.density(networkProstate)
densityBreast = nx.density(networkBreast)
print(f"Densisdad de la red de pecho: {densityBreast}. Densidad de la red de prostata: {densityProstate}")


# calculo de coeficiente de clusterizacion
clusteringCoefficient_prostate = nx.clustering(networkProstate)
avgClusteringCoeficient_prostate = nx.average_clustering(networkProstate)
clusteringCoefficient_breast = nx.clustering(networkBreast)
avgClusteringCoeficient_breast = nx.average_clustering(networkBreast)
print(f"Coeficiente de clusterizacion promedio de pecho: {avgClusteringCoeficient_breast}. Coeficiente de clusterizacion promedio de prostata: {avgClusteringCoeficient_prostate}")

# calculo de la longitud de caminos minimos -- no encuentro genes interesantes donde calcularlo

# calculo de la modularidad
modularity_prostate = nx.community.modularity(networkProstate, nx.community.label_propagation_communities(networkProstate))
modularity_Breast = nx.community.modularity(networkBreast, nx.community.label_propagation_communities(networkBreast))
print(f"La modularidad de la red de pecho es: {modularity_Breast}. La modularidad de la red de prostata es: {modularity_prostate}")

# calculo de la heterogeneidad de grado -- a partir del coeficiente de asortividad
degree_assortativity_prostate = nx.degree_assortativity_coefficient(networkProstate)
degree_assortativity_breast = nx.degree_assortativity_coefficient(networkBreast)
print(f"Coeficiente de asortatividad de pecho: {degree_assortativity_breast}. Coeficiente de asortatividad de prostata: {degree_assortativity_prostate} ")

# conectividad entre modulos -- se hace entre modulos

