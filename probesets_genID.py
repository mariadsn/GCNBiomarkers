#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 09:56:25 2024

@author: maria
"""

# MajorRevision Enero 2024 

## affy and max gen value

import pandas as pd
from collections import Counter

affyId = pd.read_csv('aafTableAnn_majorRevision.txt', sep='\t')
expressionLevel = pd.read_csv('rawAffy_expressionLevel_rmaData.csv', sep=',')


frecuencias = Counter(affyId.Symbol)
nameFilas = set(expressionLevel.index)


affyId = affyId.set_index('Probe')
expressionLevel = expressionLevel.set_index('affy_ids')
expressionLevel['Symbol'] = ''


for elem in expressionLevel.index:
    if elem in affyId.index:
        name = affyId['Symbol'][elem]
        expressionLevel['Symbol'][elem] = name
expressionLevel = expressionLevel.set_index('Symbol')    
affyId = affyId.set_index('Symbol')    
columnasDF = expressionLevel.columns

df = pd.DataFrame(columns=columnasDF)
#df['Symbol']=''

descarteGen = []
for gen in expressionLevel.index:
    if gen in affyId.index:
        if gen in frecuencias:
            if gen in descarteGen:
                pass
            else:
                descarteGen.append(gen)
                if frecuencias[gen] ==1:
                    data = expressionLevel.loc[gen]
                    df.loc[gen] = data
                else:
                    data = expressionLevel.loc[gen]
                    data1 = []
                    for x in data.columns:
                        data1.append(max(data[x]))
                    df.loc[gen] = data1
                    #print(gen, frecuencias[gen])

df.to_csv('symbol_expressionLevel_rmaDta.csv',sep='\t',index=True)
