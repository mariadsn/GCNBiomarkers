#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 10:41:11 2023

@author: maria
"""

## extraer sin duplicar de los genes normales y tumorales

import pandas as pd

activated = pd.read_csv('activatedGenes_prostate_table_010.txt', sep='\t')
repressed = pd.read_csv('repressedGenes_prostate_table_010.txt', sep='\t')

duplicated = []
f = open('tumorProstate_DEGs_table_010_symbol.csv', 'w')

data_normal = open('tumorProstate_DEGs_table_010.csv', 'r')
for line in data_normal:
    if line.startswith('attr_name'):
        f.write(line)
    spl = line.split(',')
    for i in activated.index:
        if activated.Probe[i] == spl[0]:
            if activated.Symbol[i] == "":
                print('yes')
                pass
            else:
                if not activated.Symbol[i] in duplicated:
                    print('activo')
                    duplicated.append(activated.Symbol[i])
                    f.write(str(activated.Symbol[i]) + ',' + str(spl[1]) + ',' + str(spl[2]) + ',' + str(spl[3]) + ',' + str(spl[4]) + ',' + str(spl[5]) + ',' + str(spl[6]))
    for i in repressed.index:
        if repressed.Probe[i] == spl[0]:
            if repressed.Symbol[i] == "":
                pass
            else:
                if not repressed.Symbol[i] in duplicated:
                    print('reprimido')
                    duplicated.append(repressed.Symbol[i])
                    f.write(str(repressed.Symbol[i]) + ',' + str(spl[1]) + ',' + str(spl[2]) + ',' + str(spl[3]) + ',' + str(spl[4]) + ',' + str(spl[5]) + ',' + str(spl[6]))    

f.close()
data_normal.close()        