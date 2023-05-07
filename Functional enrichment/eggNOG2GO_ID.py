# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 11:44:50 2022

@author: Marie Valerie Brasseur
"""

import os
import pandas as pd
import csv

##extract annotation data
os.chdir('C:/Users/mbras/OneDrive/Desktop/Doktorarbeit/Indoor_Genexpression/functional_annotation_eggnog/')


def filter_besthit_bitscore(k, v): ##provide k and v as lists
    replaced_b = []
    best_hit_b = {}

    for i in range(len(k)):
        if k[i].split('.p')[0] not in best_hit_b.keys():
            best_hit_b[k[i].split('.p')[0]] = [v[i], i] #get index for later filtering
            continue
        if k[i].split('.p')[0] in best_hit_b.keys():
            b1 = float(best_hit_b[k[i].split('.p')[0]][0]) 
            b2 = float(v[i])
            if b2 > b1: #higher bitscore is better
                replaced_b.append(k[i])
                best_hit_b[k[i].split('.p')[0]] = [v[i], i]
            else:                
                continue
    return best_hit_b, replaced_b #replaced should be in the final dict

##read in the annotations;
#either:
#annot = pd.read_excel('lepidostoma_trinity/lepidostoma_eggnog_mapper_results.emapper.annotations.xlsx', skiprows=range(2))  

#or:
annot = pd.read_excel('ephemera_trinity/ephemera_eggnog_mapper_results.emapper.annotations.xlsx', skiprows=range(2))  
#annot.head()
#annot.tail()

annot = annot[["query", "evalue", "score", "GOs", "KEGG_ko", "KEGG_Pathway"]]

best_hit, replaced = filter_besthit_bitscore(annot["query"], annot["score"])
#best_hit['TRINITY_DN100373_c0_g1_i2']


idx = [v[-1] for v in list(best_hit.values())]
annot = annot.filter(items=idx, axis= 0)


###gene level
##summarize isoform annotations to gene lvl

genes = [k.split('_i')[0] for k in best_hit.keys()]

len(set(genes)) #84538

iso_ID_maps = annot.set_index('query')['KEGG_Pathway'].to_dict() 
iso_ID_maps = {k:v for (k,v) in iso_ID_maps.items() if 'TRINITY' in k}

k = [k.split('_i')[0] for k in iso_ID_maps.keys()]
v = list(iso_ID_maps.values())
gene_ID_maps = {}
for i in range(len(k)):
    if k[i] not in gene_ID_maps.keys():
        gene_ID_maps[k[i]] = v[i] 
        continue
    if k[i] in gene_ID_maps.keys():
        gene_ID_maps[k[i]] += ',' + v[i] 
        
for k, v in gene_ID_maps.items():
    if ',' in v:
        new_v = list(set(v.split(',')))
    else:
        new_v = list(set(v))
    gene_ID_maps[k] = new_v

with open('Trinity_Gene2KEGG_Pathway.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for k, v in gene_ID_maps.items():
        if v != ['-']:
            IDs = ';'.join([x for x in gene_ID_maps[k] if x != '-'])
            writer.writerow([k+';'+IDs])
    csvfile.close()
        
###ko number
    
iso_ID_maps = annot.set_index('query')['KEGG_ko'].to_dict() 
iso_ID_maps = {k:v for (k,v) in iso_ID_maps.items() if 'TRINITY' in k}

k = [k.split('_i')[0] for k in iso_ID_maps.keys()]
v = list(iso_ID_maps.values())
gene_ID_maps = {}
for i in range(len(k)):
    if k[i] not in gene_ID_maps.keys():
        gene_ID_maps[k[i]] = v[i] 
        continue
    if k[i] in gene_ID_maps.keys():
        gene_ID_maps[k[i]] += ',' + v[i] 
        
for k, v in gene_ID_maps.items():
    if set(v) == {'-'} or set(v.split(',')) == {'-'}:
        new_v = ['-']        
    else:
        new_v = list(set(v.split(',')))
    gene_ID_maps[k] = new_v


with open('Trinity_Gene2KEGG_ko.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for k, v in gene_ID_maps.items():
        if v != ['-']:
            IDs = ';'.join([x for x in gene_ID_maps[k] if x != '-'])
            writer.writerow([k+';'+IDs])
    csvfile.close()

###GO
    
iso_ID_maps = annot.set_index('query')['GOs'].to_dict() 
iso_ID_maps = {k:v for (k,v) in iso_ID_maps.items() if 'TRINITY' in k}

k = [k.split('_i')[0] for k in iso_ID_maps.keys()]
v = list(iso_ID_maps.values())
gene_ID_maps = {}
for i in range(len(k)):
    if k[i] not in gene_ID_maps.keys():
        gene_ID_maps[k[i]] = v[i] 
        continue
    if k[i] in gene_ID_maps.keys():
        gene_ID_maps[k[i]] += ',' + v[i] 
        
for k, v in gene_ID_maps.items():
    if ',' in v:
        new_v = list(set(v.split(',')))
    else:
        new_v = list(set(v))
    gene_ID_maps[k] = new_v


with open('Trinity_Gene2GO.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for k, v in gene_ID_maps.items():
        if v != ['-']:
            IDs = ';'.join([x for x in gene_ID_maps[k] if x != '-'])
            writer.writerow([k+';'+IDs])
    csvfile.close()
    
