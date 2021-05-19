# -*- coding: utf-8 -*-
"""
Created on Tue May 18 15:40:27 2021

@author: floor
"""
import pandas as pd
from statistics import mean
import numpy as np
import csv
import matplotlib as plt

#genes_start_stop = fitness_eachgene(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\Processed_data_WT1.wig" )
genes_start_stop = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
genes_start_stop['fitness'] = ''
genes_start_stop['fitness_avr'] = ''
genes_start_stop['fitness_variance'] = ''
genes_start_stop['fitness_diff'] = ''

genes_start_stop['fitness']= genes_start_stop['fitness'].astype('object')

with open('fitness_to_HO (1).csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data = list(reader)
genes_start_stop['fitness'] = data

for i in range(0,len(genes_start_stop)-1):
    if (genes_start_stop.at[i,'fitness']) == []:
        genes_start_stop.at[i,'fitness_avr'] =  'no insertions'
    else:
        genes_start_stop.at[i,'fitness_avr'] = mean(genes_start_stop.at[i,'fitness'])
        genes_start_stop.at[i,'fitness_variance'] = np.var(genes_start_stop.at[i,'fitness'] )

SMFitness = pd.read_excel(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\strain_ids_and_single_mutant_fitness.xlsx")

for i in range(0, len(genes_start_stop)):
    print(i)
    for j in range(0, len(SMFitness)):
        if (genes_start_stop.at[i,'fitness']) == []:
            genes_start_stop.at[i,'fitness_diff'] =  'no insertions'
        else:
            genes_start_stop.at[i,'gene'] == SMFitness.at[j,'Systematic gene name']
            genes_start_stop.at[i,'fitness_diff'] = SMFitness.at[j,'Single mutant fitness (26°)']/genes_start_stop.at[i,'fitness_avr']
            
            
plt.figure()
plt.plot(genes_start_stop.at[i,'fitness'])   
plt.plot( SMFitness.at[j,'Single mutant fitness (26°)'])     

