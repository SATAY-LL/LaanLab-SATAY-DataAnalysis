# -*- coding: utf-8 -*-
"""
Created on Tue May 18 15:40:27 2021

@author: floor
"""
import pandas as pd
from statistics import mean
import numpy as np
import csv
import matplotlib.pyplot as plt

#STEP III:
#compare our fitness values to cellmap data.


#get data from step II and put it in a dataframe
dataframe = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
dataframe['fitness'] = ''
dataframe['fitness_avr'] = ''
dataframe['fitness_variance'] = ''
dataframe['fitness_diff'] = ''
dataframe['fitness_singlemutant'] = ''

dataframe['fitness']= dataframe['fitness'].astype('object')

with open('data_step2.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data = list(reader)
dataframe['fitness'] = data

for i in range(0,len(dataframe)):
    if (dataframe.at[i,'fitness']) == []:
        dataframe.at[i,'fitness_avr'] =  'no insertions'
    else:
        dataframe.at[i,'fitness_avr'] = mean(dataframe.at[i,'fitness'])
        dataframe.at[i,'fitness_variance'] = np.var(dataframe.at[i,'fitness'] )

#data from cellmap:
SMFitness = pd.read_excel(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\strain_ids_and_single_mutant_fitness.xlsx")

# Genes from my data and cellmap with same name are put in the same row in dataframe.
for i in range(0, len(dataframe)):
    print(i)
    for j in range(0, len(SMFitness)):
        if (dataframe.at[i,'fitness']) == []:
            dataframe.at[i,'fitness_diff'] =  'no insertions'
        elif dataframe.at[i,'gene'] == SMFitness.at[j,'Systematic gene name']:
            dataframe.at[i,'fitness_diff'] = (SMFitness.at[j,'Single mutant fitness (26°)'])-(dataframe.at[i,'fitness_avr'])
            dataframe.at[i,'fitness_singlemutant'] = SMFitness.at[j,'Single mutant fitness (26°)']

#dataSMALL =  dataframe[ 
#      (len(dataframe['fitness'])<6) 
#      ]
#print(dataSMALL)    
#smallVAR = dataframe['fitness_variance'].mean()  


#making plots >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>           

singlemutantlist = []
myfitnesslist = [] 
chromosome  = 'chrVII'
for i in range(0,len(dataframe)):
    if type(dataframe.at[i,'fitness_singlemutant']) == np.float64 and str(dataframe.at[i,'fitness_singlemutant'])!= 'nan'  and type(dataframe.at[i,'fitness_avr']) == float and len(dataframe.at[i,'fitness']) >10  and (dataframe.at[i, 'chromosome'] ==chromosome):
        singlemutantlist.append(dataframe.at[i,'fitness_singlemutant'])
        myfitnesslist.append(dataframe.at[i,'fitness_avr'])

differencelist = []
for i in range(0,len(dataframe)):
    if type(dataframe.at[i,'fitness_singlemutant']) == np.float64 and str(dataframe.at[i,'fitness_singlemutant'])!= 'nan' :
        differencelist.append(dataframe.at[i,'fitness_singlemutant']-dataframe.at[i,'fitness_avr']) 
variancelist = []
for i in range(0,len(dataframe)):
    if type(dataframe.at[i,'fitness_variance']) == np.float64:
        variancelist.append(dataframe.at[i,'fitness_variance'])

plt.style.use('seaborn-whitegrid')
plt.plot(myfitnesslist, singlemutantlist,  'x', color='black')
print(mean((differencelist)))
print(np.var(differencelist))
plt.xlabel('calculated fitness')
plt.ylabel('cellmap fitness')  
plt.title('calculated fitness compared to cellmap fitness  '+str(chromosome))
plt.legend()
plt.savefig('afb.jpg', dpi=1200)
plt.show()






