# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 15:57:01 2021

@author: floor
"""

import pandas as pd
from statistics import mean
import numpy as np
import csv
import matplotlib.pyplot as plt
import pylab

#STEP III:
#compare our fitness values to cellmap data
#get data from step II and put it in a dataframe
dataframe = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
dataframe['fitnessWT1'] = ''
dataframe['fitness_avr1'] = ''
dataframe['fitness_variance1'] = ''
dataframe['fitness_diff'] = ''
dataframe['fitnessWT2'] = ''
dataframe['fitness_avr2'] = ''


dataframe['fitnessWT1']= dataframe['fitnessWT1'].astype('object')

with open('data_step2WT1.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data = list(reader)
dataframe['fitnessWT1'] = data

dataframe['fitnessWT2']= dataframe['fitnessWT2'].astype('object')
with open('data_step2WT2.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data2 = list(reader)
dataframe['fitnessWT2'] = data2


for i in range(0,len(dataframe)):
    if (dataframe.at[i,'fitnessWT1']) == []:
        dataframe.at[i,'fitness_avr1'] =  'no insertions'
    else:
        dataframe.at[i,'fitness_avr1'] = mean(dataframe.at[i,'fitnessWT1'])
        dataframe.at[i,'fitness_variance1'] = np.var(dataframe.at[i,'fitnessWT1'] )
    if (dataframe.at[i,'fitnessWT2']) == []:
        dataframe.at[i,'fitness_avr2'] =  'no insertions'
    else:
        dataframe.at[i,'fitness_avr2'] = mean(dataframe.at[i,'fitnessWT2'])
        dataframe.at[i,'fitness_variance2'] = np.var(dataframe.at[i,'fitnessWT2'])
 
 ## PLOTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
WT1 = []
WT2 = [] 
chromosome  = 'chrX'
for i in range(0,len(dataframe)):
    if  type(dataframe.at[i,'fitness_avr1']) == float and type(dataframe.at[i,'fitness_avr2']) == float and len(dataframe.at[i,'fitnessWT1']) >5 and len(dataframe.at[i,'fitnessWT2']) >5  :
        WT1.append(dataframe.at[i,'fitness_avr1'])
        WT2.append((dataframe.at[i,'fitness_avr2']))
        
plt.style.use('Solarize_Light2')
plt.plot(WT1, WT2,  'x', color='black')
plt.xlabel('WT1')
plt.ylabel('WT2')  
plt.title('Fitness WT1 and WT2' )#( +str(chromosome))

# calculate the trendline
z = np.polyfit(WT2, WT1, 1)
p = np.poly1d(z)
pylab.plot(WT1,p(WT1),"r--")
print ("y=%.6fx+(%.6f)"%(z[0],z[1]))

plt.legend()
plt.savefig('WT1vsWT2.jpg', dpi=1200)
plt.show()
        
