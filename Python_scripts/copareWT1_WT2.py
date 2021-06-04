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
import statistics
import operator
#STEP III:
#compare our fitness values to cellmap data
#get data from step II and put it in a dataframe
dataframe = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
dataframe['readsWT1'] = ''
dataframe['fitness_avr1'] = ''
dataframe['fitness_variance1'] = ''
dataframe['fitness_variance2'] = ''
dataframe['readsWT2'] = ''
dataframe['fitness_avr2'] = ''
dataframe['error1'] = ''
dataframe['error2'] = ''


dataframe['readsWT1']= dataframe['readsWT1'].astype('object')

with open('data_step2_readsWT1.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data = list(reader)
dataframe['readsWT1'] = data

dataframe['readsWT2']= dataframe['readsWT2'].astype('object')

with open('data_step2_readsWT2.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data2 = list(reader)
dataframe['readsWT2'] = data2


z = 1.960 #for 95 percent confidence interval
for i in range(0,len(dataframe)):


    if (dataframe.at[i,'readsWT1']) == [] or len(dataframe.at[i, 'readsWT1']) < 2:
        dataframe.at[i,'fitness_avr1'] =  ''
       
    else:
        dataframe.at[i,'fitness_avr1'] = statistics.mean(dataframe.at[i,'readsWT1'])
        dataframe.at[i,'fitness_variance1'] = np.var(dataframe.at[i,'readsWT1'] )
        dataframe.at[i, 'error1'] = 2*(z*(statistics.stdev(dataframe.at[i, 'readsWT1']))/(np.sqrt(len(dataframe.at[i,'readsWT1']))))
    
    if (dataframe.at[i,'readsWT2']) == [] or len(dataframe.at[i, 'readsWT2']) < 2:
        dataframe.at[i,'fitness_avr2'] =  ''
    else:
        dataframe.at[i,'fitness_avr2'] = statistics.mean(dataframe.at[i,'readsWT2'])
        dataframe.at[i,'fitness_variance2'] = np.var(dataframe.at[i,'readsWT2'])
        dataframe.at[i, 'error2'] = 2*(z*(statistics.stdev(dataframe.at[i, 'readsWT2']))/(np.sqrt(len(dataframe.at[i,'readsWT2']))))

 ## PLOTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
WT1 = []
WT2 = [] 
err1 = []
err2 = []
chromosome  = 'chrXII'
for i in range(0,len(dataframe)):
    if  type(dataframe.at[i,'fitness_avr1']) == float and type(dataframe.at[i,'fitness_avr2']) == float and len(dataframe.at[i,'readsWT1']) >5 and len(dataframe.at[i,'readsWT2']) >5 and (dataframe.at[i,'fitness_avr1']) <500 and (dataframe.at[i,'fitness_avr2']) <500 :
        WT1.append(dataframe.at[i,'fitness_avr1'])
        WT2.append((dataframe.at[i,'fitness_avr2']))
        err1.append((dataframe.at[i,'error1']))
        err2.append((dataframe.at[i,'error2']))
        
plt.style.use('Solarize_Light2')
plt.plot(WT1, WT2,  'x', color='black')


#plt.errorbar(WT1, WT2,  xerr=err1, yerr=err2)
plt.xlabel('WT1')
plt.ylabel('WT2')  
plt.title('Fitness WT1 and WT2' )#( +str(chromosome))

# calculate the trendline
z = np.polyfit(WT2, WT1, 1)
p = np.poly1d(z)
pylab.plot(WT1,p(WT1),"r--")
print ("y=%.6fx+(%.6f)"%(z[0],z[1]))

plt.legend()
plt.savefig('afb.jpg', dpi=1200)
plt.show()
        
