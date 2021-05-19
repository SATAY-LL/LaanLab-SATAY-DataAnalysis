# -*- coding: utf-8 -*-
"""
Created on Thu May  6 16:15:48 2021

@author: floor
"""
import numpy as np
#import random
import math 
import pandas as pd
import ast
import wiggelen
import os, sys, tarfile
import matplotlib.pyplot as plt

def fitness (filepath_and_name):
    
    list = []
    for x in wiggelen.walk(open(filepath_and_name)):
        y = [x[0], x[1], x[2]]
        list.append(y)
    data = pd.DataFrame(list,columns=['chromosome','tn start position', '#reads'])
    
    #HO locus is neutral, compare to GF of this locus
    dataHO =  data[ 
      (data['chromosome']  == 'chrIV') &
      (46270 < data['tn start position']) & (data['tn start position'] <  48031) 
      ]
    readsHO = dataHO['#reads'].mean()   
    
    #remove the data of ADE2
    data = data[
        (566191<data['tn start position']) | (data['tn start position']<564476)]
    
    N=117601914802*10 #number of transposons before sequencing (())
    
        
    #use imput to find n and m
    n =data['#reads'].sum()      #total number of reads#m = len(ni) #number of different type transposons
    m =len(data) #number of different transposons
    
    ni = data['#reads'].to_numpy()
    Ni = np.zeros(m)
    Ci = np.zeros(m)
    gf = np.zeros(m)
    fitness = np.zeros(m)
    
    pcr = 15  #number of PCR cycles
    a = 2**-pcr
    
    
    for i in range (0,m):
        Ni[i] = ni[i] *N/n #calculate #tn before sequencing
        Ci[i] = Ni[i] * a  # calculate #tn before PCR
        
        if Ci[i] == 0:
            gf[i] = 0
        else: 
            gf[i] = math.log(Ci[i],(2)) #calculate the growth factor 
    CiHO = readsHO*N/m*a 
    gfHO = math.log(CiHO, (2))
    
    for i in range (0,m):
        fitness [i] = 1-(gf[i]/gfHO) #define fitness relative to HO locus
    
    data = data.join(pd.DataFrame(fitness))
    data.rename(columns = {0:'fitness'}, inplace = True)
    data = data.reset_index()
    return(data)


data = fitness(r'C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\Processed_data_WT1.wig')

plt.plot(data['fitness'],  linewidth=0.50)
plt.title('Fitness relative to the HO locus')
plt.xlabel('start position transposon')
plt.ylabel('relative fitness')

#plt.savefig('plots/afb', dpi=800)
plt.show

    
#del (a,i, list, x, y)
