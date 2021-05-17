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


 
    
list = []
for x in wiggelen.walk(open(r'C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\Processed_data_WT1.wig')):
    y = [x[0], x[1], x[2]]
    list.append(y)
data = pd.DataFrame(list,columns=['chromosome','tn start position', '#reads'])

N=117601914802*10 #number of transposons before sequencing (())
 #different transposons

#number of reads for every transposon (ni). 
# Import .txt 
    
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

mx = max (gf)
mn = min(gf)
A = (mx-mn)

for i in range (0,m):
    fitness [i] = 1- (gf[i]-mn)/A #define fitness ion scale 0-1

data = data.join(pd.DataFrame(fitness))

plt.plot(fitness)
plt.xlabel('start position transposon')
plt.ylabel('relative fitness')
#plt.savefig('plots/afb', dpi=800)
plt.show
    
del (A, mx, mn, a,i, list, x, y)
