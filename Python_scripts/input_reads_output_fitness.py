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


N=117601914802*10 #number of transposons before sequencing
 #different transposons

#number of reads for every transposon (ni). 
# Import .txt 
a = ['variableStep',0]
df = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\dataWT1.txt", sep = " ")
df.reset_index(level=0, inplace=True)
df = df[~df['index'].isin(a)]
df = df[~df['1'].isin(a)]
ins_loc = df['1'].tolist()
ni = df['2'].tolist()

for i in range(0,len(ni)):
    ni[i] = int(ni[i])
for i in range(0,len(ni)):
    ins_loc[i] = int(ins_loc[i])
    
#use imput to find n and m
n = sum(ni) #total number of reads#m = len(ni) #number of different type transposons
m = len(ins_loc) #number of different transposons

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

ins_fitness = np.column_stack((ins_loc, fitness))
    
del (A, mx, mn, a,i)
