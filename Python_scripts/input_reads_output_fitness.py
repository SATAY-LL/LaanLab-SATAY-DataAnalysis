# -*- coding: utf-8 -*-
"""
Created on Thu May  6 16:15:48 2021

@author: floor
"""
import numpy as np
import random
import math 
import pandas as pd

N=100000000 #number of transposons before sequencing
m = 50 #different transposons

#number of reads for every transposon (ni). 
#Is now created but would be a list of transposons + #reads

random.seed = 0
ni = random.sample(range(0, 501), m)

#use imput to find n and m
n = sum(ni) #total number of reads#m = len(ni) #number of different type transposons

Ni = np.zeros(m)
Ci = np.zeros(m)
gf = np.zeros(m)
fitness = np.zeros(m)

p = 10  #number of PCR cycles
a = 2**-p


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
    fitness [i] = 1- (gf[i]-mn)/A #define fitness

