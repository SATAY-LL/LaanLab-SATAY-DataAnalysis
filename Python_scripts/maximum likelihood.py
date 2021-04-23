# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 15:03:58 2021

@author: floor
"""



import matplotlib.pyplot as plt
from random import sample
from random import seed
import numpy as np
import pandas as pd


## create initial population
population = []
m=10 #total initial population size

for i in range (1, m, 2):
    population.append(1)
for i in range (1,m,2):
    population.append(2)

    
y = population.count(1) #total amount of tn1 in population

#use maximum likelihood estimator

error = []
m_over_n = []
reads1 = {} #amount of reads of type 1 present in the total population of reads
variance = []

for n in range (m//2,m):
        seed(37)
        reads = sample((population),n)
        a = m-n
        
        for i in range (1,2):
            matrix = np.zeros(shape= (2,a))
            matrix[i] = np.array(population.count(i))
            
            for j in range(1,a):
                matrix[i,j]= 1+ n/j
                    


