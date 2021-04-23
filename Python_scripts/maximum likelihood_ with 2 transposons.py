# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 10:52:22 2021

@author: floor
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from random import sample
from random import seed
import collections, numpy

## create initial population
population = []
N=1000 #total initial population size
types = 3
b=1
for i in range (1, N+1):
    population.append(b)
    b = b + 1
    if b > types:
        b = 1
    
real_Ni = population.count(1) #total amount of tn1 in population

#use maximum likelihood estimator

error = []
m_over_n = []
reads1 = {} #amount of reads of type 1 present in the total population of reads
variance = []

for n in range (50,N-N//10):
        seed(8)
        reads = sample((population),n)
        a = N-n

        ni = [reads.count(1),reads.count(2), reads.count(3)]

        j = np.arange(1,a+1).tolist()
        cols = len(j)
        rows = len(ni)   

        matrix = np.zeros((rows,cols),dtype=np.float)
        x = 0
        y = 1

        for r in range(0,rows):
            for c in range(0,cols):
                z = 1 + ni[r] / (c+1)
                position = cols * r + c       
                np.put(matrix,position, z)
        num_largest = a
        indices = matrix.argpartition(matrix.size - num_largest, axis=None)[-num_largest:]
        g , h  = np.unravel_index(indices, matrix.shape)
        y_rest = pd.Series(collections.Counter(g)).to_frame()
        estimate = ni[0]+y_rest[0][0]
        error.append(abs(estimate- real_Ni)/real_Ni*100)       
        m_over_n.append(N/n)
        a = N -y
        b = y/n
        c = N-n
        d = N-1
        variance.append(a*b*c/d)
del(a , b , c , d)
        
#plots:
# plot for the total number of reads N vs. error
plt.figure()
plt.title('Error of max. likelihood')
plt.plot(error, linewidth=.5)
plt.xlabel('total number of reads drawn (n)')
plt.ylabel('error (= |r-y|) as pecentage of m (initial tot. pop.)')  
#plt.savefig('plots/hypergeometric_M=10000_N_vs_error', dpi=800)

# plot for error vs initial total population M vs error
plt.figure()
plt.plot(m_over_n, error, linewidth=.5)
plt.xlabel('m (total initial population) / n (total reads drawn)')  
#plt.xticks(range(1,21))
plt.grid()
plt.savefig('plots/hypergeometric_M=10000_MoverN_vs_error', dpi=800)
plt.ylabel('error (= |r-y|) as pecentage of m (initial tot. pop.)')  


#plot of variance vs #reads
plt.figure()
plt.plot(variance)
plt.grid()
#plt.xticks(range(1,21))
plt.xlabel('m / n')
plt.ylabel('variance')
#plt.savefig('plots/variance vs reads', dpi=800)

plt.show

       

    
    