# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 14:31:20 2021

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
types = 4

b=1
for i in range (1, N+1):
    population.append(b)
    b = b + 1
    if b > types:
        b = 1
      
real_Ni = population.count(1) #total amount of tn1 in population

error_hyper = []
error_max =[]
m_over_n = []
reads1 = {} #amount of reads of type 1 present in the total population of reads
variance_hyper = []
variance_max=[]


#here, we create the samples. the range of total sample size can be varied here.
for n in range (300, 650):
        seed(5)
        reads = sample((population),n)
        reads1= reads.count(1)

#use reads to estimate transposons in original =hyper:
        estmateH= reads1*(N/n)  
        errhyp = abs(estmateH- real_Ni)/real_Ni*100
        error_hyper.append(errhyp) # errror hypergeometric     
        
        m_over_n.append(N/n)
#compute variance
        ah = N -real_Ni
        bh = real_Ni/n
        ch = N-n
        dh = N-1
        variance_hyper.append(ah*bh*ch/dh)
        a = N-n
        ni = []
#the maximum likelihood estimation
    #first, i make a list of all types and the number of reads.
        for count in range (1,types+1):
            ni.append(reads.count(count))

        j = np.arange(1,a+1).tolist()
        cols = len(j)
        rows = types
    #the matrix to compute the likelihood
        matrix = np.zeros((rows,cols),dtype=np.float)
        x = 0
        y = 1

        for r in range(0,rows):
            for c in range(0,cols):
                z = 1 + ni[r] / (c+1)
                position = cols * r + c       
                np.put(matrix,position, z)
    #find the N-n largest values of the matrix
        num_largest = a
        indices = matrix.argpartition(matrix.size - num_largest, axis=None)[-num_largest:]
        g , h  = np.unravel_index(indices, matrix.shape)
        y_rest = pd.Series(collections.Counter(g)).to_frame()
    #make an estimate with Yi
        estimateM = ni[0]+y_rest[0][0]
        errmax=abs((estimateM- real_Ni)/real_Ni*100)
        error_max.append(errmax)       #error max likeli
        a = N -real_Ni
        b = real_Ni/n
        c = N-n
        d = N-1
        variance_max.append(a*b*c/d)
        print(estimateM)
        print(estmateH)
        del(ni)

difference = []

zip_object = zip(error_hyper, error_max)
for list1_i, list2_i in zip_object:
    difference.append(list1_i-list2_i)
maximum_diff = max(difference)  
print('biggest difference (error%) between max.likeli & hypergeo. = ', maximum_diff)
    
del(ah , bh , ch , dh, a,b, c, d)


#plots:
# plot for the total number of reads N vs. error
plt.figure()
plt.title('error of hypergeometric and maximum likelihood, types = '+str(types))
plt.plot(error_hyper, linewidth=.5, label= 'hypergeometrix')
plt.plot(error_max, linewidth=.5, color = 'red', label = 'max. likelihood')
plt.legend()
plt.xlabel('total number of reads drawn (n)')
plt.ylabel('error (= |r-y|) as pecentage of m (initial tot. pop.)')  
#plt.savefig('plots/hypergeometric_M=10000_N_vs_error', dpi=1500)

# plot for error vs initial total population M vs error
plt.figure()
plt.plot(m_over_n, error_hyper, linewidth=.5)
plt.xlabel('m (total initial population) / n (total reads drawn), types = '+str(types))
#plt.xticks(range(1,21))
plt.grid()
#plt.savefig('plots/hypergeometric_M=10000_MoverN_vs_error', dpi=800)
plt.ylabel('error (= |r-y|) as pecentage of m (initial tot. pop.)')  


#plot of variance vs #reads
plt.figure()
plt.title('variance of hyper and max. likeli, , types = '+str(types))
plt.plot(variance_hyper, label = 'hypergeometric')
plt.plot(variance_max, label = 'maximum likelihood')
plt.grid()
#plt.xticks(range(1,21))
plt.xlabel('m / n')
plt.ylabel('variance')
#plt.savefig('plots/variance vs reads', dpi=800)

plt.figure()
plt.plot(difference)
plt.show



    