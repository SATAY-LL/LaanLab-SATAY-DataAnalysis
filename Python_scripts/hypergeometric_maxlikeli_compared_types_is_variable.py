
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 14:31:20 2021

@author: floor
"""


import numpy as np
import pandas as pd
#import statsmodels.api as sm
import matplotlib.pyplot as plt
import random
from random import seed
from random import sample
import collections, numpy

## create initial population
population = []
N=1000 #total initial population size
notypes = 5
types = []

b=1
for i in range (1, notypes+1):
    types.append(b)
    b = b + 1
# np.random.seed(4)
# randomlist = np.random.randint(1,30, size = notypes) 
randomlist = [1,20,30,10,1]
random.seed(13)
population = random.choices(types, weights=randomlist, k=N)
      
real_Ni = population.count(1) #total amount of tn1 in population
type2 = population.count(2)
type3 = population.count(3)
type4 = population.count(4)
type5 = population.count(5)

notimes = 5 #how many times to repeat the model to average over
i = 0
df = pd.DataFrame()
dfh = pd.DataFrame()
error_hyper = []
error_max = []

m_over_n = []
reads1 = {} #amount of reads of type 1 present in the total population of reads
variance_hyper = []
variance_max=[]
nn= []
#here, we create the samples. the range of total sample size can be varied here.
while i < notimes:
    error_max = []
    error_hyper = []
    for n in range (200, 900):
            random.seed(i*13)
            reads = sample((population),n)
            reads1= reads.count(1)
            if i == 0:
                nn.append(n)
            if reads1 == 0:
                print('Try another seed for the randomlist, because the transposon we want to estimate has 0 reads...')
            
            
    #use reads to estimate transposons in original =hyper:
            estmateH= reads1*(N/n)  
            errhyp = abs(estmateH- real_Ni)/real_Ni*100
            error_hyper.append(errhyp)
        
            if i == 0:
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
        #first, we make a list of all types and the number of reads.
            for count in range (1,notypes+1):
                ni.append(reads.count(count))
            
        #the matrix to compute the likelihood
            j = np.arange(1,a+1).tolist()
            cols = len(j)
            rows = notypes
    
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
            error_max.append(errmax)
            a = N -real_Ni
            b = real_Ni/n
            c = N-n
            d = N-1
            variance_max.append(a*b*c/d)
            #print(estimateM)
            #print(estmateH)
            del(ni)
    df.insert(i,i,error_max, True)  #difference 
    dfh.insert(i,i,error_hyper, True)
    i = i + 1

    
meanlist = df.mean(axis = 1)
meanlisthyper = dfh.mean(axis = 1)

difference = []
    
zip_object = zip(meanlist, meanlisthyper)
for list1_i, list2_i in zip_object:
    difference.append(list1_i-list2_i)
maximum_diff = max(difference)  
print('biggest difference (error%) between max.likeli & hypergeo. = ', maximum_diff)
        
del(ah , bh , ch , dh, a,b, c, d)

#plots:
# plot for the total number of reads N vs. error
plt.figure()
plt.title('error of method of moments and maximum likelihood, types = '+str(notypes))
plt.plot(nn, meanlisthyper, linewidth=.5, label= 'meth. of mom.')
#plt.plot(nn, error_hyper, linewidth=.5, label= 'only one time', color = 'green')
plt.plot(nn, meanlist, linewidth=.5, color = 'red', label = 'max. likelihood')
plt.legend()
plt.xlabel('total number of reads drawn (n)')
plt.ylabel('error (= |r-y|) as pecentage of m (initial tot. pop.)')  
plt.savefig('plots/afb', dpi=1500)

# # plot for error vs initial total population M vs error
# plt.figure()
# plt.plot(m_over_n, meanlisthyper, linewidth=.5,label= 'hypergeometric' )
# plt.plot(m_over_n, meanlist, linewidth =.5, label = 'maximum likelihood', color = 'red')
# plt.xlabel('m (total initial population) / n (total reads drawn), types = '+str(notypes))
# #plt.xticks(range(1,21))
# plt.grid()
# plt.legend()
# #plt.savefig('plots/hypergeometric_M=10000_MoverN_vs_error', dpi=800)
# plt.ylabel('error (= |r-y|) as pecentage of m (initial tot. pop.)')  


# #plot of variance vs #reads
# plt.figure()
# plt.title('variance of hyper and max. likeli, , types = '+str(notypes))
# plt.plot(variance_hyper, label = 'hypergeometric')
# plt.plot(variance_max, label = 'maximum likelihood')
# plt.grid()
# #plt.xticks(range(1,21))
# plt.xlabel('m / n')
# plt.ylabel('variance')
# #plt.savefig('plots/variance vs reads', dpi=800)

plt.figure()
plt.plot(difference)
plt.show
