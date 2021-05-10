
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
import time
start_time = time.time()

## create initial population
real_Ni = 4 #total amount of type 1 in initial population
N=1000000 #total initial population size
notypes = 10000 #number of different types in initial population
types = []
notimes = 1 #how many times to repeat the model to average over

population = list(np.ones(real_Ni))
b=2
for i in range (1, notypes+1):  #define the different types from 2-notypes
    types.append(b)
    b = b + 1
    
random.seed(87)
sampling = random.choices(types, k=(N-real_Ni)) #sample all types except type 1
population = population + sampling #add samples of type 2-notypes to type 1

i = 0
dfh = pd.DataFrame()
error_hyper = []
error_max = []

m_over_n = []
reads1 = {} #amount of reads of type 1 present in the total population of reads
variance_hyper = []
variance_max=[]
nn= []
#aa=[]
#esti=[]
#here, we create the samples. the range of total sample size can be varied here.
while i < notimes:
    error_max = []
    error_hyper = []

    random.seed(i*13)
    samplelist = sample((population),N//2)
    for n in range (N//2, N//2+3000):
            reads = samplelist[0:n]
            reads1= reads.count(1)
            #aa = aa + list(reads1)
            if i == 0:
                nn.append(n)
            if reads1 == 0:
                print('The transposon we want to estimate has 0 reads...')
                
                
    #use reads to estimate transposons in original =hyper:
            estmateH= reads1*(N/n) 
            #esti = esti.append(estmateH)
            errhyp = abs(estmateH- real_Ni)/real_Ni*100
            error_hyper.append(errhyp)
        
            if i == 0:
                m_over_n.append(N/n)
            print(n)
    # #compute variance
    #         ah = N -real_Ni
    #         bh = real_Ni/n
    #         ch = N-n
    #         dh = N-1
    #         variance_hyper.append(ah*bh*ch/dh)
            a = N-n
            ni = []
            del(ni)
    dfh.insert(i,i,error_hyper, True)
    i = i + 1

    
meanlisthyper = dfh.mean(axis = 1)
    
#del(ah , bh , ch , dh, a,b, c, d)

#plots:
# plot for the total number of reads N vs. error
plt.figure()
plt.title('error vs n (total reads drawn), types = '+str(notypes))
plt.plot(nn, meanlisthyper, linewidth=.5, label= 'meth. of mom.')
#plt.plot(nn, error_hyper, linewidth=.5, label= 'only one time', color = 'green')
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

# plt.figure()
# plt.plot(nn,aa, color = 'green', label = 'reads of type 1')
# plt.plot(nn,real_Ni, color = 'red', label = 'real Ni')
# plt.plot(nn,esti, color = 'blue', label = 'estimate of real Ni')


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

#plt.figure()
#plt.plot(difference)
#plt.show

print("--- %s seconds ---" % (time.time() - start_time))
