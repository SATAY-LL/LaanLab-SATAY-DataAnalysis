
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
population = []
real_Ni= 4 #how many of type 1 in the original population
N=10000 #total initial population size
notypes = 50 #number of types in initial population
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
    random.seed(i*13)
    samplelist = sample((population),N//2)
    for n in range (N//3, N//2):
            reads = samplelist[0:n]
            reads1= reads.count(1)
            if i == 0:
                nn.append(n)
            if reads1 == 0:
                print('The transposon we want to estimate has 0 reads...')
            
    #use reads to estimate transposons in original =hyper:
            estmateH= reads1*(N/n)  
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
   # the maximum likelihood estimation
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
                        
### NIEUWE CODE ###
            resultaat = [] #Lege lijst
            cols = len(j) + 1 # +1 omdat we straks een kolom met indexen gaan toevoegen
            matrix = np.zeros((rows,cols),dtype=np.float)
            #Eerst gaan we een kolom toevoegen met de originele indexen
            for r in range(0,rows):
                c = 0
                z = r
                position = cols * r + c
                np.put(matrix, position, z)
                #Daarnaast berekenen we ook alle waardes in de eerste kolom
                z = 1 + ni[r] / (c+1)
                position = cols * r + (c+1)
                np.put(matrix,position, z)
            #En we sorteren op basis van de kolom met de eerste waardes
            matrix = matrix[matrix[:, 1].argsort()]
            #Nu gaan we de rest invullen
            r = (rows - 1) #startwaarde van r
            c = 1          #startwaarde van c
            while r >= 0: #stoppen bij r = -1
                while c < cols: #stoppen als c groter is dan aantal kolommen
                    z = 1 + ni[int(matrix[r,0])] / c #bereken z
                    if len(resultaat) < a: #als de lengte van de resultaten nog kleiner zijn
                        resultaat.append(z) #voeg z toe aan de resultaten
                        position = cols * r + (c) #bereken de positie
                        np.put(matrix, position, z) #stop z in de matrix
                        c = c + 1 # ga naar de volgende kolom
                    else: #als er wel al meer dan N-n cellen berekend zijn
                        if z > min(resultaat): #als z groter is dan de kleinste waarde
                            resultaat[resultaat.index(min(resultaat))] = z #vervang de kleinste waarde
                            position = cols * r + (c) #bereken de positie
                            np.put(matrix, position, z) #stop z in de matrix
                            c = c + 1 # ga naar de volgende kolom
                        else: #als z kleiner is dan de kleinste waarde
                            break #ga dan naar de volgende rij
                r = r - 1 # ga naar de volgende rij
                c = 1 #begin weer bij c = 1
            matrix = matrix[matrix[:, 0].argsort()] #weer sorteren op de index
            matrix = np.delete(matrix,0,1) # en verwijderen van de index
    
        # #find the N-n largest values of the matrix
            num_largest = a
            indices = matrix.argpartition(matrix.size - num_largest, axis=None)[-num_largest:]
            g , h  = np.unravel_index(indices, matrix.shape)
            y_rest = pd.Series(collections.Counter(g)).to_frame()
        #make an estimate with Yi
            estimateM = ni[0]+y_rest[0][0]
            errmax=abs((estimateM- real_Ni)/real_Ni*100)
            error_max.append(errmax)
        #     a = N -real_Ni
        #     b = real_Ni/n
        #     c = N-n
        #     d = N-1
        #     variance_max.append(a*b*c/d)
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
    difference.append(abs(list1_i-list2_i))
maximum_diff = max(difference)  
print('biggest difference (error%) between max.likeli & methods of moments = ', maximum_diff)
        
#del(ah , bh , ch , dh, a,b, c, d)

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

#plt.figure()
#plt.plot(difference)
#plt.show

print("--- %s seconds ---" % (time.time() - start_time))
