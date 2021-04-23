<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 14:31:20 2021

@author: floor
"""


import matplotlib.pyplot as plt
from random import sample
from random import seed
import numpy as np


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
      
N_real = population.count(1) #total amount of tn1 in population

error_hyper = []
m_over_n = []
reads1 = {} #amount of reads of type 1 present in the total population of reads
variance_hyper = []
for n in range (50,N-N//10):
        seed(8)
        reads = sample((population),n)
        reads1= reads.count(1)

        

#use reads to estimate transposons in original =hyper:
        r= reads1*(N/n)  
        error_hyper.append(abs(r- N_real)/N_real*100)       
        
        m_over_n.append(N/n)
        ah = N -N_real
        bh = N_real/n
        ch = N-n
        dh = N-1
        variance_hyper.append(ah*bh*ch/dh)
del(ah , bh , ch , dh)


#plots:
# plot for the total number of reads N vs. error
plt.figure()
plt.title('error of hypergeometric')
plt.plot(error_hyper, linewidth=.5)
plt.xlabel('total number of reads drawn (n)')
plt.ylabel('error (= |r-y|) as pecentage of m (initial tot. pop.)')  
#plt.savefig('plots/hypergeometric_M=10000_N_vs_error', dpi=800)

# plot for error vs initial total population M vs error
plt.figure()
plt.plot(m_over_n, error_hyper, linewidth=.5)
plt.xlabel('m (total initial population) / n (total reads drawn)')  
#plt.xticks(range(1,21))
plt.grid()
plt.savefig('plots/hypergeometric_M=10000_MoverN_vs_error', dpi=800)
plt.ylabel('error (= |r-y|) as pecentage of m (initial tot. pop.)')  


#plot of variance vs #reads
plt.figure()
plt.plot(variance_hyper)
plt.grid()
#plt.xticks(range(1,21))
plt.xlabel('m / n')
plt.ylabel('variance')
#plt.savefig('plots/variance vs reads', dpi=800)

plt.show



=======
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 14:31:20 2021

@author: floor
"""


import matplotlib.pyplot as plt
from random import sample
from random import seed
import numpy as np


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
      
N_real = population.count(1) #total amount of tn1 in population

error_hyper = []
m_over_n = []
reads1 = {} #amount of reads of type 1 present in the total population of reads
variance_hyper = []
for n in range (50,N-N//10):
        seed(8)
        reads = sample((population),n)
        reads1= reads.count(1)

        

#use reads to estimate transposons in original =hyper:
        r= reads1*(N/n)  
        error_hyper.append(abs(r- N_real)/N_real*100)       
        
        m_over_n.append(N/n)
        ah = N -N_real
        bh = N_real/n
        ch = N-n
        dh = N-1
        variance_hyper.append(ah*bh*ch/dh)
del(ah , bh , ch , dh)


#plots:
# plot for the total number of reads N vs. error
plt.figure()
plt.title('error of hypergeometric')
plt.plot(error_hyper, linewidth=.5)
plt.xlabel('total number of reads drawn (n)')
plt.ylabel('error (= |r-y|) as pecentage of m (initial tot. pop.)')  
#plt.savefig('plots/hypergeometric_M=10000_N_vs_error', dpi=800)

# plot for error vs initial total population M vs error
plt.figure()
plt.plot(m_over_n, error_hyper, linewidth=.5)
plt.xlabel('m (total initial population) / n (total reads drawn)')  
#plt.xticks(range(1,21))
plt.grid()
plt.savefig('plots/hypergeometric_M=10000_MoverN_vs_error', dpi=800)
plt.ylabel('error (= |r-y|) as pecentage of m (initial tot. pop.)')  


#plot of variance vs #reads
plt.figure()
plt.plot(variance_hyper)
plt.grid()
#plt.xticks(range(1,21))
plt.xlabel('m / n')
plt.ylabel('variance')
#plt.savefig('plots/variance vs reads', dpi=800)

plt.show



>>>>>>> 7a742ccd0a602eae31a515f74c1de5b4f2edf0f9
    