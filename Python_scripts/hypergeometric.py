# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 10:54:54 2021

@author: floor
"""

import random
#import matplotlib.pyplot as plt
from random import sample
import numpy as np


#creating initial population with i types of transposons
i=70000
types = list(range(i))

m=10000000 #= total number of transposons in initial population
initial_population = random.choices(types, k=m)
DNA1_initial= initial_population.count(1)
perc_DNA1 = DNA1_initial/m
#DNA2_initial= initial_population.count(2)

del(types)
# sampling/sequencing
n=1687598 #= total number of reads
reads = sample((initial_population),n)
reads_DNA1= reads.count(1)
#reads_DNA2= reads.count(2)


# With k and n we want to make an estimation of the percentage of DNA1 in the initial population
perc_DNA1_calculated = reads_DNA1/n
var = (perc_DNA1)*(1-perc_DNA1)
print('real portion of DNA1 in population =', perc_DNA1, 'calculated portion=', perc_DNA1_calculated, 'variance=', var)