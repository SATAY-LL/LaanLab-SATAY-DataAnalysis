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
DNA1= initial_population.count(1)
perc_DNA1 = DNA1/m
DNA2= initial_population.count(2)
perc_DNA2 = DNA2/m

del(types)
# sampling/sequencing
n=1687598 #= total number of reads
reads = sample((initial_population),n)
reads_DNA1= reads.count(1)
reads_DNA2= reads.count(2)


# With k and n we want to make an estimation of the percentage of DNA1 in the initial population
DNA1_calculated = reads_DNA1*n/m
#var = n*DNA1_initial/m*(m-DNA1_initial)/m
print('real portion of DNA1 in population =', DNA1, 'calculated =', DNA1_calculated, 'variance=', 'var')

DNA2_calculated = reads_DNA2*n/m
#var2 = n*DNA2/m*(m-DNA2_initial)/m
print('real amount of DNA2 in population =', DNA2, 'calculated  DNA2=', DNA2_calculated, 'variance=', 'var2')