# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 16:14:01 2021

@author: floor
"""

import numpy as np
import pandas as pd
from ast import literal_eval
import matplotlib.pyplot as plt
import scipy.stats as stats
from collections import Counter
from reliability.Fitters import Fit_Everything

# Import curve fitting package from scipy
from scipy.optimize import curve_fit

# Import pergene_insertions.txt 
df = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\test data\ERR1533147_trimmed.sorted.bam_pergene_insertions.txt", sep = "\t")
df.columns = ["gene","chromosome", "start bp", "stop bp", "insertions", "reads"]

# Convert entire column to a list
df.loc[:,'insertions'] = df.loc[:,'insertions'].apply(lambda x: literal_eval(x))
df.loc[:,'reads'] = df.loc[:,'reads'].apply(lambda x: literal_eval(x))

# compute the average reads per insertions over all the data, this is the expected value (mu).
# so we can compare this number with all the genes

all_reads = [ ]
gene_list = [ ]
for index, row in df.iterrows():
    if row['insertions'] != [] and row['gene'] != 'ADE2':
        for x in row['reads']:
            all_reads.append(x)
            gene_list.append(row['gene'])

# There are really high values for the reads, which mess up the variance.
# With a max of 999 reads, you will still have a variance from around a few thousand.
# all_reads = [x for x in all_reads if x < 1000]
                                
totalins = len(all_reads)                  
totalreads = sum(all_reads)
print(max(all_reads))

#define mu and sigma for the normal distribution
i_min_mu_srd = []
mu = totalreads/(totalins)

for index in all_reads:
    i_min_mu_srd.append((index - mu)**2)  
    

sigma_sqrd = sum(i_min_mu_srd)/totalins
del (row, index)

#make a dataframe with genes and reads sorted so we can see to which genes the outliers belong
dictionary = {'gene':gene_list, 'reads':all_reads}
df_gene_reads = pd.DataFrame(dictionary)


# now we make a normal distribution from the expected value and sigma 
#sigma_sqrd = 20
# x = np.linspace(mu - 3*sigma_sqrd, mu + 3*sigma_sqrd, 100)
# plt.plot(x, stats.norm.pdf(x, mu, sigma_sqrd))
# plt.xlabel("#reads")
# plt.ylabel("abundance of readcount")
# plt.show()



# Now we want to see if the real distribution of read counts matches the expected distribution
abundance_reads_counter = Counter(all_reads)

abundance_reads = pd.DataFrame.from_dict(abundance_reads_counter, orient='index').reset_index()
abundance_reads.columns = ["number of reads","abundance"]
ab =  abundance_reads.sort_values(by ='number of reads' )
ab.plot.scatter(x= 'number of reads', y='abundance', s=2 ,title='Abundance of different reads counts' )
plt.xlim([0,400])
plt.ylim([0,500])
#plt.locator_params(axis='x', nbins=10) 
plt.xticks(np.arange(0, 500, 100))
#plt.xticks(np.arange(0, 10000, 200))

Fit_Everything(ab)

# labels, values = (abundance_reads.items())
# indexes = np.arange(len(labels))
# width = 1
# plt.bar(indexes, values, width)
# #plt.xticks(indexes + width * 0.5)
# plt.xlabel("#reads")
# plt.ylabel("abundance of readcount")
# plt.show()


# for index in abundance_reads:
#     xvalues = abundance_reads['Key']
#     yvalues = abundance_reads['Value']
# plt.bar(xvalues, yvalues, width=0.8, bottom=None, align='center')






     