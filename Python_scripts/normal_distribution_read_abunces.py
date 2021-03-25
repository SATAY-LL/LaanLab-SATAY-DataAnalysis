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
from scipy.stats import nbinom, poisson, geom

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
    if row['insertions'] != [] and row['gene'] != 'ADE2': #remove ADE2
        for x in row['reads']:
            all_reads.append(x)
            gene_list.append(row['gene'])

# There are really high values for the reads, which mess up the variance.
# With a max of 999 reads, you will still have a variance from around a few thousand.
# all_reads = [x for x in all_reads if x < 1000]
                                
totalins = len(all_reads)                  
totalreads = sum(all_reads)

#make a dataframe with genes and reads sorted so we can see to which genes the outliers belong
dictionary = {'gene':gene_list, 'reads':all_reads}
df_gene_reads = pd.DataFrame(dictionary)

# plot the distribution of read counts :
abundance_reads_counter = Counter(all_reads)

abundance_reads = pd.DataFrame.from_dict(abundance_reads_counter, orient='index').reset_index()
abundance_reads.columns = ["number of reads","abundance"]
ab =  abundance_reads.sort_values(by ='number of reads' )

ab.plot.scatter(x= 'number of reads', y='abundance', s=2 ,title='Abundance of different reads counts' )
plt.xlim([0,400])
plt.ylim([0,2000])
#plt.xticks(np.arange(0, 8500, 500))

# Now we are trying to fit our data to an exponential distribution

# Function to calculate the exponential with constants a and b
def exponential(x, a, b):
    return a*np.exp(b*x)

x_data = ab['number of reads'].reset_index(drop=True).truncate(0, 400, copy = False)
y_data = ab['abundance'].reset_index(drop=True).truncate(0, 2000, copy = False);

pars, cov = curve_fit(f=exponential, xdata= x_data, ydata=y_data, p0=[0, 0], bounds=(-np.inf, np.inf))

# Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
stdevs = np.sqrt(np.diag(cov))

# Calculate the residuals
res = y_data - exponential(x_data, *pars)

# Plot the fit data as an overlay on the scatter data
ax = plt.axes()
ax.plot(x_data, exponential(x_data, *pars), linestyle='--', linewidth=2, color='black')
plt.xlim([0,400])
plt.ylim([0,2000])

# plt.xscale('log')
# plt.yscale('log')

#plt.savefig('plots/fit data to exponential funtion- #reads_ 1to421-zoomed in.png', dpi=800)

del(x)

# # geometric distribution
y_data = y_data/max(ab['abundance'])

def geometric(x, p):
    return (1-p)**(x)*p

pars, cov1 = curve_fit(f=geometric, xdata= x_data, ydata=y_data, p0=[0], bounds=(-np.inf, np.inf))

# Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
stdevs1 = np.sqrt(np.diag(cov1))

# Calculate the residuals
res = y_data - geometric(x_data, *pars)

plt.figure()
plt.scatter(x_data, y_data, s=2)
plt.title('Distribution of #reads, scaled to max(abundance) ')
plt.xlabel('number of reads')
plt.ylabel('Percentage of transposons with 1 read')

ax = plt.axes()
ax.plot(x_data, geometric(x_data, *pars), linestyle='--', linewidth=2, color='black')
plt.xlim([0,400])
plt.ylim([0,0.1])
#plt.legend((line1, line2), ('label1', 'label2'))
# plt.xscale('log')
#plt.yscale('log')

#plt.savefig('plots/fit data to geometric distribution.png', dpi=800)

 


     