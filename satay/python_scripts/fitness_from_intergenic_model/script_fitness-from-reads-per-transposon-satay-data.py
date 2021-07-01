# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 16:27:36 2021

@author: linigodelacruz
"""

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import os
import seaborn as sns
import scipy 

from src.python_modules.module_intergenic_model import getting_r,pd_to_compare_libraries,filtering_significant_genes,viz_data

#%% Reading from a post processed dataset exported to an excel sheet

names_libraries={'wt':'data_wt_filtered_reads_per_tr.xlsx','dnrp1':'data_nrp1_filtered_reads_per_tr.xlsx'}
data_library=[]
for i in names_libraries.keys():
 data_library.append(pd.read_excel('datasets/'+names_libraries[i],index_col='Unnamed: 0'))

datasets=data_library

del data_library,i
#%% Removing the ADE2 and URA3 gene reads and insertions from the population

for i in np.arange(0,len(datasets)):
    datasets[i]=datasets[i][datasets[i].Standard_name != 'ADE2']
    datasets[i]=datasets[i][datasets[i].Standard_name != 'URA3']

#%%  From reads to fitness from a intergenic competition model
## asumming we start with zero reads so one copy per cell
# r=K/T(K-N)
# K=total number of reads/total of transposons (carrying capacity) per library- excluding the ade2 and URA3 gene
# T=90
# N=number of reads per gene/transposons (excluding the ade2 gene)

## Getting the rates from the function : getting_r

rates,rates_no_filter=getting_r(datasets)


#%% Preparing datasets for analyses

data_library_pd=pd.concat(datasets,keys=names_libraries.keys(),sort=True)
data_library_pd.fillna(0,inplace=True)

data_wt=data_library_pd.loc['wt']
data_nrp1=data_library_pd.loc['dnrp1']

data_wt.loc[:,'rates-intergenic']=np.abs(rates[0].values)

data_nrp1.loc[:,'rates-intergenic']=np.abs(rates[1].values)

data_wt.loc[:,'rates-intergenic-non-filter']=np.abs(rates_no_filter[0].values)

data_nrp1.loc[:,'rates-intergenic-non-filter']=np.abs(rates_no_filter[1].values)

data_wt.replace([np.inf, -np.inf], np.nan, inplace=True)
data_wt.fillna(0,inplace=True)

data_nrp1.replace([np.inf, -np.inf], np.nan, inplace=True)
data_nrp1.fillna(0,inplace=True)
#%% exporting datasets

data_wt.to_excel('datasets/data_wt_rates.xlsx')
data_nrp1.to_excel('datasets/data_dnrp1_rates.xlsx')
#%% Fitness plots
data_wt.index=data_wt['Standard_name']
data_nrp1.index=data_nrp1['Standard_name']


## To discard the non coding genes from the plot
data_wt=data_wt[data_wt['Standard_name']!= 'noncoding']
data_nrp1=data_nrp1[data_nrp1['Standard_name']!= 'noncoding']

values2HO_wt=data_wt['rates-intergenic']/data_wt.loc['HO','rates-intergenic']
values2HO_dnrp1=data_nrp1['rates-intergenic']/data_wt.loc['HO','rates-intergenic']

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
ax.plot(values2HO_wt,values2HO_dnrp1,'bo',alpha=0.2)

plt.plot(np.linspace(0.5,3),np.linspace(0.5,3)*values2HO_wt.loc['NRP1'],linewidth=1,alpha=0.4,color='k')
plt.plot(np.linspace(0.5,3),np.linspace(0.5,3),linewidth=1,alpha=0.4,color='k')

# for i in range(0,len(values2HO_wt)):
    
#     if values2HO_dnrp1[i]<0.5*(values2HO_wt[i]*values2HO_wt.loc['NRP1']):
#         ax.text(values2HO_wt[i],values2HO_dnrp1[i],values2HO_dnrp1.index[i],fontsize=6)
#     elif values2HO_dnrp1[i]/values2HO_wt[i] > 3:
#         ax.text(values2HO_wt[i],values2HO_dnrp1[i],values2HO_dnrp1.index[i],fontsize=6)

        

#ax.set_xlim(0,5)
#ax.set_ylim(0,5)
ax.set_xlabel('single knockout fitness')
ax.set_ylabel('double knockout fitness: dnrp1dgenex')
#plt.savefig('dnrp1-fitness-map-from-satay.png',dpi=300,transparent=True,format='png')

        
#%% saving fitness plots

fig.savefig('dnrp1-fitness-map-from-satay-filtered.png',dpi=300,transparent=False,format='png')


#%% reads per transposons vs fitness values 

T=90
#K=np.sum(N)
K_wt=data_wt['reads-per-tr'].sum()
N_wt=data_wt['reads-per-tr']
intergenic_model=np.log(N_wt*K_wt/(K_wt-N_wt))/T
#intergenic_model=np.abs(intergenic_model)
ref_HO_wt=intergenic_model.loc['HO']



K_dnrp1=data_nrp1['reads-per-tr'].sum()
N_dnrp1=data_nrp1['reads-per-tr']
intergenic_model_nrp1=np.log(N_dnrp1*K_dnrp1/(K_dnrp1-N_dnrp1))/T
#intergenic_model_nrp1=np.abs(intergenic_model_nrp1)
ref_HO_dnrp1=intergenic_model_nrp1.loc['HO']


fig = plt.figure(figsize=(13,5))
ax = fig.add_subplot(121)
ax.scatter(data_wt['reads-per-tr'],intergenic_model/ref_HO_wt,alpha=0.5,color='b')
ax2 = fig.add_subplot(122)
ax2.scatter(data_nrp1['reads-per-tr'],intergenic_model_nrp1/ref_HO_wt,alpha=0.5,color='b')

ax.set_title('Merged WT')
ax2.set_title('Merged dnrp1')
for axes in [ax,ax2]:
    axes.set_ylabel('fitness values intergenic model')
    axes.set_xlabel('reads per transposon per gene')
    #axes.set_xlim(0,0.12)
    axes.set_ylim(0,2)
    axes.grid()

intergenic_model.replace([np.inf, -np.inf], np.nan, inplace=True)
intergenic_model.fillna(0,inplace=True)

intergenic_model_nrp1.replace([np.inf, -np.inf], np.nan, inplace=True)
intergenic_model_nrp1.fillna(0,inplace=True)
# ## to plot genes on top
gene='MEC1'
x=data_wt.loc[gene,'reads-per-tr']
y=intergenic_model.loc[gene]/ref_HO_wt
ax.annotate(gene,(x*(1-0.02) ,y*(1+0.02)),size=10, c='green', bbox=dict(boxstyle="round", fc="w"))


x=data_nrp1.loc[gene,'reads-per-tr']
y=intergenic_model_nrp1.loc[gene]/ref_HO_wt
ax2.annotate(gene,(x*(1-0.02) ,y*(1+0.02)),size=10, c='green', bbox=dict(boxstyle="round", fc="w"))

#%%    
fig.savefig('HO_relative_merged_rates_intergenic_model_vs_reads_per_tr_Greg_'+gene+'.png',dpi=300,format='png')



