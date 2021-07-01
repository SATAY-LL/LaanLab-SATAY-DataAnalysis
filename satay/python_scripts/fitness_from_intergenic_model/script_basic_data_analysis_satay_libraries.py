# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 10:25:17 2021
This script will be used to do basic data nalaysis for the satay libraries according 
what is presenting in (Michel et al., 2017)
@author: linigodelacruz
"""
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from collections import defaultdict
import os
import seaborn as sns
import scipy 
from functools import reduce

from src.python_modules.module_analysis_transposon_sites import *
#%% Import of dataframes output from the SATAY pipeline

names_libraries={'wt':'data_wt_merged.xlsx','dnrp1':'data_dnrp1_merged.xlsx'}
data_library=[]
for i in names_libraries.keys():
 data_library.append(pd.read_excel('datasets/'+names_libraries[i],index_col='Unnamed: 0'))
#### Creating a big dataframe of the libraries

data_library_pd=pd.concat(data_library,keys=names_libraries.keys(),sort=True)
data_library_pd.fillna(0,inplace=True)

#%%

################# Computing measures ######################

freq=frequency_transposons(data_library_pd,names_libraries)
reads_per_tr=reads_per_transposon(data_library_pd,names_libraries)
tr_density=transposon_density(data_library_pd,names_libraries)
median_insertions=median_feature(data_library_pd,names_libraries,'Ninsertions')
median_insert_essentials=median_feature_essentials(data_library_pd,names_libraries,'Ninsertions')
median_insert_nonessentials=median_feature_nonessentials(data_library_pd,names_libraries,'Ninsertions')
#%% Assembling the masure into a dataframe 
analysis_libraries=defaultdict(dict)

j=0
for i in names_libraries.keys():
    
    analysis_libraries[i]['one-transposon-per-bp']=freq[j]
    analysis_libraries[i]['median-reads-per-transposons']=reads_per_tr[j]
    analysis_libraries[i]['median-tr']=median_insertions[j]
    analysis_libraries[i]['median-tr-essentials']=median_insert_essentials[j]
    analysis_libraries[i]['median-tr-non-essentials']=median_insert_nonessentials[j]
    j=j+1
    
del i,j
analysis_libraries_pd=pd.DataFrame(analysis_libraries)

#%% Defining the dataframes per type 
data_wt=data_library_pd.loc['wt'].copy()
data_nrp1=data_library_pd.loc['dnrp1'].copy()

####### Transposon density vs genes ####################

data_nrp1['tr-density']=data_nrp1['Ninsertions']/data_nrp1['Nbasepairs']
data_wt['tr-density']=data_wt['Ninsertions']/data_wt['Nbasepairs']
data_nrp1['reads-density']=data_nrp1['Nreads']/data_nrp1['Nbasepairs']
data_wt['reads-density']=data_wt['Nreads']/data_wt['Nbasepairs']


# data_wt_agnes['tr-density']=data_wt_agnes['Ninsertions']/data_wt_agnes['Nbasepairs']
# data_wt_greg2['tr-density']=data_wt_greg2['Ninsertions']/data_wt_greg2['Nbasepairs']
#%% Get the reads per transpons 
#%% Getting a better measure of the reads per transposons
#- Discard the genes that has less than 25 reads 
#- Discard the genes that has more 20% coverage of transposons (centromeres , bias genes), tr-density>0.2 (less than 3% of the data)
#- Discard the genes that has less than 1% of coverage of transposons , tr-density<0.01 (around 7% of the genes that has less than 20% coverage)

data_wt_filtered=filter_low_and_biased_reads_genes(data_wt)
data_wt=data_wt_filtered
data_wt.fillna(0,inplace=True)
data_nrp1_filtered=filter_low_and_biased_reads_genes(data_nrp1)
data_nrp1=data_nrp1_filtered
data_nrp1.fillna(0,inplace=True)# will put zeros to the discarded regions
#%% Exporting datasets 
data_nrp1.to_excel('datasets/data_nrp1_filtered_reads_per_tr.xlsx')
data_wt.to_excel('datasets/data_wt_filtered_reads_per_tr.xlsx')


#%% Plot transposon density (fig 1B Benoit) highlighting the centromere position

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
ax.plot(data_wt['tr-density'],alpha=0.5,color='b')
ax.set_ylabel('transposon density: tn/bp')
ax.set_xlabel('genes')
## annotated centromeres
for i in np.arange(0,len(data_wt)):
    
    if data_wt.loc[i,'Feature_type']=='Centromere': 
   
        ax.vlines(x=i,ymin=0,ymax=0.8,linestyles='--',alpha=0.3)
        ax.text(x=i,y=0.6,s='centromere',rotation=90,fontsize=8)
#%% Transposon density for comparing two libraries
fig=plt.figure(figsize=(10,9))
grid = plt.GridSpec(2, 1, wspace=0.0, hspace=0.0)
ax = plt.subplot(grid[0,0])
ax2 = plt.subplot(grid[1,0])   

ax.plot(data_wt['tr-density'],alpha=0.5,color='b')
ax.set_ylabel('transposond density: tn/bp')

## annotated centromeres
for i in np.arange(0,len(data_wt)):
    
    if data_wt.loc[i,'Feature_type']=='Centromere': 
   
        ax.vlines(x=i,ymin=0,ymax=0.8,linestyles='--',alpha=0.3)
        ax.text(x=i,y=0.6,s='centromere',rotation=90,fontsize=8)

ax2.plot(data_nrp1['tr-density'],alpha=0.5,color='orange')
ax2.set_ylabel('transposon density: tn/bp dnrp1 ')
ax2.set_xlabel('genes')
## annotated centromeres
for i in np.arange(0,len(data_nrp1)):
    
    if data_nrp1.loc[i,'Feature_type']=='Centromere': 
   
        ax2.vlines(x=i,ymin=0,ymax=0.8,linestyles='--',alpha=0.3)
        ax2.text(x=i,y=0.6,s='centromere',rotation=90,fontsize=8)
#%% saving the figure transposon density
fig.savefig('output_images/Transposon-density-WT-annotated-centromeres-WT-vs-dnrp1.png',dpi=300,format='png',transparent=False)
#%%  Plot reads per transposon  highlighting the centromere position
strain=data_nrp1

fig=plt.figure(figsize=(15,15))
grid = plt.GridSpec(4, 1, wspace=0.0, hspace=0.1)

ax = plt.subplot(grid[0,0])
ax2 = plt.subplot(grid[2,0])   
ax3 = plt.subplot(grid[3,0]) 
ax4 = plt.subplot(grid[1,0])  

ax.set_title('Variation along the genome for dnrp1 merged dataset')
ax.plot(strain['reads-per-tr'],alpha=0.9,color='green')
ax.set_ylabel('reads per tr filtered for high tr-density')

ax.set_ylim(0,2000)
## annotated centromeres
for i in np.arange(0,len(strain)):
    
    if strain.loc[i,'Feature_type']=='Centromere': 
   
        ax.vlines(x=i,ymin=0,ymax=2000,linestyles='--',alpha=0.3)
        #ax.text(x=i,y=4000,s='centromere',rotation=90,fontsize=8)
    elif strain.loc[i,'reads-per-tr']>500:
        ax.vlines(x=i,ymin=0,ymax=2000,linestyles='-',alpha=0.2)
        ax.text(x=i,y=1500,s=strain.loc[i,'Standard_name'],rotation=90,fontsize=8)

ax4.plot(strain['Nreadsperinsrt'],alpha=0.5,color='green')
ax4.set_ylabel('reads per tr without density filter')
ax4.set_ylim(0,2000)
#ax4.set_xticks([])
## annotated centromeres
for i in np.arange(0,len(strain)):
    
    if strain.loc[i,'Feature_type']=='Centromere': 
   
        ax4.vlines(x=i,ymin=0,ymax=2000,linestyles='--',alpha=0.3)
        #ax.text(x=i,y=4000,s='centromere',rotation=90,fontsize=8)
    elif strain.loc[i,'reads-per-tr']>500:
        ax4.vlines(x=i,ymin=0,ymax=2000,linestyles='-',alpha=0.2)
        ax4.text(x=i,y=1500,s=strain.loc[i,'Standard_name'],rotation=90,fontsize=8)       

ax2.plot(strain['tr-density'],alpha=0.7,color='black')
ax2.set_ylabel('transposon density')
ax2.set_ylim(0,1)
## annotated centromeres
for i in np.arange(0,len(data_wt)):
    
    if strain.loc[i,'Feature_type']=='Centromere': 
   
        ax2.vlines(x=i,ymin=0,ymax=1,linestyles='--',alpha=0.3)
        ax2.text(x=i,y=0.5,s='centromere',rotation=90,fontsize=8)
        
ax3.plot(strain['reads-density'],alpha=0.5,color='b')
ax3.set_ylabel('Reads density')
ax3.set_xlabel('genomic regions')
ax3.set_ylim(0,250)
## annotated centromeres
for i in np.arange(0,len(data_wt)):
    
    if strain.loc[i,'Feature_type']=='Centromere': 
   
        ax3.vlines(x=i,ymin=0,ymax=250,linestyles='--',alpha=0.3)
        ax3.text(x=i,y=150,s='centromere',rotation=90,fontsize=8)
    
#%% saving the figure reads per transposon density
fig.savefig('output_images/dnrp1-Variability-along-genome-raw-data.png',dpi=300,format='png',transparent=False)

#%% Compare the fold change of the reads per transposon per library

fold_change=data_wt['reads-per-tr']/data_nrp1['reads-per-tr']
fold_change.replace([np.inf, -np.inf], np.nan, inplace=True)
fold_change.fillna(0,inplace=True)

fold_change_nrp1=data_nrp1['reads-per-tr']/data_wt['reads-per-tr']
fold_change_nrp1.replace([np.inf, -np.inf], np.nan, inplace=True)
fold_change_nrp1.fillna(0,inplace=True)

cutoff=8

fig=plt.figure(figsize=(10,30))
grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.2)

ax = plt.subplot(grid[0,0])
ax1 = plt.subplot(grid[1,0])
ax2 = plt.subplot(grid[2,0])

ax.set_title('Potential Negative Interactors for nrp1')
ax.plot(fold_change,alpha=0.6,color='red')
ax.hlines(y=cutoff,xmin=0,xmax=14000,linestyles='--',alpha=0.3,label='cutoff')

ax.set_ylabel('fold change reads per tr')
ax.set_ylim(0,100)
## annotated centromeres
for i in np.arange(0,len(data_wt)):
    
    
    if fold_change[i]>cutoff and data_wt.loc[i,'Standard_name']!='noncoding':
        ax.vlines(x=i,ymin=0,ymax=100,linestyles='-',alpha=0.2)
        ax.text(x=i,y=60,s=data_wt.loc[i,'Standard_name'],rotation=90,fontsize=8)
    if fold_change[i]>cutoff and data_wt.loc[i,'Essentiality']==1 :
        ax.vlines(x=i,ymin=0,ymax=100,linestyles='-',alpha=0.2)
        ax.text(x=i,y=60,s=data_wt.loc[i,'Standard_name'],rotation=90,fontsize=8,color='red')

ax1.set_title('Potential Positive Interactors for nrp1')
ax1.plot(fold_change_nrp1,alpha=0.6,color='green')
ax1.hlines(y=cutoff,xmin=0,xmax=14000,linestyles='--',alpha=0.3,label='cutoff')
ax1.set_xlabel('genomic regions')
ax1.set_ylabel('fold change reads per tr')
ax1.set_ylim(0,100)
## annotated centromeres
for i in np.arange(0,len(data_wt)):
    
    
    if fold_change_nrp1[i]>cutoff and data_wt.loc[i,'Standard_name']!='noncoding' :
        ax1.vlines(x=i,ymin=0,ymax=100,linestyles='-',alpha=0.2)
        ax1.text(x=i,y=60,s=data_wt.loc[i,'Standard_name'],rotation=90,fontsize=8)
    if fold_change_nrp1[i]>cutoff and data_wt.loc[i,'Essentiality']==1 :
        ax1.vlines(x=i,ymin=0,ymax=100,linestyles='-',alpha=0.2)
        ax1.text(x=i,y=60,s=data_wt.loc[i,'Standard_name'],rotation=90,fontsize=8,color='red')
        
ax.legend()
ax1.legend()

ax2.plot(data_wt['reads-per-tr'],data_nrp1['reads-per-tr'],'o',color='gray')
ax2.set_xlim(0,3000)
ax2.set_ylim(0,3000)
ax2.set_xlabel('reads per tr WT')
ax2.set_ylabel('reads per tr dnrp1')
#%% Compare the fold change of the transposons per library

fold_change=data_wt['tr-density']/data_nrp1['tr-density']
fold_change.replace([np.inf, -np.inf], np.nan, inplace=True)
fold_change.fillna(0,inplace=True)

fold_change_nrp1=1/fold_change
fold_change_nrp1.replace([np.inf, -np.inf], np.nan, inplace=True)
fold_change_nrp1.fillna(0,inplace=True)

cutoff=1.8

fig=plt.figure(figsize=(10,30))
grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.2)

ax = plt.subplot(grid[0,0])
ax1 = plt.subplot(grid[1,0])
ax2 = plt.subplot(grid[2,0])

ax.set_title('Potential Negative Interactors for nrp1')
ax.plot(fold_change,alpha=0.6,color='red')
ax.hlines(y=cutoff,xmin=0,xmax=14000,linestyles='--',alpha=0.3,label='cutoff')

ax.set_ylabel('fold change tr density')
ax.set_ylim(0,15)
## annotated centromeres
for i in np.arange(0,len(data_wt)):
    
    
    if fold_change[i]>cutoff and data_wt.loc[i,'Standard_name']!='noncoding':
        ax.vlines(x=i,ymin=0,ymax=15,linestyles='-',alpha=0.2)
        ax.text(x=i,y=8,s=data_wt.loc[i,'Standard_name'],rotation=90,fontsize=8)
    if fold_change[i]>cutoff and data_wt.loc[i,'Essentiality']==1 :
        ax.vlines(x=i,ymin=0,ymax=15,linestyles='-',alpha=0.2)
        ax.text(x=i,y=8,s=data_wt.loc[i,'Standard_name'],rotation=90,fontsize=8,color='red')
cutoff=7
ax1.set_title('Potential Positive Interactors for nrp1')
ax1.plot(fold_change_nrp1,alpha=0.6,color='green')
ax1.hlines(y=cutoff,xmin=0,xmax=14000,linestyles='--',alpha=0.3,label='cutoff')
ax1.set_xlabel('genomic regions')
ax1.set_ylabel('fold change tr density')
ax1.set_ylim(0,30)
## annotated centromeres

for i in np.arange(0,len(data_wt)):
    
    
    if fold_change_nrp1[i]>cutoff and data_wt.loc[i,'Standard_name']!='noncoding' :
        ax1.vlines(x=i,ymin=0,ymax=30,linestyles='-',alpha=0.2)
        ax1.text(x=i,y=20,s=data_wt.loc[i,'Standard_name'],rotation=90,fontsize=8)
    if fold_change_nrp1[i]>cutoff and data_wt.loc[i,'Essentiality']==1 :
        ax1.vlines(x=i,ymin=0,ymax=30,linestyles='-',alpha=0.2)
        ax1.text(x=i,y=20,s=data_wt.loc[i,'Standard_name'],rotation=90,fontsize=8,color='red')
        
ax.legend()
ax1.legend()

ax2.plot(data_wt['tr-density'],data_nrp1['tr-density'],'o',color='gray')
# ax2.set_xlim(0,3000)
# ax2.set_ylim(0,3000)
ax2.set_xlabel('tr density WT')
ax2.set_ylabel('tr density  dnrp1')
#%% save figure
fig.savefig('output_images/fold_change_wt_vs_nrp1_tr_density.png',dpi=300,format='png',transparent=False)

#%% determine the local variation of transposons along te genome 
## Data per chromosome
#from scipy.stats import sem 
mean_wt_chrom=data_wt.groupby(by='chromosome')['Ninsertions'].mean()
std_wt_chrom=data_wt.groupby(by='chromosome')['Ninsertions'].sem()

mean_wt_chrom_trdensity=data_wt.groupby(by='chromosome')['tr-density'].mean()
std_wt_chrom_trdensity=data_wt.groupby(by='chromosome')['tr-density'].sem()

mean_wt_chrom_readspertr=data_wt.groupby(by='chromosome')['reads-per-tr'].mean()
std_wt_chrom_readspertr=data_wt.groupby(by='chromosome')['reads-per-tr'].sem()

mean_wt_chrom_reads=data_wt.groupby(by='chromosome')['Nreads'].mean()
std_wt_chrom_reads=data_wt.groupby(by='chromosome')['Nreads'].sem()

fig=plt.figure(figsize=(10,9))
grid = plt.GridSpec(4, 1, wspace=0.0, hspace=0.0)
ax = plt.subplot(grid[0,0])
ax.errorbar(data_wt.loc[:,'chromosome'].unique(), mean_wt_chrom, std_wt_chrom, marker='s', mfc='red',
         mec='green', ms=10, mew=1,capsize=4)
ax.set_ylabel('Ninsertions')

ax2 = plt.subplot(grid[1,0])
ax2.errorbar(data_wt.loc[:,'chromosome'].unique(), mean_wt_chrom_trdensity, std_wt_chrom_trdensity, marker='s', mfc='red',
         mec='green', ms=10, mew=1,capsize=4)
ax2.set_ylabel('Ninsertions per bp')

ax3 = plt.subplot(grid[2,0])
ax3.errorbar(data_wt.loc[:,'chromosome'].unique(), mean_wt_chrom_readspertr,std_wt_chrom_readspertr, marker='s', mfc='red',
         mec='green', ms=10, mew=1,capsize=4)
ax3.set_ylabel('Nreads per Ninsertions')

ax4 = plt.subplot(grid[3,0])
ax4.errorbar(data_wt.loc[:,'chromosome'].unique(), mean_wt_chrom_reads,std_wt_chrom_reads, marker='s', mfc='red',
         mec='green', ms=10, mew=1,capsize=4)
ax4.set_ylabel('Nreads')
#%% saving the figure 
fig.savefig('Variability-along-genome-sem.png',dpi=300,format='png',transparent=False)
#%% assesing local variation per chromosome. Viz per chromosome

magnitudes=['Ninsertions','tr-density','reads-per-tr']
chromosomes=data_wt.loc[:,'chromosome'].unique()

windows=5
chrom=chromosomes[5]


fig=plt.figure(figsize=(10,9))
grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.0)
ax = plt.subplot(grid[0,0])
ax.set_title('Errorbars over chrom='+str(chrom)+', every '+str(windows)+'genes')
mean_over_chromI,std_over_chromI=local_variation(chrom=chrom, windows=windows, data=data_wt,column=magnitudes[0])
chrom_data=pd.DataFrame([mean_over_chromI,std_over_chromI],index=['mean','std'])


ax.errorbar(np.arange(0,len(data_wt[data_wt.loc[:,'chromosome']==chrom]),windows), chrom_data.loc['mean',:], chrom_data.loc['std',:], marker='s', 
            mfc='red', mec='green', ms=10, mew=1,capsize=4)
ax.set_ylabel(magnitudes[0])

mean_over_chromI,std_over_chromI=local_variation(chrom=chrom, windows=windows, data=data_wt,column=magnitudes[1])
chrom_data=pd.DataFrame([mean_over_chromI,std_over_chromI],index=['mean','std'])

ax2= plt.subplot(grid[1,0])
ax2.errorbar(np.arange(0,len(data_wt[data_wt.loc[:,'chromosome']==chrom]),windows), chrom_data.loc['mean',:], chrom_data.loc['std',:], marker='s', mfc='red',
         mec='green', ms=10, mew=1,capsize=4)
ax2.set_ylabel(magnitudes[1])

mean_over_chromI,std_over_chromI=local_variation(chrom=chrom, windows=windows, data=data_wt,column=magnitudes[2])
chrom_data=pd.DataFrame([mean_over_chromI,std_over_chromI],index=['mean','std'])

ax3= plt.subplot(grid[2,0])
ax3.errorbar(np.arange(0,len(data_wt[data_wt.loc[:,'chromosome']==chrom]),windows), chrom_data.loc['mean',:], chrom_data.loc['std',:], marker='s', mfc='red',
         mec='green', ms=10, mew=1,capsize=4)
ax3.set_ylabel(magnitudes[2])
#ax.set_xlabel('genes along the windows')

#%% assesing local variation per chromosome. Visualizing more than one chromosome

magnitudes=['Ninsertions','tr-density','reads-per-tr']
chromosomes=data_wt.loc[:,'chromosome'].unique()

windows=5

chroms=[chromosomes[0],chromosomes[1],chromosomes[2],chromosomes[3]]

fig=plt.figure(figsize=(20,20))
grid = plt.GridSpec(3, len(chroms), wspace=0.0, hspace=0.0)

for i in np.arange(0,len(chroms)):
    for j in np.arange(0,len(magnitudes)): 
        ax = plt.subplot(grid[j,i])
       
        mean_over_chromI,std_over_chromI=local_variation(chrom=chroms[i], windows=windows, data=data_wt,column=magnitudes[j])
        chrom_data=pd.DataFrame([mean_over_chromI,std_over_chromI],index=['mean','std'])
        
        
        ax.errorbar(np.arange(0,len(data_wt[data_wt.loc[:,'chromosome']==chroms[i]]),windows), chrom_data.loc['mean',:], chrom_data.loc['std',:], marker='s', 
                    mfc='red', mec='green', ms=10, mew=1,capsize=4)
        
        
    ax.set_title('Errorbars over chrom='+str(chroms[i])+', every '+str(windows)+'genes')

fig.text(0.08, 0.5, magnitudes[1], va='center', rotation='vertical',fontsize=15)  
fig.text(0.08, 0.75, magnitudes[0], va='center', rotation='vertical',fontsize=15)   
fig.text(0.08, 0.25, magnitudes[2], va='center', rotation='vertical',fontsize=15)  
# axlabel.set_ylabel(magnitudes[0],fontsize=15)     
#ax.set_xlabel('genes along the windows')

#%% saving the figure
fig.savefig('local-variation-merged_wt-windows-'+str(windows)+'genes-'+ str(chroms)+'.png',dpi=300,format='png',transparent=False)
#%% Histograms of number of transposons per gene
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
ax.set_xlabel('Number of transposons per gene')
ax.set_ylabel('number of genes (CDS)')
ax.hist(data_wt['Ninsertions'],bins=300,color='gray',alpha=0.7)
ax.set_xlim(0,200)
ax.vlines(x=data_wt_greg2['Ninsertions'].median(),ymin=0,ymax=1400,linestyles='--')
ax.text(x=data_wt_greg2['Ninsertions'].median(),y=1200,s='median all genes')
#%% saving the figure
fig.savefig('Number-of-transposons-per-gene-Greg2.png',dpi=300,format='png',transparent=False)


#%% #%% Histograms of number of transposon per gene according their essentiality
### create plot
fig=plt.figure(figsize=(10,9))
grid = plt.GridSpec(2, 1, wspace=0.0, hspace=0.0)

essentials=data_wt_greg2[data_wt_greg2.loc[:,'Essentiality']==1]['Ninsertions']
nonessentials=data_wt_greg2[data_wt_greg2.loc[:,'Essentiality']==0]['Ninsertions']


ax = plt.subplot(grid[0,0])
sns.histplot(nonessentials,binwidth=2,color='gray',alpha=0.7)
max_x = ax.get_xlim()
ax.set_xlim(left=0,right=200)
ax.grid(True)
ax.set_xticklabels([])
ax.vlines(x=nonessentials.median(),ymin=0,ymax=300,linestyles='--')
ax.text(x=nonessentials.median(),y=300,s='median-non-essentials')
ax.set_ylabel('Annotated non essential genes')

ax2 = plt.subplot(grid[1,0])    
sns.histplot(essentials,binwidth=2,color='orange',alpha=0.7)
ax2.invert_yaxis()
ax2.set_xlim(0,200)
ax2.grid(True)
ax2.vlines(x=essentials.median(),ymin=0,ymax=40,linestyles='--')
ax2.text(x=essentials.median(),y=40,s='median-essentials')

ax2.set_xlabel('Number of transposons per gene')
ax2.set_ylabel('Annotated essential genes')

plt.show()
#%% saving the figure
fig.savefig('Number-of-transposons-per-gene-according-essentiality-Greg2.png',dpi=300,format='png',transparent=False)
#%% Looking for regions devoid of transposons 

difficult_genes=data_wt_agnes[data_wt_agnes['Ninsertions']==0]

difficult_genes_essentials=difficult_genes[difficult_genes['Essentiality']==1]

#%% Loooking at correlations with essential genes 
# fig = plt.figure(figsize=(11,5))
# ax = fig.add_subplot(111)
h=sns.pairplot(data=data_wt,vars=['Ninsertions','tr-density',  "reads-per-tr"],hue='Essentiality',corner=True,diag_kind="kde")
h.map_lower(sns.kdeplot, levels=5, color="b")

#%% saving the figure
h.savefig('pairplot-essentiality-WT.png',dpi=300,format='png',transparent=False)

