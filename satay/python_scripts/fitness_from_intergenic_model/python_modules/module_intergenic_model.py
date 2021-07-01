# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 11:31:37 2021

@author: linigodelacruz
"""

import numpy as np
import pandas as pd


def getting_r(datasets): 
    """
    This function computes the maximum rates per gene,effectively the fitness for 
    a strain with a knockout in gene X, given an intergenic model . The intergenic model assumes a 
    fitine population growth until the carrying capacity (K). This population growth describes a 
    logistic growth. 
    
    Parameters
    ----------
    datasets : TYPE - list of dataframes
        These are dataframes containing the information on insertions per gene and reads per gene 
        given by the transposon data sequencing.
        The columns should have the following names:
            number_of_read_per_gene
            number_of_transposon_per_gene
            

    Returns
    -------
    r : array of length of the list given as input 
        array containing the maximum rates per gene per dataset in the population according an intergenic model
   

    """
    ## add a for loop over times and compute the rates over those times and averaged them out and std 
 
    T=90
    r=[]  
    r_non_filter=[]
 
    
    for i in np.arange(0,len(datasets)):
    
    
       
        K=np.sum(datasets[i]['reads-per-tr'])
        N=datasets[i]['reads-per-tr'] # contribution of each mutant to the population density 
        # it will compute a carrying capacity per dataset 
       
                       
        r.append(np.log(N*K/(K-N))/T)
        
        K_n=np.sum(datasets[i]['Nreadsperinsrt']) # it will compute a carrying capacity per dataset 
        N_n=datasets[i]['Nreadsperinsrt']
       

                       
        r_non_filter.append(np.log(N_n*K_n/(K_n-N_n))/T)
        
              
        
        
    return r,r_non_filter



