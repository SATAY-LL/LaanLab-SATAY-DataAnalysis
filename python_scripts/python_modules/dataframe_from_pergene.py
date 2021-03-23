# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:17:11 2021

@author: gregoryvanbeek
"""

import os, sys
import numpy as np
import pandas as pd
import re


#%%INPUT
pergenefile = r""

#%%
def dataframe_from_pergenefile(pergenefile, verbose=True):
    '''
    This function creates a dataframe with the information from a pergene.txt file.
    Input is a path to a pergene.txt file
    Output is a dataframe where each row is a single gene and with the following columns:
        - gene_names
        - gene_essentiality
        - tn_per_gene
        - read_per_gene
        - Nreadsperinsrt
    
    The gene_essentiality is created based on the genes present in the Cerevisiae_EssentialGenes_List_1.txt and Cerevisiae_EssentialGenes_List_2.txt files
    The number of reads per insertion (Nreadsperinsrt) is determined by dividing the read_per_gene column by the tn_per_gene column.
    
    A more extensive version of this function is the python script genomicfeatures_dataframe.py found in the python_scripts folder.
    '''

    file_dirname = os.path.dirname(os.path.abspath('__file__'))
    sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
    from essential_genes_names import list_known_essentials #import essential_genes_names from python modules directory


# read file
    assert os.path.isfile(pergenefile), 'File not found at: %s' % pergenefile

    with open(pergenefile) as f:
        lines = f.readlines()[1:] #skip header

    genenames_list = [np.nan]*len(lines)
    tnpergene_list = [np.nan]*len(lines)
    readpergene_list = [np.nan]*len(lines) 

    line_counter = 0
    for line in lines:
        line_split = re.split(' |\t', line.strip('\n'))
        l = [x for x in line_split if x]

        if len(l) == 3:
            genenames_list[line_counter] = l[0]
            tnpergene_list[line_counter] = int(l[1])
            readpergene_list[line_counter] = int(l[2])

            line_counter += 1

    del (line, l, line_counter, pergenefile)


# determine number of reads per insertion per gene
    readperinspergene_list = [np.nan]*len(lines)
    for i in range(len(tnpergene_list)):
        if not tnpergene_list[i] == 0:
            readperinspergene_list[i] = readpergene_list[i] / tnpergene_list[i]
        else:
            readperinspergene_list[i] = 0

    del (i)


# determine essential genes
    if os.path.isfile(os.path.join(file_dirname,'..','..','data_files','Cerevisiae_EssentialGenes_List_1.txt')):
        known_essential_gene_list = list_known_essentials(input_files=[os.path.join(file_dirname,'..','..','data_files','Cerevisiae_EssentialGenes_List_1.txt'),
                                                                       os.path.join(file_dirname,'..','..','data_files','Cerevisiae_EssentialGenes_List_2.txt')],
                                                          verbose=verbose)
    elif os.path.isfile(os.path.join(file_dirname,'..','data_files','Cerevisiae_EssentialGenes_List_1.txt')):
        known_essential_gene_list = list_known_essentials(input_files=[os.path.join(file_dirname,'..','data_files','Cerevisiae_EssentialGenes_List_1.txt'),
                                                                       os.path.join(file_dirname,'..','data_files','Cerevisiae_EssentialGenes_List_2.txt')],
                                                          verbose=verbose)

    geneessentiality_list = [None]*len(lines)
    for i in range(len(genenames_list)):
        if genenames_list[i] in known_essential_gene_list:
            geneessentiality_list[i] = True
        else:
            geneessentiality_list[i] = False

    del (lines, file_dirname, known_essential_gene_list, i)

# create dataframe
    read_gene_dict = {"gene_names": genenames_list,
                      "gene_essentiality": geneessentiality_list,
                      "tn_per_gene": tnpergene_list,
                      "read_per_gene": readpergene_list,
                      "Nreadsperinsrt": readperinspergene_list}

    read_gene_df = pd.DataFrame(read_gene_dict, columns = [column_name for column_name in read_gene_dict])

    del (read_gene_dict, genenames_list, geneessentiality_list, tnpergene_list, readpergene_list, readperinspergene_list)
    
    
    return(read_gene_df)


#%%
if __name__ == '__main__':
    read_gene_df = dataframe_from_pergenefile(pergenefile=pergenefile)
