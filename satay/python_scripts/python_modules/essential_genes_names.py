# =============================================================================
# To import the functions, type:
#    import sys
#    sys.path.insert(1,r'full\path\to\this\file')
#    from essential_genes import list_known_essentials
# =============================================================================

import os
from satay.transposonmapping.importing import (
    load_default_files,
)

def list_known_essentials(input_files = None, headerlines=3, verbose=True):
    ''' Get all known essential genes from two different files and combine them in one list.
        Input is a list of of paths where files can be found with the known essential genes.
        A default list is implemented using two files present in the same folder as this file.
        It is expected that the files contain genes in a single column and nothing else.
        An option can be set for the number of headerlines, which by default is set to 3.
        The output is a list containing all the genes present in all files given in the input.
        
        Note when using the default files: The length of the output list exceed the number of known essential genes as the list sometimes contains both the standard name and the systematic name of a gene.
    '''
    
    if input_files == None:
        
        # file_dirname = os.path.dirname(os.path.abspath('__file__'))
        # if os.path.isfile(os.path.join(file_dirname,'..','data_files','Cerevisiae_EssentialGenes_List_1.txt')) and os.path.isfile(os.path.join(file_dirname,'..','data_files','Cerevisiae_EssentialGenes_List_2.txt')):
        #     essential_genes_files = [os.path.join(file_dirname,'..','data_files','Cerevisiae_EssentialGenes_List_1.txt'),
        #                              os.path.join(file_dirname,'..','data_files','Cerevisiae_EssentialGenes_List_2.txt')]
        # else:
        #     essential_genes_files = [os.path.join(file_dirname,'..','data_files','Cerevisiae_EssentialGenes_List_1.txt'),
        #                              os.path.join(file_dirname,'..','data_files','Cerevisiae_EssentialGenes_List_2.txt')]
        # If necessary, load default files
        gff_file, essential_genes_files, gene_name_file = load_default_files(
            gff_file=None, essentials_file=None, gene_names_file=None
        )

        # Verify presence of files
        data_files = {
            "gff3": gff_file,
            "essentials": essential_genes_files,
            "gene_names": gene_name_file,
        }

        for filetype, file_path in data_files.items():
            assert file_path, f"{filetype} not found at {file_path}"
        
        input_files=essential_genes_files

    else:
        essential_genes_files = input_files


    known_essential_gene_list = []
    #making sure essential_genes_files is a list
    if type(essential_genes_files)!=list:
        essential_genes_files=[input_files]
        
    for files in essential_genes_files:
        if verbose == True:
            print('Reading file :',files)
        with open(files) as f:
            for header_lines in range(headerlines):
                next(f)
            for lines in f:
                known_essential_gene_list.append(lines.rstrip('\n'))
                
    
    return(known_essential_gene_list)
#%%
if __name__ == '__main__':
    known_essential_gene_list = list_known_essentials()


