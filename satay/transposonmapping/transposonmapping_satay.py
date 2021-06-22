#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a tool developed for analysing transposon insertions for experiments using SAturated Transposon Analysis in Yeast (SATAY).
This python code contains one function called transposonmapper().
For more information about this code and the project, see github.com/Gregory94/LaanLab-SATAY-DataAnalysis

This code is based on the Matlab code created by the Kornmann lab which is available at: sites.google.com/site/satayusers/

__Author__ = Gregory van Beek. LaanLab, department of Bionanoscience, Delft University of Technology
__version__ = 1.5
__Date last update__ = 2021-01-11

Version history:
    1.1; Added code for creating two text files for storing insertion locations per gene and per essential gene [2020-07-27]
    1.2; Improved searching algorithm for essential genes [2020-08-06]
    1.3; Load file containing all essential genes so that a search for essential genes in multiple file is not needed anymore. This file is created using Create_EssentialGenes_list.py located in the same directory as this code [2020-08-07]
    1.4; Fixed bug where the gene position and transposon insertion location did not start at zero for each chromosome, causing confusing values to be stored in the _pergene_insertions.txt and _peressential_insertions.txt files [2020-08-09]
    1.5; Added functionality to handle all possible sam flags in the alignment file (bam-file) instead of only flag=0 or flag=16. This is needed for the function to handle paired-end sequencing data [2021-01-11]
"""

import os
import numpy as np
import pysam

# Local imports
from .python_modules import chromosomename_roman_to_arabic
from .python_modules import get_chromosome_names
from .python_modules import get_sequence_length
from .python_modules import get_chromosome_reads
from .Files import Files
from .python_modules import get_reads
from .python_modules import read_genes


#%%
def transposonmapper(bamfile, gfffile=None, essentialfiles=None, genenamesfile=None):
    '''
    This function is created for analysis of SATAY data using the species Saccharomyces Cerevisiae.
    It outputs the following files that store information regarding the location of all insertions:
        - .bed-file: Includes all individual basepair locations of the whole genome where at least one transposon has been mapped and the number of insertions for each locations (the number of reads) according to the Browser Extensible Data (bed) format.
                    A distinction is made between reads that had a different reading orientation during sequencing. The number of reads are stored using the equation #reads*20+100 (e.g. 2 reads is stored as 140).
        - .wig-file: Includes all individual basepair locations of the whole genome where at least one transposon has been mapped and the number of insertions for each locations (the number of reads) according to the Wiggle (wig) format.
                    In this file no distinction is made between reads that had a different reading orientation during sequencing. The number of reads are stored as the absolute count.
        - _pergene.txt-file: Includes all genes (currently 6600) with the total number of insertions and number of reads within the genomic region of the gene.
        - _peressential.txt-file: Includes all annotated essential genes (currently 1186) with the total number of insertions and number of reads within the genomic region of the gene.
        - _pergene_insertions.txt-file: Includes all genes with their genomic location (i.e. chromosome number, start and end position) and the locations of all insertions within the gene location. It also include the number number of reads per insertions.
        - _peressential_insertions.txt-file: Includes all essential genes with their genomic location (i.e. chromosome number, start and end position) and the locations of all insertions within the gene location. It also include the number number of reads per insertions.
          (note that in the latter two files, the genomic locations are continous, for example chromosome II does not start at 0, but at 'length chromosome I + 1' etc.).
    The output files are saved at the location of the input file using the same name as the input file, but with the corresponding extension.
    
    The function assumes that the reads are already aligned to a reference genome.
    The input data should be a .bam-file and the location where the .bam-file is stored should also contain an index file (.bam.bai-file, which for example can be created using sambamba).
    This function takes the following inputs:
        - bamfile [required]: Path to the bamfile. This location should also contain the .bam.bai index file (does not need to be input in this function).
        - gfffile [optional]: Path to a .gff-file including all gene information (e.g. downloaded from SGD). Default file is 'Saccharomyces_cerevisiae.R64-1-1.99.gff3'.
        - essentialfiles [optional]: Path to a .txt file containing a list all essential genes. Every line should consist of a single essential gene and the file should have one header line. Ideally this file is created using 'Create_EssentialGenes_list.py'. Default file is 'Cerevisiae_AllEssentialGenes_List.txt'.
        - genenamesfile [optional]: Path to text file that includes aliases for all genes. Default file is 'Yeast_Protein_Names.txt'.
    When the arguments for the optional files are not given, the files are used that are stored at the following location:
        "path_current_pythonscript/../data_files"
    The function uses the pysam package for handling bam files (see pysam.readthedocs.io/en/latest/index.html) and therefore this function only runs on Linux systems with SAMTools installed.
    '''
    
    filename = os.path.basename(bamfile)
   
#%% LOADING ADDITIONAL FILES
    
    # Load files paths into Files() object
    files = Files(bam_file=bamfile, gff_file=gfffile, essentials_file=essentialfiles, gene_names_file=genenamesfile)

    # Read bam file
    bam = pysam.AlignmentFile(files.bam_file, 'rb') #open bam formatted file for reading

    # Get names of all chromosomes as stored in the bam file
    ref_tid = get_chromosome_names(bam)
    ref_names = list(ref_tid.keys())
    
    # Convert chromosome names in data file to roman numerals
    ref_romannums = chromosomename_roman_to_arabic()[1]
    ref_tid_roman = {key: value for key, value in zip(ref_romannums, ref_tid)}

    # Get sequence lengths of all chromosomes
    chr_lengths, chr_lengths_cumsum = get_sequence_length(bam)

    # Get all reads within a specified genomic region
    [readnumb_array, tncoordinates_array, tncoordinatescopy_array] = get_reads(bam)

    # Read files for all genes and all essential genes
    print('Getting coordinates of all genes ...')
    [gene_coordinates, essential_coordinates, aliases_designation] = read_genes(files.gff_file, files.essentials_file, files.gene_names_file)

#%% CONCATENATE ALL CHROMOSOMES

    #FOR EACH INSERTION LOCATION, ADD THE LENGTH OF ALL PREVIOUS CHROMOSOMES.
    ll = 0
    for ii in range(1,len(ref_names)):
        ll += chr_lengths[ref_names[ii-1]]
        aa = np.where(tncoordinatescopy_array[:,0] == ii + 1)
        tncoordinatescopy_array[aa,1] = tncoordinatescopy_array[aa,1] + ll



    #FOR EACH GENE LOCATION, ADD THE LENGTH OF ALL PREVIOUS CHROMOSOMES.
    for key in gene_coordinates:
        gene_chrom = ref_tid_roman.get(gene_coordinates.get(key)[0])
        gene_coordinates[key][1] = gene_coordinates.get(key)[1] + chr_lengths_cumsum.get(gene_chrom)
        gene_coordinates[key][2] = gene_coordinates.get(key)[2] + chr_lengths_cumsum.get(gene_chrom)



    #FOR EACH ESSENTIAL GENE LOCATION, ADD THE LENGTH OF ALL PREVIOUS CHROMOSOMES.
    for key in essential_coordinates:
        gene_chrom = ref_tid_roman.get(essential_coordinates.get(key)[0])
        essential_coordinates[key][1] = essential_coordinates.get(key)[1] + chr_lengths_cumsum.get(gene_chrom)
        essential_coordinates[key][2] = essential_coordinates.get(key)[2] + chr_lengths_cumsum.get(gene_chrom)



    del (ii, ll, aa, key, gene_chrom)



# GET NUMBER OF TRANSPOSONS AND READS PER GENE
    print('Get number of insertions and reads per gene ...')
    
    #ALL GENES
    tnpergene_dict = {}
    readpergene_dict = {}
    tncoordinates_pergene_dict = {}
    # readpergenecrude_dict = {}
    for gene in gene_coordinates:
        xx = np.where(np.logical_and(tncoordinatescopy_array[:,1] >= gene_coordinates.get(gene)[1], tncoordinatescopy_array[:,1] <= gene_coordinates.get(gene)[2])) #get all insertions within range of current gene
        tnpergene_dict[gene] = np.size(xx)
        readpergene_dict[gene] = sum(readnumb_array[xx]) - max(readnumb_array[xx], default=0) #REMOVE LARGEST VALUE TO REDUCE NOISE
        # readpergenecrude_dict[gene] = sum(readnumb_array[xx])

        if np.size(xx) > 0:
            tncoordinates_pergene_dict[gene] = [gene_coordinates.get(gene)[0], gene_coordinates.get(gene)[1], gene_coordinates.get(gene)[2], list(tncoordinatescopy_array[xx[0][0]:xx[0][-1]+1, 1]), list(readnumb_array[xx])]
        else:
            tncoordinates_pergene_dict[gene] = [gene_coordinates.get(gene)[0], gene_coordinates.get(gene)[1], gene_coordinates.get(gene)[2], [], []]


    #ONLY ESSENTIAL GENES
    tnperessential_dict = {}
    readperessential_dict = {}
    tncoordinates_peressential_dict = {}
    # readperessentialcrude_dict = {}
    for gene in essential_coordinates:
        xx = np.where(np.logical_and(tncoordinatescopy_array[:,1] >= essential_coordinates.get(gene)[1], tncoordinatescopy_array[:,1] <= essential_coordinates.get(gene)[2]))
        tnperessential_dict[gene] = np.size(xx)
        readperessential_dict[gene] = sum(readnumb_array[xx]) - max(readnumb_array[xx], default=0)
        # readperessentialcrude_dict[gene] = sum(readnumb_array[xx])

        if np.size(xx) > 0:
            tncoordinates_peressential_dict[gene] = [essential_coordinates.get(gene)[0], essential_coordinates.get(gene)[1], essential_coordinates.get(gene)[2], list(tncoordinatescopy_array[xx[0][0]:xx[0][-1]+1, 1]), list(readnumb_array[xx])]
        else:
            tncoordinates_peressential_dict[gene] = [essential_coordinates.get(gene)[0], essential_coordinates.get(gene)[1], essential_coordinates.get(gene)[2], [], []]



    del (xx, gene)



# CREATE BED FILE
    bedfile = bamfile+'.bed'
    print('Writing bed file at: ', bedfile)
    print('')


    with open(bedfile, 'w') as f:
        
        f.write('track name=' + filename + ' useScore=1\n')
        
        coordinates_counter = 0
        for tn in tncoordinates_array:
            refname = [key for key, val in ref_tid.items() if val == tn[0] - 1][0]
            if refname == 'Mito':
                refname = 'M'
            f.write('chr' + refname + ' ' + str(tn[1]) + ' ' + str(tn[1] + 1) + ' . ' + str(100+readnumb_array[coordinates_counter]*20) + '\n')
            coordinates_counter += 1



    del (bedfile, coordinates_counter, refname)

# CREATE TEXT FILE WITH TRANSPOSONS AND READS PER GENE
# NOTE THAT THE TRANSPOSON WITH THE HIGHEST READ COUNT IS IGNORED.
# E.G. IF THIS FILE IS COMPARED WITH THE _PERGENE_INSERTIONS.TXT FILE THE READS DON'T ADD UP (SEE https://groups.google.com/forum/#!category-topic/satayusers/bioinformatics/uaTpKsmgU6Q)
# TOO REMOVE THIS HACK, CHANGE THE INITIALIZATION OF THE VARIABLE readpergene
    pergenefile = bamfile+'_pergene.txt'
    print('Writing pergene.txt file at: ', pergenefile)
    print('')

    with open(pergenefile, 'w') as f:

        f.write('Gene name\tNumber of transposons per gene\tNumber of reads per gene\n')

        for gene in tnpergene_dict:
            tnpergene = tnpergene_dict[gene]
            readpergene = readpergene_dict[gene]
            if gene in aliases_designation:
                gene_alias = aliases_designation.get(gene)[0]
            else:
                gene_alias = gene
            f.write(gene_alias + '\t' + str(tnpergene) + '\t' + str(readpergene) + '\n')



    del (pergenefile, gene, gene_alias, tnpergene, readpergene)

# CREATE TEXT FILE TRANSPOSONS AND READS PER ESSENTIAL GENE
    peressentialfile = bamfile+'_peressential.txt'
    print('Writing peressential.txt file at: ',peressentialfile)
    print('')

    with open(peressentialfile, 'w') as f:
        
        f.write('Gene name\tNumber of transposons per gene\tNumber of reads per gene\n')
        
        for essential in tnperessential_dict:
            tnperessential = tnperessential_dict[essential]
            readperessential = readperessential_dict[essential]
            if essential in aliases_designation:
                essential_alias = aliases_designation.get(essential)[0]
            else:
                essential_alias = essential
            f.write(essential_alias + '\t' + str(tnperessential) + '\t' + str(readperessential) + '\n')



    del (peressentialfile, essential, essential_alias, tnperessential, readperessential)


# CREATE TEXT FILE WITH LOCATION OF INSERTIONS AND READS PER GENE
    pergeneinsertionsfile = bamfile+'_pergene_insertions.txt'
    print('Witing pergene_insertions.txt file at: ',pergeneinsertionsfile)
    print('')

    with open(pergeneinsertionsfile, 'w') as f:

        f.write('Gene name\tChromosome\tStart location\tEnd location\tInsertion locations\tReads per insertion location\n')

        for gene in tncoordinates_pergene_dict:
            gene_chrom = ref_tid_roman.get(gene_coordinates.get(gene)[0])
            tncoordinates = [ins - chr_lengths_cumsum.get(gene_chrom) for ins in tncoordinates_pergene_dict[gene][3]]

            if gene in aliases_designation:
                gene_alias = aliases_designation.get(gene)[0]
            else:
                gene_alias = gene

            f.write(gene_alias + '\t' + str(tncoordinates_pergene_dict[gene][0]) + '\t' + str(tncoordinates_pergene_dict[gene][1] - chr_lengths_cumsum.get(gene_chrom)) + '\t' + str(tncoordinates_pergene_dict[gene][2] - chr_lengths_cumsum.get(gene_chrom)) + '\t' + str(tncoordinates) + '\t' + str(tncoordinates_pergene_dict[gene][4]) + '\n')



    del (gene, gene_chrom, tncoordinates, gene_alias, pergeneinsertionsfile)

# CREATE TEXT FILE WITH LOCATION OF INSERTIONS AND READS PER ESSENTIAL GENE
    peressentialinsertionsfile = bamfile+'_peressential_insertions.txt'
    print('Writing peressential_insertions.txt file at: ', peressentialinsertionsfile)
    print('')

    with open(peressentialinsertionsfile, 'w') as f:

        f.write('Essential gene name\tChromosome\tStart location\tEnd location\tInsertion locations\tReads per insertion location\n')

        for essential in tncoordinates_peressential_dict:
            gene_chrom = ref_tid_roman.get(gene_coordinates.get(essential)[0])
            tncoordinates = [ins - chr_lengths_cumsum.get(gene_chrom) for ins in tncoordinates_peressential_dict[essential][3]]

            if essential in aliases_designation:
                essential_alias = aliases_designation.get(essential)[0]
            else:
                essential_alias = essential

            f.write(essential_alias + '\t' + str(tncoordinates_peressential_dict[essential][0]) + '\t' + str(tncoordinates_peressential_dict[essential][1] - chr_lengths_cumsum.get(gene_chrom)) + '\t' + str(tncoordinates_peressential_dict[essential][2] - chr_lengths_cumsum.get(gene_chrom)) + '\t' + str(tncoordinates) + '\t' + str(tncoordinates_peressential_dict[essential][4]) + '\n')



    del (essential, gene_chrom, tncoordinates, essential_alias, peressentialinsertionsfile)

# ADD INSERTIONS AT SAME LOCATION BUT WITH DIFFERENT ORIENTATIONS TOGETHER (FOR STORING IN WIG-FILE)
    wigfile = bamfile+'.wig'
    print('Writing wig file at: ', wigfile)
    print('')


    readnumbwig_array = readnumb_array.copy()



    unique_index_array = np.array([], dtype=int) #=cc
    N_uniques_perchr_list = []
    ll = 0
    for kk in ref_names:
        index = np.where(tncoordinates_array[:,0] == int(ref_tid[kk]+1)) #get indices for current chromosome.
        unique_index = np.unique(tncoordinates_array[index][:,1], return_index=True)[1] #get all insertion locations (in tncoordinates, all rows, column 1)

        unique_index_array = np.append(unique_index_array, (unique_index+ll), axis=0)

        ll += np.count_nonzero(tncoordinates_array[:,0] == int(ref_tid[kk]+1))
        N_uniques_perchr_list.append(ll) #total amount unique indices found untill current chromosome



    del (ll, kk, unique_index)



    duplicate_list = [] #=dd
    ll = 0
    index_last_unique_previous_chromosome = 0
    for ii in N_uniques_perchr_list:
        index_last_unique = np.where(unique_index_array <= ii)[0][-1]
        for jj in range(ll,ii):
            if int(jj) not in unique_index_array[index_last_unique_previous_chromosome:index_last_unique]:
                duplicate_list.append(jj)
        index_last_unique_previous_chromosome = index_last_unique
        ll = ii

    #SUM READNUMB VALUES AT INDEX IN DUPLICATE_LIST AND DUPLICATE_LIST-1    
    for ii in duplicate_list:
        readnumbwig_array[ii-1] = readnumbwig_array[ii-1] + readnumbwig_array[ii]
    
    tncoordinateswig_duplicatesremoved_array = np.delete(tncoordinates_array, duplicate_list, axis=0)
    readnumbwig_duplicatesremoved_array = np.delete(readnumbwig_array, duplicate_list, axis=0)




    del (ll, ii, jj, N_uniques_perchr_list, index_last_unique, duplicate_list, readnumbwig_array)

# CREATING WIG FILE    
    with  open(wigfile, 'w') as f:
        f.write('track type=wiggle_0 ,maxheightPixels=60 name='+filename+'\n')
        for kk in ref_names:
            f.write('VariableStep chrom=chr' + kk + '\n')

            index = np.where(tncoordinateswig_duplicatesremoved_array[:,0] == int(ref_tid[kk]+1)) #get indices for current chromosome.
            for ii in index[0]:
                f.write(str(tncoordinateswig_duplicatesremoved_array[ii][1]) + ' ' + str(readnumbwig_duplicatesremoved_array[ii]) + '\n')


    del (wigfile, kk, ii, index)

#%%
if __name__ == '__main__':
    bamfile= 'satay/data_files/files4test/SRR062634.filt_trimmed.sorted.bam'
    transposonmapper(bamfile=bamfile)



