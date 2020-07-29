#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a tool developed for analysing transposon insertions for experiments using SAturated Transposon Analysis in Yeast (SATAY).
For more information about this code and the project, see github.com/Gregory94/LaanLab-SATAY-DataAnalysis

This code is based on the Matlab code created by the Kornmann lab and is available at: sites.google.com/site/satayusers/
It contains one function called transposonmapper().

__Author__ = Gregory van Beek. LaanLab, department of Bionanoscience, Delft University of Technology
__version__ = 1.0
__Date last update__ = 2020-07-27
"""

import os, sys
import warnings
import timeit
import numpy as np
import pysam

dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(dirname,'python_modules'))
from chromosome_and_gene_positions import chromosomename_roman_to_arabic, gene_position
from gene_names import gene_aliases
from essential_genes_names import list_known_essentials



bam_arg = sys.argv[1]



#%%
def transposonmapper(bamfile=bam_arg, gfffile=None, essentialfiles=None, genenamesfile=None):
    '''
    This function is created for analysis of SATAY data using the species Saccharomyces Cerevisiae.
    It outputs the following files that store information regarding the location of all insertions:
        - .bed-file: Include all individual basepair locations of the whole genome where at least one transposon has been mapped and the number of insertions for each locations (the number of reads) according to the Browser Extensible Data (bed) format.
                    A distinction is made between reads that had a different reading orientation during sequencing. The number of reads are stored using the equation #reads*20+100 (e.g. 2 reads is stored as 140).
        - .wig-file: Include all individual basepair locations of the whole genome where at least one transposon has been mapped and the number of insertions for each locations (the number of reads) according to the Wiggle (wig) format.
                    In this file no distinction is made between reads that had a different reading orientation during sequencing. The number of reads are stored as the absolute count.
        - _pergene.txt-file: Include all genes (currently 6600) with the total number of insertions and number of reads within the genomic region of the gene.
    The output files are saved at the location of the input file using the same name as the input file, but with the corresponding extension.
    
    The function assumes that the reads are already aligned to a reference genome.
    The input data should be a .bam-file and the location where the .bam-file is stored should also contain an index file (.bai-file, which for example can be created using sambamba).
    It takes the following inputs:
        - bamfile [required]: Path to the bamfile. This location should also contain the index file (.bai).
        - gfffile [optional]: Path to a .gff-file including all gene information (e.g. downloaded from SGD).
        - essentialfiles [optional]: List of paths with text files that includes all essential gene names.
        - genenamesfile [optional]: Path to text file that includes aliases for all genes.
    When the arguments for the optional files are not given, the files are used that are stored at the following location:
        "path_current_pythonscript/../data_files"
    The function uses the pysam package for handling bam files (see pysam.readthedocs.io/en/latest/index.html)
    '''

#%% LOADING BAM FILE
    if bamfile is None:
        path = os.path.join('/home', 'gregoryvanbeek', 'Documents', 'data_processing')
        filename = 'E-MTAB-4885.WT1.bam'
        # filename = 'SRR062634.filt_trimmed.sorted.bam'
        bamfile = os.path.join(path,filename)
    else:
        filename = os.path.basename(bamfile)
        path = bamfile.replace(filename,'')

    
    if os.path.isfile(bamfile):
        print('Running: ', bamfile)
    else:
        raise ValueError('Bam file not found at: ', bamfile)



#%% LOADING ADDITIONAL FILES
    files_path = os.path.join(dirname,'..','..','data_files')
    
    #LOADING GFF-FILE
    if gfffile is None:
        gfffile = os.path.join(files_path,'Saccharomyces_cerevisiae.R64-1-1.99.gff3')
    if not os.path.isfile(gfffile):
        raise ValueError('Path to GFF-file does not exist.')
    
    #LOADING TEXT FILES WITH ESSENTIAL GENES
    if essentialfiles is None:
        essentialfiles = [os.path.join(files_path,'Cervisiae_EssentialGenes_List_1.txt'), os.path.join(files_path,'Cervisiae_EssentialGenes_List_2.txt')]
    for essential_file in essentialfiles:
        if not os.path.isfile(essential_file):
            raise ValueError('Following path does not exist: ' + essential_file)
    del essential_file
    
    #LOADING TEXT FILE WITH GENE NAME ALIASES
    if genenamesfile is None:
        genenamesfile = os.path.join(files_path,'Yeast_Protein_Names.txt')
    if not os.path.isfile(genenamesfile):
        raise ValueError('Following path does not exist: ' + genenamesfile)



#%% READ BAM FILE
    bam = pysam.AlignmentFile(bamfile, 'rb') #open bam formatted file for reading



#%% GET NAMES OF ALL CHROMOSOMES AS STORED IN THE BAM FILE
    ref_tid_dict = {} # 'I' | 0, 'II' | 1, ...
    ref_name_list = [] # 'I', 'II', ...
    for i in range(bam.nreferences): #if bam.nreferences does not work, use range(17) #16 chromosomes and the mitochondrial chromosome
        ref_name = bam.get_reference_name(i)
        ref_tid_dict[ref_name] = bam.get_tid(ref_name)
        ref_name_list.append(ref_name)



    del (ref_name, i)



#%% CONVERT CHROMOSOME NAMES IN DATA FILE TO ROMAN NUMERALS
    ref_romannums = chromosomename_roman_to_arabic()[0]
    ref_tid_roman_dict = {}
    for key,val in ref_tid_dict.items():
        ref_tid_roman_dict[ref_romannums[int(val)+1]] = key



    del (key, val, ref_romannums)



#%% GET SEQUENCE LENGTHS OF ALL CHROMOSOMES
    chr_length_dict = {} # 'I' | 230218, 'II' | 813184, ...
    chr_summedlength_dict = {} # 'I' | 0, 'II' | 230218, 'III' |  1043402, ...
    ref_summedlength = 0
    for key in ref_tid_dict:
        ref_length = bam.get_reference_length(key)
        chr_length_dict[key] = ref_length
        chr_summedlength_dict[key] = ref_summedlength
        ref_summedlength += ref_length



    del (key, ref_length, ref_summedlength)



#%% GET NUMBER OF MAPPED, UNMAPPED AND TOTAL AMOUNT OF READS PER CHROMOSOME
    # total_reads = bam.mapped
    stats = bam.get_index_statistics()
    chr_mappedreads_dict = {} # 'I' | [mapped, unmapped, total reads]
    for stat in stats:
        chr_mappedreads_dict[stat[0]] = [stat[1], stat[2], stat[3]]
        if stat[2] != 0:
            warnings.warn('Unmapped reads found in chromosome ' + stat[0])
    
    
    
    del (stat, stats)



#%% GET ALL READS WITHIN A SPECIFIED GENOMIC REGION
    tnnumber_dict = {}
    ll = 0 #Number of unique insertions in entire genome
    for kk in ref_name_list:
        timer_start = timeit.default_timer()
        read_counter = 0

        N_reads_kk = chr_mappedreads_dict[kk][2]
        start_array = np.empty(shape=(N_reads_kk), dtype=int)
        flag_array = np.empty(shape=(N_reads_kk), dtype=int)
        readlength_array = np.empty(shape=(N_reads_kk), dtype=int)
        
        
        #RETREIVING ALL THE READS FROM THE CURRENT CHROMOSOME.
        print('Getting reads for chromosme ', kk ,' ...')
        for reads in bam.fetch(kk, 0, chr_length_dict[kk], until_eof=True):
            read = str(reads).split('\t')

            start_array[read_counter] = int(read[3]) + 1
            flag_array[read_counter] = int(read[1])
            readlength_array[read_counter] = int(len(read[9]))
    
            read_counter += 1



        #CORRECT STARTING POSITION FOR READS WITH REVERSED ORIENTATION
        flag0coor_array = np.where(flag_array==0) #coordinates reads 5' -> 3'
        flag16coor_array = np.where(flag_array==16) # coordinates reads 3' -> 5'

        startdirect_array = start_array[flag0coor_array]
        flagdirect_array = flag_array[flag0coor_array]

        startindirect_array = start_array[flag16coor_array] + readlength_array[flag16coor_array]
        flagindirect_array = flag_array[flag16coor_array]

        start2_array = np.concatenate((startdirect_array, startindirect_array), axis=0)
        flag2_array = np.concatenate((flagdirect_array, flagindirect_array), axis=0)

        del (flag0coor_array, flag16coor_array, startdirect_array, flagdirect_array, startindirect_array, flagindirect_array)



        start2_sortindices = start2_array.argsort(kind='mergesort') #use mergesort for stable sorting
        start2_array = start2_array[start2_sortindices]
        flag2_array = flag2_array[start2_sortindices]

        del start2_sortindices



        #CREATE ARRAY OF START POSITION AND FLAGS OF ALL READS IN GENOME
        ref_tid_kk = int(ref_tid_dict[kk]+1)
        if ll == 0:
            tncoordinates_array = np.array([])

        mm = 0 # Number of unique reads per insertion
        jj = 1 # Number of unique reads in current chromosome (Number of transposons in current chromosome)
        for ii in range(1,len(start2_array)):
            if abs(start2_array[ii]-start2_array[ii-1]) <= 2 and flag2_array[ii] == flag2_array[ii-1]: #If two subsequent reads are within two basepairs and have the same orientation, add them together.
                mm += 1
            else:
                avg_start_pos = abs(round(np.mean(start2_array[ii-mm-1 : ii])))
                if tncoordinates_array.size == 0: #include first read
                    tncoordinates_array = np.array([ref_tid_kk, int(avg_start_pos), int(flag2_array[ii-1])])
                    readnumb_list = [mm+1]
                else:
                    tncoordinates_array = np.vstack((tncoordinates_array, [ref_tid_kk, int(avg_start_pos), int(flag2_array[ii-1])]))    
                    readnumb_list.append(mm+1)
                mm = 0
                jj += 1
                ll += 1

            if ii == len(start2_array) - 1: #include last read
                avg_start_pos = abs(round(np.mean(start2_array[ii-mm-1 : ii])))
                tncoordinates_array = np.vstack((tncoordinates_array, [ref_tid_kk, int(avg_start_pos), int(flag2_array[ii-1])]))
                readnumb_list.append(mm+1)

        tnnumber_dict[kk] = jj

        del (jj, start_array, flag_array, readlength_array, flag2_array, start2_array, ref_tid_kk)



        timer_end = timeit.default_timer()
        print('Chromosome ', kk, ' completed in ', (timer_end - timer_start), ' seconds.')



    readnumb_array = np.array(readnumb_list)
    del readnumb_list

    tncoordinatescopy_array = np.array(tncoordinates_array, copy=True)




#%% GET LIST OF ALL GENES AND ALL ESSENTIAL GENES
    print('Getting coordinates of all genes ...')

    # GET POSITION GENES
    gff_path = os.path.join(files_path,'Saccharomyces_cerevisiae.R64-1-1.99.gff3')
    genecoordinates_dict = gene_position(gff_path) #'YAL069W' | ['I', 335, 649], ...



    # GET ALL ANNOTATED ESSENTIAL GENES
    essential_path1 = os.path.join(files_path,'Cervisiae_EssentialGenes_List_1.txt')
    essential_path2 = os.path.join(files_path,'Cervisiae_EssentialGenes_List_2.txt')
    essentialnames_list = list_known_essentials([essential_path1, essential_path2])


    
    # GET ALIASES OF ALL GENES
    names_path = os.path.join(files_path,'Yeast_Protein_Names.txt')
    aliases_designation_dict = gene_aliases(names_path)[0] #'YMR056C' | ['AAC1'], ...


    
    # FOR ALL GENE IN GENECOORDINATES_DICT, CHECK IF THEY ARE ANNOTATED AS ESSENTIAL
    essentialcoordinates_dict = {}
    for gene in genecoordinates_dict:
        if gene in essentialnames_list:
            essentialcoordinates_dict[gene] = genecoordinates_dict.get(gene).copy()
        else:
            gene_aliases_list = []
            for key, val in aliases_designation_dict.items():
                if gene == key or gene in val: #If gene occurs as key or in the values list in aliases_designation_dict, put all its aliases in a single list.
                    gene_aliases_list.append(key)
                    for aliases in aliases_designation_dict.get(key):
                        gene_aliases_list.append(aliases)
    
            for gene_alias in gene_aliases_list:
                if gene_alias in essentialnames_list:
                    essentialcoordinates_dict[gene_alias] = genecoordinates_dict.get(gene).copy()
                    break



    del (gff_path, essential_path1, essential_path2, names_path, gene, gene_alias, essentialnames_list, gene_aliases_list, key, val)



#%% CONCATENATE ALL CHROMOSOMES

    #FOR EACH INSERTION LOCATION, ADD THE LENGTH OF ALL PREVIOUS CHROMOSOMES.
    ll = 0
    for ii in range(1,len(ref_name_list)):
        ll += chr_length_dict[ref_name_list[ii-1]]
        aa = np.where(tncoordinatescopy_array[:,0] == ii + 1)
        tncoordinatescopy_array[aa,1] = tncoordinatescopy_array[aa,1] + ll



    #FOR EACH GENE LOCATION, ADD THE LENGTH OF ALL PREVIOUS CHROMOSOMES.
    for key in genecoordinates_dict:
        gene_chrom = ref_tid_roman_dict.get(genecoordinates_dict.get(key)[0])
        genecoordinates_dict[key][1] = genecoordinates_dict.get(key)[1] + chr_summedlength_dict.get(gene_chrom)
        genecoordinates_dict[key][2] = genecoordinates_dict.get(key)[2] + chr_summedlength_dict.get(gene_chrom)



    #FOR EACH ESSENTIAL GENE LOCATION, ADD THE LENGTH OF ALL PREVIOUS CHROMOSOMES.
    for key in essentialcoordinates_dict:
        gene_chrom = ref_tid_roman_dict.get(essentialcoordinates_dict.get(key)[0])
        essentialcoordinates_dict[key][1] = essentialcoordinates_dict.get(key)[1] + chr_summedlength_dict.get(gene_chrom)
        essentialcoordinates_dict[key][2] = essentialcoordinates_dict.get(key)[2] + chr_summedlength_dict.get(gene_chrom)



    del (ii, ll, aa, key, gene_chrom)



#%% GET NUMBER OF TRANSPOSONS AND READS PER GENE
    print('Get number of insertions and reads per gene ...')
    
    #ALL GENES
    tnpergene_dict = {}
    readpergene_dict = {}
    # readpergenecrude_dict = {}
    for gene in genecoordinates_dict:
        xx = np.where(np.logical_and(tncoordinatescopy_array[:,1] >= genecoordinates_dict.get(gene)[1], tncoordinatescopy_array[:,1] <= genecoordinates_dict.get(gene)[2])) #get all insertions within range of current gene
        tnpergene_dict[gene] = np.size(xx)
        readpergene_dict[gene] = sum(readnumb_array[xx]) - max(readnumb_array[xx], default=0)
        # readpergenecrude_dict[gene] = sum(readnumb_array[xx])



    #ONLY ESSENTIAL GENES
    tnperessential_dict = {}
    readperessential_dict = {}
    # readperessentialcrude_dict = {}
    for gene in essentialcoordinates_dict:
        xx = np.where(np.logical_and(tncoordinatescopy_array[:,1] >= essentialcoordinates_dict.get(gene)[1], tncoordinatescopy_array[:,1] <= essentialcoordinates_dict.get(gene)[2]))
        tnperessential_dict[gene] = np.size(xx)
        readperessential_dict[gene] = sum(readnumb_array[xx]) - max(readnumb_array[xx], default=0)
        # readperessentialcrude_dict[gene] = sum(readnumb_array[xx])



    del (xx, gene)



#%% CREATE BED FILE
    bedfile = bamfile+'.bed'
    print('Writing bed file at: ', bedfile)


    with open(bedfile, 'w') as f:
        
        f.write('track name=' + filename + ' useScore=1\n')
        
        coordinates_counter = 0
        for tn in tncoordinates_array:
            refname = [key for key, val in ref_tid_dict.items() if val == tn[0] - 1][0]
            if refname == 'Mito':
                refname = 'M'
            f.write('chr' + refname + ' ' + str(tn[1]) + ' ' + str(tn[1] + 1) + ' . ' + str(100+readnumb_array[coordinates_counter]*20) + '\n')
            coordinates_counter += 1



    del (bedfile, coordinates_counter, refname)

#%% CREATE TEXT FILE WITH TRANSPOSONS AND READS PER GENE
    pergenefile = bamfile+'_pergene.txt'
    print('Writing pergene.txt file at: ', pergenefile)


    with open(pergenefile, 'w') as f:

        f.write('Gene name\tNumber of transposons per gene\tNumber of reads per gene\n')

        for gene in tnpergene_dict:
            tnpergene = tnpergene_dict[gene]
            readpergene = readpergene_dict[gene]
            if gene in aliases_designation_dict:
                gene_alias = aliases_designation_dict.get(gene)[0]
            else:
                gene_alias = gene
            f.write(gene_alias + '\t' + str(tnpergene) + '\t' + str(readpergene) + '\n')



    del (pergenefile, gene, gene_alias, tnpergene, readpergene)

#%% CREATE TEXT FILE TRANSPOONS AND READS PER ESSENTIAL GENE
    peressentialfile = bamfile+'_peressential.txt'
    print('Writing peressential.txt file at: ',peressentialfile)
    
    
    with open(peressentialfile, 'w') as f:
        
        f.write('Gene name\tNumber of transposons per gene\tNumber of reads per gene\n')
        
        for essential in tnperessential_dict:
            tnperessential = tnperessential_dict[essential]
            readperessential = readperessential_dict[essential]
            if essential in aliases_designation_dict:
                essential_alias = aliases_designation_dict.get(essential)[0]
            else:
                essential_alias = essential
            f.write(essential_alias + '\t' + str(tnperessential) + '\t' + str(readperessential) + '\n')



    del (peressentialfile, essential, essential_alias, tnperessential, readperessential)
                
#%% ADD INSERTIONS AT SAME LOCATION BUT WITH DIFFERENT ORIENTATIONS TOGETHER (FOR STORING IN WIG-FILE)
    wigfile = bamfile+'.wig'
    print('Writing wig file at: ', wigfile)



    readnumbwig_array = readnumb_array.copy()



    unique_index_array = np.array([], dtype=int) #=cc
    N_uniques_perchr_list = []
    ll = 0
    for kk in ref_name_list:
        index = np.where(tncoordinates_array[:,0] == int(ref_tid_dict[kk]+1)) #get indices for current chromosome.
        unique_index = np.unique(tncoordinates_array[index][:,1], return_index=True)[1] #get all insertion locations (in tncoordinates, all rows, column 1)

        unique_index_array = np.append(unique_index_array, (unique_index+ll), axis=0)

        ll += np.count_nonzero(tncoordinates_array[:,0] == int(ref_tid_dict[kk]+1))
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

#%% CREATING WIG FILE
    with  open(wigfile, 'w') as f:
        f.write('track type=wiggle_0 ,maxheightPixels=60 name='+filename+'\n')
        for kk in ref_name_list:
            f.write('VariableStep chrom=chr' + kk + '\n')

            index = np.where(tncoordinateswig_duplicatesremoved_array[:,0] == int(ref_tid_dict[kk]+1)) #get indices for current chromosome.
            for ii in index[0]:
                f.write(str(tncoordinateswig_duplicatesremoved_array[ii][1]) + ' ' + str(readnumbwig_duplicatesremoved_array[ii]) + '\n')


    del (wigfile, kk, ii, index)

#%%
if __name__ == '__main__':
    transposonmapper()


