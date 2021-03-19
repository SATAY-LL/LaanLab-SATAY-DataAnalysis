# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 17:04:02 2021

@author: gregoryvanbeek
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'..','python_modules'))
from chromosome_and_gene_positions import chromosome_position, chromosomename_roman_to_arabic, gene_position
from essential_genes_names import list_known_essentials
from gene_names import gene_aliases
from chromosome_names_in_files import chromosome_name_bedfile


#%%INPUT
bed_files=[r"",
           r""]
variable = "insertions" # "insertions", "reads"
chromosome='i'#['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']
set_barwidth=None
set_logscale = False
savefig=False


#%%
def compareplot(bed_files=None, variable="insertions", chromosome=None, set_barwidth=None, set_logscale=False, savefig=False):
    '''This function creates a bar plot along a specified chromosome for the number of transposons.
    The height of each bar represents the number of transposons at the genomic position indicated on the x-axis.
    The input is as follows:
        -The bed-files ('bed_files', a list containing two paths, each refering to a bed-file [mandatory]),
        -Which chromosome ('chromosome', indicated by roman numeral or list of roman numerals [optional]),
        -The width of the bars ('bar_width-user_set', indicated by an integer [optional]),
        -Path to where to save the figures ('savefigure_path', string containing an existing path [optional]),
        -Name of the figures ('savefigure_name', string containing a single name, the name will be automatically extended with the chromosomal number [optional]).
    
    The bed_file is one of the files created by the Matlab code from the kornmann-lab.
    The figure shows two graphs, the top one represents the first bed-file given in the list, the bottom plot the second bed-file in the list.
    If the chromosome number is not set by the user, it automatically loops over all chromosomes and determines the figures for each of them.    
    The bar_width determines how many basepairs are put in one bin. Little basepairs per bin may be slow. Too many basepairs in one bin and possible low transposon areas might be obscured.
    When either the savefigure_path and/or the savefigure_name is left empty, the figure won't be saved.
    If the both these variables are given, the figures are saved using the path/figurename_chromX where the _chromX extension is automatically added.
    
    The background of the graph is color coded to indicate areas that code for genes.
    For this a list for essential genes is needed (used in 'list_known_essentials' function) and a .gff file is required (for the functions in 'chromosome_and_gene_positions.py') and a list for gene aliases (used in the function 'gene_aliases').
    '''
#%% USED FILES
    gff_file = os.path.join(file_dirname,'..','data_files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
    essential_genes_files = [os.path.join(file_dirname,'..','data_files','Cerevisiae_EssentialGenes_List_1.txt'),
                            os.path.join(file_dirname,'..','data_files','Cerevisiae_EssentialGenes_List_2.txt')]
    gene_information_file = os.path.join(file_dirname,'..','data_files','Yeast_Protein_Names.txt')
#%% GET CHROMOSOME LENGTHS AND POSITIONS
    chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file)
    
#%% GET ALL GENES IN CURRENT CHROMOSOME
    gene_pos_dict = gene_position(gff_file)
    genes_essential_list = list_known_essentials(essential_genes_files, verbose=False)
    gene_alias_list = gene_aliases(gene_information_file)[0]
    
#%% DETERMINE WHICH CHROMOSOME NEEDS TO BE ANALYZED AND LOOP OVER THE CHROMOSOMES
    if type(chromosome) is list:
        chrom_list = chromosome
    elif type(chromosome) is str:
        chrom_list = [chromosome.upper()]
    else:
        chrom_list = []
        roman_to_arabic_numerals = chromosomename_roman_to_arabic()[1]
        for keys in roman_to_arabic_numerals:
            chrom_list.append(keys)
    
    for chrom in chrom_list:
        print('')
        print('Analyzing chromosome: ', chrom)
        genes_currentchrom_pos_list = [k for k, v in gene_pos_dict.items() if chrom in v]

#%% READ BED FILE
        allinsertionsites_allfiles_list = []
        alltransposoncounts_allfiles_binnedlist = []
        for bed_file in bed_files:
            print("Processing file: %s" % bed_file)
            with open(bed_file) as f:
                lines = f.readlines()
    
#%% GET NAMES FOR THE CHROMOSOMES IN THE BED FILE
            chrom_start_index_dict, chrom_end_index_dict= chromosome_name_bedfile(bed_file)[1:3]

#%% GET ALL TRANSPOSON COUNTS
            allcounts_list = np.zeros(chr_length_dict.get(chrom)+2)
            if variable == "insertions":
                for line in lines[chrom_start_index_dict.get(chrom):chrom_end_index_dict.get(chrom)+1]:
                    line = line.strip('\n').split()
                    allcounts_list[int(line[1])] += 1

            elif variable == "reads":
                for line in lines[chrom_start_index_dict.get(chrom):chrom_end_index_dict.get(chrom)+1]:
                        line = line.strip('\n').split()
                        allcounts_list[int(line[1])] += int(line[4])

#%% BINNING OF THE READS
            if set_barwidth == None:
                bar_width = int(chr_length_dict.get(chrom)/500)
            else:
                bar_width = set_barwidth
    
            allcounts_binnedlist = []
            val_counter = 0
            sum_values = 0
            if bar_width == 1:
                allcounts_binnedlist = allcounts_list
                allinsertionsites_list = np.linspace(0,chr_length_dict.get(chrom),int(chr_length_dict.get(chrom)/float(bar_width)))
            else:
                for n in range(len(allcounts_list)):
                    if val_counter % bar_width != 0:
                        sum_values += allcounts_list[n]
                    elif val_counter % bar_width == 0:
                        allcounts_binnedlist.append(sum_values)
                        sum_values = 0
                    val_counter += 1
                    
                allinsertionsites_list = np.linspace(0,chr_length_dict.get(chrom),int(chr_length_dict.get(chrom)/bar_width)+1)
        
            allinsertionsites_allfiles_list.append(allinsertionsites_list)
            alltransposoncounts_allfiles_binnedlist.append(allcounts_binnedlist)

#%% DETERMINE DIFFERENCE BETWEEN DATASETS TRANSPOSONCOUNTS
        transposoncounts_positivedifference_list = [0]*len(alltransposoncounts_allfiles_binnedlist[0])
        transposoncounts_negativedifference_list = [0]*len(alltransposoncounts_allfiles_binnedlist[0])
        for i in range(0,len(alltransposoncounts_allfiles_binnedlist[0])):
            difference = alltransposoncounts_allfiles_binnedlist[0][i]-alltransposoncounts_allfiles_binnedlist[1][i]
            if difference >= 0:
                transposoncounts_positivedifference_list[i] = difference
            elif difference < 0:
                transposoncounts_negativedifference_list[i] = -difference

#%% PLOTTING
        print('Plotting chromosome ', chrom, '...')
        print('bar width for plotting is ',bar_width)
        binsize = bar_width
        font_size = 12
        max_ylim = max([item for sublist in alltransposoncounts_allfiles_binnedlist for item in sublist]) #GET MAXIMUM VALUE FOR SETTING THE Y AXIS LIMIT EQUAL FOR BOTH GRAPHS
        max_ylim = max_ylim + 0.1*max_ylim
        
        
        plt.figure(figsize=(19,9))
        grid = plt.GridSpec(2, 1, wspace=0.0, hspace=0.0)
    
    
        ax1 = plt.subplot(grid[0,0])
        for gene in genes_currentchrom_pos_list:
            gene_start_pos = int(gene_pos_dict.get(gene)[1])
            gene_end_pos = int(gene_pos_dict.get(gene)[2])
            if gene in genes_essential_list:
                ax1.axvspan(gene_start_pos,gene_end_pos,facecolor='g',alpha=0.3)
                ax1.text(gene_start_pos,max_ylim,gene_alias_list.get(gene)[0], rotation=45)
            else:
                ax1.axvspan(gene_start_pos,gene_end_pos,facecolor='r',alpha=0.3)
    
        ax1.bar(allinsertionsites_allfiles_list[0],alltransposoncounts_allfiles_binnedlist[0],width=binsize,color=(0.2,0.2,0.2,0.8))
        ax1.bar(allinsertionsites_allfiles_list[0],transposoncounts_positivedifference_list,width=binsize,color=(0.52,0.71,0.90,0.8))

        if set_logscale == True:
            ax1.set_yscale('log')
        else:
            ax1.set_ylim(0,max_ylim)
        ax1.set_axisbelow(True)
        ax1.grid(True)
        if variable == "insertions":
            ax1.set_ylabel('Aboslute insertion count', fontsize=font_size)
        elif variable == "reads":
            ax1.set_ylabel('Aboslute read count', fontsize=font_size)
        ax1.set_xlim(0,chr_length_dict.get(chrom))
    
    
        ax2 = plt.subplot(grid[1,0])
        for gene in genes_currentchrom_pos_list:
            gene_start_pos = int(gene_pos_dict.get(gene)[1])
            gene_end_pos = int(gene_pos_dict.get(gene)[2])
            if gene in genes_essential_list:
                ax2.axvspan(gene_start_pos,gene_end_pos,facecolor='g',alpha=0.3)
            else:
                ax2.axvspan(gene_start_pos,gene_end_pos,facecolor='r',alpha=0.3)

        if variable == "insertions":
            ax2.bar(allinsertionsites_allfiles_list[1],alltransposoncounts_allfiles_binnedlist[1],width=binsize,color=(0.2,0.2,0.2,0.8), label='Number of transposons')
        elif variable == "reads":
            ax2.bar(allinsertionsites_allfiles_list[1],alltransposoncounts_allfiles_binnedlist[1],width=binsize,color=(0.2,0.2,0.2,0.8), label='Number of reads')
        ax2.bar(allinsertionsites_allfiles_list[1],transposoncounts_negativedifference_list,width=binsize,color=(0.52,0.71,0.90,0.8), label='Absolute difference datasets (set1-set2)')

        if set_logscale == True:
            ax2.set_yscale('log')
        else:
            ax2.set_ylim(0,max_ylim)
        ax2.set_axisbelow(True)
        ax2.grid(True)
        if variable == "insertions":
            ax2.set_ylabel('Aboslute insertion count', fontsize=font_size)
        elif variable == "reads":
            ax2.set_ylabel('Aboslute read count', fontsize=font_size)
        ax2.set_xlabel('Basepair position on chromosome '+chrom, fontsize=font_size)
        ax2.set_xlim(0,chr_length_dict.get(chrom))
        ax2.invert_yaxis()
        ax2.legend(loc='lower left', fontsize=font_size)
        
        plt.tight_layout()
        
        
        if savefig == True:
            saving_name = os.path.join(os.path.dirname(bed_files[0]), os.path.basename(bed_files[0]).strip(".bed")+"_compareplot_chrom" + chrom + ".png")
            plt.savefig(saving_name)
            plt.close()



#%%
if __name__ == '__main__':
    compareplot(bed_files=bed_files,
                variable=variable,
                chromosome=chromosome,
                set_barwidth=set_barwidth,
                set_logscale=set_logscale,
                savefig=savefig)