######THIS FILE LOADS A TEXT FILE WITH INFORMATION ABOUT THE NUMBER OF READS AND TRANSPOSONS PER GENE.######
######PLOTS A HISTOGRAM OR VIOLINPLOT FOR THE NUMBER OF READS AND TRANSPOSONS DIFFERENTIATED PER ESSENTIALITY.######

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def stats_pergene(normalize,filepath, filename):

    #GET ALL KNOWN ESSENTIAL GENES FROM TWO DIFFERENT FILES AND COMBINE THEM IN ONE LIST
    essential_genes_files = ['/Users/gregory/Documents/LaanLab/EssentialGenes_Database/Cervisiae_EssentialGenes_List_1.txt',
                            '/Users/gregory/Documents/LaanLab/EssentialGenes_Database/Cervisiae_EssentialGenes_List_2.txt']
    known_essential_gene_list = []

    for files in essential_genes_files:
        with open(files) as f:
            for header_lines in range(3):
                next(f)
            for lines in f:
                known_essential_gene_list.append(lines.rstrip('\n'))




    #CREATE A DICTIONARY WITH ALL GENES (BOTH THE COMMON NAME AND THEIR SYSTEMATIC NAME) AND SAVE THEM WITH THEIR RESPECTIVE LENGHTS (IN TERMS OF BP WHICH IS DEFINED AS bp=aa*3)
    gene_information_file = '/Users/gregory/Documents/LaanLab/EssentialGenes_Database/Yeast_Protein_Names.txt'
    gene_length_dict = {}
    with open(gene_information_file) as f:
        lines = f.readlines()
        for i in range(58,len(lines)-6):
            n=0
            l = lines[i]
            extra_columns = l.count(';')    #count how many times ';' occurs in a line. This is needed to get the right columns as sometimes aliases of gene names are presented in extra columns
            l_short = ' '.join(l.split())
            l_list = l_short.split(' ')
            if l_list[1+extra_columns] == 'GAG' or l_list[1+extra_columns] == 'POL':    #these are sequences that sometimes occurs and have to be ignored
                extra_columns = extra_columns + 1
            gene_length_dict[l_list[0].strip(';')] = int(l_list[5+extra_columns])*3
            gene_length_dict[l_list[1+extra_columns]] = int(l_list[5+extra_columns])*3
            if extra_columns > 0:
                for n in range(extra_columns+1):
                    gene_length_dict[l_list[0+n]] = int(l_list[5+extra_columns])*3




    gene_name_list = []
    reads_pergene_list = []
    tn_pergene_list = []
    gene_essential_boolean_list = []
    reads_pergene_norm_list = []
    tn_pergene_norm_list = []

    file = os.path.join(filepath,filename)

    with open(file) as f:
        next(f)
        gene_counter = 0
        fail_counter = 0
        for lines in f:
            split_line = lines.split('\t')

            gene_name = split_line[0].rstrip()



            if len(split_line) > 1:
                tn_pergene = int(split_line[1])
            else:
                tn_pergene = np.nan

            if len(split_line) > 2:
                reads_pergene = int(split_line[2])
            else:
                reads_pergene = np.nan

            reads_pergene_list.append(reads_pergene)
            tn_pergene_list.append(tn_pergene)


            if normalize == 'False':
                if gene_name in known_essential_gene_list:
                    gene_essential_boolean_list.append('True') 
                else:
                    gene_essential_boolean_list.append('False')
                gene_name_list.append(gene_name)

            else:
                gene_length = gene_length_dict.get(gene_name)
                try:
                    reads_pergene_norm_list.append(reads_pergene/gene_length)
                    tn_pergene_norm_list.append(tn_pergene/gene_length)

                    if gene_name in known_essential_gene_list:
                        gene_essential_boolean_list.append('True') 
                    else:
                        gene_essential_boolean_list.append('False')

                    gene_name_list.append(gene_name)
                    gene_counter = gene_counter + 1
                except:
                    fail_counter = fail_counter + 1

    print('Number of genes analyzed: ',gene_counter)
    print('Number of genes not found: ',fail_counter)
    print('')
    print('')




    #CREATE DATAFRAME
    if normalize == 'False':
        genes = {'Gene_name': gene_name_list,
                'Transposons_per_gene': tn_pergene_list,
                'Reads_per_gene': reads_pergene_list,
                'Essential_gene': gene_essential_boolean_list
                }
        df = pd.DataFrame(genes,columns=['Gene_name','Transposons_per_gene','Reads_per_gene', 'Essential_gene'])
    else:
        genes = {'Gene_name': gene_name_list,
                'Transposon_density_per_gene': tn_pergene_norm_list,
                'Read_density_per_gene': reads_pergene_norm_list,
                'Essential_gene': gene_essential_boolean_list
                }
        df = pd.DataFrame(genes,columns=['Gene_name','Transposon_density_per_gene','Read_density_per_gene', 'Essential_gene'])




    print('Average number of reads all genes: ',np.nanmean(reads_pergene_list))
    print('Median number of reads all genes: ',np.nanmedian(reads_pergene_list))
    print('Standard Deviation number of reads all genes: ',np.nanstd(reads_pergene_list))
    print('')

    print('Average number of transposons all genes: ',np.nanmean(tn_pergene_list))
    print('Median number of transposons all genes: ',np.nanmedian(tn_pergene_list))
    print('Standard Deviation number of transposons all genes: ',np.nanstd(tn_pergene_list))
    print('')
    print('')



    indices_essential_genes_list = [i for i, x in enumerate(gene_essential_boolean_list) if x == 'True']
    reads_peressentialgene_list = [reads_pergene_list[i] for i in indices_essential_genes_list]
    print('Average number of reads of essential genes: ',np.nanmean(reads_peressentialgene_list))
    print('Median number of reads of essential genes: ',np.nanmedian(reads_peressentialgene_list))
    print('Maximum number of reads of essential genes: ',np.nanmax(reads_peressentialgene_list))
    print('Minimum number of reads of essential genes: ',np.nanmin(reads_peressentialgene_list))
    print('')

    indices_nonessential_genes_list = [i for i, x in enumerate(gene_essential_boolean_list) if x == 'False']
    reads_pernonessentialgene_list = [reads_pergene_list[i] for i in indices_nonessential_genes_list]
    print('Average number of reads of nonessential genes: ',np.nanmean(reads_pernonessentialgene_list))
    print('Median number of reads of nonessential genes: ',np.nanmedian(reads_pernonessentialgene_list))
    print('Maximum number of reads of nonessential genes: ',np.nanmax(reads_pernonessentialgene_list))
    print('Minimum number of reads of nonessential genes: ',np.nanmin(reads_pernonessentialgene_list))
    print('')
    print('')



    tn_peressentialgene_list = [tn_pergene_list[i] for i in indices_essential_genes_list]
    print('Average number of transposons of essential genes: ',np.nanmean(tn_peressentialgene_list))
    print('Median number of transposons of essential genes: ',np.nanmedian(tn_peressentialgene_list))
    print('Standard deviation number of transposons of essential genes: ',np.nanstd(tn_peressentialgene_list))
    print('Maximum number of transposons of essential genes: ',np.nanmax(tn_peressentialgene_list))
    print('Minimum number of transposons of essential genes: ',np.nanmin(tn_peressentialgene_list))
    print('')

    tn_pernonessentialgene_list = [tn_pergene_list[i] for i in indices_nonessential_genes_list]
    print('Average number of transposons of nonessential genes: ',np.nanmean(tn_pernonessentialgene_list))
    print('Median number of transposons of nonessential genes: ',np.nanmedian(tn_pernonessentialgene_list))
    print('Standard deviation number of transposons of nonessential genes: ',np.nanstd(tn_pernonessentialgene_list))
    print('Maximum number of transposons of nonessential genes: ',np.nanmax(tn_pernonessentialgene_list))
    print('Minimum number of transposons of nonessential genes: ',np.nanmin(tn_pernonessentialgene_list))
    print('')




    fig, (ax1,ax2) = plt.subplots(1,2)
    n_reads,bins_reads,patches_reads = ax1.hist(reads_pergene_list, density=False, bins=range(min(reads_pergene_list), max(reads_pergene_list) + 1000, 1000), facecolor='blue',alpha=0.5, label='Reads')
    ax1.axvline(x=np.nanmedian(reads_pergene_list), linestyle='--', label='Median Reads', color='k')
    ax1.set_xlabel('Reads per gene (all genes)')
    ax1.set_ylabel('Occurance')
    ax1.legend()
    ax1.grid()
    n_tn,bins_tn,patches_tn = ax2.hist(tn_pergene_list, density=False, bins=range(min(tn_pergene_list), max(tn_pergene_list) + 50, 50), facecolor='red',alpha=0.5, label='Transposons')
    ax2.axvline(x=np.nanmedian(tn_pergene_list), linestyle='--', label='Median Transposons', color='k')
    ax2.set_xlabel('Counts per gene (all genes)')
    ax2.set_ylabel('Occurance')
    ax2.legend()
    ax2.grid()
    plt.show()




    fig, (ax1, ax2) = plt.subplots(1,2)
    sns.set(style='whitegrid', palette='pastel', color_codes=True)
    df['Reads'] = ''
    if normalize == 'False':
        sns.violinplot(x='Reads',y="Reads_per_gene", hue='Essential_gene', inner='quarter', gridsize=400, palette={'False':'r', 'True':'g'}, split=True, cut=0, orient='v', data=df, ax=ax1)
    else:
        sns.violinplot(x='Reads',y="Read_density_per_gene", hue='Essential_gene', inner='quarter', gridsize=400, palette={'False':'r', 'True':'g'}, split=True, cut=0, orient='v', data=df, ax=ax1)
    ax1.grid()

    df['Transposons'] = ''
    if normalize == 'False':
        sns.violinplot(x='Transposons',y="Transposons_per_gene", hue='Essential_gene', inner='quarter', gridsize=400, palette={'False':'r', 'True':'g'}, split=True, cut=0, orient='v', data=df, ax=ax2)
    else:
        sns.violinplot(x='Transposons',y="Transposon_density_per_gene", hue='Essential_gene', inner='quarter', gridsize=400, palette={'False':'r', 'True':'g'}, split=True, cut=0, orient='v', data=df, ax=ax2)
    ax2.grid()
    plt.show()




if __name__ == '__main__':
    stats_pergene('False','/Users/gregory/Documents/LaanLab/LaanLab_Data/Michel2017_WT1/','E-MTAB-4885.WT1.bam_pergene.txt')#'Cerevisiae_WT1_Michel2017_Trimmed_Aligned.sorted.bam_pergene.txt')#
