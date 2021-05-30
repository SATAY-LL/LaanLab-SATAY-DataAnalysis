# -*- coding: utf-8 -*-
"""
"""

import os
import numpy as np


#%%INPUT

file = r""
verbose = True

#%%
def total_mapped_reads(file, verbose=True):
    '''
    This function gets the total number of reads that are present in the dataset.
    Input is either a .bed-file or a .wig-file.
    Output is an dict containing the median read count per insertion, number of insertions and number of reads.
    '''
    if not os.path.isfile(file):
        print('WARNING: Following path not found: ', file)
        return()

    else:

        extension = os.path.splitext(file)



        if extension[1] == '.wig':
            
            with open(file, 'r') as f:
                lines = f.readlines()

            total_reads = 0
            total_ins = 0
            reads_list = []
            chrom_line = 'variablestep'
            for line in lines[1:]:
                if not line.lower().startswith(chrom_line):
                    line_split = line.strip('\n').split(' ')
                    reads_list.append(int(line_split[1]))
                    total_reads += int(line_split[1])
                    total_ins += 1

            if verbose == True:
                print('Number of mapped reads found = %i' % total_reads)
                print('Number of insertions found = %i' % total_ins)
                print('Median read per transposon = %i' % np.median(reads_list))

            count_dict = {"Ninsertions": total_ins,
                          "Nreads": total_reads,
                          "Median": np.median(reads_list)}

            return(count_dict)



        elif extension[1] == '.bed':

            with open(file,'r') as f:
                lines = f.readlines()

            total_reads = 0
            total_ins = 0
            reads_list = []
            for line in lines[1:]:
                line_split = line.strip('\n').split(' ')
                reads_list.append((int(line_split[4]) - 100)/20)
                total_reads += (int(line_split[4]) - 100)/20
                total_ins += 1

            if verbose == True:
                print('Number of mapped reads found = %i' % total_reads)
                print('Number of insertions found = %i' % total_ins)
                print('Median read per transposon = %i' % np.median(reads_list))

            count_dict = {"Ninsertions": total_ins,
                          "Nreads": total_reads,
                          "Median": np.median(reads_list)}

            return(count_dict)


        else:
            print('WARNING: No valid file extension. Files must have either of the following extension; .wig or .bed')



#%%
if __name__ == '__main__':
    count_dict = total_mapped_reads(file=file, verbose=verbose)