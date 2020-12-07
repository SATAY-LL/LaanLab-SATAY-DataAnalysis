# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 10:44:55 2020

@author: gregoryvanbeek

This file reads a wig_file and creates a folder at the same directory where the
wigfile is seperated in wigfiles for individual chromosomes (including Mito)

REQUIRES chromosome_names_in_files (https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/python_modules/chromosome_names_in_files.py)
"""

#inputfile=r"C:\Users\gregoryvanbeek\Documents\Data_Sets\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig"


import os, sys
dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(dirname,'python_modules'))
from chromosome_names_in_files import chromosome_name_wigfile

def split_wigfile(inputfile='', verbose=True):



    if not os.path.isfile(inputfile) and not inputfile == '':
        print('WARNING: inputfile does not exists')
        exit()



    filepath = os.path.dirname(inputfile)
    filename = os.path.splitext(os.path.basename(inputfile))[0]

    directoryname = os.path.join(filepath, filename + '_chromosomesplit')
    
    if not os.path.exists(directoryname):
        os.mkdir(directoryname)



    chromosome_names = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','M']



    chrom_names_dict, chrom_start_line_dict, chrom_end_line_dict = chromosome_name_wigfile(inputfile)



    with open(inputfile, 'r') as f:
        lines = f.readlines()
    header = lines[0]

    for chrom in chromosome_names:

        outputfile = os.path.join(directoryname, filename + '_' + str(chrom) + '.wig')


        with open(outputfile, 'w+') as f:
            f.write(header)
            for l in range(chrom_start_line_dict.get(chrom)-1, chrom_end_line_dict.get(chrom)+1):
                f.write(lines[l])



if __name__ == '__main__':
    split_wigfile(r"C:\Users\gregoryvanbeek\Documents\Data_Sets\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig")