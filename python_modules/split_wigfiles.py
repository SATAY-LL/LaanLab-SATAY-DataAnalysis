# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 10:44:55 2020

@author: gregoryvanbeek

This file reads a wig_file and creates a folder at the same directory where the
wigfile is seperated in wigfiles for individual chromosomes (including Mito)
"""

inputfile=r"C:\Users\gregoryvanbeek\Documents\Data_Sets\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig"


import os

def split_wigfile(inputfile='', verbose=True):



    if not os.path.isfile(inputfile) and not inputfile == '':
        print('WARNING: inputfile does not exists')
        exit()



    filepath = os.path.dirname(inputfile)
    filename = os.path.basename(inputfile)
    
    directoryname = os.path.join(filepath, os.path.splitext(filename)[0] + '_chromosomesplit')
    
    if not os.path.exists(directoryname):
        os.mkdir(directoryname)



    chromosome_names = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','Mito']



    with open(inputfile, 'r') as f:
        lines = f.readlines()





if __name__ == '__main__':
    split_wigfile(r"C:\Users\gregoryvanbeek\Documents\Data_Sets\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig")