# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a python script for demultiplexing reads from transposon sequencing data.
"""

#import os, sys


file1 = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq"
file2 = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"



barcode_forward_sample1 = "GCCACATA"
barcode_reverse_sample1 = "GCGAGTAA"
barcode_forward_sample2 = "GAGCTGAA"
barcode_reverse_sample2 = "GATAGACA"

barcode_dict = {barcode_forward_sample1: 'Sample1_for',
                barcode_reverse_sample1: 'Sample1_rev',
                barcode_forward_sample2: 'Sample2_for',
                barcode_reverse_sample2: 'Sample2_rev'}



f1 = open(file1)
f2 = open(file2)


for line1 in f1:
    line2 = f2.readline()
    if line1.startswith('@') and line2.startswith('@') and line1.split(' ')[0] == line2.split(' ')[0]:
        line1 = f1.readline()
        bc1 = [barcode for barcode in barcode_dict if barcode in line1]
    
        line2 = f2.readline()
        bc2 = [barcode for barcode in barcode_dict if barcode in line1]


f1.close()
f2.close()