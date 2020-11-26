# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 07:43:46 2020

@author: gregoryvanbeek

This is a code to test the output of demulitplexing.py
"""

input_file = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1_sample1_forward.fq"
#input_file = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1_sample1_reverse.fq"
#input_file = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1_sample2_forward.fq"
#input_file = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1_sample2_reverse.fq"
#input_file = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1_unmatched.fq"

#barcode_list = ["GCCACATA",
#                "GCGAGTAA",
#                "GAGCTGAA",
#                "GATAGACA"]

f = open(input_file,'r')


barcode_counter_dict = {"GCCACATA": 0,
                        "GCGAGTAA": 0,
                        "GAGCTGAA": 0,
                        "GATAGACA": 0}

line_counter = 0
for line in f:
    if line.startswith('@'):
        line_counter += 1
        line_seq = f.readline()
        for barcode in barcode_counter_dict:
            if barcode in line_seq:
                barcode_counter_dict[barcode] += 1
                break
            else:
                pass
                
f.close()

print("Total number of reads found: %i" % line_counter)
print(barcode_counter_dict)

