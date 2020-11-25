# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a python script for demultiplexing reads from transposon sequencing data.
"""

import os, sys


inputfile_list = [r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq",
                  r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"]


#barcode_forward_sample1 = "GCCACATA"
#barcode_reverse_sample1 = "GCGAGTAA"
#barcode_forward_sample2 = "GAGCTGAA"
#barcode_reverse_sample2 = "GATAGACA"

barcode_list = ["GCCACATA",
                "GCGAGTAA",
                "GAGCTGAA",
                "GATAGACA"]

barcode_dict = {barcode_list[0]: 'sample1_forward',
                barcode_list[1]: 'sample1_reverse',
                barcode_list[2]: 'sample2_forward',
                barcode_list[3]: 'sample2_reverse'}



print("Creating files for storing output data...")
print('')
outputfile_list = []
for i in range(int(len(barcode_list)/2)):
    outputfile_list.append(os.path.splitext(inputfile_list[0])[0] + '_' + "sample" + str(i+1) + "_forward" + ".fq")
    outputfile_list.append(os.path.splitext(inputfile_list[0])[0] + '_' + "sample" + str(i+1) + "_reverse" + ".fq")
print(outputfile_list)
print('')

#CREATE FILES FOR EACH OF THE SAMPLES AND ADD THE 4 LINES CORRESPONDING TO A READ TO THE SAMPLE FILE OF THE READ.
#ALLOW INPUT AS MANY BARCODES (AS A LIST) AS WANTED BY CREATING A FOR LOOP OVER ALL INPUTTED BARCODES THAT ARE STORED IN A DICT WITH VALUES ARE THE NAMES CORRESPONDING TO THE NUMBER OF SAMPLES
#!!!CHECK IF FORWARD AND REVERSE READ COME FROM THE SAME SAMPLE



f1 = open(inputfile_list[0])
f2 = open(inputfile_list[1])


temp=0
for line1 in f1:
    if temp < 10:
        line2 = f2.readline()
        if line1.startswith('@') and line2.startswith('@') and line1.split(' ')[0] == line2.split(' ')[0]:
            print(line1)
            print(line2)
            line1 = f1.readline()
            bc1 = [barcode for barcode in barcode_dict if barcode in line1]
            if not bc1 == []:
                print(bc1[0], barcode_dict.get(bc1[0]))
            else:
                print(bc1)
    
            line2 = f2.readline()
            bc2 = [barcode for barcode in barcode_dict if barcode in line2]
            if not bc2 == []:
                print(bc2[0], barcode_dict.get(bc2[0]))
            else:
                print(bc2)
    
            print('')
            print('')
            temp+=1
    elif temp >= 10:
        break

f1.close()
f2.close()