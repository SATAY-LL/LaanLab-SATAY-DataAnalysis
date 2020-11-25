# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a python script for demultiplexing reads from transposon sequencing data.
"""

import os, sys


inputfile_list = [r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq",
                  r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"]


#KEP THE SAME NAMING FORMAT FOR THE SAMPLES (I.E. sample#_...)
barcode_dict = {"GCCACATA": 'sample1_forward',
                "GCGAGTAA": 'sample1_reverse',
                "GAGCTGAA": 'sample2_forward',
                "GATAGACA": 'sample2_reverse'}



print("Creating files for storing output data...")
print('')
outputfile_dict = {}
for barcode, samplename in barcode_dict.items():
    outputfile_dict[samplename] = os.path.splitext(inputfile_list[0])[0] + '_' + samplename + ".fq"
outputfile_dict['unmatched'] = os.path.splitext(inputfile_list[0])[0] + '_' + samplename + ".fq"
print(outputfile_dict)
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

            
            if not bc1 == [] and not bc2 == []:
                if barcode_dict.get(bc1[0]).split('_')[0] == barcode_dict.get(bc2[0]).split('_')[0]:
                    print(barcode_dict.get(bc1[0]).split('_')[0])
                else:
                    print("Read pair do not belong to same sample")
            else:
                print("Barcode not found in one of the reads")

            print('')
            print('')
            temp+=1
    elif temp >= 10:
        break

f1.close()
f2.close()