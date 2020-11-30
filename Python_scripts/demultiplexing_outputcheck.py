# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 07:43:46 2020

@author: gregoryvanbeek

This is a code to test the output of demulitplexing.py
"""

input_file1 = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq"
input_file2 = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"


f1 = open(input_file1,'r')
f2 = open(input_file2,'r')


barcode_sample_dict = {"GCCACATA": 1,
                       "GCGAGTAA": 1,
                       "GAGCTGAA": 2,
                       "GATAGACA": 2}


barcode_counter_dict = {"GCCACATA": 0,
                        "GCGAGTAA": 0,
                        "GAGCTGAA": 0,
                        "GATAGACA": 0}
#barcode_counter_dict = {"GCTTAGGTAGTGCTGGTCGTAGTGAGCCACATAtccgtcccgcaagttaaata":0,
#                        "AGGTCAGTCACATGGTTAGGACGCAGCGAGTAAacgaaaacgaacgggataaa":0,
#                        "GCTTAGGTAGTGCTGGTCGTAGTGAGAGCTGAAtccgtcccgcaagttaaata":0,
#                        "AGGTCAGTCACATGGTTAGGACGCAGATAGACAacgaaaacgaacgggataaa":0}

line_counter = 0
file_counter = 1
read_dict = {}
for f in [f1,f2]:
    print("Analyzing file number %i" % file_counter)
    
    for line in f:
        if line.startswith('@'):
            line_counter += 1
            
            line_header = line
            line_seq = f.readline()
            for barcode in barcode_counter_dict:
                if barcode.upper() in line_seq:
                    barcode_counter_dict[barcode] += 1
                    read_dict[line_header] = barcode
                    break
                else:
                    read_dict[line_header] = 'N'

    file_counter += 1
                
f1.close()
f2.close()


for read, barcode in read_dict.items():
    if read.split(' ')[1][0] == '1':
        barcode1 = barcode
        barcode2 = read_dict.get(read.split(' ')[0] + ' 2' + read.split(' ')[1][1:])
        if not barcode_sample_dict.get(barcode1) == barcode_sample_dict.get(barcode2) and not barcode2 == None:
            print(read)
            print(read.split(' ')[0] + ' 2' + read.split(' ')[1][1:])
            print(barcode1)
            print(barcode2)
#            break
#            pass
        


print("Total number of reads found: %i" % line_counter)
print(barcode_counter_dict)
print("Number of lines with no match found: %i" % (abs(sum([hits for key, hits in barcode_counter_dict.items()]) - line_counter)))

