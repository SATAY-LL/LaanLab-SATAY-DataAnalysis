# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 07:43:46 2020

@author: gregoryvanbeek

This is a code to test the output of demulitplexing.py
"""

#"GCCACATA"
#"GCGAGTAA"
#"GAGCTGAA"
#"GATAGACA"

#"GCTTAGGTAGTGCTGGTCGTAGTGAGCCACATAtccgtcccgcaagttaaata"
#"AGGTCAGTCACATGGTTAGGACGCAGCGAGTAAacgaaaacgaacgggataaa"
#"GCTTAGGTAGTGCTGGTCGTAGTGAGAGCTGAAtccgtcccgcaagttaaata"
#"AGGTCAGTCACATGGTTAGGACGCAGATAGACAacgaaaacgaacgggataaa"
#
#files_list = [r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq",
#              r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"]
#
#barcode_list = ["GCCACATA",
#                "GCGAGTAA",
#                "GAGCTGAA",
#                "GATAGACA"]
#
#sample_list = [1,1,2,2]

import os
import pandas as pd

def which_barcode_in_read(files_list=[], barcode_list=[], sample_list=[], saving_dataframe=False):
    '''
    Get the barcode in all reads in the files and counts how often each barcode is present in all files.
    
    files_list: list of paths to fastq files.
    barcode_list: list of all barcodes.
    sample_list: list of samples to which the barcodes belong. This must be in the same order as the barcodes.
    '''

    if len(files_list) == 0:
        print("No input files are given. Required input: files_list[file1.fq, file2.fq, ...]")
        return
    else:
        for file in files_list:
            if not os.path.isfile(file):
                print("File %s does not exist." % file)
                return

    if len(barcode_list) == 0:
        print("No barcodes are given. Required input: barcode_list[barcode1, barcode2, ...]")
        return
    if len(sample_list) == 0:
        print("No samples are given. Required input: sample_list[sample_number1, sample_number2, ...]")
        return
    elif not len(barcode_list) == len(sample_list):
        print("barcode_list and sample_list do not have the same length. Enter just as many sample numbers as there are barcode sequences.")
        return
    


    barcode_counter_dict = {}
    for barcode in barcode_list:
        barcode_counter_dict[barcode] = 0
    del (barcode)



    barcode_sample_dict = {}
    ind = 0
    for barcode in barcode_list:
        barcode_sample_dict[barcode] = sample_list[ind]
        ind += 1
    del (barcode, sample_list, ind)



    line_counter = 0
    file_counter = 1
    read_header_list = []
    read_seq_list = []
    read_phred_list = []
    read_barcode_list = []
    read_sample_list = []
#    read_barcode_dict = {}
#    read_sample_dict = {}
    for file in files_list:
        with open(file, 'r') as f:
            print("Analyzing file number %i" % file_counter)

            for line in f:
                if line.startswith('@'):
                    line_counter += 1

                    line_header = line
                    read_header_list.append(line_header)
                    line_seq = f.readline()
                    read_seq_list.append(line_seq)
                    line_dummy = f.readline()
                    if not line_dummy[0] == '+':
                        print('WARNING: dummy line in sequence is not +.')
                        print(line_header)
                    line_phred = f.readline()
                    read_phred_list.append(line_phred)

                    barcode_counter = 0
                    for barcode in barcode_counter_dict:
                        barcode_counter += 1

                        if barcode_counter <= len(barcode_list) and barcode.upper() in line_seq:
                            barcode_counter_dict[barcode] += 1

                            read_barcode_list.append(barcode)
                            read_sample_list.append(barcode_sample_dict.get(barcode))
#                            read_barcode_dict[line_header] = barcode
#                            read_sample_dict[line_header] = barcode_sample_dict.get(barcode)
                            break
                        elif barcode_counter >= len(barcode_list):
                            read_barcode_list.append('N')
                            read_sample_list.append(0)
#                            read_barcode_dict[line_header] = 'N'

#                    if line_counter > 10:
#                        break
#                    else:
#                        print(line_header)
#                        print(line_seq)

        file_counter += 1

    reads_dict = {"header": read_header_list,
                    "sequence": read_seq_list,
                    "phred": read_phred_list,
                    "barcode": read_barcode_list,
                    "sample": read_sample_list}

    reads_df = pd.DataFrame(reads_dict, columns = [column_name for column_name in reads_dict])

    return(reads_df)
    if saving_dataframe == True:
        reads_df.to_csv(os.path.join(os.path.dirname(files_list[0]), 'reads_dataframe.csv'), index=False)
    

#    unmatched_reads = abs(sum([hits for key, hits in barcode_counter_dict.items()]) - line_counter)
#    barcode_counter_dict['unmatched_reads'] = unmatched_reads
#
#    print('')
#    print("Total number of reads found: %i" % line_counter)
#    print(barcode_counter_dict)
#
#    print('First reads:')
#    i = 0
#    for k,v in read_barcode_dict.items():
#        if i <= 5:
#            print(k)
#            print(v)
#            print('')
#        else:
#            break
#        i += 1




#barcode_sample_dict = {"GCCACATA": 1,
#                       "GCGAGTAA": 1,
#                       "GAGCTGAA": 2,
#                       "GATAGACA": 2}
#
#for read, barcode in read_barcode_dict.items():
#    if read.split(' ')[1][0] == '1':
#        barcode1 = barcode
#        barcode2 = read_barcode_dict.get(read.split(' ')[0] + ' 2' + read.split(' ')[1][1:])
#        if not barcode_sample_dict.get(barcode1) == barcode_sample_dict.get(barcode2) and not barcode2 == None:
#            print(read)
#            print(read.split(' ')[0] + ' 2' + read.split(' ')[1][1:])
#            print(barcode1)
#            print(barcode2)
#            break
#            pass
        


if __name__ == '__main__':
    reads_df = which_barcode_in_read(files_list = [r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq",
                                        r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"],
                          barcode_list = ["GCCACATA",
                                          "GCGAGTAA",
                                          "GAGCTGAA",
                                          "GATAGACA"],
                          sample_list = [1,1,2,2],
                          saving_dataframe=True)
